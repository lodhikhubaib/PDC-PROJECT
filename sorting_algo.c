// 		mpirun -n <number_of_nodes> <executive_file_name> <sorting_technique{Q,M,R}>
// 		<printing/not printing data{P,N}> <size_of_array>

#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#define max_rows 100000000
#define send_data_tag 560
#define recv_data_tag 1030

int array_c[max_rows];
int array_1[max_rows];
int array_m[max_rows];
void quicksort(int array_c[],int low, int high);
void mergeSort(int *alist,int n);
void radixsort(int array_c[], int n);
int  getMax(int array_c[], int n);
void countSort(int array_c[], int n, int exp);


void main(int argc, char** argv)
{

int procs0,s,send,prev,next,small_n,start,my_id,my_id0,procs,root,extra,rows,i,avg_rows,rows_to_send,row_start,row_end,row_end0,iteration_number,info[2];
char *print = argv[2];
char *p = "P";
char *n = "N";
clock_t begin0,end0;
MPI_Status status;
MPI_Request request;

srand(0);

/////// INITIALISATION OF PARALLEL NODES  /////////////
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &procs);
root=0;

////// GETTING THE PROCESSORS IN THE FORM OF 2^N //////
procs0=procs;
i=1;
while(procs0!=0)
{
extra=procs%i;
procs0=procs0/2;
i=i*2;
}
procs=procs-extra;
if(my_id<procs){////  to avoid non-working processes to get into loop

///////   ROOT NODE  /////
if(my_id==root)
{
begin0 = clock();

/////// Unused Processors //////
if (strcmp(print, p) == 0)
printf("\n%i processor(s) were unused : Provide processors of form 2^n to completely utilize them.",extra);


rows=atoi(argv[3]);
if(rows>max_rows){
printf("The number of rows exceed the permitted number of rows\n");
exit(1);
}

printf(" \n");
///////  CREATING AN ARRAY   /////
for(i=0;i<rows;i++){
array_c[i]=rand()% (100000 + 1 - 0) + 0; //rand() % (max_number + 1 - minimum_number) + minimum_number
if (strcmp(print, p) == 0)
printf("%i ",array_c[i]);
}
if (strcmp(print, p) == 0)
printf(" \n");
//////   DIVIDING AND GIVING THE ARRAY TO SEPERATE NODES /////
avg_rows=(rows)/procs;
row_end=avg_rows-1;row_end0=row_end;
extra=rows%procs;
if(extra!=0)
{row_end+=1;extra=extra-1;row_end0=row_end;}

for(i=1;i<procs;i++){
row_start=row_end+1;
row_end=row_start+avg_rows-1;
if(extra!=0)
{row_end += 1; extra=extra-1;}
rows_to_send=row_end-row_start+1;
info[0]=rows_to_send;info[1]=rows;
MPI_Send(&info, 2, MPI_INT, i, send_data_tag, MPI_COMM_WORLD);
MPI_Send(&array_c[row_start],rows_to_send,MPI_INT,i, send_data_tag, MPI_COMM_WORLD);
}
info[0]=row_end0+1;

}

else{  /////// CHILD NODES RECEIVING THEIR PART OF THE ARRAY   //////////////////////////////

MPI_Recv(&info,2,MPI_INT,root,send_data_tag,MPI_COMM_WORLD,&status);
MPI_Recv(&array_c,info[0],MPI_INT,root,send_data_tag,MPI_COMM_WORLD,&status);
}

///////////   SORTING WITHIN ALL ACTIVE NODES   ///////////////////

char *sort = argv[1];
char *q = "Q";
char *m = "M";
char *r = "R";
if (strcmp(sort, q) == 0) 
quicksort( array_c,0,info[0]-1);

if (strcmp(sort, m) == 0)
mergeSort(array_c,info[0]);

if (strcmp(sort, r) == 0)
radixsort(array_c,info[0]);


////////// MERGING OF SORTED ARRAYS ///////
my_id0=my_id;
send=1;
procs0=procs;
while(procs0!=1)
{
if(my_id0%2==1) 
{
MPI_Ssend(&info[0],1,MPI_INT,my_id-send, recv_data_tag, MPI_COMM_WORLD);
MPI_Ssend(&array_c[0],info[0],MPI_INT,my_id-send,recv_data_tag, MPI_COMM_WORLD);
procs0=1;
}
if(my_id0%2==0)
{
MPI_Recv(&info[1],1,MPI_INT,my_id+send,recv_data_tag,MPI_COMM_WORLD,&status);
MPI_Recv(&array_1[0],info[1],MPI_INT,my_id+send,recv_data_tag,MPI_COMM_WORLD,&status);
procs0=procs0/2;
send=send*2;
my_id0=my_id0/2;

int i=0;
int j=0;
int k=0;
while (i<info[0] && j<info[1])
{
if (array_c[i]<array_1[j]){
array_m[k++]=array_c[i++];
}
else{
array_m[k++]=array_1[j++];
}
}
while (i<info[0]){
array_m[k++]=array_c[i++];
}
while (j<info[1]){
array_m[k++]=array_1[j++];
}

info[0] += info[1];
for(i=0;i<info[0];i++)
{array_c[i]=array_m[i];}
}}

///////  FINAL PRINTING OF ARRAY IN THE ROOT NODE   ///////////////

if(my_id==root)
{

for(i=0;i<rows;i++)
{
if (strcmp(print, p) == 0)
printf("%i " ,array_c[i]);
}

end0 = clock();
printf("\ntime: %Lf sort:%s procs:%i  rows:%i\n\n", (long double)(end0 - begin0)/(long double)CLOCKS_PER_SEC *1000.0f,sort,procs,rows);
}

}
MPI_Finalize();
}
//////////////    QUICKSORT CALL FUNCTION    //////////////////////

void quicksort(int array_c[],int low, int high)
{
int temp,pivot,i,j,pi;
if (low<high)
{
pivot = array_c[high];
i=low-1;
for(j=low; j<high; j++)
{
if(array_c[j]<=pivot)
{
i++;
temp=array_c[i];
array_c[i]=array_c[j];
array_c[j]=temp;
}
}
temp=array_c[i+1];
array_c[i+1]=array_c[high];
array_c[high]=temp;
pi=i+1;
quicksort(array_c,low,pi-1);
quicksort(array_c,pi+1,high);
}}
//////////////////////////////    MERGESORT CALL FUNCTION    //////////////////////
void  mergeSort(int *alist,int n)
{

if (n>1){
  int mid;
      int *lefthalf;
      int *righthalf;
      mid=n/2;

      lefthalf=(int*)malloc(mid*sizeof(int));
      righthalf=(int*)malloc((n-mid)*sizeof(int));

      for (int m=0;m<mid;m++){
        lefthalf[m]=alist[m];

          }
      for (int i=mid;i<n;i++){
        righthalf[i-mid]=alist[i];
      }
      mergeSort(lefthalf,mid);
      mergeSort(righthalf,n-mid);

        int i=0;
        int j=0;
        int k=0;

        int left_size,right_size;
        left_size=mid;
        right_size=n-mid;

        while (i<left_size && j<right_size)
           {
             if (lefthalf[i]<righthalf[j]){
                alist[k++]=lefthalf[i++];
 
           }
             else{
               alist[k++]=righthalf[j++];
               }

       }
        while (i<left_size){
            alist[k++]=lefthalf[i++];

          }
        while (j<right_size){
            alist[k++]=righthalf[j++];

          }

   }

}
//////////////////////////////    RADIXSORT CALL FUNCTION    //////////////////////

void radixsort(int array_c[], int n)
{

    int m = getMax(array_c, n);

    for (int exp = 1; m/exp > 0; exp *= 10)
        countSort(array_c, n, exp);
}

///////////////  Counting Sort   ////////////

void countSort(int array_c[], int n, int exp)
{
    int output[n]; 
    int i, count[10] = {0};

    
    for (i = 0; i < n; i++)
        count[ (array_c[i]/exp)%10 ]++;

  
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (i = n - 1; i >= 0; i--)
    {
        output[count[ (array_c[i]/exp)%10 ] - 1] = array_c[i];
        count[ (array_c[i]/exp)%10 ]--;
    }

    // Copy the output array to arr[], so that arr[] now
    for (i = 0; i < n; i++)
        array_c[i] = output[i];
}


// A function to get maximum value in arr[]
int getMax(int array_c[], int n)
{
    int mx = array_c[0];
    for (int i = 1; i < n; i++)
        if (array_c[i] > mx)
            mx = array_c[i];
    return mx;
}

