# PDC-PROJECT
# Performance Comparison of Sorting Algorithms

# Group Members
Khuzaima Ahsan 21K-3328

Aahil Ashiq Ali 21K-4549

Khubaib Khan Lodhi 21K-4596

# Teacher Name: Sir Nadeem Kafi Khan

# OPEN MP:

# Introduction
In this project, we implemented parallel programming on sorting random numbers using three different sorting algorithms and compared their performance. 
The performance of each algorithm is subjective of the time taken to sort the data. 
We tested count sort, merge sort and selection sort algorithms.
We changed the input size on all 3 sorting algorithms and displayed the time taken on single thread and multiple threads. 

# COUNT SORT 

Counting sort is an algorithm for sorting a collection of objects according to keys.
This code can help in performance analysis of merge sort algorithm and trends from small dataset of 1024 to the length of data set the user provides. 
The counting sort algorithm does a better job when the range of the n items is within n. However, when the range is increased to n2, the performance of counting sort is worse than the comparison-based sorting algorithms. Either you can use radix sort to solve the problem or parallelize the existing code by not making count array but instead determining the index of every element in count variable and placing in output array.



