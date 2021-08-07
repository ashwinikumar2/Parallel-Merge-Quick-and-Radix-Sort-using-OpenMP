#include <iostream>
#include <omp.h>
#include "sort.h"
#include <vector>
#include <cmath>
// #include <unistd.h>
// #define THREAD_NUM 4
using namespace std;
#define THREAD_NUM 200
void pSort::init()
{}
   
void pSort::close()
{}

///////////////////////////////////////////////////////////////////////////////////////MERGE SORT////////////////////////////////////////////////////////////////////
pSort::dataType* merge(pSort::dataType *arr, pSort::dataType *arr1 , int l, int m, int r)                       //used to merge two arrays by comparing the first elements(in iteration) of left portion as compared to right portion
{                                                                                                               //or we can say the usual merge used in serial merge sort algorithm    

    for(int i = 0; i < m-l+1; i++)
    {    
        arr1[i] = arr[l + i];
    }
    for(int j = 0; j < r-m; j++)
    {
        arr1[j+m-l+1] = arr[m + 1 + j];
    }
    int i = 0,j=0,k=l; 
    
    while (i < m - l + 1 && j < r - m)                  //while we are iterating in the left portion of arr i.e. from 0 to m-l and iterating in the right portion of arr i.e. m-l+1 to r-m
    {
        if (arr1[i].key <= arr1[j+m-l+1].key)                       
        {
            arr[k] = arr1[i];
            i=i+1;
        }
        else
        {
            arr[k] = arr1[j+m-l+1];
            j=j+1;
        }
        k=k+1;
    }

    while (i < m - l + 1) 
    {
        arr[k] = arr1[i];
        i=i+1;
        k=k+1;
    }
    while (j < r - m)
    {
        arr[k] = arr1[j+m-l+1];
        j=j+1;
        k=k+1;
    }

    return arr;                                         
}
void mergeSort(pSort::dataType *arr, int l, int r)                                                  //serial merge sort used for doing merge sort serially on an array by the usual merge sort algorithm
{
    if (l < r)
    {   pSort:: dataType *temp=arr;

        int m = (l+r) / 2;
        pSort::dataType *arr1=new pSort::dataType[r-l+1];
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        pSort::dataType *output=merge(arr,arr1, l, m, r);
        for (int i = l; i < r; i++)                                                                 //used for setting the pointer values of what was obtained as data at initialization of this function
        {
            /* code */
            temp[i]=output[i];
        }
    }
}

void parallel_merge(pSort::dataType *data, int ndata, int initial, int last)                                               //parallel version of merge sort algorithm 
{
    

    if(last-initial+1==2)                                       //if size of input data is 2
    {
        pSort::dataType *temp=data;
        if(data[initial].key>data[last].key)                //if second element is less than first element in the array then swap them
        {

            temp[initial]=data[last];
            temp[last]=data[initial];
            data[last]=temp[last];
            data[initial]=temp[initial];
            return;
        }
        return;
    }
    if(last<=initial){                                          //if size of input data is 1 then return input array as it is
        
        return;
    }
    omp_set_num_threads(2);                                     // create 2 threads for this recursive call of parallel merge
    #pragma omp parallel                                        //start parallel computation using threads
    {
        #pragma omp for                                         //create parallel for loop
        for (int i = 0; i < 2; i++)
        {
            if(i==0)                                            //if it is the first thread then call recursively on left half
            {
                parallel_merge(data, ndata, initial, (last+initial)/2); 
            }
            else                                                //else call recursively on right half
            {
                parallel_merge(data, ndata, (last+initial)/2+1, last);
            }
        }
    }

    pSort::dataType *arr1=new pSort::dataType[last-initial+1];          //a dummy array for passing in merge function
    data=merge(data, arr1,initial, (last+initial)/2, last);             //merge left and right sorted data by comparing each element of left and right subarray and then merging them into a big array
    
}
////////////////////////////////////////////////////////////////////////////////////////MERGE SORT///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////QUICK SORT////////////////////////////////////////////////////////////////////
pSort::dataType* sequential_quick_sort(pSort::dataType *data, int ndata, int low, int high)                         //serial quick sort implementing the usual quick sort algorithm by finding the partition index and 
{                                                                                                       //then recursively use this on left and right portions seperated by partition index
    if(low<high)
    {
        // pSort::dataType *temp=data;

        pSort::dataType pivot = data[high];                         //taking pivot to be the value of key of last element in the recieved array  
        int i = (low - 1); 
      
        for (int j=low;j<high;j++)  
        {  
            if (data[j].key < pivot.key)  
            {  
                i++; 
                pSort::dataType x=data[i];
                data[i]=data[j];                                    //swapping values of array elements at indexes i and j
                data[j]=x; 
            }  
        }  
        pSort::dataType x=data[i+1];
        data[i+1]=data[high];                                       //swapping values of array elements at indexes i+1 and high
        data[high]=x;

        int pi=i+1;
        data= sequential_quick_sort(data, ndata, low, pi-1);
        data=sequential_quick_sort(data, ndata, pi+1, high);
        return data;
        
    }
    return data;
}


pSort::dataType* quick_sort(pSort::dataType *data, int ndata, int low, int high)
{
    if(high-low+1<=1<<25){                                 //if size of quick sort is less than 10 then use sequential quick sort
        data=sequential_quick_sort(data, ndata, low, high);
        return data;
    }
    if(low<high)                                        //if low is less than high then only perform the computation for quick sort
    {
        // pSort::dataType *temp=data;

        pSort::dataType pivot = data[high];                         //taking pivot to be the value of key of last element in the recieved array  
        int i = (low - 1); 
      
        for (int j=low;j<high;j++)  
        {  
            if (data[j].key < pivot.key)  
            {  
                i++; 
                pSort::dataType x=data[i];
                data[i]=data[j];                                    //swapping values of array elements at indexes i and j
                data[j]=x; 
            }  
        }  
        pSort::dataType x=data[i+1];
        data[i+1]=data[high];                                       //swapping values of array elements at indexes i+1 and high
        data[high]=x;

        int pi=i+1;
        #pragma omp parallel sections                                      //start parallel computation using sections
        {
                #pragma omp section                                                //if it is the first section then call recursively on left subarray
                {
                    data=quick_sort(data, ndata, low, pi-1);
                }
                #pragma omp section                                                    //else call recursively on right subarray
                {
                    data=quick_sort(data, ndata, pi+1, high);
                }
        }
        return data;
        
    }
    return data;


}

///////////////////////////////////////////////////////////////////////////////////////QUICK SORT////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////RADIX SORT////////////////////////////////////////////////////////////////////
int number_of_digits(int n)                                     //returns the number of digits in an integer(without sign)
{
    if(n>0)
    return floor(log10(n)+1);
    else if (n==0)
        return 1;
    else{                       //if integer is negative then also return the digit without any sign
        int x= floor(log10(-n)+1);
        if(x<0)
            return (-x);
        else
            return x;
    }   
}
pSort::dataType* radix_sort(pSort::dataType *data, int ndata)                               //parallel version of radix sort
{
    for(int i=1;i<=15;i++)                                                                  //for each number of digit
    {   

        vector<vector<vector<pSort::dataType>>> buckets;        //a 3-D matrix such that digit->thread_bucket->bucket
        for(int j=0;j<10;j++)                                   //for each digit from 0 to 9 initialize a vector into the bucket matrix
        {
            vector<vector<pSort::dataType>> v;
            buckets.push_back(v);
        }
        int processes=omp_get_num_procs();                      //to get number of processors available
        if(processes<1)
            processes=1;
        int num=1;
        for(int k=0;k<10;k++)                                   //for each digit from 0 to 9
        {
            for(int j=0;j<processes;j++)                        //for each thread from 0 to processes-1
            {
                vector<pSort::dataType> v;
                ((buckets[k])).push_back(v);


            }
        }
        omp_set_num_threads(processes);                         // set number of threads in "parallel" blocks
        #pragma omp parallel                                    //start parallel computation using threads
        {
            #pragma omp for                                     //create parallel for loop
            for (int thread = 0; thread < processes; thread++)
            {
                vector<vector<pSort::dataType>> m;              //creating local buckets for this thread
                for(int j=0;j<10;j++)                           //for each digit from 0 to 9 initialize a vector into the local buckets matrix
                {
                    vector<pSort::dataType> v;
                    m.push_back(v);
                }


                if(thread!=(processes-1))                       //if it is not the last thread 
                {   
                    int initial=thread*(ndata/processes);       //initial index for the share of memory for this thread from shared memory
                    int last=(thread+1)*(ndata/processes)-1;    //last index for the share of memory for this thread from shared memory
                    for(int j=initial;j<=last;j++)              //for all the data between initial and last indexes i.e. the share of memory for this thread from shared memory
                    {
                        int y=data[j].key;          
                        if(number_of_digits(y)>=i)              //if number of digits greater than current digit iteration i.e. whether it is more than i (e.g. LSD)
                        {
                            y=y/((int)(pow(10, i-1)));
                            if(y<0)
                                y=-y;
                            int x=(y)%(10);
                            (m[x]).push_back(data[j]);
                        }
                        else                                    //if not then push the number into bucket of 0 in local buckets matrix
                        {
                            int x=0;
                            (m[x]).push_back(data[j]);
                        }

                    }
                    for(int i1=0;i1<10;i1++){                   //copying all the local buckets into shared memory buckets collection
                    
                        (buckets[i1])[thread]=m[i1];            
                    }
                    

                }
                else                                            //if it is last thread
                {
                    int initial=thread*(ndata/processes);       //initial index for the share of memory for this thread from shared memory
                    int last=ndata-1;                           //last index for the share of memory for this thread from shared memory
                    for(int j=initial;j<=last;j++)              //for all the data between initial and last indexes i.e. the share of memory for this thread from shared memory
                    {

                        
                        int y=data[j].key;
                        if(number_of_digits(y)>=i)              //if number of digits greater than current digit iteration i.e. whether it is more than i (e.g. LSD)
                        {
                            y=y/((int)(pow(10, i-1)));
                            if(y<0)
                                y=-y;
                            int x=(y)%(10);
                            (m[x]).push_back(data[j]);
                        }
                        else                                    //if not then push the number into bucket of 0 in local buckets matrix
                        {
                            int x=0;
                            (m[x]).push_back(data[j]);
                        }
                    }
                    for(int i1=0;i1<10;i1++){                   //copying all the local buckets into shared memory buckets collection
                    
                        (buckets[i1])[thread]=m[i1];
                        
                    }

                }
                /* code */
            }
        }
        pSort::dataType *data_new=new pSort::dataType[ndata];           //creating a dummy array to keep it as the new data array
        int index_of_data_new=0;                                        //index of data new= 0;
        for(int k=0;k<buckets.size();k++)                               //for each digit bucket in the 3-D matrix
        {
            for(int j=0;j<(buckets[k]).size();j++)                      //for each thread
            {
                for(int p=0;p<buckets[k][j].size();p++){                //for each thread's buckets of this digit 'k'
                
                    data_new[index_of_data_new]=((buckets[k])[j])[p];
                    index_of_data_new++;
                }

            }
        }
        data=data_new;                                                  //data becomes data new
    }

    //now for seperating negative and positive numbers

    vector<vector<vector<pSort::dataType>>> buckets;        //a 3-D matrix such that digit->thread_bucket->bucket
    for(int j=0;j<2;j++)                                    //0 for negative and 1 for positive
    {
        vector<vector<pSort::dataType>> v;
        buckets.push_back(v);
    }
    int processes=omp_get_num_procs();                      //to get number of processors available
    if(processes<1)
        processes=1;
    int num=1;
    for(int k=0;k<2;k++)                                    //for both 0 and 1
    {
        for(int j=0;j<processes;j++)                        //for each thread from 0 to processes-1
        {
            vector<pSort::dataType> v;
            ((buckets[k])).push_back(v);


        }
    }
    omp_set_num_threads(processes);                         // set number of threads in "parallel" blocks
    int negative=0,positive=0;                              //number of negative numbers and number of positive numbers
    #pragma omp parallel
    {
        #pragma omp for                                     //create parallel for loop
        for (int thread = 0; thread < processes; thread++)
        {

            vector<vector<pSort::dataType>> m;              //creating local buckets for this thread
            for(int j=0;j<2;j++)                            //for both 0 and 1 initialize a vector into the local buckets matrix
            {
                vector<pSort::dataType> v;
                m.push_back(v);
            }


            if(thread!=(processes-1))                       //if it is not the last thread 
            {   
                int initial=thread*(ndata/processes);       //initial index for the share of memory for this thread from shared memory
                int last=(thread+1)*(ndata/processes)-1;    //last index for the share of memory for this thread from shared memory
                for(int j=initial;j<=last;j++)              //for all the data between initial and last indexes i.e. the share of memory for this thread from shared memory
                {

                    pSort::dataType y=data[j];
                    if(y.key<0){                            //if the number is negative then put in 0th bucket
                        (m[0]).push_back(data[j]);
                    }
                    else                                    //else put it in 1st bucket
                    {
                        m[1].push_back(data[j]);
                    }

                }
                negative+=m[0].size();                      //increase the number of negative ints by 1
                positive+=m[1].size();                      //increase the number of positive ints by 1
                for(int i1=0;i1<2;i1++){                    //copy local buckets to shared memory buckets
                
                (buckets[i1])[thread]=m[i1];  
                }
                

            }
            else                                             //if it is last thread
            {
                int initial=thread*(ndata/processes);       //initial index for the share of memory for this thread from shared memory
                int last=ndata-1;                           //last index for the share of memory for this thread from shared memory
                for(int j=initial;j<=last;j++)              //for all the data between initial and last indexes i.e. the share of memory for this thread from shared memory
                {

                    
                    pSort::dataType y=data[j];
                    if(y.key<0){                            //if the number is negative then put in 0th bucket
                        (m[0]).push_back(data[j]);
                    }
                    else                                    //else put it in 1st bucket
                    {
                        m[1].push_back(data[j]);
                    }
                }
                negative+=m[0].size();                      //increase the number of negative ints by 1
                positive+=m[1].size();                      //increase the number of positive ints by 1
                for(int i1=0;i1<2;i1++){                    //copy local buckets to shared memory buckets
                
                (buckets[i1])[thread]=m[i1];
                    
                }

            }
            /* code */
        }
    }
    pSort::dataType *data_new=new pSort::dataType[ndata];       //creating a dummy array to keep it as the new data array
    int index_of_data_new=0;                                    //index of data new= 0;
    for(int k=0;k<buckets.size();k++)                           //for each digit bucket in the 3-D matrix
    {
        for(int j=0;j<(buckets[k]).size();j++)                   //for each thread
        {
            for(int p=0;p<buckets[k][j].size();p++){            //for each thread's buckets of this digit 'k'
            
                data_new[index_of_data_new]=((buckets[k])[j])[p];
                index_of_data_new++;
            }

        }
    }
    for(int j=0;j<negative;j++)                             //copy the negative ints obtained in reverse order
    {
        data[j]=data_new[negative-1-j];
    }
    for(int j=0;j<negative;j++)                             //copy the positive ints in the same order
    {
        data[j]=data_new[negative+j];
    }
    return data;
}

void pSort::sort(dataType *data, int ndata, SortType sorter)
{
    pSort::dataType *temp=data;
    if(sorter== BEST)
        data=quick_sort(data, ndata, 0, ndata-1);
    else if (sorter== QUICK)
        data=quick_sort(data, ndata, 0, ndata-1);
    else if(sorter == MERGE)
        parallel_merge(data, ndata, 0, ndata-1);
    else 
    data=radix_sort(data, ndata);  

    for(int i=0;i<ndata;i++)
    {
        temp[i]=data[i];
    }
    // for(int i=0;i<ndata;i++)
    // {
        // cout<<"data["<<i<<"]= "<<data[i].key<<endl;
    // }
}

///////////////////////////////////////////////////////////////////////////////////////RADIX SORT////////////////////////////////////////////////////////////////////