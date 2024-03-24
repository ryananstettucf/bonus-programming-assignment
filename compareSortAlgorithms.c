/* 
 * Bonus assignment
 * Implement parseData, selectionSort, insertionSort, 
 * bubbleSort, mergeSort and heapSort functions. Each sorting function 
 * will count the number of extra memory allocated using the global 
 * variable extraMemoryAllocated.
 *
 * Author: Ryan Anstett
 */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// global variable to count extra memory allocated
int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

void someSortFunction(int arr[], int n) {
    int* tempArray = (int*)Alloc(n * sizeof(int));
    // Use tempArray for sorting...

    // Sorting logic here

    // Once done with tempArray, free the allocated memory
    DeAlloc(tempArray);
}


// worker function to implement heap sort
// recursively calls itself to sort the subtree
void heapWorker(int arr[], int n, int i) {
    // Initialize largest as root
    int largest = i; 

    // left = 2*i + 1
    int left = 2 * i + 1; 

    // right = 2*i + 2
    int right = 2 * i + 2; 

    // If left child is larger than root, set largest
    if (left < n && arr[left] > arr[largest])
        largest = left;

    // If right child is larger than largest so far, set largest
    if (right < n && arr[right] > arr[largest])
        largest = right;

    // If largest is not root, update root, recursively call heapWorker
    if (largest != i) {
        // Swap arr[i] with arr[largest]
        int swap = arr[i];
        arr[i] = arr[largest];
        arr[largest] = swap;

        // recursion to sort the subtree
        heapWorker(arr, n, largest);
    }
}


// implement heap sort, more efficient sort using recursion
void heapSort(int arr[], int n) {
  
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapWorker(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i > 0; i--) {
        // Move current root to end
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        // call worker function to perform heap logic
        heapWorker(arr, i, 0);
    }
}

// mergeWorker function to implement merge sort
void mergeWorker(int arr[], int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    // temporary arrays for merge
    int L[n1], R[n2];

    // Temporary arrays for merge, dynamically allocated
    //int *L = (int *)Alloc(n1 * sizeof(int));
    // int *R = (int *)Alloc(n2 * sizeof(int));


    // Copy data to temporary arrays L[] and R[]
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    // Merge the temporary arrays back into arr[l..r]
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of temp L[] array
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of temp R[] array
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

// implement mergeSort, more efficient sort using recursion
// Divide the array into two halves, sort each half, then merge them
void mergeSort(int arr[], int l, int r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        mergeWorker(arr, l, m, r);
    }
}


// implement insertion sort
// This sort does not require extra memory
void insertionSort(int* pData, int n) {
  
    int i, key, j;
  
    for (i = 1; i < n; i++) {
      
        // The element to be inserted in the sorted sequence
        key = pData[i]; 
      
      j = i - 1;

        /* Move elements of pData[0..i-1], that are
           greater than key, to one position ahead
           of their current position. This loop will
           shift elements to the right, making room
           to insert the key value in its correct
           position */
        while (j >= 0 && pData[j] > key) {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        // Insert the key value into correct position 
        pData[j + 1] = key; 
    }
}


// implement bubble sort, less efficient, slower sort algorithm
// sorts in place with no need for extra memory allocation
void bubbleSort(int* pData, int n)
{
    int i, j;

    // loop through entire array
    for (i = 0; i < n-1; i++)    
    {
        // Last i elements are already in place
        for (j = 0; j < n-i-1; j++)
        {
            if (pData[j] > pData[j+1])
            {
                // swap pData[j] and pData[j+1]
                int temp = pData[j];
                pData[j] = pData[j+1];
                pData[j+1] = temp;
            }
        }
    }
}


// implement selection sort
// selectionSort does not requre dynamic memory
void selectionSort(int* pData, int n)
{
    int i, j, minIdx;

    printf("In selectionSort \n");
  
    // One by one move boundary of unsorted subarray
    for (i = 0; i < n-1; i++)
    {    
        // Find the minimum element in unsorted array
        minIdx = i;
        for (j = i+1; j < n; j++)
        {
            if (pData[j] < pData[minIdx])
                minIdx = j;
        }

        // Swap the found minimum element with the first element
        if (minIdx != i)
        {
            int temp = pData[minIdx];
            pData[minIdx] = pData[i];
            pData[i] = temp;
        }
    }
}


int parseData(char *inputFileName, int **ppData)
{
    FILE* inFile = fopen(inputFileName, "r");
    int dataSz = 0;
    *ppData = NULL;

    if (inFile)
    {
        // Read the size of data
        if(!fscanf(inFile, "%d\n", &dataSz))
        {
          printf("Error reading data size\n");
          fclose(inFile);
          return 0;
        }
        
        // Allocate memory for the data array
        *ppData = (int *)Alloc(sizeof(int) * dataSz);

        // Check if memory allocation was successful
        if (*ppData == NULL)
        {
            fclose(inFile); // Close the file if memory allocation fails
            return 0; // Return 0 to indicate failure
        }

        // Read the actual data into the array
        for (int i = 0; i < dataSz; i++)
        {
            if (fscanf(inFile, "%d", (*ppData) + i) != 1)
            {
                // Handle the error if any data can't be read
                DeAlloc(*ppData); // Deallocate the allocated memory
                *ppData = NULL; // Set the pointer to NULL
                fclose(inFile); // Close the file
                return 0; // Return 0 to indicate failure
            }
        }

        fclose(inFile); // Close the file after reading all data
    }

    return dataSz; // Return the size of the data read
}


// prints first 100 and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	
  printf("\tData:\n\t");
  
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

// main function provided for the program
int main(void)
{
    clock_t start, end;
    int i;
    double cpu_time_used;
    char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};

    // Loop for each of the 3 input files
    for (i=0;i<3;++i)
    {
		  int *pDataSrc, *pDataCopy;
		  int dataSz = parseData(fileNames[i], &pDataSrc);

      // if empty file continue to next item in loop
		  if (dataSz <= 0)
		    continue;
		
		    pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		    printf("---LOOP %d ------------------------\n", i+1);
		
        printf("Dataset Size : %d\n",dataSz);
		    printf("---------------------------\n");

        // Start sorting algorithms
		    printf("\nSelection Sort:\n");
		    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		    extraMemoryAllocated = 0;
		    start = clock();
		    selectionSort(pDataCopy, dataSz);
		    end = clock();
		    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		    
        printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		    printArray(pDataCopy, dataSz);
  
	 	    // Start the Insertion sort
        printf("Insertion Sort:\n");
		    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		    extraMemoryAllocated = 0;
		    start = clock();
		    insertionSort(pDataCopy, dataSz);
		    end = clock();
		    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		    printArray(pDataCopy, dataSz);

        // Start the Bubble sort
		    printf("Bubble Sort:\n");
		    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		    extraMemoryAllocated = 0;
		    start = clock();
		    bubbleSort(pDataCopy, dataSz);
		    end = clock();
		    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		    printArray(pDataCopy, dataSz);
      
        // Start the Merge sort
  		  printf("Merge Sort:\n");
		    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		    extraMemoryAllocated = 0;
		    start = clock();
		    
        mergeSort(pDataCopy, 0, dataSz - 1);
		    end = clock();
		    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		    printArray(pDataCopy, dataSz);


        // Start the Heap sort
        printf("Heap Sort:\n");
		    memcpy(pDataCopy, pDataSrc, (dataSz*sizeof(int)));
		    extraMemoryAllocated = 0;
	      start = clock();
        
        // orig    heapSort(pDataCopy, 0, dataSz - 1);
		    heapSort(pDataCopy, dataSz);
        end = clock();
		    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		    printArray(pDataCopy, dataSz);
		
		    DeAlloc(pDataCopy);
		    DeAlloc(pDataSrc);
    }
    printf("done processing all 3 files\n");
    
  return 0;
}