#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);
   /* Stop the timer */

   /* Add you code here  */
    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1) * (n - 1) / p;
    low_value = low_value + (low_value + 1) % 2;
    high_value = high_value - (high_value + 1) % 2;
    size = (high_value - low_value) / 2 + 1;
    
    /**
     * process 0 must holds all primes used
     */
    proc0_size = (n - 1) / p;
    if ((2 + proc0_size) < (int) sqrt((double) n))
    {
        if (!id)
            printf("Too many processes.\n");
        MPI_Finalize();
        exit(1);
    }

    /**
     * Allocation
     */
    marked = (char*) malloc(size);
    if (marked == NULL)
    {
        printf(" Cannot allocate enough memory.\n");
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size; i++)
        marked[i] = 0;

    /**
     * Core Function
     */
    if (!id)
        index = 0;
    prime = 3;
    
    do
    {
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2;
        else
        {
            if ((low_value % prime) == 0)
                first = 0;
            else
                first = (prime - (low_value % prime) + low_value / prime % 2 * prime) / 2;
        }
        for (i = first; i < size; i += prime)
            marked[i] = 1;
        if (!id)
        {
            while(marked[++index] == 1);
            prime = 3 + index * 2;
        }
        if (p > 1)
            MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++)
        if (marked[i] == 0)
            count++;
    if (id == 0)
        count++;    
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   










   ///* Stop the timer */
   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

