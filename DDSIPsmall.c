/*  Authors:            Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
	Description:
	A few small functions for signal handling, cpu time, minimum, maximum,
	feasibility test, and error handling w.r.t. cplex errors

	License:
	This file is part of DDSIP.

    DDSIP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    DDSIP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDSIP; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#define THRESHOLD 5
//int DDSIP_Partition (double *, int *, int, int, int);
//void DDSIP_QuickSortRe (double *, int *, int, int, int);

// DDSIP_Quicksort (descending) combined with insertion sort - for an indexed array, sorting just the index
void DDSIP_qsort_ins_D (const double * a, int * ind, int l, int r)
{
    int i, j, k;
    double tmp;
    if (r-l > THRESHOLD)  //Quicksort
    {
        i=l-1;
        j=r;
        for(;;)
        {
            k = ind[r];
            while (i<r && a[ind[++i]]>a[k]);
            while (j>l && a[ind[--j]]<a[k]);
            if (i>=j) break;
            DDSIP_SWAP (ind[i], ind[j]);
        }
        DDSIP_SWAP (ind[i], ind[r]);

        DDSIP_qsort_ins_D (a, ind, l, i-1);
        DDSIP_qsort_ins_D (a, ind, i+1, r);
    }
    else  //insertion sort
    {
        for(i=l+1; i<=r; ++i)
        {
            k = ind[i];
            tmp=a[k];
            for(j=i-1; j>=l && tmp>a[ind[j]]; --j)
                ind[j+1]=ind[j];
            ind[j+1]=k;
        }
    }
}
// Quicksort (ascending) combined with insertion sort - for an indexed array, sorting just the index
void DDSIP_qsort_ins_A (const double * a, int * ind, int l, int r)
{
    int i, j, k;
    double tmp;
    if (r-l > THRESHOLD)  //Quicksort
    {
        i=l-1;
        j=r;
        for(;;)
        {
            k = ind[r];
            while (i<r && a[ind[++i]]<a[k]);
            while (j>l && a[ind[--j]]>a[k]);
            if (i>=j) break;
            DDSIP_SWAP (ind[i], ind[j]);
        }
        DDSIP_SWAP (ind[i], ind[r]);

        DDSIP_qsort_ins_A (a, ind, l, i-1);
        DDSIP_qsort_ins_A (a, ind, i+1, r);
    }
    else  //insertion sort
    {
        for(i=l+1; i<=r; ++i)
        {
            k = ind[i];
            tmp=a[k];
            for(j=i-1; j>=l && tmp<a[ind[j]]; --j)
                ind[j+1]=ind[j];
            ind[j+1]=k;
        }
    }
}
// Quicksort and Partition from
// http://ciips.ee.uwa.edu.au/~morris/Year2/PLDS210/qsort.html
// Modified to sort array and indices
//==========================================================================
//int
//Partition (double *a, int *ind, int low, int high, int length)
//{
//  int left, right, pivot;
//  double pivot_item, tmp;
//  pivot_item = a[low];
//  pivot = low;
//  left = low;
//  right = high;
//  while (left < right) {
//    /* Move left while item < pivot */
//    while (left < length && a[left] <= pivot_item)
//      left++;
//    /* Move right while item > pivot */
//    while (right && a[right] > pivot_item)
//      right--;
//    if (left < right) {
//      DDSIP_SWAP (a[left], a[right]);
//      DDSIP_SWAP (ind[left], ind[right]);
//    }
//  }
//  /* right is final position for the pivot */
//  a[low] = a[right];
//  a[right] = pivot_item;
//  DDSIP_SWAP (ind[low], ind[right]);
//  return right;
//}
//
////==========================================================================
//void
//QuickSortRe (double *val, int *ind, int low, int high, int length)
//{
//  int pivot;
//  /* Termination condition! */
//  if (high > low) {
//    pivot = DDSIP_Partition (val, ind, low, high, length);
//    DDSIP_QuickSortRe (val, ind, low, pivot - 1, length);
//    DDSIP_QuickSortRe (val, ind, pivot + 1, high, length);
//  }
//}
//
////==========================================================================
//// For sorting the whole arrays, this one is shorter
//void
//QuickSort (double *val, int *ind, int length)
//{
//  DDSIP_QuickSortRe (val, ind, 0, length - 1, length);
//}

//==========================================================================
// Signal handling
void
DDSIP_HandleSignal (int signal_number)
{
    DDSIP_killsignal = signal_number;
    printf ("\n");
}

//==========================================================================
// Register for signal handling
void
DDSIP_RegisterSignalHandlers (void)
{
    void (*error) (int);

    error = signal (SIGINT, DDSIP_HandleSignal);
    if (error == SIG_ERR)
        fprintf (stderr, "*Warning: Failed to register handler for 'SIGINT'!");
    if (error == SIG_IGN)
        signal (SIGINT, SIG_IGN);
}

//==========================================================================
// returns cpu-time in seconds
double
DDSIP_GetCpuTime (void)
{
#ifdef _WIN32
    clock_t now;
#else
    struct tms now;
#endif
    int sec, hun;
    double total;

#ifndef CLK_TCK
    long clocks_per_second = 0;

    clocks_per_second = sysconf (_SC_CLK_TCK);
#ifdef _WIN32
    now = clock ();
#else
    times (&now);
    hun = ((now.tms_utime % clocks_per_second) * 100) / clocks_per_second;
    sec = (now.tms_utime / clocks_per_second);
#endif

#else
#ifdef _WIN32
    now = clock ();
#else
    times (&now);
    hun = ((now.tms_utime % CLK_TCK) * 100) / CLK_TCK;
    sec = (now.tms_utime / CLK_TCK);
#endif

#endif

#ifdef _WIN32
    total = now / CLK_TCK;
#else
    total = (double) sec + ((double) (hun) / 100.0);
#endif

    return total;
}

//==========================================================================
// The CPLEX errors below indicate infeasibility during presolve. We take care of those separately
int
DDSIP_Error (int stat)
{
    if ((stat == 0) || (stat == CPXERR_PRESLV_INForUNBD) || (stat == CPXMIP_OPTIMAL))
        return 0;
    else
        return 1;
}

//==========================================================================
// Returns 1 if status implies some kind of infeasibility, 0 otherwise
// CPLEX 8.0 (103, 106, 108, 110, 112, 114, 117, 119) // 115, 118
int
DDSIP_NoSolution (int stat)
{
//war auch mal drin: (stat == CPXMIP_OPTIMAL_INFEAS) || und das auskommentiert: (stat == CPXMIP_UNBOUNDED) || (stat == CPXMIP_INForUNBD)
    if ((stat == CPXMIP_INFEASIBLE) || (stat == CPXERR_NO_SOLN) || (stat == CPXERR_NO_SOLNPOOL) || (stat == CPXERR_PRESLV_INF)
            || (stat == CPXMIP_NODE_LIM_INFEAS) || (stat == CPXMIP_TIME_LIM_INFEAS) || (stat == CPXMIP_FAIL_INFEAS)
            || (stat == CPXMIP_MEM_LIM_INFEAS) || (stat == CPXMIP_ABORT_INFEAS) || (stat == CPXMIP_FAIL_INFEAS_NO_TREE)
            || (stat == CPXMIP_UNBOUNDED) || (stat == CPXMIP_INForUNBD))
        return 1;
    else
        return 0;
}

//==========================================================================
// The CPLEX errors below indicate infeasibility during presolve. We take care of those separately
int
DDSIP_Infeasible (int stat)
{
    if ((stat == CPX_STAT_INFEASIBLE) || (stat == CPXERR_PRESLV_INForUNBD)
            || (stat == CPXERR_SUBPROB_SOLVE) || (stat == CPXMIP_INForUNBD) || (stat == CPXMIP_INFEASIBLE))
        return 1;
    else
        return 0;
}
//==========================================================================
// translate time to more readable form
void
DDSIP_translate_time (double total, int * hours, int * mins, double * secs)
{
    *hours = floor (total/3600.0);
    total -= 3600.0* (*hours);
    *mins  = floor (total/60.0);
    *secs  = total - 60.0* (*mins);
    return;
}
