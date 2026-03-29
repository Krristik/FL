/***************************************************************************
 *            dppopulation.h
 *
 *  Fri Mar 23 15:30:20 2012
 *  Copyright  2012  Konstantin Kozlov
 *  <kozlov@spbcas.ru>
 ****************************************************************************/

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor Boston, MA 02110-1301,  USA
 */

// Implimentation of population as a set of individs

#ifndef _DP_POPULATION_H
#define _DP_POPULATION_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include "dpindivid.h"

typedef struct DpPopulation {
	int size;
    int cur_size;
	int iter;
    int imin;  // i - individ
    int target_index;
    double dmin; // target func data
	double fmean;
	double fmax;
	int ifmax;
	int iage;
	int aage;
    int ind_size;   // individ size
    double *mean;   // M[param]
    double *variance;  // D[param]
	DpIndivid**individ;
    int *ages_descending;   // Sorted data of individ ages
    int *cost_ascending;    // Sorted data of target func values
    int nfronts;    // Pareto front number (?)
    GArray*fronts;  // Pareto front array
    int slice_a;    // Slice left border
    int slice_b;    // Slice right border

    //modifications
    double p_param;
    int xpbest_idx;
} DpPopulation;

DpPopulation*dp_population_new(int size, int initial_size, int ind_size, int targets_size, int precond_size);

void dp_population_delete(DpPopulation*pop);

void dp_population_update(DpPopulation*pop, int start_index, int end_index);

int dp_individ_ages_descending(void *p1, void *p2, void *user_data);

int dp_individ_failures_descending(void *p1, void *p2, void *user_data);

int dp_individ_cost_ascending(void *p1, void *p2, void *user_data);

void dp_population_mpi_distribute(DpPopulation*pop, int mpi_id, int mpi_nnodes);    // mpi - message passing interface

void dp_population_mpi_gather(DpPopulation*pop, int mpi_id, int mpi_nnodes);

void dp_population_save(FILE*fp, DpPopulation*pop);

void dp_population_load(FILE*fp, DpPopulation*pop);

void dp_population_cr2cost(DpPopulation*pop);   // For sophisticated target func, cr - crowd radius (?)

DpPopulation*dp_population_union(DpPopulation*population, DpPopulation*trial);

#ifdef __cplusplus
}
#endif

#endif /* _DP_POPULATION_H */
