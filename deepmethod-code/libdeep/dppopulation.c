/***************************************************************************
 *            dppopulation.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef MPIZE
#include <mpi.h>
#endif

#include "dpindivid.h"
#include "dppopulation.h"

DpPopulation*dp_population_new(int size, int initial_size, int ind_size, int targets_size, int precond_size)
{
	int i;
	DpPopulation*pop;
	pop = (DpPopulation*)malloc(sizeof(DpPopulation));
	pop->size = size;
	pop->cur_size = initial_size;
	pop->individ = (DpIndivid**)calloc(size, sizeof(DpIndivid*));
	pop->mean = (double*)calloc(ind_size, sizeof(double));
	pop->variance = (double*)calloc(ind_size, sizeof(double));
	pop->cost_ascending = (int*)calloc(size, sizeof(int));
	pop->ages_descending = (int*)calloc(size, sizeof(int));
	pop->ind_size = ind_size;
	for ( i = 0; i < size; i++ ) {
		pop->individ[i] = dp_individ_new(ind_size, targets_size, precond_size);
		pop->cost_ascending[i] = i;
		pop->ages_descending[i] = i;
	}
	pop->imin = 0;
	pop->dmin = -1;
	pop->iage = 0;
	pop->aage = 0;
	pop->iter = 0;
	pop->nfronts = 0;
	pop->fronts = NULL;
	pop->target_index = 0;
	pop->slice_a = 0;
	pop->slice_b = initial_size;
	pop->p_param = 1;
	pop->xpbest_idx = 0;
#ifdef MPIZE
/* MPI initialization steps */
	int world_id = 0, world_count = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &world_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_id);
	int ind_per_node = (int)ceil(pop->size / world_count);
	int ind_per_last_node = pop->size - ind_per_node * (world_count - 1);
	pop->slice_a = world_id * ind_per_node;
	pop->slice_b = (world_id == world_count - 1) ? size : pop->slice_a + ind_per_node;
#endif
	return pop;
}

DpPopulation*dp_population_union(DpPopulation*population, DpPopulation*trial)
{
	int i;
	DpPopulation*pop;
	pop = (DpPopulation*)malloc(sizeof(DpPopulation));
	pop->size = population->size + trial->size;
	pop->ind_size = population->ind_size;
	pop->individ = (DpIndivid**)calloc(pop->size, sizeof(DpIndivid*));
	pop->mean = (double*)calloc(pop->ind_size, sizeof(double));
	pop->variance = (double*)calloc(pop->ind_size, sizeof(double));
	pop->cost_ascending = (int*)calloc(pop->size, sizeof(int));
	pop->ages_descending = (int*)calloc(pop->size, sizeof(int));
	for ( i = 0; i < population->size; i++ ) {
		pop->individ[i] = population->individ[i];
		pop->cost_ascending[i] = i;
		pop->ages_descending[i] = i;
	}
	for ( i = 0; i < trial->size; i++ ) {
		pop->individ[population->size + i] = trial->individ[i];
		pop->cost_ascending[i] = population->size + i;
		pop->ages_descending[i] = population->size + i;
	}
	pop->imin = 0;
	pop->dmin = -1;
	pop->iage = 0;
	pop->aage = 0;
	pop->iter = 0;
	pop->nfronts = 0;
	pop->fronts = NULL;
	pop->target_index = 0;
	pop->p_param = 1;
	pop->xpbest_idx = 0;
	return pop;
}

void dp_population_save(FILE*fp, DpPopulation*pop)
{
	int i;
	for ( i = 0; i < pop->size; i++ ) {
	  dp_individ_save(fp, pop->individ[i]);
	  fprintf(fp, "%d\n", pop->cost_ascending[i]);
	  fprintf(fp, "%d\n", pop->ages_descending[i]);
	}
	for ( i = 0; i < pop->ind_size; i++ ) {
		fprintf(fp, "%13.9f ", pop->mean[i]);
	}
	fprintf(fp, "\n");
	for ( i = 0; i < pop->ind_size; i++ ) {
		fprintf(fp, "%13.9f ", pop->variance[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "%d\n", pop->imin);
	fprintf(fp, "%13.9f\n", pop->dmin);
	fprintf(fp, "%d\n", pop->iage);
	fprintf(fp, "%d\n", pop->aage);
	fprintf(fp, "%d\n", pop->iter);
}

void dp_population_load(FILE*fp, DpPopulation*pop)
{
	int i;
	for ( i = 0; i < pop->size; i++ ) {
	  dp_individ_load(fp, pop->individ[i]);
	  fscanf(fp, "%d", &(pop->cost_ascending[i]));
	  fscanf(fp, "%d", &(pop->ages_descending[i]));
	}
	for ( i = 0; i < pop->ind_size; i++ ) {
		fscanf(fp, "%lf", &(pop->mean[i]));
	}
	for ( i = 0; i < pop->ind_size; i++ ) {
		fscanf(fp, "%lf", &(pop->variance[i]));
	}
	fscanf(fp, "%d", &(pop->imin));
	fscanf(fp, "%lf", &(pop->dmin));
	fscanf(fp, "%d", &(pop->iage));
	fscanf(fp, "%d", &(pop->aage));
	fscanf(fp, "%d", &(pop->iter));
}

void dp_population_delete(DpPopulation*pop)
{
	int i;
	for ( i = 0; i < pop->size; i++ ) {
		dp_individ_delete(pop->individ[i]);
	}
	free(pop->individ);
	free(pop->mean);
	free(pop->variance);
	free(pop);
}

void dp_population_update(DpPopulation*pop, int start_index, int end_index)
{
	int i, j, n, kounter;
	double smean, smean2, s;
	double fmin, fmean, fmax;
	int amax, imax, imin, ifmax;
	int pbest_idx;

	//xpbest - p -> p*n -> cost_ascending[p*n]

	g_qsort_with_data(pop->ages_descending, pop->cur_size, sizeof(pop->ages_descending[0]), (GCompareDataFunc)dp_individ_ages_descending, pop);
	g_qsort_with_data(pop->cost_ascending, pop->cur_size, sizeof(pop->cost_ascending[0]), (GCompareDataFunc)dp_individ_cost_ascending, pop);

	kounter = 0;
	fmean = 0;
	fmax = 0;
	ifmax = -1;
	if (pop->individ[start_index]->cost < G_MAXDOUBLE) {
		fmean += pop->individ[start_index]->cost;
		kounter++;
		if ( pop->individ[start_index]->cost > fmax ) {
			fmax = pop->individ[start_index]->cost;
			ifmax = start_index;
		}
	}
	fmin = pop->individ[start_index]->cost;
	pop->individ[start_index]->cost_ind = pop->cost_ascending[start_index];
	imin = start_index;
	amax = pop->individ[start_index]->age;
	imax = start_index;
	n = end_index - start_index + 1;
	for ( j = start_index + 1; j < end_index; j++) {
		pop->individ[pop->cost_ascending[j]]->cost_ind = j;
		if (pop->individ[j]->cost < G_MAXDOUBLE) {
			fmean += pop->individ[j]->cost;
			kounter++;
			if ( pop->individ[j]->cost > fmax ) {
				fmax = pop->individ[j]->cost;
				ifmax = j;
			}
		}
		if ( pop->individ[j]->cost < fmin ) {
			fmin = pop->individ[j]->cost;
			imin = j;
		}
	}
	if ( imax == imin )
		imax = end_index - 1;
	for ( j = start_index + 1; j < end_index; j++) {
		if ( j != imin && pop->individ[j]->age > amax ) {
			amax = pop->individ[j]->age;
			imax = j;
		}
	}

	//xamean - аналогично, только добавить веса
	for ( i = 0; i < pop->ind_size; i++) {
		smean = 0;
		smean2 = 0;
		for ( j = start_index; j < end_index; j++) {
			s = pop->individ[j]->x[i];
			smean += s;
			smean2 += s * s;
		}
		smean /= n;
		smean2 /= n;
		pop->mean[i] = smean;
		pop->variance[i] = ( smean2 - smean * smean ) / ( n - 1 ) ;
	}
	pop->imin = imin;
	pop->dmin = fmin;
	pop->fmean = fmean / kounter;
	pop->iage = imax;
	pop->aage = amax;
	pop->fmax = fmax;
	pop->ifmax = ifmax;
	//xpbest modification
	pbest_idx = floor(pop->p_param * pop->size);
	pop->xpbest_idx = pop->cost_ascending[pbest_idx];
}

int dp_population_target_compare(gconstpointer a, gconstpointer b, gpointer user_data)
{
	int need_swap, target_index;
	int *ia = (int*)a;
	int *ib = (int*)b;
	DpPopulation*popunion = (DpPopulation*)user_data;
	DpIndivid*i1 = popunion->individ[(*ia)];
	DpIndivid*i2 = popunion->individ[(*ib)];
	target_index = popunion->target_index;
	need_swap = 0;
	if ( i1->targets[target_index] > i2->targets[target_index] ) {
		need_swap = 1;
	} else if ( i1->targets[target_index] < i2->targets[target_index] ) {
		need_swap = -1;
	}
	return need_swap;
}

void dp_population_cr2cost(DpPopulation*pop)
{
    GArray*curr_front;
    int i, curr_ind, j, k;
    int ntargets = pop->individ[0]->ntargets;
    double tmax, tmin, trange;
    DpIndivid* individ, *individ_l, *individ_r;
    for ( k = 0; k < pop->fronts->len; k++ ) {
        curr_front = g_array_index (pop->fronts, GArray*, k);
        for ( j = 0; j < curr_front->len; j++ ) {
            curr_ind = g_array_index (curr_front, int, j);
            individ = pop->individ[curr_ind];
            individ->crdist = 0;
        }
        for( i = 0; i < ntargets; i++) {
            pop->target_index = i;
            g_array_sort_with_data(curr_front, (GCompareDataFunc) dp_population_target_compare, (gpointer) pop);
            curr_ind = g_array_index (curr_front, int, curr_front->len - 1);
            individ = pop->individ[curr_ind];
            tmax = individ->targets[i];
            individ->crdist = G_MAXDOUBLE;
            curr_ind = g_array_index (curr_front, int, 0);
            individ = pop->individ[curr_ind];
            tmin = individ->targets[i];
            individ->crdist = G_MAXDOUBLE;
            trange = tmax - tmin;
            if ( trange < G_MINDOUBLE ) {
                continue;
            }
            curr_ind = g_array_index (curr_front, int, 1);
            individ_r = pop->individ[curr_ind];
            for ( j = 1; j < curr_front->len - 1; j++ ) {
                individ_l = individ;
                individ = individ_r;
                curr_ind = g_array_index (curr_front, int, j + 1);
                individ_r = pop->individ[curr_ind];
                individ->crdist += (individ_r->targets[i] - individ_l->targets[i])/trange;
            }
        }
    }
    for ( j = 0; j < pop->size; j++) {
		pop->individ[j]->cost = ( pop->individ[j]->crdist > G_MINDOUBLE ) ? 1 / pop->individ[j]->crdist : G_MAXDOUBLE;
	}
}

int dp_individ_ages_descending(void *p1, void *p2, void *user_data)
{
	DpPopulation*population = (DpPopulation*)user_data;
	int*i1 = (int*)p1;
	int*i2 = (int*)p2;
	DpIndivid*individ_1 = population->individ[(*i1)];
	DpIndivid*individ_2 = population->individ[(*i2)];
	if ( individ_1->age < individ_2->age ) {
		return 1;
	} else if ( individ_1->age > individ_2->age ) {
		return -1;
	}
	return 0;
}

int dp_individ_failures_descending(void *p1, void *p2, void *user_data)
{
	DpPopulation*population = (DpPopulation*)user_data;
	int*i1 = (int*)p1;
	int*i2 = (int*)p2;
	DpIndivid*individ_1 = population->individ[(*i1)];
	DpIndivid*individ_2 = population->individ[(*i2)];
	if ( individ_1->failures < individ_2->failures ) {
		return 1;
	} else if ( individ_1->failures > individ_2->failures ) {
		return -1;
	}
	return 0;
}

int dp_individ_cost_ascending(void *p1, void *p2, void *user_data)
{
	DpPopulation*population = (DpPopulation*)user_data;
	int*i1 = (int*)p1;
	int*i2 = (int*)p2;
	DpIndivid*individ_1 = population->individ[(*i1)];
	DpIndivid*individ_2 = population->individ[(*i2)];
	if ( individ_1->cost > individ_2->cost ) {
		return 1;
	} else if ( individ_1->cost < individ_2->cost ) {
		return -1;
	}
	return 0;
}

void dp_population_mpi_distribute(DpPopulation*pop, int mpi_id, int mpi_nnodes)
{
#ifdef MPIZE
	int dest;
	int source = 0;
	int ind_index, dest_index, ind_per_node, ind_per_last_node, number_of_comm_lines;
	MPI_Status *status;       /* status array of the MPI_Waitall function */
	MPI_Request *request;           /* handle array for receiving messages */
	ind_per_node = (int)ceil(pop->size / mpi_nnodes);
	ind_per_last_node = pop->size - ind_per_node * (mpi_nnodes - 1);
	if (mpi_id == 0) { /* master - sending*/
		number_of_comm_lines = pop->size - ind_per_node;
		request = (MPI_Request *)calloc(number_of_comm_lines, sizeof(MPI_Request));
		status =  (MPI_Status *)calloc(number_of_comm_lines, sizeof(MPI_Status));
		dest_index = 0;
		for (dest = 1; dest < mpi_nnodes - 1; dest++) { /* my_id = 0 doesn't send anything */
			dest_index += ind_per_node;
			for (ind_index = dest_index; ind_index < dest_index + ind_per_node; ind_index++) {
				MPI_Isend(pop->individ[ind_index]->x, pop->ind_size, MPI_DOUBLE, dest, ind_index, MPI_COMM_WORLD, &request[ind_index - ind_per_node]);
			}
		}
		dest_index += ind_per_node;
		dest = mpi_nnodes - 1;
		for (ind_index = dest_index; ind_index < dest_index + ind_per_last_node; ind_index++) {
			MPI_Isend(pop->individ[ind_index]->x, pop->ind_size, MPI_DOUBLE, dest, ind_index, MPI_COMM_WORLD, &request[ind_index - ind_per_node]);
		}
		MPI_Waitall(number_of_comm_lines, request, status);
		free(request);
		free(status);
	} else { /* slave - recieveing */
		number_of_comm_lines = (mpi_id == (mpi_nnodes - 1)) ? ind_per_last_node : ind_per_node;
		request = (MPI_Request *)calloc(number_of_comm_lines, sizeof(MPI_Request));
		status =  (MPI_Status *)calloc(number_of_comm_lines, sizeof(MPI_Status));
		dest_index = ind_per_node * mpi_id;
		for (ind_index = dest_index; ind_index < dest_index + number_of_comm_lines; ind_index++) {
			MPI_Irecv(pop->individ[ind_index]->x, pop->ind_size, MPI_DOUBLE, source, ind_index, MPI_COMM_WORLD, &request[ind_index - dest_index]);
		}
		MPI_Waitall(number_of_comm_lines, request, status);
		free(request);
		free(status);
	}
#endif
}

void dp_population_mpi_gather(DpPopulation*pop, int mpi_id, int mpi_nnodes)
{
#ifdef MPIZE
    int dest;
	int source = 0;
	int ind_index, dest_index, ind_per_node, ind_per_last_node, number_of_comm_lines, comm_index;
	MPI_Status *status;       /* status array of the MPI_Waitall function */
	MPI_Request *request;           /* handle array for receiving messages */
	ind_per_node = (int)ceil(pop->size / mpi_nnodes);
	ind_per_last_node = pop->size - ind_per_node * (mpi_nnodes - 1);
	if (mpi_id == 0) { /* master - gathering*/
		number_of_comm_lines = 3 * (pop->size - ind_per_node);
		request = (MPI_Request *)calloc(number_of_comm_lines, sizeof(MPI_Request));
		status =  (MPI_Status *)calloc(number_of_comm_lines, sizeof(MPI_Status));
		dest_index = 0;
		comm_index = 0;
		for (dest = 1; dest < mpi_nnodes - 1; dest++) { /* my_id = 0 doesn't send anything */
			dest_index += ind_per_node;
			for (ind_index = dest_index; ind_index < dest_index + ind_per_node; ind_index++) {
				MPI_Irecv(pop->individ[ind_index]->targets, pop->individ[ind_index]->ntargets, MPI_DOUBLE, dest, 3 * ind_index, MPI_COMM_WORLD, &request[comm_index]);
				comm_index++;
				MPI_Irecv(pop->individ[ind_index]->precond, pop->individ[ind_index]->nprecond, MPI_DOUBLE, dest, 3 * ind_index + 1, MPI_COMM_WORLD, &request[comm_index]);
				comm_index++;
				MPI_Irecv(&(pop->individ[ind_index]->cost), 1, MPI_DOUBLE, dest, 3 * ind_index + 2, MPI_COMM_WORLD, &request[comm_index]);
				comm_index++;
			}
		}
		dest_index += ind_per_node;
		dest = mpi_nnodes - 1;
		for (ind_index = dest_index; ind_index < dest_index + ind_per_last_node; ind_index++) {
			MPI_Irecv(pop->individ[ind_index]->targets, pop->individ[ind_index]->ntargets, MPI_DOUBLE, dest, 3 * ind_index, MPI_COMM_WORLD, &request[comm_index]);
			comm_index++;
			MPI_Irecv(pop->individ[ind_index]->precond, pop->individ[ind_index]->nprecond, MPI_DOUBLE, dest, 3 * ind_index + 1, MPI_COMM_WORLD, &request[comm_index]);
			comm_index++;
			MPI_Irecv(&(pop->individ[ind_index]->cost), 1, MPI_DOUBLE, dest, 3 * ind_index + 2, MPI_COMM_WORLD, &request[comm_index]);
			comm_index++;
		}
		MPI_Waitall(number_of_comm_lines, request, status);
		free(request);
		free(status);
	} else { /* slave - sending */
		number_of_comm_lines = (mpi_id == (mpi_nnodes - 1)) ? ind_per_last_node : ind_per_node;
		request = (MPI_Request *)calloc(3 * number_of_comm_lines, sizeof(MPI_Request));
		status =  (MPI_Status *)calloc(3 * number_of_comm_lines, sizeof(MPI_Status));
		dest_index = ind_per_node * mpi_id;
		for (ind_index = dest_index; ind_index < dest_index + number_of_comm_lines; ind_index++) {
			MPI_Isend(pop->individ[ind_index]->targets, pop->individ[ind_index]->ntargets, MPI_DOUBLE, dest, 3 * ind_index, MPI_COMM_WORLD, &request[3 * (ind_index - dest_index)]);
			MPI_Isend(pop->individ[ind_index]->precond, pop->individ[ind_index]->nprecond, MPI_DOUBLE, dest, 3 * ind_index + 1, MPI_COMM_WORLD, &request[3 * (ind_index - dest_index) + 1]);
			MPI_Isend(&(pop->individ[ind_index]->cost), 1, MPI_DOUBLE, dest, 3 * ind_index + 2, MPI_COMM_WORLD, &request[3 * (ind_index - dest_index) + 2]);
		}
		MPI_Waitall(number_of_comm_lines, request, status);
		free(request);
		free(status);
	}
#endif
}
