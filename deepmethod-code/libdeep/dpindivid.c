/***************************************************************************
 *            dpindivid.c
 *
 *  Fri Mar 23 15:25:40 2012
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
#include "dpindivid.h"

DpIndivid*dp_individ_new(int size, int targets_size, int precond_size)
{
	int i;
    DpIndivid*individ;
	individ = (DpIndivid*)malloc(sizeof(DpIndivid));
	individ->size = size;
	individ->cost = 0;
	individ->x = (double*)calloc(size, sizeof(double));
	individ->y = (double*)calloc(size, sizeof(double));
	individ->z = (double*)calloc(size, sizeof(double));
	individ->targets = (double*)calloc(targets_size, sizeof(double));
	individ->ntargets = targets_size;
	for ( i = 0; i < targets_size; i++ ) {
        individ->targets[i] = G_MAXDOUBLE;
	}
	individ->precond = (double*)calloc(precond_size, sizeof(double));
	individ->nprecond = precond_size;
	individ->age = 0;
	individ->r1 = -1;
	individ->r2 = -1;
	individ->r3 = -1;
	individ->r4 = -1;
	individ->moves = 0;
	individ->failures = 0;
	individ->grads = 0;
	individ->pareto_front = -1;
	individ->dom_count = 0;
	individ->user_data = NULL;
	individ->dominated = NULL;
	individ->crdist = 0;
	g_mutex_init( &(individ->m) );
	individ->status = 0;
    individ->cost_ind = 0;
    individ->useWriteVector = FALSE;
	return individ;
}

void dp_individ_delete(DpIndivid*individ)
{
	free(individ->x);
	free(individ->y);
	free(individ->z);
	free(individ->targets);
	free(individ->precond);
	g_mutex_clear (&(individ->m));
	free(individ);
}

void dp_individ_copy_values(DpIndivid*individ, DpIndivid*trial)
{
	int i;
	individ->cost = trial->cost;
	individ->age = trial->age;
	individ->r1 = trial->r1;
	individ->r2 = trial->r2;
	individ->r3 = trial->r3;
	individ->r4 = trial->r4;
	individ->moves = trial->moves;
	individ->failures = trial->failures;
	individ->grads = trial->grads;
	individ->pareto_front = trial->pareto_front;
	individ->dom_count = trial->dom_count;
	individ->crdist = trial->crdist;
	for ( i = 0; i < individ->size; i++ ) {
		individ->x[i] = trial->x[i];
		individ->y[i] = trial->y[i];
		individ->z[i] = trial->z[i];
	}
	for ( i = 0; i < individ->ntargets && i < trial->ntargets; i++ ) {
		individ->targets[i] = trial->targets[i];
	}
	for ( i = 0; i < individ->nprecond && i < trial->nprecond; i++ ) {
		individ->precond[i] = trial->precond[i];
	}
	if (individ->user_data == NULL) individ->user_data = trial->user_data;
	individ->cost_ind = trial->cost_ind;

    for ( i = 0; i < TRIANGULAR; i++ ) {
        individ->p[i] = trial->p[i];
        individ->F[i] = trial->F[i];
    }
}

void dp_individ_save(FILE*fp, DpIndivid*individ)
{
	int i;
	fprintf(fp, "%13.9f\n", individ->cost);
	fprintf(fp, "%d\n", individ->age);
	fprintf(fp, "%d\n", individ->r1);
	fprintf(fp, "%d\n", individ->r2);
	fprintf(fp, "%d\n", individ->r3);
	fprintf(fp, "%d\n", individ->r4);
	fprintf(fp, "%d\n", individ->moves);
	fprintf(fp, "%d\n", individ->failures);
	fprintf(fp, "%d\n", individ->grads);
	for ( i = 0; i < individ->size; i++ ) {
	  fprintf(fp, "%13.9f ", individ->x[i]);
	}
	fprintf(fp, "\n");
	for ( i = 0; i < individ->size; i++ ) {
	  fprintf(fp, "%13.9f ", individ->y[i]);
	}
	fprintf(fp, "\n");
	for ( i = 0; i < individ->size; i++ ) {
	  fprintf(fp, "%13.9f ", individ->z[i]);
	}
	fprintf(fp, "\n");
	for ( i = 0; i < individ->ntargets; i++ ) {
	  fprintf(fp, "%13.9f ", individ->targets[i]);
	}
	fprintf(fp, "\n");
	for ( i = 0; i < individ->nprecond; i++ ) {
	  fprintf(fp, "%13.9f ", individ->precond[i]);
	}
    for ( i = 0; i < TRIANGULAR; i++ ) {
      fprintf(fp, "%13.9f ", &(individ->p[i]));
    }
    for ( i = 0; i < TRIANGULAR; i++ ) {
      fprintf(fp, "%13.9f ", &(individ->F[i]));
    }
	fprintf(fp, "\n");
}

void dp_individ_load(FILE*fp, DpIndivid*individ)
{
	int i;
	fscanf(fp, "%lf", &(individ->cost));
	fscanf(fp, "%d", &(individ->age));
	fscanf(fp, "%d", &(individ->r1));
	fscanf(fp, "%d", &(individ->r2));
	fscanf(fp, "%d", &(individ->r3));
	fscanf(fp, "%d", &(individ->r4));
	fscanf(fp, "%d", &(individ->moves));
	fscanf(fp, "%d", &(individ->failures));
	fscanf(fp, "%d", &(individ->grads));
	for ( i = 0; i < individ->size; i++ ) {
	  fscanf(fp, "%lf", &(individ->x[i]));
	}
	for ( i = 0; i < individ->size; i++ ) {
	  fscanf(fp, "%lf", &(individ->y[i]));
	}
	for ( i = 0; i < individ->size; i++ ) {
	  fscanf(fp, "%lf", &(individ->z[i]));
	}
	for ( i = 0; i < individ->ntargets; i++ ) {
	  fscanf(fp, "%lf", &(individ->targets[i]));
	}
	for ( i = 0; i < individ->nprecond; i++ ) {
	  fscanf(fp, "%lf", &(individ->precond[i]));
	}
    for ( i = 0; i < TRIANGULAR; i++ ) {
      fscanf(fp, "%lf", &(individ->p[i]));
    }
    for ( i = 0; i < TRIANGULAR; i++ ) {
      fscanf(fp, "%lf", &(individ->F[i]));
    }
}

void dp_individ_pack(DpIndivid*individ, double**buffer2send, int*bufferDim)
{
	int i;
	int ind_size = individ->size;
	*bufferDim = ind_size + 1 + ind_size + 1;
	*buffer2send = (double*)calloc(*bufferDim, sizeof(double));
	for( i = 0; i < ind_size; i++ ) {
		(*buffer2send)[i] = individ->x[i];
		(*buffer2send)[i + ind_size + 1] = individ->y[i];
	}
	(*buffer2send)[ind_size] = individ->cost;
}

void dp_individ_unpack(DpIndivid*individ, double*buffer2recv, int bufferDim)
{
	int i;
	int ind_size = individ->size;
	for( i = 0; i < ind_size; i++ ) {
		individ->x[i] = buffer2recv[i];
		individ->y[i] = buffer2recv[i + ind_size + 1];
	}
	individ->cost = buffer2recv[ ind_size ];
	individ->age = 0;
}

void dp_individ_dump(DpIndivid*individ)
{
	int i;
	int ind_size = individ->size;
	fprintf(stdout, "\n%% Individ dump (%d):\n", ind_size);
	for( i = 0; i < ind_size; i++ ) {
		fprintf(stdout, "%%%d %e %e %e\n", i, individ->x[i], individ->y[i], individ->z[i]);
	}
}

