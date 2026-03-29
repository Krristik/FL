/***************************************************************************
 *            dpindivid.h
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

// Implementation of DE individ

#ifndef _DP_INDIVID_H
#define _DP_INDIVID_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#define TRIANGULAR 3

typedef struct DpIndivid { // структура данных для индивида
	int size;
	double cost;
	int cost_ind;
	int age;
	int pareto_front;
	int dom_count;
	double crdist;
    GArray *dominated; // GArray хранит и данные и длину массива
    double *x; // параматры безусловные - не ограничены
    double *y; // только если градиент
    double *z; // параметры, которые лежат внутри границ(ограничений)
	double *targets;
	int ntargets;
	double *precond;
	int nprecond;
	int invalid;
    //coefs for triangular rec
    double p[TRIANGULAR];
    double F[TRIANGULAR];
    //
    int r1; // случайные индексы для дипа, при мутации
    int r2; //
    int r3; //
    int r4; //
	int r5; //
	int moves;
	int failures;
	int grads;
	gpointer user_data;
    GMutex m; // переменная для потокобезопасности
	int status;
    gboolean useWriteVector;

} DpIndivid;

DpIndivid*dp_individ_new(int size, int targets_size, int precond_size);

// Print individ info into console
void dp_individ_dump(DpIndivid*individ);

void dp_individ_unpack(DpIndivid*individ, double*buffer2recv, int bufferDim);

void dp_individ_pack(DpIndivid*individ, double**buffer2send, int*bufferDim);

void dp_individ_copy_values(DpIndivid*individ, DpIndivid*trial);

void dp_individ_delete(DpIndivid*individ);

void dp_individ_save(FILE*fp, DpIndivid*individ);

void dp_individ_load(FILE*fp, DpIndivid*individ);

#ifdef __cplusplus
}
#endif

#endif /* _DP_INDIVID_H */
