/***************************************************************************
 *            dpmutation.h
 *
 *  Fri Mar 23 15:34:26 2012
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
 
#ifndef _DP_MUTATION_H
#define _DP_MUTATION_H

#ifdef __cplusplus
extern "C"
{
#endif
//#include "dppopulation.h"
//#include "dpdeep.h"
#include <glib.h>

typedef enum DpMutationStrategy {
	ES_0,
	ES_1,
	ES_2
} DpMutationStrategy;

typedef struct DpMutationControl {
	int size;
	double tau_1;
	double tau_2;
	double alpha;
	double phi;
	double *f;
	double *p;
	DpMutationStrategy mutation_strategy;
} DpMutationControl;

void dp_individ_mutation(DpMutationControl*control, GRand*hrand, DpIndivid*individ, DpIndivid*input_1, int start_index, int end_index);

DpMutationControl*dp_mutation_control_init(DpEvaluationCtrl*hevalctrl, DpPopulation*pop, double alpha, double phi);

void dp_mutation_control_update(DpMutationControl*mc, DpPopulation*pop, int start_index, int end_index);

#ifdef __cplusplus
}
#endif

#endif /* _DP_MUTATION_H */
