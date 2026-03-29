/***************************************************************************
 *            dpmutation.c
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
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include "dpindivid.h"
#include "dppopulation.h"
#include "dpevaluation.h"
#include "dpmutation.h"

void dp_individ_mutation(DpMutationControl*control, GRand*hrand, DpIndivid*individ, DpIndivid*input_1, int start_index, int end_index)
{
	int i;
	double x1, x2, y;
	switch (control->mutation_strategy) {
		case ES_0:
			for ( i = start_index; i < end_index; i++ ) {
				y = individ->y[i];
				individ->z[i] = y;
/* FIXME */
//				x1 = g_rand_gauss(hrand, 0, 1);
//				x2 = g_rand_gauss(hrand, 0, 1);
				y = y * exp( control->tau_2 * x1 + control->tau_1 * x2 );
				if ( y > control->f[i] )
					y = control->f[i];
//				individ->x[i] = input_1->x[i] + y * g_rand_gauss(hrand, 0, 1);
				individ->y[i] = y;
			}
		break;
		default:
			fprintf(stderr, "\ndp_individ_mutation: unknown type %d\n", control->mutation_strategy);
			exit(-1);
	}
}

DpMutationControl*dp_mutation_control_init(DpEvaluationCtrl*hevalctrl, DpPopulation*pop, double alpha, double phi)
{
	DpMutationControl*mc;
	int i, j;
	double chi, y;
	double factor = 1.0 / sqrt((double)pop->ind_size);
	mc = (DpMutationControl*)malloc(sizeof(DpMutationControl));
	mc->size = pop->ind_size;
	mc->f = (double*)calloc(mc->size, sizeof(double));
	mc->p = (double*)calloc(mc->size, sizeof(double));
	for ( i = 0; i < mc->size; i++ ) {
		y =  factor * (hevalctrl->eval->points[i]->upper - hevalctrl->eval->points[i]->lower);
		mc->f[i] = y;
		mc->p[i] = y;
		for ( j = 0; j < pop->size; j++) {
			pop->individ[j]->y[i] = y;
		}
	}
	chi = 1.0 / ( 2.0 * mc->size ) + 1.0 / ( 2.0 * sqrt( (double)mc->size ) );
	mc->phi = sqrt( ( 2.0 / chi ) * log( ( 1.0 / alpha ) * ( exp( phi * phi * chi / 2.0 ) - ( 1 - alpha ) ) ) );
	mc->alpha = alpha;
	mc->tau_1 = sqrt( (double)mc->size );
	mc->tau_1 = sqrt( 2.0 * mc->tau_1 );
	mc->tau_1 = mc->phi / mc->tau_1;
	mc->tau_2 = sqrt( 2.0 * (double)mc->size );
	mc->tau_2 = mc->phi / mc->tau_2;
	return mc;
}

void dp_mutation_control_update(DpMutationControl*mc, DpPopulation*pop, int start_index, int end_index)
{
	int i, j;
	double y;
	for ( i = 0; i < mc->size; i++ ) {
		for ( j = start_index; j < end_index; j++) {
			y = pop->individ[j]->y[i];
			y = pop->individ[j]->z[i] + ( 1 - mc->alpha ) * ( y - pop->individ[j]->z[i] );
/*			y = pop->individ[j]->z[i] + mc->alpha * ( y - pop->individ[j]->z[i] );*/
			pop->individ[j]->y[i] = y;
		}
	}
}
