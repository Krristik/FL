/***************************************************************************
 *            dposda.c
 *
 *  Tue Aug 21 17:37:57 2012
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
#include "dposda.h"

DpOsdaInfo *dp_osda_info_init(DpEvaluation*heval, DpTarget*htarget, int worldid, int seed, double gamma_init, double roundoff_error, DpEvaluationStrategy eval_strategy, int number_of_trials, double step_parameter, double step_decrement, double derivative_step)
{
	DpOsdaInfo*hosdainfo = dp_osda_info_new(number_of_trials, step_parameter, step_decrement, derivative_step);
//	hosdainfo->hevalctrl = dp_evaluation_init(heval, htarget, worldid, seed, gamma_init, roundoff_error, 0, eval_strategy);
	hosdainfo->trial = dp_evaluation_individ_init(hosdainfo->hevalctrl);
	hosdainfo->current = dp_evaluation_individ_init(hosdainfo->hevalctrl);
	return hosdainfo;
}

void dp_osda_info_save(FILE*fp, DpOsdaInfo*hosdainfo)
{

}

void dp_osda_info_load(FILE*fp, DpOsdaInfo*hosdainfo)
{

}

DpOsdaInfo *dp_osda_info_new (int number_of_trials, double step_parameter, double step_decrement, double derivative_step)
{
	DpOsdaInfo*hosdainfo;
	hosdainfo = (DpOsdaInfo*)malloc(sizeof(DpOsdaInfo));
	hosdainfo->number_of_trials = number_of_trials;
	hosdainfo->step_parameter = step_parameter;
	hosdainfo->step_decrement = step_decrement;
	hosdainfo->derivative_step = derivative_step;
	hosdainfo->current_index = 0;
	return hosdainfo;
}

void dp_osda_update_step(DpOsdaInfo*hosdainfo)
{

}

void dp_osda_accept_step(DpOsdaInfo*hosdainfo, double*value)
{
	*value = hosdainfo->current->cost;
}

void dp_osda_step_func (gpointer data, gpointer user_data)
{
	DpOsdaInfo*hosdainfo = (DpOsdaInfo*)user_data;
	DpIndivid*my_trial = hosdainfo->trial;
	DpIndivid*my_individ = hosdainfo->current;
	double *pdelta = (double*)data;
	int iter, index = hosdainfo->current_index;
	double delta, cost = my_individ->cost;
	dp_evaluation_individ_evaluate_prime(hosdainfo->hevalctrl, my_individ, my_individ, index, cost);
	iter = 0;
	delta = hosdainfo->step_parameter * my_individ->y[index];
	my_individ->age++;
	do {
		iter++;
		dp_individ_copy_values(my_trial, my_individ);
		my_trial->age = 0;
		my_trial->x[hosdainfo->current_index] = my_individ->x[hosdainfo->current_index] - delta;
		dp_evaluation_individ_evaluate(hosdainfo->hevalctrl, my_trial, my_individ, index, cost);
		delta *= hosdainfo->step_decrement;
	} while ( my_trial->cost >= my_individ->cost && iter < hosdainfo->number_of_trials );
	if ( my_trial->cost >= my_individ->cost ) {
		dp_individ_copy_values(my_trial, my_individ);
	}
	(*pdelta) = delta;
}

void dp_osda_step(DpOsdaInfo*hosdainfo)
{
	DpIndivid*trial = hosdainfo->trial;
	DpIndivid*individ = hosdainfo->current;
	int index = hosdainfo->current_index;
	double delta;
	dp_osda_step_func((gpointer) &delta, (gpointer) hosdainfo);
	individ->y[index] = delta;
	hosdainfo->current_index++;
	if ( hosdainfo->current_index >= hosdainfo->current->size ) {
		hosdainfo->current_index = 0;
	}
	hosdainfo->current = trial;
	hosdainfo->trial = individ;
}

void dp_osda_post(DpOsdaInfo*hosdainfo)
{
	DpIndivid*current = hosdainfo->current;
	dp_evaluation_individ_copy(hosdainfo->hevalctrl, current, current, hosdainfo->current_index, current->cost);
}

void dp_osda_post_evaluate(DpOsdaInfo*hosdainfo)
{
	DpIndivid*current = hosdainfo->current;
	dp_evaluation_individ_evaluate_precond(hosdainfo->hevalctrl, current, current, hosdainfo->current_index, current->cost);
}

