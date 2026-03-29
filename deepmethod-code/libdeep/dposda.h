/***************************************************************************
 *            dposda.h
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
 
// implemantation of optmal steepest descent algorithm

#ifndef _DP_OSDA_H
#define _DP_OSDA_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include "dpindivid.h"
#include "dpevaluation.h"

typedef struct DpOsdaInfo {
	double step_parameter;
	double step_decrement;
	double derivative_step;
	int number_of_trials;
	DpIndivid*trial;
	DpIndivid*current;
	DpEvaluationCtrl *hevalctrl;
	int current_index;
} DpOsdaInfo;

DpOsdaInfo *dp_osda_info_init(DpEvaluation*heval, DpTarget*htarget, int worldid, int seed, double gamma_init, double roundoff_error, DpEvaluationStrategy eval_strategy, int number_of_trials, double step_parameter, double step_decrement, double derivative_step);

DpOsdaInfo *dp_osda_info_new (int number_of_trials, double step_parameter, double step_decrement, double derivative_step);

void dp_osda_update_step(DpOsdaInfo*hosdainfo);

void dp_osda_accept_step(DpOsdaInfo*hosdainfo, double*value);

void dp_osda_step_func (gpointer data, gpointer user_data);

void dp_osda_step(DpOsdaInfo*hosdainfo);

void dp_osda_post(DpOsdaInfo*hosdainfo);

void dp_osda_post_evaluate(DpOsdaInfo*hosdainfo);

void dp_osda_info_save(FILE*fp, DpOsdaInfo*hosdainfo);

void dp_osda_info_load(FILE*fp, DpOsdaInfo*hosdainfo);

#ifdef __cplusplus
}
#endif

#endif /* _DP_OSDA_H */
