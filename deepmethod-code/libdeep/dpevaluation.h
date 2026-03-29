/***************************************************************************
 *            dpevaluation.h
 *
 *  Fri Mar 23 15:42:14 2012
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

// Target func evaluation #1

#ifndef _DP_EVALUATION_H
#define _DP_EVALUATION_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include "dptarget.h"
#include "dpindivid.h"
#include "dppopulation.h"

typedef enum DpEvaluationStrategy {
	tanh_trans_flag,
	sin_trans_flag,
	alg_trans_flag,
	rand_trans_flag
} DpEvaluationStrategy;

typedef struct DpEvaluationPoint { // точка эволюции
    double *param;  // Target func parameters
	int index;
	double upper;
	double lower;
	int limited;
    int dedupl;
    double *grad;   // Params for gradient method func
	double alpha;
	double beta;
	double gamma;
	double scale;
} DpEvaluationPoint;

typedef struct DpEvaluation { // структура для эволюции дипа
	int size;
	DpEvaluationPoint**points;
} DpEvaluation;

typedef struct DpEvaluationCtrl { //
	DpEvaluationStrategy eval_strategy;
    GThreadPool*gthreadpool; // пул потоков
	gboolean exclusive;
	int eval_max_threads;
	int kount;
	DpEvaluation*eval;
	int seed;
	int yoffset;
	double gamma_init;
	double roundoff_error;
	DpTarget*eval_target;
	GRand*hrand;
	GMutex m;
	int kounter;
} DpEvaluationCtrl;

void dp_evaluation_individ_evaluate_prime(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost);

void dp_evaluation_individ_copy(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost);

int dp_evaluation_individ_prepare(DpEvaluationCtrl*hevalctrl, DpIndivid*individ);

int dp_evaluation_individ_compare(const void *p1, const void *p2, void *user_data);

int dp_evaluation_individ_dominates_compare(const void *p1, const void *p2, void *user_data);

void dp_evaluation_individ_evaluate(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost);

void dp_evaluation_individ_evaluate_precond(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost);

DpEvaluationCtrl*dp_evaluation_init(DpEvaluation*heval, DpTarget*htarget, int worldid, GKeyFile*gkf, gchar*groupname);

DpEvaluationCtrl*dp_evaluation_ctrl_new();

void dp_evaluation_copy_to_grad(DpEvaluationCtrl*hevalctrl, double*x);

void dp_evaluation_copy_to_param(DpEvaluationCtrl*hevalctrl, double*x);

void dp_evaluation_individ_scramble(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, double eps);

void dp_evaluation_individ_set(DpEvaluationCtrl*hevalctrl, DpIndivid*individ);

void dp_evaluation_individ_set_grad(DpEvaluationCtrl*hevalctrl, DpIndivid*individ);

void dp_evaluation_individ_transform_grad(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, int index);

DpPopulation*dp_evaluation_population_init(DpEvaluationCtrl*hevalctrl, int size, int max_size, double noglobal_eps);

DpEvaluationCtrl*dp_evaluation_ctrl_init(int worldid, GKeyFile*gkf, gchar*groupname);

DpIndivid*dp_evaluation_individ_init(DpEvaluationCtrl*hevalctrl);

void dp_evaluation_population_init_func (gpointer data, gpointer user_data);

DpPopulation*dp_evaluation_population_init_serial(DpEvaluationCtrl*hevalctrl, int size, int max_size, double noglobal_eps);

int dp_evaluation_cr_compare(gconstpointer a, gconstpointer b, gpointer user_data);

#ifdef __cplusplus
}
#endif

#endif /* _DP_EVALUATION_H */
