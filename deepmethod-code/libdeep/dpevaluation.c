/***************************************************************************
 *            dpevaluation.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef MPIZE
#include <mpi.h>
#endif
#include "dpevaluation.h"
#include "dpindivid.h"
#include "dppopulation.h"

DpEvaluationCtrl*dp_evaluation_ctrl_new()
{
	DpEvaluationCtrl*ec;
	ec = (DpEvaluationCtrl*)malloc(sizeof(DpEvaluationCtrl));
	ec->kount = 0;
	ec->roundoff_error = 1e-8;
	ec->eval = NULL;
	ec->eval_target = NULL;
	ec->eval_max_threads = 0;
	ec->exclusive = TRUE;
	return ec;
}

DpEvaluationCtrl*dp_evaluation_ctrl_init(int worldid, GKeyFile*gkf, gchar*groupname)
{
	DpEvaluationCtrl*hevalctrl;
	hevalctrl = dp_evaluation_ctrl_new();
	GError *gerror = NULL;
	gchar*str,**strlist;
	int retval = 0;
	if ( ( str = g_key_file_get_string(gkf, groupname, "gamma_init", &gerror) ) != NULL ) {
		hevalctrl->gamma_init = g_strtod( str , NULL);
		g_free(str);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "roundoff_error", &gerror) ) != NULL ) {
		hevalctrl->roundoff_error = g_strtod( str , NULL);
		g_free(str);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "seed", &gerror) ) != NULL ) {
		hevalctrl->seed = g_strtod( str , NULL);
		g_free(str);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "max_threads", &gerror) ) != NULL ) {
		hevalctrl->eval_max_threads = g_strtod( str , NULL);
		g_free(str);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "transform", &gerror) ) != NULL ) {
		if ( !g_strcmp0(str, "tanh") ) {
			hevalctrl->eval_strategy = tanh_trans_flag;
		} else if ( !g_strcmp0(str, "sin") ) {
			hevalctrl->eval_strategy = sin_trans_flag;
		} else if ( !g_strcmp0(str, "alg") ) {
			hevalctrl->eval_strategy = alg_trans_flag;
		} else if ( !g_strcmp0(str, "rand") ) {
			hevalctrl->eval_strategy = rand_trans_flag;
		}
		g_free(str);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	hevalctrl->yoffset = worldid;
	hevalctrl->hrand = g_rand_new_with_seed ((guint32)(hevalctrl->seed + worldid));
	g_mutex_init( &(hevalctrl->m) );
	hevalctrl->kounter = 0;
	return hevalctrl;
}

DpEvaluationCtrl*dp_evaluation_init(DpEvaluation*heval, DpTarget*htarget, int worldid, GKeyFile*gkf, gchar*groupname)
{
	int i;
	DpEvaluationCtrl*hevalctrl;
	hevalctrl = dp_evaluation_ctrl_init(worldid, gkf, groupname);
	hevalctrl->eval = heval;
	hevalctrl->eval_target = htarget;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		if ( hevalctrl->eval->points[i]->limited ) {
			hevalctrl->eval->points[i]->upper -= hevalctrl->roundoff_error;
			hevalctrl->eval->points[i]->lower += hevalctrl->roundoff_error;
			hevalctrl->eval->points[i]->alpha = 0.5 * (hevalctrl->eval->points[i]->upper + hevalctrl->eval->points[i]->lower);
			hevalctrl->eval->points[i]->beta = 0.5 * (hevalctrl->eval->points[i]->upper - hevalctrl->eval->points[i]->lower);
			hevalctrl->eval->points[i]->gamma = hevalctrl->gamma_init;
		}
	}
	return hevalctrl;
}

void dp_evaluation_individ_evaluate_prime(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost)
{
	int max_value_flag;
	individ->invalid = 1;
	dp_evaluation_individ_prepare(hevalctrl, individ);
	dp_target_eval_update_user_data(hevalctrl->eval_target, individ->user_data, tabu->z, index, cost);
	max_value_flag = dp_target_eval_prime (hevalctrl->eval_target, individ->z, &(individ->invalid), individ->y, individ->user_data, index, cost);
	dp_evaluation_individ_transform_grad(hevalctrl, individ, index);
}

void dp_evaluation_individ_copy(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost)
{
	dp_evaluation_individ_prepare(hevalctrl, individ);
	dp_target_eval_update_user_data(hevalctrl->eval_target, individ->user_data, tabu->z, index, cost);
	dp_evaluation_copy_to_param(hevalctrl, individ->z);
}

void dp_evaluation_individ_evaluate(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost)
{
	int max_value_flag = 0, i;
	individ->invalid = 1;
	max_value_flag = dp_evaluation_individ_prepare(hevalctrl, individ);
	if (max_value_flag == 1) {
		individ->cost = G_MAXDOUBLE;
		for ( i = 0; i < individ->ntargets; i++ ) {
			individ->targets[i] = G_MAXDOUBLE;
		}
		return;
	}
	dp_target_eval_update_user_data(hevalctrl->eval_target, individ->user_data, tabu->z, index, cost);
	max_value_flag = dp_target_eval (hevalctrl->eval_target, individ->z, &(individ->invalid), &(individ->cost), individ->targets, individ->precond, individ->user_data, index, cost);
	g_mutex_lock( &(hevalctrl->m) );
	hevalctrl->kounter++;
	g_mutex_unlock( &(hevalctrl->m) );
}

void dp_evaluation_individ_evaluate_precond(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, DpIndivid*tabu, int index, double cost)
{
	int max_value_flag = 0, i;
	individ->invalid = 1;
	max_value_flag = dp_evaluation_individ_prepare(hevalctrl, individ);
	if (max_value_flag == 1) {
/*		individ->cost = G_MAXDOUBLE;
		for ( i = 0; i < individ->ntargets; i++ ) {
			individ->targets[i] = G_MAXDOUBLE;
		}*/
		return;
	}
	dp_target_eval_update_user_data(hevalctrl->eval_target, individ->user_data, tabu->z, index, cost);
	max_value_flag = dp_target_eval_precond (hevalctrl->eval_target, individ->z, &(individ->invalid), individ->precond, individ->user_data, index, cost);
}

/*
 * if trial is better return 1
 *
 */

int dp_evaluation_individ_compare(const void *p1, const void *p2, void *user_data)
{
	int i, select_trial;
	DpIndivid**i1 = (DpIndivid**)p1;
	DpIndivid**i2 = (DpIndivid**)p2;
	DpIndivid*individ = *i1;
	DpIndivid*trial = *i2;
	DpEvaluationCtrl*hevalctrl = (DpEvaluationCtrl*)user_data;
	int ignore_cost = hevalctrl->eval_target->ignore_cost;
	double use_crdist = hevalctrl->eval_target->use_crdist;
	GRand*hrand = hevalctrl->hrand;
	select_trial = 0;
/* individ is accepted if -1 is returned */
	if ( individ->invalid == 1 && trial->invalid == 0 ) {
		select_trial = 1;
	} else if ( ignore_cost == 0 && trial->cost > individ->cost ) {
		select_trial = -1;
	} else if ( use_crdist > 0 && ( trial->pareto_front < individ->pareto_front || (( trial->pareto_front == individ->pareto_front ) && trial->crdist > individ->crdist) ) && g_rand_double(hrand) < use_crdist ) {
		select_trial = 1;
	} else {
		select_trial = -1;
		for ( i = 0; i < hevalctrl->eval_target->array_size; i++ ) {
			if ( trial->targets[i] < individ->targets[i] && g_rand_double(hrand) < hevalctrl->eval_target->penalty[i]->rank ) {
				select_trial = 1;
				break;
			}
		}
	}
	return select_trial;
}

int dp_evaluation_cr_compare(gconstpointer a, gconstpointer b, gpointer user_data)
{
	int need_swap;
	int *ia = (int*)a;
	int *ib = (int*)b;
	DpPopulation*popunion = (DpPopulation*)user_data;
	DpIndivid*i1 = popunion->individ[(*ia)];
	DpIndivid*i2 = popunion->individ[(*ib)];
	need_swap = 0;
	if ( i1->crdist < i2->crdist ) {
		need_swap = 1;
	} else if ( i1->crdist > i2->crdist ) {
		need_swap = -1;
	}
	return need_swap;
}

/* return 1 if p1 nondominates p2 *
   that is p2 dominates p1 and we *
   need to swap them to get non-  *
   dominated sorting              *
***********************************/

int dp_evaluation_individ_dominates_compare(const void *p1, const void *p2, void *user_data)
{
	int i, dominates;
	DpIndivid**i1 = (DpIndivid**)p1;
	DpIndivid**i2 = (DpIndivid**)p2;
	DpIndivid*individ = *i1;
	DpIndivid*trial = *i2;
	DpEvaluationCtrl*hevalctrl = (DpEvaluationCtrl*)user_data;
	dominates = 1;
	if ( individ->invalid == 1 && trial->invalid == 0 ) {
		dominates = 0;
	} else {
		dominates = 1;
		for ( i = 0; i < hevalctrl->eval_target->array_size; i++ ) {
			if ( trial->targets[i] < individ->targets[i] ) {
				dominates = 0;
				break;
			}
		}
	}
	return dominates;
}

/*********************************************************************
 ** gaussian random: normal distribution                            **
 ** N(mean, dev)                                                    **
 **                                                                 **
 ** void ShareNormalRandVec(s, n, mean, dev)                        **
 ** return s=vec(n)                                                 **
 **                                                                 **
 ** source code from                                                **
 **        http://remus.rutgers.edu/~rhoads/Code/code.html          **
 ** taken from Numerical Recipes in C p.289                         **
 *********************************************************************/

double dp_rand_gauss(GRand*hrand, double mean, double dev)
{
	static double V2, fac;
	static int phase = 0;
	double S, Z, U1, U2, V1;
	if ( phase == 1 ) {
		Z = V2 * fac;
	} else {
		do {
			U1 = g_rand_double(hrand);
			U2 = g_rand_double(hrand);
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while ( S >= 1 || S == 0.0 );
		fac = sqrt (-2 * log(S) / S);
		Z = V1 * fac;
	}
	phase = 1 - phase;
	Z = mean + dev*Z;
	return Z;
}

int dp_evaluation_individ_prepare(DpEvaluationCtrl*hevalctrl, DpIndivid*individ)
{
	int i;
	double y;
	GRand*hrand = hevalctrl->hrand;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		if ( hevalctrl->eval->points[i]->limited ) {
			if ( hevalctrl->eval_strategy == sin_trans_flag ) {
				y = hevalctrl->eval->points[i]->alpha + hevalctrl->eval->points[i]->beta * sin(hevalctrl->eval->points[i]->gamma * individ->x[i]);
				y /= hevalctrl->eval->points[i]->scale;
				individ->z[i] = y;
			} else if ( hevalctrl->eval_strategy == tanh_trans_flag ) {
				y = hevalctrl->eval->points[i]->alpha + hevalctrl->eval->points[i]->beta * tanh(hevalctrl->eval->points[i]->gamma * individ->x[i]);
				y /= hevalctrl->eval->points[i]->scale;
				individ->z[i] = y;
			} else if ( hevalctrl->eval_strategy == alg_trans_flag ) {
				if ( individ->x[i] > hevalctrl->eval->points[i]->lower && individ->x[i] < hevalctrl->eval->points[i]->upper ) {
					individ->z[i] = individ->x[i] / hevalctrl->eval->points[i]->scale;
				} else {
					individ->x[i] = hevalctrl->eval->points[i]->lower + 2.0 * g_rand_double(hrand) * hevalctrl->eval->points[i]->beta;
					individ->z[i] = individ->x[i] / hevalctrl->eval->points[i]->scale;
				}
			} else if ( hevalctrl->eval_strategy == rand_trans_flag ) {
				while ( individ->x[i] < hevalctrl->eval->points[i]->lower || individ->x[i] > hevalctrl->eval->points[i]->upper ) {
					individ->x[i] = hevalctrl->eval->points[i]->alpha + dp_rand_gauss(hrand, 0, 1) * hevalctrl->eval->points[i]->beta;
				}
				individ->z[i] = individ->x[i] / hevalctrl->eval->points[i]->scale;
			}
		} else {
			individ->z[i] = individ->x[i] / hevalctrl->eval->points[i]->scale;
		}
		if (individ->z[i] != individ->z[i]){ // is NaN
			return 1;
		}
	}
	return 0;
}

void dp_evaluation_copy_to_grad(DpEvaluationCtrl*hevalctrl, double*x)
{
	int i;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		(*(hevalctrl->eval->points[i]->grad)) = x[i];
	}
}

void dp_evaluation_copy_to_param(DpEvaluationCtrl*hevalctrl, double*x)
{
	int i;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		(*(hevalctrl->eval->points[i]->param)) = x[i];
	}
}

void dp_evaluation_individ_scramble(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, double eps)
{
	int i;
	double y, z;
	GRand*hrand = hevalctrl->hrand;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		if ( eps > 0 ) {
			z = (*(hevalctrl->eval->points[i]->param)) * hevalctrl->eval->points[i]->scale;
			y = g_rand_double_range(hrand, (1.0 - eps) * z, (1.0 + eps) * z);
		} else {
			y = g_rand_double_range(hrand, hevalctrl->eval->points[i]->lower, hevalctrl->eval->points[i]->upper);
		}
		if ( hevalctrl->eval->points[i]->limited ) {
			if ( hevalctrl->eval_strategy == sin_trans_flag )
				individ->x[i] = asin ( ( y - hevalctrl->eval->points[i]->alpha ) / hevalctrl->eval->points[i]->beta) / hevalctrl->eval->points[i]->gamma;
			else if ( hevalctrl->eval_strategy == tanh_trans_flag )
				individ->x[i] = atanh ( ( y - hevalctrl->eval->points[i]->alpha ) / hevalctrl->eval->points[i]->beta) / hevalctrl->eval->points[i]->gamma;
			else if ( hevalctrl->eval_strategy == alg_trans_flag )
				individ->x[i] = y;
			else if ( hevalctrl->eval_strategy == rand_trans_flag )
				individ->x[i] = y;
		} else {
			individ->x[i] = y;
		}
	}
}


void dp_evaluation_individ_set(DpEvaluationCtrl*hevalctrl, DpIndivid*individ)
{
	int i;
	double y;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		y = (*(hevalctrl->eval->points[i]->param)) * hevalctrl->eval->points[i]->scale;
		if ( hevalctrl->eval->points[i]->limited ) {
			if ( hevalctrl->eval_strategy == sin_trans_flag )
				individ->x[i] = asin ( ( y - hevalctrl->eval->points[i]->alpha ) / hevalctrl->eval->points[i]->beta) / hevalctrl->eval->points[i]->gamma;
			else if ( hevalctrl->eval_strategy == tanh_trans_flag )
				individ->x[i] = atanh ( ( y - hevalctrl->eval->points[i]->alpha ) / hevalctrl->eval->points[i]->beta) / hevalctrl->eval->points[i]->gamma;
			else if ( hevalctrl->eval_strategy == alg_trans_flag )
				individ->x[i] = y;
			else if ( hevalctrl->eval_strategy == rand_trans_flag )
				individ->x[i] = y;
		} else {
			individ->x[i] = y;
		}
	}
}

void dp_evaluation_individ_set_grad(DpEvaluationCtrl*hevalctrl, DpIndivid*individ)
{
	int i;
	double z, w;
	for ( i = 0; i < hevalctrl->eval->size; i++) {
		z = (*(hevalctrl->eval->points[i]->grad)) * hevalctrl->eval->points[i]->scale;
		if ( hevalctrl->eval->points[i]->limited ) {
			if ( hevalctrl->eval_strategy == sin_trans_flag )
				individ->y[i] = z * hevalctrl->eval->points[i]->gamma * hevalctrl->eval->points[i]->beta * cos ( individ->x[i] * hevalctrl->eval->points[i]->gamma );
			else if ( hevalctrl->eval_strategy == tanh_trans_flag ) {
				w = tanh( individ->x[i] * hevalctrl->eval->points[i]->gamma );
				individ->y[i] = z * hevalctrl->eval->points[i]->gamma * hevalctrl->eval->points[i]->beta * ( 1 - w * w );
			} else if ( hevalctrl->eval_strategy == alg_trans_flag )
				individ->y[i] = z;
			else if ( hevalctrl->eval_strategy == rand_trans_flag )
				individ->y[i] = z;
		} else {
			individ->y[i] = z;
		}
	}
}

void dp_evaluation_individ_transform_grad(DpEvaluationCtrl*hevalctrl, DpIndivid*individ, int index)
{
	double orig_grad = individ->y[index];
	double z = orig_grad * hevalctrl->eval->points[index]->scale;
	double w;
	if ( hevalctrl->eval->points[index]->limited ) {
		if ( hevalctrl->eval_strategy == sin_trans_flag ) {
			individ->y[index] = z * hevalctrl->eval->points[index]->gamma * hevalctrl->eval->points[index]->beta * cos ( individ->x[index] * hevalctrl->eval->points[index]->gamma );
		} else if ( hevalctrl->eval_strategy == tanh_trans_flag ) {
			w = tanh( individ->x[index] * hevalctrl->eval->points[index]->gamma );
			individ->y[index] = z * hevalctrl->eval->points[index]->gamma * hevalctrl->eval->points[index]->beta * ( 1 - w * w );
		} else if ( hevalctrl->eval_strategy == alg_trans_flag ) {
			individ->y[index] = z;
		} else if ( hevalctrl->eval_strategy == rand_trans_flag ) {
			individ->y[index] = z;
		}
	} else {
		individ->y[index] = z;
	}
}

DpPopulation*dp_evaluation_population_init_serial(DpEvaluationCtrl*hevalctrl, int size, int max_size, double noglobal_eps)
{
	DpPopulation*pop;
	int i;
	pop = dp_population_new(max_size, size, hevalctrl->eval->size, hevalctrl->eval_target->array_size, hevalctrl->eval_target->precond_size);
	dp_evaluation_individ_set(hevalctrl, pop->individ[0]);
	pop->individ[0]->user_data = dp_target_eval_get_user_data(hevalctrl->eval_target);
	dp_evaluation_individ_evaluate(hevalctrl, pop->individ[0], pop->individ[0], 0, 0);
	for ( i = 1; i < size; i++) {
		dp_evaluation_individ_scramble(hevalctrl, pop->individ[i], noglobal_eps);
		pop->individ[i]->user_data = dp_target_eval_get_user_data(hevalctrl->eval_target);
		dp_evaluation_individ_evaluate(hevalctrl, pop->individ[i], pop->individ[i], i, 0);
	}
	dp_population_update(pop, 0, pop->size);
	return pop;
}

void dp_evaluation_population_init_func (gpointer data, gpointer user_data)
{
	DpEvaluationCtrl*hevalctrl = (DpEvaluationCtrl*)user_data;
	DpIndivid*individ = (DpIndivid*)data;
	g_mutex_lock( &(individ->m) );
	individ->status = 1;
	g_mutex_unlock( &(individ->m) );
	dp_evaluation_individ_evaluate(hevalctrl, individ, individ, 0, 0);
	g_mutex_lock( &(individ->m) );
	individ->status = 0;
	g_mutex_unlock( &(individ->m) );
}

DpPopulation*dp_evaluation_population_init(DpEvaluationCtrl*hevalctrl, int size, int max_size, double noglobal_eps)
{
	DpPopulation*pop;
	int i, istart = 0;
	gboolean immediate_stop = TRUE;
	gboolean wait_finish = FALSE;
	GError *gerror = NULL;
	GMainContext *gcontext = g_main_context_default();
	gulong microseconds = G_USEC_PER_SEC / 1000;
	int notdone, currstatus;
	pop = dp_population_new(max_size, size, hevalctrl->eval->size, hevalctrl->eval_target->array_size, hevalctrl->eval_target->precond_size);
	if ( noglobal_eps == 0 ) {
		dp_evaluation_individ_set(hevalctrl, pop->individ[0]);
		pop->individ[0]->user_data = dp_target_eval_get_user_data(hevalctrl->eval_target);
		istart = 1;
		pop->individ[0]->cost = G_MAXDOUBLE;
	}
	for ( i = istart; i < max_size; i++) {
		dp_evaluation_individ_scramble(hevalctrl, pop->individ[i], noglobal_eps);
		pop->individ[i]->user_data = dp_target_eval_get_user_data(hevalctrl->eval_target);
		pop->individ[i]->cost = G_MAXDOUBLE;
	}
#ifdef MPIZE
/* MPI initialization steps */
	int world_id = 0, world_count = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &world_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_id);
	int ind_per_node = (int)ceil(pop->size / world_count);
	int ind_per_last_node = pop->size - ind_per_node * (world_count - 1);
	dp_population_mpi_distribute(pop, world_id, world_count);
#endif
	if ( hevalctrl->eval_max_threads > 0 ) {
		hevalctrl->gthreadpool = g_thread_pool_new ((GFunc) dp_evaluation_population_init_func, (gpointer) hevalctrl, hevalctrl->eval_max_threads, hevalctrl->exclusive, &gerror);
		if ( gerror != NULL ) {
			g_error("%s", gerror->message);
		}
		for ( i = pop->slice_a; i < pop->slice_b; i++) {
			g_thread_pool_push (hevalctrl->gthreadpool, (gpointer)(pop->individ[i]), &gerror);
			if ( gerror != NULL ) {
				g_error("%s", gerror->message);
			}
		}
		notdone = 1;
		while(notdone == 1 || g_thread_pool_unprocessed (hevalctrl->gthreadpool) > 0) {
			g_main_context_iteration(gcontext, FALSE);
			g_usleep (microseconds);
			notdone = 0;
			for ( i = pop->slice_a; i < pop->slice_b; i++ ) {
				g_mutex_lock( &(pop->individ[i]->m) );
				currstatus = pop->individ[i]->status;
				g_mutex_unlock( &(pop->individ[i]->m) );
/*                                currstatus = dp_deep_evaluate_status (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);*/
				if ( currstatus == 1 ) {
					notdone = 1;
				}
			}
		}
		g_thread_pool_free (hevalctrl->gthreadpool, immediate_stop, wait_finish);
	} else {
		for ( i = pop->slice_a; i < pop->slice_b; i++) {
			dp_evaluation_population_init_func ((gpointer)(pop->individ[i]), (gpointer) hevalctrl);
		}
	}
#ifdef MPIZE
	dp_population_mpi_gather(pop, world_id, world_count);
#endif
	dp_population_update(pop, 0, pop->cur_size);
	return pop;
}

DpIndivid*dp_evaluation_individ_init(DpEvaluationCtrl*hevalctrl)
{
	DpIndivid*individ;
	individ = dp_individ_new(hevalctrl->eval->size, hevalctrl->eval_target->array_size, hevalctrl->eval_target->precond_size);
	dp_evaluation_individ_set(hevalctrl, individ);
	individ->user_data = dp_target_eval_get_user_data(hevalctrl->eval_target);
	dp_evaluation_individ_evaluate(hevalctrl, individ, individ, 0, 0);
	return individ;
}

