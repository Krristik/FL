/***************************************************************************
 *            dpdeep.h
 *
 *  Fri Mar 23 15:45:54 2012
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

// Control deep during execution

#ifndef _DP_DEEP_H
#define _DP_DEEP_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include "dpindivid.h"
#include "dppopulation.h"
#include "dprecombination.h"
#include "dpevaluation.h"
// #include "dparchive.h"

typedef enum DpSelector { // выбор состояния
	DpSelectorGenerated = (1 << 0),
	DpSelectorEvaluated  = (1 << 1),
	DpSelectorSelected  = (1 << 2),
	DpSelectorNone = (1 << 3)
} DpSelector;

typedef struct DpDeepInfo { // информация из .ini файла
	int debug;
	int population_size;
	double fmin; // значение функции на оптимуме предыдущего шага
	int population_max_size;
    int lower_bound;
	int lb; //счетчики
	int b;  //
	int st; //
	int w;  //
    int R;  //допустимое количество операций
    int s;  //часть популяции, которую мы отбираем для увеличения/уменьшения популяции
	int es_lambda;
	int es_cutoff;
	int es_kind;
	int ca_flag;
	double p_param; //параметр p -> передать в population
	double e_param; //параметр e -> передать в population
	DpSelector selector;
	double noglobal_eps;
	double substeps;
	double substieps;
	int mean_cost;
	int last_method;
	DpPopulation*trial;
	DpPopulation*population;
	DpPopulation*popunion;
	DpRecombinationControl *recombination_control;
	DpEvaluationCtrl *hevalctrl;
    GThreadPool*gthreadpool;
	gint max_threads;
	gboolean exclusive;
    DpArchive* archive;
} DpDeepInfo;

DpDeepInfo *dp_deep_info_init(DpEvaluation*heval, DpTarget*htarget, int worldid, GKeyFile*gkf, gchar*groupname);
// инициализация переменных дипа

DpDeepInfo *dp_deep_info_new(GKeyFile*gkf, gchar*groupname);
// входит в dp_deep_info_init

void dp_deep_update_step(DpDeepInfo*hdeepinfo);

void dp_deep_accept_step(DpDeepInfo*hdeepinfo, double*value);

void dp_deep_step_func (gpointer data, gpointer user_data);

void dp_deep_step(DpDeepInfo*hdeepinfo);

void dp_deep_post(DpDeepInfo*hdeepinfo);    // Selection of the best individ after DEEP step

void dp_deep_post_evaluate(DpDeepInfo*hdeepinfo);

void dp_deep_info_save(FILE*fp, DpDeepInfo *hdeepinfo);

void dp_deep_generate_step(DpDeepInfo*hdeepinfo);

void dp_deep_generate_dd_step(DpDeepInfo*hdeepinfo);

void dp_deep_generate_ca_step(DpDeepInfo*hdeepinfo);

void dp_deep_generate_ac_step(DpDeepInfo*hdeepinfo);

void dp_deep_evaluate_step(DpDeepInfo*hdeepinfo);

void dp_deep_select_step(DpDeepInfo*hdeepinfo);

void dp_de_select_step(DpDeepInfo*hdeepinfo);

void dp_deep_generate_func (gpointer data, gpointer user_data);

void dp_deep_generate_dd_func (gpointer data, gpointer user_data);

void dp_deep_generate_ca_func (gpointer data, gpointer user_data);

void dp_deep_generate_ac_func (gpointer data, gpointer user_data);

void dp_deep_evaluate_func (gpointer data, gpointer user_data);

void dp_deep_select_func (gpointer data, gpointer user_data);

void dp_de_select_func (gpointer data, gpointer user_data);

void dp_deep_info_load(FILE*fp, DpDeepInfo *hdeepinfo);

void change_pointers(DpPopulation*population, const int reduce_number);

void dp_deep_pop_size_change(DpDeepInfo*hdeepinfo);

void dp_deep_redstrat(DpDeepInfo*hdeepinfo);

void dp_deep_augstrat_1(DpDeepInfo*hdeepinfo);

void dp_deep_augstrat_2(DpDeepInfo*hdeepinfo);

int* generate_rand_ind(DpDeepInfo* hdeepinfo, int min);

void dp_init_xpbest_param(DpDeepInfo*hdeepinfo);

void dp_init_xamean_param(DpDeepInfo*hdeepinfo);

#ifdef __cplusplus
}
#endif

#endif /* _DP_DEEP_H */
