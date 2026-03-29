/***************************************************************************
 *            dpdeep.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dpdeep.h"


DpDeepInfo *dp_deep_info_new(GKeyFile*gkf, gchar*groupname)
{
    DpDeepInfo*hdeepinfo;
    hdeepinfo = (DpDeepInfo*)malloc(sizeof(DpDeepInfo));
    GError *gerror = NULL;
    gchar*str,**strlist;
    int retval = 0;
    hdeepinfo->fmin = 0;
    hdeepinfo->population_size = 5;
    hdeepinfo->population_max_size = 10;
    hdeepinfo->lower_bound = 5;
    hdeepinfo->lb = 0;
    hdeepinfo->w = 0;
    hdeepinfo->b = 0;
    hdeepinfo->st = 0;
    hdeepinfo->R = 4;
    hdeepinfo->s = 0.05;
    hdeepinfo->es_lambda = 2;
    hdeepinfo->es_cutoff = 2;
    hdeepinfo->es_kind = 0;
    hdeepinfo->noglobal_eps = 0;
    hdeepinfo->substeps = 0;
    hdeepinfo->substieps = 0;
    hdeepinfo->mean_cost = 0;
    hdeepinfo->ca_flag = 0;
    hdeepinfo->p_param = 1;
    hdeepinfo->e_param = 1;
    if ( ( str = g_key_file_get_string(gkf, groupname, "population_size", &gerror) ) != NULL ) {
        hdeepinfo->population_size = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "max_size", &gerror) ) != NULL ) {
        hdeepinfo->population_max_size = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if (( str = g_key_file_get_string(gkf, groupname, "lower_bound", &gerror)) != NULL){
        hdeepinfo->lower_bound = g_strtod(str, NULL);
        g_free(str);
    } else{
        g_debug("%s", gerror->message);
        g_clear_error(&gerror);
    }
    if (( str = g_key_file_get_string(gkf, groupname, "lb_parameter", &gerror)) != NULL){
        hdeepinfo->lb = g_strtod(str, NULL);
        g_free(str);
    } else{
        g_debug("%s", gerror->message);
        g_clear_error(&gerror);
    }
    if (( str = g_key_file_get_string(gkf, groupname, "B_parameter", &gerror)) != NULL){
        hdeepinfo->b = g_strtod(str, NULL);
        g_free(str);
    } else{
        g_debug("%s", gerror->message);
        g_clear_error(&gerror);
    }
    if (( str = g_key_file_get_string(gkf, groupname, "W_parameter", &gerror)) != NULL){
        hdeepinfo->w = g_strtod(str, NULL);
        g_free(str);
    } else{
        g_debug("%s", gerror->message);
        g_clear_error(&gerror);
    }
    if (( str = g_key_file_get_string(gkf, groupname, "St_parameter", &gerror)) != NULL){
        hdeepinfo->st = g_strtod(str, NULL);
        g_free(str);
    } else{
        g_debug("%s", gerror->message);
        g_clear_error(&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "operation_number", &gerror) ) != NULL ) {
        hdeepinfo->R = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "selection", &gerror) ) != NULL ) {
        hdeepinfo->s = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "es_lambda", &gerror) ) != NULL ) {
        hdeepinfo->es_lambda = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "es_cutoff", &gerror) ) != NULL ) {
        hdeepinfo->es_cutoff = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "es_kind", &gerror) ) != NULL ) {
        hdeepinfo->es_kind = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "ca_flag", &gerror) ) != NULL ) {
        hdeepinfo->ca_flag = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "noglobal_eps", &gerror) ) != NULL ) {
        hdeepinfo->noglobal_eps = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "substeps", &gerror) ) != NULL ) {
        hdeepinfo->substeps = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "substieps", &gerror) ) != NULL ) {
        hdeepinfo->substieps = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "mean_cost", &gerror) ) != NULL ) {
        hdeepinfo->mean_cost = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "max_threads", &gerror) ) != NULL ) {
        hdeepinfo->max_threads = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    //p параметр для xpbest стратегии
    if ( ( str = g_key_file_get_string(gkf, groupname, "p_param", &gerror) ) != NULL ) {
        hdeepinfo->p_param = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    //e параметр для вычисления веса в xamean стратегии
    if ( ( str = g_key_file_get_string(gkf, groupname, "e_param", &gerror) ) != NULL ) {
        hdeepinfo->e_param = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    hdeepinfo->exclusive = TRUE;
    hdeepinfo->selector = DpSelectorNone;
    hdeepinfo->debug = 0;
    hdeepinfo->gthreadpool = NULL;
    hdeepinfo->last_method = -1;
    return hdeepinfo;
}

DpDeepInfo *dp_deep_info_init(DpEvaluation*heval, DpTarget*htarget, int worldid, GKeyFile*gkf, gchar*groupname)
{
    GError *gerror = NULL;
    int i;
    int pbest_idx;
    DpDeepInfo*hdeepinfo = dp_deep_info_new(gkf, groupname);
    DpRecombinationStrategy strategy;
    hdeepinfo->hevalctrl = dp_evaluation_init(heval, htarget, worldid, gkf, groupname);
    hdeepinfo->trial = dp_population_new(hdeepinfo->population_max_size, hdeepinfo->population_size, hdeepinfo->hevalctrl->eval->size, hdeepinfo->hevalctrl->eval_target->array_size, hdeepinfo->hevalctrl->eval_target->precond_size);
    hdeepinfo->population = dp_evaluation_population_init(hdeepinfo->hevalctrl, hdeepinfo->population_size, hdeepinfo->population_max_size, hdeepinfo->noglobal_eps);
    hdeepinfo->fmin = hdeepinfo->population->dmin;
    hdeepinfo->recombination_control = dp_recombination_control_init(hdeepinfo->hevalctrl->hrand, hdeepinfo->population, gkf, groupname);
    hdeepinfo->archive = NULL;
    DpArchiveMode mode;
    if (hdeepinfo->recombination_control->strategy == DE_3_pBest ||
        hdeepinfo->recombination_control->strategy == DE_3_Amean ||
        hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
            mode = INDIVIDS;
            hdeepinfo->archive = dp_archive_init(hdeepinfo->population, mode);
        } else {
            mode = DIFF;
        }
	if (hdeepinfo->recombination_control->use_archive_prob > 0) {
	    hdeepinfo->archive = dp_archive_init(hdeepinfo->population, mode); // initializing archvie 2NP size
	}
    hdeepinfo->popunion = dp_population_union(hdeepinfo->population, hdeepinfo->trial);
    //передача параметров p, e, m
    if (hdeepinfo->recombination_control->strategy == DE_3_pBest || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        dp_init_xpbest_param(hdeepinfo);
    }

    if (hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        dp_init_xamean_param(hdeepinfo);
    }

    // Init ind for trifngular
    ind_triand_init(hdeepinfo->population, hdeepinfo->recombination_control, hdeepinfo->hevalctrl->hrand);

    if ( hdeepinfo->max_threads > 0 ) {
        hdeepinfo->gthreadpool = g_thread_pool_new ((GFunc) dp_deep_evaluate_func, (gpointer) hdeepinfo, hdeepinfo->max_threads, hdeepinfo->exclusive, &gerror);
        g_thread_pool_set_max_unused_threads (-1);
        if ( gerror != NULL ) {
            g_error("%s", gerror->message);
        }
    }
    for (i = 0; i < hdeepinfo->population->size; i++) {
        hdeepinfo->trial->individ[i]->user_data = hdeepinfo->population->individ[i]->user_data;
    }
    return hdeepinfo;
}

void dp_deep_info_save(FILE*fp, DpDeepInfo *hdeepinfo)
{
    dp_population_save(fp, hdeepinfo->trial);
    dp_population_save(fp, hdeepinfo->population);
    dp_recombination_control_save(fp, hdeepinfo->recombination_control);
}

void dp_deep_info_load(FILE*fp, DpDeepInfo *hdeepinfo)
{
    dp_population_load(fp, hdeepinfo->trial);
    dp_population_load(fp, hdeepinfo->population);
    hdeepinfo->popunion = dp_population_union(hdeepinfo->population, hdeepinfo->trial);
    dp_recombination_control_load(fp, hdeepinfo->recombination_control);
}

void dp_deep_step_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4, r5;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    recombination_control->current_individ = my_id;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    r1 = ( ignore_cost == 0 ) ? population->imin : -1;
    do {
        r2 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r2 == my_id || r2 == r1 );
    do {
        r3 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r3 == my_id || r3 == r1 || r3 == r2 );
    do {
        r4 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r4 == my_id || r4 == r1 || r4 == r2 || r4 == r3 );
    do {
        r5 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r5 == my_id || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4);
    r1 = population->imin;
    my_tabu = population->individ[r1];
    start_index = g_rand_int_range (hrand, 0, population->ind_size);
    end_index = population->ind_size;
    dp_individ_copy_values(my_trial, my_individ);
    my_trial->age = 0;
    my_trial->r1 = r1;
    my_trial->r2 = r2;
    my_trial->r3 = r3;
    my_trial->r4 = r4;
    my_trial->r5 = r5;
    my_trial->moves = 0;
    my_trial->failures = 0;
    my_trial->grads = 0;

    //dp_individ_recombination(recombination_control, hrand, my_trial, population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], start_index, end_index);
    DifferenceVector* vectorRead;
    DifferenceVector* vectorWrite;

	if (recombination_control->use_archive_prob > 0 &&
        hdeepinfo->recombination_control->strategy != DE_3_pBest &&
        hdeepinfo->recombination_control->strategy != DE_3_Amean &&
        hdeepinfo->recombination_control->strategy != DE_3_dual_adaptive) {
		vectorRead = hdeepinfo->archive->difference_vectors[my_id];
    	vectorWrite = hdeepinfo->archive->difference_vectors[hdeepinfo->population_size + my_id];
	}

    if (hdeepinfo->recombination_control->strategy == DE_3_pBest) {
        dp_individ_recombination(recombination_control, hrand, my_trial,
                                population->individ[r1], population->individ[r2], population->individ[r3], population->individ[population->xpbest_idx], population->individ[r5],
                                start_index, end_index, vectorWrite, vectorRead);
    } else if (hdeepinfo->recombination_control->strategy == DE_3_Amean) {
        vectorRead = hdeepinfo->archive->xamean;
        dp_individ_recombination(recombination_control, hrand, my_trial,
                                population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], population->individ[r5],
                                start_index, end_index, vectorWrite, vectorRead);
    } else if (hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        vectorRead = hdeepinfo->archive->xamean;
        dp_individ_recombination(recombination_control, hrand, my_trial,
            population->individ[r1], population->individ[r2], population->individ[r3], population->individ[population->xpbest_idx], population->individ[r5],
            start_index, end_index, vectorWrite, vectorRead);
    } else {
        dp_individ_recombination(recombination_control, hrand, my_trial,
            population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], population->individ[r5],
            start_index, end_index, vectorWrite, vectorRead);
    }

    dp_evaluation_individ_evaluate(hdeepinfo->hevalctrl, my_trial, my_tabu, my_id, my_tabu->cost);
    recombination_control->FEs++;
    if ( ignore_cost == 0 && my_id == population->imin ) {
        if ( my_trial->cost >= my_individ->cost ) {
            //добавлять в архив здесь
            if (hdeepinfo->recombination_control->strategy == DE_3_pBest || hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
                hdeepinfo->recombination_control->Scr[hdeepinfo->recombination_control->Scr_size] = hdeepinfo->recombination_control->mod_CR[my_id];
                hdeepinfo->recombination_control->Scr_size++;
                hdeepinfo->recombination_control->Sf[hdeepinfo->recombination_control->Sf_size] = hdeepinfo->recombination_control->mod_F[my_id];
                hdeepinfo->recombination_control->Sf_size++;
                add_individ(hdeepinfo->archive, my_trial, hdeepinfo->hevalctrl->eval_target->array_size, hdeepinfo->hevalctrl->eval_target->precond_size);
            }
            dp_individ_copy_values(my_trial, my_individ);
            my_trial->age++;

        }
    } else if ( -1 != dp_evaluation_individ_compare((const void *)(&my_individ), (const void *)(&my_trial), (void*)(hdeepinfo->hevalctrl)) ) {
        //добавлять в архив здесь
        if (hdeepinfo->recombination_control->strategy == DE_3_pBest || hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
            hdeepinfo->recombination_control->Scr[hdeepinfo->recombination_control->Scr_size] = hdeepinfo->recombination_control->mod_CR[my_id];
            hdeepinfo->recombination_control->Scr_size++;
            hdeepinfo->recombination_control->Sf[hdeepinfo->recombination_control->Sf_size] = hdeepinfo->recombination_control->mod_F[my_id];
            hdeepinfo->recombination_control->Sf_size++;
            add_individ(hdeepinfo->archive, my_trial, hdeepinfo->hevalctrl->eval_target->array_size, hdeepinfo->hevalctrl->eval_target->precond_size);
        }
        dp_individ_copy_values(my_trial, my_individ);
        my_trial->age++;
    }
}

void dp_deep_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    if ( hdeepinfo->max_threads > 0 ) {
        if (hdeepinfo->gthreadpool == NULL) {
            hdeepinfo->gthreadpool = g_thread_pool_new ((GFunc) dp_deep_step_func, (gpointer) hdeepinfo, hdeepinfo->max_threads, hdeepinfo->exclusive, &gerror);
        }
        if ( gerror != NULL ) {
            g_error("%s", gerror->message);
        }
        for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
            g_thread_pool_push (hdeepinfo->gthreadpool, GINT_TO_POINTER(individ_id + 1), &gerror);
            if ( gerror != NULL ) {
                g_error("%s", gerror->message);
            }
        }
        g_thread_pool_free (hdeepinfo->gthreadpool, immediate_stop, wait_finish);
        hdeepinfo->gthreadpool = NULL;
    } else {
        for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
            dp_deep_step_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
        }
    }

    dp_population_update(trial, 0, trial->size);
    if (hdeepinfo->recombination_control->strategy == DE_3_pBest || hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        dp_recombination_control_update_Mcr_and_Mf(hdeepinfo->recombination_control, population, trial);
        dp_recombination_control_renew_Scr_and_Sf(hdeepinfo->recombination_control);
    }
    if (hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        calculate_xamean(hdeepinfo->archive);
    }
    trial->iter = population->iter + 1;
    hdeepinfo->recombination_control->iter = trial->iter;
    hdeepinfo->population = trial;
    hdeepinfo->trial = population;
}

typedef struct duple{
    int el;
    int ind;
} duple;

void cycle(int mean, int begin, int end, duple* tmp, int* unique, int *unique_counter){
    for(int i = begin; i < end; i++){
        for (int j = 0; j < *unique_counter; j++){
            if(tmp[i].el == unique[j])
                tmp[i].ind = 3;
        }
        if(tmp[i].ind != 3){
            tmp[i].ind = mean;
            unique[*unique_counter] = tmp[i].el;
            (*unique_counter)++;
        }
    }
}

int compare(const duple *i, const duple *j){
    return i->ind > j->ind;
}

void delete_dupl(double* ind1, double* ind2, int N){
    int ind1_t[N], ind2_t[N];
    duple tmp[2 * N];
    int unique[2 * N];
    int unique_counter = 0;
    int d = 500, h = 999;

    for(int i = 0; i < N; i++){
        ind1_t[i] = (int)(ind1[i] * h / d);
        ind2_t[i] = (int)(ind2[i] * h / d);
    }

    for(int i = 0; i < N; i++){
        tmp[i].el = ind1_t[i];
        tmp[i + N].el = ind2_t[i];
        tmp[i].ind = tmp[i + N].ind = 0;
    }

    cycle(1, 0, N, tmp, unique, &unique_counter);
    cycle(2, N, 2 * N, tmp, unique, &unique_counter);

    qsort(tmp, 2 * N, sizeof (duple),\
          (int(*)(const void *, const void *)) compare);

    for(int i = 0; i < N; i++)
        ind1[i] = (double)((double)(tmp[i].el) * d / h);
}

void main_dedupl(DpIndivid* my_trial, DpPopulation* population, DpDeepInfo* hdeepinfo, int ind){
    DpEvaluationPoint** dedupl_ind = hdeepinfo->hevalctrl->eval->points;
    double tmp1[my_trial->size], tmp2[my_trial->size];
    int N = 0;
    for(int i = 0; i < my_trial->size; i++){
        if(dedupl_ind[i]->dedupl == 0){
            tmp1[N] = my_trial->x[i];
            tmp2[N] = population->individ[ind]->x[i];
            N++;
        }
    }

    delete_dupl(tmp1, tmp2, N);
    int counter = 0;

    for(int i = 0; i < my_trial->size; i++){
        if(dedupl_ind[i]->dedupl == 0){
            my_trial->x[i] = tmp1[counter];
            counter++;
        }
    }
}

void dp_deep_generate_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4, r5;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpArchive* archive = hdeepinfo->archive;
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    recombination_control->current_individ = my_id;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    r1 = ( ignore_cost == 0 ) ? population->imin : -1;
    do {
        r2 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r2 == my_id || r2 == r1 );
    do {
        r3 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r3 == my_id || r3 == r1 || r3 == r2 );
    do {
        r4 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r4 == my_id || r4 == r1 || r4 == r2 || r4 == r3 );
    do {
        r5 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r5 == my_id || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4);
    r1 = population->imin; // best individ
    // my_tabu = population->individ[r1];
    start_index = g_rand_int_range (hrand, 0, population->ind_size);
    end_index = population->ind_size;
    dp_individ_copy_values(my_trial, my_individ);
    my_trial->age = 0;
    my_trial->r1 = r1;
    my_trial->r2 = r2;
    my_trial->r3 = r3;
    my_trial->r4 = r4;
    my_trial->r5 = r5;
    my_trial->moves = 0;
    my_trial->failures = 0;
    my_trial->grads = 0;

    DifferenceVector* vectorRead;
    DifferenceVector* vectorWrite;

	if (recombination_control->use_archive_prob > 0 &&
        recombination_control->strategy != DE_3_pBest &&
        recombination_control->strategy != DE_3_Amean &&
        recombination_control->strategy != DE_3_dual_adaptive) {
		vectorRead = archive->difference_vectors[my_id];
    	vectorWrite = archive->difference_vectors[hdeepinfo->population_size + my_id];
	}

    if (recombination_control->strategy == DE_3_pBest) {
        dp_individ_recombination(recombination_control, hrand, my_trial,
                                population->individ[r1], population->individ[r2], population->individ[r3], population->individ[population->xpbest_idx], population->individ[r5],
                                start_index, end_index, vectorWrite, vectorRead);
    } else if (recombination_control->strategy == DE_3_Amean) {
        vectorRead = archive->xamean;
        dp_individ_recombination(recombination_control, hrand, my_trial,
                                population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], population->individ[r5],
                                start_index, end_index, vectorWrite, vectorRead);
    } else if (recombination_control->strategy == DE_3_dual_adaptive) {
        vectorRead = archive->xamean;
        dp_individ_recombination(recombination_control, hrand, my_trial,
            population->individ[r1], population->individ[r2], population->individ[r3], population->individ[population->xpbest_idx], population->individ[r5],
            start_index, end_index, vectorWrite, vectorRead);
    } else {
        dp_individ_recombination(recombination_control, hrand, my_trial,
            population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], population->individ[r5],
            start_index, end_index, vectorWrite, vectorRead);
    }

    my_trial->cost = G_MAXDOUBLE;
}

void dp_deep_generate_dd_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4, r5;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    r1 = ( ignore_cost == 0 ) ? population->imin : -1;
    do {
        r2 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r2 == my_id || r2 == r1 );
    do {
        r3 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r3 == my_id || r3 == r1 || r3 == r2 );
    do {
        r4 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r4 == my_id || r4 == r1 || r4 == r2 || r4 == r3 );
    do {
        r5 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r5 == my_id || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4);
    r1 = population->imin;
    my_tabu = population->individ[r1];
    start_index = g_rand_int_range (hrand, 0, population->ind_size);
    end_index = population->ind_size;
    dp_individ_copy_values(my_trial, my_individ);
    my_trial->age = 0;
    my_trial->r1 = r1;
    my_trial->r2 = r2;
    my_trial->r3 = r3;
    my_trial->r4 = r4;
    my_trial->r5 = r5;
    my_trial->moves = 0;
    my_trial->failures = 0;
    my_trial->grads = 0;

    // dp_individ_recombination(recombination_control, hrand, my_trial, population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], start_index, end_index);
    DifferenceVector* vectorRead;
    DifferenceVector* vectorWrite;

	if (recombination_control->use_archive_prob > 0) {
		vectorRead = hdeepinfo->archive->difference_vectors[my_id];
    	vectorWrite = hdeepinfo->archive->difference_vectors[hdeepinfo->population_size + my_id];
	}

    dp_individ_recombination(recombination_control, hrand, my_trial,
                             population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], population->individ[r5],
                             start_index, end_index, vectorWrite, vectorRead);

    //delete duple

    int ind = g_rand_int_range(hrand, 0, population->size);
    while(ind == my_id || ind == my_trial->r1 || ind == my_trial->r2 || \
          ind == my_trial->r3 || ind == my_trial->r4) {
        ind = g_rand_int_range(hrand, 0, population->size);
    }

    main_dedupl(my_trial, population, hdeepinfo, ind);

    //end del duple

    my_trial->cost = G_MAXDOUBLE;
}

void dp_deep_generate_ca_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4;
    int start_index, end_index;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];

    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;

    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    r1 = population->ifmax;
    r2 = my_id;
    r3 = my_id;
    r4 = my_id;
    start_index = my_id;
    end_index = hdeepinfo->es_lambda;
    my_trial->age = 0;
    my_trial->r1 = r1;
    my_trial->r2 = r2;
    my_trial->r3 = r3;
    my_trial->r4 = r4;
    my_trial->moves = 0;
    my_trial->failures = 0;
    my_trial->grads = 0;
    dp_individ_recombination_ca(recombination_control, hrand, my_trial, population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], start_index, end_index);
    my_trial->cost = G_MAXDOUBLE;
}

void dp_deep_generate_ac_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4;
    int start_index, end_index;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];

    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;

    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    /*	r1 = population->ifmax;
    r2 = my_id;
    r3 = my_id;
    r4 = my_id;*/
    r1 = (my_individ->r1 < 0) ? population->imin : my_individ->r1;
    r2 = (my_individ->r2 < 0) ? population->imin : my_individ->r2;
    r3 = (my_individ->r3 < 0) ? population->imin : my_individ->r3;
    r4 = (my_individ->r4 < 0) ? population->imin : my_individ->r4;
    start_index = my_id;
    end_index = hdeepinfo->es_lambda;
    /* The following code checks if
 * the inidvid has parents
 */
    my_trial->age = 0;
    my_trial->r1 = r1;
    my_trial->r2 = r2;
    my_trial->r3 = r3;
    my_trial->r4 = r4;
    my_trial->moves = 0;
    my_trial->failures = 0;
    my_trial->grads = 0;
    dp_individ_recombination_ac(recombination_control, hrand, my_trial, population->individ[r1], population->individ[r2], population->individ[r3], population->individ[r4], start_index, end_index);
    my_trial->cost = G_MAXDOUBLE;
}

void dp_deep_evaluate_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    r1 = population->imin;
    my_tabu = population->individ[r1];
    my_trial->age = 0;
    g_mutex_lock( &(my_trial->m) );
    my_trial->status = 1;
    g_mutex_unlock( &(my_trial->m) );
    dp_evaluation_individ_evaluate(hdeepinfo->hevalctrl, my_trial, my_tabu, my_id, my_tabu->cost);
    g_mutex_lock( &(my_trial->m) );
    my_trial->status = 0;
    g_mutex_unlock( &(my_trial->m) );
}

int dp_deep_evaluate_status (gpointer data, gpointer user_data)
{
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1, s;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    g_mutex_lock( &(my_trial->m) );
    s = my_trial->status;
    g_mutex_unlock( &(my_trial->m) );
    return s;
}

void dp_deep_select_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    r1 = population->imin;
    my_tabu = population->individ[r1];
    if ( ignore_cost == 0 && my_id == population->imin ) {
        if ( my_trial->cost >= my_individ->cost ) {
            population->individ[my_id] = my_trial;
            trial->individ[my_id] = my_individ;
            my_individ->age++;
        }
    } else if ( -1 != dp_evaluation_individ_compare((const void *)(&my_individ), (const void *)(&my_trial), (void*)(hdeepinfo->hevalctrl)) ) {
        /* individ is accepted if -1 is returned */
        population->individ[my_id] = my_trial;
        trial->individ[my_id] = my_individ; /* accepted */
        my_individ->age++;
    }
}

void dp_de_select_func (gpointer data, gpointer user_data)
{
    int r1, r2, r3, r4;
    int start_index, end_index;
    DpIndivid*my_tabu;
    DpDeepInfo*hdeepinfo = (DpDeepInfo*)user_data;
    int my_id = GPOINTER_TO_INT(data) - 1;
    DpPopulation*trial = hdeepinfo->trial;
    DpIndivid*my_trial = trial->individ[my_id];
    DpPopulation*population = hdeepinfo->population;
    DpIndivid*my_individ = population->individ[my_id];
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    int ignore_cost = hdeepinfo->hevalctrl->eval_target->ignore_cost;
    r1 = population->imin;
    my_tabu = population->individ[r1];
    if ( my_trial->cost >= my_individ->cost ) {
        population->individ[my_id] = my_trial;
        trial->individ[my_id] = my_individ;
        my_individ->age++;
        if (my_trial->useWriteVector)
            hdeepinfo->archive->difference_vectors[my_id] = hdeepinfo->archive->difference_vectors[my_id + hdeepinfo->population_size];
    }
    else{
        if (my_trial->useWriteVector)
            hdeepinfo->archive->difference_vectors[my_id] = hdeepinfo->archive->difference_vectors[my_id + hdeepinfo->population_size];
    }
}

void dp_deep_generate_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;

    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_deep_generate_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }


}

void dp_deep_generate_dd_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_deep_generate_dd_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }
}

void dp_deep_generate_ca_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_deep_generate_ca_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }
}

void dp_deep_generate_ac_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_deep_generate_ac_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }
}

void dp_deep_evaluate_step(DpDeepInfo*hdeepinfo)
{
    GMainContext *gcontext = g_main_context_default();
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    gulong microseconds = G_USEC_PER_SEC / 1000;
    int notdone, currstatus;
    if ( hdeepinfo->max_threads > 0 ) {
        for ( individ_id = population->slice_a; individ_id < population->slice_b; individ_id++ ) {
            g_thread_pool_push (hdeepinfo->gthreadpool, GINT_TO_POINTER(individ_id + 1), &gerror);
            if ( gerror != NULL ) {
                g_error("%s", gerror->message);
            }
        }
        notdone = 1;
        //        while(g_thread_pool_unprocessed (hdeepinfo->gthreadpool) > 0) {
        while(notdone == 1 || g_thread_pool_unprocessed (hdeepinfo->gthreadpool) > 0) {
            g_main_context_iteration(gcontext, FALSE);
            g_usleep (microseconds);
            notdone = 0;
            for ( individ_id = population->slice_a; individ_id < population->slice_b; individ_id++ ) {
                currstatus = dp_deep_evaluate_status (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
                if ( currstatus == 1 ) {
                    notdone = 1;
                }
            }
        }
    } else {
        for ( individ_id = population->slice_a; individ_id < population->slice_b; individ_id++ ) {
            dp_deep_evaluate_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
        }
    }
}

/* can be moved to ppopulaiton.c */
void dp_deep_select_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_deep_select_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }
    dp_population_update(trial, 0, trial->size);
    trial->iter = population->iter + 1;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        int r4 = trial->individ[individ_id]->r4;
        if (r4 > -1 && trial->individ[individ_id]->age > 0 && trial->individ[r4]->age > 0) {
            trial->individ[r4]->failures++;
        }
    }
    hdeepinfo->population = trial; /* accepted */
    hdeepinfo->trial = population;
}

void dp_de_select_step(DpDeepInfo*hdeepinfo)
{
    int individ_id;
    gboolean immediate_stop = FALSE;
    gboolean wait_finish = TRUE;
    DpPopulation*population = hdeepinfo->population;
    DpPopulation*trial = hdeepinfo->trial;
    GError *gerror = NULL;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        dp_de_select_func (GINT_TO_POINTER(individ_id + 1), (gpointer) hdeepinfo);
    }
    dp_population_update(trial, 0, trial->size);
    trial->iter = population->iter + 1;
    for ( individ_id = 0; individ_id < population->size; individ_id++ ) {
        int r4 = trial->individ[individ_id]->r4;
        if (r4 > -1 && trial->individ[individ_id]->age > 0 && trial->individ[r4]->age > 0) {
            trial->individ[r4]->failures++;
        }
    }
    hdeepinfo->population = trial; /* accepted */
    hdeepinfo->trial = population;
}

void dp_deep_accept_step(DpDeepInfo*hdeepinfo, double*value)
{
    if (hdeepinfo->mean_cost == 0) {
        *value = hdeepinfo->population->dmin;
    } else 	if (hdeepinfo->mean_cost == 1) {
        *value = hdeepinfo->population->fmean;
    }
}

void dp_deep_post(DpDeepInfo*hdeepinfo)
{
    DpPopulation*pop = hdeepinfo->population;
    double *x = pop->individ[pop->imin]->x;
    dp_evaluation_individ_copy(hdeepinfo->hevalctrl, pop->individ[pop->imin], pop->individ[pop->imin], pop->imin, pop->dmin);
}

void dp_deep_post_evaluate(DpDeepInfo*hdeepinfo)
{
    DpPopulation*pop = hdeepinfo->population;
    dp_evaluation_individ_evaluate_precond(hdeepinfo->hevalctrl, pop->individ[pop->imin], pop->individ[pop->imin], pop->imin, pop->dmin);
}

void dp_deep_update_step(DpDeepInfo*hdeepinfo)
{
    dp_population_update(hdeepinfo->population, 0, hdeepinfo->population->size);
    if (hdeepinfo->recombination_control->strategy == DE_3_pBest || hdeepinfo->recombination_control->strategy == DE_3_Amean || hdeepinfo->recombination_control->strategy == DE_3_dual_adaptive) {
        dp_recombination_control_update_modified(hdeepinfo->hevalctrl->hrand, hdeepinfo->recombination_control, hdeepinfo->population);
    }
    dp_recombination_control_update(hdeepinfo->recombination_control, hdeepinfo->hevalctrl->hrand, hdeepinfo->population, 0, hdeepinfo->population->size);
}

int cmpfunc(const void * a, const void * b) {
        return (*(int*)a - *(int*)b);
    }

void change_pointers(DpPopulation*population, const int reduce_number){
    int *last_indices;
    int len = population->cur_size;
    int i = 0;
    int k;
    int j = reduce_number - 1;
    DpIndivid* tmp;

    last_indices = &population->cost_ascending[len - reduce_number];
    qsort(last_indices, reduce_number, sizeof(int), cmpfunc);

    // Перемещение последних reduce_number индексов в конец списка individ
    for (k = 0; k < reduce_number; k++){
        if (last_indices[i] >= len - reduce_number)
            break;

        if (last_indices[j] != len - i - 1){
            tmp = population->individ[len - i - 1];
            population->individ[len - i -  1] = population->individ[last_indices[i]];
            population->individ[last_indices[i]] = tmp;
            i++;
        }
        else j--;
    }
}

int* generate_rand_ind(DpDeepInfo* hdeepinfo, int min){
    int r2, r3, r4, r5;
    int* rand_ind = (int*)malloc(4*sizeof(int));
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    do {
        r2 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r2 == min );
    do {
        r3 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r3 == min || r3 == r2 );
    do {
        r4 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r4 == min || r4 == r2 || r4 == r3 );
    do {
        r5 = g_rand_int_range (hrand, 0, hdeepinfo->population_size);
    } while ( r5 == min || r5 == r2 || r5 == r3 || r5 == r4);
    rand_ind[0] = r2;
    rand_ind[1] = r3;
    rand_ind[2] = r4;
    rand_ind[3] = r5;
    return rand_ind;
}

void dp_deep_pop_size_change(DpDeepInfo*hdeepinfo){

    if (hdeepinfo->population->dmin < hdeepinfo->fmin){
        hdeepinfo->b = 1;
        hdeepinfo->fmin = hdeepinfo->population->dmin;
    }
    else{
        hdeepinfo->w = 1;
        hdeepinfo->st += 1;
    }
    if (hdeepinfo->population->cur_size <= hdeepinfo->lower_bound){
        hdeepinfo->lb += 1;
    }
    if ((hdeepinfo->w == 1) && (hdeepinfo->st <= hdeepinfo->R)){
        dp_deep_redstrat(hdeepinfo);
        hdeepinfo->w = 0;
    }
    if ((hdeepinfo->b == 1) && (hdeepinfo->population->cur_size != hdeepinfo->population_max_size)){
        dp_deep_augstrat_1(hdeepinfo);
        hdeepinfo->b = 0;
        hdeepinfo->st = 0;
    }
    if (((hdeepinfo->st > hdeepinfo->R) || (hdeepinfo->lb > hdeepinfo->R)) && (hdeepinfo->population->cur_size != hdeepinfo->population_max_size)){
        dp_deep_augstrat_2(hdeepinfo);
        hdeepinfo->st = 0;
        hdeepinfo->lb = 0;
    }
}

// стратегии изменения размера популяции
void dp_deep_redstrat(DpDeepInfo*hdeepinfo){
    const int reduce_number = (int)floor(hdeepinfo->population->cur_size * hdeepinfo->s);
    // перестановка индексов - худших в конец
    change_pointers(hdeepinfo->population, reduce_number);
    for (int i = 0; i < hdeepinfo->population->cur_size; i++) {
        hdeepinfo->trial->individ[i]->user_data = hdeepinfo->population->individ[i]->user_data;
    }
    hdeepinfo->population_size -= reduce_number;
    hdeepinfo->population->cur_size -= reduce_number;
    hdeepinfo->population->slice_b -= reduce_number;
    dp_population_update(hdeepinfo->population, 0, hdeepinfo->population->cur_size);
}

void dp_deep_augstrat_1(DpDeepInfo*hdeepinfo){
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    DpPopulation*population = hdeepinfo->population;
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    DifferenceVector* vectorRead;
    DifferenceVector* vectorWrite;
    int start_index = g_rand_int_range (hrand, 0, hdeepinfo->population->ind_size);
    int end_index = population->ind_size;
    int* indices = generate_rand_ind(hdeepinfo, population->imin);
    DpRecombinationStrategy prev_strat = recombination_control->strategy;
    recombination_control->strategy = DE_4_bin;
    dp_individ_recombination(recombination_control, hrand, population->individ[population->cur_size],
                             population->individ[population->imin], population->individ[indices[0]], population->individ[indices[1]], population->individ[indices[2]], population->individ[indices[3]],
                             start_index, end_index, vectorWrite, vectorRead);
    hdeepinfo->population_size += 1;
    population->cur_size += 1;
    population->slice_b += 1;
    recombination_control->strategy = prev_strat;
    dp_population_update(population, 0, population->cur_size);
    free(indices);
}

void dp_deep_augstrat_2(DpDeepInfo*hdeepinfo){
    const int add_size = (int)floor(hdeepinfo->population->cur_size * hdeepinfo->s);
    GRand*hrand = hdeepinfo->hevalctrl->hrand;
    DpPopulation*population = hdeepinfo->population;
    DpRecombinationControl *recombination_control = hdeepinfo->recombination_control;
    DifferenceVector* vectorRead;
    DifferenceVector* vectorWrite;
    int start_index = g_rand_int_range (hrand, 0, hdeepinfo->population->ind_size);
    int end_index = population->ind_size;
    int* indices = generate_rand_ind(hdeepinfo, population->imin);
    DpRecombinationStrategy prev_strat = recombination_control->strategy;
    recombination_control->strategy = DE_3_bin;
    for (int i = 0; i < add_size; i++){
        dp_individ_recombination(recombination_control, hrand, population->individ[population->cur_size+i],
                                population->individ[population->imin], population->individ[indices[0]], population->individ[indices[1]], population->individ[indices[2]], population->individ[indices[3]],
                                start_index, end_index, vectorWrite, vectorRead);
    }
    hdeepinfo->population_size += add_size;
    population->cur_size += add_size;
    population->slice_b += add_size;
    recombination_control->strategy = prev_strat;
    dp_population_update(population, 0, population->cur_size);
    free(indices);
}

void dp_init_xpbest_param(DpDeepInfo*hdeepinfo) {
    int pbest_idx;

    hdeepinfo->population->p_param = hdeepinfo->p_param;
    pbest_idx = floor(hdeepinfo->population->p_param * hdeepinfo->population->size);
    hdeepinfo->population->xpbest_idx = hdeepinfo->population->cost_ascending[pbest_idx];

    hdeepinfo->trial->p_param = hdeepinfo->p_param;
    pbest_idx = floor(hdeepinfo->trial->p_param * hdeepinfo->trial->size);
    hdeepinfo->trial->xpbest_idx = hdeepinfo->trial->cost_ascending[pbest_idx];

    hdeepinfo->popunion->p_param = hdeepinfo->p_param;
    pbest_idx = floor(hdeepinfo->popunion->p_param * hdeepinfo->popunion->size);
    hdeepinfo->popunion->xpbest_idx = hdeepinfo->popunion->cost_ascending[pbest_idx];
}

void dp_init_xamean_param(DpDeepInfo*hdeepinfo) {
    if (hdeepinfo->recombination_control->use_archive_prob > 0) {
        hdeepinfo->archive->e = hdeepinfo->e_param;
        hdeepinfo->archive->m = round(hdeepinfo->archive->e * hdeepinfo->archive->last_index);
    }
}