/***************************************************************************
 *            dprecombination.h
 *
 *  Fri Mar 23 15:39:19 2012
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

/*  Recombination-mutation functions for DEEP
 *
 *
*/

#ifndef _DP_RECOMBINATION_H
#define _DP_RECOMBINATION_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <glib.h>
#include "dparchive.h"

#define FW_PARAM 0.7

typedef enum DpRecombinationStrategy {
	Simple,
	DE_3_bin,
    DE_4_bin,
	DE_3_exp,
	DE_3_bin_T,
	DE_3_bin_A,
	DE_3_exp_T,
	DE_3_bin_rand,
	DE_3_exp_rand,
	DE_3_bin_self,
	DE_3_exp_self,
	DE_3_bin_rand_T,
	DE_3_exp_rand_T,
	DE_3_bin_rand_phi,
	DE_3_exp_rand_phi,
	DE_3_bin_self_phi,
	DE_3_exp_self_phi,
	DE_3_bin_rand_D,
	DE_3_exp_rand_D,
	DE_3_bin_self_D,
	DE_3_exp_self_D,
	DE_3_bin_self_T,
    DE_3_exp_self_T,
    DE_3_triangular,
	DE_3_pBest,
	DE_3_Amean,
	DE_3_dual_adaptive,
	MaxStrategy
} DpRecombinationStrategy;

typedef struct DpRecombinationControl {
	DpRecombinationStrategy strategy;
	int size;
	int pop_size;
	double gamma;
    double *f;  // Array of params
    double *p;  //
    double *c;  //
    double *v;  //
	double f_inf;
	double f_supp;
	double p_inf;
	double p_supp;
	// параметры модификаций
	double a_param;
	int iter;
	int max_iters;
	//для модификации параметров контроля
	double *mod_CR; //массив индивидуальных пар-ров
	double *mod_F;  //контроля для каждого ИНДИВИДА
	int current_individ; //индекс индивида находящегося на рекомбинации
	double *Mf;
	double *Mcr;
	double *Sf;
	int Sf_size;
	double *Scr;
	int Scr_size;
	//
    int adjust;  // Adaptaion flag
    int toggle; // Adaptation of params on even steps
	double use_archive_prob; // setting this probability to 0 will switch archive off
	//переработаны FEs (functions evaluations) и MaxFEs (maximum functions evaluations)
	int FEs;
	int MaxFEs;
} DpRecombinationControl;
// те, которые используют input_1 - best-ы, методы без self & rand

void ind_triand_init(DpPopulation *population, DpRecombinationControl *control, GRand *hrand);

// void dp_individ_recombination(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index);

void dp_individ_recombination(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, DpIndivid*input_5, int start_index, int end_index, DifferenceVector* vectorWrite, DifferenceVector* vectorRead);

void dp_individ_recombination_ca(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index);

void dp_individ_recombination_ac(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index);

DpRecombinationControl*dp_recombination_control_init(GRand*hrand,  DpPopulation*pop, GKeyFile*gkf, gchar*groupname);

void dp_recombination_control_renew(DpRecombinationControl*rc, GRand*hrand, DpPopulation*pop);

void dp_recombination_control_update(DpRecombinationControl*rc, GRand*hrand, DpPopulation*pop, int start_index, int end_index);

void dp_recombination_control_save(FILE*fp, DpRecombinationControl*rc);

void dp_recombination_control_load(FILE*fp, DpRecombinationControl*rc);

double dp_rand_normal(GRand*hrand, double mean, double stddev);

double dp_rand_Cauchy(GRand*hrand, double loc, double scale);

void dp_recombination_control_update_Mcr_and_Mf(DpRecombinationControl*rc, DpPopulation*pop, DpPopulation*trial);

void dp_recombination_control_renew_Scr_and_Sf(DpRecombinationControl*rc);

// #include "dparchive.h"

#ifdef __cplusplus
}
#endif

#endif /* _DP_RECOMBINATION_H */
