/***************************************************************************
 *            dpopt.c
 *
 *  Fri Mar 23 15:49:34 2012
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
#include <glib/gprintf.h>
#include "dploop.h"
#include "dpopt.h"
#include "dpdeep.h"
#include "dposda.h"

#ifdef GIO_STANDALONE_SOURCE
#define g_strcmp0(str1, str2) strcmp(str1, str2)
#endif

DpOpt *dp_opt_init(DpEvaluation*heval, DpTarget*htarget, int world_id, int world_count,char*filename, DpOptStopType stop_type, double criterion, int stop_count, int pareto_all, int precision)
{
	DpOpt *hopt;
	GList *funcs = NULL;
	GList *run_once_before = NULL;
	GList *run_once_after = NULL;
	hopt = (DpOpt*)malloc(sizeof(DpOpt));
	hopt->world_id = world_id;
	hopt->world_count = world_count;
	hopt->filename = g_strdup (filename);
	hopt->stop_type = stop_type;
	hopt->criterion = criterion;
	hopt->stop_count = stop_count;
	hopt->logname = g_strdup_printf( "%s.log", hopt->filename, hopt->world_id);
	hopt->chkname = g_strdup_printf( "%s.chk", hopt->filename, hopt->world_id);
	hopt->tstname = g_strdup_printf( "%s.time", hopt->filename);
	hopt->htarget = htarget;
	hopt->heval = heval;
	hopt->hloop = dp_loop_new(NULL, NULL, NULL);
	hopt->opt_type = H_OPT_NONE;
	hopt->method_info = NULL;
	hopt->pareto_all = pareto_all;
	hopt->precision = precision;
	hopt->debug = 0;
	hopt->delay = 0;
	hopt->logdepth = 0;
	hopt->monitor = -1;
	hopt->keep_order = 0;
	hopt->cancel_flags = 0;
	hopt->cancel_count = 0;
	hopt->cancel_iter = 0;
	hopt->cancel_time = 0;
	hopt->cancel_prop = 0;
	hopt->cancel_abs = 0;
	hopt->cancel_score = 0;
	hopt->cancel_counter = 0;
	hopt->advance_kounter = 0;
	return hopt;
}

void dp_opt_add_func(DpOpt *hopt, DpLoopFunc func, int tau_flag, DpOptType opt_type, int order, gpointer method_info)
{
	DpLoopRunFunc *run_func;
	DpLoop *hloop = hopt->hloop;
	run_func = dp_loop_run_func_new(tau_flag, func, (gpointer) hopt);
	if ( hopt->opt_type == H_OPT_NONE ) {
		hopt->opt_type = opt_type;
	}
	if ( hopt->method_info == NULL ) {
		hopt->method_info = method_info;
	}
	if ( order < 0 ) {
		dp_loop_add_func_before(hloop, run_func);
	} else if ( order == 0 ) {
		dp_loop_add_func(hloop, run_func);
	} else {
		dp_loop_add_func_after(hloop, run_func);
	}
}

void dp_opt_add_from_func_list(gchar**list, DpOpt *hopt, int order, GKeyFile*gkf, gchar*groupname, int world_id, DpEvaluation*heval, DpTarget*htarget, GError**err)
{
	DpDeepInfo*hdeepinfo = NULL;
	DpOsdaInfo*hosdainfo = NULL;
    gpointer method_info; // Make 3 separeate fields for each opt_type instacnce (if not init - init, else continue processing)
	int i, tau_flag, incr;
	DpOptType opt_type;
	incr = (order == 0) ? 2 : 1;
	tau_flag = 1;
	for ( i = 0; list[i]; i += incr ) {
		if (incr == 2) {
			tau_flag = g_strtod(list[i + 1], NULL);
		}
		if (hopt->world_id == 0) {
			if ( !g_strcmp0(list[i], "writelog") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_write_log, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "printlog") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_print_log, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "readlog") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_read_log, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "writestate") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_write_state, tau_flag, opt_type, order, method_info);
			} else if (!g_strcmp0(list[i], "readstate")) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_read_state, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "writepareto") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_write_pareto, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "evalpareto") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_evaluate_pareto_front, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "selpareto") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_select_pareto_front, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "sortpareto") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_sort_pareto_front, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "permutepop") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_permute, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "cr2cost") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_cr2cost, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "checkstop") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_check_stop, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "initstop") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_init_stop, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "checkcancel") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_check_cancel, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "initcancel") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_init_cancel, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "substitute") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_substitute, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "dpupdate") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_deep_update, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "optpost") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_post, tau_flag, opt_type, order, method_info);
			} else if ( !g_strcmp0(list[i], "optposteval") ) {
				opt_type = H_OPT_NONE;
				method_info = NULL;
				dp_opt_add_func(hopt, dp_opt_post_evaluate, tau_flag, opt_type, order, method_info);
			}
		}
		if ( !g_strcmp0(list[i], "rotatetarget") ) {
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_rotate_target, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "mpidistribute") ) {
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_mpi_distribute, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "mpigather") ) {
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_mpi_gather, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "deep") ) {
			opt_type = H_OPT_DEEP;
			if (!hdeepinfo) {
				hdeepinfo = dp_deep_info_init(heval, htarget, world_id, gkf, groupname);
			}
			method_info = (gpointer) hdeepinfo;
			dp_opt_add_func(hopt, dp_opt_deep, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "gdeep") ) {
			opt_type = H_OPT_DEEP;
			if (!hdeepinfo) {
				hdeepinfo = dp_deep_info_init(heval, htarget, world_id, gkf, groupname);
			}
			method_info = (gpointer) hdeepinfo;
			dp_opt_add_func(hopt, dp_opt_deep_generate, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "gdddeep") ) {
			opt_type = H_OPT_DEEP;
			if (!hdeepinfo) {
				hdeepinfo = dp_deep_info_init(heval, htarget, world_id, gkf, groupname);
			}
			method_info = (gpointer) hdeepinfo;
			dp_opt_add_func(hopt, dp_opt_deep_generate_dd, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "gcadeep") ) {
			opt_type = H_OPT_DEEP;
			if (!hdeepinfo) {
				hdeepinfo = dp_deep_info_init(heval, htarget, world_id, gkf, groupname);
			}
			method_info = (gpointer) hdeepinfo;
			dp_opt_add_func(hopt, dp_opt_deep_generate_ca, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "gacdeep") ) {
			opt_type = H_OPT_DEEP;
			if (!hdeepinfo) {
				hdeepinfo = dp_deep_info_init(heval, htarget, world_id, gkf, groupname);
			}
			method_info = (gpointer) hdeepinfo;
			dp_opt_add_func(hopt, dp_opt_deep_generate_ac, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "edeep") ) {
            opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_deep_evaluate, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "sdeep") ) {
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_deep_select, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "selde") ) {
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_de_select, tau_flag, opt_type, order, method_info);
		} else if ( !g_strcmp0(list[i], "osda") ) {
			opt_type = H_OPT_OSDA;
			method_info = (gpointer) hosdainfo;
			dp_opt_add_func(hopt, dp_opt_osda, tau_flag, opt_type, order, method_info);
		} else if (!g_strcmp0(list[i], "pdeep"))
		{
			opt_type = H_OPT_NONE;
			method_info = NULL;
			dp_opt_add_func(hopt, dp_opt_popsize, tau_flag, opt_type, order, method_info);
		}
	}
}


void dp_opt_run(DpOpt *hopt)
{
	DpLoop *hloop = hopt->hloop;
	dp_loop_zero_counters(hopt->hloop);
	dp_loop_run (hloop);
}

void dp_opt_monitor(DpOpt *hopt, int monitor, GError**gerror)
{
	DpLoop *hloop;
	GList *funcs = NULL;
	GList *run_once_before = NULL;
	GList *run_once_after = NULL;
	DpLoopRunFunc *func;

	func = (DpLoopRunFunc *)malloc(sizeof(DpLoopRunFunc));
	func->tau_flag = 1;
	func->user_data = (void*)hopt;
	func->func = (DpLoopFunc)dp_read_log;
	funcs = g_list_append (funcs, (void*)func);
	func = (DpLoopRunFunc *)malloc(sizeof(DpLoopRunFunc));
	func->tau_flag = 1;
	func->user_data = (void*)hopt;
	func->func = (DpLoopFunc)dp_write_tst;
	funcs = g_list_append (funcs, (void*)func);

	hloop = dp_loop_new(run_once_before, funcs, run_once_after);
	dp_loop_zero_counters(hloop);
	hopt->hloop = hloop;
	hopt->monitor = monitor;
	dp_loop_run (hloop);
}

DpLoopExitCode dp_opt_popsize(DpLoop*hloop, gpointer user_data)
{
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	DpLoopExitCode exit_code = DP_LOOP_EXIT_NOEXIT;
	dp_deep_pop_size_change(hdeepinfo);
	return exit_code;
}

DpLoopExitCode dp_opt_init_stop(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_deep_accept_step(hdeepinfo, &(hopt->cost));
			hopt->old_cost = hopt->cost;
			hopt->stop_counter = 0;
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
			dp_osda_accept_step(hosdainfo, &(hopt->cost));
			hopt->old_cost = hopt->cost;
			hopt->stop_counter = 0;
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_check_stop(DpLoop*hloop, gpointer user_data) // don't use
{
	DpOpt*hopt = (DpOpt*)user_data;
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	double delta;
	if ( hopt->stop_type != H_OPT_ABSOLUTE_TIME && hopt->stop_type != H_OPT_ABSOLUTE_ITER ) {
		switch( hopt->stop_type ) {
			case H_OPT_PROPORTIONAL_STOP:
				delta = ( hopt->old_cost - hopt->cost );
				if ( hopt->old_cost != 0 ) {
					delta /= hopt->old_cost;
				}
			break;
			case H_OPT_ABSOLUTE_STOP:
				delta = hopt->old_cost - hopt->cost;
			break;
			case H_OPT_ABSOLUTE_SCORE:
				delta = hopt->cost;
			break;
			default:
				hloop->exit_str = g_strdup_printf ( "dp_opt_checkstop: unknown %d", (int)(hopt->stop_type));
				return DP_LOOP_EXIT_ERROR;
		}
		if ( fabs(delta) < hopt->criterion ) {
			hopt->stop_counter++;
		} else {
			hopt->stop_counter = 0;
			hopt->old_cost = hopt->cost;
		}
		if ( hopt->stop_counter == hopt->stop_count) {
			hopt->stop_counter = 0;
			hopt->old_cost = hopt->cost;
			ret_val = DP_LOOP_EXIT_SUCCESS;
		}
	} else if ( hopt->stop_type == H_OPT_ABSOLUTE_TIME && hloop->w_time >= hopt->criterion ) {
		ret_val = DP_LOOP_EXIT_SUCCESS;
	} else if ( hopt->stop_type == H_OPT_ABSOLUTE_ITER && hloop->tau_counter >= hopt->criterion ) {
		ret_val = DP_LOOP_EXIT_SUCCESS;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_init_cancel(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_deep_accept_step(hdeepinfo, &(hopt->cost));
			hopt->old_cost = hopt->cost;
			hopt->cancel_counter = 0;
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
			dp_osda_accept_step(hosdainfo, &(hopt->cost));
			hopt->old_cost = hopt->cost;
			hopt->cancel_counter = 0;
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_check_cancel(DpLoop*hloop, gpointer user_data)
{
	DpOpt*hopt = (DpOpt*)user_data;
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	double delta = 0;
	double criterion = -1;
	if ( hopt->cancel_flags & DP_PROPORTIONAL_CANCEL ) {
		delta = ( hopt->old_cost - hopt->cost );
		if ( hopt->old_cost != 0 ) {
			delta /= hopt->old_cost;
		}
		criterion = hopt->cancel_prop;
		delta = (delta > 0) ? delta : -delta;
		if ( delta < criterion ) {
			hopt->cancel_counter++;
		} else {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
		}
		if ( hopt->cancel_counter == hopt->cancel_count) {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
			ret_val = DP_LOOP_EXIT_SUCCESS;
		}
	}
	if ( hopt->cancel_flags & DP_ABSOLUTE_CANCEL ) {
		delta = hopt->old_cost - hopt->cost;
		criterion = hopt->cancel_abs;
		delta = (delta > 0) ? delta : -delta;
		if ( delta < criterion ) {
			hopt->cancel_counter++;
		} else {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
		}
		if ( hopt->cancel_counter == hopt->cancel_count) {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
			ret_val = DP_LOOP_EXIT_SUCCESS;
		}
	}
	if ( hopt->cancel_flags & DP_ABSOLUTE_SCORE ) {
		delta = hopt->cost;
		criterion = hopt->cancel_score;
		delta = (delta > 0) ? delta : -delta;
		if ( delta < criterion ) {
			hopt->cancel_counter++;
		} else {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
		}
		if ( hopt->cancel_counter == hopt->cancel_count) {
			hopt->cancel_counter = 0;
			hopt->old_cost = hopt->cost;
			ret_val = DP_LOOP_EXIT_SUCCESS;
		}
	}
	if ( (hopt->cancel_flags & DP_ABSOLUTE_TIME) && (hloop->w_time >= hopt->cancel_time) ) {
		ret_val = DP_LOOP_EXIT_SUCCESS;
	}
	if ( (hopt->cancel_flags & DP_ABSOLUTE_ITER) && (hloop->tau_counter >= hopt->cancel_iter) ) {
		ret_val = DP_LOOP_EXIT_SUCCESS;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_deep(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	double energy_start;
	energy_start = hopt->cost;
	dp_deep_step(hdeepinfo);
	dp_deep_accept_step(hdeepinfo, &(hopt->cost));
	return ret_val;
}

DpLoopExitCode dp_opt_deep_generate(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	if ( hdeepinfo->selector == DpSelectorGenerated ) {
        return ret_val;
    }
	dp_deep_generate_step(hdeepinfo);
	hdeepinfo->selector = DpSelectorGenerated;
	hdeepinfo->last_method = 0;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_generate_dd(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	if ( hdeepinfo->selector == DpSelectorGenerated ) {
        return ret_val;
    }
	dp_deep_generate_dd_step(hdeepinfo);
	hdeepinfo->selector = DpSelectorGenerated;
	hdeepinfo->last_method = 0;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_generate_ca(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	DpPopulation*population = hdeepinfo->population;
	DpPopulation*trial = hdeepinfo->trial;
	DpIndivid*tmp_ind;
	int b1, i, j, k, ind;
	if ( hdeepinfo->selector == DpSelectorGenerated ) {
        return ret_val;
    }
    g_debug("dp_opt_deep_generate_ca tau:%d selector:%d counter:%d delay:%d stop_counter:%d", hloop->tau_counter, hdeepinfo->selector, hdeepinfo->hevalctrl->kounter, hopt->delay, hopt->stop_counter);
	b1 = (hdeepinfo->ca_flag == 0) ? hdeepinfo->es_lambda : (int)(population->size / hdeepinfo->es_lambda + 0.5);
	k = 0;
	for (i = 0; i < hdeepinfo->es_lambda; i++) {
		ind = population->cost_ascending[i];
		tmp_ind = population->individ[ind];
		for (j = 0; j < b1 && k < population->size; j++, k++) {
			dp_individ_copy_values(trial->individ[k], tmp_ind);
		}
	}
	ind = b1 - 1;
	trial->fmax = population->individ[ind]->cost;
	trial->ifmax = ind;
	trial->imin = 0;
	for (i = 0; i < b1; i++) {
		if (i * b1 > population->size - 1) break;
		tmp_ind = trial->individ[i * b1];
		for (j = i; j < population->size; j += b1) {
			dp_individ_copy_values(population->individ[j], tmp_ind);
		}
	}
	hdeepinfo->population = trial;
	hdeepinfo->trial = population;
	dp_deep_generate_ca_step(hdeepinfo);
	hdeepinfo->selector = DpSelectorGenerated;
	hdeepinfo->last_method = 1;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_generate_ac(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	if ( hdeepinfo->selector == DpSelectorGenerated ) {
        return ret_val;
    }
    if (hopt->keep_order == 1 && hdeepinfo->last_method != 1) {
        return ret_val;
	}
	dp_deep_generate_ac_step(hdeepinfo);
	hdeepinfo->selector = DpSelectorGenerated;
	hdeepinfo->last_method = 2;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_evaluate(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	dp_deep_evaluate_step(hdeepinfo);
	hdeepinfo->selector = DpSelectorEvaluated;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_select(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	if ( hdeepinfo->selector == DpSelectorSelected ) {
        return ret_val;
    }
	dp_deep_select_step(hdeepinfo);
	dp_deep_accept_step(hdeepinfo, &(hopt->cost));
	hdeepinfo->selector = DpSelectorSelected;
	return ret_val;
}

DpLoopExitCode dp_opt_de_select(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	if ( hdeepinfo->selector == DpSelectorSelected ) {
        return ret_val;
    }
	dp_de_select_step(hdeepinfo);
	dp_deep_accept_step(hdeepinfo, &(hopt->cost));
	hdeepinfo->selector = DpSelectorSelected;
	return ret_val;
}

DpLoopExitCode dp_opt_deep_update(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	dp_deep_update_step(hdeepinfo);
	return ret_val;
}

DpLoopExitCode dp_write_log(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	DpTarget*htarget = (DpTarget*)(hopt->htarget);
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	FILE*fp;
	int i, kounter;
	int precision = hopt->precision;
	fp = fopen(hopt->logname, "a");
	if ( !fp ) {
		hloop->exit_str = g_strdup_printf ( "dp_write_log: can't open %s", hopt->logname);
		return DP_LOOP_EXIT_ERROR;
	}
	fprintf(fp, "wtime:%e tau:%d freeze:%d score:%.*f", hloop->w_time, hloop->tau_counter, hopt->cancel_counter, precision, hopt->cost);
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			g_mutex_lock( &(hdeepinfo->hevalctrl->m) );
			kounter = hdeepinfo->hevalctrl->kounter;
			g_mutex_unlock( &(hdeepinfo->hevalctrl->m) );
			fprintf(fp, " fmean:%.*f fmax:%.*f", precision, hdeepinfo->population->fmean, precision, hdeepinfo->population->fmax);
			fprintf(fp, " L:%d", hdeepinfo->last_method);
			fprintf(fp, " kounter:%d", kounter);
            DpIndivid* individ = hdeepinfo->population->individ[hdeepinfo->population->imin];
            fprintf(fp, " cost:%.*f", precision, individ->cost);
            if (hopt->logdepth > 1) {
				for( i = 0; i < individ->ntargets; i++) {
					fprintf(fp, " target[%i]:%.*f", i, precision, individ->targets[i]);
				}
			}
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
		break;
	}
	if (hopt->logdepth > 0) {
		for ( i = 0; i < heval->size; i++) {
			fprintf(fp, " p[%d]:%.*f", heval->points[i]->index, precision, *(heval->points[i]->param));
		}
	}
	fprintf(fp, "\n");
	fclose(fp);
	return ret_val;
}

DpLoopExitCode dp_print_log(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	DpTarget*htarget = (DpTarget*)(hopt->htarget);
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	char*msg = g_path_get_basename(hopt->logname);
	FILE*fp;
	int i, kounter;
	int precision = hopt->precision;
	fprintf(stdout, "%s wtime:%e tau:%d freeze:%d score:%.*f", msg, hloop->w_time, hloop->tau_counter, hopt->cancel_counter, precision, hopt->cost);
		g_free(msg);
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			g_mutex_lock( &(hdeepinfo->hevalctrl->m) );
			kounter = hdeepinfo->hevalctrl->kounter;
			g_mutex_unlock( &(hdeepinfo->hevalctrl->m) );
			fprintf(stdout, " fmean:%.*f fmax:%.*f", precision, hdeepinfo->population->fmean, precision, hdeepinfo->population->fmax);
			fprintf(stdout, " L:%d", hdeepinfo->last_method);
			fprintf(stdout, " kounter:%d", kounter);
            DpIndivid* individ = hdeepinfo->population->individ[hdeepinfo->population->imin];
            fprintf(stdout, " cost:%.*f", precision, individ->cost);
            if (hopt->logdepth > 1) {
				for( i = 0; i < individ->ntargets; i++) {
					fprintf(stdout, " target[%i]:%.*f", i, precision, individ->targets[i]);
				}
			}
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
		break;
	}
	if (hopt->logdepth > 0) {
		for ( i = 0; i < heval->size; i++) {
			fprintf(stdout, " p[%d]:%.*f", heval->points[i]->index, precision, *(heval->points[i]->param));
		}
	}
	fprintf(stdout, "\n");
	return ret_val;
}

DpLoopExitCode dp_read_log(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	DpTarget*htarget = (DpTarget*)(hopt->htarget);
	DpDeepInfo*hdeepinfo;
	DpPopulation*pop;
	DpEvaluationCtrl*hevalctrl;
	hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	double fmean, fmax, cost, *targets;
	int last_method, kounter;
	char *base;
	gchar **base_tokens;
	int n_base_tokens, n_extra_tokens, len_log, len_par;
	base = (char *)calloc(4 * MAX_RECORD, sizeof(char *));
	int monitor = hopt->monitor;
	FILE*fp;
	int i;
	fp = fopen(hopt->logname, "r");
	if ( !fp && hopt->monitor > 0 ) {
		hloop->exit_str = g_strdup_printf ( "dp_read_log: can't open %s", hopt->logname);
		return DP_LOOP_EXIT_ERROR;
	} else if ( !fp && hopt->monitor < 1 ) {
		return ret_val;
	}
    base = fgets(base, 16 * MAX_RECORD, fp);
    base_tokens = g_strsplit_set(base, ": ", -1);
	len_log = strlen(base_tokens[0]);
    n_base_tokens = g_strv_length(base_tokens) - 1;
	n_extra_tokens = 0;
	len_par = 0;
	for ( i = 0; i < n_base_tokens; i++) {
		if (g_str_has_prefix (base_tokens[i], "target")) {
			n_extra_tokens++;
		}
		if (g_str_has_prefix (base_tokens[i], "p[")) {
			len_par++;
		}
	}
	if (heval->size != len_par) {
		g_error("Number of parameters %d != number of p %d", heval->size, len_par);
	}
    g_strfreev(base_tokens);
    free(base);
    targets = g_new(double, n_extra_tokens);
	do {
		for ( i = 0; i < len_log; i++) {
			fscanf(fp, "%*c");
		}
		fscanf(fp, " wtime:%lf tau:%d freeze:%d score:%lf", &(hloop->w_time), &(hloop->tau_counter), &(hopt->stop_counter), &(hopt->cost));
		fscanf(fp, " fmean:%lf fmax:%lf", &fmean, &fmax);
		fscanf(fp, " L:%d", &last_method);
		fscanf(fp, " kounter:%d", &kounter);
		fscanf(fp, " cost:%lf", &cost);
/*		if (hopt->logdepth > 1) {*/
			for ( i = 0; i < n_extra_tokens; i++) {
				fscanf(fp, " target[%*i]:%lf", &targets[i]);
			}
/*		}*/
/*		if (hopt->logdepth > 0) {*/
			for ( i = 0; i < heval->size; i++) {
				fscanf(fp, " p[%d]:%lf", &(heval->points[i]->index), heval->points[i]->param);
			}
/*		}*/
	} while ( !feof(fp) && monitor != hloop->tau_counter );
	fclose(fp);
	if ( hopt->monitor < 1 ) {
		pop = hdeepinfo->population;
		hevalctrl = hdeepinfo->hevalctrl;
		dp_evaluation_individ_set(hevalctrl, pop->individ[pop->imin]);
		for ( i = 0; i < n_extra_tokens; i++) {
			pop->individ[pop->imin]->targets[i] = targets[i];
		}
		pop->individ[pop->imin]->cost = cost;
	}
	g_free(targets);
	return ret_val;
}

DpLoopExitCode dp_write_state(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	FILE*fp;
	fp = fopen(hopt->chkname, "w");
	if ( !fp ) {
		hloop->exit_str = g_strdup_printf ( "dp_write_state: can't open %s", hopt->chkname);
		return DP_LOOP_EXIT_ERROR;
	}
	fprintf(fp, "%e %e %e %e %d %e %d %e\n", hloop->w_time, hloop->u_time, hloop->start, hloop->finish, hloop->tau_counter, hopt->cost, hopt->stop_counter, hopt->old_cost);
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_deep_info_save(fp, hdeepinfo);
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
			dp_osda_info_save(fp, hosdainfo);
		break;
	}
	fclose(fp);
	return ret_val;
}

DpLoopExitCode dp_rotate_target(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpTarget*htarget = (DpTarget*)(hopt->htarget);
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	int *indices_of_targets;
	int n_indices;
	int offset;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_target_shift_penalty_ranks (htarget);
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_read_state(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	FILE*fp;
	if ( g_file_test (hopt->chkname, G_FILE_TEST_EXISTS) ) {
	  fp = fopen(hopt->chkname, "r");
	  if ( !fp ) {
	    hloop->exit_str = g_strdup_printf ( "dp_read_state: can't open %s", hopt->chkname);
	    return DP_LOOP_EXIT_ERROR;
	  }
	  fscanf(fp, "%e %e %e %e %d %e %d %e", &(hloop->w_time), &(hloop->u_time), &(hloop->start), &(hloop->finish), &(hloop->tau_counter), &(hopt->cost), &(hopt->stop_counter), &(hopt->old_cost));
	  if (hopt->advance_kounter == 1) {
		hopt->cancel_iter += hloop->tau_counter;
	  }
	  switch (hopt->opt_type) {
	  case H_OPT_DEEP:
	    hdeepinfo = (DpDeepInfo*)(hopt->method_info);
	    dp_deep_info_load(fp, hdeepinfo);
	    break;
	  case H_OPT_OSDA:
	    hosdainfo = (DpOsdaInfo*)(hopt->method_info);
	    dp_osda_info_load(fp, hosdainfo);
	    break;
	  }
	  fclose(fp);
	}
	/*	dp_loop_del_func(hloop, dp_read_state);*/
	return ret_val;
}

DpLoopExitCode dp_write_tst(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_SUCCESS;
	DpOpt*hopt = (DpOpt*)user_data;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	DpTarget*htarget = (DpTarget*)(hopt->htarget);
	FILE*fp;
	int i;
	fp = fopen(hopt->tstname, "w");
	if ( !fp ) {
		hloop->exit_str = g_strdup_printf ( "dp_write_tst: can't open %s", hopt->tstname);
		return DP_LOOP_EXIT_ERROR;
	}
	fprintf(fp, "user       time:%e\n", hloop->u_time);
	fprintf(fp, "wall clock time:%e\n", hloop->w_time);
	fprintf(fp, "tau     counter:%d\n", hloop->tau_counter);
	fprintf(fp, "final      cost:%e\n", hopt->cost);
//	fprintf(fp, "cost       runs:%d\n", hopt->htarget->target->kount);
	for ( i = 0; i < hopt->htarget->array_size; i++ ) {
//		fprintf(fp, "penalty[%d]cost:%e\n", i, hopt->htarget->penalty[i]->retval);
//		fprintf(fp, "penalty[%d]runs:%d\n", i, hopt->htarget->penalty[i]->kount);
	}
	fprintf(fp, "\n");
	fclose(fp);
	return ret_val;
}

DpLoopExitCode dp_opt_osda(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpOsdaInfo*hosdainfo = (DpOsdaInfo*)(hopt->method_info);
	double energy_start;
	energy_start = hopt->cost;
	dp_osda_step(hosdainfo);
	dp_osda_accept_step(hosdainfo, &(hopt->cost));
	dp_osda_update_step(hosdainfo);
	return ret_val;
}

DpLoopExitCode dp_opt_mpi_distribute(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	int stop_flag = ( hloop->stop_flag == DP_LOOP_EXIT_SUCCESS ) ? 1 : 0;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_population_mpi_distribute(hdeepinfo->trial, hopt->world_id, hopt->world_count);
		break;
	}
	hloop->stop_flag = ( stop_flag == 1 ) ? DP_LOOP_EXIT_SUCCESS : DP_LOOP_EXIT_NOEXIT;
	return ret_val;
}

DpLoopExitCode dp_opt_mpi_gather(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	int stop_flag = ( hloop->stop_flag == DP_LOOP_EXIT_SUCCESS ) ? 1 : 0;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_population_mpi_gather(hdeepinfo->trial, hopt->world_id, hopt->world_count);
		break;
	}
	hloop->stop_flag = ( stop_flag == 1 ) ? DP_LOOP_EXIT_SUCCESS : DP_LOOP_EXIT_NOEXIT;
	return ret_val;
}

/*
 * Copy the best several times
 *
 */

DpLoopExitCode dp_opt_substitute(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	int stop_flag = ( hloop->stop_flag == DP_LOOP_EXIT_SUCCESS ) ? 1 : 0;
	int src_ind, dst_ind, i, k, substieps, substeps;
	GMainContext *gcontext = g_main_context_default();
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			DpEvaluationCtrl*hevalctrl = hdeepinfo->hevalctrl;
			GError*gerror = NULL;
			gboolean immediate_stop = TRUE;
			gboolean wait_finish = FALSE;
			gulong microseconds = G_USEC_PER_SEC / 1000;
			int notdone, currstatus;
			hevalctrl->gthreadpool = g_thread_pool_new ((GFunc) dp_evaluation_population_init_func, (gpointer) hevalctrl, hevalctrl->eval_max_threads, hevalctrl->exclusive, &gerror);
			if ( gerror != NULL ) {
				g_error("%s", gerror->message);
			}
			if (hdeepinfo->es_kind == 0) {
				src_ind = hdeepinfo->population->cost_ascending[0];
			}
            for (i = 0; i < hdeepinfo->es_lambda; i++) {
				if (hdeepinfo->es_kind == 1) {
					src_ind = hdeepinfo->population->cost_ascending[i];
				}
                dst_ind = hdeepinfo->population->ages_descending[i];
				g_debug ("dp_opt_substitute: trying count=%d; iter=%d; src=%d; age=%d; failures=%d; cost=%f; dst=%d; age=%d; failures=%d; cost=%f;", i, hdeepinfo->population->iter, src_ind, hdeepinfo->population->individ[src_ind]->age, hdeepinfo->population->individ[src_ind]->failures, hdeepinfo->population->individ[src_ind]->cost, dst_ind, hdeepinfo->population->individ[dst_ind]->age, hdeepinfo->population->individ[dst_ind]->failures, hdeepinfo->population->individ[dst_ind]->cost);
                if (hdeepinfo->population->individ[dst_ind]->age < hdeepinfo->es_cutoff) break;
                if ( dst_ind != hdeepinfo->population->imin && dst_ind != src_ind ) {
					g_debug ("dp_opt_substitute: changed count=%d; iter=%d; src=%d; age=%d; failures=%d; cost=%f; dst=%d; age=%d; failures=%d; cost=%f;", i, hdeepinfo->population->iter, src_ind, hdeepinfo->population->individ[src_ind]->age, hdeepinfo->population->individ[src_ind]->failures, hdeepinfo->population->individ[src_ind]->cost, dst_ind, hdeepinfo->population->individ[dst_ind]->age, hdeepinfo->population->individ[dst_ind]->failures, hdeepinfo->population->individ[dst_ind]->cost);
                    dp_individ_copy_values(hdeepinfo->population->individ[dst_ind], hdeepinfo->population->individ[src_ind]);
                    substeps = hdeepinfo->substeps;
                    if (hopt->delay > 0 && hopt->stop_counter >= hopt->delay) {
						substeps = 0;
					}
                    if (substeps >= 0) {
						dp_evaluation_individ_scramble(hdeepinfo->hevalctrl, hdeepinfo->population->individ[dst_ind], substeps);
						g_thread_pool_push (hevalctrl->gthreadpool, (gpointer)(hdeepinfo->population->individ[dst_ind]), &gerror);
						if ( gerror != NULL ) {
							g_error("%s", gerror->message);
						}
						hdeepinfo->population->individ[dst_ind]->age = 0;
						hdeepinfo->population->individ[dst_ind]->moves = 0;
						hdeepinfo->population->individ[dst_ind]->failures = 0;
						hdeepinfo->population->individ[dst_ind]->grads = 0;
					}
                }
			}
			src_ind = hdeepinfo->population->cost_ascending[0];
            for (i = 0; i < hdeepinfo->population->size; i++) {
                dst_ind = hdeepinfo->population->ages_descending[i];
				g_debug ("dp_opt_substitute: trying count=%d; iter=%d; src=%d; age=%d; failures=%d; cost=%f; dst=%d; age=%d; failures=%d; cost=%f;", i, hdeepinfo->population->iter, src_ind, hdeepinfo->population->individ[src_ind]->age, hdeepinfo->population->individ[src_ind]->failures, hdeepinfo->population->individ[src_ind]->cost, dst_ind, hdeepinfo->population->individ[dst_ind]->age, hdeepinfo->population->individ[dst_ind]->failures, hdeepinfo->population->individ[dst_ind]->cost);
                if (hdeepinfo->population->individ[dst_ind]->invalid == 1) {
	                if ( dst_ind != hdeepinfo->population->imin && dst_ind != src_ind ) {
						g_debug ("dp_opt_substitute: changed count=%d; iter=%d; src=%d; age=%d; failures=%d; cost=%f; dst=%d; age=%d; failures=%d; cost=%f;", i, hdeepinfo->population->iter, src_ind, hdeepinfo->population->individ[src_ind]->age, hdeepinfo->population->individ[src_ind]->failures, hdeepinfo->population->individ[src_ind]->cost, dst_ind, hdeepinfo->population->individ[dst_ind]->age, hdeepinfo->population->individ[dst_ind]->failures, hdeepinfo->population->individ[dst_ind]->cost);
	                    dp_individ_copy_values(hdeepinfo->population->individ[dst_ind], hdeepinfo->population->individ[src_ind]);
	                    substieps = hdeepinfo->substieps;
	                    if (hopt->delay > 0 && hopt->stop_counter >= hopt->delay) {
							substieps = 0;
						}
	                    if (substieps >= 0) {
							dp_evaluation_individ_scramble(hdeepinfo->hevalctrl, hdeepinfo->population->individ[dst_ind], substieps);
							g_thread_pool_push (hevalctrl->gthreadpool, (gpointer)(hdeepinfo->population->individ[dst_ind]), &gerror);
							if ( gerror != NULL ) {
								g_error("%s", gerror->message);
							}
							hdeepinfo->population->individ[dst_ind]->age = 0;
							hdeepinfo->population->individ[dst_ind]->moves = 0;
							hdeepinfo->population->individ[dst_ind]->failures = 0;
							hdeepinfo->population->individ[dst_ind]->grads = 0;
						}
	                }
				}
			}
			notdone = 1;
			while(notdone == 1 || g_thread_pool_unprocessed (hevalctrl->gthreadpool) > 0) {
				g_main_context_iteration(gcontext, FALSE);
	            g_usleep (microseconds);
				notdone = 0;
				for (i = 0; i < hdeepinfo->population->size; i++) {
					g_mutex_lock( &(hdeepinfo->population->individ[i]->m) );
					currstatus = hdeepinfo->population->individ[i]->status;
					g_mutex_unlock( &(hdeepinfo->population->individ[i]->m) );
					if ( currstatus == 1 ) {
						notdone = 1;
					}
				}
	        }
			g_thread_pool_free (hevalctrl->gthreadpool, immediate_stop, wait_finish);
			dp_population_update(hdeepinfo->population, 0, hdeepinfo->population->size);
		break;
	}
	hloop->stop_flag = ( stop_flag == 1 ) ? DP_LOOP_EXIT_SUCCESS : DP_LOOP_EXIT_NOEXIT;
	return ret_val;
}

DpLoopExitCode dp_opt_permute(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
    int i, src_ind, dst_ind, seed_ind;
    DpIndivid*tmp_ind;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			DpPopulation*pop = hdeepinfo->population;
			GRand*hrand = hdeepinfo->hevalctrl->hrand;
			seed_ind = hdeepinfo->es_lambda % pop->size;
			for ( i = 0; i < hdeepinfo->es_lambda; i++ ) {
                src_ind = g_rand_int_range (hrand, 0, pop->size);
                dst_ind = g_rand_int_range (hrand, 0, pop->size);
                tmp_ind = pop->individ[dst_ind];
                pop->individ[dst_ind] = pop->individ[src_ind];
                pop->individ[src_ind] = tmp_ind;
			}
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_post(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_deep_post(hdeepinfo);
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
			dp_osda_post(hosdainfo);
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_post_evaluate(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_deep_post_evaluate(hdeepinfo);
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
			dp_osda_post_evaluate(hosdainfo);
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_opt_evaluate_pareto_front(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpPopulation*pop;
	GArray*fronts;
	GArray*curr_front, *next_front, *tmp_front;
	int i, curr_front_index, curr_ind, j, flag, k;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			pop = hdeepinfo->popunion;
			for ( i = 0; i < pop->size; i++ ) {
                pop->individ[i]->pareto_front = -1;
                pop->individ[i]->dom_count = 0;
                if ( pop->individ[i]->dominated != NULL ) {
                    g_array_free (pop->individ[i]->dominated, TRUE);
                }
                pop->individ[i]->dominated = g_array_new (FALSE, FALSE, sizeof (int));
			}
			fronts = g_array_new (FALSE, FALSE, sizeof (GArray*));
			curr_front = g_array_new (FALSE, FALSE, sizeof (int));
            for ( i = 0; i < pop->size; i++ ) {
                for ( j = 0; j < pop->size; j++ ) {
                    if ( 1 == dp_evaluation_individ_dominates_compare((const void *)(&(pop->individ[i])), (const void *)(&(pop->individ[j])), (void*)(hdeepinfo->hevalctrl)) ) {
                        g_array_append_val (pop->individ[i]->dominated, j);
                    } else if ( 1 == dp_evaluation_individ_dominates_compare((const void *)(&(pop->individ[j])), (const void *)(&(pop->individ[i])), (void*)(hdeepinfo->hevalctrl)) ) {
                        pop->individ[i]->dom_count++;
                    }
                }
                if ( pop->individ[i]->dom_count == 0 ) {
                    pop->individ[i]->pareto_front = 0;
                    g_array_append_val (curr_front, i);
                }
            }
            g_array_append_val (fronts, curr_front);
            curr_front_index = 0;
            do {
                next_front = g_array_new (FALSE, FALSE, sizeof (int));
                flag = 0;
                for ( k = 0; k < curr_front->len; k++ ) {
                    i = g_array_index (curr_front, int, k);
                    for ( j = 0; j < pop->individ[i]->dominated->len; j++ ) {
                        curr_ind = g_array_index (pop->individ[i]->dominated, int, j);
                        pop->individ[curr_ind]->dom_count--;
                        if ( pop->individ[curr_ind]->dom_count == 0 ) {
                            pop->individ[curr_ind]->pareto_front = curr_front_index + 1;
                            g_array_append_val (next_front, curr_ind);
                            flag = 1;
                        }
                    }
                }
                if ( flag == 1 ) {
                    curr_front_index++;
                    g_array_append_val (fronts, next_front);
                } else {
                    g_array_free ( next_front, FALSE );
                }
                tmp_front = curr_front;
                curr_front = next_front;
            } while ( flag == 1 );
            pop->nfronts = curr_front_index;
            if ( pop->fronts != NULL ) {
                for ( k = 0; k < pop->fronts->len; k++ ) {
                    tmp_front = g_array_index (pop->fronts, GArray*, k);
                    g_array_free ( tmp_front, TRUE );
                }
                g_array_free ( pop->fronts, TRUE );
            }
            pop->fronts = fronts;
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_sort_pareto_front(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpPopulation*pop, *trial, *population;
	GArray*curr_front;
	int i, curr_ind, j, k;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			pop = hdeepinfo->popunion;
            for ( k = 0; k < pop->fronts->len; k++ ) {
                curr_front = g_array_index (pop->fronts, GArray*, k);
                g_array_sort_with_data(curr_front, (GCompareDataFunc) dp_evaluation_cr_compare, (gpointer) pop);
            }
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_select_pareto_front(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpPopulation*pop, *trial, *population;
	GArray*curr_front;
	int i, curr_ind, j, k;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
            if ( hdeepinfo->selector == DpSelectorSelected ) {
                return ret_val;
            }
			pop = hdeepinfo->popunion;
			population = hdeepinfo->population;
			trial = hdeepinfo->trial;
			i = 0;
            for ( k = 0; k < pop->fronts->len; k++ ) {
                curr_front = g_array_index (pop->fronts, GArray*, k);
                for ( j = 0; j < curr_front->len; j++ ) {
                    curr_ind = g_array_index (curr_front, int, j);
                    if ( i < population->size ) {
                        population->individ[i] = pop->individ[curr_ind];
                        population->individ[i]->age = 0;
                    } else {
                        trial->individ[i - population->size] = pop->individ[curr_ind];
                        trial->individ[i - population->size]->age++;
                    }
                    i++;
                }
            }
            trial->iter = population->iter;
            population->iter++;
            dp_population_update(population, 0, population->size);
            dp_deep_accept_step(hdeepinfo, &(hopt->cost));
            hdeepinfo->selector = DpSelectorSelected;
		break;
	}
	return ret_val;
}

DpLoopExitCode dp_write_pareto(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpPopulation*pop;
	GArray*curr_front;
	GString *params;
	int i, curr_ind, j, k, lim_ind;
	DpEvaluation*heval = (DpEvaluation*)(hopt->heval);
	FILE*fp;
	gchar*pareto_logname;
	pareto_logname = g_strdup_printf("%s_pareto", hopt->logname);
	fp = fopen(pareto_logname, "a");
	if ( !fp ) {
		hloop->exit_str = g_strdup_printf ( "dp_write_pareto: can't open %s", pareto_logname);
		return DP_LOOP_EXIT_ERROR;
	}
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
		    hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			DpTarget*htarget = hopt->htarget;
			pop = hdeepinfo->popunion;
			lim_ind = ( hopt->pareto_all < 0 ) ? pop->fronts->len : hopt->pareto_all + 1;
            for ( k = 0; k < lim_ind; k++ ) {
                curr_front = g_array_index (pop->fronts, GArray*, k);
                for ( j = 0; j < curr_front->len; j++ ) {
                    curr_ind = g_array_index (curr_front, int, j);
                    DpIndivid* individ = pop->individ[curr_ind];
                    params = htarget->params_to_string(individ->user_data, individ->z);
                    fprintf(fp, "wtime:%e tau:%d cost:%f front:%d:%d", hloop->w_time, hloop->tau_counter, individ->cost, k, j);
                    for( i = 0; i < individ->ntargets; i++) {
                        fprintf(fp, " target[%i]:%f", i, individ->targets[i]);
                    }
/*                    for ( i = 0; i < heval->size; i++) {
                        fprintf(fp, " p[%d]:%f", heval->points[i]->index, individ->z[i]);
                    }*/
                    fprintf(fp, " %s", params->str);
                    fprintf(fp, "\n");
                    g_string_free(params, TRUE);
                }
                fprintf(fp, "\n\n");
            }
            fprintf(fp, "\n\n\n\n");
		break;
	}
	fclose(fp);
	g_free(pareto_logname);
	return ret_val;
}

DpLoopExitCode dp_opt_cr2cost(DpLoop*hloop, gpointer user_data)
{
	DpLoopExitCode ret_val = DP_LOOP_EXIT_NOEXIT;
	DpOpt*hopt = (DpOpt*)user_data;
	DpDeepInfo*hdeepinfo;
	DpOsdaInfo*hosdainfo;
	switch (hopt->opt_type) {
		case H_OPT_DEEP:
			hdeepinfo = (DpDeepInfo*)(hopt->method_info);
			dp_population_cr2cost(hdeepinfo->popunion);
		break;
		case H_OPT_OSDA:
			hosdainfo = (DpOsdaInfo*)(hopt->method_info);
		break;
	}
	return ret_val;
}
