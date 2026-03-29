/***************************************************************************
 *            dpopt.h
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

#ifndef _DP_OPT_H
#define _DP_OPT_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include "dpevaluation.h"
#include "dprecombination.h"
#include "dploop.h"
#include "dptarget.h"

#define MAX_RECORD 4096

typedef enum DpOptType {
	H_OPT_FIRST,
	H_OPT_NONE,
	H_OPT_DEEP,
	H_OPT_OSDA,
	H_OPT_LAST
} DpOptType;

typedef enum DpOptStopType {
    H_OPT_PROPORTIONAL_STOP, //
    H_OPT_ABSOLUTE_STOP, //
    H_OPT_ABSOLUTE_SCORE, //
    H_OPT_ABSOLUTE_TIME, //
    H_OPT_ABSOLUTE_ITER //
} DpOptStopType;

typedef enum DpCancel {
	DP_PROPORTIONAL_CANCEL = (1 << 0),
	DP_ABSOLUTE_CANCEL = (1 << 1),
	DP_ABSOLUTE_SCORE = (1 << 2),
	DP_ABSOLUTE_TIME = (1 << 3),
	DP_ABSOLUTE_ITER = (1 << 4)
} DpCancel;

typedef struct DpOpt {
	DpOptType opt_type;
	DpEvaluation*heval;
	DpTarget*htarget;
	void *method_info;
	double cost;
	double new_cost;
	double old_cost;
	double cost_start;
	int delay;
	int logdepth;
/* deprecated */
	DpOptStopType stop_type;
	int stop_count;
	double criterion;
	int stop_counter;
/* end deprecated */
	guint cancel_flags;
	int cancel_count;
	int cancel_iter;
	double cancel_time;
	double cancel_prop;
	double cancel_abs;
	double cancel_score;
	int cancel_counter;
/* end new api */
	char*filename;
	char*logname;
	char*chkname;
	char*tstname;
	DpLoop*hloop;
	int world_id;
	int world_count;
	int monitor;
	int pareto_all;
	int precision;
	int debug;
	int keep_order;
	int advance_kounter;
} DpOpt;

//DpOptSettings*dp_opt_settings_new();

DpOpt *dp_opt_init(DpEvaluation*heval, DpTarget*htarget, int world_id, int world_count,char*filename, DpOptStopType stop_type, double criterion, int stop_count, int pareto_all, int precision);

void dp_opt_add_func(DpOpt *hopt, DpLoopFunc func, int tau_flag, DpOptType opt_type, int order, gpointer method_info);

void dp_opt_add_from_func_list(gchar**list, DpOpt *hopt, int order, GKeyFile*gkf, gchar*groupname, int world_id, DpEvaluation*heval, DpTarget*htarget, GError**err);

void dp_opt_run(DpOpt *hopt);

void dp_opt_monitor(DpOpt *hopt, int monitor, GError**gerror);

DpLoopExitCode dp_opt_popsize(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_write_log(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_print_log(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_write_state(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_read_state(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_write_tst(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_check_stop(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_check_cancel(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_generate(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_generate_dd(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_generate_ca(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_generate_ac(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_evaluate(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_select(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_de_select(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_osda(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_mpi_distribute(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_mpi_gather(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_post(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_post_evaluate(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_init_stop(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_init_cancel(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_read_log(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_rotate_target(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_write_pareto(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_evaluate_pareto_front(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_cr2cost(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_select_pareto_front(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_sort_pareto_front(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_permute(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_substitute(DpLoop*hloop, gpointer user_data);

DpLoopExitCode dp_opt_deep_update(DpLoop*hloop, gpointer user_data);

#ifdef __cplusplus
}
#endif

#endif /* _DP_OPT_H */
