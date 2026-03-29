/***************************************************************************
 *            dploop.h
 *
 *  Fri Mar 23 16:25:49 2012
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
 
// Implementation of eternal cycle
// Not related wth mathematical part of DEEP

#ifndef _DP_LOOP_H
#define _DP_LOOP_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

typedef enum DpLoopExitCode {
	DP_LOOP_EXIT_SUCCESS,
	DP_LOOP_EXIT_ERROR,
	DP_LOOP_EXIT_NOEXIT
} DpLoopExitCode;

typedef struct DpLoop {
	struct tms *cpu_start;
	struct tms *cpu_finish;
	double start;
	double finish;
	double w_time;
	double u_time;
	int tau_counter;
	DpLoopExitCode stop_flag;
	gchar *exit_str;
	GList *run_once_before;
	GList *funcs;
	GList *run_once_after;
	double clk_tck; /* system-specific clock tick duration */
} DpLoop;

typedef DpLoopExitCode (*DpLoopFunc)(DpLoop*hloop, void *user_data);

typedef struct DpLoopRunFunc {
	DpLoopFunc func;
	gpointer user_data;
	int tau_flag;
} DpLoopRunFunc;

void dp_loop_run_func(void *arg, gpointer user_data);

void dp_loop_run_once_func(void *arg, gpointer user_data);

DpLoopRunFunc *dp_loop_run_func_new(int tau_flag, DpLoopFunc loop_func, gpointer user_data);

DpLoop *dp_loop_new(GList *before_funcs, GList *funcs, GList *after_funcs);

void dp_loop_zero_counters(DpLoop*hloop);

void dp_loop_run(DpLoop*hloop);

void dp_loop_del_func(DpLoop*hloop, DpLoopRunFunc *func);

void dp_loop_add_func_before(DpLoop*hloop, DpLoopRunFunc *func);

void dp_loop_add_func(DpLoop*hloop, DpLoopRunFunc *func);

void dp_loop_add_func_after(DpLoop*hloop, DpLoopRunFunc *func);

#ifdef __cplusplus
}
#endif

#endif /* _DP_LOOP_H */
