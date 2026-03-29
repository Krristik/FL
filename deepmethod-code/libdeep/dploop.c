/***************************************************************************
 *            dploop.c
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
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#if defined(__MINGW32__) || defined(_MSC_VER)

#else
#include <sys/times.h>
#endif
#if defined(_MSC_VER)

#else
#include <unistd.h>
#endif
#include <time.h>
#include "dploop.h"

DpLoopRunFunc *dp_loop_run_func_new(int tau_flag, DpLoopFunc loop_func, gpointer user_data)
{
	DpLoopRunFunc *func = (DpLoopRunFunc *)malloc(sizeof(DpLoopRunFunc));
	func->tau_flag = tau_flag;
	func->user_data = user_data;
	func->func = loop_func;
	return func;
}

DpLoop *dp_loop_new(GList *before_funcs, GList *funcs, GList *after_funcs)
{
	DpLoop *hloop = (DpLoop*)malloc(sizeof(DpLoop));
	hloop->funcs = funcs;
	hloop->run_once_before = before_funcs;
	hloop->run_once_after = after_funcs;
#if defined(__MINGW32__) || defined(_MSC_VER)

#else
	hloop->clk_tck = (double)sysconf(_SC_CLK_TCK);
	hloop->cpu_start  = (struct tms *)malloc(sizeof(struct tms));
	hloop->cpu_finish = (struct tms *)malloc(sizeof(struct tms));
#endif
	return hloop;
}

void dp_loop_add_func(DpLoop*hloop, DpLoopRunFunc *func)
{
	hloop->funcs = g_list_append(hloop->funcs, (gpointer) func);
}

void dp_loop_del_func(DpLoop*hloop, DpLoopRunFunc *func)
{
	hloop->funcs = g_list_remove_all(hloop->funcs, (gconstpointer) func);
}

void dp_loop_add_func_after(DpLoop*hloop, DpLoopRunFunc *func)
{
	hloop->run_once_after = g_list_append(hloop->run_once_after, (gpointer) func);
}

void dp_loop_add_func_before(DpLoop*hloop, DpLoopRunFunc *func)
{
	hloop->run_once_before = g_list_append(hloop->run_once_before, (gpointer) func);
}

void dp_loop_run_func(void *arg, gpointer user_data)
{
	DpLoop*hloop = (DpLoop*)user_data;
	DpLoopRunFunc*func = (DpLoopRunFunc*)arg;
	if ( ( ( ( hloop->tau_counter ) % ( func->tau_flag ) ) == 0 ) && ( hloop->stop_flag == DP_LOOP_EXIT_NOEXIT ) ) {
		hloop->stop_flag = func->func (hloop, func->user_data);
	}
}

void dp_loop_run_once_func(void *arg, gpointer user_data)
{
	DpLoop*hloop = (DpLoop*)user_data;
	DpLoopRunFunc*func = (DpLoopRunFunc*)arg;
	if ( hloop->stop_flag != DP_LOOP_EXIT_ERROR ) {
		hloop->stop_flag = func->func (hloop, func->user_data);
	}
}

void dp_loop_zero_counters(DpLoop*hloop)
{
#if defined(__MINGW32__) || defined(_MSC_VER)

#else
	times(hloop->cpu_start);
	hloop->start = time(NULL);
#endif
	hloop->stop_flag = DP_LOOP_EXIT_NOEXIT;
	hloop->tau_counter = 0;
	hloop->w_time = 0;
	hloop->u_time = 0;
}

void dp_loop_run(DpLoop*hloop)
{
	g_list_foreach(hloop->run_once_before, dp_loop_run_once_func, (gpointer)hloop);
	while ( hloop->stop_flag == DP_LOOP_EXIT_NOEXIT ) {
		hloop->tau_counter++;
		g_list_foreach(hloop->funcs, dp_loop_run_func, (gpointer)hloop);
#if defined(__MINGW32__) || defined(_MSC_VER)

#else
		times(hloop->cpu_finish);
		hloop->finish = time(NULL);
		hloop->w_time += hloop->finish - hloop->start;
		hloop->u_time += (hloop->cpu_finish->tms_utime - hloop->cpu_start->tms_utime) / ( hloop->clk_tck );
		hloop->start = hloop->finish;
		hloop->cpu_start->tms_utime = hloop->cpu_finish->tms_utime;
#endif
	}
	if ( hloop->stop_flag == DP_LOOP_EXIT_SUCCESS ) {
		g_list_foreach(hloop->run_once_after, dp_loop_run_once_func, (gpointer)hloop);
	}
}

