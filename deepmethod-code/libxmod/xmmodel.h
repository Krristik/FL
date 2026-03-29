/***************************************************************************
 *            xmmodel.h
 *
 *  Tue Mar 27 12:06:21 2012
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


// Here you can find model-related data
// COnverts extern param into proccessable model

#ifndef _XM_MODEL_H
#define _XM_MODEL_H

#ifdef __cplusplus
extern "C"
{
#endif
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

typedef gchar* (*XmParmsConvert)(gpointer *user_data, GError **err);

typedef enum XmModelType {
	XmModelCommand = (1 << 0),
	XmModelIntprt  = (1 << 1),
	XmModelNone = (1 << 2)
} XmModelType;

typedef struct XmModelConn {
	int source;
	int source_conn;
	int dest;
	int dest_conn;
	int rank;
} XmModelConn;

typedef struct XmModelPart {
	gchar*name;
	int num_parms;
	int *index;
	int checked_flag;
	int checked_num_conn;
} XmModelPart;

typedef struct XmModel { //
	XmModelType type;
	int debug;
	int size;
	int *parms;
	int *iparms;
	int *mask;
	int *param_type; /* 0 - double, 1 - index, 2 - round, 3 - fixed */
	int *has_params;
	int *subsubset;
	int **type_index;
	int **index_type;
	int *limited;
	int *index;
	int *tweak_index;
	int index_size;
	int *tweak;
    int *dedupl;
	double *dparms;
	double *bparms;
	double *lbound;
	double *hbound;
	double *scale;
	XmModelPart *part;
	int *lookup;
	int num_parts;
	gchar*command;
	gchar*prime_command;
	gchar*convert;
	XmParmsConvert converter;
	gchar*delimiters;
	int num_keys; /* length of keys, array */
	int *keys;
	int num_values; /* length of mapping */
	double *array;
	int *mapping;
	gchar*prime_delimiters;
	int *prime_keys;
	int num_prime_values;
	double *prime_array;
	int *prime_mapping;
	int current_penalty_index;
	double current_functional_value;
	double functional_value;
	int prime_index;
	int a_precision;
	int b_precision;
	int ref_counter;
	int copy_val_parms;
	int copy_counter;
	// for synch parallel exec xmlmodel
	GMutex m;
	GAsyncQueue *queue_intprts;
	int num_threads;
	int timeoutsec;
} XmModel;

XmModel*xm_model_new();

int xm_model_init(gchar*filename, gchar*groupname, XmModel*xmmodel, GError **err);

double xm_model_score(gpointer user_data);

double xm_model_parms_double_to_int(gpointer user_data, double*x);

double xm_model_dparms_to_int(gpointer user_data);

double xm_model_dparms_to_index(gpointer user_data);

double xm_model_parms_double_to_index(gpointer user_data, double*x);

double xm_model_score_double(gpointer user_data, double*x);

double xm_model_objfunc(gpointer user_data, double*x);

double xm_model_score_int(gpointer user_data, double*x);

double xm_model_read_penalty(gpointer user_data, double*x);

double xm_model_barrier_penalty(gpointer user_data, double*x);

gchar *xm_model_convert_parms_to_sdf(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_gcdm(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_gemstat(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_dgemstat(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_subset(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_subsubset(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_ini(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_octave(gpointer *user_data, GError **err);

gchar *xm_model_convert_parms_to_r(gpointer *user_data, GError **err);

XmModelConn*xm_model_conn_new(int source, int source_conn, int dest, int dest_conn, int rank);

int xm_model_conn_comp (gconstpointer a, gconstpointer b);

int xm_model_conn_comp_ports (gconstpointer a, gconstpointer b);

void xm_model_conn_free(XmModelConn*conn);

GString*xm_model_sdf_contents(XmModel*xmmodel);

GString*xm_model_gcdm_contents(XmModel*xmmodel);

GString*xm_model_gemstat_contents(XmModel*xmmodel);

GString*xm_model_dgemstat_contents(XmModel*xmmodel);

GString*xm_model_subset_contents(XmModel*xmmodel);

GString*xm_model_subsubset_contents(XmModel*xmmodel);

GString*xm_model_ini_contents(XmModel*xmmodel);

GString*xm_model_octave_contents(XmModel*xmmodel);

GString*xm_model_r_contents(XmModel*xmmodel);

void xm_model_save(XmModel*xmmodel, gchar*filename);

int xm_model_get_next_free_conn(int connected_to_conn, XmModel*xmmodel);

gint xm_model_sort_atoms (gconstpointer a, gconstpointer b);

gpointer xm_model_copy_values(gpointer psrc);

void xm_model_set_dparms(XmModel *xmmodel, double*x);

void xm_model_update_values(gpointer psrc, double*x, int prime_index, double functional_value);

int xm_model_run(XmModel *xmmodel);

int xm_model_run_prime(XmModel *xmmodel);

double xm_model_prime_double(gpointer user_data, double*x);

double xm_model_prime_int(gpointer user_data, double*x);

GString *xm_model_score_int_to_string(gpointer user_data, double*x);

GString *xm_model_score_double_to_string(gpointer user_data, double*x);

void xm_model_convert_param_type(XmModel *xmmodel);

#ifdef __cplusplus
}
#endif

#endif /* _XM_MODEL_H */
