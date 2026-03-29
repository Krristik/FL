/***************************************************************************
 *            xmmodel.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <glib/gprintf.h>
#include <glib/gstdio.h>
#ifdef GIO_STANDALONE_SOURCE
#include <gio/gfile.h>
#else
#include <gio/gio.h>
#endif
#include "xmmodel.h"

#ifdef GIO_STANDALONE_SOURCE
#define g_strcmp0(str1, str2) strcmp(str1, str2)
#endif

int check_debug (gchar*debug_str) {
	int verbose = 0;
	const gchar*p = g_getenv("G_MESSAGES_DEBUG");
	if (p == NULL) {
		return verbose;
	}
	gchar*debug_env = g_strdup(p);
	if (debug_env == NULL) {
		return verbose;
	}
	gchar**debug_envv = g_strsplit(debug_env, " ", -1);
	if (debug_envv == NULL) {
		g_free(debug_env);
		return verbose;
	}
#ifndef OLDGLIB
	verbose = g_strv_contains ((const gchar* const*)debug_envv, debug_str);
#endif
	g_free(debug_env);
	g_strfreev(debug_envv);
	return verbose;
}

XmModel*xm_model_new()
{
	XmModel *xmmodel = (XmModel*)g_malloc( sizeof(XmModel) );
	xmmodel->debug = check_debug(G_LOG_DOMAIN);
	xmmodel->size = 0;
	xmmodel->tweak = NULL;
	xmmodel->command = NULL;
	xmmodel->delimiters = NULL;
	xmmodel->keys = NULL;
	xmmodel->num_values = 0;
	xmmodel->num_keys = 0;
	xmmodel->array = NULL;
	xmmodel->mapping = NULL;
	xmmodel->prime_command = NULL;
	xmmodel->prime_delimiters = NULL;
	xmmodel->prime_keys = NULL;
	xmmodel->num_prime_values = 0;
	xmmodel->prime_array = NULL;
	xmmodel->prime_mapping = NULL;
	xmmodel->current_penalty_index = 0;
	xmmodel->current_functional_value = G_MAXDOUBLE;
	xmmodel->prime_index = 0;
	xmmodel->functional_value = G_MAXDOUBLE;
	xmmodel->parms = NULL;
	xmmodel->limited = NULL;
	xmmodel->index = NULL;
	xmmodel->dparms = NULL;
	xmmodel->bparms = NULL;
	xmmodel->lbound = NULL;
	xmmodel->hbound = NULL;
	xmmodel->scale = NULL;
	xmmodel->convert = NULL;
	xmmodel->a_precision = 16;
	xmmodel->b_precision = 9;
	xmmodel->ref_counter = 0;
	xmmodel->copy_val_parms = 0;
	xmmodel->copy_counter = 0;
	xmmodel->type = XmModelCommand;
	xmmodel->timeoutsec = 600;
	xmmodel->param_type = NULL;
	xmmodel->has_params = NULL;
	xmmodel->mask = NULL;
    xmmodel->dedupl = NULL;
	return xmmodel;
}

gpointer xm_model_copy_values(gpointer psrc)
{
	int j, k;
	XmModel *src = (XmModel *)psrc;
	XmModel *xmmodel = xm_model_new();
	xmmodel->index_size = src->index_size;
	xmmodel->index = g_new0 ( int, xmmodel->index_size );
	xmmodel->tweak_index = g_new0 ( int, xmmodel->index_size );
	for ( j = 0; j < xmmodel->index_size; j++ ) {
		xmmodel->index[j] = src->index[j];
		xmmodel->tweak_index[j] = src->tweak_index[j];
	}
	xmmodel->num_parts = src->num_parts;
	xmmodel->part = g_new ( XmModelPart, xmmodel->num_parts );
	for ( j = 0; j < xmmodel->num_parts; j++ ) {
		xmmodel->part[j].name = g_strdup(src->part[j].name);
		xmmodel->part[j].num_parms = src->part[j].num_parms;
		xmmodel->part[j].index = g_new( int, xmmodel->part[j].num_parms );
		for ( k = 0; k < xmmodel->part[j].num_parms; k++ ) {
			xmmodel->part[j].index[k] = src->part[j].index[k];
		}
	}
	xmmodel->converter = src->converter;
	xmmodel->size = src->size;
	xmmodel->debug = src->debug;
	xmmodel->parms = g_new0 ( int, xmmodel->size );
	xmmodel->iparms = g_new0 ( int, xmmodel->size );
	xmmodel->lookup = g_new0 ( int, xmmodel->size );
	xmmodel->tweak = g_new0 ( int, xmmodel->size );
	xmmodel->limited = g_new0 ( int, xmmodel->size );
	xmmodel->mask = g_new0 ( int, xmmodel->size );
	xmmodel->param_type = g_new0 ( int, xmmodel->size );
    xmmodel->dedupl = g_new0 ( int, xmmodel->size );
    xmmodel->dparms = g_new0 ( double, xmmodel->size );
	xmmodel->bparms = g_new0 ( double, xmmodel->size );
	xmmodel->lbound = g_new0 ( double, xmmodel->size );
	xmmodel->hbound = g_new0 ( double, xmmodel->size );
	xmmodel->scale = g_new0 ( double, xmmodel->size );
	for ( j = 0; j < xmmodel->size; j ++ ) {
		xmmodel->parms[j] = src->parms[j];
		xmmodel->iparms[j] = src->iparms[j];
		xmmodel->lookup[j] = src->lookup[j];
		xmmodel->mask[j] = src->mask[j];
		xmmodel->tweak[j] = src->tweak[j];
		xmmodel->limited[j] = src->limited[j];
		xmmodel->dparms[j] = src->dparms[j];
		xmmodel->bparms[j] = src->bparms[j];
		xmmodel->lbound[j] = src->lbound[j];
		xmmodel->hbound[j] = src->hbound[j];
		xmmodel->scale[j] = src->scale[j];
		xmmodel->param_type[j] = src->param_type[j];
        xmmodel->dedupl[j] = src->dedupl[j];
	}
	xmmodel->command = g_strdup(src->command);
	if ( src->prime_command != NULL) {
		xmmodel->prime_command = g_strdup(src->prime_command);
		xmmodel->prime_delimiters = g_strdup(src->prime_delimiters);
		xmmodel->num_prime_values = src->num_prime_values;
		xmmodel->prime_array = g_new0 ( double, xmmodel->num_prime_values );
		xmmodel->prime_keys = g_new0 ( int, xmmodel->num_prime_values );
		xmmodel->prime_mapping = g_new0 ( int, xmmodel->num_prime_values );
		for ( j = 0; j < xmmodel->num_prime_values; j++ ) {
			xmmodel->prime_array[j] = src->prime_array[j];
			xmmodel->prime_keys[j] = src->prime_keys[j];
			xmmodel->prime_mapping[j] = src->prime_mapping[j];
		}
		xmmodel->prime_index = src->prime_index;
	}
	xmmodel->delimiters = g_strdup(src->delimiters);
	xmmodel->num_values = src->num_values;
	xmmodel->num_keys = src->num_keys;
	xmmodel->array = g_new0 ( double, xmmodel->num_keys );
	xmmodel->keys = g_new0 ( int, xmmodel->num_keys );
	xmmodel->mapping = g_new0 ( int, xmmodel->num_values );
	for ( j = 0; j < xmmodel->num_keys; j++ ) {
		xmmodel->array[j] = src->array[j];
		xmmodel->keys[j] = src->keys[j];
	}
	for ( j = 0; j < xmmodel->num_values; j++ ) {
		xmmodel->mapping[j] = src->mapping[j];
	}
	xmmodel->current_penalty_index = src->current_penalty_index;
	xmmodel->current_functional_value = src->current_functional_value;
	xmmodel->convert = g_strdup(src->convert);
	xmmodel->functional_value = src->functional_value;
	xmmodel->a_precision = src->a_precision;
	xmmodel->b_precision = src->b_precision;
	xmmodel->ref_counter = src->copy_counter + 1;
	xmmodel->copy_val_parms = src->copy_val_parms;
	xmmodel->copy_counter = 0;
	src->copy_counter++;
	xmmodel->queue_intprts = src->queue_intprts;
	xmmodel->type = src->type;
	xmmodel->timeoutsec = src->timeoutsec;
    xmmodel->has_params = g_new0(int, 5);
	xmmodel->has_params[0] = src->has_params[0];
	xmmodel->has_params[1] = src->has_params[1];
	xmmodel->has_params[2] = src->has_params[2];
	xmmodel->has_params[3] = src->has_params[3];
    xmmodel->has_params[4] = src->has_params[4];
    xmmodel->index_type = g_new0(int*, 5);
    xmmodel->type_index = g_new0(int*, 5);
	xmmodel->index_type[0] = g_new0(int, xmmodel->has_params[0]);
	xmmodel->type_index[0] = g_new0(int, xmmodel->has_params[0]);
	xmmodel->index_type[1] = g_new0(int, xmmodel->has_params[1]);
	xmmodel->type_index[1] = g_new0(int, xmmodel->has_params[1]);
	xmmodel->index_type[2] = g_new0(int, xmmodel->has_params[2]);
	xmmodel->type_index[2] = g_new0(int, xmmodel->has_params[2]);
	xmmodel->index_type[3] = g_new0(int, xmmodel->has_params[3]);
	xmmodel->type_index[3] = g_new0(int, xmmodel->has_params[3]);
    xmmodel->index_type[4] = g_new0(int, xmmodel->has_params[4]);
    xmmodel->type_index[4] = g_new0(int, xmmodel->has_params[4]);
    for ( j = 0; j < xmmodel->has_params[0]; j++ ) {
		xmmodel->type_index[0][j] = src->type_index[0][j];
	}
	for ( j = 0; j < xmmodel->has_params[1]; j++ ) {
		xmmodel->type_index[1][j] = src->type_index[1][j];
	}
	for ( j = 0; j < xmmodel->has_params[2]; j++ ) {
		xmmodel->type_index[2][j] = src->type_index[2][j];
	}
	for ( j = 0; j < xmmodel->has_params[3]; j++ ) {
		xmmodel->type_index[3][j] = src->type_index[3][j];
    }
    for ( j = 0; j < xmmodel->has_params[4]; j++ ) {
        xmmodel->type_index[4][j] = src->type_index[4][j];
    }
	return xmmodel;
}

void xm_model_update_values(gpointer psrc, double*x, int prime_index, double functional_value)
{
	XmModel *xmmodel = (XmModel *)psrc;
	int i, index;
	for ( i = 0; i < xmmodel->index_size; i++) {
		index = xmmodel->tweak_index[i];
		xmmodel->dparms[index] = x[i];
	}
	xmmodel->prime_index = prime_index;
	xmmodel->functional_value = functional_value;
}

void xm_model_set_dparms(XmModel *xmmodel, double*x)
{
	int i, index, j;
	for ( i = 0; i < xmmodel->index_size; i++) {
		index = xmmodel->tweak_index[i];
		xmmodel->dparms[index] = x[i];
	}
}

gchar*param2str(XmModel*xmmodel, int i)
{
	gchar*str = NULL;
	if (xmmodel->dparms[i] != xmmodel->dparms[i]) { // is NaN
		g_debug("p[%d] = %.13f is nan", i, xmmodel->dparms[i]);
		return NULL;
	} else {
		if (xmmodel->param_type[i] == 0) {
			str = g_strdup_printf("%.*f", xmmodel->b_precision, xmmodel->dparms[i]);
		} else if (xmmodel->param_type[i] == 1) {
			str = g_strdup_printf("%d", xmmodel->parms[i]);
		} else if (xmmodel->param_type[i] == 2) {
			str = g_strdup_printf("%d", xmmodel->parms[i]);
		} else if (xmmodel->param_type[i] == 3) {
			str = g_strdup_printf("%d", xmmodel->parms[i]);
        } else if (xmmodel->param_type[i] == 4) {
            str = g_strdup_printf("%d", xmmodel->parms[i]);
        }
	}
	return str;
}

gchar*dparam2str(XmModel*xmmodel, int i)
{
	gchar*str = NULL;
	if (xmmodel->dparms[i] != xmmodel->dparms[i]) { // is NaN
		g_debug("p[%d] = %.13f is nan", i, xmmodel->dparms[i]);
		return NULL;
	} else {
		str = g_strdup_printf("%.*f", xmmodel->b_precision, xmmodel->dparms[i]);
	}
	return str;
}

typedef struct API {
	gchar*file_cmd;
	gchar*source_cmd;
	gchar*fflush_cmd;
	gchar*func_cmd;
} API;

// interpreter with io channels
typedef struct Interpreter{
	GIOChannel *in,
               *out,
               *err;
	GMutex m;
	GCond cond;
	gchar *response;
	GPid child_pid;
	int debug;
	GMutex mstatus;
	int status;
	API *api;
} Interpreter;

GString *interpreter_recieve_response(GIOChannel*out, int debug, API*api)
{
	GString *response = g_string_new("");
	GIOStatus status = G_IO_STATUS_NORMAL;
    gchar *response_line = NULL, *suffix;
    gboolean is_done = FALSE;
    gsize size;
    GError *gerror = NULL;
	int max_reponse_size = G_MAXINT16;
    int j = 0;
    suffix = g_strdup_printf("%s\n", api->fflush_cmd);
    do {
        status = g_io_channel_read_line( out, &response_line, &size, NULL, &gerror );
		if (status != G_IO_STATUS_NORMAL) {
			g_debug("Status is not NORMAL %d", status);
			if (gerror != NULL) {
				g_print("%s", gerror->message);
				g_error_free(gerror);
			}
			break;
		}
        response = g_string_append(response, response_line);
        is_done = g_str_has_suffix(response->str, suffix);
		g_debug("line [%d] '%s' last=%d '%s' %d\n", j, response_line, is_done, response->str, status);
        j++;
		if (j > max_reponse_size) is_done = TRUE;
        g_free(response_line);
    } while(!is_done);
    g_free(suffix);
    return response;
}

static gboolean out_watch( GIOChannel   *channel,
        GIOCondition  cond, gpointer data)
{
    GString *response;
    gsize size;
    int currstatus, wrongstart;
    Interpreter *intprt = (Interpreter*)data;
	wrongstart = 0;
    if( cond == G_IO_HUP ){
        g_debug("Cond %d is G_IO_HUP %d", intprt->child_pid, (int)cond);
/*        g_io_channel_unref( channel );*/
        return( FALSE );
    } else if( cond == G_IO_IN ) {
		g_debug("Cond %d is G_IO_IN %d", intprt->child_pid, (int)cond);
    } else if( cond == G_IO_OUT ) {
		g_debug("Cond %d is G_IO_OUT %d", intprt->child_pid, (int)cond);
		wrongstart = 1;
    }
/* status:
 * 0 - interpreter waits for command
 * 1 - command was sent to interpreter
 * 2 - we want to read the command output
 */
	g_mutex_lock( &(intprt->mstatus) );
	currstatus = intprt->status;
	g_mutex_unlock( &(intprt->mstatus) );
	if (currstatus == 0 || currstatus == 2) {
		wrongstart = 1;
	}
	if (wrongstart == 1) {
		g_debug("Cond %d is G_IO_IN %d but status %d", intprt->child_pid, (int)cond, currstatus);
/*		g_string_free(response, TRUE);*/
		return TRUE;
	}
/* status was 1 - let's change it to 2 */
	g_mutex_lock( &(intprt->mstatus) );
	intprt->status = 2;
	g_mutex_unlock( &(intprt->mstatus) );
	g_debug("Interpreter out [%d]:", intprt->child_pid);
    response = interpreter_recieve_response(channel, intprt->debug, intprt->api);
/* read info - "free" interpreter */
    g_mutex_lock( &(intprt->mstatus) );
	intprt->status = 0;
	g_mutex_unlock( &(intprt->mstatus) );
    g_mutex_lock( &(intprt->m) );
    intprt->response = g_string_free(response, FALSE);
    g_cond_signal( &(intprt->cond) );
    g_mutex_unlock( &(intprt->m) );
    return TRUE;
}

static gboolean err_watch( GIOChannel   *channel,
        GIOCondition  cond, gpointer data)
{
    gchar *string;
    gsize  size;

    if( cond == G_IO_HUP ) {
/*        g_io_channel_unref( channel );*/
		g_debug("Err G_IO_HUP");
        return (TRUE);
    }

    g_io_channel_read_line( channel, &string, &size, NULL, NULL );
    g_print("Err: %s", string);
    g_free( string );

    return TRUE;
}

void write_to_interpreter(GIOChannel* in, GString* msg){
    GError *err = NULL;
    GIOStatus gs = g_io_channel_write_chars( in,
            msg->str,
            msg->len,
            NULL,
            &err);
    if (gs != G_IO_STATUS_NORMAL) {
	g_debug("Write chars status bad %d", gs);
	if (err != NULL) {
		g_debug("%s", err->message);
	}
	return;
    }
    g_io_channel_flush(in, NULL);
}

void init_source_to_interpreter(GIOChannel * intprt_in, API*api){
	if (strlen(api->source_cmd) > 0) {
		GString * command = g_string_new(api->source_cmd);
		g_string_append_printf(command, "(\'%s\')\r\n%s\r\n", api->file_cmd, api->fflush_cmd);
//		g_string_append_printf(command, "(\'%s\');\r\n", api->file_cmd);
		write_to_interpreter(intprt_in, command);
	}
}

Interpreter* init_interpreter(XmModel * xmmodel){
/*    gchar      *argvr[] = { "R", "--no-save", "--silent", "--vanilla", NULL };
 * { "R", "--no-save", "--silent", "--vanilla", "script.R", "source", "func", "flush.console()", NULL } */
	gchar      **args;
    gint        in,
                out,
                err;
    gboolean    ret;
	GPid child_pid;
	GError*gerror = NULL;
	int length;
	if ( !g_shell_parse_argv(xmmodel->command, &length, &args, &gerror) ) {
		if ( gerror ) {
			g_error("g_shell_parse failed for %s\nwith %s", xmmodel->command, gerror->message);
		}
	}
/* remove and parse extra args */
    API* api = g_new(API, 1);
    api->file_cmd = g_strdup(args[length - 4]);
    g_free(args[length - 4]);
    args[length - 4] = NULL;
    api->source_cmd = g_strdup(args[length - 3]);
    g_free(args[length - 3]);
    args[length - 3] = NULL;
    api->func_cmd = g_strdup(args[length - 2]);
    g_free(args[length - 2]);
    args[length - 2] = NULL;
    api->fflush_cmd = g_strdup(args[length - 1]);
    g_free(args[length - 1]);
    args[length - 1] = NULL;
    /* Spawn child process */
    ret = g_spawn_async_with_pipes( NULL, args, NULL,
            G_SPAWN_SEARCH_PATH | G_SPAWN_DO_NOT_REAP_CHILD, NULL,
            NULL, &child_pid,
            &in, &out, &err, NULL );
    if ( !ret ) {
        g_printf("g_spawn_async_with_pipes in init interpreter\nfailed for %s\nwith %s", args[0], gerror->message);
        return NULL;
    }
    Interpreter* intprt = g_new(Interpreter, 1);
    intprt->in = g_io_channel_unix_new( in );
    intprt->out = g_io_channel_unix_new( out );
    intprt->err = g_io_channel_unix_new( err );
	intprt->child_pid = child_pid;
	intprt->status = 0;
	g_mutex_init( &(intprt->mstatus) );
    g_mutex_init( &(intprt->m) );
    g_cond_init( &(intprt->cond) );
    intprt->response = NULL;
	intprt->debug = xmmodel->debug;
	intprt->api = api;
    /* Add watches to channels */
    g_io_add_watch( intprt->out, G_IO_IN | G_IO_OUT | G_IO_HUP, (GIOFunc)out_watch, intprt);
    g_io_add_watch( intprt->err, G_IO_IN | G_IO_OUT | G_IO_HUP, (GIOFunc)err_watch, intprt);
	if (strlen(api->source_cmd) > 0) {
		init_source_to_interpreter(intprt->in, intprt->api);
		g_debug("Interpreter init [%d]:", intprt->child_pid);
		GString *response = interpreter_recieve_response(intprt->out, intprt->debug, intprt->api);
		g_string_free(response, TRUE);
	}
    return intprt;
}

void api_free(API*api)
{
	g_free(api->file_cmd);
	g_free(api->source_cmd);
	g_free(api->func_cmd);
	g_free(api->fflush_cmd);
	g_free(api);
}

void kill_interpreter(Interpreter* intprt)
{
	if (g_source_remove_by_user_data ((gpointer) intprt)) {
		g_debug("Source 1 removed [%d]:", intprt->child_pid);
	} else {
		g_debug("Source not 1 removed [%d]:", intprt->child_pid);
	}
	if (g_source_remove_by_user_data ((gpointer) intprt)) {
		g_debug("Source 2 removed [%d]:", intprt->child_pid);
	} else {
		g_debug("Source not 2 removed [%d]:", intprt->child_pid);
	}
	g_io_channel_unref (intprt->out);
	g_io_channel_unref (intprt->err);
	g_mutex_clear (&(intprt->m));
	g_cond_clear (&(intprt->cond));
	api_free(intprt->api);
#ifndef G_OS_WIN32
	kill(intprt->child_pid, SIGKILL);
#endif
	g_spawn_close_pid (intprt->child_pid);
	g_debug("Interp killed [%d]:", intprt->child_pid);
	g_free(intprt);
}

GString *create_command(XmModel*xmmodel, API*api)
{
	GString *command = g_string_new(api->func_cmd);
	gchar*buf;
	int i = 0;
	command = g_string_append_c(command, '(');
	for(i = 0; i < xmmodel->size - 1; i++) {
		if ( (buf = param2str(xmmodel, i)) == NULL) {
			g_string_free(command, TRUE);
			return NULL;
		}
		command = g_string_append(command, buf);
		g_free(buf);
		command = g_string_append_c(command, ',');
	}
	i = xmmodel->size - 1;
	if ( (buf = param2str(xmmodel, i)) == NULL) {
		g_string_free(command, TRUE);
		return NULL;
	}
	command = g_string_append(command, buf);
	g_free(buf);
	g_string_append_printf(command, ");\r\n%s\r\n", api->fflush_cmd);
//	g_string_append_printf(command, ");\r\n");
	return command;
}

int xm_model_run_interpreter(XmModel *xmmodel)
{
	gchar**result;
	gchar**margv;
	int argcp, i, j;
	GString *command;
	GError*gerror = NULL;
	// посмотреть откуда взять статус
	int flaggs, child_exit_status = 1;
	gchar *standard_output;
	gchar*conversion;
	double max_value = G_MAXDOUBLE;
	gint64 end_time;
	int nsec = xmmodel->timeoutsec;/* wait for command completion that number of seconds*/
	gboolean failed = FALSE;
	GAsyncQueue * queue = (GAsyncQueue*)xmmodel->queue_intprts;
	end_time = g_get_monotonic_time () + nsec * G_TIME_SPAN_SECOND;
	Interpreter * intprt = NULL;
	intprt = (Interpreter*)g_async_queue_timeout_pop(queue, end_time);
	if (intprt == NULL) {
		g_debug ( "Couldn't get interpreter" );
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			xmmodel->array[i] = max_value;
		}
		xmmodel->copy_val_parms = 0;
		return child_exit_status;
	}
	// create command
	command = create_command(xmmodel, intprt->api);
	if (command == NULL) {
		g_debug ( "Couldn't create command with NaN" );
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			xmmodel->array[i] = max_value;
		}
		xmmodel->copy_val_parms = 0;
		// return interpreter to queue
		g_async_queue_push(queue, (gpointer)intprt);
		return child_exit_status;
	}
	g_mutex_lock( &(intprt->mstatus) );
	intprt->status = 1;
	g_mutex_unlock( &(intprt->mstatus) );
	g_debug("Run starting [%d]:", intprt->child_pid);
    write_to_interpreter(intprt->in, command);
	failed = FALSE;
    g_mutex_lock( &(intprt->m) );
    while( !(intprt->response) && !failed ){
		if (!g_cond_wait_until (&(intprt->cond), &(intprt->m), end_time)) {
			      // timeout has passed.
//			g_mutex_unlock (&data_mutex);
//			return NULL;
			g_debug("Run timeouted [%d]:", intprt->child_pid);
			failed = TRUE;
		}
    }
//	if (!(intprt->response) || failed) {
	if (intprt->response == NULL || (intprt->response != NULL && strlen(intprt->response) < 2) || failed) {
		g_mutex_unlock (&(intprt->m));
		g_debug ( "Interpreter failed [%d]:", intprt->child_pid);
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			xmmodel->array[i] = max_value;
		}
		xmmodel->copy_val_parms = 0;
		kill_interpreter(intprt);
		intprt = init_interpreter(xmmodel);
		if (!intprt) {
			g_debug("Couldn't re-add interpreter!");
			return child_exit_status;
		}
		g_async_queue_push(queue, (gpointer)intprt);
		g_debug("Re-added intprt [%d]:", intprt->child_pid);
		return child_exit_status;
	}
    standard_output = g_strdup( intprt->response );
    g_free( intprt->response );
    intprt->response = NULL;
    g_mutex_unlock (&(intprt->m));
	g_mutex_lock( &(intprt->mstatus) );
	intprt->status = 0;
	g_mutex_unlock( &(intprt->mstatus) );
	result = g_strsplit_set(standard_output, xmmodel->delimiters, -1);
	int result_length = g_strv_length(result);
	int parsing_failed = 0;
	int need_intprt_restart = 0;
	if ( strlen(standard_output) > 0 && result != NULL ) {
		if ( xmmodel->debug == 1 ) {
			g_printf("result_length = %d;\n", result_length);
			for ( j = 0; j < result_length; j++ ) {
				g_printf("result[%d] = %s;\n", j, result[j]);
			}
		}
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			int ik = (xmmodel->keys[i] >= 0) ? xmmodel->keys[i] : result_length + xmmodel->keys[i];
			if (ik < 0 || ik > result_length - 1) {
				g_debug ( "Couldn't parse interpreter [%d] result: %s", intprt->child_pid, standard_output);
				parsing_failed = 1;
				need_intprt_restart = 1;
			} else if ( result[ik] != NULL ) {
				xmmodel->array[i] = g_strtod(result[ik], NULL);
			} else {
				g_debug ( "result[%d] doesn't exist", xmmodel->keys[i]);
			}
		}
	} else {
		parsing_failed = 1;
	}
	if (need_intprt_restart == 1) {
		kill_interpreter(intprt);
		intprt = init_interpreter(xmmodel);
		if (intprt) {
			g_debug("Re-added new intprt [%d]:", intprt->child_pid);
		}
	}
	if (parsing_failed == 1) {
		g_debug ( "Couldn't parse interpreter [%d] output: %s", intprt->child_pid, standard_output);
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			xmmodel->array[i] = max_value;
		}
	}
	g_strfreev(result);
	g_free(standard_output);
    xmmodel->copy_val_parms = 0;

	// return interpreter to queue
	if (intprt) {
		g_async_queue_push(queue, (gpointer)intprt);
	}
	child_exit_status = 0;
	return child_exit_status;
}

int xm_model_run_command(XmModel *xmmodel)
{
	gchar**result;
	gchar**margv;
	int argcp, i, j;
	GString *command;
	GError*gerror = NULL;
	int flaggs, child_exit_status = 1;
	gchar *standard_output;
	gchar *standard_error;
	gchar*conversion;
	double max_value = G_MAXDOUBLE;
	command = g_string_new(xmmodel->command);
	if ( xmmodel->convert != NULL ) {
		if ( ( conversion = xmmodel->converter((gpointer)xmmodel, &gerror) ) == NULL ) {
			g_string_free(command, TRUE);
			for ( i = 0; i < xmmodel->num_values; i++ ) {
				xmmodel->array[i] = max_value;
			}
			xmmodel->copy_val_parms = 0;
			return child_exit_status;
		}
		if ( gerror ) {
			g_error("converter failed for %s", gerror->message);
		}
		gchar*buf;
		buf = g_shell_quote(conversion);
		command = g_string_append_c(command, ' ');
		command = g_string_append(command, buf);
		g_free(buf);
	} else {
        GString *params = g_string_new("");
        gchar*buf;
        for ( i = 0; i < xmmodel->size; i++ ) {
			if ( (buf = param2str(xmmodel, i)) == NULL) {
				g_string_free(params, TRUE);
				for ( i = 0; i < xmmodel->num_values; i++ ) {
					xmmodel->array[i] = max_value;
				}
				xmmodel->copy_val_parms = 0;
				return child_exit_status;
			}
			params = g_string_append(params, buf);
			g_free(buf);
			params = g_string_append_c(params, ' ');
		}
		g_string_append_printf(command, " %s", params->str);
        g_string_free(params, TRUE);
	}
	if ( !g_shell_parse_argv(command->str, &argcp, &margv, &gerror) ) {
		if ( gerror ) {
			g_error("g_shell_parse failed for %s\nwith %s", command->str, gerror->message);
		}
	}
	g_string_free(command, TRUE);
	int g_spawn_failed = 0;
	standard_output = NULL;
	standard_error = NULL;
/*	flaggs = G_SPAWN_SEARCH_PATH | G_SPAWN_STDOUT_TO_DEV_NULL; */
/*	flaggs |= G_SPAWN_STDERR_TO_DEV_NULL; */
	flaggs = G_SPAWN_SEARCH_PATH;
	if ( !g_spawn_sync (NULL, margv, NULL, (GSpawnFlags)flaggs, NULL, NULL, &standard_output, &standard_error, &child_exit_status, &gerror) ) {
		if ( gerror ) {
			g_debug("g_spawn_sync failed for %s\nwith %s", margv[0], gerror->message);
			g_spawn_failed = 1;
		}
	}
	g_strfreev(margv);
	result = (standard_output == NULL) ? NULL:g_strsplit_set(standard_output, xmmodel->delimiters, -1);
	int parsing_failed = 0;
	if ( g_spawn_failed == 0 && standard_output != NULL && strlen(standard_output) > 0 && result != NULL ) {
		int result_length = g_strv_length(result);
		if ( xmmodel->debug == 1 ) {
            g_printf("result_length = %d;\n", result_length);
			for ( j = 0; j < result_length; j++ ) {
                g_printf("result[%d] = %s;\n", j, result[j]);
			}
		}
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			int ik = (xmmodel->keys[i] >= 0) ? xmmodel->keys[i] : result_length + xmmodel->keys[i];
			if (ik < 0 || ik > result_length - 1) {
				g_debug ( "Couldn't parse command result: %s", standard_output);
				parsing_failed = 1;
			} else if (result[ik] == NULL) {
				g_debug ( "result[%d] doesn't exist", xmmodel->keys[i]);
				parsing_failed = 1;
			} else if (!g_ascii_isdigit(result[ik][0]) && result[ik][0] != '-' && result[ik][0] != '.' && result[ik][0] != '+') {
				g_debug ( "result[%d] starts with not a digit, point, plus or minus: %c", xmmodel->keys[i], result[ik][0]);
				parsing_failed = 1;
			} else if ((result[ik][0] == '-' || result[ik][0] == '.' || result[ik][0] == '+') && !g_ascii_isdigit(result[ik][1])) {
				g_debug ( "result[%d] starts with point, plus or minus and no second digit: %c, %c", xmmodel->keys[i], result[ik][0], result[ik][1]);
				parsing_failed = 1;
			} else {
				xmmodel->array[i] = g_strtod(result[ik], NULL);
			}
		}
	} else {
		parsing_failed = 1;
	}
	if (parsing_failed == 1) {
		g_debug ( "Couldn't parse output: %s", standard_output);
		g_debug ( "Standard error: %s", standard_error);
		for ( i = 0; i < xmmodel->num_keys; i++ ) {
			xmmodel->array[i] = max_value;
		}
	}
	if (result != NULL) g_strfreev(result);
	if (standard_output != NULL) g_free(standard_output);
	if (standard_error != NULL) g_free(standard_error);
	if ( xmmodel->convert != NULL ) {
		g_unlink(conversion);
	}
    xmmodel->copy_val_parms = 0;
	return child_exit_status;
}

int xm_model_run(XmModel *xmmodel) {
	int retval;
	if (xmmodel->type == XmModelCommand) {
		retval = xm_model_run_command(xmmodel);
	} else if (xmmodel->type == XmModelIntprt) {
		retval = xm_model_run_interpreter (xmmodel);
	}
	return (retval);
}

int xm_model_run_prime(XmModel *xmmodel)
{
	gchar**result;
	gchar**margv;
	int argcp, i;
	GString *command;
	GError*gerror = NULL;
	int flaggs, child_exit_status;
	gchar *standard_output;
	gchar *standard_error;
	gchar*conversion;
	command = g_string_new(xmmodel->prime_command);
	if ( xmmodel->convert != NULL ) {
		if ( ( conversion = xmmodel->converter((gpointer)xmmodel, &gerror) ) == NULL ) {
			g_string_free(command, TRUE);
			for ( i = 0; i < xmmodel->num_values; i++ ) {
				xmmodel->array[i] = G_MAXDOUBLE;
			}
			return child_exit_status;
		}
		if ( gerror ) {
			g_error("converter failed for %s", gerror->message);
		}
		g_string_append_printf(command, " %s", conversion);
	} else {
        GString *params = g_string_new("");
        for ( i = 0; i < xmmodel->size; i++ ) {
            g_string_append_printf(params, "%f ", xmmodel->dparms[i]);
        }
		g_string_append_printf(command, " '%s'", params->str);
        g_string_free(params, TRUE);
	}
	g_string_append_printf(command, " %d", xmmodel->prime_index);
	if ( !g_shell_parse_argv(command->str, &argcp, &margv, &gerror) ) {
		if ( gerror ) {
			g_error("g_shell_parse failed for %s\nwith %s", command->str, gerror->message);
		}
	}
	g_string_free(command, TRUE);
/*	flaggs = G_SPAWN_SEARCH_PATH | G_SPAWN_STDOUT_TO_DEV_NULL; */
/*	flaggs |= G_SPAWN_STDERR_TO_DEV_NULL; */
	flaggs = G_SPAWN_SEARCH_PATH;
	if ( !g_spawn_sync (NULL, margv, NULL, (GSpawnFlags)flaggs, NULL, NULL, &standard_output, &standard_error, &child_exit_status, &gerror) ) {
		if ( gerror ) {
			g_error("g_spawn_sync failed for %s\nwith %s", margv[0], gerror->message);
		}
	}
	g_strfreev(margv);
	result = g_strsplit_set(standard_output, xmmodel->prime_delimiters, -1);
	for ( i = 0; i < xmmodel->num_prime_values; i++ ) {
		xmmodel->prime_array[i] = g_strtod(result[xmmodel->prime_keys[i]], NULL);
	}
	g_strfreev(result);
	g_free(standard_output);
	g_free(standard_error);
	if ( xmmodel->convert != NULL ) {
		g_unlink(conversion);
	}
	return child_exit_status;
}

int xm_model_parms_sort_index_ascending(gpointer p1, gpointer p2, gpointer user_data)
{
	XmModel *xmmodel = (XmModel *)user_data;
	int*i1 = (int*)p1;
	int*i2 = (int*)p2;
	int ii1 = xmmodel->tweak_index[(*i1)];
	int ii2 = xmmodel->tweak_index[(*i2)];
	int a = xmmodel->dparms[ii1];
	int b = xmmodel->dparms[ii2];
	if ( a < b ) {
		return -1;
	} else if ( a > b ) {
		return 1;
	}
	return 0;
}

double xm_model_dparms_to_int(gpointer user_data)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = -1;
	int i, j;
	if (xmmodel->copy_val_parms > 0 && xmmodel->ref_counter == 1) {
        for ( i = 0; i < xmmodel->size; i++ ) {
			xmmodel->parms[i] = xmmodel->iparms[i];
		}
        return val;
	}
	for ( i = 0; i < xmmodel->index_size; i++ ) {
		xmmodel->index[i] = i;
	}
	g_qsort_with_data(xmmodel->index, xmmodel->index_size, sizeof(xmmodel->index[0]), (GCompareDataFunc)xm_model_parms_sort_index_ascending, xmmodel);
	j = 0;
	for ( i = 0; i < xmmodel->size; i++ ) {
		if ( xmmodel->tweak[i] == 1 ) {
			xmmodel->parms[ xmmodel->tweak_index[ j ] ] = xmmodel->tweak_index[ xmmodel->index[ j ] ];
			j++;
		} else {
			xmmodel->parms[i] = xmmodel->iparms[i];
		}
	}
	return val;
}

double xm_model_parms_double_to_int(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = -1;
	int i, j;
	xm_model_set_dparms(xmmodel, x);
	val = xm_model_dparms_to_int(user_data);
	return val;
}

double xm_model_dparms_to_index(gpointer user_data)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = -1, z, beta;
	double f1 = 500;
	double f2 = 999;
	int alpha;
	int i, j, k;
	j = 0;
	for ( i = 0; i < xmmodel->size; i++ ) {
		if ( xmmodel->tweak[i] == 1 ) {
			k = xmmodel->tweak_index[ j ];
/*			z = f2 * ( 1 + xmmodel->dparms[ k ] ) / f1;*/
			z = xmmodel->dparms[ k ];
			alpha = (int)(z + 0.5);
			beta = alpha - z;
			xmmodel->parms[ k ] = ( beta > 0.5 ) ? alpha - 1 : alpha;
			j++;
		} else {
			xmmodel->parms[i] = xmmodel->iparms[i];
		}
	}
	return val;
}

double xm_model_parms_double_to_index(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = -1;
	int i, j;
	xm_model_set_dparms(xmmodel, x);
	val = xm_model_dparms_to_index(user_data);
	if ( xmmodel->debug == 1 ) {
		for ( i = 0; i < xmmodel->size; i++ ) {
			g_printf("%f;", xmmodel->dparms[i]);
		}
		g_printf("\n");
		for ( i = 0; i < xmmodel->size; i++ ) {
			g_printf("%d;", xmmodel->parms[i]);
		}
	}
	return val;
}

double xm_model_score_double(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	xm_model_set_dparms(xmmodel, x);
    xm_model_run(xmmodel);
    val = xmmodel->array[xmmodel->mapping[0]];
	xmmodel->current_penalty_index = 1;
	xmmodel->current_functional_value = val;
	return val;
}

double xm_model_score_int(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	xm_model_run(xmmodel);
	val = xmmodel->array[xmmodel->mapping[0]];
	xmmodel->current_penalty_index = 1;
	xmmodel->current_functional_value = val;
	return val;
}

int xm_model_type_sort_index_ascending(gpointer p1, gpointer p2, gpointer user_data)
{
	XmModel *xmmodel = (XmModel *)user_data;
	int*i1 = (int*)p1;
	int*i2 = (int*)p2;
	int ii1 = xmmodel->type_index[2][(*i1)];
	int ii2 = xmmodel->type_index[2][(*i2)];
	int a = xmmodel->dparms[ii1];
	int b = xmmodel->dparms[ii2];
	if ( a < b ) {
		return -1;
	} else if ( a > b ) {
		return 1;
	}
	return 0;

}

void xm_model_convert_param_type(XmModel *xmmodel)
{
	double z, beta;
	int alpha;
    int t1, t2, t3, i, j;
/* rounded */
	if (xmmodel->has_params[1] > 0) {
		for ( t1 = 0; t1 < xmmodel->has_params[1]; t1++ ) {
			i = xmmodel->type_index[1][t1];
			z = xmmodel->dparms[i];
			alpha = (int)(z + 0.5);
			beta = alpha - z;
			xmmodel->parms[i] = ( beta > 0.5 ) ? alpha - 1 : alpha;
		}
	}
/*rounded 2 */
    int t4, h = 500, d = 999;
    double x_d;
    if (xmmodel->has_params[4] > 0){
        for ( t4 = 0; t4 < xmmodel->has_params[4]; t4++ ) {
            i = xmmodel->type_index[4][t4];
            x_d = xmmodel->dparms[i];
            xmmodel->parms[i] = (int)(((1 + x_d) * d) / h);
        }
    }
/* index converted */
	if (xmmodel->has_params[2] > 0) {
		for ( t2 = 0; t2 < xmmodel->has_params[2]; t2++ ) {
			xmmodel->index_type[2][t2] = t2;
		}
		g_qsort_with_data(xmmodel->index_type[2], xmmodel->has_params[2], sizeof(xmmodel->index_type[2][0]), (GCompareDataFunc)xm_model_type_sort_index_ascending, xmmodel);
		for ( t2 = 0; t2 < xmmodel->has_params[2]; t2++ ) {
			i = xmmodel->type_index[2][t2];
			j = xmmodel->index_type[2][t2];
			xmmodel->parms[i] = xmmodel->type_index[2][j];
		}
	}

/* fixed */
	if (xmmodel->has_params[3] > 0) {
		for ( t3 = 0; t3 < xmmodel->has_params[3]; t3++ ) {
			i = xmmodel->type_index[3][t3];
			xmmodel->parms[i] = xmmodel->iparms[i];
		}
	}
}

double xm_model_objfunc(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	int ret_val;
	xm_model_set_dparms(xmmodel, x);
	xm_model_convert_param_type(xmmodel);
    ret_val = xm_model_run(xmmodel);
    val = xmmodel->array[xmmodel->mapping[0]];
	xmmodel->current_penalty_index = 1;
	xmmodel->current_functional_value = val;
	if ( xmmodel->debug == 1 ) {
		if (ret_val != 0 && val == G_MAXDOUBLE) {
			g_message("Ret val %d val G_MAXDOUBLE", ret_val);
		} else if (ret_val != 0) {
			g_message("Ret val %d val %f", ret_val, val);
		}
	}
	g_debug("Objfunc: %f", val);
	return val;
}

double xm_model_prime_double(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	xm_model_set_dparms(xmmodel, x);
	xm_model_run_prime(xmmodel);
	val = xmmodel->prime_array[xmmodel->prime_mapping[0]];
	return val;
}

double xm_model_prime_int(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	xm_model_run_prime(xmmodel);
	val = xmmodel->prime_array[xmmodel->prime_mapping[0]];
	return val;
}

double xm_model_read_penalty(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	if ( xmmodel->current_penalty_index > xmmodel->num_values - 1 ) {
        xmmodel->current_penalty_index = 1;
	}
	val = xmmodel->array[xmmodel->mapping[ xmmodel->current_penalty_index ]];
	xmmodel->current_penalty_index++;
	return val;
}

double xm_model_barrier_penalty(gpointer user_data, double*x)
{
	double sqr_dist, dist;
	int i, index;
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	sqr_dist = 0;
	for ( i = 0; i < xmmodel->index_size; i++ ) {
		index = xmmodel->tweak_index[i];
		dist = xmmodel->bparms[index] - x[i];
		sqr_dist += dist * dist;
	}
	val = ( sqr_dist > 0 ) ? ( xmmodel->current_functional_value - xmmodel->functional_value ) / sqr_dist : xmmodel->current_functional_value;
	return val;
}

gchar*xm_model_read(gchar*filename, gsize*size, GError**err)
{
	GFile*gfile;
	char *etag;
	gchar*data = NULL;
	GError *gerror = NULL;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
    gfile = g_file_new_for_commandline_arg(filename);
	if ( !g_file_load_contents(gfile, (GCancellable*)NULL, &data, size, &etag, &gerror) ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
    g_object_unref(gfile);
	return data;
}

int xm_model_load(gchar*data, gsize size, gchar*groupname, XmModel *xmmodel, GError**err)
{
	GKeyFile*gkf;
	GError *gerror = NULL;
	gchar*str, **strlist;
	gsize length;
	int *ilist, ii, j, jj, k, kk, nn, n;
	double *dlist, dval;
	int retval = 0;
	g_return_val_if_fail (err == NULL || *err == NULL, 1);
	gkf = g_key_file_new();
	if ( !g_key_file_load_from_data(gkf, (const gchar*)data, size, G_KEY_FILE_NONE, &gerror ) ) {
		g_propagate_error (err, gerror);
		return 1;
	}
    if ( ( ii = g_key_file_get_integer(gkf, groupname, "a_precision", &gerror) ) != 0  || gerror == NULL ) {
		xmmodel->a_precision = ii;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
    if ( ( ii = g_key_file_get_integer(gkf, groupname, "b_precision", &gerror) ) != 0  || gerror == NULL ) {
		xmmodel->b_precision = ii;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	nn = 0;
	if ( ( strlist = g_key_file_get_string_list(gkf, groupname, "parts", &length, NULL) ) != NULL ) {
		xmmodel->num_parts = length / 2;
		xmmodel->part = g_new ( XmModelPart, xmmodel->num_parts );
		for ( j = 0, jj = 0; j < length; j += 2, jj++ ) {
			xmmodel->part[jj].name = g_strdup(strlist[j]);
			n = g_strtod (strlist[j + 1], NULL);
			nn += n;
			xmmodel->part[jj].num_parms = n;
			xmmodel->part[jj].index = g_new( int, xmmodel->part[jj].num_parms );
		}
		xmmodel->lookup = g_new0 ( int, nn );
		kk = 0;
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				xmmodel->part[jj].index[k] = kk;
				xmmodel->lookup[kk] = jj;
				kk++;
			}
		}
		g_free(strlist);
	} else if ( ( strlist = g_key_file_get_string_list(gkf, groupname, "partnames", &length, &gerror) ) != NULL ) {
		xmmodel->num_parts = length;
		xmmodel->part = g_new ( XmModelPart, xmmodel->num_parts );
		for ( j = 0; j < length; j++ ) {
			xmmodel->part[j].name = g_strdup(strlist[j]);
		}
		if (gerror != NULL) g_clear_error (&gerror);
		if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "partsizes", &length, &gerror) ) == NULL ) {
			g_error("No part sizes!");
		}
		for ( j = 0; j < length; j++ ) {
			n = ilist[j];
			nn += n;
			xmmodel->part[j].num_parms = n;
			xmmodel->part[j].index = g_new( int, xmmodel->part[j].num_parms );
		}
		xmmodel->lookup = g_new0 ( int, nn );
		kk = 0;
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				xmmodel->part[jj].index[k] = kk;
				xmmodel->lookup[kk] = jj;
				kk++;
			}
		}
		g_free(ilist);
		g_free(strlist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if (nn > 0) {
		xmmodel->size = nn;
		xmmodel->parms = g_new0 ( int, xmmodel->size );
		xmmodel->iparms = g_new0 ( int, xmmodel->size );
		xmmodel->param_type = g_new0 ( int, xmmodel->size );
		xmmodel->limited = g_new0 ( int, xmmodel->size );
        xmmodel->dedupl = g_new0 ( int, xmmodel->size );
        xmmodel->dparms = g_new0 ( double, xmmodel->size );
		xmmodel->bparms = g_new0 ( double, xmmodel->size );
		xmmodel->scale = g_new0 ( double, xmmodel->size );
		xmmodel->mask = g_new0 ( int, xmmodel->size );
		xmmodel->tweak = g_new0 ( int, xmmodel->size );
		xmmodel->lbound = g_new0 ( double, xmmodel->size );
		xmmodel->hbound = g_new0 ( double, xmmodel->size );
	} else {
		g_error ("parts absent");
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "parms", &length, NULL) ) != NULL ) {
		if (length != xmmodel->size) {
			g_error ("model size %d != num parms %d", xmmodel->size, length);
		}
		for ( j = 0; j < xmmodel->size; j++ ) {
			xmmodel->parms[j] = ilist[j];
			xmmodel->iparms[j] = xmmodel->parms[j];
			xmmodel->limited[j] = 1;
			xmmodel->param_type[j] = 0;
            xmmodel->dedupl[j] = 0;
        }
		g_free(ilist);
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "valparms", &length, NULL) ) != NULL ) {
		k = 0;
		for ( j = 0; j < xmmodel->size; j++ ) {
			xmmodel->limited[j] = 1;
            xmmodel->dedupl[j] = 0;
            xmmodel->param_type[j] = 0;
			xmmodel->parms[j] = ilist[k];
			xmmodel->iparms[j] = ilist[k];
			xmmodel->dparms[j] = j;
			xmmodel->bparms[j] = xmmodel->dparms[j];
			xmmodel->scale[j] = 1.0;
            k++;
            if ( k >= length ) {
                k = 0;
            }
		}
		xmmodel->copy_val_parms = 1;
		g_free(ilist);
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "parmvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num parms %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->parms[kk] = ilist[jj];
				xmmodel->iparms[kk] = ilist[jj];
			}
		}
		xmmodel->copy_val_parms = 1;
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "partype", &length, NULL) ) != NULL ) {
		for ( j = 0; j < xmmodel->size; j++) {
			xmmodel->param_type[j] = ilist[j];
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "parmtypes", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num types %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->param_type[kk] = ilist[jj];
			}
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
    if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "subsubset", &length, &gerror) ) != NULL ) {
		xmmodel->subsubset = ilist;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "hasparams", &length, &gerror) ) != NULL ) {
		xmmodel->has_params = ilist;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "mask", &length, NULL) ) != NULL ) {
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->mask[j] = ilist[j];
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else if ( ( length = g_key_file_get_integer(gkf, groupname, "nummask", NULL) ) != 0  || gerror == NULL ) {
		for ( j = 0; j < length; j ++ ) {
			xmmodel->mask[j] = 1;
		}
		for ( j = length; j < xmmodel->size; j ++ ) {
			xmmodel->mask[j] = 0;
		}
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "maskvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num mask %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->mask[kk] = ilist[jj];
			}
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "tweak", &length, NULL) ) != NULL ) {
		xmmodel->index_size = 0;
		for ( j = 0; j < xmmodel->size; j++) {
			xmmodel->tweak[j] = ilist[j];
			if ( xmmodel->tweak[j] == 1 ) {
				xmmodel->index_size++;
			}
		}
		xmmodel->index = g_new0 ( int, xmmodel->index_size );
		xmmodel->tweak_index = g_new0 ( int, xmmodel->index_size );
		k = 0;
		for ( j = 0; j < xmmodel->size; j++) {
			if ( xmmodel->tweak[j] == 1 ) {
				xmmodel->tweak_index[k] = j;
				k++;
			}
		}
		g_free(ilist);
	} else if ( ( length = g_key_file_get_integer(gkf, groupname, "numtweak", NULL) ) != 0 ) {
		xmmodel->index_size = length;
		xmmodel->index = g_new0 ( int, xmmodel->index_size );
		xmmodel->tweak_index = g_new0 ( int, xmmodel->index_size );
		for ( j = 0; j < length; j ++ ) {
			xmmodel->tweak[j] = 1;
			xmmodel->tweak_index[j] = j;
		}
		for ( j = length; j < xmmodel->size; j ++ ) {
			xmmodel->tweak[j] = 0;
		}
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "valnotweak", &length, NULL) ) != NULL ) {
		ii = xmmodel->size - length;
		xmmodel->index_size = ii;
		xmmodel->index = g_new0 ( int, ii );
		xmmodel->tweak_index = g_new0 ( int, ii );
		for ( j = 0; j < length; j ++ ) {
			xmmodel->tweak[ ilist[j] ] = -1;
		}
		k = 0;
		for ( j = 0; j < xmmodel->size; j ++ ) {
			if ( xmmodel->tweak[j] == 0 ) {
				xmmodel->tweak[j] = 1;
				xmmodel->tweak_index[k] = j;
				k++;
			}
		}
		for ( j = 0; j < length; j ++ ) {
			xmmodel->tweak[ ilist[j] ] = 0;
		}
		g_free(ilist);
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "tweakvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num tweak %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->tweak[kk] = ilist[jj];
			}
		}
		xmmodel->index_size = 0;
		for ( j = 0; j < xmmodel->size; j++) {
			if ( xmmodel->tweak[j] == 1 ) {
				xmmodel->index_size++;
			}
		}
		xmmodel->index = g_new0 ( int, xmmodel->index_size );
		xmmodel->tweak_index = g_new0 ( int, xmmodel->index_size );
		k = 0;
		for ( j = 0; j < xmmodel->size; j++) {
			if ( xmmodel->tweak[j] == 1 ) {
				xmmodel->tweak_index[k] = j;
				k++;
			}
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "dparms", &length, NULL) ) != NULL ) {
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->dparms[j] = dlist[j];
			xmmodel->bparms[j] = xmmodel->dparms[j];
			xmmodel->scale[j] = 1.0;
		}
		g_free(dlist);
	} else if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "dparmvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num dparm %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->dparms[kk] = dlist[jj];
				xmmodel->bparms[kk] = xmmodel->dparms[kk];
			}
		}
		g_free(dlist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "lbound", &length, NULL) ) != NULL ) {
		if (length != xmmodel->size) {
			g_error ("parts size %d != num lbound %d", xmmodel->size, length);
		}
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->lbound[j] = dlist[j];
		}
		g_free(dlist);
	} else if ( ( dval = g_key_file_get_double(gkf, groupname, "vallbound", NULL) ) != 0 ) {
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->lbound[j] = dval;
		}
	} else if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "lboundvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num lbound %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->lbound[kk] = dlist[jj];
			}
		}
		g_free(dlist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "hbound", &length, NULL) ) != NULL ) {
		if (length != xmmodel->size) {
			g_error ("parts size %d != num hbound %d", xmmodel->size, length);
		}
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->hbound[j] = dlist[j];
		}
		g_free(dlist);
	} else if ( ( dval = g_key_file_get_double(gkf, groupname, "valhbound", NULL) ) != 0 ) {
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->hbound[j] = dval;
		}
	} else if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "hboundvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num hbound %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->hbound[kk] = dlist[jj];
			}
		}
		g_free(dlist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "limited", &length, NULL) ) != NULL ) {
		if (length != xmmodel->size) {
			g_error ("parts size %d != num limit %d", xmmodel->size, length);
		}
		for ( j = 0; j < xmmodel->size; j++) {
			xmmodel->limited[j] = ilist[j];
		}
	} else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "limitvals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num limit %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->limited[kk] = ilist[jj];
			}
		}
		g_free(ilist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "scale", &length, NULL) ) != NULL ) {
		if (length != xmmodel->size) {
			g_error ("parts size %d != num scale %d", xmmodel->size, length);
		}
		for ( j = 0; j < xmmodel->size; j ++ ) {
			xmmodel->scale[j] = dlist[j];
		}
		g_free(dlist);
	} else if ( ( dlist = g_key_file_get_double_list(gkf, groupname, "scalevals", &length, &gerror) ) != NULL ) {
		if (length != xmmodel->num_parts) {
			g_error ("parts number %d != num scale %d", xmmodel->num_parts, length);
		}
		for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
			for ( k = 0; k < xmmodel->part[jj].num_parms; k++ ) {
				kk = xmmodel->part[jj].index[k];
				xmmodel->scale[kk] = dlist[jj];
			}
		}
		g_free(dlist);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "keys", &length, &gerror) ) != NULL ) {
		xmmodel->keys = ilist;
		xmmodel->num_keys = length;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "mapping", &length, &gerror) ) != NULL ) {
		xmmodel->mapping = ilist;
		xmmodel->num_values = length;
		xmmodel->array = g_new0 ( double, length );
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "delimiters", &gerror) ) != NULL ) {
		xmmodel->delimiters = str;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "command", &gerror) ) != NULL ) {
		xmmodel->command = str;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "prime_keys", &length, &gerror) ) != NULL ) {
		xmmodel->prime_keys = ilist;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "prime_mapping", &length, &gerror) ) != NULL ) {
		xmmodel->prime_mapping = ilist;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ii = g_key_file_get_integer(gkf, groupname, "num_prime_values", &gerror) ) != 0  || gerror == NULL ) {
		xmmodel->num_prime_values = ii;
		xmmodel->prime_array = g_new0 ( double, ii );
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ii = g_key_file_get_integer(gkf, groupname, "num_threads", &gerror) ) != 0  || gerror == NULL ) {
		xmmodel->num_threads = ii;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "prime_delimiters", &gerror) ) != NULL ) {
		xmmodel->prime_delimiters = str;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "prime_command", &gerror) ) != NULL ) {
		xmmodel->prime_command = str;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "convert", &gerror) ) != NULL ) {
		xmmodel->convert = str;
		if ( !g_strcmp0 ( xmmodel->convert, "sdf" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_sdf;
		} else if ( !g_strcmp0 ( xmmodel->convert, "gcdm" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_gcdm;
		} else if ( !g_strcmp0 ( xmmodel->convert, "gemstat" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_gemstat;
		} else if ( !g_strcmp0 ( xmmodel->convert, "dgemstat" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_dgemstat;
		} else if ( !g_strcmp0 ( xmmodel->convert, "subset" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_subset;
		} else if ( !g_strcmp0 ( xmmodel->convert, "subsubset" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_subsubset;
		} else if ( !g_strcmp0 ( xmmodel->convert, "octave" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_octave;
		} else if ( !g_strcmp0 ( xmmodel->convert, "ini" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_ini;
		} else if ( !g_strcmp0 ( xmmodel->convert, "r" ) ) {
			xmmodel->converter = xm_model_convert_parms_to_r;
		}
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( str = g_key_file_get_string(gkf, groupname, "type", &gerror) ) != NULL ) {
		if ( !g_strcmp0 ( str, "command" ) ) {
			xmmodel->type = XmModelCommand;
		} else if ( !g_strcmp0 ( str, "interpreter" ) ) {
			xmmodel->type = XmModelIntprt;
		}
		g_free(str);
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}
	if ( ( ii = g_key_file_get_integer(gkf, groupname, "timeoutsec", &gerror) ) != 0  || gerror == NULL ) {
		xmmodel->timeoutsec = ii;
		if (gerror != NULL) g_clear_error (&gerror);
	} else {
		g_debug ("%s", gerror->message );
		g_clear_error (&gerror);
	}

    if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "dedupled", &length, NULL) ) != NULL ) {
        if (length != xmmodel->size) {
            g_error ("parts size %d != num to deduplicate %d", xmmodel->size, length);
        }
        for ( j = 0; j < xmmodel->size; j++) {
            xmmodel->dedupl[j] = ilist[j];
        }
    } else if ( ( ilist = g_key_file_get_integer_list(gkf, groupname, "deduplvals", &length, &gerror) ) != NULL ) {
        if (length != xmmodel->num_parts) { // корректно обработать ошибку
            g_error ("parts number %d != num to deduplicate %d", xmmodel->num_parts, length);
        }

        for ( jj = 0; jj < xmmodel->num_parts; jj++ ) {
            if (ilist[jj] > xmmodel->part[jj].num_parms){
                g_error ("num to deduplicate %d > part size %d", ilist[jj], xmmodel->part[jj].num_parms);
            }
            for ( k = 0; k < ilist[jj]; k++ ) {
                kk = xmmodel->part[jj].index[k];
                xmmodel->dedupl[kk] = 1;
            }

            for ( k = ilist[jj]; k < xmmodel->part[jj].num_parms; k++ ) {
                kk = xmmodel->part[jj].index[k];
                xmmodel->dedupl[kk] = 0;
            }
        }
        g_free(ilist);
        if (gerror != NULL) g_clear_error (&gerror);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }

	g_key_file_free(gkf);
	return retval;
}

int xm_model_init(gchar*filename, gchar*groupname, XmModel*xmmodel, GError **err)
{
	int ret_val = 0;
	gchar*data = NULL;
	gsize size;
	GError *gerror = NULL;
    int t0, t1, t2, t3, t4, j, k;
	g_return_val_if_fail (err == NULL || *err == NULL, 1);
	if ( ( data = xm_model_read(filename, &size, &gerror) ) == NULL ) {
		g_propagate_error (err, gerror);
		return 1;
	}
	ret_val = xm_model_load(data, size, groupname, xmmodel, &gerror);
	if ( ret_val == 1 ) {
		g_propagate_error (err, gerror);
		ret_val = 1;
	}
	g_free(data);
	if (!xmmodel->has_params) {
        xmmodel->has_params = g_new0(int, 5);
		for ( j = 0; j < xmmodel->size; j++ ) {
			if (xmmodel->tweak[j] == 1) {
				xmmodel->has_params[xmmodel->param_type[j]]++;
			} else if (xmmodel->param_type[j] == 1 || xmmodel->param_type[j] == 2 ) {
				xmmodel->has_params[3]++;
				xmmodel->param_type[j] = 3;
			}
		}
        xmmodel->index_type = g_new0(int*, 5);
        xmmodel->type_index = g_new0(int*, 5);
		xmmodel->index_type[0] = g_new0(int, xmmodel->has_params[0]);
		xmmodel->type_index[0] = g_new0(int, xmmodel->has_params[0]);
		xmmodel->index_type[1] = g_new0(int, xmmodel->has_params[1]);
		xmmodel->type_index[1] = g_new0(int, xmmodel->has_params[1]);
		xmmodel->index_type[2] = g_new0(int, xmmodel->has_params[2]);
		xmmodel->type_index[2] = g_new0(int, xmmodel->has_params[2]);
		xmmodel->index_type[3] = g_new0(int, xmmodel->has_params[3]);
		xmmodel->type_index[3] = g_new0(int, xmmodel->has_params[3]);
        xmmodel->index_type[4] = g_new0(int, xmmodel->has_params[4]);
        xmmodel->type_index[4] = g_new0(int, xmmodel->has_params[4]);

		t0 = 0;
		t1 = 0;
		t2 = 0;
		t3 = 0;
        t4 = 0;
		for ( j = 0; j < xmmodel->size; j++ ) {
			if (xmmodel->param_type[j] == 0 && xmmodel->tweak[j] == 1) {
				xmmodel->type_index[0][t0] = j;
				t0++;
			}
			if (xmmodel->param_type[j] == 1) {
				xmmodel->type_index[1][t1] = j;
				t1++;
			}
			if (xmmodel->param_type[j] == 2) {
				xmmodel->type_index[2][t2] = j;
				t2++;
			}
			if (xmmodel->param_type[j] == 3) {
				xmmodel->type_index[3][t3] = j;
				t3++;
			}
            if (xmmodel->param_type[j] == 4) {
                xmmodel->type_index[4][t4] = j;
                t4++;
            }
		}
	} else {
        xmmodel->index_type = g_new0(int*, 5);
        xmmodel->type_index = g_new0(int*, 5);
		xmmodel->index_type[0] = g_new0(int, xmmodel->has_params[0]);
		xmmodel->type_index[0] = g_new0(int, xmmodel->has_params[0]);
		xmmodel->index_type[1] = g_new0(int, xmmodel->has_params[1]);
		xmmodel->type_index[1] = g_new0(int, xmmodel->has_params[1]);
		xmmodel->index_type[2] = g_new0(int, xmmodel->has_params[2]);
		xmmodel->type_index[2] = g_new0(int, xmmodel->has_params[2]);
        xmmodel->index_type[3] = g_new0(int, xmmodel->has_params[3]);
        xmmodel->type_index[3] = g_new0(int, xmmodel->has_params[3]);
        xmmodel->index_type[4] = g_new0(int, xmmodel->has_params[4]);
        xmmodel->type_index[4] = g_new0(int, xmmodel->has_params[4]);
        j = 0;
		for ( t0 = 0; t0 < xmmodel->has_params[0]; t0++ ) {
			xmmodel->type_index[0][t0] = j;
			xmmodel->param_type[j] = 0;
			xmmodel->tweak[j] = 1;
			j++;
		}
		for ( t1 = 0; t1 < xmmodel->has_params[1]; t1++ ) {
			xmmodel->type_index[1][t1] = j;
			xmmodel->param_type[j] = 1;
			xmmodel->tweak[j] = 1;
			j++;
		}
		for ( t2 = 0; t2 < xmmodel->has_params[2]; t2++ ) {
			xmmodel->type_index[2][t2] = j;
			xmmodel->param_type[j] = 2;
			xmmodel->tweak[j] = 1;
			j++;
		}
		for ( t3 = 0; t3 < xmmodel->has_params[3]; t3++ ) {
			xmmodel->type_index[3][t3] = j;
			xmmodel->param_type[j] = 3;
			xmmodel->tweak[j] = 0;
			j++;
		}
        for ( t4 = 0; t4 < xmmodel->has_params[4]; t4++ ) {
            xmmodel->type_index[4][t4] = j;
            xmmodel->param_type[j] = 4;
            xmmodel->tweak[j] = 1;
            j++;
        }
		k = 0;
		for ( j = 0; j < xmmodel->size; j++) {
			if ( xmmodel->tweak[j] == 1 ) {
				xmmodel->tweak_index[k] = j;
				k++;
			}
		}
	}
	if (xmmodel->type == XmModelIntprt) {
	// create queue with interpreters
		xmmodel->queue_intprts = g_async_queue_new();
		GAsyncQueue *q = xmmodel->queue_intprts;
		int i;
		Interpreter* intprt;
		for(i = 0; i < xmmodel->num_threads; i++){
	        intprt = init_interpreter(xmmodel);
    		g_async_queue_push(q, (gpointer)intprt);
		}
	}
	return ret_val;
}

int xm_model_get_next_free_conn(int connected_to_conn, XmModel*xmmodel)
{
	int curr, i;
	int ret_val = connected_to_conn;
	int connected_to;
	curr = connected_to_conn;
	for ( i = 0; i < xmmodel->size; i++ ) {
		curr++;
		if ( curr > xmmodel->size - 1 ) {
			curr = 0;
		}
		connected_to = xmmodel->lookup[curr];
		if ( xmmodel->parms[curr] != -1 && xmmodel->tweak[curr] == 1 && xmmodel->part[connected_to].checked_num_conn < xmmodel->part[connected_to].num_parms ) {
			ret_val = curr;
			break;
		}
	}
	return ret_val;
}

gint xm_model_sort_atoms (gconstpointer a, gconstpointer b)
{
	int i = GPOINTER_TO_INT(a);
	int j = GPOINTER_TO_INT(b);
	if ( i > j ) {
		return 1;
	} else if ( i < j ) {
		return -1;
	}
	return 0;
}

GString*xm_model_sdf_contents(XmModel*xmmodel)
{
	GQueue *gqueue;
	GList*bonds = NULL, *curr_bond;
	GList*atoms = NULL, *curr_atom;
	GString*file_contents, *file_contents_atoms, *file_contents_bonds;
	int i, j, k;
	int ABORT = 0;
	XmModelConn*conn, *curr_conn;
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		xmmodel->part[i].checked_flag = 0;
		xmmodel->part[i].checked_num_conn = 0;
	}
/*
for ( i = 0; i < xmmodel->size; i++ ) {
	if ( xmmodel->tweak[i] == 1 ) {
		fprintf(stderr, "%d[%d] ", i, xmmodel->parms[i]);
	}
}
fprintf(stderr, "\n");
*/
	int num_atoms = 0;
	int num_bonds = 0;
	file_contents_atoms = g_string_new("");
	file_contents_bonds = g_string_new("");
	gqueue = g_queue_new ();
/*	g_queue_init(gqueue);*/
	for ( i = 0; i < xmmodel->size; i++ ) {
		if ( xmmodel->mask[i] == 1 ) {
			int index_of_part, index_of_part_conn;
			int connected_to, connected_to_conn;
			g_queue_push_head(gqueue, GINT_TO_POINTER(i));
			while ( !g_queue_is_empty(gqueue) ) {
				index_of_part_conn = GPOINTER_TO_INT(g_queue_pop_head(gqueue));
				if ( ( connected_to_conn = xmmodel->parms[index_of_part_conn] ) == -1 ) {
					continue;
				}
				index_of_part = xmmodel->lookup[index_of_part_conn];
				if ( xmmodel->tweak[connected_to_conn] == 0 && index_of_part_conn != xmmodel->parms[connected_to_conn] ) {
					connected_to_conn = xm_model_get_next_free_conn(connected_to_conn, xmmodel);
					xmmodel->parms[index_of_part_conn] = connected_to_conn;
				}
				connected_to = xmmodel->lookup[connected_to_conn];
				if ( xmmodel->part[connected_to].checked_num_conn > xmmodel->part[connected_to].num_parms - 1 || connected_to_conn == index_of_part_conn ) {
					ABORT = 1;
					break;
				}
				xmmodel->parms[connected_to_conn] = -1;
				if ( xmmodel->part[index_of_part].checked_flag == 0 ) {
					xmmodel->part[index_of_part].checked_flag = 1;
					atoms = g_list_prepend(atoms, GINT_TO_POINTER(index_of_part));
					num_atoms++;
					for ( k = 0; k < xmmodel->part[index_of_part].num_parms; k++ ) {
						if ( xmmodel->part[index_of_part].index[k] != index_of_part_conn ) {
							g_queue_push_head(gqueue, GINT_TO_POINTER(xmmodel->part[index_of_part].index[k]));
						}
					}
				}
				if ( xmmodel->part[connected_to].checked_flag == 0 ) {
					xmmodel->part[connected_to].checked_flag = 1;
					atoms = g_list_prepend(atoms, GINT_TO_POINTER(connected_to));
					num_atoms++;
					for ( k = 0; k < xmmodel->part[connected_to].num_parms; k++ ) {
						if ( xmmodel->part[connected_to].index[k] != connected_to_conn ) {
							g_queue_push_head(gqueue, GINT_TO_POINTER(xmmodel->part[connected_to].index[k]));
						}
					}
				}
				conn = xm_model_conn_new(index_of_part, index_of_part_conn, connected_to, connected_to_conn, 0);
				if ( ( curr_bond = g_list_find_custom ( bonds, (gconstpointer) conn, (GCompareFunc)xm_model_conn_comp ) ) != NULL ) {
					curr_conn = (XmModelConn*)(curr_bond->data);
					if ( conn->dest_conn != curr_conn->dest_conn && conn->source_conn != curr_conn->source_conn ) {
						xm_model_conn_free(conn);
						curr_conn->rank++;
						xmmodel->part[curr_conn->dest].checked_num_conn++;
						xmmodel->part[curr_conn->source].checked_num_conn++;
					}
				} else {
					xmmodel->part[conn->dest].checked_num_conn++;
					xmmodel->part[conn->source].checked_num_conn++;
					bonds = g_list_prepend(bonds, (gpointer)conn);
					num_bonds++;
				}
			}
		}
	}
	g_queue_free (gqueue);
	if ( ABORT == 1 ) {
		g_list_free(bonds);
		g_list_free(atoms);
		g_string_free(file_contents_atoms, TRUE);
		g_string_free(file_contents_bonds, TRUE);
		return NULL;
	}
	atoms = g_list_sort (atoms, (GCompareFunc) xm_model_sort_atoms);
	for (i = 0, curr_atom = g_list_first(atoms); i < num_atoms; i++, curr_atom = g_list_next(curr_atom)) {
		j = GPOINTER_TO_INT(curr_atom->data);
		g_string_append_printf(file_contents_atoms, "        0          0    0      %s   0  0  0  0  0  0  0  0  0  0  0  0\n", xmmodel->part[j].name);
	}
	for (i = 0, curr_bond = g_list_first(bonds); i < num_bonds; i++, curr_bond = g_list_next(curr_bond)) {
		conn = (XmModelConn*)(curr_bond->data);
		int source = g_list_index ( atoms, GINT_TO_POINTER(conn->source) );
		int dest = g_list_index ( atoms, GINT_TO_POINTER(conn->dest) );
		g_string_append_printf(file_contents_bonds, "%3d%3d%3d  0  0  0  0\n", source + 1, dest + 1, conn->rank + 1);
		xm_model_conn_free(conn);
	}
	g_list_free(bonds);
	g_list_free(atoms);
	file_contents = g_string_new("Name\n  Soft\n\n");
	g_string_append_printf(file_contents, "%3d%3d  0  0  0  0  0  0  0  0  1 V2000\n", num_atoms, num_bonds);
	file_contents = g_string_append(file_contents, (const gchar*)(file_contents_atoms->str));
	g_string_free(file_contents_atoms, TRUE);
	file_contents = g_string_append(file_contents, (const gchar*)(file_contents_bonds->str));
	g_string_free(file_contents_bonds, TRUE);
	file_contents = g_string_append(file_contents, "M END\n$$$$");
/*
for ( i = 0; i < xmmodel->size; i++ ) {
	if ( xmmodel->tweak[i] == 1 ) {
		fprintf(stderr, "%d[%d] ", i, xmmodel->parms[i]);
	}
}
fprintf(stderr, "\n");
fprintf (stderr, "%s", file_contents->str);*/
	return file_contents;
}

GString*xm_model_gcdm_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k, target;
	int num_targets, num_regulators;
	file_contents = g_string_new("");
	g_string_append_printf(file_contents, "$input\n");
	i = 0;
	k = 0;
	g_string_append_printf(file_contents, "%s\n", xmmodel->part[i].name);
	num_targets = xmmodel->part[i].num_parms;
	for ( j = 0; j < num_targets; j++ ) {
		if ( (buf = param2str(xmmodel, k)) == NULL) {
			g_string_free(file_contents, TRUE);
			return NULL;
		}
		file_contents = g_string_append(file_contents, buf);
		g_free(buf);
		file_contents = g_string_append_c(file_contents, ' ');
		k++;
	}
	file_contents = g_string_append_c(file_contents, '\n');
	for ( i = 1; i < 5; i++ ) {
		g_string_append_printf(file_contents, "%s\n", xmmodel->part[i].name);
		if ( xmmodel->part[i].num_parms > 0 ) {
			num_regulators = xmmodel->part[i].num_parms / num_targets;
			for ( target = 0; target < num_targets; target++ ) {
				for ( j = 0; j < num_regulators; j++ ) {
					if ( (buf = param2str(xmmodel, k)) == NULL) {
						g_string_free(file_contents, TRUE);
						return NULL;
					}
					file_contents = g_string_append(file_contents, buf);
					g_free(buf);
					file_contents = g_string_append_c(file_contents, ' ');
					k++;
				}
				file_contents = g_string_append_c(file_contents, '\n');
			}
		}
	}
	for ( i = 5; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s\n", xmmodel->part[i].name);
		for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, ' ');
			k++;
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	g_string_append_printf(file_contents, "$$\n");
	return file_contents;
}

GString*xm_model_gemstat_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k;
	file_contents = g_string_new("");
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s ", xmmodel->part[i].name);
		for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
			k = xmmodel->part[i].index[j];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, ' ');
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	return file_contents;
}

GString*xm_model_dgemstat_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k;
	file_contents = g_string_new("");
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s ", xmmodel->part[i].name);
		for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
			k = xmmodel->part[i].index[j];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, ' ');
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s ", xmmodel->part[i].name);
		for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
			k = xmmodel->part[i].index[j];
			if ( (buf = dparam2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, ' ');
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	return file_contents;
}

GString*xm_model_subset_contents(XmModel*xmmodel)
{
	GString*file_contents;
	int i, j, k;
	file_contents = g_string_new("");
	k = 0;
	for ( i = 0; i < xmmodel->size; i++ ) {
		if ( xmmodel->tweak[i] == 1 ) {
			if ( xmmodel->mask[k] == 1 ) {
				g_string_append_printf(file_contents, "%d\n", xmmodel->parms[i]);
				k++;
			}
		} else {
			g_string_append_printf(file_contents, "%d\n", xmmodel->parms[i]);
		}
	}
	return file_contents;
}

GString*xm_model_subsubset_contents(XmModel*xmmodel)
{
	GString*file_contents;
	int i, j, k;
	int main_index = 0;
	int main_index_parm = xmmodel->parms[main_index];
	if ( main_index_parm == main_index ) {
		main_index_parm = xmmodel->parms[xmmodel->size - 1];
		xmmodel->parms[xmmodel->size - 1] = main_index;
		xmmodel->parms[main_index] = main_index_parm;
	}
    file_contents = g_string_new("");
	for ( i = 1; i < main_index_parm + 1; i++ ) {
        g_string_append_printf(file_contents, "%d\n", xmmodel->parms[i]);
	}
	return file_contents;
}

GString*xm_model_ini_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k;
	file_contents = g_string_new("[DEEP]\n");
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		if ( xmmodel->part[i].num_parms == 1 ) {
			g_string_append_printf(file_contents, "%s=", xmmodel->part[i].name);
			k = xmmodel->part[i].index[0];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, '\n');
		} else if ( xmmodel->part[i].num_parms > 1 ) {
			g_string_append_printf(file_contents, "%s/size=%d\n", xmmodel->part[i].name, xmmodel->part[i].num_parms);
			for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
				k = xmmodel->part[i].index[j];
				if ( (buf = param2str(xmmodel, k)) == NULL) {
					g_string_free(file_contents, TRUE);
					return NULL;
				}
				g_string_append_printf(file_contents, "%s/%d/value=", xmmodel->part[i].name, j + 1);
				file_contents = g_string_append(file_contents, buf);
				g_free(buf);
				file_contents = g_string_append_c(file_contents, '\n');
			}
		}
	}
	return file_contents;
}

GString*xm_model_octave_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k;
	file_contents = g_string_new("");
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s =", xmmodel->part[i].name);
		if ( xmmodel->part[i].num_parms == 1 ) {
			k = xmmodel->part[i].index[0];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			file_contents = g_string_append_c(file_contents, ';');
		} else if ( xmmodel->part[i].num_parms > 1 ) {
			g_string_append_printf(file_contents, "[");
			for ( j = 0; j < xmmodel->part[i].num_parms; j++ ) {
				k = xmmodel->part[i].index[j];
				if ( (buf = param2str(xmmodel, k)) == NULL) {
					g_string_free(file_contents, TRUE);
					return NULL;
				}
				file_contents = g_string_append(file_contents, buf);
				g_free(buf);
				file_contents = g_string_append_c(file_contents, ' ');
			}
			g_string_append_printf(file_contents, "];");
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	return file_contents;
}

GString*xm_model_r_contents(XmModel*xmmodel)
{
	GString*file_contents;
	gchar*buf = NULL;
	int i, j, k;
	file_contents = g_string_new("");
	for ( i = 0; i < xmmodel->num_parts; i++ ) {
		g_string_append_printf(file_contents, "%s =", xmmodel->part[i].name);
		if ( xmmodel->part[i].num_parms == 1 ) {
			k = xmmodel->part[i].index[0];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
		} else if ( xmmodel->part[i].num_parms > 1 ) {
			g_string_append_printf(file_contents, "c(");
			for (j = 0; j < xmmodel->part[i].num_parms - 1; j++) {
				k = xmmodel->part[i].index[j];
				if ( (buf = param2str(xmmodel, k)) == NULL) {
					g_string_free(file_contents, TRUE);
					return NULL;
				}
				file_contents = g_string_append(file_contents, buf);
				g_free(buf);
				file_contents = g_string_append_c(file_contents, ',');
			}
			j = xmmodel->part[i].num_parms - 1;
            k = xmmodel->part[i].index[j];
			if ( (buf = param2str(xmmodel, k)) == NULL) {
				g_string_free(file_contents, TRUE);
				return NULL;
			}
			file_contents = g_string_append(file_contents, buf);
			g_free(buf);
			g_string_append_printf(file_contents, ")");
		}
		file_contents = g_string_append_c(file_contents, '\n');
	}
	return file_contents;
}

gchar *xm_model_convert_parms_to_sdf(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_sdf_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_gcdm(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_gcdm_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_gemstat(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_gemstat_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_dgemstat(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_dgemstat_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_subset(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_subset_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_subsubset(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_subsubset_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_ini(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_ini_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_octave(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_octave_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

gchar *xm_model_convert_parms_to_r(gpointer *user_data, GError **err)
{
	XmModel*xmmodel = (XmModel*)user_data;
	gint fhandle;
	const gchar *tmpl = NULL;
	gchar *name_used = NULL;
	GError *gerror = NULL;
	GString*file_contents;
	g_return_val_if_fail (err == NULL || *err == NULL, NULL);
	if ( ( fhandle = g_file_open_tmp (tmpl, &name_used, &gerror) ) == -1 ) {
		g_propagate_error (err, gerror);
		return NULL;
	}
	close(fhandle);
	if ( ( file_contents = xm_model_r_contents(xmmodel) ) == NULL ) {
		g_unlink(name_used);
		return NULL;
	}
	if ( !g_file_set_contents ((const gchar *)name_used, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
		g_string_free(file_contents, TRUE);
		g_propagate_error (err, gerror);
		g_unlink(name_used);
		return NULL;
	}
	g_string_free(file_contents, TRUE);
	return name_used;
}

XmModelConn*xm_model_conn_new(int source, int source_conn, int dest, int dest_conn, int rank)
{
	XmModelConn*conn;
	conn = (XmModelConn*)malloc(sizeof(XmModelConn));
	conn->source = source;
	conn->source_conn = source_conn;
	conn->dest = dest;
	conn->dest_conn = dest_conn;
	conn->rank = rank;
	return conn;
}

int xm_model_conn_comp (gconstpointer a, gconstpointer b)
{
	XmModelConn*conn1 = (XmModelConn*)a;
	XmModelConn*conn2 = (XmModelConn*)b;
	if ( conn1->source < conn2->source && conn1->dest < conn2->dest ) {
		return -1;
	} else if ( conn1->source == conn2->source && conn1->dest == conn2->dest ) {
		return 0;
	} else if ( conn1->source == conn2->dest && conn1->dest == conn2->source ) {
		return 0;
	}
	return 1;
}

int xm_model_conn_comp_ports (gconstpointer a, gconstpointer b)
{
	XmModelConn*conn1 = (XmModelConn*)a;
	XmModelConn*conn2 = (XmModelConn*)b;
	if ( conn1->source < conn2->source && conn1->dest < conn2->dest ) {
		return -1;
	} else if ( conn1->source == conn2->source && conn1->dest == conn2->dest && conn1->source_conn == conn2->source_conn && conn1->dest_conn == conn2->dest_conn ) {
		return 0;
	} else if ( conn1->source == conn2->dest && conn1->dest == conn2->source && conn1->source_conn == conn2->dest_conn && conn1->dest_conn == conn2->source_conn) {
		return 0;
	}
	return 1;
}

void xm_model_conn_free(XmModelConn*conn)
{
	free(conn);
}

void xm_model_save(XmModel*xmmodel, gchar*filename)
{
	GString*file_contents = NULL;
	int i;
	GError*gerror = NULL;
	xm_model_convert_param_type(xmmodel);
	if ( !g_strcmp0 ( xmmodel->convert, "sdf" ) ) {
		if ( ( file_contents = xm_model_sdf_contents(xmmodel) ) == NULL ) {
			g_error("Can't get sdf contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "gcdm" ) ) {
		if ( ( file_contents = xm_model_gcdm_contents(xmmodel) ) == NULL ) {
			g_error("Can't get gcdm contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "gemstat" ) ) {
		if ( ( file_contents = xm_model_gemstat_contents(xmmodel) ) == NULL ) {
			g_error("Can't get gemstat contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "dgemstat" ) ) {
		if ( ( file_contents = xm_model_dgemstat_contents(xmmodel) ) == NULL ) {
			g_error("Can't get dgemstat contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "subset" ) ) {
		if ( ( file_contents = xm_model_subset_contents(xmmodel) ) == NULL ) {
			g_error("Can't get subset contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "subsubset" ) ) {
		if ( ( file_contents = xm_model_subsubset_contents(xmmodel) ) == NULL ) {
			g_error("Can't get subsubset contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "octave" ) ) {
		if ( ( file_contents = xm_model_octave_contents(xmmodel) ) == NULL ) {
			g_error("Can't get octave contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "ini" ) ) {
		if ( ( file_contents = xm_model_ini_contents(xmmodel) ) == NULL ) {
			g_error("Can't get ini contents");
		}
	} else if ( !g_strcmp0 ( xmmodel->convert, "r" ) ) {
		if ( ( file_contents = xm_model_r_contents(xmmodel) ) == NULL ) {
			g_error("Can't get r contents");
		}
	} else {
		file_contents = g_string_new("");
		for ( i = 0; i < xmmodel->size; i++ ) {
			g_string_append_printf(file_contents, "%f ", xmmodel->dparms[i]);
		}
	}
	if ( file_contents != NULL ) {
		if ( !g_file_set_contents ((const gchar *)filename, (const gchar *)(file_contents->str), (gssize)(file_contents->len), &gerror) ) {
			g_string_free(file_contents, TRUE);
			g_error ("%s", gerror->message);
		}
		g_string_free(file_contents, TRUE);
	}
}

GString *xm_model_score_double_to_string(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	GString *params = g_string_new("");
	int i;
	xm_model_set_dparms(xmmodel, x);
	for ( i = 0; i < xmmodel->size; i++ ) {
		g_string_append_printf(params, "x[%d]:%.*f ", i, xmmodel->b_precision, xmmodel->dparms[i]);
	}
	return params;
}

GString *xm_model_score_int_to_string(gpointer user_data, double*x)
{
	XmModel *xmmodel = (XmModel *)user_data;
	double val = G_MAXDOUBLE;
	GString *params = g_string_new("");
	int i;
	if ( !g_strcmp0 ( xmmodel->convert, "subsubset" ) ) {
	    int main_index = 0;
        int main_index_parm = xmmodel->parms[main_index];
        for ( i = 1; i < main_index_parm + 1; i++ ) {
            g_string_append_printf(params, "x[%d]:%d ", i, xmmodel->parms[i]);
        }
	} else {
        for ( i= 0; i < xmmodel->size; i++ ) {
            g_string_append_printf(params, "x[%d]:%d ", i, xmmodel->parms[i]);
        }
	}
	return params;
}
