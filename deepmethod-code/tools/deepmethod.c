/***************************************************************************
 *            dpdeepctl.c
 *
 *  Fri Mar 23 16:45:04 2012
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

#if defined(_MSC_VER)

#else
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#ifdef GIO_STANDALONE_SOURCE
#include <gio/gfile.h>
#else
#include <gio/gio.h>
#endif
#ifdef MPIZE
#include <mpi.h>
#endif
#include "dpdeep.h"
#include "dprecombination.h"
#include "dpopt.h"
#include "dpsettings.h"
#include "xmmodel.h"
#include "xmtranslate.h"

static char*default_name;
static char*hopt_filename;
static char*model_file;
static char*model_group;
static char*settings_file;
static char*settings_group;
static char*target_file;
static char*target_group;
static char*output_file;
static char*operation;
static int monitor;
static int hoptdebug;
static gboolean dont_write_params = FALSE;

/*
 * Standard gettext macros.
*/
#ifdef ENABLE_NLS
#  include <libintl.h>
#  undef _
#  define _(String) dgettext (PACKAGE, String)
#  ifdef gettext_noop
#    define N_(String) gettext_noop (String)
#  else
#    define N_(String) (String)
#  endif
#else
#  define textdomain(String) (String)
#  define gettext(String) (String)
#  define dgettext(Domain,Message) (Message)
#  define dcgettext(Domain,Message,Type) (Message)
#  define bindtextdomain(Domain,Directory) (Domain)
#  define _(String) (String)
#  define N_(String) (String)
#endif

static gboolean
option_version_cb (const gchar *option_name,
                   const gchar *value,
                   gpointer     data,
                   GError     **error)
{
#if defined(_MSC_VER)
  g_print ("%s %s\n", _("DEEP optimizer"));
#else
  g_print ("%s %s\n", _("DEEP optimizer, version "), VERSION);
#endif
  exit (0);
  return FALSE;
}

static GOptionEntry entries[] =
{
	{ "default-name", 0, 0, G_OPTION_ARG_STRING, &default_name, N_("File name for everything and default group names"), N_("FILENAME") },
	{ "model-file", 0, 0, G_OPTION_ARG_STRING, &model_file, N_("File name where model is described"), N_("FILENAME") },
	{ "model-group", 0, 0, G_OPTION_ARG_STRING, &model_group, N_("Group name where model is described"), N_("GROUP") },
	{ "settings-file", 0, 0, G_OPTION_ARG_STRING, &settings_file, N_("File name where settings are described"), N_("FILENAME") },
	{ "settings-group", 0, 0, G_OPTION_ARG_STRING, &settings_group, N_("Group name where settings are described"), N_("GROUP") },
	{ "target-file", 0, 0, G_OPTION_ARG_STRING, &target_file, N_("File name where target is described"), N_("FILENAME") },
	{ "target-group", 0, 0, G_OPTION_ARG_STRING, &target_group, N_("Group name where target is described"), N_("GROUP") },
	{ "output-file", 0, 0, G_OPTION_ARG_STRING, &output_file, N_("File name to write"), N_("FILENAME") },
	{ "operation", 0, 0, G_OPTION_ARG_STRING, &operation, N_("What to do"), N_("OPERATION") },
	{ "monitor", 0, 0, G_OPTION_ARG_INT, &monitor, N_("Extract line from log"), N_("NUMBER") },
	{ "debug", 0, 0, G_OPTION_ARG_INT, &hoptdebug, N_("Print out run-time info"), N_("NUMBER") },
	{ "version", 0, G_OPTION_FLAG_NO_ARG | G_OPTION_FLAG_HIDDEN, G_OPTION_ARG_CALLBACK, option_version_cb, NULL, NULL },
	{ "dont-write-params", 'd', 0, G_OPTION_ARG_NONE, &dont_write_params, N_("Suppress parameters file"), NULL },
	{ NULL }
};

#ifdef GIO_STANDALONE_SOURCE
#define g_strcmp0(str1, str2) strcmp(str1, str2)
#endif

int main(int argc, char **argv)
{
	XmModel*xmmodel;
    DpOpt *hopt;
	DpSettings*dpsettings = dp_settings_new();
	DpEvaluation*heval;
	DpTarget*htarget;
	int world_id = 0, world_count = 1;
	double dval;
	GOptionContext *context;
	GKeyFile*gkf;
	GError *gerror = NULL;
#ifdef ENABLE_NLS
	bindtextdomain (GETTEXT_PACKAGE, DEEPMETHOD_LOCALEDIR);
	bind_textdomain_codeset (GETTEXT_PACKAGE, "UTF-8");
	textdomain (GETTEXT_PACKAGE);
#endif
	hopt_filename = NULL;
	context = g_option_context_new (_("- DEEP optimizer"));
	g_option_context_add_main_entries(context, (const GOptionEntry *)entries, NULL);
	g_option_context_set_ignore_unknown_options(context, TRUE);
#if defined(_MSC_VER)

#else
	g_option_context_set_translation_domain(context, GETTEXT_PACKAGE);
#endif
	if (!g_option_context_parse (context, &argc, &argv, &gerror)) {
		g_error (_("option parsing failed: %s\n"), gerror->message);
	}
	g_option_context_free (context);
	if ( model_file == NULL || model_group == NULL ) {
		if ( default_name == NULL ) {
			g_error(_("%s called with wrong options for model"), g_get_prgname());
		} else {
			model_file = g_strdup(default_name);
			model_group = g_strdup("default_model");
		}
	}
	if ( settings_file == NULL || settings_group == NULL ) {
		if ( default_name == NULL ) {
			g_error(_("%s called with wrong options for settings"), g_get_prgname());
		} else {
			settings_file = g_strdup(default_name);
			settings_group = g_strdup("default_settings");
		}
	}
	if ( target_file == NULL || target_group == NULL ) {
		if ( default_name == NULL ) {
			g_error(_("%s called with wrong options for target"), g_get_prgname());
		} else {
			target_file = g_strdup(default_name);
			target_group = g_strdup("default_target");
		}
	}
	if ( output_file == NULL ) {
		if ( default_name == NULL ) {
			g_error(_("%s called with wrong options for output"), g_get_prgname());
		} else {
			output_file = g_strdup_printf("%s-deep-output", default_name);
			hopt_filename = g_strdup(default_name);
		}
	} else {
		hopt_filename = g_strdup(output_file);
	}
	if ( operation == NULL ) {
		if ( default_name == NULL ) {
			g_error(_("%s called with wrong options for operation"), g_get_prgname());
		} else {
			operation = g_strdup("optimize");
		}
	}
	gkf = dp_settings_init(settings_file, settings_group, dpsettings, &gerror);
	if ( gerror != NULL ) {
		g_error(_("Settings init:%s"), gerror->message);
	}
	xmmodel = xm_model_new();
	xm_model_init(model_file, model_group, xmmodel, &gerror);
#ifdef MPIZE
/* MPI initialization steps */
	MPI_Init(&argc, &argv);     /* initializes the MPI execution environment */
	MPI_Comm_size(MPI_COMM_WORLD, &world_count);         /* number of processors? */
	MPI_Comm_rank(MPI_COMM_WORLD, &world_id);          /* ID of local processor? */
#endif
	if ( !g_strcmp0 ( operation, "optimize" ) ) {
		heval = xm_translate_parms(xmmodel);
		htarget = dp_target_new();
		dp_target_init(target_file, target_group, htarget, &gerror);
		if ( gerror != NULL ) {
			g_error("%s", gerror->message);
		}
		xm_translate_score(htarget, xmmodel);
		hopt = dp_opt_init(heval, htarget, world_id, world_count, hopt_filename, dpsettings->stop_type, dpsettings->criterion, dpsettings->stop_count, dpsettings->pareto_all, dpsettings->precision);
		hopt->debug = hoptdebug;
		dp_settings_process_run(dpsettings, gkf, settings_group, hopt, world_id, heval, htarget, &gerror);
		if ( gerror != NULL ) {
			g_error(_("Settings process init:%s"), gerror->message);
		}
		dp_opt_run(hopt);
	} else if ( !g_strcmp0 ( operation, "monitor" ) ) {
		heval = xm_translate_parms(xmmodel);
		htarget = dp_target_new();
		dp_target_init(target_file, target_group, htarget, &gerror);
		if ( gerror != NULL ) {
			g_error("%s", gerror->message);
		}
		xm_translate_score(htarget, xmmodel);
		hopt = dp_opt_init(heval, htarget, world_id, world_count, hopt_filename, dpsettings->stop_type, dpsettings->criterion, dpsettings->stop_count, dpsettings->pareto_all, dpsettings->precision);
		hopt->debug = hoptdebug;
		dp_opt_monitor(hopt, monitor, &gerror);
		if ( gerror != NULL ) {
			g_error(_("Monitor error:%s"), gerror->message);
		}
	}
#ifdef MPIZE
	if ( world_id == 0 ) { /* master process */
#endif
		if (!dont_write_params) {
			if ( hopt->hloop->stop_flag != DP_LOOP_EXIT_ERROR ) {
				xm_retranslate_precond(htarget, xmmodel);
				xm_model_save(xmmodel, output_file);

			} else {
				g_error(_("Loop finished with an unknown error.\nOutput not produced."));
			}
        }
#ifdef MPIZE
	}
	MPI_Finalize();
#endif
	return 0;
}

