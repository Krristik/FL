/***************************************************************************
 *            xmtranslate.c
 *
 *  Tue Mar 27 12:24:31 2012
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
#include <glib.h>
#include "xmmodel.h"
#include "xmtranslate.h"
#include "dpevaluation.h"

#ifdef GIO_STANDALONE_SOURCE
#define g_strcmp0(str1, str2) strcmp(str1, str2)
#endif

DpEvaluation *xm_translate_parms(XmModel *xmmodel)
{
	int i, j;
	int index = 0;
	DpEvaluation*plist;
	plist = (DpEvaluation*)malloc(sizeof(DpEvaluation));
	plist->size = 0;
	plist->points = NULL;
	for ( i = 0; i < xmmodel->size; i++ ) {
		if ( xmmodel->tweak[i] == 1 ) {
			plist->points = (DpEvaluationPoint**)realloc( plist->points, (index + 1) * sizeof(DpEvaluationPoint*) );
			plist->points[index] = (DpEvaluationPoint*)malloc( sizeof(DpEvaluationPoint) );
			plist->points[index]->param = &(xmmodel->dparms[i]);
			plist->points[index]->index = i;
			plist->points[index]->upper = xmmodel->hbound[ i ];
			plist->points[index]->lower = xmmodel->lbound[ i ];
			plist->points[index]->scale = xmmodel->scale[ i ];
			plist->points[index]->limited = xmmodel->limited[ i ];
            plist->points[index]->dedupl = xmmodel->dedupl[ i ];
            index++;
		}
	}
	plist->size  = index;
	return plist;
}

DpFunc xm_translate_get_score_func_by_name ( char*name )
{
	if ( !g_strcmp0(name, "doubletoint") ) {
		return (DpFunc)xm_model_parms_double_to_int;
	} else if ( !g_strcmp0(name, "doubletoindex") ) {
		return (DpFunc)xm_model_parms_double_to_index;
	} else if ( !g_strcmp0(name, "scoreint") ) {
		return (DpFunc)xm_model_score_int;
	} else if ( !g_strcmp0(name, "scoredouble") ) {
		return (DpFunc)xm_model_score_double;
	} else if ( !g_strcmp0(name, "objfunc") ) {
		return (DpFunc)xm_model_objfunc;
	} else if ( !g_strcmp0(name, "readpenalty") ) {
		return (DpFunc)xm_model_read_penalty;
	} else if ( !g_strcmp0(name, "barrier") ) {
		return (DpFunc)xm_model_barrier_penalty;
	} else if ( !g_strcmp0(name, "primedouble") ) {
		return (DpFunc)xm_model_prime_double;
	} else if ( !g_strcmp0(name, "primeint") ) {
		return (DpFunc)xm_model_prime_int;
	}
}

void xm_translate_score(DpTarget*htarget, XmModel *xmmodel)
{
	int i;
	htarget->target->f = xm_translate_get_score_func_by_name(htarget->target->name);
	htarget->user_data = (gpointer)(xmmodel);
	htarget->copy_model = xm_model_copy_values;
	htarget->update_model = xm_model_update_values;
	for ( i = 0; i < htarget->penalty_size; i++ ) {
		htarget->penalty[i]->f = xm_translate_get_score_func_by_name(htarget->penalty[i]->name);
	}
	htarget->params_to_string = xm_model_score_double_to_string;
	for ( i = 0; i < htarget->precond_size; i++ ) {
		htarget->precond[i]->f = xm_translate_get_score_func_by_name(htarget->precond[i]->name);
        if ( !g_strcmp0(htarget->precond[i]->name, "doubletoint") ) {
            htarget->params_to_string = xm_model_score_int_to_string;
        } else if ( !g_strcmp0(htarget->precond[i]->name, "doubletoindex") ) {
            htarget->params_to_string = xm_model_score_int_to_string;
        }
	}
	if ( htarget->prime != NULL ) {
		htarget->prime->f = xm_translate_get_score_func_by_name(htarget->prime->name);
	}
	for ( i = 0; i < htarget->constreq_size; i++ ) {
		htarget->constreq[i]->f = xm_translate_get_score_func_by_name(htarget->constreq[i]->name);
	}
	for ( i = 0; i < htarget->constrneq_size; i++ ) {
		htarget->constrneq[i]->f = xm_translate_get_score_func_by_name(htarget->constrneq[i]->name);
	}
}

void xm_retranslate_precond(DpTarget*htarget, XmModel *xmmodel)
{
	int i;
	for ( i = 0; i < htarget->precond_size; i++ ) {
        if ( !g_strcmp0(htarget->precond[i]->name, "doubletoint") ) {
            xm_model_dparms_to_int((gpointer) xmmodel);
        } else if ( !g_strcmp0(htarget->precond[i]->name, "doubletoindex") ) {
            xm_model_dparms_to_index((gpointer) xmmodel);
        }
	}
}
