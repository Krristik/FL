/***************************************************************************
 *            xmtranslate.h
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

#ifndef _XM_TRANSLATE_H
#define _XM_TRANSLATE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include "xmmodel.h"
#include "dpevaluation.h"

DpEvaluation *xm_translate_parms(XmModel *xmmodel);

DpFunc xm_translate_get_score_func_by_name ( char*name );

void xm_translate_score(DpTarget*htarget, XmModel *xmmodel);

void xm_retranslate_precond(DpTarget*htarget, XmModel *xmmodel);

#ifdef __cplusplus
}
#endif

#endif /* _XM_TRANSLATE_H */
