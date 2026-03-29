/***************************************************************************
 *            dparchive.c
 *
 *  Copyright  2022  Andrei Ivanov
 *  Andrey Ivanov <ezhva2010@gmail.com>
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

#ifndef DPARCHIVE_H
#define DPARCHIVE_H

#ifdef __cplusplus
extern "C"
{
#endif


#include <stdio.h>
#include <stdlib.h>
#include <dppopulation.h>

//режим работы архива - векторы разностей или индивиды
typedef enum DpArchiveMode {
    DIFF,
    INDIVIDS
} DpArchiveMode;

typedef struct {
    double* value;
    int generation;
} DifferenceVector;

typedef struct DpArchive{
    DpArchiveMode mode;
    int max_size;
    int last_index;
    int iter;
    int generation;
    DifferenceVector* xamean;
    double e;
    int m;
    int indSize;
    DifferenceVector** difference_vectors;
    DpIndivid** individ;
} DpArchive;

DpArchive* dp_archive_new(int size, int individSize, DpArchiveMode mode);

DpArchive* dp_archive_init(DpPopulation* population, DpArchiveMode mode);

void dp_archive_delete(DpArchive* archive);

void add_difference_vector(DpArchive* archive, DifferenceVector* new_vector);

void add_individ(DpArchive* archive, DpIndivid* new_individ, int targets_size, int precond_size);

void shuffle_archive(DpArchive* oldArchive);

void calculate_xamean(DpArchive* archive);

#ifdef __cplusplus
}
#endif

#endif /* _DP_RECOMBINATION_H */
