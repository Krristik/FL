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

#include "dparchive.h"
#include <time.h>
#include <math.h>

DpArchive* dp_archive_new(int size, int individSize, DpArchiveMode mode){
    DpArchive* archive = malloc(sizeof(DpArchive));
    archive->mode = mode;
    if (archive->mode == DIFF) {
        archive->difference_vectors = (DifferenceVector**)malloc(sizeof(DifferenceVector) * 2 * size);
        for (int i = 0; i < size * 2; i++){
            archive->difference_vectors[i] = (DifferenceVector*)malloc(sizeof(DifferenceVector));
            archive->difference_vectors[i]->value = (double*)malloc(sizeof(double) * individSize);
        }
        archive->last_index = size;
        archive->max_size = size * 2;

    } else if (archive->mode == INDIVIDS) {
        archive->individ = NULL;
        archive->xamean = (DifferenceVector*)malloc(sizeof(DifferenceVector));
        archive->xamean->value = (double*)malloc(sizeof(double) * individSize);
        archive->last_index = 0;
        archive->max_size = round(size * 2.6);
    }
    archive->iter = 0;
    archive->generation = 0;
    archive->indSize = individSize;
    return archive;
}

DpArchive* dp_archive_init(DpPopulation* population, DpArchiveMode mode){
    DpArchive* newArchive = dp_archive_new(population->size, population->individ[0]->size, mode);
    if (newArchive->mode == DIFF) {
        for (int i = 0; i < population->size * 2; i++){
            int i1, i2;
            i1 = i < population->size ? i : i - population->size;
            i2 = i1 == population->size - 1 ? 0 : i1 + 1;
            for (int j = 0; j < population->individ[0]->size; j++){
                newArchive->difference_vectors[i]->value[j] = population->individ[i1]->x[j] - population->individ[i2]->x[j];
            }
            newArchive->difference_vectors[i]->generation = 0;
        }
    } else if (newArchive->mode == INDIVIDS) {
        for (int j = 0; j < newArchive->indSize; j++) {
            newArchive->xamean->value[j] = 0;
        }
    }
    return newArchive;
}

void dp_archive_delete(DpArchive* archive){
    if (archive->mode == DIFF) {
        free(archive->difference_vectors);
    } else if (archive->mode == INDIVIDS) {
        for (int i=0; i < archive->last_index; i++) {
            dp_individ_delete(archive->individ[i]);
        }
        free(archive->individ);
        free(archive->xamean);
    }
}

void add_difference_vector(DpArchive* archive, DifferenceVector* new_vector){
    archive->difference_vectors[archive->last_index] = new_vector;
    archive->last_index++;
}

void insert_individ_in_archive(DpIndivid*** individs, int* size, int* max_size, DpIndivid* to_archive) {
    //проверка на заполненность архива: если массив заполнен полностью, случайный элемент исключается из архива
    if (*size == *max_size) {
        static int seeded = 0;
        if (!seeded) {
            srand(time(NULL));
            seeded = 1;
        }

        int remove_idx = rand() % *size;
        dp_individ_delete((*individs)[remove_idx]);
        for (int i = remove_idx; i < *size - 1; i++) {
            (*individs)[i] = (*individs)[i + 1];
        }
        (*size)--;
    }
    //бинарный поиск индекса вставки
    int low = 0, high = *size - 1, pos = *size;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if ((*individs)[mid]->cost == to_archive->cost) {
            pos = mid;
            break;
        } else if ((*individs)[mid]->cost < to_archive->cost) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    if (pos == *size) pos = low;
    //перевыделение памяти для нового элемента, если архив не заполнен до конца
    if (*size < *max_size) {
        DpIndivid **tmp = realloc(*individs, (*size + 1) * sizeof(DpIndivid*));
        if (!tmp) return NULL;
        *individs = tmp;
    }
    //вставка нового индивида
    for (int i = *size; i > pos; i--) {
        (*individs)[i] = (*individs)[i - 1];
    }
    (*individs)[pos] = to_archive;
    (*size)++;
}

void add_individ(DpArchive* archive, DpIndivid* new_individ, int targets_size, int precond_size){
    DpIndivid* to_archive = dp_individ_new(archive->indSize, targets_size, precond_size);
    dp_individ_copy_values(to_archive, new_individ);
    insert_individ_in_archive(&(archive->individ), &(archive->last_index), &(archive->max_size), to_archive);
}

void shuffle_archive(DpArchive *oldArchive) {

    srand(time(NULL));
    for (int i = oldArchive->max_size - 1; i >= 1; i--) {
        int j = rand() % (i + 1);
        DifferenceVector* tmp = oldArchive->difference_vectors[j];
        oldArchive->difference_vectors[j] = oldArchive->difference_vectors[i];
        oldArchive->difference_vectors[i] = tmp;
    }
}

void calculate_xamean(DpArchive* archive) {
    double w, smean, sum_of_logs = 0;
    archive->m = round(archive->e * archive->last_index);
    if (archive->m < 1) {
        archive->m = 1;
    }
    for (int k = 0; k < archive->m; k++) {
        sum_of_logs += (log(archive->m + 0.5) - log(k + 1));
    }
    for (int j = 0; j < archive->indSize; j++) {
        smean = 0;
        for (int i = 0; i < archive->m; i++) {
            w = (log(archive->m + 0.5) - log(i + 1))/sum_of_logs;
            smean += w * archive->individ[i]->x[j];
        }
        archive->xamean->value[j] = smean;
    }
}
