/***************************************************************************
 *            dprecombination.c
 *
 *  Fri Mar 23 15:39:19 2012
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
#include <glib.h>
#include <time.h>
#include "dpindivid.h"
#include "dppopulation.h"
#include "dprecombination.h"

void dp_triangular_rec(double*x,  double*x_1,  double*x_2,  double*x_3, double*F, double*p, int end_index){
    //всегда преобразуем F по тангенсу и используем уже преобразованный

    double norm = p[0] + p[1] + p[2];
    double w[3] = {p[0]/ norm, p[1]/ norm, p[2]/ norm};
    for(int i = 0; i < end_index; i++) {
        x[i] = w[0] * x_1[i] + w[1] * x_2[i] + w[2] * x_3[i] +\
                F[0] * (x_1[i] - x_2[i]) + F[1] * (x_1[i] - x_3[i]) +\
                F[2] * (x_2[i] - x_3[i]);
    }
}

int cmp_double(const double*i, const double*j){
    return *i < *j;
}

void ind_triand_init(DpPopulation *population, DpRecombinationControl *control, GRand *hrand){
    if (control->strategy != DE_3_triangular)
        return;

    for(int i = 0; i < population->size; i++){
        population->individ[i]->p[0] = 1;
        population->individ[i]->p[1] = g_rand_double_range(hrand, 0.75, 1);
        population->individ[i]->p[2] = g_rand_double_range(hrand, 0.5, population->individ[i]->p[1]);

        if(control->gamma != 0){
            for(int j = 0; j < TRIANGULAR; j++) {
                population->individ[i]->F[j] = g_rand_double_range(hrand, control->f_inf, control->f_supp);
            }
        } else {
            for(int j = 0; j < TRIANGULAR; j++){
                population->individ[i]->F[j] = control->f[0];
            }
        }
    }
}

//Transfer archive vector for writing and trial vector for reading
//void dp_individ_recombination(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index)
void dp_individ_recombination(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,
                              DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, DpIndivid*input_5,
                              int start_index, int end_index,
                              DifferenceVector* vectorWrite, DifferenceVector* vectorRead)
{
    int i;
    int L;
    int flag;
    double u, phi, alpha;
    int FEs = control->FEs, MaxFEs = control->MaxFEs;
    double Fw, SPg;

    DpIndivid* tmp[3] = {input_1, input_2, input_3};

    DifferenceVector *archiveVector;

	if (control->use_archive_prob > 0) {
	    for (int i = 0; i < input_1->size; i++) {
    	    vectorWrite->value[i] = input_2->x[i] - input_3->x[i];
		}
    	vectorWrite->generation = 0;
    // предусмотреть количество использований
    	double rand = g_rand_double(hrand);

    	if (rand > control->use_archive_prob) { /* setting this probability to 0 will switch the archive off */
        	archiveVector = vectorWrite;
        	individ->useWriteVector = TRUE;
    	} else {
        	archiveVector = vectorRead;
    	}
	}

    switch (control->strategy) {
    case Simple:
        i = start_index;

        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = input_1->x[i] + control->f[i] * input_2->x[i];
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->moves++;
        break;
    case DE_3_bin:  // !

        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = input_1->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
		break;
    case DE_4_bin:

        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = input_1->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] ) + control->f[i] * ( input_4->x[i] - input_5->x[i] );
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
        input_4->grads++;
        input_5->grads++;
		break;
    case DE_3_bin_A:  // !

        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = input_1->x[i] + control->f[i] * archiveVector->value[i];
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
        break;

	case DE_3_exp:  // !
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            individ->x[i] = input_1->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_rand: // !
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = input_4->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_rand: // !
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            individ->x[i] = input_4->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_self:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->p[i] || L == 0 ) {
                individ->x[i] = individ->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_self:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            individ->x[i] = individ->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_T:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] ) {
                individ->x[i] = input_1->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );

            } else if ( u < 1 - control->p[i] ) {
                phi = input_1->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - input_1->cost ) * ( input_1->x[i] - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( input_1->cost - input_3->cost ) * ( input_3->x[i] - input_1->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( input_1->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_rand_T:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] ) {
                individ->x[i] = input_4->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            } else if ( u < 1 - control->p[i] ) {
                phi = input_4->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - input_4->cost ) * ( input_4->x[i] - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( input_4->cost - input_3->cost ) * ( input_3->x[i] - input_4->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( input_4->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_rand_phi:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] || L == 0 ) {
                phi = input_4->cost + input_2->cost + input_3->cost;
                individ->x[i] = control->f[i] * ( input_2->cost - input_4->cost ) * ( input_4->x[i] - input_2->x[i] );
                individ->x[i] += control->f[i] * ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += control->f[i] * ( input_4->cost - input_3->cost ) * ( input_3->x[i] - input_4->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( input_4->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_self_phi:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] || L == 0 ) {
                phi = individ->cost + input_2->cost + input_3->cost;
                individ->x[i] = control->f[i] * ( input_2->cost - individ->cost ) * ( individ->x[i] - input_2->x[i] );
                individ->x[i] += control->f[i] * ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += control->f[i] * ( individ->cost - input_3->cost ) * ( input_3->x[i] - individ->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( individ->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_self_T:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] ) {
                individ->x[i] = individ->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
            } else if ( u < 1 - control->p[i] ) {
                double indx = individ->x[i];
                phi = individ->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - individ->cost ) * ( indx - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( individ->cost - input_3->cost ) * ( input_3->x[i] - indx );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( indx + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_T:
        i = start_index;
        flag = 1;
        for ( L = 0; L < end_index; L++ ) {
            if ( flag == 1 ) {
                individ->x[i] = input_1->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 0;
            } else if ( flag == 0 ) {
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 2;
                phi = input_1->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - input_1->cost ) * ( input_1->x[i] - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( input_1->cost - input_3->cost ) * ( input_3->x[i] - input_1->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( input_1->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_1->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_rand_phi:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            phi = input_4->cost + input_2->cost + input_3->cost;
            individ->x[i] = control->f[i] * ( input_2->cost - input_4->cost ) * ( input_4->x[i] - input_2->x[i] );
            individ->x[i] += control->f[i] * ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
            individ->x[i] += control->f[i] * ( input_4->cost - input_3->cost ) * ( input_3->x[i] - input_4->x[i] );
            phi = ( phi <  0 ) ? -phi : phi;
            individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
            individ->x[i] += ( input_4->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_self_phi:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            phi = individ->cost + input_2->cost + input_3->cost;
            individ->x[i] = control->f[i] * ( input_2->cost - individ->cost ) * ( individ->x[i] - input_2->x[i] );
            individ->x[i] += control->f[i] * ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
            individ->x[i] += control->f[i] * ( individ->cost - input_3->cost ) * ( input_3->x[i] - individ->x[i] );
            phi = ( phi <  0 ) ? -phi : phi;
            individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
            individ->x[i] += ( individ->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_rand_T:
        i = start_index;
        flag = 1;
        for ( L = 0; L < end_index; L++ ) {
            if ( flag == 1 ) {
                individ->x[i] = input_4->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 0;
            } else if ( flag == 0 ) {
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 2;
                phi = input_4->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - input_4->cost ) * ( input_4->x[i] - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( input_4->cost - input_3->cost ) * ( input_3->x[i] - input_4->x[i] );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( input_4->x[i] + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_4->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_self_T:
        i = start_index;
        flag = 1;
        for ( L = 0; L < end_index; L++ ) {
            if ( flag == 1 ) {
                individ->x[i] = individ->x[i] + control->f[i] * ( input_2->x[i] - input_3->x[i] );
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 0;
            } else if ( flag == 0 ) {
                if ( g_rand_double(hrand) > control->p[i] )
                    flag = 2;
                double indx = individ->x[i];
                phi = individ->cost + input_2->cost + input_3->cost;
                individ->x[i] = ( input_2->cost - individ->cost ) * ( indx - input_2->x[i] );
                individ->x[i] += ( input_3->cost - input_2->cost ) * ( input_2->x[i] - input_3->x[i] );
                individ->x[i] += ( individ->cost - input_3->cost ) * ( input_3->x[i] - indx );
                phi = ( phi <  0 ) ? -phi : phi;
                individ->x[i] = ( phi ==  0 ) ? individ->x[i] : individ->x[i]/phi;
                individ->x[i] += ( indx + input_2->x[i] + input_3->x[i] ) / 3.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_rand_D:
        i = start_index;
        phi = input_2->cost + input_3->cost;
        phi = ( phi <  0 ) ? -phi : phi;
        alpha = input_3->cost - input_2->cost;
        alpha = ( alpha <  0 ) ? -alpha : alpha;
        phi = phi / alpha;
        alpha = (input_2->cost < input_3->cost) ? 1 : -1;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] || L == 0 ) {
                individ->x[i] = ( input_2->x[i] + input_3->x[i] ) / 2.0 - control->f[i] * alpha * ( input_3->x[i] - input_2->x[i] ) / 2.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_2->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_bin_self_D:
        i = start_index;
        phi = individ->cost + input_3->cost;
        phi = ( phi <  0 ) ? -phi : phi;
        alpha = input_3->cost - individ->cost;
        alpha = ( alpha <  0 ) ? -alpha : alpha;
        phi = phi / alpha;
        alpha = (individ->cost < input_3->cost) ? 1 : -1;
        for ( L = 0; L < end_index; L++ ) {
            u = g_rand_double(hrand);
            if ( u < control->p[i] || L == 0 ) {
                individ->x[i] = ( individ->x[i] + input_3->x[i] ) / 2.0 - control->f[i] * alpha * ( input_3->x[i] - individ->x[i] ) / 2.0;
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_3->grads++;
        break;
    case DE_3_exp_rand_D:
        i = start_index;
        phi = input_2->cost + input_3->cost;
        phi = ( phi <  0 ) ? -phi : phi;
        alpha = input_3->cost - input_2->cost;
        alpha = ( alpha <  0 ) ? -alpha : alpha;
        phi = phi / alpha;
        alpha = (input_2->cost < input_3->cost) ? 1 : -1;
        for ( L = 0; L < end_index; L++ ) {
            individ->x[i] = ( input_2->x[i] + input_3->x[i] ) / 2.0 - control->f[i] * alpha * ( input_3->x[i] - input_2->x[i] ) / 2.0;
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        input_2->moves++;
        input_2->grads++;
        input_3->grads++;
        break;
    case DE_3_exp_self_D:
        i = start_index;
        phi = individ->cost + input_3->cost;
        phi = ( phi <  0 ) ? -phi : phi;
        alpha = input_3->cost - individ->cost;
        alpha = ( alpha <  0 ) ? -alpha : alpha;
        phi = phi / alpha;
        alpha = (individ->cost < input_3->cost) ? 1 : -1;
        for ( L = 0; L < end_index; L++ ) {
            individ->x[i] = ( individ->x[i] + input_3->x[i] ) / 2.0 - control->f[i] * alpha * ( input_3->x[i] - individ->x[i] ) / 2.0;
            if ( g_rand_double(hrand) > control->p[i] )
                break;
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_2->grads++;
        input_3->grads++;
        break;

    case DE_3_triangular:   // !
        for(int i = 0; i < TRIANGULAR - 1; i++){
            for(int j = i + 1; j < TRIANGULAR; j++){
                if(tmp[i]->cost > tmp[j]->cost){
                    DpIndivid* tmp_ind = tmp[i];
                    tmp[i] = tmp[j];
                    tmp[j] = tmp_ind;
                }
            }
        }

        dp_triangular_rec(individ->x, tmp[0]->x, tmp[1]->x, tmp[2]->x,\
                individ->F, individ->p, individ->size);

        dp_triangular_rec(individ->F, tmp[0]->F, tmp[1]->F, tmp[2]->F,\
                individ->F, individ->p, TRIANGULAR);

        dp_triangular_rec(individ->p, tmp[0]->p, tmp[1]->p, tmp[2]->p,\
                individ->F, individ->p, TRIANGULAR);

        qsort(individ->p, TRIANGULAR, sizeof (double),\
              (int(*)(const void *, const void *)) cmp_double);

        for(int i = 0; i < TRIANGULAR; i++){
            individ->F[i] = (control->f_inf + control->f_supp)/ 2 + ((control->f_supp - control->f_inf)/ 2) * tanh(individ->F[i]);
        }

        individ->moves++;
        input_1->grads++;
        input_2->grads++;
        input_3->grads++;
        break;

    case DE_3_pBest:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->mod_CR[control->current_individ] || L == 0 ) {
                Fw = (FW_PARAM + FEs * (control->a_param - FW_PARAM) / MaxFEs) * control->mod_F[control->current_individ];
                //input4 = xpbest
                individ->x[i] = individ->x[i] + Fw * (input_4->x[i] - individ->x[i]) + control->mod_F[control->current_individ] * (input_1->x[i] - input_2->x[i]);
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_1->grads++;
        input_2->grads++;
        break;

    case DE_3_Amean:
        i = start_index;
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->mod_CR[control->current_individ] || L == 0 ) {
                Fw = (FW_PARAM + FEs * (control->a_param - FW_PARAM) / MaxFEs) * control->mod_F[control->current_individ];
                //xamean - vectorRead
                individ->x[i] = individ->x[i] + Fw * (vectorRead->value[i] - individ->x[i]) + control->mod_F[control->current_individ] * (input_1->x[i] - input_2->x[i]);
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_1->grads++;
        input_2->grads++;
        break;

    case DE_3_dual_adaptive:
        i = start_index;
        SPg = 1 / (1 + exp(1 - pow(FEs / MaxFEs, 2)));
        for ( L = 0; L < end_index; L++ ) {
            if ( g_rand_double(hrand) < control->mod_CR[control->current_individ] || L == 0 ) {
                Fw = (FW_PARAM + FEs * (control->a_param - FW_PARAM) / MaxFEs) * control->mod_F[control->current_individ];
                if (g_rand_double(hrand) < SPg) {
                    //xpbest - input_4
                    individ->x[i] = individ->x[i] + Fw * (input_4->x[i] - individ->x[i]) + control->mod_F[control->current_individ] * (input_1->x[i] - input_2->x[i]);
                } else {
                    //xamean - vectorRead
                    individ->x[i] = individ->x[i] + Fw * (vectorRead->value[i] - individ->x[i]) + control->mod_F[control->current_individ] * (input_1->x[i] - input_2->x[i]);
                }
            }
            i++;
            if ( i > end_index - 1 ) {
                i = 0;
            }
        }
        individ->moves++;
        input_1->grads++;
        input_2->grads++;
        break;
    default:
        g_error("dp_individ_recombination: unknown type %d", control->strategy);
    }
}

void dp_individ_recombination_ca(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index)
{
    int i, b;
    double cost_max;
    cost_max = input_1->cost;
    b = end_index;
    double alpha, beta, c1, c2, d, u;
    //	u = g_rand_double(hrand);
    if (individ->cost_ind < input_2->cost_ind) {
        alpha = 1;
        beta = (double)(input_2->cost_ind - individ->cost_ind) / (b - 1);
    } else {
        alpha = -1;
        beta = (double)(individ->cost_ind - input_2->cost_ind) / (b - 1);
    }
    for ( i = 0; i < individ->size; i++ ) {
        d = (input_2->x[i] - individ->x[i]) / 2.0;
        c1 = individ->x[i] - d * (1 + alpha * beta);
        c2 = individ->x[i] + d * (1 - alpha * beta);
        individ->x[i] = c1 + (c2 - c1) * g_rand_double(hrand);
        //		individ->x[i] = c1 + (c2 - c1) * u;
    }
    input_2->moves++;
}

void dp_individ_recombination_ac(DpRecombinationControl *control, GRand*hrand, DpIndivid*individ,  DpIndivid*input_1,  DpIndivid*input_2,  DpIndivid*input_3,  DpIndivid*input_4, int start_index, int end_index)
{
    int i;
    double d;
    for ( i = 0; i < individ->size; i++ ) {
        d = (input_4->x[i] - individ->x[i]) * 2.0;
        individ->x[i] += d * g_rand_double(hrand);
    }
    input_4->moves++;
}

DpRecombinationControl*dp_recombination_control_init(GRand*hrand,  DpPopulation*pop, GKeyFile*gkf, gchar*groupname)
{
    DpRecombinationControl*rc;
    int i;
    rc = (DpRecombinationControl*)malloc(sizeof(DpRecombinationControl));
    GError *gerror = NULL;
    gchar*str,**strlist;
    int retval = 0;
    rc->strategy = DE_3_bin_rand;
    double weight = 0.3;
    double prob = 0.3;
    rc->gamma = 0.9;
    rc->max_iters = 0;
    if ( ( str = g_key_file_get_string(gkf, groupname, "recombination_weight", &gerror) ) != NULL ) {
        weight = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "recombination_prob", &gerror) ) != NULL ) {
        prob = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "recombination_gamma", &gerror) ) != NULL ) {
        rc->gamma = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "a_param", &gerror) ) != NULL ) {
        rc->a_param = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "absolute_iter", &gerror) ) != NULL ) {
        rc->max_iters = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    if ( ( str = g_key_file_get_string(gkf, groupname, "recombination_strategy", &gerror) ) != NULL ) {
        if ( !g_strcmp0(str, "simple") ) {
            rc->strategy = Simple;
        } else if ( !g_strcmp0(str, "de_3_bin") ) {
            rc->strategy = DE_3_bin;
        } else if ( !g_strcmp0(str, "de_3_exp") ) {
            rc->strategy = DE_3_exp;
        } else if ( !g_strcmp0(str, "de_3_exp_T") ) {
            rc->strategy = DE_3_exp_T;
        } else if ( !g_strcmp0(str, "de_3_bin_T") ) {
            rc->strategy = DE_3_bin_T;
        } else if ( !g_strcmp0(str, "de_3_bin_A") ) {
            rc->strategy = DE_3_bin_A;
        } else if ( !g_strcmp0(str, "de_3_bin_rand") ) {
            rc->strategy = DE_3_bin_rand;
        } else if ( !g_strcmp0(str, "de_3_exp_rand") ) {
            rc->strategy = DE_3_exp_rand;
        } else if ( !g_strcmp0(str, "de_3_bin_self") ) {
            rc->strategy = DE_3_bin_self;
        } else if ( !g_strcmp0(str, "de_3_exp_self") ) {
            rc->strategy = DE_3_exp_self;
        } else if ( !g_strcmp0(str, "de_3_exp_rand_T") ) {
            rc->strategy = DE_3_exp_rand_T;
        } else if ( !g_strcmp0(str, "de_3_bin_rand_T") ) {
            rc->strategy = DE_3_bin_rand_T;
        } else if ( !g_strcmp0(str, "de_3_exp_rand_phi") ) {
            rc->strategy = DE_3_exp_rand_phi;
        } else if ( !g_strcmp0(str, "de_3_bin_rand_phi") ) {
            rc->strategy = DE_3_bin_rand_phi;
        } else if ( !g_strcmp0(str, "de_3_exp_self_phi") ) {
            rc->strategy = DE_3_exp_self_phi;
        } else if ( !g_strcmp0(str, "de_3_bin_self_phi") ) {
            rc->strategy = DE_3_bin_self_phi;
        } else if ( !g_strcmp0(str, "de_3_exp_rand_D") ) {
            rc->strategy = DE_3_exp_rand_D;
        } else if ( !g_strcmp0(str, "de_3_bin_rand_D") ) {
            rc->strategy = DE_3_bin_rand_D;
        } else if ( !g_strcmp0(str, "de_3_exp_self_D") ) {
            rc->strategy = DE_3_exp_self_D;
        } else if ( !g_strcmp0(str, "de_3_bin_self_D") ) {
            rc->strategy = DE_3_bin_self_D;
        } else if ( !g_strcmp0(str, "de_3_exp_self_T") ) {
            rc->strategy = DE_3_exp_self_T;
        } else if ( !g_strcmp0(str, "de_3_bin_self_T") ) {
            rc->strategy = DE_3_bin_self_T;
        } else if( !g_strcmp0(str, "de_3_triangular") ) {
            rc->strategy = DE_3_triangular;
        } else if( !g_strcmp0(str, "de_3_pbest") ) {
            rc->strategy = DE_3_pBest;
        } else if( !g_strcmp0(str, "de_3_amean") ) {
            rc->strategy = DE_3_Amean;
        } else if( !g_strcmp0(str, "de_3_dual_adaptive") ) {
            rc->strategy = DE_3_dual_adaptive;
        }
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    rc->size = pop->ind_size;
    rc->pop_size = pop->cur_size;
    rc->FEs = rc->pop_size;
    rc->MaxFEs = rc->pop_size * rc->max_iters;
    rc->toggle = g_rand_double(hrand) > 0.5 ? 1:0;
    rc->p_inf = G_MINDOUBLE;
    rc->p_supp = 1.0;
    rc->f_inf = 1.0 / sqrt((double)pop->ind_size);
    rc->f_supp = 2.0;
    rc->iter = pop->iter;
    if (rc->strategy == DE_3_pBest || rc->strategy == DE_3_Amean || rc->strategy == DE_3_dual_adaptive)
    {
        rc->mod_CR = (double*)calloc(pop->size, sizeof(double));
        rc->mod_F = (double*)calloc(pop->size, sizeof(double));
        rc->Mf = (double*)calloc(pop->size, sizeof(double));
        rc->Mcr = (double*)calloc(pop->size, sizeof(double));
        for (int i = 0; i < pop->size; i++) {
            rc->Mf[i] = 0.5;
            rc->Mcr[i] = 0.5;
        }
        rc->Sf = (double*)calloc(pop->size, sizeof(double));
        rc->Sf_size = 0;
        rc->Scr = (double*)calloc(pop->size, sizeof(double));
        rc->Scr_size = 0;
    }
    rc->f = (double*)calloc(rc->size, sizeof(double));
    rc->p = (double*)calloc(rc->size, sizeof(double));
    rc->c = (double*)calloc(rc->size, sizeof(double));
    rc->v = (double*)calloc(rc->size, sizeof(double));
    if (rc->strategy == DE_3_pBest || rc->strategy == DE_3_Amean || rc->strategy == DE_3_dual_adaptive) {
        dp_recombination_control_update_modified(hrand, rc, pop);
    }
    if ( weight == 0 && prob == 0 && rc->gamma > 0 ) {
        rc->adjust = 1;
        for ( i = 0; i < rc->size; i++ ) {
            rc->f[i] = g_rand_double_range(hrand, rc->f_inf, rc->f_supp);
            rc->p[i] = g_rand_double_range(hrand, rc->p_inf, rc->p_supp);
            rc->v[i] = pop->variance[i];
        }
    } else {
        rc->adjust = 0;
        for ( i = 0; i < rc->size; i++ ) {
            rc->f[i] = weight;
            rc->p[i] = prob;
            rc->v[i] = pop->variance[i];
        }
    }
	rc->use_archive_prob = 0;
    if ( ( str = g_key_file_get_string(gkf, groupname, "use_archive_prob", &gerror) ) != NULL ) {
        rc->use_archive_prob = g_strtod( str , NULL);
        g_free(str);
    } else {
        g_debug ("%s", gerror->message );
        g_clear_error (&gerror);
    }
    return rc;
}

void dp_recombination_control_save(FILE*fp, DpRecombinationControl*rc)
{
    int i;
    for ( i = 0; i < rc->size; i++ ) {
        fprintf(fp, "%13.9f ", rc->f[i]);
    }
    fprintf(fp, "\n");
    for ( i = 0; i < rc->size; i++ ) {
        fprintf(fp, "%13.9f ", rc->p[i]);
    }
    fprintf(fp, "\n");
    for ( i = 0; i < rc->size; i++ ) {
        fprintf(fp, "%13.9f ", rc->v[i]);
    }
    fprintf(fp, "\n");
    for ( i = 0; i < rc->size; i++ ) {
        fprintf(fp, "%13.9f ", rc->c[i]);
    }
}

void dp_recombination_control_load(FILE*fp, DpRecombinationControl*rc)
{
    int i;
    for ( i = 0; i < rc->size; i++ ) {
        fscanf(fp, "%lf", &(rc->f[i]));
    }
    for ( i = 0; i < rc->size; i++ ) {
        fscanf(fp, "%lf", &(rc->p[i]));
    }
    for ( i = 0; i < rc->size; i++ ) {
        fscanf(fp, "%lf", &(rc->v[i]));
    }
    for ( i = 0; i < rc->size; i++ ) {
        fscanf(fp, "%lf", &(rc->c[i]));
    }
}

void dp_recombination_control_renew(DpRecombinationControl*rc, GRand*hrand, DpPopulation*pop)
{
    int i;
    if ( rc->adjust == 1 ) {
        for ( i = 0; i < rc->size; i++ ) {
            rc->f[i] = g_rand_double_range(hrand, rc->f_inf, rc->f_supp);
            rc->p[i] = g_rand_double_range(hrand, rc->p_inf, rc->p_supp);
            rc->v[i] = pop->variance[i];
        }
    }
}

void dp_recombination_control_update(DpRecombinationControl*rc, GRand*hrand, DpPopulation*pop, int start_index, int end_index)
{
    int i, j;
    double flag;
    if ( rc->adjust == 0 )
        return;
    for ( i = 0; i < rc->size; i++ ) {
        if ( pop->variance[i] > 0 )
            rc->c[i] = rc->gamma * rc->v[i] / pop->variance[i];
        else
            rc->c[i] = 0.0;
        rc->v[i] = pop->variance[i];
    }
    if ( rc->toggle == 1 ) {
        rc->toggle = 0;
        for ( j = 0; j < rc->size; j++ ) {
            flag = rc->pop_size * rc->f[j] * rc->f[j] - 1.0;
            rc->p[j] = ( rc->c[j] < 1 ) ? rc->p_inf : -flag + sqrt ( flag * flag - rc->pop_size * ( 1.0 - rc->c[j] ) );
            if ( ( rc->p[j] >= rc->p_supp ) || ( rc->p[j] <= rc->p_inf ) )
                rc->p[j] = g_rand_double_range(hrand, rc->p_inf, rc->p_supp);
        }
    } else {
        rc->toggle = 1;
        for ( j = 0; j < rc->size; j++ ) {
            flag = rc->pop_size * ( rc->c[j] - 1.0 ) + rc->p[j] * ( rc->p_supp + 1 - rc->p[j] );
            rc->f[j] = ( flag < 0 ) ? rc->f_inf : sqrt ( flag / ( 2.0 * rc->pop_size * rc->p[j] ) );
            if ( ( rc->f[j] >= rc->f_supp ) || ( rc->f[j] <= rc->f_inf ) )
                rc->f[j] = g_rand_double_range(hrand, rc->f_inf, rc->f_supp);
        }
    }
}

double dp_rand_normal(GRand*hrand, double mean, double stddev) {
    static double V2, fac;
	static int phase = 0;
	double S, Z, U1, U2, V1;
	if ( phase == 1 ) {
		Z = V2 * fac;
	} else {
		do {
			U1 = g_rand_double(hrand);
			U2 = g_rand_double(hrand);
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while ( S >= 1 || S == 0.0 );
		fac = sqrt (-2 * log(S) / S);
		Z = V1 * fac;
	}
	phase = 1 - phase;
	Z = mean + stddev*Z;
	return Z;
}

double dp_rand_Cauchy(GRand*hrand, double loc, double scale) {
    double u;
    do {
        u = g_rand_double(hrand);
    } while (u == 0.0 || u == 1.0);
    return loc + scale * tan(acos(-1.0) * (u - 0.5));
}

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void dp_recombination_control_update_modified(GRand*hrand, DpRecombinationControl*rc, DpPopulation*pop) {
    for (int i = 0; i < pop->size; i++) {
        if (rc->Mcr[i] < 0) {
            rc->mod_CR[i] = 0;
        } else {
            rc->mod_CR[i] = dp_rand_normal(hrand, rc->Mcr[i], 0.1);
        }

        rc->mod_CR[i] = MAX(rc->mod_CR[i], 0);
        rc->mod_CR[i] = MIN(rc->mod_CR[i], 1);
        rc->mod_F[i] = dp_rand_Cauchy(hrand, rc->Mf[i], 0.1);
    }
}

void dp_recombination_control_update_Mcr_and_Mf(DpRecombinationControl*rc, DpPopulation*pop, DpPopulation*trial) {
    double wm, upper_sum_wm, lower_sum_wm;
    double x_squared, u_squared;
    double upper_sum_CR, upper_sum_F, lower_sum_CR, lower_sum_F;
    if (rc->Scr_size == 0) {
        return;
    }

    for (int i = 0; i < pop->cur_size; i++) {
        upper_sum_CR = 0;
        upper_sum_F = 0;
        lower_sum_CR = 0;
        lower_sum_F = 0;
        for (int m = 0; m < rc->Scr_size; m++) {
            wm = 0;
            upper_sum_wm = 0;
            lower_sum_wm = 0;
            for (int j = 0; j < rc->size; j++) {
                upper_sum_wm += pop->individ[m]->x[j] * trial->individ[m]->x[j];
            }
            for (int k = 0; k < rc->Scr_size; k++) {
                x_squared = 0;
                u_squared = 0;
                for (int in_j; in_j < rc->size; in_j++) {
                    x_squared += pow(pop->individ[k]->x[in_j], 2);
                    u_squared += pow(trial->individ[k]->x[in_j], 2);
                }
                lower_sum_wm += sqrt(x_squared) * sqrt(u_squared);
            }
            wm = upper_sum_wm / lower_sum_wm;
            upper_sum_CR += wm * pow(rc->Scr[m], 2);
            lower_sum_CR += wm * rc->Scr[m];
            upper_sum_F += wm * pow(rc->Sf[m], 2);
            lower_sum_F += wm * rc->Sf[m];
        }
        rc->Mcr[i] = upper_sum_CR / lower_sum_CR;
        rc->Mf[i] = upper_sum_F / lower_sum_F;
    }
}

void dp_recombination_control_renew_Scr_and_Sf(DpRecombinationControl*rc) {
    for (int i=0; i < rc->Scr_size; i++) {
        rc->Scr[i] = 0;
    }
    rc->Scr_size = 0;
    for (int j=0; j < rc->Sf_size; j++) {
        rc->Sf[j] = 0;
    }
    rc->Sf_size = 0;
}