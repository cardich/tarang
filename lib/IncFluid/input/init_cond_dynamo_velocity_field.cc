/* Tarang-4.0
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-4.0 .
 *
 * Tarang-4.0 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-4.0 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-4.0; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */


/*! \file  init_cond_TG.cc
 * 
 * @brief Initial conditions as Taylor Green flow (TG).
 *
 * @note The parameters are read from parameter file. 
 *
 *		Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x)sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "../IncFluid.h"


//*********************************************************************************************

//
// Vector
//

void IncFluid::Init_cond_dynamo_full_velocity_field(IncVF& W)
{
	Init_cond();  // Read the velocity field..
	
// Init_cond_double_para(1,2,3) = totalEb,	kkmin, kkmax
// Distribute totalEb among the modes in the shell
	
	DP totalEb = (*init_cond_double_para)(1);
	DP kkmin = (*init_cond_double_para)(2);
	DP kkmax = (*init_cond_double_para)(3);
	
	int total_no_modes = (4*M_PI/3)*(my_pow(kkmax,3) - my_pow(kkmin,3));
	DP ekW = totalEb/total_no_modes;   // Energy per mode
	
	DP kkmag, ampW, phase1W, phase2W, phase3W;
	int index, maxN3;
	
	if (N[3] > 2)
		maxN3 = N[3]/2;
	
	else	// 2D
		maxN3 = 0;
	
	for (int l1=0; l1<local_N1; l1++)
		for (int l2=0; l2<N[2]; l2++) 
			for (int l3=0; l3<=maxN3; l3++)  {
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if ((kkmag > kkmin) && (kkmag <= kkmax)) {
					index = (int) ceil(kkmag);
					
					ampW = sqrt(2*ekW);
					
					phase1W = 2*M_PI * SPECrand.random();
					phase2W = 2*M_PI * SPECrand.random();
					phase3W = 2*M_PI * SPECrand.random();
					
					W.Put_vector_amp_phase_comp_conj(l1, l2, l3, N, ampW, 
													 phase1W, phase2W, phase3W);
				}	
			}
	
	if (my_id == master_id) 
		(*W.V1)(0,0,0) = (*W.V2)(0,0,0) = (*W.V2)(0,0,0) = 0.0;
	
	
	if (N[3] == 2) {		
		(*W.V1)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V2)(Range::all(), Range::all(), 1) = 0.0;
		(*W.V3)(Range::all(), Range::all(), Range(0,1)) = 0.0;
	}
	
	if (N[2] == 1) 
		(*W.V2)(Range::all(), 0, Range::all()) = 0.0;
	
	if (alias_switch == "DEALIAS")		Dealias(W);
	
	Satisfy_reality_condition_field(W);
	
}



//******************************** End of Init_cond_TG.cc *************************************







