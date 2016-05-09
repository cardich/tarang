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

/*! \file  compute_force_ek_hk_supply.cc
 * 
 * @brief Compute force when ek and hk supply rate is given
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug sk=1, -1 needs to be handled separately.
 */

#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_const_energy_supply()
{
	if (force_switch == 1)  {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_supply;
		
		static int nf; 
		static DP energy_supply_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_supply = (*force_double_para)(3);		
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_supply_per_mode = energy_supply / nf;	
			
			kx_max = (int) ceil(outer_radius/kfactor[1]);
			
			if (N[2] > 1)
				ky_max = (int) ceil(outer_radius/kfactor[2]);
			else
				ky_max = 0;		
			
			if (N[3] > 2)
				kz_max = (int) ceil(outer_radius/kfactor[3]);
			else	
				kz_max = 0;
			
			
			
			if (basis_type == "FOUR")
				kx_min = -kx_max;
			
			else if (basis_type == "SCFT")
				kx_min = 0;
			
			is_force_field_para_read = 1;	
			
		}
	
		int lx, ly, lz;
		DP kkmag, alpha_k, beta_k, modal_energy;
		
		*Force1 = 0.0; 
		*Force2 = 0.0; 
		*Force3 = 0.0;
					
			
		for (int kx = kx_min; kx <= kx_max; kx++)
			for (int ky = -ky_max; ky <= ky_max; ky++)  
				for (int kz = 0; kz <= kz_max; kz++)
				{
					lx = Get_lx(basis_type, kx, N);
					ly = Get_ly3D(basis_type, ky, N);
					lz = kz;
					
					if ( (lx >= 0) && (lx < local_N1) ) {
						
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
							
							if (energy_supply_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = energy_supply_per_mode/ (2*modal_energy);
									beta_k = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}	
							}
						}
					}	//  of if ((lx >= 0)...)						
				}		// of for
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1,	outer_radius, kfactor) == 1) )
			Dealias_force();
		
		Satisfy_reality_condition_force_field();
	}	
}

//*********************************************************************************************
//
//	Scalar
//

void IncFluid::Compute_force_const_energy_supply(IncSF& T)
{

	if (T.force_switch == 0)
		Compute_force_const_energy_supply();
	
	else if (T.force_switch == 1) {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_supply;
		static DP energy_supply_scalar;
		
		
		static int nf; 
		static DP energy_supply_per_mode;
		static DP energy_supply_scalar_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_supply = (*force_double_para)(3);
			energy_supply_scalar = (*force_double_para)(4);
			
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_supply_per_mode = energy_supply / nf;
			energy_supply_scalar_per_mode = energy_supply_scalar / nf;	
			

			kx_max = (int) ceil(outer_radius/kfactor[1]);
					
			if (N[2] > 1)
				ky_max = (int) ceil(outer_radius/kfactor[2]);
			else
				ky_max = 0;		
			
			if (N[3] > 2)
				kz_max = (int) ceil(outer_radius/kfactor[3]);
			else	
				kz_max = 0;
			
			if (basis_type == "FOUR")
				kx_min = -kx_max;
			
			else if (basis_type == "SCFT")
				kx_min = 0;
			
			is_force_field_para_read = 1;	
			
		}
		
		
		int lx, ly, lz;	
		DP kkmag, alpha_k, beta_k, alpha_k_scalar, modal_energy;
	
		if (force_switch == 1) {
			*Force1 = 0.0; 
			*Force2 = 0.0; 
			*Force3 = 0.0;
		}
		
		if (T.force_switch == 1) 
			*T.Force = 0.0;
							
			
		for (int kx = kx_min; kx <= kx_max; kx++)
			for (int ky = -ky_max; ky <= ky_max; ky++)  
				for (int kz = 0; kz <= kz_max; kz++)
				{
					lx = Get_lx(basis_type, kx, N);
					ly = Get_ly3D(basis_type, ky, N);
					lz = kz;
					
					if ( (lx >= 0) && (lx < local_N1) ) {
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
							
							if (energy_supply_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = energy_supply_per_mode/ (2*modal_energy);
									beta_k = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}	
							}
							
							
							if (energy_supply_scalar_per_mode > MYEPS) {
								
								modal_energy = T.CS_Modal_energy(lx, ly, lz); 
								if (modal_energy > MYEPS) {
									alpha_k_scalar = energy_supply_scalar_per_mode 
													/ (2*modal_energy); 
									
									T.Force_scalar_alpha(kx, ky, kz, alpha_k_scalar);
								}
							}
						}
					}   // of ( (lx >= 0)..)						
				}		// of for
			
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(T);
		
		Satisfy_reality_condition_force_field(T);
	}	
						
}

//*********************************************************************************************
// 
// Vector
//
void IncFluid::Compute_force_const_energy_supply(IncVF& W)
{

	if (W.force_switch == 0)
		Compute_force_const_energy_supply();
	
	else if (W.force_switch == 1) {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_supply;
		static DP energy_supply_W;
		
		static int nf; 
		static DP energy_supply_per_mode;
		static DP energy_supply_W_per_mode;	
			
		static int kx_max, ky_max, kz_max, kx_min;
			
		if (is_force_field_para_read == 0)
		{
			
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_supply = (*force_double_para)(3);
			energy_supply_W = (*force_double_para)(4);
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_supply_per_mode = energy_supply / nf;
			energy_supply_W_per_mode = energy_supply_W / nf;	
			
			
			kx_max = (int) ceil(outer_radius/kfactor[1]);
			
			if (N[2] > 1)
				ky_max = (int) ceil(outer_radius/kfactor[2]);
			else
				ky_max = 0;		
			
			if (N[3] > 2)
				kz_max = (int) ceil(outer_radius/kfactor[3]);
			else	
				kz_max = 0;
			

			if (basis_type == "FOUR")
				kx_min = -kx_max;
			
			else if (basis_type == "SCFT")
				kx_min = 0;
			
			
			is_force_field_para_read = 1;	
		
		}
		
		int lx, ly, lz;
		DP kkmag;
		DP alpha_k, beta_k, modal_energy;
		DP alpha_k_W, beta_k_W;
		
		if (force_switch == 1) {
			*Force1 = 0.0; 
			*Force2 = 0.0; 
			*Force3 = 0.0;
		}
		
		if (W.force_switch == 1) {
			*W.Force1 = 0.0; 
			*W.Force2 = 0.0; 
			*W.Force3 = 0.0;
		}
				
		for (int kx = kx_min; kx <= kx_max; kx++)
			for (int ky = -ky_max; ky <= ky_max; ky++)  
				for (int kz = 0; kz <= kz_max; kz++)
				{
					lx = Get_lx(basis_type, kx, N);
					ly = Get_ly3D(basis_type, ky, N);
					lz = kz;
					
					if ( (lx >= 0) && (lx < local_N1) ) {
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
							
							// for U
							if (energy_supply_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = energy_supply_per_mode/ (2*modal_energy);
									beta_k = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}
							}
										
							// for W
							if (energy_supply_W_per_mode > MYEPS) {
								
								modal_energy = W.CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k_W = energy_supply_W_per_mode/ (2*modal_energy);
									beta_k_W = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);	
								}
							}
						}  
					}	// of ( (lx >= 0)...)					
				}		// of for
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(W);
		
		Satisfy_reality_condition_force_field(W);
	}	
}


//*********************************************************************************************
//
// Vector + Scalar
//

void IncFluid::Compute_force_const_energy_supply(IncVF& W, IncSF& T)
{

	if ((W.force_switch == 0) && (T.force_switch == 0))
		Compute_force_const_energy_supply();
	
	else if ((W.force_switch == 1) && (T.force_switch == 0))
		Compute_force_const_energy_supply(W);
	
	else if ((W.force_switch == 1) && (T.force_switch == 1)) {
		
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_supply;
		static DP energy_supply_W;
		static DP energy_supply_scalar;
		
		static int nf; 
		static DP energy_supply_per_mode;
		static DP energy_supply_W_per_mode;
		static DP energy_supply_scalar_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_supply = (*force_double_para)(3);
			energy_supply_W = (*force_double_para)(4);
			energy_supply_scalar = (*force_double_para)(5);
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_supply_per_mode = energy_supply / nf;
			energy_supply_W_per_mode = energy_supply_W / nf;
			energy_supply_scalar_per_mode = energy_supply_scalar / nf;
			

			kx_max = (int) ceil(outer_radius/kfactor[1]);
			
			if (N[2] > 1)
				ky_max = (int) ceil(outer_radius/kfactor[2]);
			else
				ky_max = 0;
					
			if (N[3] > 2)
				kz_max = (int) ceil(outer_radius/kfactor[3]);
			else	
				kz_max = 0;
			
			if (basis_type == "FOUR")
				kx_min = -kx_max;
			
			else if (basis_type == "SCFT")
				kx_min = 0;
			
			
			is_force_field_para_read = 1;	
			
		}
		
	
		int lx, ly, lz;
		DP kkmag;
		DP alpha_k, beta_k, modal_energy;
		DP alpha_k_W, beta_k_W;
		DP alpha_k_scalar;
		
		if (force_switch == 1) {
			*Force1 = 0.0; 
			*Force2 = 0.0; 
			*Force3 = 0.0;
		}
		
		if (W.force_switch == 1) {
			*W.Force1 = 0.0; 
			*W.Force2 = 0.0; 
			*W.Force3 = 0.0;
		}
		
		if (T.force_switch == 1) 
			*T.Force = 0.0;
			
				
		for (int kx = kx_min; kx <= kx_max; kx++)
			for (int ky = -ky_max; ky <= ky_max; ky++)  
				for (int kz = 0; kz <= kz_max; kz++)
				{
					lx = Get_lx(basis_type, kx, N);
					ly = Get_ly3D(basis_type, ky, N);
					lz = kz;
					
					if ( (lx >= 0) && (lx < local_N1) ) 
					{
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
			
							// for U
							if (energy_supply_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = energy_supply_per_mode/ (2*modal_energy);
									beta_k = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}
							}
							
							// for W
							if (energy_supply_W_per_mode > MYEPS) {
								
								modal_energy = W.CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k_W = energy_supply_W_per_mode/ (2*modal_energy);
									beta_k_W = 0.0;
									
									Force_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);	
								}
							}
								
							
							if (energy_supply_scalar_per_mode > MYEPS) {

								modal_energy = T.CS_Modal_energy(lx, ly, lz); 
								if (modal_energy > MYEPS) {
									alpha_k_scalar = energy_supply_scalar_per_mode 
														/ (2*modal_energy); 
									
									T.Force_scalar_alpha(kx, ky, kz, alpha_k_scalar);
								}
							}
						}
					}	// of if ( (lx >= 0) ..)	
				}		// of for loop
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(W, T);
		
		Satisfy_reality_condition_force_field(W, T);
	}	
					
}


//***********************  End of compute_force_ek_hk_supply.cc *******************************


