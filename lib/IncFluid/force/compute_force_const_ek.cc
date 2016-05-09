/* TarangMPI-4.0
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

/*! \file  compute_force_const_ek_hk.cc
 * 
 * @brief Compute force when ek and hk are held constant.
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates.
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 * 
 * @bug sk=1, -1 needs to be handled separately. Maximum helicity state requires 
 *		V field to be maximum helical, which is zero probability state.
 */ 


#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_const_energy()
{
	if (force_switch == 1)  {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_level;
		
		
		static int nf; 
		static DP energy_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_level = (*force_double_para)(3);
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_per_mode = energy_level / nf;	
			
			
			kx_max = (int) ceil(outer_radius/kfactor[1]);
			
			if (N[2] > 1)
				ky_max = (int) ceil(outer_radius/kfactor[2]);
			else
				ky_max = 0;
			
			if (N[3] >= 2)
				kz_max = (int) ceil(outer_radius/kfactor[3]);
			else	
				kz_max = 0;
			
			if (basis_type == "FOUR")
				kx_min = -kx_max;
			
			else if (basis_type == "SCFT")
				kx_min = 0;
			
			is_force_field_para_read = 1;	
			
		}
		
		*Force1 = 0.0; 
		*Force2 = 0.0; 
		*Force3 = 0.0;
	
		int lx, ly, lz;
		DP kkmag, alpha_k, beta_k, modal_energy;
		
			
		for (int kx = kx_min; kx <= kx_max; kx++)
			for (int ky = -ky_max; ky <= ky_max; ky++)  
				for (int kz = 0; kz <= kz_max; kz++)
				{
					lx = Get_lx(basis_type, kx, N);
					ly = Get_ly3D(basis_type, ky, N);
					lz = kz;
					
					if ( (lx >= 0) && (lx < local_N1) )  {
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
				
							if (energy_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = sqrt(energy_per_mode / modal_energy); 
									beta_k = 0.0;
												   
									Const_vector_field_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}			   
							}
						}			
					}
				}						
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force();
		
		Satisfy_reality_condition_force_field();
	}	

}


//*********************************************************************************************
//
// scalar
//

void IncFluid::Compute_force_const_energy(IncSF& T)
{
	if (T.force_switch == 0)
		Compute_force_const_energy();
	
	else if (T.force_switch == 1) {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_level;
		static DP energy_scalar_level;
		
		static int nf; 
		static DP energy_per_mode;
		static DP energy_scalar_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_level = (*force_double_para)(3);
			energy_scalar_level = (*force_double_para)(4);
			
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_per_mode = energy_level / nf;	
			energy_scalar_per_mode = energy_scalar_level / nf;
			

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
		DP kkmag, alpha_k, beta_k, alpha_k_scalar, sk, modal_energy;
		
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
						
						if ((kkmag > inner_radius) && (kkmag <= outer_radius))  {
							if (energy_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = sqrt(energy_per_mode / modal_energy); 
									beta_k = 0.0;
									
									Const_vector_field_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}
							}

							if (energy_scalar_per_mode > MYEPS) {
								
								modal_energy = T.CS_Modal_energy(lx, ly, lz); 
								if (modal_energy > MYEPS) {
									alpha_k_scalar = sqrt(energy_scalar_per_mode / modal_energy);
									
									T.Const_scalar_field_alpha(kx, ky, kz, alpha_k_scalar);
								}
							}
						}			
					}
				}
				
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


void IncFluid::Compute_force_const_energy(IncVF& W)
{
	if (W.force_switch == 0)
		Compute_force_const_energy();
	
	else if (W.force_switch == 1) {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_level;
		static DP energyW_level;
		
		
		static int nf; 
		static DP energy_per_mode;
		static DP energyW_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_level = (*force_double_para)(3);
			energyW_level = (*force_double_para)(4);	
		
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_per_mode = energy_level / nf;	
			energyW_per_mode = energyW_level / nf;
			
			
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
							
							if (energy_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = sqrt(energy_per_mode / modal_energy); 
									beta_k = 0.0;
									
									Const_vector_field_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}
							}
							
							// For W 
							
							if (energyW_per_mode > MYEPS) {
								
								modal_energy = W.CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k_W = sqrt(energyW_per_mode / modal_energy); 
									beta_k_W = 0.0;
									
									W.Const_vector_field_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);	
								}
							}
							
						}		
					}					
				}
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
			Dealias_force(W);
		
		Satisfy_reality_condition_force_field(W);
	}				
}


//*********************************************************************************************
//
// Vector + scalar
//

void IncFluid::Compute_force_const_energy(IncVF& W, IncSF& T)
{
	if ((W.force_switch == 0) && (T.force_switch == 0))
		Compute_force_const_energy();
	
	else if ((W.force_switch == 1) && (T.force_switch == 0))
		Compute_force_const_energy(W);
	
	else if ((W.force_switch == 1) && (T.force_switch == 1)) {
		static DP inner_radius;
		static DP outer_radius;
		static DP energy_level;
		static DP h_by_k_E;				// Hk/(k*e)
		static DP energyW_level;
		static DP h_by_k_E_W;	
		static DP energy_scalar_level;
		
		
		static int nf; 
		static DP energy_per_mode;
		static DP energyW_per_mode;
		static DP energy_scalar_per_mode;
		
		static int kx_max, ky_max, kz_max, kx_min;
		
		if (is_force_field_para_read == 0)
		{
			inner_radius = (*force_double_para)(1);
			outer_radius = (*force_double_para)(2);
			energy_level = (*force_double_para)(3);
			energyW_level = (*force_double_para)(4);
			energy_scalar_level = (*force_double_para)(5);
			
			
			nf = Get_number_modes_in_shell(basis_type, N, inner_radius, outer_radius, kfactor);
			
			energy_per_mode = energy_level / nf;	
			energyW_per_mode = energyW_level / nf;
			energy_scalar_per_mode = energy_scalar_level / nf;
			
			if (N[1] > 1)
				kx_max = (int) ceil(outer_radius/kfactor[1]);
			else 
				kx_max = 0;
			
			
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
					
					if ( (lx >= 0) && (lx < local_N1) ) {
						
						kkmag = Kmagnitude(basis_type, lx, ly, lz, N, kfactor);
					
						if ((kkmag > inner_radius) && (kkmag <= outer_radius)) {
							
							// For U
							if (energy_per_mode > MYEPS) {
								
								modal_energy = CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k = sqrt(energy_per_mode / modal_energy); 
									beta_k = 0.0;
									
									Const_vector_field_alpha_beta(kx, ky, kz, alpha_k, beta_k);	
								}
							}
							
							// For W 
							
							if (energyW_per_mode > MYEPS) {
								
								modal_energy = W.CV_Modal_energy(lx, ly, lz);
								if (modal_energy > MYEPS) {
									alpha_k_W = sqrt(energyW_per_mode / modal_energy); 
									beta_k_W = 0.0;
									
									W.Const_vector_field_alpha_beta(kx, ky, kz, alpha_k_W, beta_k_W);	
								}
							}

							// scalar
							
							if (energy_scalar_per_mode > MYEPS) {
								modal_energy = T.CS_Modal_energy(lx, ly, lz); 
								if (modal_energy > MYEPS) {
									alpha_k_scalar = sqrt(energy_scalar_per_mode / modal_energy);
									
									T.Const_scalar_field_alpha(kx, ky, kz, alpha_k_scalar);
								}
							}
							
						} // of if ((kmag > ..) 	
					}	// of if (lx ..)				
				} // of for ()
				
		if ( (alias_switch == "DEALIAS") 
				&& (Is_alias_array(basis_type, N, *V1, outer_radius, kfactor) == 1) )
				Dealias_force(W, T);
		
		Satisfy_reality_condition_force_field(W, T);
	}	
}

//*****************************  End of compute_force_const_ek_hk.cc **************************


