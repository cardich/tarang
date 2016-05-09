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

/*! \file  compute_force_LM.cc
 * 
 * @brief  Set up force using Liquid metal for under strong mean magnetic field.
 *
 * @note F =  (hat(B0).hat(K))^2 V
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_Liquid_metal()
{	
	
	DP B0x = (*force_double_para)(1);
	DP B0y = (*force_double_para)(2);
	DP B0z = (*force_double_para)(3);
	
	TinyVector <DP,3> B0;
	B0 = B0x, B0y, B0z;
	
	*Force1 = *V1;	
	*Force2 = *V2; 
	*Force3 = *V3;
	
	Array_mult_V0_khat_sqr(basis_type, N, *Force1, B0, kfactor);
	Array_mult_V0_khat_sqr(basis_type, N, *Force2, B0, kfactor);
	Array_mult_V0_khat_sqr(basis_type, N, *Force3, B0, kfactor);
	
	*Force1 = complx(-1,0) * (*Force1);
	*Force2 = complx(-1,0) * (*Force2);
	*Force3 = complx(-1,0) * (*Force3);
	
	if (alias_switch == "DEALIAS")   Dealias_force();

}	



void IncFluid::Compute_force_Liquid_metal_const_energy_supply()
{	
	
	DP B0x = (*force_double_para)(4);
	DP B0y = (*force_double_para)(5);
	DP B0z = (*force_double_para)(6);
	
	DP inner_radius = (*force_double_para)(1);
	DP outer_radius = (*force_double_para)(2);
	DP energy_supply = (*force_double_para)(3);
	
	// first feed const eps force;
	Compute_force_const_energy_supply();
	
	// LM force
	TinyVector <DP,3> B0;
	B0 = B0x, B0y, B0z;
	
	*VF_temp = *V1;	
	Array_mult_V0_khat_sqr(basis_type, N, *VF_temp, B0, kfactor);
	*Force1 = *Force1 + complx(-1,0) * (*VF_temp);
	
	
	*VF_temp = *V2;	
	Array_mult_V0_khat_sqr(basis_type, N, *VF_temp, B0, kfactor);
	*Force2 = *Force2 + complx(-1,0) * (*VF_temp);
	
	*VF_temp = *V3;	
	Array_mult_V0_khat_sqr(basis_type, N, *VF_temp, B0, kfactor);
	*Force3 = *Force3 + complx(-1,0) * (*VF_temp);
	
	if (alias_switch == "DEALIAS")   Dealias_force();
	
}	


//*********************************************************************************************

void IncFluid::Compute_force_Liquid_metal(IncSF& T)
{
	
	Compute_force_Liquid_metal();
	
	// For the velocity field
	if (globalvar_Pr_switch == "PRLARGE") {
		
		if (globalvar_RB_Uscaling == "USMALL") 		
			*Force1 = *Force1 + (globalvar_Ra*globalvar_Pr)*(*T.F);			// (u.grad)u-Ra*Pr*theta	
		
		else if (globalvar_RB_Uscaling  == "ULARGE") 		
			*Force1 = *Force1 + (*T.F);					// (u.grad)u-theta			
	}
	
	else if ((globalvar_Pr_switch == "PRSMALL") || (globalvar_Pr_switch == "PRZERO")) {
		
		if (globalvar_RB_Uscaling == "USMALL") 
			*Force1 = *Force1 + (globalvar_Ra)*(*T.F);				// (u.grad)u-theta
		
		else if (globalvar_RB_Uscaling == "ULARGE") 		
			*Force1 = *Force1 + (globalvar_Pr)*(*T.F);				// (u.grad)u-theta
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") {
		*Force1 = *Force1 + (globalvar_Ra)*(*T.F);
	}
	
	
	*Force2 = 0.0; 
	*Force3 = 0.0;
	
	
	// For the temperature field
	
	if (globalvar_Pr_switch == "PRZERO") 
		*T.Force = 0.0;  
	// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
	
	else if (globalvar_Pr_switch == "PRLARGE")	
		*T.Force =   globalvar_temperature_grad*(*V1);										
	// F(T) = globalvar_temperature_grad * ux(k)	
	// globalvar_temperature_grad = +1 for RB, and -1 for stratified flows
	
	else if (globalvar_Pr_switch == "PRSMALL") 	
		*T.Force =   complex<DP>(globalvar_temperature_grad/globalvar_Pr, 0)*(*V1); 			
	//F(T) =  globalvar_temperature_grad * ux(k)/Pr
	
	else if (globalvar_Pr_switch == "PRINFTY") {
		*T.Force =   globalvar_temperature_grad*(*V1);
	}
	
	
	if (alias_switch == "DEALIAS")		Dealias_force(T);	
}


//*********************************************************************************************

// Eqn 2D NS with additional forcing: -u/Rh

/** @brief Kolmogorov flow forcing
 * 
 * @note Force_x = amp*sin(k0 x) cos(k0 y) cos(k0 z) - V1/Rh
 *		Force_y = 0
 *		Force_z = -amp*cos(k0 x)sin(k0 y) cos(k0 z) -V3/Rh
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncFluid::Compute_force_Kolmogorov_flow()
{	
	int k0 = ((int) (*force_double_para)(1));
	DP force_amp = (*force_double_para)(2);
	DP Rh = (*force_double_para)(3);
	
	if (force_switch == 1) {
		
		// -u/Rh
		*Force1 = complx(-1/Rh,0) * (*V1);
		*Force2 = 0.0; 
		*Force3 = complx(-1/Rh,0) * (*V3);
		
		// Add the current forces
		if (basis_type == "FOUR") {
			
			int lx_k0 = Get_lx("FOUR", k0, N);
			int lx_minus_k0 = Get_lx("FOUR", -k0, N);
			
			if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) {
				(*Force1)(lx_k0, 0, k0)  = (*Force1)(lx_k0, 0, k0) - I*(force_amp/4);
				(*Force3)(lx_k0, 0, k0)  = (*Force3)(lx_k0, 0, k0) + I*(force_amp/4);
			}
			
			if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) {
				(*Force1)(lx_minus_k0, 0, k0) = (*Force1)(lx_minus_k0, 0, k0) + I*(force_amp/4);
				(*Force3)(lx_minus_k0, 0, k0) = (*Force3)(lx_minus_k0, 0, k0) + I*(force_amp/4);
			}
		}
		
		else if (basis_type == "SCFT") {
			
			int lx_k0 = Get_lx("SCFT", k0, N);
			
			if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) {
				(*Force1)(lx_k0, 0, k0) = (*Force1)(lx_k0, 0, k0) + complx(force_amp/4,0.0);
				(*Force3)(lx_k0, 0, k0) = (*Force3)(lx_k0, 0, k0) + I*(force_amp/4);
			}
		}
	}
	
	if (alias_switch == "DEALIAS")   Dealias_force();
}

//*********************************************************************************************


void IncFluid::Compute_force_Ekman_friction()
{	
	
	DP alpha = (*force_double_para)(1);
	
	*Force1 = *V1;	
	*Force2 = *V2; 
	*Force3 = *V3;
	
	*Force1 = complx(-alpha,0) * (*Force1);
	*Force2 = complx(-alpha,0) * (*Force2);
	*Force3 = complx(-alpha,0) * (*Force3);
	
	if (alias_switch == "DEALIAS")   Dealias_force();
	
}	



void IncFluid::Compute_force_Ekman_friction_const_energy_supply()
{	
	
	DP alpha = (*force_double_para)(4);
	
	DP inner_radius = (*force_double_para)(1);
	DP outer_radius = (*force_double_para)(2);
	DP energy_supply = (*force_double_para)(3);
	
	// first feed const eps force;
	Compute_force_const_energy_supply();
	
	// Ekman friction
	*VF_temp = *V1;	
	*Force1 = *Force1 + complx(-alpha,0) * (*VF_temp);
	
	*VF_temp = *V2;	
	*Force2 = *Force2 + complx(-alpha,0) * (*VF_temp);
	
	*VF_temp = *V3;	
	*Force3 = *Force3 + complx(-alpha,0) * (*VF_temp);
	
	if (alias_switch == "DEALIAS")   Dealias_force();
	
}	


//****************************** End of compute_force_LM.cc ***********************************


