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
/*! \file  time_advance_SF.cc
 * 
 * @brief Time advances velocity and passive scalar fields in Incompressible flow by one unit
 *
 * Status before entering the fn   <BR>
 *		nlin = FT[U.grad U]; U.F = p(k);  <BR>
 *		T.nlin = FT[U.grad T]
 *		(computed by Compute_nlin and compute_pressure)  <BR>
 *	
 *	Procedure (EULER)	 <BR>											
 *		1. nlin[i] = -nlin[i] - grad(pressure) in fn compute_rhs()     <BR>          
 *		2. V,T(t+dt) = V,T(t) + dt*nlin(V,T) in fn single_time_step() <BR>
 *		
 *	Procedure (RK2)	 <BR>
 *		1. Save V[i] in Vcopy[i] <BR>
 *		2. Compute_rhs(T): nlin[i] = -nlin[i] - grad(pressure) and T.nlin = -T.nlin <BR>
 *		3. V,T(t+dt) = V,T(t) + (dt/2)*nlin(V,T) using fn single_time_step(): MID POINT <BR>
 *		4. Compute_nlin(V(t+dt/2),T(t+dt/2)) <BR>
 *		5. compute p(t+dt/2) using compute_pressure() <BR>
 *		6. V,T[i] = V,Tcopy[i] (copy back) <BR>
 *		7. V,T(t+dt) = V,T(t) + dt*nlinU,T(t+dt/2) using fn single_time_step()  <BR>
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
                
   
#include "../IncFluid.h"


//********************************************************************************************* 


void IncFluid::Time_advance_NonBoussinesq(IncSF& T)
{

	
//*********************************************************************************************
//		EULER

	if (integ_scheme == "EULER") 
	{ 
		
		Compute_rhs(T);  
	
		//  Temperature fluctuations..
		//	Single_time_step(T, Tdt, 1, 1, 1)::: 
		//  Add_nlin_dt(dt) and Mult_field_exp_ksqr_dt(1, dt);
		T.Add_nlin_dt(Tdt);
		T.Mult_field_exp_ksqr_dt(Tdt, 1.0);
		
		
		// total_rho = rho'+ rho(x) + mean density(=1); rho(x) = -x
		// delta_rho/rho = alpha DT * (T_cond + theta)
		// Density at the bottom plate is rho_0 (mean density). Here normalized to 1.
		// *T.Fr contains total density rho/rho_0 after this.
		Compute_density(T); 
		
		*V1r = *V1;
		*V2r = *V2;
		*V3r = *V3;
		
		RV_Inverse_transform(*VF_temp_r);
		
		Array_real_mult(N, *V1r, *T.Fr, *V1r);						
		Array_real_mult(N, *V2r, *T.Fr, *V2r);	
		Array_real_mult(N, *V3r, *T.Fr, *V3r);	
		
		RV_Forward_transform(*VF_temp_r);
		
		*V1 = *V1r;
		*V2 = *V2r;
		*V3 = *V3r; 
		
		Add_nlin_dt(Tdt);
		Mult_field_exp_ksqr_dt(Tdt, 1.0);				// \hat(rho u_i) time advanced..
		
		// Now get u_i
		CV_Inverse_transform(*VF_temp_r);
		
		// Vi = Vi/rho
		Array_real_divide(N, *V1, *T.Fr, *V1);						
		Array_real_divide(N, *V2, *T.Fr, *V2);	
		Array_real_divide(N, *V3, *T.Fr, *V3); 
		
		CV_Forward_transform(*VF_temp_r); 
		// *Vi now has the final field in the Fourier space.
	 
	}
		
		
//*********************************************************************************************		
//		RK2			
	else if (integ_scheme == "RK2") 
	{
		
		// Allocated once and for all- Vcopy,Scopy  (static)
		static PlainCVF Vcopy(N);							
		static PlainCSF Scopy(N); 

		
		Copy_field_to(Vcopy); 
		T.Copy_field_to(Scopy);								// Vcopy[i] <- V[i] ; Scopy <- T.F     
  
		Compute_rhs(T);  
		  
		// Temp fluctuations
		// Single_time_step(T, Tdt, 0.5, 0.5, 0.5);		:Goto the mid point
		T.Add_nlin_dt(0.5*Tdt);
		T.Mult_field_exp_ksqr_dt(Tdt, 0.5);
		
		Compute_density(T); 
		
		*V1r = *V1;
		*V2r = *V2;
		*V3r = *V3;
		
		RV_Inverse_transform(*VF_temp_r);
		
		Array_real_mult(N, *V1r, *T.Fr, *V1r);						
		Array_real_mult(N, *V2r, *T.Fr, *V2r);	
		Array_real_mult(N, *V3r, *T.Fr, *V3r);	
		
		RV_Forward_transform(*VF_temp_r);
		
		*V1 = *V1r;
		*V2 = *V2r;
		*V3 = *V3r; 
		
		Add_nlin_dt(0.5*Tdt);
		Mult_field_exp_ksqr_dt(Tdt, 0.5);				// T(rho u_i) time advanced..
		
		// Now get u_i
		CV_Inverse_transform(*VF_temp_r);
		
		// Vi = Vi/rho
		Array_real_divide(N, *V1, *T.Fr, *V1);						
		Array_real_divide(N, *V2, *T.Fr, *V2);	
		Array_real_divide(N, *V3, *T.Fr, *V3);
		
		CV_Forward_transform(*VF_temp_r);
		// *Vi now has the final field in the Fourier space.
		
		
		Compute_force_TO_rhs_NonBoussinesq(T);
  

		Copy_field_from(Vcopy); 
		T.Copy_field_from(Scopy);							// Copy back into V[i],T

		// Single_time_step(T, Tdt, 1, 0.5, 1);							
		// Time-step by Tdt now using the mid-point slopes 
		
		T.Mult_field_exp_ksqr_dt(Tdt, 0.5);
		T.Add_nlin_dt(Tdt);	
		T.Mult_field_exp_ksqr_dt(Tdt, 0.5);
		
		Compute_density(T); 
		
		*V1r = *V1;
		*V2r = *V2;
		*V3r = *V3;
		
		RV_Inverse_transform(*VF_temp_r);
		
		Array_real_mult(N, *V1r, *T.Fr, *V1r);						
		Array_real_mult(N, *V2r, *T.Fr, *V2r);	
		Array_real_mult(N, *V3r, *T.Fr, *V3r);	
		
		RV_Forward_transform(*VF_temp_r);
		
		*V1 = *V1r;
		*V2 = *V2r;
		*V3 = *V3r; 
		
		Mult_field_exp_ksqr_dt(Tdt, 0.5);
		Add_nlin_dt(Tdt);
		Mult_field_exp_ksqr_dt(Tdt, 0.5);				// \hat(rho u_i) time advanced..
		
		// Now get u_i
		CV_Inverse_transform(*VF_temp_r);
		
		// Vi = Vi/rho
		Array_real_divide(N, *V1, *T.Fr, *V1);						
		Array_real_divide(N, *V2, *T.Fr, *V2);	
		Array_real_divide(N, *V3, *T.Fr, *V3);
		
		CV_Forward_transform(*VF_temp_r);
		// *Vi now has the final field in the Fourier space.
	}
	
	// if free_slip_verticalwall condition is initialized as IC, it should
	// be satisfied at all times, yet we set it again just to make sure everytime step.
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field(T);
	
	// For all schemes
	if (alias_switch == "DEALIAS")		Dealias(T);					// Keep V(k) dealiased	

}

//**********************************   End of time_advance_SF.cc  *****************************


