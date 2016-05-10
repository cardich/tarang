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


/*! \file  single_time_step.cc
 * 
 * @brief single time step fields using Euler, RK2, and semi-implicit scheme.
 * @sa void IncFluid::Single_time_step_EULER(DP dt)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"

//*********************************************************************************************
//*********************************************************************************************
/** @brief Single time step
 *
 * @return \f$ \vec{V}(t+dt) = [\vec{V}(t) *\exp(-K^2*\nu*dt*a) 
 *								+ c dt \vec{N}(t+dt) *\exp(-K^2*\nu*dt*b) ] 
 *
 * @note For hyperviscosity:  \f$ \exp(-K^2*\nu*dt) \f$ is replaced by
 *					 $ \exp(-K^2*\nu*dt - K^4*\nu_h dt) \f$
 */
void IncFluid::Single_time_step(DP dt, DP a, DP b, DP c)
{	
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
	{
		Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS))
	{
		Mult_nlin_exp_ksqr_dt(dt, b);
		Add_nlin_dt(c*dt);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, a);
		Add_nlin_dt(c*dt);	
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS))
	{
		Add_nlin_dt(c*dt);									// V = V + c*dt*nlin
		Mult_field_exp_ksqr_dt(dt, a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, (a-b));
		Add_nlin_dt(c*dt);	
		Mult_field_exp_ksqr_dt(dt, b);
	}
}


//*********************************************************************************************	
//	Passive scalar  
void IncFluid::Single_time_step(IncSF& T, DP dt, DP a, DP b, DP c)
{	
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Single_time_step_scalar(T,dt, a, b, c);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq")) 
		Single_time_step_RB(T,dt, a, b, c);
	
}	

//*********************************************************************************
void IncFluid::Single_time_step_scalar(IncSF& T, DP dt, DP a, DP b, DP c)
{	
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
	{
		Add_nlin_dt(c*dt);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS))
	{
		Mult_nlin_exp_ksqr_dt(dt, b);
		Add_nlin_dt(c*dt);
		
		T.Mult_nlin_exp_ksqr_dt(dt, b);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, a);
		Add_nlin_dt(c*dt);	
		
		T.Mult_field_exp_ksqr_dt(dt, a);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS))
	{
		Add_nlin_dt(c*dt);									// V = V + c*dt*nlin
		Mult_field_exp_ksqr_dt(dt, a);
		
		T.Add_nlin_dt(c*dt);								
		T.Mult_field_exp_ksqr_dt(dt, a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, (a-b));
		Add_nlin_dt(c*dt);	
		Mult_field_exp_ksqr_dt(dt, b);
		
		T.Mult_field_exp_ksqr_dt(dt, (a-b));
		T.Add_nlin_dt(c*dt);	
		T.Mult_field_exp_ksqr_dt(dt, b);
	}
	
}


//*********************************************************************************
void IncFluid::Single_time_step_RB(IncSF& T, DP dt, DP a, DP b, DP c)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{	
		Single_time_step(dt, a, b, c);						// V advance
		
		*T.F = *V1;
		Array_divide_ksqr_SCFT(N, *T.F, kfactor);				// theta(k) = v1(k)/k^2
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
	
		*V1 = *nlin1;
		*V2 = *nlin2;
		*V3 = *nlin3;
		
		Array_divide_ksqr_SCFT(N, *V1, kfactor);				// u(k) = nlin(k)/k^2	
		Array_divide_ksqr_SCFT(N, *V2, kfactor);
		Array_divide_ksqr_SCFT(N, *V3, kfactor);
		
		if (globalvar_RB_Uscaling == "ULARGE") {
			(*V1) = (*V1)/dissipation_coefficient;
			(*V2) = (*V2)/dissipation_coefficient;
			(*V3) = (*V3)/dissipation_coefficient;
		}
		
		// Time advance T
		if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
		{
			T.Add_nlin_dt(c*dt);
		}
		
		else if ((abs(a) < MYEPS) && (abs(b) > MYEPS))
		{
			T.Mult_nlin_exp_ksqr_dt(dt, b);
			T.Add_nlin_dt(c*dt);
		}
		
		else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))
		{
			T.Mult_field_exp_ksqr_dt(dt, a);
			T.Add_nlin_dt(c*dt);
		}
		
		else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS))
		{
			T.Add_nlin_dt(c*dt);								
			T.Mult_field_exp_ksqr_dt(dt, a);
		}
		
		else if ((abs(a) > MYEPS) && (abs(b) > MYEPS))
		{
			T.Mult_field_exp_ksqr_dt(dt, (a-b));
			T.Add_nlin_dt(c*dt);	
			T.Mult_field_exp_ksqr_dt(dt, b);
		}
	}
	
	else 
		Single_time_step_scalar(T, dt, a, b, c);
	
}  


//*********************************************************************************************
//	MHD		

void IncFluid::Single_time_step(IncVF& W, DP dt, DP a, DP b, DP c)
{  
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
	{
		Add_nlin_dt(c*dt);
		W.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS))
	{
		Mult_nlin_exp_ksqr_dt(dt, b);
		Add_nlin_dt(c*dt);
		
		W.Mult_nlin_exp_ksqr_dt(dt, b);
		W.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, a);
		Add_nlin_dt(c*dt);	
		
		W.Mult_field_exp_ksqr_dt(dt, a);
		W.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS))
	{
		Add_nlin_dt(c*dt);									// V = V + c*dt*nlin
		Mult_field_exp_ksqr_dt(dt, a);
		
		W.Add_nlin_dt(c*dt);								
		W.Mult_field_exp_ksqr_dt(dt, a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, (a-b));
		Add_nlin_dt(c*dt);	
		Mult_field_exp_ksqr_dt(dt, b);
		
		W.Mult_field_exp_ksqr_dt(dt, (a-b));
		W.Add_nlin_dt(c*dt);	
		W.Mult_field_exp_ksqr_dt(dt, b);
	}
}  


//*********************************************************************************************
//	MHD + passive scalar  

void IncFluid::Single_time_step(IncVF& W, IncSF& T,  DP dt, DP a, DP b, DP c)
{
	
	if ((abs(a) < MYEPS) && (abs(b) < MYEPS))
	{
		Add_nlin_dt(c*dt);
		W.Add_nlin_dt(c*dt);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a) < MYEPS) && (abs(b) > MYEPS))
	{
		Mult_nlin_exp_ksqr_dt(dt, b);
		Add_nlin_dt(c*dt);
		
		W.Mult_nlin_exp_ksqr_dt(dt, b);
		W.Add_nlin_dt(c*dt);
		
		T.Mult_nlin_exp_ksqr_dt(dt, b);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(b) < MYEPS) && (abs(a) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, a);
		Add_nlin_dt(c*dt);	
		
		W.Mult_field_exp_ksqr_dt(dt, a);
		W.Add_nlin_dt(c*dt);
		
		T.Mult_field_exp_ksqr_dt(dt, a);
		T.Add_nlin_dt(c*dt);
	}
	
	else if ((abs(a-b) < MYEPS) && (abs(a) > MYEPS))
	{
		Add_nlin_dt(c*dt);									// V = V + c*dt*nlin
		Mult_field_exp_ksqr_dt(dt, a);
		
		W.Add_nlin_dt(c*dt);								
		W.Mult_field_exp_ksqr_dt(dt, a);
		
		T.Add_nlin_dt(c*dt);								
		T.Mult_field_exp_ksqr_dt(dt, a);
	}
	
	else if ((abs(a) > MYEPS) && (abs(b) > MYEPS))
	{
		Mult_field_exp_ksqr_dt(dt, (a-b));
		Add_nlin_dt(c*dt);	
		Mult_field_exp_ksqr_dt(dt, b);
		
		W.Mult_field_exp_ksqr_dt(dt, (a-b));
		W.Add_nlin_dt(c*dt);	
		W.Mult_field_exp_ksqr_dt(dt, b);
		
		T.Mult_field_exp_ksqr_dt(dt, (a-b));
		T.Add_nlin_dt(c*dt);	
		T.Mult_field_exp_ksqr_dt(dt, b);
	}
	
}

//*********************************************************************************************
//*********************************************************************************************
/** @brief Compute_force_TO_rhs()
 *
 * @return Computes force and nlin; add them; compute pressure; combine them to compute rhs.
 *
 */
void IncFluid::Compute_force_TO_rhs()
{
	Compute_force();
	
	Compute_nlin();										// Compute nlin using V(t+dt/2)
	
//	Satisfy_reality_condition_nlin();
	
	Add_force();										// nlin = nlin - f
	
	Compute_pressure();									// Compute pressure using V(t+dt/2)
	
	Compute_rhs();
	
}

void IncFluid::Compute_force_TO_rhs(IncSF& T)
{
	Compute_force(T);
	
	Compute_nlin(T);									// Compute nlin using V(t+dt/2)
	
//	Satisfy_reality_condition_nlin(T);
	
	Add_force(T);										// nlin = nlin - f
	
	Compute_pressure();									// Compute pressure using V(t+dt/2)
	
	Compute_rhs(T);		
}


void IncFluid::Compute_force_TO_rhs(IncVF& W)
{
	Compute_force(W);
	
	Compute_nlin(W);										// Compute nlin using V(t+dt/2)
	
//	Satisfy_reality_condition_nlin(W);
	
	Add_force(W);											// nlin = nlin - f
	
	Compute_pressure();										// Compute pressure using V(t+dt/2)
	
	Compute_rhs(W);
}


void IncFluid::Compute_force_TO_rhs(IncVF& W, IncSF& T)
{
	Compute_force(W, T);
	
	Compute_nlin(W, T);										// Compute nlin using V(t+dt/2)
	
//	Satisfy_reality_condition_nlin(W, T);
	
	Add_force(W, T);										// nlin = nlin - f
	
	Compute_pressure();										// Compute pressure using V(t+dt/2)
	
	Compute_rhs(W, T);
}


void IncFluid::Compute_force_TO_rhs_NonBoussinesq(IncSF& T)
{
	Compute_force(T);
	
	Compute_nlin_NonBoussinesq(T);									
	// Compute nlin using V(t+dt/2)
	
	//	Satisfy_reality_condition_nlin(T);
	
	Add_force(T);										// nlin = nlin - f
	
	Compute_pressure();									// Compute pressure using V(t+dt/2)
	
	Compute_rhs(T);		
}


//*********************************************************************************************
//*********************************************************************************************
/** @brief Compute_RK4_parts(PlainCVF& tot_Vrhs, DP dt, DP b, DP factor) conmputes Ci's for 
 *		for computing fields at t+dt.
 *
 * @return tot_Vrhs += factor*RHS(t)*exp(-nu k^2 b*dt)
 *			RHS is contained in *nlin
 */
void IncFluid::Compute_RK4_parts(PlainCVF& tot_Vrhs, DP dt, DP b, DP factor)
{	
	Mult_nlin_exp_ksqr_dt(Tdt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
	Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
}


void IncFluid::Compute_RK4_parts(IncSF& T, PlainCVF& tot_Vrhs, PlainCSF& tot_Srhs, 
								 DP dt, DP b, DP factor)
{	
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_RK4_parts_scalar(T,tot_Vrhs, tot_Srhs, dt, b, factor);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Compute_RK4_parts_RB(T,tot_Vrhs, tot_Srhs, dt, b, factor);
}


void IncFluid::Compute_RK4_parts_scalar(IncSF& T, PlainCVF& tot_Vrhs, 
										PlainCSF& tot_Srhs, DP dt, DP b, DP factor)
{	
	Mult_nlin_exp_ksqr_dt(dt, b);								// *nlin = *nlin x exp(-nu k^2 dt/d)
	Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
	
	T.Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
	T.Add_nlin_to_field(tot_Srhs, factor);					// rhs(t0, u0) in nlin
}


void IncFluid::Compute_RK4_parts_RB(IncSF& T, PlainCVF& tot_Vrhs, 
									PlainCSF& tot_Srhs, DP dt, DP b, DP factor)
{	
	if ((globalvar_Pr_switch != "PRZERO") && (globalvar_Pr_switch != "PRINFTY"))
	{
		Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
		Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
		
		T.Mult_nlin_exp_ksqr_dt(dt, b);							
		T.Add_nlin_to_field(tot_Srhs, factor);					
	}
	
	else if (globalvar_Pr_switch == "PRZERO")
	{
		Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
		Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
	}
	
	else if (globalvar_Pr_switch == "PRINFTY")
	{
		T.Mult_nlin_exp_ksqr_dt(dt, b);							
		T.Add_nlin_to_field(tot_Srhs, factor);	
	}
	
	
}


void IncFluid::Compute_RK4_parts(IncVF& W, PlainCVF& tot_Vrhs, PlainCVF& tot_Wrhs, 
								 DP dt, DP b, DP factor)
{	
	Mult_nlin_exp_ksqr_dt(dt, b);								// *nlin = *nlin x exp(-nu k^2 dt/d)
	Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
	
	W.Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
	W.Add_nlin_to_field(tot_Wrhs, factor);					// rhs(t0, u0) in nlin
}


void IncFluid::Compute_RK4_parts(IncVF& W, IncSF& T, PlainCVF& tot_Vrhs, PlainCVF& tot_Wrhs, 
								 PlainCSF& tot_Srhs, DP dt, DP b, DP factor)
{	
	Mult_nlin_exp_ksqr_dt(dt, b);								// *nlin = *nlin x exp(-nu k^2 dt/d)
	Add_nlin_to_field(tot_Vrhs, factor);					// rhs(t0, u0) in nlin
	
	T.Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
	T.Add_nlin_to_field(tot_Srhs, factor);					// rhs(t0, u0) in nlin
	
	W.Mult_nlin_exp_ksqr_dt(dt, b);							// *nlin = *nlin x exp(-nu k^2 dt/d)
	W.Add_nlin_to_field(tot_Wrhs, factor);					// rhs(t0, u0) in nlin
}


//*********************************************************************************************
// total_rho = rho'+ rho(x) + mean density(=1); rho(x) = -x
// z=0 plane is the reference plane where the total density = 1.
// delta_rho/rho = alpha DT * (T_cond + theta)
// Density at the bottom plate is rho_0 (mean density). Here normalized to 1.

void IncFluid::Compute_density(IncSF& T)
{
	*T.Fr = *T.F;
	T.RS_Inverse_transform(*VF_temp_r);	
	
	int kx;
	DP x;
	for (int l1=0; l1<local_N1; l1++) {
		kx = Get_kx("SCFT", l1, N);
		x = (kx+0.5)/N[1];
		
		for (int l2=0; l2<N[2]; l2++)
			for (int l3=0; l3<(N[3]/2); l3++)
				(*T.Fr)(l1, l2, l3) += (0.5-x);
	}
	
	real(*T.Fr) = globalvar_alpha_DT*real(*T.Fr) + 1.0; 
	imag(*T.Fr) = globalvar_alpha_DT*imag(*T.Fr) + 1.0; // mean density = 1.0
	
	
}

//**********************************   End of Single_time_step.cc  ****************************



