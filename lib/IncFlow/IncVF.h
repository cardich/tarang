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

/*! \file  IncVF.h
 * 
 * @brief  Class declaration of IncVF, Incompressible Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   No known bugs
 */
 
//*********************************************************************************************

#ifndef _IncVF
#define _IncVF

#include "../fields/fields.h"			// blitz, fft etc. declared
#include "../plainfields/plainfields.h"
#include "Nlin.h"
#include "IncSF.h"
#include "EnergyTr.h"

#include <iostream>
#include <string>
#include <fstream>



//! @brief Incompressible vector field IncVF 
/*!
 *  Inherits CVF that contains the complex vector field. <BR>
 *  RVF that contains real vector field, typically Inverse tranform of CVF. <BR>
 *  NLIN that contains the nonlinear term \f$ N_i = \mathcal{F} (D_j V_j V_i) \f$.<BR>
 *  EnergyTr that contains energy transfer vars like flux, shell-to-shell transfers. <BR>
 * 
 *	Compute nonlinear terms <BR>
 *  Compute energy transfer functions: <BR>
 *	Isotropic: flux, shell-to-shell <BR>
 *  Anisotropic: ring-to-ring in spherical shell and in cylinderical shells.
 *
 *	@sa IncSF.h
 *	@sa Nlin.h
 *  @sa EnergyTr.h
 */
 
 
//*********************************************************************************************	

class IncVF: public CVF, public RVF, public NLIN 
{ 
 public:

	//!  Force along x \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<complx,3> *Force1;						
	
	//!  Force along y \f$ F_y(local_{N1}, N_2, N_3/2+1) \f$.
	Array<complx,3> *Force2;						
	
	//!  Force along z \f$ F_x(local_{N1}, N_2, N_3/2+1) \f$.
	Array<complx,3> *Force3;						
	
	//!  temp array \f$ (local_{N1}, N_2, N_3/2+1) \f$.
	Array<complx,3> *VF_temp;					
	
	//!  Another temp array \f$ (local_{N1}, N_2, N_3/2+1) \f$.			
	Array<complx,3> *VF_temp2;	
	
	//!  temp array \f$ (local_{N2}, N_1, N_3/2+1) \f$.
	Array<complx,3> *VF_temp_r;
																
	//!  Nusselt spectrum \f$ \sum \Re(V_1(\vec{K}) \theta^*(\vec{K})) \f$ over the shell.					
	Array<DP,1>		*shell_ek_cross_V1T;
	Array<DP,1>		*shell_ek_cross_V2T;
	Array<DP,1>		*shell_ek_cross_V3T;	
	
	//!  \f$ \sum \Re(V_1(\vec{K}) \theta^*(\vec{K})) \f$ over the rings.  
	//! Allocated only for anisotropic cases.
	Array<DP,2>		*ring_ek_cross_V1T;
	Array<DP,2>		*ring_ek_cross_V2T;
	Array<DP,2>		*ring_ek_cross_V3T;
	
	//!  \f$ \sum \Re(V_1(\vec{K}) \theta^*(\vec{K})) \f$ over the cylindrical rings.  
	//! Allocated only for anisotropic cases.
	Array<DP,2>		*cylinder_ring_ek_cross_V1T;
	Array<DP,2>		*cylinder_ring_ek_cross_V2T;
	Array<DP,2>		*cylinder_ring_ek_cross_V3T;					

	//!  spectrum of cross helicity \f$ \sum \Re(\vec{V}(\vec{K}) \cdot 
	/// \vec{W}^*(\vec{K})) \f$.		
	Array<DP,1>		*shell_ek_cross;
	
	//!  spectrum of energy feed by force 
	/// \f$ \sum \Re(\vec{F}(\vec{K}) \cdot \vec{V}^*(\vec{K})) \f$.
	Array<DP,1>		*shell_spectrum_force_Vk;
	
	
	//! Size of complex array of CVF, force etc.
	int		N[4];
	
	//! Basis type: FOUR or SCFT
	string	basis_type;
	
	//! Alias switch: ALIAS or DEALIAS
	string	alias_switch;
	
	int force_switch;		// 1 if on. 0 otherwise
	
	//!  Conversion factor for grid to actual wavenumber: \f$ f_i \f$.
	DP		kfactor[4];
	
	//! Dissipation coefficient appearing before laplacian \f$ \nu \f$.
	DP		dissipation_coefficient;
	int		hyper_dissipation_switch;
	
	//! Hyper_dissipation coefficient appearing before \f$\nabla^exponent \f$: 
	DP		hyper_dissipation_coefficient;
	int		hyper_dissipation_exponent;
  
	
	//************************* Anisotropic vars  *********************************************
	
	//! anisotropic spectrum
	int		sector_spectrum_input_scheme;
	
	//!  Number of sectors for ring spectrum.
	int		no_sectors_spectrum;					
	
	//!  sector_angles array; ring_n = [ sector_angle(n-1),sector_angle(n) );
	Array<DP,1>  *sector_angle_array_spectrum;			
	
	int		no_cylinder_slabs_spectrum;
	
	//! linearly spaced kz array
	Array<DP,1>  *cylinder_kpll_array_spectrum;		
	
	// anisotropic Hc
	Array<DP,2>		*ring_ek_cross;
	Array<DP,2>		*cylinder_ring_ek_cross;
	
	// Anisotropic force feed
	Array<DP,2>		*ring_spectrum_force_Vk;
	Array<DP,2>		*cylinder_ring_spectrum_force_Vk;
	

	//*****************************************************************************************

 public:
							
	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */						
	IncVF
	(
		int NN[], 
		string string_switches[], 
		Array<int,1> switches,
		int force_switch_para,
		DP *prog_kfactor, 
		DP diss_coefficient, 
		DP hyper_diss_coefficient, 
		int hyper_diss_exponent,	
		Array<int,1> misc_output_para
	);
	
	void Mult_field_exp_ksqr_dt(DP dt, DP a);
	void Mult_nlin_exp_ksqr_dt(DP dt, DP a);
	void Add_nlin_dt(DP dt);
	
	void Copy_field_to(CVF& W);
	void Copy_field_to(PlainCVF& W);
	void Copy_field_from(CVF& W);
	void Copy_field_from(PlainCVF& W);
	
	void Add_nlin_to_field(CVF& W, DP factor);
	void Add_nlin_to_field(PlainCVF& W, DP factor);
	
	void Dealias();
	void Dealias(IncVF& W);
	void Dealias(IncSF& T);
	void Dealias(IncVF& W, IncSF& T);
	
	void Dealias_force();
	void Dealias_force(IncVF& W);
	void Dealias_force(IncSF& T);
	void Dealias_force(IncVF& W, IncSF& T);
	
	void Dealias_nlin();
	void Dealias_nlin(IncVF& W);
	void Dealias_nlin(IncSF& T);
	void Dealias_nlin(IncVF& W, IncSF& T);
	
	/// Elsasser.cc file
			
	void UB_to_Elsasser_field(IncVF& W);
	void Elsasser_to_UB_field(IncVF& W);
	
	void UB_to_Elsasser_force(IncVF& W);
	void Elsasser_to_UB_force(IncVF& W);
	
	void UB_to_Elsasser_nlin(IncVF& W);
	void Elsasser_to_UB_nlin(IncVF& W);
	
	
	DP Get_cross_helicity(IncVF& W);
	DP Get_Nusselt_no(IncSF& T);
	
	void Compute_cross_vT_shell_spectrum(IncSF& T);
	void Compute_cross_vT_ring_spectrum(IncSF& T);
	void Compute_cross_vT_cylinder_ring_spectrum(IncSF& T);
	
	void Compute_cross_helicity_shell_spectrum(IncVF& W);
	void Compute_cross_helicity_ring_spectrum(IncVF& W);
	void Compute_cross_helicity_cylinder_ring_spectrum(IncVF& W);
	
	void Compute_force_shell_spectrum();
	void Compute_force_ring_spectrum();
	void Compute_force_cylinder_ring_spectrum();
	
	void Compute_divergence_nlin(); 
	void Compute_divergence_field();
  
	/// For computing FT[nabla.(V V)]
	
	void Compute_RSprod_diag();						///  Vr[i] -> Vr[i]^2 stored in nlin[i]
	void Compute_RSprod_offdiag();					/// Computes Vr[i]*Vr[j]  
	void Derivative_RSprod_VV();					/// Computes kj*FT[Vr[j]*Vr[i]
	
	void Compute_RSprod_diag(IncVF& W);				/// nlin[i]= Vr[i]*Wr[i]
	void Compute_RSprod_offdiag(IncVF& W);			/// Computes Vr[i]*Wr[j] 
	
	/// Computes kj*FT[Vr[j]*Wr[i]] and completes W.nlin[i]
	void Derivative_RSprod_VV(IncVF& W);		
	
	void Compute_RSprod_diag_ft_derivative();				// For transpose_order
	void Compute_RSprod_diag_ft_derivative(IncVF& W);
	void Forward_transform_derivative_RSprod_offdiag();
	void Forward_transform_derivative_RSprod_offdiag(IncVF& W);
	
		
	void Xderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B);
	void Yderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B);
	void Zderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B);
	
	void Xderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B);
	void Yderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B);
	void Zderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B);
	
	
	void Compute_nlin();							/// nlin[i] = FT[v.grad(v)][i]
	void Compute_nlin(IncSF& T);					///  nlin[i] = FT[v.grad(v)][i]; T.nlin= FT[v.grad(T)]
	
	void Compute_nlin_scalar(IncSF& T);				// scalar convection
	void Compute_nlin_RB(IncSF& T);					// RB convection
	void Compute_nlin_NonBoussinesq(IncSF& T);					// NonBoussinesq flows convection
	void Compute_nlin(IncVF& W);
	void Compute_nlin(IncVF& W, IncSF& T);  


	/// U.nlin[i] = FT[U.grad U]; W.nlin[i] = FT[W.grad W]   
	/// nlinWdU contains FT[W.grad U]; nlinUdW contains FT[U.grad W]	
	
	void Compute_nlin(IncVF& W, CVF& nlinWdU, CVF& nlinUdW); 
  
	void Compute_pressure();									// For Fluid, Passive Scalar, and MHD
	void Compute_shell_spectrum_pressure();
	void Compute_ring_spectrum_pressure();
	void Compute_cylinder_ring_spectrum_pressure();
	
	void Satisfy_reality_condition_field();
	void Satisfy_reality_condition_field(IncSF& T);
	void Satisfy_reality_condition_field(IncVF& W);
	void Satisfy_reality_condition_field(IncVF& W, IncSF& T);
	
	void Satisfy_reality_condition_force_field();
	void Satisfy_reality_condition_force_field(IncSF& T);
	void Satisfy_reality_condition_force_field(IncVF& W);
	void Satisfy_reality_condition_force_field(IncVF& W, IncSF& T);
	
	void Satisfy_reality_condition_nlin();
	void Satisfy_reality_condition_nlin(IncSF& T);
	void Satisfy_reality_condition_nlin(IncVF& W);
	void Satisfy_reality_condition_nlin(IncVF& W, IncSF& T);
	
	void Satisfy_reality_condition();
	void Satisfy_reality_condition(IncSF& T);
	void Satisfy_reality_condition(IncVF& W);
	void Satisfy_reality_condition(IncVF& W, IncSF& T);
	
	// Useful for init.cond
	
	void Last_component(int kx, int ky, complx &Vx, complx &Vy);
	void Last_component(int kx, int ky, int kz, complx &Vx, complx &Vy, complx &Vz);
	void Compute_VyVz(int kx, int ky, int kz, complx& Vx, complx Omega, complx &Vy, complx &Vz);
	
	void Add_complex_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz);
	void Assign_field_comp_conj(int kx, int ky, int kz, complx Vx, complx Vy, complx Vz);

	void Add_complex_conj_force(int kx, int ky, int kz, complx Fx, complx Fy, complx Fz);
	void Assign_force_and_comp_conj(int kx, int ky, int kz, complx Fx, complx Fy, complx Fz);
	
	void Zero_modes_RB_slip(IncSF& T);
	void Zero_modes_RB_slip(IncVF& W, IncSF& T);
	
	DP Get_Tk(int kx, int ky, int kz);

	void Get_random_vector(int kx, int ky, int kz, DP rand_range, complx& Vx, complx& Vy, complx& Vz);
	void Get_random_scalar(int kx, int ky, int kz, DP rand_range, complx& F);
	

	
//	void Put_random_vector_add_conj(int l1, int l2, int l3, DP amp);
	
	void Put_vector_amp_phase_comp_conj
	(
		int lx, int ly, int lz,
		int N[],
		DP amp, 
		DP phase1, DP phase2, DP phase3
	);
	
		//	void Put_vector_amp_phase_comp_conj_2Din3Dgrid(int lx, int ly, DP amp, DP phase);

	void Force_alpha(int kx, int ky, DP alpha);
	void Force_alpha_beta(int kx, int ky, int kz, DP alpha, DP beta);
	void Force_alpha_2Din3Dgrid(int kx, int ky, DP alpha);
	
	void Const_vector_field_alpha_beta(int kx, int ky, int kz, DP alpha, DP beta);
	
	void Setup_Taylor_Green_field(int k0, DP amp);
	void Setup_Taylor_Green_field_mag_field(int k0, DP amp);
	void Setup_ABC_field(int k0, DP amp, DP A, DP B, DP C);
	
	void Setup_Taylor_Green_force_field(int k0, DP amp);
	void Setup_ABC_force_field(int k0, DP amp, DP A, DP B, DP C);
	
	void Setup_SIX_MODE_field(int k0, DP amp101, DP amp011, DP amp112, DP h);
	void Setup_SIX_MODE_force_field(int k0, DP amp101, DP amp011, DP amp112, DP h);
	
	// Misc functions
	void free_slip_verticalwall_field();
	void free_slip_verticalwall_field(IncSF& T);
	void free_slip_verticalwall_field(IncVF& W);	
	void free_slip_verticalwall_field(IncVF& W, IncSF& T);
	
	void free_slip_verticalwall_force_field();
	void free_slip_verticalwall_force_field(IncSF& T);
	void free_slip_verticalwall_force_field(IncVF& W);	
	void free_slip_verticalwall_force_field(IncVF& W, IncSF& T);
	
	DP Get_mag_V0();
	
	void Model_initial_shell_spectrum_Pope(Array<DP,1> Sk, DP epsilon);
			
};

#endif

//************************************  End of IncVF.h  ***************************************


