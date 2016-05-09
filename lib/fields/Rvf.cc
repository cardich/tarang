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



/*! \file  Rvf.cc
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @sa	Rvf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 

#include "Rvf.h"

/*
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
 */

//*********************************************************************************************

RVF::RVF
(
	int NN[], 
	string string_switches[],
	Array<int,1> switches,
	DP *kfactor, 
	Array<int,1> misc_output_para
) 
{ 

	RV_basis_type = string_switches[1];
	
	RV_structure_fn_switch			= switches(6);
	RV_planar_structure_fn_switch	= switches(7);
	RV_structure_fn_approx_switch	= switches(24);
	
	RV_structurefn_q_min			= misc_output_para(4);
	RV_structurefn_q_max			= misc_output_para(5);
	RV_structurefn_r_arraysize		= misc_output_para(6);
	RV_rmax_div_st_fn_r_max			= misc_output_para(7);   // rmax/str_fn_rmax
	RV_structurefn_neighbour_index	= misc_output_para(8);   // neighbour number
	
    for (int i=1; i<=3; i++) {
		Nrv[i]=NN[i];
		RV_xfactor[i] = 1/kfactor[i];
	}
     
	
	if (RV_basis_type == "FOUR") 
		for (int i=1; i<=3; i++) {
			RV_L[i] = 2*M_PI/kfactor[i];
			RV_Delta_x[i] = RV_L[i]/Nrv[i];
		}
	
	else if (RV_basis_type == "SCFT") {
		RV_L[1] = M_PI/kfactor[1];
		RV_Delta_x[1] = RV_L[1]/Nrv[1];
		
		for (int i=2; i<=3; i++)	{
			RV_L[i] = 2*M_PI/kfactor[i];
			RV_Delta_x[i] = RV_L[i]/Nrv[i];
		}
	}
	  
	// If Transpose switch on, FT(*Vtr) = *V
	
#ifdef TRANSPOSE
	V1r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	V2r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	V3r = new Array<complx,3>(local_N2, Nrv[1], Nrv[3]/2+1);  
	
#else 
	V1r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1);			
	V2r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1);  
	V3r = new Array<complx,3>(local_N1, Nrv[2], Nrv[3]/2+1); 
#endif	

	(*V1r) = 0.0; 
	(*V2r) = 0.0; 
	(*V3r) = 0.0;
	
	// Structure function
	
	if ((RV_structure_fn_switch == 1) || (RV_planar_structure_fn_switch == 1))
	{
		RV_structurefn_r_max = Max_radius_inside_real_space(RV_basis_type, Nrv, RV_Delta_x)	
									/RV_rmax_div_st_fn_r_max;
		
		int q_range = RV_structurefn_q_max - RV_structurefn_q_min;
		
		RV_St = new Array<DP,3>(RV_structurefn_r_arraysize+1, q_range+1, 2);
		(*RV_St) = 0.0;
		
		RV_St_count = new Array<DP,1>(RV_structurefn_r_arraysize+1);
		(*RV_St_count) = 0.0;
		
		// Contains averaged over all procs in master proc.
		RV_St_final = new Array<DP,3>(RV_structurefn_r_arraysize+1, q_range+1, 2);
		(*RV_St_final) = 0.0;
		
		RV_St_count_final = new Array<DP,1>(RV_structurefn_r_arraysize+1);
		(*RV_St_count_final) = 0.0;
		

		// To optimize, we take max distance to be N[i]/RV_rmax_div_st_fn_r_max-1,
	/*	if (RV_structure_fn_approx_switch == 0)
			dist_farthest_proc = numprocs-1;
		
		else if (RV_structure_fn_approx_switch == 1) */
			dist_farthest_proc = numprocs/RV_rmax_div_st_fn_r_max-1; 
	}	
	
	if (RV_planar_structure_fn_switch == 1) {
	
		if (globalvar_anisotropy_switch == 1)
			RV_structurefn_pll_ind_max = Nrv[1];
		
		else if (globalvar_anisotropy_switch == 2)
			RV_structurefn_pll_ind_max = Nrv[2];
		
		else if (globalvar_anisotropy_switch == 3)
			RV_structurefn_pll_ind_max = Nrv[3];
	}
	
	
	
}

/**********************************************************************************************

				Forward transform: Inplace 
   
**********************************************************************************************/

void RVF::RV_Forward_transform(Array<complx,3> temp_r)
{
	Forward_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// SFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// CFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 0);			// CFT	
}

/**********************************************************************************************

			Forward_transform(*Vir) = Vi 
				Vir unchanged.
   
**********************************************************************************************/


void RVF::RV_Forward_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
)
{
	temp_r = *V1r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V1, 1);			// SFT
	
	temp_r = *V2r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V2, 0);			// CFT
	
	temp_r = *V3r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V3, 0);			// CFT	
}


/**********************************************************************************************

		Inplace Inverse Fourier transform.
   
**********************************************************************************************/

void RVF::RV_Inverse_transform(Array<complx,3> temp_r)
{

	Inverse_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// ISFT

	Inverse_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// ICFT

	Inverse_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 0);			// ICFT

}


/**********************************************************************************************

			Inverse_transform(Vi) = *Vir 
				Keeping Vi unchanged....
				temp = N1 x N2 x N3
   
**********************************************************************************************/

void RVF::RV_Inverse_transform_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp
)
{
	temp = V1;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V1r, 1);	// ISFT
	
	temp = V2;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V2r, 0);	// ICFT
	
	temp = V3;
	Inverse_transform_array_transpose_order(RV_basis_type, Nrv, temp, *V3r, 0);	// ICFT

}


/**********************************************************************************************

			 Forward  transform RSprod 
			  3D: *Vr[i] With SCFT -- 1,3 SFT;  2 CFT
   
**********************************************************************************************/


void RVF::RV_Forward_transform_RSprod(Array<complx,3> temp_r)
{	
	Forward_transform_array(RV_basis_type, Nrv, *V1r, temp_r, 1);			// SFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V2r, temp_r, 0);			// CFT
	
	Forward_transform_array(RV_basis_type, Nrv, *V3r, temp_r, 1);			// SFT	
}
	

// *Vir unchanged..
void RVF::RV_Forward_transform_RSprod_transpose_order
(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
)
{
	temp_r = *V1r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V1, 1);			// SFT
	
	temp_r = *V2r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V2, 0);			// CFT
	
	temp_r = *V3r;
	Forward_transform_array_transpose_order(RV_basis_type, Nrv, temp_r, V3, 1);			// SFT	
}


/**********************************************************************************************

			Dealiase RVF
   
**********************************************************************************************/

void RVF::RV_Dealias()
{
	Dealias_array(RV_basis_type, Nrv, *V1r);
	
	Dealias_array(RV_basis_type, Nrv, *V2r);
	
	Dealias_array(RV_basis_type, Nrv, *V3r);
}


/**********************************************************************************************

			 Outputs Vi in real space
   
**********************************************************************************************/

void RVF::RV_Output(ofstream& fileout, Array<complx,3> temp, string format)
{		
	Output_asreal(fileout, Nrv, *V1r, temp, format);
	Output_asreal(fileout, Nrv, *V2r, temp, format); 
	Output_asreal(fileout, Nrv, *V3r, temp, format);
}

void RVF::RV_Output_transpose_order(ofstream& fileout, Array<complx,3> temp_array, string format)
{	
	Output_asreal_transpose_order(fileout, Nrv, *V1r, temp_array, format);
	Output_asreal_transpose_order(fileout, Nrv, *V2r, temp_array, format); 
	Output_asreal_transpose_order(fileout, Nrv, *V3r, temp_array, format);
}

void RVF::RV_input(ifstream& file_in, Array<complx,3> temp_array, string format)
{
	Read_real_field_data_MPI(RV_basis_type, file_in, Nrv, *V1r, temp_array, format);
	Read_real_field_data_MPI(RV_basis_type, file_in, Nrv, *V2r, temp_array, format);
	Read_real_field_data_MPI(RV_basis_type, file_in, Nrv, *V3r, temp_array, format);
}


void RVF::RV_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3)
{
	
	if (RV_structure_fn_approx_switch == 0) {
		Qi1_start = Qi2_start = Qi3_start = 0;
		
		Qi1_end = local_N1-1;
		Qi2_end = Nrv[2]-1;
		Qi3_end = Nrv[3]/2-1;
	}
	
	else if (RV_structure_fn_approx_switch == 1) {
		Qi1_start = Pi1;
		Qi2_start = Pi2;
		Qi3_start = Pi3;
		
		Qi1_end = min(((int) (Pi1 + Nrv[1]/RV_rmax_div_st_fn_r_max)), local_N1-1);
		Qi2_end = min(((int) (Pi2 + Nrv[2]/RV_rmax_div_st_fn_r_max)), Nrv[2]-1);
		Qi3_end = min(((int) (Pi3 + Nrv[3]/(2*RV_rmax_div_st_fn_r_max))), Nrv[3]/2-1);
		// To optimize, we take max distance to be min of the above two.
	}
}

//**************************  End of RVF class definitions ************************************






