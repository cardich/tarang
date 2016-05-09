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


/*! \file  Rsf.cc
 * 
 * @brief  Class declaration of Rsf, a Real Scalar Field 
 *
 * @sa	Rsf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug No known bugs
 */
 

#include "Rsf.h"


//*********************************************************************************************

//  Class constructor; Allocation for F, Ek etc. Initialize the arrays.

RSF::RSF
(
	int NN[], 
	string string_switches[],
	Array<int,1> switches,
	DP *kfactor, 
	Array<int,1> misc_output_para
) 
{

	RS_basis_type = string_switches[1];
	
	RS_structure_fn_switch			= switches(6);
	RS_planar_structure_fn_switch	= switches(7);
	RS_structure_fn_approx_switch	= switches(23);
	
	RS_structurefn_q_min			= misc_output_para(4);
	RS_structurefn_q_max			= misc_output_para(5);
	RS_structurefn_r_arraysize		= misc_output_para(6);
	RS_rmax_div_st_fn_r_max			= misc_output_para(7);   // rmax/str_fn_rmax
	RS_structurefn_neighbour_index	= misc_output_para(8);   // points to skip in each dirn
	
    for (int i=1; i<=3; i++) {
		Nrs[i] = NN[i];
		RS_xfactor[i] = 1/kfactor[i];
	}
	
	if (RS_basis_type == "FOUR") 
		for (int i=1; i<=3; i++) {
			RS_L[i] = 2*M_PI/kfactor[i];
			RS_Delta_x[i] = RS_L[i]/Nrs[i];
		}
	
	else if (RS_basis_type == "SCFT") {
		RS_L[1] = M_PI/kfactor[1];
		RS_Delta_x[1] = RS_L[1]/Nrs[1];
		
		for (int i=2; i<=3; i++)	{
			RS_L[i] = 2*M_PI/kfactor[i];
			RS_Delta_x[i] = RS_L[i]/Nrs[i];
		}
	}
	
	// Structure function
	
	
	if ((RS_structure_fn_switch == 1) || (RS_planar_structure_fn_switch == 1))
	{
		RS_structurefn_r_max = Max_radius_inside_real_space(RS_basis_type, Nrs, RS_Delta_x)
								/RS_rmax_div_st_fn_r_max;
		
		int q_range = RS_structurefn_q_max - RS_structurefn_q_min;
		
		RS_St = new Array<DP,2>(RS_structurefn_r_arraysize+1, q_range+1);
		(*RS_St) = 0.0;
		
		RS_St_count = new Array<DP,1>(RS_structurefn_r_arraysize+1);
		(*RS_St_count) = 0.0;
		
		// Contains averaged over all procs in master proc.
		RS_St_final = new Array<DP,2>(RS_structurefn_r_arraysize+1, q_range+1);
		(*RS_St_final) = 0.0;
		
		RS_St_count_final = new Array<DP,1>(RS_structurefn_r_arraysize+1);
		(*RS_St_count_final) = 0.0;
		
		// To optimize, we take max distance to be N[i]/RS_rmax_div_st_fn_r_max,
		/*	if (RS_structure_fn_approx_switch == 0)
			dist_farthest_proc = numprocs-1;
		
		else if (RS_structure_fn_approx_switch == 1) */
		dist_farthest_proc = numprocs/RS_rmax_div_st_fn_r_max-1;
	}
	
	
	if (RS_planar_structure_fn_switch == 1)
	{
		
		if (globalvar_anisotropy_switch == 1)
			RS_structure_fn_rpll_max = Nrs[1];
		
		else if (globalvar_anisotropy_switch == 2)
			RS_structure_fn_rpll_max = Nrs[2];
		
		else if (globalvar_anisotropy_switch == 3)
			RS_structure_fn_rpll_max = Nrs[3];
	}
	

#ifdef TRANSPOSE
	Fr = new Array<complx,3>(local_N2, Nrs[1], Nrs[3]/2+1); 
#else    
    Fr = new Array<complx,3>(local_N1, Nrs[2], Nrs[3]/2+1);  
#endif	

    *Fr = 0.0;   // initialize
	
}

/**********************************************************************************************

					Inplace Forward Fourier transform
   
**********************************************************************************************/

void RSF::RS_Forward_transform(Array<complx,3> temp_r)
{
	Forward_transform_array(RS_basis_type, Nrs, *Fr, temp_r, 1);			// SFT
}



/**********************************************************************************************

					Forward_transform_transpose_order(*Fr) = F 
						*Fr unchanged
   
**********************************************************************************************/

void RSF::RS_Forward_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp_r)
{
	temp_r = *Fr;
	
	Forward_transform_array_transpose_order(RS_basis_type, Nrs, temp_r, F, 1);		// SFT
}



/**********************************************************************************************

					Inplace Inverse Fourier transform
   
**********************************************************************************************/

void RSF::RS_Inverse_transform(Array<complx,3> temp_r)
{
	Inverse_transform_array(RS_basis_type, Nrs, *Fr, temp_r, 1);			// ISFT
}

/**********************************************************************************************

					Inverse_transform_transpose_order(F) = *Fr 
						F unchanged
   
**********************************************************************************************/

void RSF::RS_Inverse_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp)
{
	temp = F;
	
	Inverse_transform_array_transpose_order(RS_basis_type, Nrs, temp, *Fr, 1);		// ISFT
}


/**********************************************************************************************

				Dealias RSF
   
**********************************************************************************************/

void RSF::RS_Dealias()
{
	Dealias_array(RS_basis_type, Nrs, *Fr);
}




/**********************************************************************************************

    Outputs F in real space
   
**********************************************************************************************/

void RSF::RS_Output(ofstream& fileout, Array<complx,3> temp, string format)
{
	Output_asreal(fileout, Nrs, *Fr, temp, format);
}


void RSF::RS_Output_transpose_order(ofstream& fileout, Array<complx,3> temp_r, string format)
{
	Output_asreal_transpose_order(fileout, Nrs, *Fr, temp_r, format);
}

void RSF::RS_input(ifstream& file_in, Array<complx,3> temp_array, string format)
{
	Read_real_field_data_MPI(RS_basis_type, file_in, Nrs, *Fr, temp_array, format);
}


void RSF::RS_structure_fn_Q_ends(int Pi1, int Pi2, int Pi3)
{
	// To optimize, we take max distance to be min(N/5,
	if (RS_structure_fn_approx_switch == 0) {
		Qi1_start = Qi2_start = Qi3_start = 0;
		
		Qi1_end = local_N1-1;
		Qi2_end = Nrs[2]-1;
		Qi3_end = Nrs[3]/2-1;
	}
	
	else if (RS_structure_fn_approx_switch == 1) {
		Qi1_start = Pi1;
		Qi2_start = Pi2;
		Qi3_start = Pi3;
		
		Qi1_end = min(((int) (Pi1 + Nrs[1]/RS_rmax_div_st_fn_r_max)), local_N1-1);
		Qi2_end = min(((int) (Pi2 + Nrs[2]/RS_rmax_div_st_fn_r_max)), Nrs[2]-1);
		Qi3_end = min(((int) (Pi3 + Nrs[3]/(2*RS_rmax_div_st_fn_r_max))), Nrs[3]/2-1);
		// To optimize, we take max distance to be min of the above two.
	}
}

//************************ END of RSF class Definitions ***************************************



