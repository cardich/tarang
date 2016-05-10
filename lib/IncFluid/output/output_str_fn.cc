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

/*! \file  Output_str_fn.cc
 * 
 * @brief  Output Structure function, Planar structure function.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"   

//*********************************************************************************************


void IncFluid::Output_structure_fn()
{
	if (my_id == master_id) {
		structure_fn_file << "%% Time = " << Tnow << endl; 
		structure_fn_file << "%% q_min, qmax, r_arraysize = " 
			<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
			<< RV_structurefn_r_arraysize << endl; 
		structure_fn_file << "%% structurefn_r_max = " << RV_structurefn_r_max << endl << endl; 
	}
		
	
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	//Now *Vir contains real space vars
	
	// if 2D, then ignore V2 transfer.
	if (N[2] == 1)   *V2r = 0.0;
	
	RV_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
			structure_fn_file	<< r				<< "  ";
			
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
				structure_fn_file	<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
									<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
								
			structure_fn_file << endl;
		}
	
//	if (my_id == 0)
//		cout << "OUTPUT " << (*RV_St_final)(Range::all(), Range::all(), 0) << endl;
}
  
  
//
//

void IncFluid::Output_structure_fn(IncSF& T)
{
	if (my_id == master_id) {
		structure_fn_file << "%% Time = " << Tnow << endl; 
		structure_fn_file << "%% q_min, qmax, r_arraysize = " 
		<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
		<< RV_structurefn_r_arraysize << endl; 
		structure_fn_file << "%% structurefn_r_max = " << RV_structurefn_r_max << endl << endl; 
	}
	
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
	//Now *Vir, *T.Fr contains real space vars
	
	// if 2D, then ignore V2 transfer.
	if (N[2] == 1)   *V2r = 0.0;
	
	
	RV_Compute_structure_function();
	T.RS_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
			structure_fn_file	<< r					<< "  "; 
			
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
				structure_fn_file	<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
									<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " "
									<< (*T.RS_St_final)(r,q-RV_structurefn_q_min)	<< "   ";
								
			structure_fn_file << endl;
		}	
}
    
//
//

void IncFluid::Output_structure_fn(IncVF& W)
{
	if (my_id == master_id) {
		structure_fn_file << "%% Time = " << Tnow << endl; 
		structure_fn_file << "%% q_min, qmax, r_arraysize = " 
		<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
		<< RV_structurefn_r_arraysize << endl; 
		structure_fn_file << "%% structurefn_r_max = " << RV_structurefn_r_max << endl << endl; 
	} 
		
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3; 
	
	W.RV_Inverse_transform(*VF_temp_r);
	//Now *Vir, *W.Vir contains real space vars
	
	// if 2D, then ignore V2 transfer.
	if (N[2] > 1)   {
		*V2r = 0.0;
		*W.V2r = 0.0;
	}
	
	
	RV_Compute_structure_function();
	W.RV_Compute_structure_function();
	
	if (my_id == master_id)
		for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
			structure_fn_file	<< r					<< "  "; 
			
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
				structure_fn_file	<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
									<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " "
									<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
									<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,1)	<< "    ";
								
			structure_fn_file << endl;
		}		
} 

//
//

void IncFluid::Output_structure_fn(IncVF& W, IncSF& T)
{

	/*
	if (my_id == master_id) {
		structure_fn_file << "%% Time = " << Tnow << endl; 
		structure_fn_file << "%% q_min, qmax, r_arraysize = " 
		<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
		<< RV_structurefn_r_arraysize << endl; 
		structure_fn_file << "%% structurefn_r_max = " << RV_structurefn_r_max << endl << endl; 
	}
		
	 // Inverse transform and get real space
	 *V1r = *V1;  
	 *V2r = *V2; 
	 *V3r = *V3; 
	 
	 RV_Inverse_transform(*VF_temp_r);
	 
	 *W.V1r = *W.V1;  
	 *W.V2r = *W.V2;  
	 *W.V3r = *W.V3; 
	 W.RV_Inverse_transform(*VF_temp_r);
 
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
	 //Now *Vir, *W.Vir contains real space vars
 
	RV_Compute_structure_function();
	W.RV_Compute_structure_function();
	T.RS_Compute_structure_function();
	
	
	if (my_id == master_id)
		for (int r=1; r<=RV_structurefn_r_max; r++) {
			structure_fn_file	<< r					<< "  " 
 
			for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
				structure_fn_file	<< (*RV_St)(r,q-RV_structurefn_q_min,0)		<< " " 
									<< (*RV_St)(r,q-RV_structurefn_q_min,1)		<< "    "
									<< (*W.RV_St)(r,q-RV_structurefn_q_min,0)	<< " " 
									<< (*W.RV_St)(r,q-RV_structurefn_q_min,1)	<< "   "
									<< (*T.RS_St)(r,q-RV_structurefn_q_min);
								
			structure_fn_file << endl;
		}	
	*/
} 



//*********************************************************************************************

// Planar structure function
//
void IncFluid::Output_planar_structure_fn()
{
	
	if (my_id == master_id) {
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		planar_structure_fn_file << "%% q_min, qmax, r_arraysize = " 
				<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
				<< RV_structurefn_r_arraysize << endl; 
		planar_structure_fn_file << "%% planar_structurefn_r_max = " << RV_structurefn_r_max << endl;
		planar_structure_fn_file << "%% Planar struncture function: V" << endl << endl;
	}
	
	
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	//Now *Vir contains real space vars

	
//	cout << (*V1r)(Range::all(),0,Range::all()) << (*V2r)(Range::all(),1,Range::all())  << endl;
		
	
	if (globalvar_anisotropy_switch == 1) {
		for (int pll_index=0; pll_index<local_N1; pll_index++) 				
			RV_Compute_planar_structure_function_anis1(planar_structure_fn_file, pll_index);
			// Computes Structure function for each y-z plane (indexed by pll_index), and 
			// print in file_out. No need of any data transfer across processors.
	}
		
	
	
	else if (globalvar_anisotropy_switch == 2) 
		for (int pll_index=0; pll_index<N[2]; pll_index++) {
			RV_Compute_planar_structure_function_anis2(pll_index);
			// Compute Structure function for each x-z plane (indexed by pll_index).  
			// Each processor computes for a subpart of the plane; 
			// these subparts are summed in the master node.
			
			planar_structure_fn_file  << "%% pll_index = " << pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		} 
	
	// Compute Structure function for each x-y plane (indexed by pll_index).  
	// Each processor computes for a subpart of the plane; 
	// these subparts are summed in the master node.
	// Real and imag parts planes are done according to the imag_switch.		
	else if (globalvar_anisotropy_switch == 3) 
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
		
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index << endl;
			RV_Compute_planar_structure_function_anis3(pll_index, 0);
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
			
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index+1 << endl;	
			RV_Compute_planar_structure_function_anis3(pll_index, 1);
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}	
		}
		
}
  
//
//


void IncFluid::Output_planar_structure_fn(IncSF& T)
{
	if (my_id == master_id) {
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		planar_structure_fn_file << "%% q_min, qmax, r_arraysize = " 
			<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
			<< RV_structurefn_r_arraysize << endl; 
		planar_structure_fn_file << "%% planar_structurefn_r_max = " << RV_structurefn_r_max << endl;
		planar_structure_fn_file	<< "%% Planar struncture function" << endl << endl; 
	}
	
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	
	*T.Fr = *T.F;   
	T.RS_Inverse_transform(*VF_temp_r);
	//Now *Vir, *T.Fr contains real space vars
	
	if (globalvar_anisotropy_switch == 1) 
		for (int pll_index=0; pll_index<local_N1; pll_index++) {
			RV_Compute_planar_structure_function_anis1(planar_structure_fn_file, pll_index);
			T.RS_Compute_planar_structure_function_anis1(planar_structure_fn_file, pll_index);	
			// Computes Structure function for each y-z plane (indexed by pll_index), and 
			// print in file_out. No need of any data transfer across processors.
		}
	
	else if (globalvar_anisotropy_switch == 2)	{
		
		for (int pll_index=0; pll_index<N[2]; pll_index++) {
			RV_Compute_planar_structure_function_anis2(pll_index);
			// Compute Structure function for each x-z plane (indexed by pll_index).  
			// Each processor computes for a subpart of the plane; 
			// these subparts are summed in the master node.
			
			planar_structure_fn_file  << "%% pll_index = " << pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
						<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
						<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		} 
		
		
		// Same as above but for T field
		for (int pll_index=0; pll_index<N[2]; pll_index++) {
			T.RS_Compute_planar_structure_function_anis2(pll_index);
			
			planar_structure_fn_file  << "%% T: pll_index = " << pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r << "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file 
							<< (*T.RS_St_final)(r,q-RV_structurefn_q_min) << " ";
					
					planar_structure_fn_file << endl;
				}
		}
	}
	
	// Compute Structure function for each x-y plane (indexed by pll_index).  
	// Each processor computes for a subpart of the plane; 
	// these subparts are summed in the master node.
	// Real and imag parts planes are done according to the imag_switch.
	else if (globalvar_anisotropy_switch == 3) {
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			RV_Compute_planar_structure_function_anis3(pll_index, 0);
			
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index << endl;
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			RV_Compute_planar_structure_function_anis3(pll_index, 1);
			
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index+1 << endl;
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			T.RS_Compute_planar_structure_function_anis3(pll_index, 0);
			
			planar_structure_fn_file  << "%% T: pll_index = " << 2*pll_index << endl;
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file 
						<< (*T.RS_St_final)(r,q-RV_structurefn_q_min)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			T.RS_Compute_planar_structure_function_anis3(pll_index, 1);
			
			planar_structure_fn_file  << "%% T: pll_index = " << 2*pll_index+1 << endl;
			
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*T.RS_St_final)(r,q-RV_structurefn_q_min)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
	}										
}


// Planar structure function
//
void IncFluid::Output_planar_structure_fn(IncVF& W)
{
	if (my_id == master_id) {
		planar_structure_fn_file << "%% Time = " << Tnow << endl; 
		planar_structure_fn_file << "%% q_min, qmax, r_arraysize = " 
				<< RV_structurefn_q_min << " " << RV_structurefn_q_max  << " " 
				<< RV_structurefn_r_arraysize << endl; 
		planar_structure_fn_file << "%% planar_structurefn_r_max = " << RV_structurefn_r_max << endl;
		planar_structure_fn_file	<< "%% Planar struncture function" << endl << endl; 
	}
	
	// Inverse transform and get real space
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	RV_Inverse_transform(*VF_temp_r);
	
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2;  
	*W.V3r = *W.V3; 
	
	W.RV_Inverse_transform(*VF_temp_r);
	//Now *Vir, *W.Vir contains real space vars
		
	if (globalvar_anisotropy_switch == 1) {
		for (int pll_index=0; pll_index<local_N1; pll_index++) 
			RV_Compute_planar_structure_function_anis1(planar_structure_fn_file, pll_index);
		
		planar_structure_fn_file << endl << "W fields: " << endl;
		
		for (int pll_index=0; pll_index<local_N1; pll_index++)
			W.RV_Compute_planar_structure_function_anis1(planar_structure_fn_file, pll_index);	
	}	
	
	else if (globalvar_anisotropy_switch == 2)	{
		
		for (int pll_index=0; pll_index<N[2]; pll_index++) {
			RV_Compute_planar_structure_function_anis2(pll_index);
			
			planar_structure_fn_file  << "%% pll_index = " << pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[2]; pll_index++) {
			W.RV_Compute_planar_structure_function_anis2(pll_index);
			
			planar_structure_fn_file  << "%% W: pll_index = " << pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " "
							<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
	}
		
	// Compute Structure function for each x-y plane (indexed by pll_index).  
	// Each processor computes for a subpart of the plane; 
	// these subparts are summed in the master node.
	// Real and imag parts planes are done according to the imag_switch.
	else if (globalvar_anisotropy_switch == 3) {
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			RV_Compute_planar_structure_function_anis3(pll_index, 0);
			
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
							<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}	
		
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			RV_Compute_planar_structure_function_anis3(pll_index, 1);
			
			planar_structure_fn_file  << "%% pll_index = " << 2*pll_index+1 << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
						<< (*RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " " 
						<< (*RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			W.RV_Compute_planar_structure_function_anis3(pll_index, 0);
			
			planar_structure_fn_file  << "%% W: pll_index = " << 2*pll_index << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
						<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " "
						<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
		
		for (int pll_index=0; pll_index<N[3]/2; pll_index++) {
			W.RV_Compute_planar_structure_function_anis3(pll_index, 1);
			
			planar_structure_fn_file  << "%% W: pll_index = " << 2*pll_index+1 << endl;
			
			// master_id has the planar structure fn for plane= pll_index.
			if (my_id == master_id)	
				for (int r=1; r<=RV_structurefn_r_arraysize; r++) {
					planar_structure_fn_file	<< r	<< "  ";
					
					for (int q=RV_structurefn_q_min; q<=RV_structurefn_q_max; q++)
						planar_structure_fn_file	
						<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,0)	<< " "
						<< (*W.RV_St_final)(r,q-RV_structurefn_q_min,1)	<< " ";
					
					planar_structure_fn_file << endl;
				}
		}
				
	}
		
}
  
//
//

void IncFluid::Output_planar_structure_fn(IncVF& W, IncSF& T)
{											
}


//********************************* End of output_str_fn **************************************


