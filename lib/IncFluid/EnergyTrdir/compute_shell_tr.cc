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

/*! \file  compute_shell_tr.cc
 * 
 * @brief  Computes shell-to-shell transfer between V/W to V/W.
 *
 *	The Giver field is filled in shell  (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 *  
 *  The shell-from-index = 1:no-shells-1.
 *   shell-mult-all() provides us shell-to-shell energy transfer from these shells to
 *	 1:no_shells.  We output shell-to-shell(1:no-shells-1, 1:no-shells).
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 

#include "../IncFluid.h"


//*********************************************************************************************

// Fluid
//
void IncFluid::Compute_shell_tr()
{
	
	(*shelltoshell_self) = 0.0;
	
	
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{	
		Fill_shell(shell_from_index);	
		
		EnergyTr_Compute_nlin();												
		// nlin = U.grad Um	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, *V1, *V2, *V3, 
						*shell_radius, *temp_shell_tr, kfactor);
		// results(shell_index) in (*temp_shell_tr)(index)
				
		(*shelltoshell_self)(shell_from_index, Range::all()) = -*temp_shell_tr;
			
	}
}

//
//	SCALAR   
//

void IncFluid::Compute_shell_tr(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Compute_shell_tr_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Compute_shell_tr_RB(T);
}


void IncFluid::Compute_shell_tr_scalar(IncSF& T)
{

	// U to U
	Compute_shell_tr();
	
	(*shelltoshell_SF) = 0.0;
	
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{	
		Fill_shell( shell_from_index, T );	
		
		EnergyTr_Compute_nlin(T);												
		// nlin1 = U.grad Tm	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *shell_radius, 
						*temp_shell_tr, kfactor);	
		
		(*shelltoshell_SF)(shell_from_index, Range::all()) = -*temp_shell_tr;	
			
	}
}

// RB
void IncFluid::Compute_shell_tr_RB(IncSF& T)
{
	
	(*shelltoshell_self) = 0.0;
	(*shelltoshell_SF) = 0.0;
	
	if (globalvar_Pr_switch == "PRZERO") 
		Compute_shell_tr();
	
	else if (globalvar_Pr_switch == "PRINFTY")		// fill only Temperature flux
	{
		(*shelltoshell_SF) = 0.0;
		
		for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
		{
			Fill_shell( shell_from_index, T );	
			
			EnergyTr_Compute_nlin(T);										
			// nlin1 = U.grad Tm	
			
			Shell_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *shell_radius, 
						   *temp_shell_tr, kfactor);		
			
			(*shelltoshell_SF)(shell_from_index, Range::all()) = -*temp_shell_tr;	
			
		}
	}
	
	else
		Compute_shell_tr_scalar(T);
}

//*********************************************************************************************
//	 MHD   
//
void IncFluid::Compute_shell_tr(IncVF& W)
{
	// U to U
	Compute_shell_tr();
	
	(*shelltoshell_VF_WtoW) = 0.0;
	
	
	// W to W
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{
		Fill_shell(shell_from_index, W);	
		
		EnergyTr_Compute_nlin();											
		// nlin = U.grad Wm	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
						*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);		

		(*shelltoshell_VF_WtoW)(shell_from_index, Range::all()) = -*temp_shell_tr;	
		
	}


	// U to W
	
	(*shelltoshell_VF_UtoW) = 0.0;
	

	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{	
		Fill_shell(shell_from_index);	
		
		EnergyTr_Compute_nlin(W);											
		// nlin = W.grad Um	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
						*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);		
						
		(*shelltoshell_VF_UtoW)(shell_from_index, Range::all()) = *temp_shell_tr;
			

	}
	
	
	// Shell_to_shell transfers for Elsasser vars
	//
	
	(*shelltoshell_Elsasser_plus) = 0.0;
	(*shelltoshell_Elsasser_minus) = 0.0;



	IncFluid::UB_to_Elsasser_field(W);											
	// U=Zp=(U+B); B=Zm=(U-B);
	
	// Shell_to_shell: Zp to Zp
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{	
		Fill_shell(shell_from_index);										
		// G = Zp<
		
		EnergyTr_Compute_nlin(W);											
		// nlin = Zm.grad Zp<	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, *V1, *V2, *V3, 
						*shell_radius, *temp_shell_tr, kfactor);	
							
		(*shelltoshell_Elsasser_plus)(shell_from_index, Range::all()) = -*temp_shell_tr;	
	}

	
	// Shell_to_shell: Zm to Zm
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{
		Fill_shell(shell_from_index, W);									
		// G = Zm<
		
		EnergyTr_Compute_nlin();											
		// nlin = Zp.grad Zm<	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *nlin2, *nlin3, 
						*W.V1, *W.V2, *W.V3, *shell_radius, *temp_shell_tr, kfactor);		
						
		(*shelltoshell_Elsasser_minus)(shell_from_index, Range::all()) = -*temp_shell_tr;	
	}

	IncFluid::Elsasser_to_UB_field(W);													
	// Back to U, B vars
	
}



//*********************************************************************************************
// Vector + Scalar
//
void IncFluid::Compute_shell_tr(IncVF& W, IncSF& T)
{
	
	// U/W to U/W
	Compute_shell_tr(W);
	
	(*shelltoshell_SF) = 0.0;

	
	
	// T to T
	for (int shell_from_index = 1; shell_from_index < no_shells; shell_from_index++) 
	{	
		Fill_shell( shell_from_index, T);	
		
		EnergyTr_Compute_nlin(T);												
		// nlin = U.grad Tm	
		
		Shell_mult_all(basis_type, alias_switch, N, *nlin1, *T.F, *shell_radius, 
						*temp_shell_tr, kfactor);	
							
		(*shelltoshell_SF)(shell_from_index, Range::all()) = -*temp_shell_tr;	
		
	}
		
}	


//******************************  End of Compute_shell_tr.cc  *********************************


