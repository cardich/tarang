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

/*! \file  add_pressure_grad.cc
 * 
 * @brief  Add pressure gradient to nlin.   
 * 
 * @note nlin[i] contains Dj T[Vr[j]*Vr[i]];
 *		F.p contains pressure;
 *		*p is Cos transformed
 *
 * @note Result => nlin[i] = nlin[i] +  Di [p(k)] (grad(pressure)) 
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "IncFluid.h"


//*********************************************************************************************

void IncFluid::Add_pressure_gradient()       
{	

	DP q_omega_t = globalvar_Keplerian_q_shear* globalvar_Keplerian_omega* globalvar_Tnow;
	
	// parity  =0; -dp/dx= f (sin)
	Xderiv(basis_type, N, *F, *VF_temp, kfactor, 0);   
	*nlin1 = *nlin1 + *VF_temp; 
	if (globalvar_prog_kind == "KEPLERIAN") {
		Yderiv(basis_type, N, *F, *VF_temp, kfactor, 0); 
		*nlin1 = *nlin1 + complex<DP>(q_omega_t,0)*(*VF_temp);
	}
		
	
	Yderiv(basis_type, N, *F, *VF_temp, kfactor, 0);   
	*nlin2 = *nlin2 + *VF_temp;
	
	Zderiv(basis_type, N, *F, *VF_temp, kfactor, 0);   
	*nlin3 = *nlin3 + *VF_temp;
		
}

//******************************* End of add_pressure_grad.cc  ********************************

