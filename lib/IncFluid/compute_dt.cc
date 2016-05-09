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

/*! \file  compute_dt.cc
 * 
 * @brief Compute dt using CFL criteria.
 *			Returns min(dx/Urms, Tdt_fixed) to avoid higher dt.
 *
 * @sa DP IncFluid::Get_dt()
 * @note	Reference Pope for CFL condition
 *
 * @author  M. K. Verma
 * @version 4.0   MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "IncFluid.h"


//*********************************************************************************************

DP IncFluid::Get_dt()
{

		//	CV_Compute_totalenergy();
	
	DP delta_u = sqrt(Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V1));
	delta_u += sqrt(Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V2));
	delta_u += sqrt(Get_total_energy(CV_basis_type, CV_alias_switch, Ncv, *V3));
	
	DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
	DP dx = (2*M_PI)/kmax;
	DP dt = dx / (delta_u *20);

	MPI_Bcast( &dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	return min(dt, Tdt_fixed);
}	


//
// Scalar: Same as fluid
//

DP IncFluid::Get_dt(IncSF& T)
{
	return Get_dt();
}	

//
// Vector
//

DP IncFluid::Get_dt(IncVF& W)
{
	
	DP dt = Get_dt();
	
	DP kmax = Max_radius_inside(basis_type, alias_switch, N, kfactor);
	DP dx = (2*M_PI)/kmax;
	
	if (my_id == master_id) {
		DP abs_B0 = abs((*W.V1)(0,0,0)) + abs((*W.V2)(0,0,0)) + abs((*W.V3)(0,0,0));
		
		if (abs_B0 > MYEPS)
			dt = min(dt, dx/(abs_B0 * 20));
	}		
		
	MPI_Bcast( &dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
	
	return min(dt, Tdt_fixed);
}	

//
// Vector + scalar
//

DP IncFluid::Get_dt(IncVF& W, IncSF& T)
{

	return Get_dt(W);
}



//**********************************   End of compute_dt.cc  **********************************


