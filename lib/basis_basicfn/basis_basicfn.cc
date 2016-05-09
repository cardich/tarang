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


/*! \file	basis_basicfn.cc
 * 
 * @brief basic functions definitions that are common to all basis functions.
 *
 * @sa basis_basicfn.h
 *
 * @version 4.0 Parallel version
 * @author  M. K. Verma
 * @date	Sept 2008
 * @bug		The Transpose for 2D is not working at present.. it must transport 
 *				real array.
 */ 


#include "basis_basicfn.h"

using namespace blitz ;

//
//
//*********************************************************************************************
void Transpose_array(int N[], Array<complx,3> A, Array<complx,3> Atr)
{
	
	if (N[2] > 1) {

		if (numprocs == 1)
			Atr = A.transpose(1,0,2);
		
		else {
			static Array<complx,2> Axy(local_N1, N[2]);	
			static Array<complx,2> Axy_tr(N[2], local_N1);
			
			static Array<complx,2> temp2(N[2], local_N1);
			
			int data_size = 2* local_N1 * local_N2;								
			// 2 for complex to double
			
			for (int i3=0; i3<=N[3]/2; i3++) {
				
				Axy = A(Range::all(), Range::all(), i3);
				Axy_tr = Axy.transpose(1,0);
				
				// The jth block sent from process i is received by process j and 
				// is placed in the ith block of recvbuf. 
				
				MPI_Alltoall(reinterpret_cast<DP*>(Axy_tr.data()), data_size,  MPI_DP, 
							 reinterpret_cast<DP*>(temp2.data()),
							 data_size,   MPI_DP, MPI_COMM_WORLD);	
				
				for (int source = 0; source < numprocs; source++) 
					Atr(Range::all(), Range(source*local_N1,(source+1)*local_N1-1), i3)
					= temp2(Range(source*local_N2,(source+1)*local_N2-1), Range::all());
			}
		}
	}
	
	// 2D case
	else if (N[2] == 1) {
		
		if (numprocs == 1) 
			Atr = A.transpose(2,1,0);
		
		else {
			static Array<complx,2> Apart(local_N1, N[3]/2);	
			static Array<complx,2> Apart_tr(N[3]/2, local_N1);
			
			static Array<complx,2> temp2(N[3]/2, local_N1);
			
			int data_size = 2* local_N1 * local_N3;								
			// 2 for complex to double
			
			Apart = A(Range::all(), 0, Range(0,N[3]/2-1));
			Apart_tr = Apart.transpose(1,0);
			
			MPI_Alltoall(reinterpret_cast<DP*>(Apart_tr.data()), data_size, MPI_DP, 
						 reinterpret_cast<DP*>(temp2.data()),
						 data_size,  MPI_DP, MPI_COMM_WORLD);
			
			for (int source = 0; source < numprocs; source++)									
				Atr(Range(0,local_N3-1), 0, Range(source*local_N1,(source+1)*local_N1-1))
				= temp2(Range(source*local_N3,(source+1)*local_N3-1), Range::all());
			
			
			// Taking care of the kz=N[3]/2 row.  It is in the last processor.
			static Array<complx,1>	temp_array(local_N1);
			data_size = 2*local_N1;
			int tag = 999;
			
			if (my_id == (numprocs-1)) {
				Atr(local_N3, 0, Range(my_id*local_N1,(my_id+1)*local_N1-1)) 
				= A(Range::all(), 0, N[3]/2);		// my_id = last proc
				
				for (int source = 0; source < (numprocs-1); source++) {
					MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size, 
							 MPI_DP, source, tag, MPI_COMM_WORLD, &status );
					
					Atr(local_N3, 0, Range(source*local_N1,(source+1)*local_N1-1))
					= temp_array(Range::all());
				}	
			}
			
			else {
				temp_array = A(Range::all(),0,N[3]/2);
				
				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size, 
						 MPI_DP, numprocs-1, tag, MPI_COMM_WORLD );
				// send to last processor		 
			}
			
			MPI_Barrier(MPI_COMM_WORLD);

		}
	}
}

//
//
//*********************************************************************************************


void Inverse_transpose_array(int N[], Array<complx,3> Atr, Array<complx,3> A)
{

	if (N[2] > 1) {
		
		if (numprocs == 1)
			A = Atr.transpose(1,0,2);
		
		else{
			static Array<complx,2> Axy(N[1], local_N2);
			static Array<complx,2> Atr_xy(local_N2, N[1]);	
			
			static Array<complx,2> temp1(N[1], local_N2);
			
			
			
			
			int data_size = 2* local_N1 * local_N2;								
			// 2 for complex to double
			
			for (int i3=0; i3<=N[3]/2; i3++)
			{
				
				Atr_xy = Atr(Range::all(), Range::all(), i3);
				Axy = Atr_xy.transpose(1,0);
				
				// The jth block sent from process i is received by process j and 
				// is placed in the ith block of recvbuf. 
				
				MPI_Alltoall(reinterpret_cast<DP*>(Axy.data()), data_size, MPI_DP, 
							 reinterpret_cast<DP*>(temp1.data()),
							 data_size,  MPI_DP, MPI_COMM_WORLD);	
				
				for (int source = 0; source < numprocs; source++) 
					A(Range::all(), Range(source*local_N2,(source+1)*local_N2-1), i3)
					= temp1(Range(source*local_N1,(source+1)*local_N1-1), Range::all());
			}	
		}
	}
	
	// 2D case
	else if (N[2] == 1) {
		
		if (numprocs == 1) 
			A = Atr.transpose(2,1,0);
		
		else {
			static Array<complx,2> Apart_tr(local_N3, N[1]);	
			static Array<complx,2> Apart(N[1], local_N3);
			
			static Array<complx,2> temp1(N[1], local_N3);
			
			int data_size = 2* local_N1 * local_N3;								
			// 2 for complex to double
			
			Apart_tr = Atr(Range(0,local_N3-1), 0, Range::all());
			Apart = Apart_tr.transpose(1,0);
			
			MPI_Alltoall(reinterpret_cast<DP*>(Apart.data()), data_size, MPI_DP, 
						 reinterpret_cast<DP*>(temp1.data()),
						 data_size,  MPI_DP, MPI_COMM_WORLD);
			
			for (int source = 0; source < numprocs; source++)												
				A(Range::all(), 0, Range(source*local_N3,(source+1)*local_N3-1))
				= temp1(Range(source*local_N1,(source+1)*local_N1-1), Range::all());			 
			
			// Taking care of the kz=N[3]/2 row.  It is in the last processor.
			static Array<complx,1>	temp_array(local_N1);
			data_size = 2*local_N1;
			int tag = 999;
			
			if (my_id == (numprocs-1)) {
				A(Range::all(), 0, N[3]/2) 
				= Atr(local_N3, 0, Range(my_id*local_N1, (my_id+1)*local_N1-1));		
				// my_id = last proc,
				
				for (int dest = 0; dest < (numprocs-1); dest++) {
					temp_array = Atr(local_N3, 0, Range(dest*local_N1, (dest+1)*local_N1-1));
					
					MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size, 
							 MPI_DP, dest, tag, MPI_COMM_WORLD);
					
				}		 
			}
			
			else {
				MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size, 
						 MPI_DP, numprocs-1, tag, MPI_COMM_WORLD, &status);
				
				A(Range::all(), 0, N[3]/2) = temp_array;	
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
		
}

//
//
//*********************************************************************************************


int Get_shell_index(DP kkmag, Array<DP, 1> shell_radius_array)
{
	int index;
	
	for (index=0; index <= (shell_radius_array.length())(0); index++)
	{
		if (kkmag <= shell_radius_array(index))
			break;
    }
		
	return index;
}

//
//*********************************************************************************************

int Get_sector_index(DP theta, Array<DP, 1> sector_angle_array)
{
	int index;
	
	if (theta > M_PI/2)
		theta = M_PI - theta;		// 0 <= theta <= M_PI/2 
	
	if (abs(theta) < MYEPS)			// theta = 0 line included in sector=1
		return 1;
		
	else
	{	
		for (index=0; index <= (sector_angle_array.length())(0); index++)
			if (theta <= sector_angle_array(index))
			break;
		
		return index;
	}

	// The special case theta=pi belongs to the last sector.
}

//
//*********************************************************************************************

void Compute_ring_index
(
	DP kkmag, 
	DP theta, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	int& shell_index, 
	int& sector_index
)
{
	int ind1, ind2;
	
	for (ind1=0; ind1 <= (shell_radius_array.length())(0); ind1++)
	{
		if (kkmag <= shell_radius_array(ind1))
			break;
    }
	
	shell_index = ind1;
	
	
	// sector index 
	
	if (theta > M_PI/2)
		theta = M_PI - theta;		// 0 <= theta <= M_PI/2
	
	if (abs(theta) < MYEPS)			// theta = 0 line included in sector=1
		sector_index = 1;
	
	else	
	{
		for (ind2=0; ind2 <= (sector_angle_array.length())(0); ind2++)
			if (theta <= sector_angle_array(ind2))
				break;
	
		sector_index = ind2;
		// The special case theta=pi belongs to the last sector.
	}	
}


//
//*********************************************************************************************

int Get_slab_index(DP kkpll, DP kkperp, Array<DP,1> cylinder_kpll_array)
{
	int index;
	
	kkpll = abs(kkpll);			// if kkpll < 0, include it in the positive slab.
	
	// kz = 0 included in the first slab
	if (kkpll < MYEPS )
		return 1;
	
	else
	{
		for (index=1; index <= (cylinder_kpll_array.length())(0); index++)
			if (kkpll <= cylinder_kpll_array(index))
				break;
		
		return index;
	}

	
	/*
		// kz = 0 included in the first slab
	if ( (kkpll >= cylinder_kpll_array(0)) && (kkpll <= cylinder_kpll_array(1)) )
		return 1;
		
	else
	{
		for (index=2; index <= (cylinder_kpll_array.length())(0); index++)
			if (kkpll <= cylinder_kpll_array(index))
				break;
		
		return index;
	}
	 */
	
}

//
//*********************************************************************************************

void Compute_cylinder_ring_index
(
	DP kkpll, 
	DP kkperp, 
	Array<DP,1> cylinder_shell_radius_array,
	Array<DP,1> cylinder_kpll_array, 
	int& shell_index, 
	int& slab_index
)	
{

	int ind1, ind2;
	
	for (ind1=0; ind1 <= (cylinder_shell_radius_array.length())(0); ind1++)
		if (kkperp <= cylinder_shell_radius_array(ind1))
			break;
	
	shell_index = ind1;
	
	// for slab index
	
	kkpll = abs(kkpll);			// if kkpll < 0, include it in the positive slab.
	
	
	// kz = 0 included in the first slab
	if (kkpll < MYEPS )
		slab_index = 1;
	
	else
	{
		for (ind2=1; ind2 <= (cylinder_kpll_array.length())(0); ind2++)
			if (kkpll <= cylinder_kpll_array(ind2))
				break;
		
		slab_index = ind2;
	}
	
	/*
	// kz = 0 included in the first slab
	if ( (kkpll >= cylinder_kpll_array(0)) && (kkpll <= cylinder_kpll_array(1)) )
		slab_index = 1;
		
	else
	{																						
		for (ind2=2; ind2 <= (cylinder_kpll_array.length())(0); ind2++)
			if (kkpll <= cylinder_kpll_array(ind2))
				break;
		
		slab_index = ind2;
	}	
	 */
}


//*********************************************************************************************

// Given i_last index returns k_last index (Fourier space). 
// If i_last is even, then take the real part of the Fourier space value.
// else take the imaginary part of the value.
void Real_to_fourier_index(int i_last, int& k_last, int& is_real)
{
	k_last = i_last/2;
	
	if ((i_last % 2) == 0)
		is_real = 1;
	
	else
		is_real = 0;

}


//***********************************  End of basis_basicfn.cc ********************************



