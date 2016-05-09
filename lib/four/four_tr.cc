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
/*! \file fourier.cc 
 * 
 * @sa fourier.h
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		none known
 */ 


#include "four_tr.h"

using namespace blitz ;


//*********************************************************************************************

void Init_fftw_plan_FOUR(int N[], Array<complx,3> A)
{
	
	if (globalvar_fftw_original_switch == 1) {
		
		if (N[2] > 1) {
			r2c_plan_FOUR = FFTW_MPI_PLAN_DFT_R2C_3D_DP(N[1], N[2], N[3], 
								reinterpret_cast<DP*>(A.data()), 
								reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()),  
								MPI_COMM_WORLD, FFTW_PLAN_FLAG);
									
			c2r_plan_FOUR = FFTW_MPI_PLAN_DFT_C2R_3D_DP(N[1], N[2], N[3], 
								reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), 
								reinterpret_cast<DP*>(A.data()), 
								MPI_COMM_WORLD, FFTW_PLAN_FLAG);	
		}	
		
		else if (N[2] ==1) {
			static Array<complx,2> A2d(N[1], (N[3]/2)+1);
			
			r2c_plan_FOUR = FFTW_MPI_PLAN_DFT_R2C_2D_DP(N[1], N[3], 
										 reinterpret_cast<DP*>(A2d.data()), 
										 reinterpret_cast<FFTW_COMPLEX_DP*>(A2d.data()),  
										 MPI_COMM_WORLD, FFTW_PLAN_FLAG);
			
			c2r_plan_FOUR = FFTW_MPI_PLAN_DFT_C2R_2D_DP(N[1], N[3], 
										 reinterpret_cast<FFTW_COMPLEX_DP*>(A2d.data()), 
										 reinterpret_cast<DP*>(A2d.data()), 
										 MPI_COMM_WORLD, FFTW_PLAN_FLAG);
		}
	}
	// 3D using 2D & 1D transforms
	else {
		if (N[2] > 1) {
			
			Array<complx,1> temp_row(N[2]);
			Array<complx,2> temp_plane(N[1],(N[3]/2)+1);
			
			r2c_2d_plan_FOUR = FFTW_PLAN_DFT_R2C_2D_DP(N[1], N[3], 
										reinterpret_cast<DP*>(temp_plane.data()),
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_plane.data()),
										FFTW_PLAN_FLAG);
			
			c2r_2d_plan_FOUR = FFTW_PLAN_DFT_C2R_2D_DP(N[1], N[3], 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_plane.data()),
										reinterpret_cast<DP*>(temp_plane.data()),
										FFTW_PLAN_FLAG);
			
			c2c_1d_forward_plan_FOUR = FFTW_PLAN_DFT_1D_DP(N[2], 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row.data()), 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row.data()), 
										FFTW_FORWARD, FFTW_PLAN_FLAG);
			
			c2c_1d_inverse_plan_FOUR = FFTW_PLAN_DFT_1D_DP(N[2], 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row.data()), 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row.data()), 
										FFTW_BACKWARD, FFTW_PLAN_FLAG);
		}	
		
		else if (N[2] ==1) {
			Array<complx,1> temp_row2(N[1]);
			Array<complx,1> temp_col2(N[3]/2+1);
			
			c2c_1d_forward_plan_FOUR = FFTW_PLAN_DFT_1D_DP(N[1], 
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row2.data()),
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row2.data()),
										FFTW_FORWARD, FFTW_PLAN_FLAG);
			
			c2c_1d_inverse_plan_FOUR = FFTW_PLAN_DFT_1D_DP(N[1],  
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row2.data()),
										reinterpret_cast<FFTW_COMPLEX_DP*>(temp_row2.data()),
										FFTW_BACKWARD, FFTW_PLAN_FLAG);
			
			r2c_1d_plan_FOUR = FFTW_PLAN_DFT_R2C_1D_DP(N[3], 
									   reinterpret_cast<DP*>(temp_col2.data()), 
									   reinterpret_cast<FFTW_COMPLEX_DP*>(temp_col2.data()), 
									   FFTW_PLAN_FLAG);
			
			c2r_1d_plan_FOUR = FFTW_PLAN_DFT_C2R_1D_DP(N[3], 
									   reinterpret_cast<FFTW_COMPLEX_DP*>(temp_col2.data()), 
									   reinterpret_cast<DP*>(temp_col2.data()), 
									   FFTW_PLAN_FLAG);
		}
	}
}


//*********************************************************************************************

void FTr2c_plane_FOUR(int N[], Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_2d_plan_FOUR, reinterpret_cast<DP*>(Plane.data()), 
						 reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()));
}

void FTc2r_plane_FOUR(int N[], Array<complx,2> Plane)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_2d_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(Plane.data()), 
						 reinterpret_cast<DP*>(Plane.data()));
}

void FT_col_FOUR(int N[], Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_1d_forward_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), 
					 reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void IFT_col_FOUR(int N[], Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_DP(c2c_1d_inverse_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), 
					 reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void FTr2c_col_FOUR(int N[], Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_R2C_DP(r2c_1d_plan_FOUR, reinterpret_cast<DP*>(Col.data()), 
						reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()));
}

void FTc2r_col_FOUR(int N[], Array<complx,1> Col)
{
	FFTW_EXECUTE_DFT_C2R_DP(c2r_1d_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(Col.data()), 
						reinterpret_cast<DP*>(Col.data()));
}

//*********************************************************************************************

void Zero_pad_last_plane_FOUR(int N[],  Array<complx,3> A)
{
	A(Range(0,local_N1-1),Range(0,N[2]-1),N[3]/2) = 0.0;
}

//*********************************************************************************************

void Norm_FOUR(int N[], Array<complx,3> A) 
{
	A = A/(DP(N[1]) *  DP(N[2]) * DP(N[3]));
}

//**************************************************************************************************

void ArrayFFT_FOUR_transpose_order
(
	int N[], 
	Array<complx,3> Atr, 
	Array<complx,3> A
)
{
	if (N[2] > 1)  {
		
		Atr(Range(0,local_N2-1),Range(0,N[1]-1),N[3]/2) = 0.0;	// zero_pad k=N[3]/2
		
		static Array<complx,1> temp_col(N[2]);					// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[1], (N[3]/2)+1);	// Temp location to hold *Vi(i,:)
		
		for (int l1=0; l1<local_N2; l1++) {
			temp_plane = Atr(l1, Range::all(), Range::all());
			FTr2c_plane_FOUR(N, temp_plane);
			Atr(l1, Range::all(), Range::all()) = temp_plane;
		}
		
		Inverse_transpose_array(N, Atr, A);
		
		for (int l1=0; l1<local_N1; l1++) 
			for (int l3=0; l3<=N[3]/2; l3++) {
				temp_col = A(l1, Range::all(), l3);
				FT_col_FOUR(N, temp_col);
				A(l1, Range::all(), l3) = temp_col;
			}
		
		Norm_FOUR(N, A);
	}
	
	else if (N[2] == 1) {
		cout << "ERROR:  ArrayFFT_FOUR_transpose_order for 2D array not" 
			 <<		"allowed in FOUR basis " << endl;
		exit(1);
	}
	
}


//**************************************************************************************************

void ArrayFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
) 
{
	
	if (globalvar_fftw_original_switch == 1) {
		
		Zero_pad_last_plane_FOUR(N, A);
		
		if (N[2] > 1) {	
			FFTW_EXECUTE_DFT_R2C_DP(r2c_plan_FOUR, reinterpret_cast<DP*>(A.data()), 
								 reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()));
		}
		
		else if (N[2] == 1) {
			
			static Array<complx,2> A2d(N[1], (N[3]/2)+1); 
	
			A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
			
			FFTW_EXECUTE_DFT_R2C_DP(r2c_plan_FOUR, reinterpret_cast<DP*>(A2d.data()), 
				 reinterpret_cast<FFTW_COMPLEX_DP*>(A2d.data())); 
			
			A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
			
		}
		
		Norm_FOUR(N, A);
	}
	
	// not original switch.. split in 2D and 1D ffts.
	else {
		if (N[2] > 1) {
			Transpose_array(N, A, temp_r);
			ArrayFFT_FOUR_transpose_order(N, temp_r, A);	// Norm done here
		}
		
		else if (N[2] == 1) {
			static Array<complx,1> temp_row(N[1]);			// Temp location to hold *Vi(:,j,k)
			static Array<complx,1> temp_col(N[3]/2+1);	
			
			A(Range::all(),0,N[3]/2) = 0.0;					// zero_pad k=N[3]/2
			
			for (int l1=0; l1 < local_N1; l1++)	{
				temp_col = A(l1, 0, Range::all());
				FTr2c_col_FOUR(N, temp_col);				// r2c
				A(l1, 0, Range::all()) = temp_col;
			}
			
			Transpose_array(N, A, temp_r);				
			
			for (int l3=0; l3<local_N3; l3++) {
				temp_row = temp_r(l3, 0, Range::all());
				FT_col_FOUR(N, temp_row);	//c2c
				temp_r(l3, 0, Range::all()) = temp_row;
			}
			
			if (my_id == numprocs-1) {
				temp_row = temp_r(local_N3, 0, Range::all());
				FT_col_FOUR(N, temp_row);
				temp_r(local_N3, 0, Range::all()) = temp_row;
			}
			
			Inverse_transpose_array(N, temp_r, A);		
			
			Norm_FOUR(N, A);
		}
		
	}
	
}




//*********************************************************************************************

/*
void ArrayFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
)
{
	Zero_pad_last_plane_FOUR(N, A); 
	ArrayFFTW_FOUR(N, A, temp_r); 
	Norm_FOUR(N, A); 
}


void ArrayFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
)
{
	Zero_pad_last_plane_FOUR(N, A); 
	ArrayFFTW_FOUR_transpose_order(N, Atr, A); 
	Norm_FOUR(N, A); 
}
 */

//*********************************************************************************************


void ArrayIFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> A, 
 Array<complx,3> Atr
 )
{
	if (N[2] > 1)
	{
		
		static Array<complx,1> temp_col(N[2]);					// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[1], (N[3]/2)+1);    // Temp location to hold *Vi(i,:)
		
		for (int l1=0; l1<local_N1; l1++) 
			for (int l3=0; l3<=N[3]/2; l3++) {
				temp_col = A(l1, Range::all(), l3);
				IFT_col_FOUR(N, temp_col);
				A(l1, Range::all(), l3) = temp_col;
			}
		
		Transpose_array(N, A, Atr);
		
		for (int l1=0; l1 < local_N2; l1++)	{
			temp_plane = Atr(l1, Range::all(), Range::all());
			FTc2r_plane_FOUR(N, temp_plane);
			Atr(l1, Range::all(), Range::all()) = temp_plane;
		}
		
		Zero_pad_last_plane_FOUR(N, Atr);
	}
	
	else if (N[2] == 1) {
		cout << "ERROR:  ArrayIFFT_FOUR_transpose_order for 2D array not" 
			 <<		"allowed in FOUR basis " << endl;
		exit(1);
	}	
}


//**************************************************************************************************


void ArrayIFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
)
{
	
	if (globalvar_fftw_original_switch == 1) {
		
		if (N[2] > 1)
			FFTW_EXECUTE_DFT_C2R_DP(c2r_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(A.data()), 
								 reinterpret_cast<DP*>(A.data()));
		
					
		else if (N[2] == 1) {
			
			static Array<complx,2> A2d(N[1], (N[3]/2)+1);
			
			A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
			

			FFTW_EXECUTE_DFT_C2R_DP(c2r_plan_FOUR, reinterpret_cast<FFTW_COMPLEX_DP*>(A2d.data()), 
								 reinterpret_cast<DP*>(A2d.data()));
			
			A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
		}
	}
	
	// using 2d and 1d ffts
	else {
		if (N[2] > 1) {
			ArrayIFFT_FOUR_transpose_order(N, A, temp_r);
			Inverse_transpose_array(N, temp_r, A);
		}
		
		else if (N[2] == 1) {
			static Array<complx,1> temp_row(N[1]);			// Temp location to hold *Vi(:,j,k)
			static Array<complx,1> temp_col(N[3]/2+1);	
			
			Transpose_array(N, A, temp_r);  
			
			for (int l3=0; l3<local_N3; l3++) {
				temp_row = temp_r(l3, 0, Range::all());
				IFT_col_FOUR(N, temp_row);	//c2c
				temp_r(l3, 0, Range::all()) = temp_row;
			}
			
			if (my_id == numprocs-1) {
				temp_row = temp_r(local_N3, 0, Range::all());
				IFT_col_FOUR(N, temp_row);
				temp_r(local_N3, 0, Range::all()) = temp_row;
			}
			
			Inverse_transpose_array(N, temp_r, A);  
			
			for (int l1=0; l1 < local_N1; l1++)	{
				temp_col = A(l1, 0, Range::all());
				FTc2r_col_FOUR(N, temp_col);	// c2r
				A(l1, 0, Range::all()) = temp_col;
			}
		}
	}
	
}


//*********************************************************************************************

void Xderiv_FOUR(int N[],Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 			// l1 is the local array-index along x
	{
		k1 = Get_kx_FOUR(l1, N);
		
		B(l1,Range::all(),Range::all()) = 
		    complex<DP>(0, kfactor[1])*((DP) 1.0*k1)*( A(l1,Range::all(),Range::all()) ); 	
	}
}	


//*********************************************************************************************


void Yderiv_FOUR(int N[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	secondIndex l2;
	B(Range::all(),Range(0,N[2]/2),Range::all()) = 
		complex<DP>(0, kfactor[2])*((DP) 1.0*l2)
		* ( A(Range::all(),Range(0,N[2]/2),Range::all()) ); 
		// l2 = 0:N2/2; k2 = l2.
	
	if (N[2] > 1)
		B(Range::all(),Range(N[2]/2+1,N[2]-1),Range::all()) = 
			complex<DP>(0, kfactor[2])*((DP) 1.0*(l2+1-N[2]/2))*
			( A(Range::all(),Range(N[2]/2+1,N[2]-1),Range::all()) );
			// l2 = 0:N2/2;  k2 = l2-N2/2+1.
}

//*********************************************************************************************


void Zderiv_FOUR(int N[],Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	thirdIndex l3;
	B = complex<DP>(0, kfactor[3])*((DP) 1.0*l3)*(A);
}

//********************************	End of four_tr.cc *****************************************


