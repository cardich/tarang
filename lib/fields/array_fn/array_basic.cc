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


/*! \file array_basic.cc
 *
 * @sa field_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */


#include "array_basic.h"

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
/**********************************************************************************************

     		C = A.B (assuming A and B are real)
			Term by term multiplication.

***********************************************************************************************/

void Array_real_mult
(
	int N[],
	Array<complx,3> A, Array<complx,3> B,
	Array<complx,3> C
)
{
	real(C) = real(A) * real(B);
	imag(C) = imag(A) * imag(B);
}

void Array_real_divide
(
 int N[],
 Array<complx,3> A, Array<complx,3> B,
 Array<complx,3> C
)
{
	real(C) = real(A) / real(B);
	imag(C) = imag(A) / imag(B);
}


/**********************************************************************************************

     		Output A assuiming it to be real.

***********************************************************************************************/

void Output_asreal
(
	ofstream& file_out,
	int N[],
	Array<complx,3> A,
	Array<complx,3> temp_array,
	string format
)
{
	int data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
	int tag = 123;

	if (my_id == master_id) {
//		file_out << "%% Array " <<  N[1] << " x "  << N[2] << " x " << N[3] << endl;

		for (int i=0; i<local_N1; i++) {
			for (int j=0; j<N[2]; j++)
				for (int k=0; k<N[3]/2; k++)
					file_out << real(A(i,j,k)) << " " << imag(A(i,j,k)) << " " ;
			file_out << endl;
		}

		for (int source = 1; source <= numprocs-1; source++) {
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
						MPI_DP, source, tag, MPI_COMM_WORLD, &status );

			for (int i=0; i<local_N1; i++) {
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<N[3]/2; k++)
						file_out << real(temp_array(i,j,k)) << " " << imag(temp_array(i,j,k))
								 << " " ;
				file_out << endl;
			}
		}

		file_out << endl << endl;
	}

	else		// process is not master
		MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
						MPI_DP, master_id, tag, MPI_COMM_WORLD );


	MPI_Barrier(MPI_COMM_WORLD);


}


void Output_asreal_transpose_order
(
	 ofstream& file_out,
	 int N[],
	 Array<complx,3> A,
	 Array<complx,3> temp_r,
	string format
)
{

	if (N[2] > 1) {
		int data_size = 2* local_N2 * N[1] * (N[3]/2 + 1);
		int tag = 123;

		if (my_id == master_id) {
		//	file_out << "%% Array " <<  N[2] << " x "  << N[1] << " x " << N[3] << endl;

			for (int i=0; i<local_N2; i++) {
				for (int j=0; j<N[1]; j++)
					for (int k=0; k<N[3]/2; k++)
						file_out << real(A(i,j,k)) << " " << imag(A(i,j,k)) << " " ;
				file_out << endl;
			}

			for (int source = 1; source <= numprocs-1; source++) {

				MPI_Recv( reinterpret_cast<DP*>((temp_r).data()), data_size,
						 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				for (int i=0; i<local_N2; i++) {
					for (int j=0; j<N[1]; j++)
						for (int k=0; k<N[3]/2; k++)
							file_out << real(temp_r(i,j,k)) << " " << imag(temp_r(i,j,k)) << " " ;

					file_out << endl;
				}
			}

			file_out << endl << endl;
		}

		else		// process is not master
			MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
					 MPI_DP, master_id, tag, MPI_COMM_WORLD );


		MPI_Barrier(MPI_COMM_WORLD);
	}

	// 2D arrays
	else if (N[2] ==1) {
		int data_size = 2* local_N3 *  N[1];
		int tag = 123;

		if (my_id == master_id) {
			file_out << "%% Array " <<  N[3] << " x "  << N[1] << endl;

			for (int i=0; i<local_N3; i++) {
				for (int j=0; j<N[1]; j++)
					file_out << real(A(i,0,j)) << " " << imag(A(i,0,j)) << " " ;

				file_out << endl;
			}

			for (int source = 1; source <= numprocs-1; source++) {

				MPI_Recv( reinterpret_cast<DP*>((temp_r).data()), data_size,
						 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				for (int i=0; i<local_N2; i++) {
					for (int j=0; j<N[1]; j++)
						file_out << real(temp_r(i,0,j)) << " " << imag(temp_r(i,0,j)) << " " ;

					file_out << endl;
				}
			}

			file_out << endl << endl;
		}

		else		// process is not master
			MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
					 MPI_DP, master_id, tag, MPI_COMM_WORLD );


		MPI_Barrier(MPI_COMM_WORLD);
	}


}


/**********************************************************************************************

		Master reads data from file_in and distributes to A in all procs

***********************************************************************************************/

void Read_data_MPI
(
	string basis_type,
	ifstream& file_in,
	int N[],
	Array<complx,3> A,
	Array<complx,3> temp_array,
	string format
)
{

	int data_size =  2* local_N1 * N[2] * (N[3]/2 + 1);
	int tag1 = 111;

	int count = 0;

	complex<DP> temp_x;

	if (my_id == master_id) {
		for (int i=0; i<local_N1; i++) {
			for (int j=0; j<N[2]; j++)
				for (int k=0; k<=N[3]/2; k++) {
					if (!file_in.eof())	{
						file_in >> temp_x;
						A(i,j,k) = (complx)temp_x;
					//	cout << i << "  " <<  j << "  "<<  k  << "  "<<  temp_x <<endl;
					}
					else {
						cout << "ERROR: INPUT DATA < REQD SIZE " << endl;
						exit(1);
					}

				}
			// cout << "index i " << i << endl;
		}	

		for (int dest = 1; dest <= numprocs-1; dest++) {

			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<=N[3]/2; k++) {
						if (!file_in.eof()) {
							file_in >> temp_x;
							temp_array(i,j,k) = (complx)temp_x;
						}
						else {
							cout << "ERROR: INPUT DATA < REQD SIZE " << endl;
							exit(1);
						}
					}

			MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						MPI_DP, dest, tag1, MPI_COMM_WORLD );
		}
	}

	else		// process is not master
		MPI_Recv( reinterpret_cast<DP*>(A.data()), data_size,
						MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

	MPI_Barrier(MPI_COMM_WORLD);


/*	if ((my_id == master_id) && (count != (N[1] * N[2] * (N[3]/2 + 1)))) {
		cout << "ERROR:  THE DATA READ DOES NOT MATCH WITH THE GRID SIZE " << endl;
		exit(1);
	} */
}


//***********

void Read_data_kz0plane_MPI
(
	string basis_type,
	ifstream& file_in,
	int N[],
	Array<complx,3> A,
	string format
)
{
	static Array<complx,2>	temp_array(local_N1, N[2]);
	temp_array = 0.0;

	int data_size =  2* local_N1 * N[2];
	int tag1 = 111;

	complex<DP> temp_x;

	if (my_id == master_id) {

		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++)
			{
				file_in >> temp_x;
				A(i,j,0) = (complx)temp_x;
			}

		for (int dest = 1; dest <= numprocs-1; dest++) {

			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
				{
					file_in >> temp_x;
					temp_array(i,j) = (complx)temp_x;
				}

			//	data_size = dest_data_size[dest];
			MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
					 MPI_DP, dest, tag1, MPI_COMM_WORLD );
		}
	}

	else {  // process is not master
		MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
				 MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

		A(Range::all(),Range::all(),0) = temp_array;
	}

	MPI_Barrier(MPI_COMM_WORLD);

}


/**********************************************************************************************

		Master reads data of REDUCED SIZE from file_in and distributes to A in all procs

***********************************************************************************************/

void Read_data_MPI
(
	string basis_type,
	ifstream& file_in,
	int N[],
	int Nreduced[],
	Array<complx,3> A,
	string format
)
{
	static Array<complx,3>	temp_array(local_N1,Nreduced[2],Nreduced[3]/2+1);
	temp_array = 0.0;

	static Array<complx,2>  plane_kxNr1by2(Nreduced[2],Nreduced[3]/2+1);
	complex<DP> temp_x;

	int tag1 = 111;

	int data_size = 2* local_N1 * Nreduced[2] * (Nreduced[3]/2 + 1);

	int last_dest;					// The last processor that will contain reducedA.
	int last_i;
	int	remaining_i;
	int last_left_dest, first_right_dest;
	int first_i;

	// SCFT
	if (basis_type == "SCFT") {

		if (Nreduced[1] > local_N1)
			last_dest = newdivision(Nreduced[1], local_N1) - 1;
		else
			last_dest = 0;


		A = 0.0;		// Initialize for every proc

		if (my_id == master_id) {

			last_i = minimum(local_N1, Nreduced[1]) - 1;
			// lower value of local_N1-1 or Nreducec[1]-1.

			temp_array = 0.0;
			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
					{
						file_in >> temp_x;
						temp_array(i,j,k) = (complx)temp_x;
					}

			A(Range::all(),Range(0,Nreduced[2]/2),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(0,Nreduced[2]/2),Range(0,Nreduced[3]/2));

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1),Range(0,Nreduced[3]/2));


			for (int dest = 1; dest <= last_dest; dest++)
			{
				remaining_i = Nreduced[1] - dest*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				temp_array = 0.0;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
						{
							file_in >> temp_x;
							temp_array(i,j,k) = (complx)temp_x;
						}

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
								MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}
		}

		else if ((my_id > 0) && (my_id <= last_dest))			// process is not master
		{
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
							MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

			A(Range::all(),Range(0,Nreduced[2]/2),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(0,Nreduced[2]/2),Range::all());

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1),Range::all());
		}

		else
			A = 0.0;

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// FOUR BASIS: Note kx=\pm Nr[1]/2 are included both sides in A.
	else if (basis_type == "FOUR") {

		if (local_N1 >= (Nreduced[1]/2+1)) {
			last_left_dest = 0;
			first_right_dest = numprocs-1;
		}
		else if (local_N1 == Nreduced[1]/2) {
			first_right_dest = numprocs-1;
			last_left_dest = newdivision(Nreduced[1]/2+1, local_N1) - 1;
		}
		else if (local_N1 < Nreduced[1]/2)  {
			last_left_dest = newdivision(Nreduced[1]/2+1, local_N1) - 1;
			first_right_dest = numprocs-newdivision(Nreduced[1]/2, local_N1);
		}

		A = 0.0;

		if (my_id == master_id) {
			last_i = minimum(local_N1, Nreduced[1]/2+1) - 1;
			// lower value of local_N1-1 or Nreducec[1]-1.

			temp_array = 0.0;

			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
					{
						file_in >> temp_x;
						temp_array(i,j,k) = (complx)temp_x;
					}

			A(Range::all(),Range(0,Nreduced[2]/2),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(0,Nreduced[2]/2),Range::all());

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1),Range::all());

			for (int dest = 1; dest <= last_left_dest; dest++) {
				remaining_i = Nreduced[1]/2+1 - dest*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				temp_array = 0.0;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
						{
							file_in >> temp_x;
							temp_array(i,j,k) = (complx)temp_x;
						}

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
								MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}


			//****** for dest = first_right_dest

			plane_kxNr1by2 = temp_array(last_i,Range::all(),Range::all());
			// kx=Nr[1]/2 to be saved for the right side of A.

			remaining_i = (Nreduced[1]/2) - (numprocs-1-first_right_dest)*local_N1;
			// remaining i's from the right edge of the proc (i=Nr/2 to Nr-1 or kx=-Nr/2+1:-1

			if (remaining_i >= local_N1) // present proc inside neg k regime
				first_i = 0;
			else // present proc contains part of neg k regme
				first_i = local_N1-remaining_i;

			temp_array = 0.0;
			temp_array(first_i,Range::all(),Range::all()) = plane_kxNr1by2;  // kx=Nr/2

			for (int i=first_i+1; i<local_N1; i++)
				for (int j=0; j<Nreduced[2]; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
					{
						file_in >> temp_x;
						temp_array(i,j,k) = (complx)temp_x;
					}

			if (first_right_dest == 0) {
				// for numproc=1, fill k=-Nr/2 to -1 in master node only
				A(Range(first_i,local_N1),Range(0,Nreduced[2]/2),Range(0,Nreduced[3]/2))
				= temp_array(Range(first_i,local_N1),Range(0,Nreduced[2]/2),Range::all());

				A(Range(first_i,local_N1),Range(N[2]-Nreduced[2]/2,N[2]-1),Range(0,Nreduced[3]/2))
				= temp_array(Range(first_i,local_N1),Range(Nreduced[2]/2,Nreduced[2]-1),Range::all());
			}

			else
				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, first_right_dest, tag1, MPI_COMM_WORLD );


			//****************** Later destinations

			for (int dest = first_right_dest+1; dest <= numprocs-1; dest++) {

				temp_array = 0.0;

				first_i = 0;	// all the procs will have full parts of Nreduced A array.

				for (int i=first_i; i<local_N1; i++)
					for (int j=0; j<Nreduced[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
						{
							file_in >> temp_x;
							temp_array(i,j,k) = (complx)temp_x;
						}

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
							MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}
		}

		// process is not master
		else if ((my_id > 0) && ((my_id <= last_left_dest) || (my_id >= first_right_dest))) {
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
							MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

			A(Range::all(),Range(0,Nreduced[2]/2), Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(0,Nreduced[2]/2),Range::all());

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),Range(0,Nreduced[3]/2))
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1),Range::all());
		}

		else // A(Nreduced) not contained in these procs
			A = 0.0;

		MPI_Barrier(MPI_COMM_WORLD);
	}

}

//***************

void Read_data_kz0plane_MPI
(
 string basis_type,
 ifstream& file_in,
 int N[],
 int Nreduced[],
 Array<complx,3> A,
 string format
 )
{
	static Array<complx,2>	temp_array(local_N1,Nreduced[2]);
	temp_array = 0.0;

	static Array<complx,1>  column_kxNr1by2(Nreduced[2]);
	complex<DP> temp_x;
	int tag1 = 111;

	int data_size = 2* local_N1 * Nreduced[2];

	int last_dest;					// The last processor that will contain reducedA.
	int last_i;
	int	remaining_i;
	int last_left_dest, first_right_dest;
	int first_i;

	// SCFT
	if (basis_type == "SCFT") {

		if (Nreduced[1] > local_N1)
			last_dest = newdivision(Nreduced[1], local_N1) - 1;
		else
			last_dest = 0;


		A = 0.0;		// Initialize for every proc

		if (my_id == master_id) {

			last_i = minimum(local_N1, Nreduced[1]) - 1;
			// lower value of local_N1-1 or Nreducec[1]-1.

			temp_array = 0.0;
			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++)
					file_in >> temp_array(i,j);

			A(Range::all(),Range(0,Nreduced[2]/2),0)
			= temp_array(Range::all(),Range(0,Nreduced[2]/2));

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),0)
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1));

			for (int dest = 1; dest <= last_dest; dest++) {
				remaining_i = Nreduced[1] - dest*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				temp_array = 0.0;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						file_in >> temp_array(i,j);

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}
		}

		else if ((my_id > 0) && (my_id <= last_dest))	{// process is not master
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
					 MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

			A(Range::all(),Range(0,Nreduced[2]/2),0)
			= temp_array(Range::all(),Range(0,Nreduced[2]/2));

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),0)
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1));
		}

		else
			A = 0.0;

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// FOUR BASIS: Note kx=\pm Nr[1]/2 are included both sides in A.
	else if (basis_type == "FOUR") {

		if (local_N1 >= (Nreduced[1]/2+1)) {
			last_left_dest = 0;
			first_right_dest = numprocs-1;
		}
		else if (local_N1 == Nreduced[1]/2) {
			first_right_dest = numprocs-1;
			last_left_dest = newdivision(Nreduced[1]/2+1, local_N1) - 1;
		}
		else if (local_N1 < Nreduced[1]/2)  {
			last_left_dest = newdivision(Nreduced[1]/2+1, local_N1) - 1;
			first_right_dest = numprocs-newdivision(Nreduced[1]/2, local_N1);
		}

		A = 0.0;

		if (my_id == master_id) {
			last_i = minimum(local_N1, Nreduced[1]/2+1) - 1;
			// lower value of local_N1-1 or Nreducec[1]-1.

			temp_array = 0.0;

			for (int i=0; i<=last_i; i++)
				for (int j=0; j<Nreduced[2]; j++)
				{
					file_in >> temp_x;
					temp_array(i,j) = (complx)temp_x;
				}

			A(Range::all(),Range(0,Nreduced[2]/2),0)
			= temp_array(Range::all(),Range(0,Nreduced[2]/2));

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),0)
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1));

			for (int dest = 1; dest <= last_left_dest; dest++) {
				remaining_i = Nreduced[1]/2+1 - dest*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				temp_array = 0.0;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
					{
						file_in >> temp_x;
						temp_array(i,j) = (complx)temp_x;
					}

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}


			//****** for dest = first_right_dest

			column_kxNr1by2 = temp_array(last_i,Range::all());
			// kx=Nr[1]/2 to be saved for the right side of A.

			remaining_i = (Nreduced[1]/2) - (numprocs-1-first_right_dest)*local_N1;
			// remaining i's from the right edge of the proc (i=Nr/2 to Nr-1 or kx=-Nr/2+1:-1

			if (remaining_i >= local_N1) // present proc inside neg k regime
				first_i = 0;
			else // present proc contains part of neg k regime
				first_i = local_N1-remaining_i;

			temp_array = 0.0;
			temp_array(first_i,Range::all()) = column_kxNr1by2;  // kx=Nr/2

			for (int i=first_i+1; i<local_N1; i++)
				for (int j=0; j<Nreduced[2]; j++)
				{
					file_in >> temp_x;
					temp_array(i,j) = (complx)temp_x;
				}

			if (first_right_dest == 0) {
				// for numproc=1, fill k=-Nr/2 to -1 in master node only
				A(Range(first_i,local_N1),Range(0,Nreduced[2]/2),0)
				= temp_array(Range(first_i,local_N1),Range(0,Nreduced[2]/2));

				A(Range(first_i,local_N1),Range(N[2]-Nreduced[2]/2,N[2]-1),0)
				= temp_array(Range(first_i,local_N1),Range(Nreduced[2]/2,Nreduced[2]-1));
			}

			else
				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, first_right_dest, tag1, MPI_COMM_WORLD );


			//****************** Later destinations

			for (int dest = first_right_dest+1; dest <= numprocs-1; dest++) {
				temp_array = 0.0;
				first_i = 0;	// all the procs will have full parts of Nreduced A array.

				for (int i=first_i; i<local_N1; i++)
					for (int j=0; j<Nreduced[2]; j++)
					{
							file_in >> temp_x;
							temp_array(i,j) = (complx)temp_x;
					}

				MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, dest, tag1, MPI_COMM_WORLD );
			}
		}

		// process is not master
		else if ((my_id > 0) && ((my_id <= last_left_dest) || (my_id >= first_right_dest))) {

			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
					 MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);

			A(Range::all(),Range(0,Nreduced[2]/2),0)
			= temp_array(Range::all(),Range(0,Nreduced[2]/2));

			A(Range::all(),Range(N[2]-Nreduced[2]/2,N[2]-1),0)
			= temp_array(Range::all(),Range(Nreduced[2]/2,Nreduced[2]-1));
		}

		else // A(Nreduced) not contained in these procs
			A = 0.0;

		MPI_Barrier(MPI_COMM_WORLD);
	}

}


/**********************************************************************************************

		Master colllects  A  from all procs and writes in file_out

***********************************************************************************************/

/*
void Write_data_MPI_ASCII
(
	string basis_type,
	ofstream& file_out,
	int N[],
	Array<complx,3> A,
	Array<complx,3> temp_array
)
{

//	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
	int tag = 123;

	if (my_id == master_id) {
		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++)
				for (int k=0; k<=N[3]/2; k++)
					file_out << A(i,j,k) << " " ;

		for (int source = 1; source <= numprocs-1; source++) {
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
						MPI_DP, source, tag, MPI_COMM_WORLD, &status );

			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<=N[3]/2; k++)
						file_out << temp_array(i,j,k) << " " ;
		}

		file_out << endl << endl;
	}

	else		// process is not master
		MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
						MPI_DP, master_id, tag, MPI_COMM_WORLD );


	MPI_Barrier(MPI_COMM_WORLD);

}

*/

void Write_data_MPI
(
	string basis_type,
	ofstream& file_out,
	int N[],
	Array<complx,3> A,
	Array<complx,3> temp_array,
	string format
)
{

//	int source_local_N1, source_local_N1_start;			// that of worker
	int data_size = 2* local_N1 * N[2] * (N[3]/2 + 1);
	int tag = 123;

	if (my_id == master_id) {
		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++)
				for (int k=0; k<=N[3]/2; k++)
					file_out << A(i,j,k) << " " ;

		for (int source = 1; source <= numprocs-1; source++) {
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
						MPI_DP, source, tag, MPI_COMM_WORLD, &status );

			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<=N[3]/2; k++)
						file_out << temp_array(i,j,k) << " " ;
		}

		file_out << endl << endl;
	}

	else		// process is not master
		MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
						MPI_DP, master_id, tag, MPI_COMM_WORLD );


	MPI_Barrier(MPI_COMM_WORLD);

}

//**********

void Write_data_kz0plane_MPI
(
	string basis_type,
	ofstream& file_out,
	int N[],
	Array<complx,3> A,
	string format
)
{
	static Array<complx,2>	temp_array(local_N1, N[2]);
	int data_size = 2* local_N1 * N[2];
	int tag = 123;

	if (my_id == master_id) {
		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++)
				file_out << A(i,j,0) << " " ;

		for (int source = 1; source <= numprocs-1; source++) {
			MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
					 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					file_out << temp_array(i,j,0) << " " ;
		}

		file_out << endl << endl;
	}

	else {													// processor: not master
		temp_array = A(Range::all(),Range::all(),0);

		MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
				 MPI_DP, master_id, tag, MPI_COMM_WORLD );
	}

	MPI_Barrier(MPI_COMM_WORLD);

}



/**********************************************************************************************

		Master outputs data of REDUCED SIZE into file_out
		-- We assume that the Nreduced[1] > local_N1, so the Reduced data is to be
				fit in the master node.

***********************************************************************************************/


void Write_data_MPI
(
	string basis_type,
	ofstream& file_out,
	int N[],
	int Nreduced[],
	Array<complx,3> A,
	string format
)
{
	// Problems for 2D implementation...
	if (basis_type == "FOUR") {
		cout << "ERROR: Write field in Nreduced has bugs in FOUR basis.." << endl;
		return;
	}
	
	static Array<complx,3>	temp_array(local_N1,Nreduced[2],Nreduced[3]/2+1);

	int data_size = 2* local_N1 * Nreduced[2] * (Nreduced[3]/2 + 1);  // to be transferred

	int source_local_N1, source_local_N1_start;			// that of worker
	int last_source, last_i;
	int remaining_i;
	int tag = 123;

	if (basis_type == "SCFT") {

		if (Nreduced[1] > local_N1)
			last_source = newdivision(Nreduced[1], local_N1) - 1;
		else
			last_source = 0;

		if (my_id == master_id) {

			last_i = minimum(local_N1, Nreduced[1]) - 1;

			for (int i=0; i<=last_i; i++) {
				for (int j=0; j<=Nreduced[2]/2; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
						file_out << A(i,j,k);

				for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
						file_out << A(i,j,k);
			}


			for (int source = 1; source <= last_source; source++) {
				MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
							MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				remaining_i = Nreduced[1] - source*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
							file_out << temp_array(i,j,k);
			}

			file_out << endl << endl;
		}

		else if ((my_id > 0) && (my_id <= last_source)) {

				temp_array(Range::all(),Range(0,Nreduced[2]/2),Range::all())
				= A(Range::all(),Range(0,Nreduced[2]/2), Range(0,Nreduced[3]/2));

				temp_array(Range::all(),Range(Nreduced[2]/2+1,Nreduced[2]-1),Range::all())
				= A(Range::all(),Range(N[2]-Nreduced[2]/2+1,N[2]-1),Range(0,Nreduced[3]/2));

				MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
								MPI_DP, master_id, tag, MPI_COMM_WORLD );
			}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// FOUR BASIS TYPE
	else if (basis_type == "FOUR") {

		int last_left_source, first_right_source;
		int first_i;

		if (local_N1 >= (Nreduced[1]/2+1)) {
			last_left_source = 0;
			first_right_source = numprocs-1;
		}
		else if ((local_N1 <= Nreduced[1]/2) && (local_N1 >= (Nreduced[1]/2-1)))  {
			first_right_source = numprocs-1;
			last_left_source = newdivision(Nreduced[1]/2+1, local_N1) - 1;
		}
		else if (local_N1 < (Nreduced[1]/2-1))  {
			last_left_source = newdivision(Nreduced[1]/2+1, local_N1) - 1;
			first_right_source = numprocs-newdivision(Nreduced[1]/2-1, local_N1);
		}


		if (my_id == master_id) {
			last_i = minimum(local_N1, Nreduced[1]/2+1) - 1;

			for (int i=0; i<=last_i; i++) {
				for (int j=0; j<=Nreduced[2]/2; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
						file_out << A(i,j,k);

				for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
					for (int k=0; k<=Nreduced[3]/2; k++)
						file_out << A(i,j,k);
			}

			for (int source = 1; source <= last_left_source; source++) {
			/*	MPI_Recv( &source_local_N1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );
				MPI_Recv( &source_local_N1_start, 1, MPI_INT, source, tag, MPI_COMM_WORLD,&status);
				MPI_Recv( &data_size, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status );*/

				MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
							MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				remaining_i = Nreduced[1]/2 - source*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
							file_out << temp_array(i,j,k);
			}

			//**** right sources now..
			if (first_right_source == 0) {// numproc =1, so master node has the full data
				for (int i=N[1]-Nreduced[1]/2+1; i<=N[1]-1; i++) {
					for (int j=0; j<=Nreduced[2]/2; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
							file_out << A(i,j,k);

					for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
						for (int k=0; k<=Nreduced[3]/2; k++)
							file_out << A(i,j,k);
				}
			}

			else
				for (int source = first_right_source; source < numprocs; source++) {
					MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
								MPI_DP, source, tag, MPI_COMM_WORLD, &status );

					remaining_i = (Nreduced[1]/2-1) - (numprocs-1-source)*local_N1;
					// Remaining i's from the right edge of the present proc.
					// Note that negative k are from -Nr/2+1 to -1.. so (Nreduced[1]/2-1)

					if (remaining_i >= local_N1) // present proc inside neg k regime
						first_i = 0;
					else // present proc contains part of neg k regme
						first_i = local_N1-remaining_i;

					for (int i=first_i; i<local_N1; i++)
						for (int j=0; j<Nreduced[2]; j++)
							for (int k=0; k<=Nreduced[3]/2; k++)
								file_out << temp_array(i,j,k);
				}

			file_out << endl << endl;
		}

		// Not a master: compute nodes
		else if ((my_id > 0) && ((my_id <= last_left_source) || (my_id >= first_right_source))) {
				temp_array = 0.0;

				temp_array(Range::all(),Range(0,Nreduced[2]/2),Range::all())
				= A(Range::all(),Range(0,Nreduced[2]/2), Range(0,Nreduced[3]/2));

				temp_array(Range::all(),Range(Nreduced[2]/2+1,Nreduced[2]-1),Range::all())
				= A(Range::all(),Range(N[2]-Nreduced[2]/2+1,N[2]-1),Range(0,Nreduced[3]/2));

				MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
								MPI_DP, master_id, tag, MPI_COMM_WORLD );
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

}

//**********

void Write_data_kz0plane_MPI
(
 string basis_type,
 ofstream& file_out,
 int N[],
 int Nreduced[],
 Array<complx,3> A,
 string format
)
{
	static Array<complx,2>	temp_array(local_N1,Nreduced[2]);

	int data_size = 2* local_N1 * Nreduced[2];  // to be transferred

	int source_local_N1, source_local_N1_start;			// that of worker
	int last_source, last_i;
	int remaining_i;
	int tag = 123;

	if (basis_type == "SCFT") {

		if (Nreduced[1] > local_N1)
			last_source = newdivision(Nreduced[1], local_N1) - 1;
		else
			last_source = 0;

		if (my_id == master_id) {

			last_i = minimum(local_N1, Nreduced[1]) - 1;

			for (int i=0; i<=last_i; i++) {
				for (int j=0; j<=Nreduced[2]/2; j++)
					file_out << A(i,j,0);

				for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
					file_out << A(i,j,0);
			}


			for (int source = 1; source <= last_source; source++) {
				MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				remaining_i = Nreduced[1] - source*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
						file_out << temp_array(i,j);
			}

			file_out << endl << endl;
		}

		else if ((my_id > 0) && (my_id <= last_source)) {
				temp_array(Range::all(),Range(0,Nreduced[2]/2))
				= A(Range::all(),Range(0,Nreduced[2]/2),0);

				temp_array(Range::all(),Range(Nreduced[2]/2+1,Nreduced[2]-1))
				= A(Range::all(),Range(N[2]-Nreduced[2]/2+1,N[2]-1),0);

				MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
						 MPI_DP, master_id, tag, MPI_COMM_WORLD );
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// FOUR BASIS TYPE
	else if (basis_type == "FOUR") {

		int last_left_source, first_right_source;
		int first_i;

		if (local_N1 >= (Nreduced[1]/2+1)) {
			last_left_source = 0;
			first_right_source = numprocs-1;
		}
		else if ((local_N1 <= Nreduced[1]/2) && (local_N1 >= (Nreduced[1]/2-1)))  {
			first_right_source = numprocs-1;
			last_left_source = newdivision(Nreduced[1]/2+1, local_N1) - 1;
		}
		else if (local_N1 < (Nreduced[1]/2-1))  {
			last_left_source = newdivision(Nreduced[1]/2+1, local_N1) - 1;
			first_right_source = numprocs-newdivision(Nreduced[1]/2-1, local_N1);
		}


		if (my_id == master_id) {
			last_i = minimum(local_N1, Nreduced[1]/2+1) - 1;

			for (int i=0; i<=last_i; i++) {
				for (int j=0; j<=Nreduced[2]/2; j++)
					file_out << A(i,j,0);

				for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
					file_out << A(i,j,0);
			}

			for (int source = 1; source <= last_left_source; source++) {
				MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
						 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

				remaining_i = Nreduced[1]/2 - source*local_N1;
				last_i = minimum(local_N1, remaining_i) - 1;

				for (int i=0; i<=last_i; i++)
					for (int j=0; j<Nreduced[2]; j++)
							file_out << temp_array(i,j);
			}

			//**** right sources now..
			if (first_right_source == 0) {// numproc =1, so master node has the full data
				for (int i=N[1]-Nreduced[1]/2+1; i<=N[1]-1; i++) {
					for (int j=0; j<=Nreduced[2]/2; j++)
						file_out << A(i,j,0);

					for (int j=N[2]-Nreduced[2]/2+1; j<N[2]; j++)
						file_out << A(i,j,0);
				}
			}

			else
				for (int source = first_right_source; source < numprocs; source++) {
					MPI_Recv( reinterpret_cast<DP*>((temp_array).data()), data_size,
							 MPI_DP, source, tag, MPI_COMM_WORLD, &status );

					remaining_i = (Nreduced[1]/2-1) - (numprocs-1-source)*local_N1;
					// Remaining i's from the right edge of the present proc.
					// Note that negative k are from -Nr/2+1 to -1.. so (Nreduced[1]/2-1)

					if (remaining_i >= local_N1) // present proc inside neg k regime
						first_i = 0;
					else // present proc contains part of neg k regme
						first_i = local_N1-remaining_i;

					for (int i=first_i; i<local_N1; i++)
						for (int j=0; j<Nreduced[2]; j++)
							file_out << temp_array(i,j);
				}

			file_out << endl << endl;
		}

		// Not a master: compute nodes
		else if ((my_id > 0) && ((my_id <= last_left_source) || (my_id >= first_right_source))) {
			temp_array = 0.0;

			temp_array(Range::all(),Range(0,Nreduced[2]/2))
			= A(Range::all(),Range(0,Nreduced[2]/2),0);

			temp_array(Range::all(),Range(Nreduced[2]/2+1,Nreduced[2]-1))
			= A(Range::all(),Range(N[2]-Nreduced[2]/2+1,N[2]-1),0);

			MPI_Send( reinterpret_cast<DP*>(A.data()), data_size,
					 MPI_DP, master_id, tag, MPI_COMM_WORLD );
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

}



/**********************************************************************************************

		Collect pieced of kz planes from all the procs and put them in plane

***********************************************************************************************/


void Get_XY_plane(int N[], Array<complx,3> A, Array<complx,2> plane, int kz)
{

	int source_local_N1_start;			// that of worker
	int data_size = 2* local_N1 * N[2];
	int tag = 999;

	static Array<complx,2> plane_part(local_N1, N[2]);

	if (my_id == master_id) 	{
		plane( Range(0,local_N1-1), Range::all() ) = A(Range::all(), Range::all(), kz);

		for (int source = 1; source <= numprocs-1; source++) {
			MPI_Recv( reinterpret_cast<DP*>(plane_part.data()), data_size,
						MPI_DP, source, tag, MPI_COMM_WORLD, &status );

			source_local_N1_start = source*local_N1;

			plane( Range(source_local_N1_start, source_local_N1_start+local_N1-1),
							Range::all(), kz )  = plane_part;
		}
	}
	else { // Not a master: compute nodes
		plane_part = A(Range::all(), Range::all(), kz);
		MPI_Send( reinterpret_cast<DP*>(plane_part.data()), data_size,
						MPI_DP, master_id, tag, MPI_COMM_WORLD );
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Now master sends the plane to all the nodes

	tag = 899;
	data_size = 2 * N[1] * N[2];

	if (my_id == master_id)
		for (int dest = 1; dest <= numprocs-1; dest++)
			MPI_Send( reinterpret_cast<DP*>(plane.data()), data_size,
						MPI_DP, dest, tag, MPI_COMM_WORLD );

	else
		MPI_Recv( reinterpret_cast<DP*>(plane.data()), data_size,
						MPI_DP, master_id, tag, MPI_COMM_WORLD, &status );


	MPI_Barrier(MPI_COMM_WORLD);

}

/**********************************************************************************************

		Model energy spectrum for initial condition
		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(1+2.8/12)
		with q = 1.5, b = 0.02

***********************************************************************************************/


void Model_initial_shell_spectrum(int N[],  Array<DP,1> Sk, DP a, DP q, DP b, DP c)
{

//	DP q = 1.5;
//	DP b = 1.0;
//	DP c = 1+2.8/12;
	DP k;	// radius

//	TinyVector<int, 1> len = Sk.length();

	Sk = 0.0;

	for (int i=0; i < (Sk.length())(0); i++)
	{
		k = 1.0*i;
		Sk(i) = a * my_pow(k,4) * exp(- b * pow(k,1.1)) /pow( (my_pow(k,4)+my_pow(q,4)), (DP) c);
    }

}


/**********************************************************************************************
 
Real field data..
 
 ***********************************************************************************************/



void Read_real_field_data_MPI
(
 string basis_type,
 ifstream& file_in,
 int N[],
 Array<complx,3> A,
 Array<complx,3> temp_array,
 string format
 )
{
	
	int data_size =  2* local_N1 * N[2] * (N[3]/2 + 1);
	int tag1 = 111;
	
	int count = 0;
	
	DP temp_x;
	
	if (my_id == master_id) {
		
		// read for master_id
		for (int i=0; i<local_N1; i++)
			for (int j=0; j<N[2]; j++) 
				for (int k=0; k<N[3]/2; k++) {
					if (!file_in.eof())	{
						file_in >> real(A(i,j,k)) >> imag(A(i,j,k));
					}
					else {
						cout << "ERROR: INPUT DATA < REQD SIZE " << endl;
						exit(1);
					}
					
				}
		A(Range::all(), Range::all(), N[3]/2) = 0.0;
		
		// read for other procs..
		for (int dest = 1; dest <= numprocs-1; dest++) {
			
			for (int i=0; i<local_N1; i++)
				for (int j=0; j<N[2]; j++)
					for (int k=0; k<N[3]/2; k++) {
						if (!file_in.eof()) {
							file_in >> real(temp_array(i,j,k)) >> imag(temp_array(i,j,k));
						}
						else {
							cout << "ERROR: INPUT DATA < REQD SIZE " << endl;
							exit(1);
						}
					}
			temp_array(Range::all(), Range::all(), N[3]/2) = 0.0;
			
			MPI_Send( reinterpret_cast<DP*>((temp_array).data()), data_size,
					 MPI_DP, dest, tag1, MPI_COMM_WORLD );
		}
	}
	
	else		// process is not master
		MPI_Recv( reinterpret_cast<DP*>(A.data()), data_size,
				 MPI_DP, master_id, tag1, MPI_COMM_WORLD, &status);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	/*	if ((my_id == master_id) && (count != (N[1] * N[2] * (N[3]/2 + 1)))) {
	 cout << "ERROR:  THE DATA READ DOES NOT MATCH WITH THE GRID SIZE " << endl;
	 exit(1);
	 } */
}

//*********************************   End of array_basic.cc ***********************************
