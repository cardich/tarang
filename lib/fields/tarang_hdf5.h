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

/*! \file tarang_hdf5.cc
 * @author  Anando Gopal Chatterjee
 * @version 4.0
 * @date June 2011
 * @bug  No known bugs
 */


#ifndef _TARANG_HDF5
#define _TARANG_HDF5

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


#include "../basis_basicfn/basis_basicfn_inline.h"
#include "../basis_basicfn/basis_basicfn.h"

using namespace blitz;

//---------------------------------HDF5--------------------------------------------

struct h5file
{
	string filename;			//name template e.g.      "out/CV.hdf5"
	int insert_info_at;			//position of insertion
};

string DP_to_string(DP time);

class Tarang_HDF5_IO
{
	hsize_t dimens_3d[3];
	hid_t err;

	hsize_t start_3d[3];
	hsize_t stride_3d[3];
	hsize_t count_3d[3];

	hid_t dataspace, memspace, dataset;
	hid_t acc_template;

	hsize_t      dimsr[3];


	/* define an info object to store MPI-IO information */
	MPI_Info FILE_INFO_TEMPLATE;
	hid_t file_identifier;
	int ierr;

	int dim1, dim2, dim3;
	DP *dataArray;

	stringstream info;			//information to be inserted

	h5file input, output;

	string prefix_prinfinity[1];
	string prefix_fluid[3];
	string prefix_fluid_T[4];
	string prefix_mhd[6];
	string prefix_mhd_T[7];
	map<int,string*> dataset_name;

	int i;

	Array<complx,3> *A3;
	Array<complx,2>	*temp_array;//(local_N1, N[2]);

	bool exist_7z;
	ofstream frequent_time_log;
	ofstream output_time_log;
	bool frequent;
	string data_dir;

	void Read_DP_Array(string dataset_name, int rank, DP *dataArray,...);
	void Read_DP_Array_reduced(string dataset_name, int rank, DP *dataArray,...);
	void Write_DP_Array(string dataset_name, int rank, DP *dataArray, ...);

	void init();
	void find7z();

public:
	Tarang_HDF5_IO(string solver, int N[]);
	~Tarang_HDF5_IO();


	void CV_Input_HDF5(int N[], int number_of_fields, ...);
	void CV_Input_plane_HDF5(int N[], int number_of_fields, ...);
	void CV_Input_reduced_HDF5(int N[], int number_of_fields, ...);

	void CV_Output_HDF5(string info, int N[], int number_of_fields, ...);
	void CV_Output_plane_HDF5(string info, int N[], int number_of_fields, ...);
	void CV_Output_reduced_HDF5(string info, int N[], int number_of_fields, ...);

};


/** @brief Write given complex array to given file
 *  @param file_out_name Name of output file,
 *  @param N[] The size of the array A,
 *  @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$.
 */
// void Write_data_MPI_HDF5(string file_out_name, int N[], Array<complx,3> *A);
// void Read_data_MPI_HDF5 (string file_out_name,          int N[], Array<complx,3> *A);



#endif


//*********************************   End of tarang_hdf5.h  ***********************************

