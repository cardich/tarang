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

//**************************************************************************************************

#include "tarang_hdf5.h"

string DP_to_string(DP time){
	stringstream ss;
	ss.precision(6);
	ss<<time;
	return ss.str();
}

// void Tarang_HDF5_IO::find7z(){
//   system( "env 7z &> /tmp/tarang_7z" );
//   ifstream t7z;
//   t7z.open("/tmp/tarang_7z");
//   char line[100];
//   t7z.getline(line,100);

//   if (string(line) == "env: 7z: No such file or directory"){
//     exist_7z=false;
//     cout<<"7z compression is disabled"<<endl;
//   }
//   else{
//     exist_7z=true;
//     cout<<"7z compression is enabled"<<endl;
//   }
// }

Tarang_HDF5_IO::Tarang_HDF5_IO(string data_dir_name ,int N[])
{
  	/* set the file access template for parallel IO access */
	acc_template = H5Pcreate(H5P_FILE_ACCESS);

	/* create an MPI_INFO object -- on some platforms it is useful to
	 * pass some information onto the underlying MPI_File_open call
	 */
	ierr = MPI_Info_create(&FILE_INFO_TEMPLATE);
	ierr = H5Pset_sieve_buf_size(acc_template, SIEVE_BUF_SIZE);
	ierr = H5Pset_alignment(acc_template, 524288, 262144);

	/*
	 * visit http://www.mpi-forum.org/docs/mpi21-report/node267.htm for information on access_style, ..
	 */
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"access_style", (char*)"write_once");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"collective_buffering", (char*)"true");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_block_size",  (char*)"1048576");
	ierr = MPI_Info_set(FILE_INFO_TEMPLATE, (char*)"cb_buffer_size", (char*)"4194304");

	/* tell the HDF5 library that we want to use MPI-IO to do the writing */
	ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, FILE_INFO_TEMPLATE);

  dim1=0;
  dim2=0;
  dim3=0;


  //input is always from datadir/in/input.hdf5
  input.filename=data_dir_name+"/in/input.h5";


  //output is always to datadir/out/output.hdf5
  output.filename=data_dir_name+"/out/output.h5";
  output.insert_info_at=data_dir_name.length()+11;


  //initialize prefixes
  prefix_prinfinity[0]="T";

  prefix_fluid[0]="CV1";
  prefix_fluid[1]="CV2";
  prefix_fluid[2]="CV3";

  prefix_fluid_T[0]="CV1";
  prefix_fluid_T[1]="CV2";
  prefix_fluid_T[2]="CV3";
  prefix_fluid_T[3]="T";

  prefix_mhd[0]="CV1";
  prefix_mhd[1]="CV2";
  prefix_mhd[2]="CV3";
  prefix_mhd[3]="CW1";
  prefix_mhd[4]="CW2";
  prefix_mhd[5]="CW3";

  prefix_mhd_T[0]="CV1";
  prefix_mhd_T[1]="CV2";
  prefix_mhd_T[2]="CV3";
  prefix_mhd_T[3]="CW1";
  prefix_mhd_T[4]="CW2";
  prefix_mhd_T[5]="CW3";
  prefix_mhd_T[6]="T";


  //assign them to map
  dataset_name[1]=prefix_prinfinity;
  dataset_name[3]=prefix_fluid;
  dataset_name[4]=prefix_fluid_T;
  dataset_name[6]=prefix_mhd;
  dataset_name[7]=prefix_mhd_T;

  data_dir=data_dir_name;

  if (my_id == master_id){
//     find7z();
    output_time_log.open((data_dir+"/out/output_time.d").c_str());

    frequent_time_log.open((data_dir+"/out/output_frequent_time.d").c_str());
    frequent_time_log<<"#Times at which output_frequent were written."<<endl;
  }

  temp_array = new  Array<complx,2>(local_N1, N[2]);

}

Tarang_HDF5_IO::~Tarang_HDF5_IO(){
      /* release the file access template */
  ierr = H5Pclose(acc_template);
  ierr = MPI_Info_free(&FILE_INFO_TEMPLATE);

  if (my_id == master_id)//{
    frequent_time_log.close();
//     system( "rm /tmp/tarang_7z" );
//  }

}


//**************************************************************************************************

void Tarang_HDF5_IO::Read_DP_Array(string dataset_name, int rank, DP *dataArray,...)
{


  va_list vl;
  va_start( vl, dataArray );

  if ( rank>3 )
    return;


  if ( rank >=1 )
    dim1 = va_arg( vl, int );

  if ( rank >=2 )
    dim2 = va_arg( vl, int );

  if ( rank == 3 )
    dim3 = va_arg( vl, int );



    dimens_3d[0] = (dim1)*numprocs;
    dimens_3d[1] = (dim2);
    dimens_3d[2] = (dim3);

    dataspace = H5Screate_simple(rank, dimens_3d, NULL);

    /* figure out the offset into the dataspace for the current processor */
    start_3d[0] = (dim1)*my_id;
    start_3d[1] = 0;
    start_3d[2] = 0;

    stride_3d[0] = 1;
    stride_3d[1] = 1;
    stride_3d[2] = 1;

    count_3d[0] = (dim1);
    count_3d[1] = (dim2);
    count_3d[2] = (dim3);

    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d, stride_3d, count_3d, NULL);


    /* now we need to create a memory space.  In this case, we want to
     *     extract the interior of the patch out when we write, so we'll need
     *     to use a hyperslab */
    dimens_3d[0] = (dim1);
    dimens_3d[1] = (dim2);
    dimens_3d[2] = (dim3);

    memspace = H5Screate_simple(rank, dimens_3d, NULL);

    /* now, use the HDF5 hyperslab function to pick out the portion of the data
     *     array that we intend to store */

    start_3d[0] = 0;
    start_3d[1] = 0;
    start_3d[2] = 0;

    stride_3d[0] = 1;
    stride_3d[1] = 1;
    stride_3d[2] = 1;

    count_3d[0] = (dim1);
    count_3d[1] = (dim2);
    count_3d[2] = (dim3);

    err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			      start_3d, NULL, count_3d, NULL);

  /* now create the dataset -- if it is the first time through, otherwise
   *     open the exisiting dataset */
  dataset = H5Dopen2(file_identifier, dataset_name.c_str(), H5P_DEFAULT);

  /* and finally write the data to disk -- in this call, we need to
   *     include both the memory space and the data space. */
  err = H5Dread(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, dataArray);

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}


//**************************************************************************************************

void Tarang_HDF5_IO::Write_DP_Array(string dataset_name, int rank, DP *dataArray,...)
{


        va_list vl;
        va_start( vl, dataArray );

        if ( rank>3 )
                return;


        if ( rank >=1 )
                dim1 = va_arg( vl, int );

        if ( rank >=2 )
                dim2 = va_arg( vl, int );

        if ( rank == 3 )
                dim3 = va_arg( vl, int );



        dimens_3d[0] = (dim1)*numprocs;
        dimens_3d[1] = (dim2);
        dimens_3d[2] = (dim3);

        dataspace = H5Screate_simple(rank, dimens_3d, NULL);

        /* figure out the offset into the dataspace for the current processor */
        start_3d[0] = (dim1)*my_id;
        start_3d[1] = 0;
        start_3d[2] = 0;

        stride_3d[0] = 1;
        stride_3d[1] = 1;
        stride_3d[2] = 1;

        count_3d[0] = (dim1);
        count_3d[1] = (dim2);
        count_3d[2] = (dim3);

        err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d, stride_3d, count_3d, NULL);


        /* now we need to create a memory space.  In this case, we want to
         *     extract the interior of the patch out when we write, so we'll need
         *     to use a hyperslab */
        dimens_3d[0] = (dim1);
        dimens_3d[1] = (dim2);
        dimens_3d[2] = (dim3);

        memspace = H5Screate_simple(rank, dimens_3d, NULL);

        /* now, use the HDF5 hyperslab function to pick out the portion of the data
         *     array that we intend to store */

        start_3d[0] = 0;
        start_3d[1] = 0;
        start_3d[2] = 0;

        stride_3d[0] = 1;
        stride_3d[1] = 1;
        stride_3d[2] = 1;

        count_3d[0] = (dim1);
        count_3d[1] = (dim2);
        count_3d[2] = (dim3);

        err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                                  start_3d, NULL, count_3d, NULL);

        /* now create the dataset -- if it is the first time through, otherwise
         *     open the exisiting dataset */
        dataset = H5Dcreate2(file_identifier, dataset_name.c_str(), H5_DATATYPE,
                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        /* and finally write the data to disk -- in this call, we need to
         *     include both the memory space and the data space. */
        err = H5Dwrite(dataset, H5_DATATYPE, memspace, dataspace, H5P_DEFAULT, dataArray);

        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);
}


//**************************************************************************************************
void Tarang_HDF5_IO::CV_Input_HDF5(int N[], int number_of_fields, ...)
{
	va_list vl;
	va_start( vl, number_of_fields);

//        if ( my_id == master_id && exist_7z)
// 	  system(("env 7z x "+input.filename+".7z -aoa -o"+data_dir+"/in &> /dev/null").c_str());

	file_identifier = H5Fopen(input.filename.c_str(), H5F_ACC_RDONLY, acc_template);

	for (i=0; i<number_of_fields; i++){
	    A3 = reinterpret_cast<Array<complx, 3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    Read_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], N[3]+2);
	}

	H5Fclose(file_identifier);

// 	 if ( my_id == master_id && exist_7z)
// 	  system(("env 7z t "+input.filename+".7z &> /dev/null && rm "+input.filename+" &> /dev/null").c_str());

}

//**************************************************************************************************
void Tarang_HDF5_IO::CV_Input_plane_HDF5(int N[], int number_of_fields, ...)
{
	va_list vl;
	va_start( vl, number_of_fields);

// 	if ( my_id == master_id && exist_7z)
// 	  system(("env 7z x "+input.filename+".7z -aoa  -o"+data_dir+"/in &> /dev/null").c_str());


	file_identifier = H5Fopen(input.filename.c_str(), H5F_ACC_RDONLY, acc_template);


	for (i=0; i<number_of_fields; i++){
	  if (i==2 || i==5){
	    dataArray = reinterpret_cast<DP *>((*temp_array).data());
	    Read_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], 2);
	    A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    (*A3)(Range::all(),Range::all(),0)=*temp_array;
	  }
	  else
	  {
	    A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    Read_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], N[3]+2);
	  }
	}

	H5Fclose(file_identifier);

// 	if ( my_id == master_id && exist_7z)
// 	  system(("env 7z t "+input.filename+".7z &> /dev/null && rm "+input.filename+" &> /dev/null").c_str());
}


//**************************************************************************************************

void Tarang_HDF5_IO::CV_Output_HDF5(string info, int N[], int number_of_fields, ...)
{

	va_list vl;
	va_start( vl, number_of_fields);

	string filename;

	filename=output.filename;
	if (info.substr(0,8) == "frequent"){
	  frequent=true;
	  filename.insert(output.insert_info_at,"_frequent");
	}
	else{
	  frequent=false;
	  filename.insert(output.insert_info_at,info);
	}

	file_identifier = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);

	for (i=0; i<number_of_fields; i++){
	    A3= reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    Write_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], N[3]+2);
	}



      H5Fclose(file_identifier);

      if ( my_id == master_id )
	if (frequent)
	    frequent_time_log<<info.substr(9)<<endl;
	else
	    output_time_log<<filename<<endl;

// 	if ( my_id == master_id && !frequent && exist_7z){
// 	  system(("env 7za a -t7z "+filename+".7z "+filename+" -aoa &> /dev/null").c_str());
// 	  system(("env 7z t "+filename+".7z &> /dev/null && rm "+filename+" &> /dev/null").c_str());
// 	}

}

//**************************************************************************************************

void Tarang_HDF5_IO::CV_Output_plane_HDF5(string info, int N[], int number_of_fields, ...)
{

	va_list vl;
	va_start( vl, number_of_fields);


	string filename;
	filename=output.filename;

        if (info.substr(0,8) == "frequent"){
                frequent=true;
                filename.insert(output.insert_info_at,"_frequent");
        }
        else{
                frequent=false;
                filename.insert(output.insert_info_at,info);
        }


	file_identifier = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);


	for (i=0; i<number_of_fields; i++){
	  if (i==2 || i==5){
	    A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    *temp_array = (*A3)(Range::all(),Range::all(),0);
	    dataArray = reinterpret_cast<DP *>((*temp_array).data());
	    Write_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], 2);
	  }
	  else
	  {
	    A3 = reinterpret_cast<Array<complx,3>*>(va_arg( vl, int*));
	    dataArray = reinterpret_cast<DP *>((*A3).data());
	    Write_DP_Array(dataset_name[number_of_fields][i], 3, dataArray, local_N1, N[2], N[3]+2);
	  }
	}


	H5Fclose(file_identifier);

      if ( my_id == master_id )
	if (frequent)
	    frequent_time_log<<info.substr(9)<<endl;
	else
	    output_time_log<<filename<<endl;

// 	if ( my_id == master_id && !frequent && exist_7z){
// 	  system(("env 7za a -t7z "+filename+".7z "+filename+" -aoa &> /dev/null").c_str());
// 	  system(("env 7z t "+filename+".7z &> /dev/null && rm "+filename+" &> /dev/null").c_str());
// 	}

}


//*********************************   End of tarang_hdf5.cc ***********************************
