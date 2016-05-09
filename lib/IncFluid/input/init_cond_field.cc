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

/*! \file  init_cond_field.cc
 *
 * @brief   Read field vars from file field_in_file.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */


#include "../IncFluid.h"


/**********************************************************************************************

						Input from a file: field_in_file

***********************************************************************************************/


void  IncFluid::Init_cond()
{
	int kx, ky, kz;
	complx vz;

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;

	if ( input_field_format == "ASCII" ){
		Input_prefix(field_in_file);

		if (input_vx_vy_switch == 0)
			CV_input(field_in_file, *VF_temp, input_field_format);
		else {
			Read_data_MPI(CV_basis_type, field_in_file, N, *V1, *VF_temp, input_field_format);
			Read_data_MPI(CV_basis_type, field_in_file, N, *V2, *VF_temp, input_field_format);

			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *V3, input_field_format);
		}
	}
	else if ( input_field_format == "HDF5" ) {
		if (input_vx_vy_switch == 0)
			CV_Input_HDF5(N, 3, V1, V2, V3);
		else
			CV_Input_plane_HDF5(N, 3, V1, V2, V3);
	}

	if (input_vx_vy_switch == 1) {

		// Vz(lx,ly,lz=0) already read; construct for lz>=1.
		int kx, ky, kz;
		complx vz;

		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);

					(*V3)(lx,ly,lz) = vz;
				}
	}

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field();

	 if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;


	if (alias_switch == "DEALIAS")		Dealias();

}

//*********************************************************************************************

void  IncFluid::Init_cond(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_scalar(T);

	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Init_cond_RB(T);
}


void  IncFluid::Init_cond_scalar(IncSF& T)
{
	int kx, ky, kz;
	complx vz;

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*T.F) = 0.0;

	if ( input_field_format == "ASCII" ){
		Input_prefix(field_in_file);

		if (input_vx_vy_switch == 0)
			CV_input(field_in_file, *VF_temp, input_field_format);

		else {
			Read_data_MPI(CV_basis_type, field_in_file, N, *V1, *VF_temp, input_field_format);
			
			Read_data_MPI(CV_basis_type, field_in_file, N, *V2, *VF_temp, input_field_format);
			
			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *V3, input_field_format);
		}

		T.CS_input(field_in_file, *VF_temp, input_field_format);
	}


	else if ( input_field_format == "HDF5" ) {
		if (input_vx_vy_switch == 0)
			CV_Input_HDF5(N, 4, V1, V2, V3, T.F);
		else
			CV_Input_plane_HDF5(N, 4, V1, V2, V3, T.F);
	}

	if (input_vx_vy_switch == 1) {
		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);

					(*V3)(lx,ly,lz) = vz;
				}
	}

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(T);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias(T);
}


void  IncFluid::Init_cond_RB(IncSF& T)
{
	if ( input_field_format == "ASCII" ){
		if (globalvar_Pr_switch == "PRZERO")
		{
			Init_cond();

			*T.F = *V1;
			Array_divide_ksqr(basis_type, N, *T.F, kfactor);
		}

		else if (globalvar_Pr_switch == "PRINFTY")
		{
			T.CS_input(field_in_file, *VF_temp, input_field_format);

			if (apply_realitycond_IC_switch == 1)
				Satisfy_reality_array(basis_type, N, *T.F);

			Init_cond_Prinfty(T);
		}

		else
			Init_cond_scalar(T);
	}

	else if ( input_field_format == "HDF5" ) {
		if (globalvar_Pr_switch == "PRZERO")
		{
			Init_cond();

			*T.F = *V1;
			Array_divide_ksqr(basis_type, N, *T.F, kfactor);
		}

		else if (globalvar_Pr_switch == "PRINFTY")
		{
			CV_Input_HDF5(N, 1, T.F);

			if (apply_realitycond_IC_switch == 1)
				Satisfy_reality_array(basis_type, N, *T.F);

			Init_cond_Prinfty(T);
		}

		else
			Init_cond_scalar(T);
	}

	if (basis_type == "SCFT")
		Zero_modes_RB_slip(T);

}


//
//
//*********************************************************************************************
void  IncFluid::Init_cond(IncVF& W)
{
	int kx, ky, kz;
	complx vz;

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;

	if ( input_field_format == "ASCII" ){

		Input_prefix(field_in_file);

		if (input_vx_vy_switch == 0) {
			CV_input(field_in_file, *VF_temp, input_field_format);
			W.CV_input(field_in_file, *VF_temp, input_field_format);
		}

		else {
			Read_data_MPI(CV_basis_type, field_in_file, N, *V1, *VF_temp, input_field_format);
			Read_data_MPI(CV_basis_type, field_in_file, N, *V2, *VF_temp, input_field_format);

			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *V3, input_field_format);

			Read_data_MPI(CV_basis_type, field_in_file, N, *W.V1, *VF_temp, input_field_format);
			Read_data_MPI(CV_basis_type, field_in_file, N, *W.V2, *VF_temp, input_field_format);

			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *W.V3, input_field_format);
		}
	}

	else if ( input_field_format == "HDF5" ) {
		if (input_vx_vy_switch == 0)
			CV_Input_HDF5(N, 6, V1, V2, V3, W.V1, W.V2, W.V3);
		else
			CV_Input_plane_HDF5(N, 6, V1, V2, V3, W.V1, W.V2, W.V3);
	}


	if (input_vx_vy_switch == 1) {
		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);
					(*V3)(lx,ly,lz) = vz;

					Last_component(kx, ky, kz, (*W.V1)(lx,ly,lz), (*W.V2)(lx,ly,lz), vz);
					(*W.V3)(lx,ly,lz) = vz;
				}
	}

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias(W);
}

//
//*********************************************************************************************
void  IncFluid::Init_cond(IncVF& W, IncSF& T)
{
	int kx, ky, kz;
	complx vz;

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;
	(*T.F) = 0.0;

	if ( input_field_format == "ASCII" ){

		Input_prefix(field_in_file);

		if (input_vx_vy_switch == 0) {
			CV_input(field_in_file, *VF_temp, input_field_format);
			W.CV_input(field_in_file, *VF_temp, input_field_format);
		}

		else {
			Read_data_MPI(CV_basis_type, field_in_file, N, *V1, *VF_temp, input_field_format);
			Read_data_MPI(CV_basis_type, field_in_file, N, *V2, *VF_temp, input_field_format);

			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *V3, input_field_format);

			Read_data_MPI(CV_basis_type, field_in_file, N, *W.V1, *VF_temp, input_field_format);
			Read_data_MPI(CV_basis_type, field_in_file, N, *W.V2, *VF_temp, input_field_format);

			Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, *W.V3, input_field_format);
		}

		T.CS_input(field_in_file, *VF_temp, input_field_format);

	}

	else if ( input_field_format == "HDF5" ) {
		if (input_vx_vy_switch == 0)
			CV_Input_HDF5(N, 7, V1, V2, V3, W.V1, W.V2, W.V3, T.F);
		else
			CV_Input_plane_HDF5(N, 7, V1, V2, V3, W.V1, W.V2, W.V3, T.F);
	}

	if (input_vx_vy_switch == 0) {
		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);
					(*V3)(lx,ly,lz) = vz;

					Last_component(kx, ky, kz, (*W.V1)(lx,ly,lz), (*W.V2)(lx,ly,lz), vz);
					(*W.V3)(lx,ly,lz) = vz;
				}
	}


	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W, T);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;


	if (alias_switch == "DEALIAS")		Dealias(W, T);
}


/**********************************************************************************************

						Input from a file: field_in_file(N_in_reduced[])

***********************************************************************************************/


// Fluid
void  IncFluid::Init_cond_reduced()
{

	Input_prefix(field_in_file);

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;

	if (input_vx_vy_switch == 0)
		CV_input(field_in_file, N_in_reduced, input_field_format);

	else {
		int kx, ky, kz;
		complx vz;

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V2, input_field_format);

		*V3 = 0.0;
		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V3, input_field_format);

		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);

					(*V3)(lx,ly,lz) = vz;
				}
	}

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field();

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias();
}

//*********************************************************************************************
// Passive scalar + RB convection


void  IncFluid::Init_cond_reduced(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_reduced_scalar(T);

	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Init_cond_reduced_RB(T);
}


void  IncFluid::Init_cond_reduced_scalar(IncSF& T)
{
	Input_prefix(field_in_file);

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*T.F) = 0.0;

	if (input_vx_vy_switch == 0)
		CV_input(field_in_file, N_in_reduced, input_field_format);

	else {
		int kx, ky, kz;
		complx vz;

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V2, input_field_format);

		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V3, input_field_format);

		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(Ncv[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);

					(*V3)(lx,ly,lz) = vz;
				}
	}

	T.CS_input(field_in_file, N_in_reduced, input_field_format);

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(T);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias(T);
}

void  IncFluid::Init_cond_reduced_RB(IncSF& T)
{

	if (globalvar_Pr_switch == "PRZERO")
	{
		Init_cond_reduced();

		*T.F = *V1;
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);
	}

	if (globalvar_Pr_switch == "PRINFTY")
	{

		T.CS_input(field_in_file, N_in_reduced, input_field_format);

		if (apply_realitycond_IC_switch == 1)
			Satisfy_reality_array(basis_type, N, *T.F);

		Init_cond_Prinfty(T);
	}

	else
		Init_cond_reduced_scalar(T);

	Zero_modes_RB_slip(T);
}

//*********************************************************************************************
// MHD

void  IncFluid::Init_cond_reduced(IncVF& W)
{
	Input_prefix(field_in_file);

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;

	if (input_vx_vy_switch == 0) {
		CV_input(field_in_file, N_in_reduced, input_field_format);
		W.CV_input(field_in_file, N_in_reduced, input_field_format);
	}

	else {
		int kx, ky, kz;
		complx vz;

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V2, input_field_format);

		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V3, input_field_format);

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V2, input_field_format);

		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V3, input_field_format);

		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(Ncv[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);
					(*V3)(lx,ly,lz) = vz;

					Last_component(kx, ky, kz, (*W.V1)(lx,ly,lz), (*W.V2)(lx,ly,lz), vz);
					(*W.V3)(lx,ly,lz) = vz;
				}
	}

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias(W);
}


//*********************************************************************************************
// Magnetoconvection

void  IncFluid::Init_cond_reduced(IncVF& W, IncSF& T)
{
	Input_prefix(field_in_file);

	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;
	(*T.F) = 0.0;

	if (input_vx_vy_switch == 0) {
		CV_input(field_in_file, N_in_reduced, input_field_format);
		W.CV_input(field_in_file, N_in_reduced, input_field_format);
	}

	else {
		int kx, ky, kz;
		complx vz;

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V2, input_field_format);

		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *V3, input_field_format);

		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V1, input_field_format);
		Read_data_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V2, input_field_format);

		Read_data_kz0plane_MPI(CV_basis_type, field_in_file, N, N_in_reduced, *W.V3, input_field_format);

		for (int lx=0; lx<local_N1; lx++)
			for (int ly=0; ly<N[2]; ly++)
				for (int lz=1; lz<=(N[3]/2); lz++) {
					kx = Get_kx(basis_type, lx, N);
					ky = Get_ky3D(basis_type, ly, N);
					kz = lz;

					Last_component(kx, ky, kz, (*V1)(lx,ly,lz), (*V2)(lx,ly,lz), vz);
					(*V3)(lx,ly,lz) = vz;

					Last_component(kx, ky, kz, (*W.V1)(lx,ly,lz), (*W.V2)(lx,ly,lz), vz);
					(*W.V3)(lx,ly,lz) = vz;
				}
	}

	T.CS_input(field_in_file, N_in_reduced, input_field_format);

	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W, T);

	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;

	if (alias_switch == "DEALIAS")		Dealias(W, T);
}


/**********************************************************************************************
 
 Input real-field from a file: field_in_file
 
 ***********************************************************************************************/


void  IncFluid::Init_cond_real_field()
{
	int kx, ky, kz;
	complx vz;
	
	(*V1r) = 0.0;  (*V2r) = 0.0;  (*V3r) = 0.0;
	
	if ( input_field_format == "ASCII" ) {
		Input_prefix(field_in_file);
		RV_input(field_in_file, *VF_temp, input_field_format);
	}
/*	else if ( input_field_format == "HDF5" ) 
			CV_Input_HDF5(N, 3, V1, V2, V3);  */
	
	*V1 = *V1r;
	*V2 = *V2r;
	*V3 = *V3r;
	CV_Forward_transform(*VF_temp_r); 
	// *Vi contains the transformed data now
	
	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field();
	
	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;
	
	//	No need for dealiasing.. Dealiasing is done before product (ui*uj) computation.
	
}



void  IncFluid::Init_cond_real_field(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_real_field_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
			 || (globalvar_prog_kind == "NonBoussinesq"))
		Init_cond_real_field_RB(T);
}


void  IncFluid::Init_cond_real_field_scalar(IncSF& T)
{
	int kx, ky, kz;
	complx vz;
	
	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*T.F) = 0.0;
	
	if ( input_field_format == "ASCII" ){
		Input_prefix(field_in_file);
		RV_input(field_in_file, *VF_temp, input_field_format);
		T.RS_input(field_in_file, *VF_temp, input_field_format);
	}
	
/*	else if ( input_field_format == "HDF5" ) 
		CV_Input_HDF5(N, 4, V1, V2, V3, T.F); */
	
	*V1 = *V1r;
	*V2 = *V2r;
	*V3 = *V3r;
	*T.F = *T.Fr;
	
	CV_Forward_transform(*VF_temp_r);
	T.CS_Forward_transform(*VF_temp_r);
	// *Vi contains the transformed data now
	
	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(T);
	
	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;
	
	//	No need for dealiasing.. Dealiasing is done before product (ui*uj) computation.
}

void  IncFluid::Init_cond_real_field_RB(IncSF& T)
{
	if ( input_field_format == "ASCII" ){
		if (globalvar_Pr_switch == "PRZERO") {
			Init_cond_real_field();
			// *Vi contains the required field now (in Fourier space)
			
			*T.F = *V1;
			Array_divide_ksqr(basis_type, N, *T.F, kfactor);
		}
		
		else if (globalvar_Pr_switch == "PRINFTY") {
			
			T.RS_input(field_in_file, *VF_temp, input_field_format);  // in *T.Fr
			T.CS_Forward_transform(*VF_temp_r);
			// *T.F contains the required field now (in Fourier space)
			
			if (apply_realitycond_IC_switch == 1)
				Satisfy_reality_array(basis_type, N, *T.F);
			
			Init_cond_Prinfty(T);
		}
		
		else{
			Init_cond_real_field_scalar(T);
		}
			
	}
	
/*	else if ( input_field_format == "HDF5" ) {
		if (globalvar_Pr_switch == "PRZERO") {
			Init_cond_real_field();
			RV_Forward_transform(*VF_temp);
			
			*T.F = *V1;
			Array_divide_ksqr(basis_type, N, *T.F, kfactor);
		}
		
		else if (globalvar_Pr_switch == "PRINFTY")	{
			CV_Input_HDF5(N, 1, T.F);
			RS_Forward_transform(*VF_temp);
			
			if (apply_realitycond_IC_switch == 1)
				Satisfy_reality_array(basis_type, N, *T.F);
			
			Init_cond_Prinfty(T);
		}
		
		else
			Init_cond_real_field_scalar(T);
	}  */
	
	if (basis_type == "SCFT")
		Zero_modes_RB_slip(T);
	
}


//
//
//*********************************************************************************************
void  IncFluid::Init_cond_real_field(IncVF& W)
{
	int kx, ky, kz;
	complx vz;
	
	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;
	
	if ( input_field_format == "ASCII" ) {
		Input_prefix(field_in_file);
		RV_input(field_in_file, *VF_temp, input_field_format);
		W.RV_input(field_in_file, *VF_temp, input_field_format);
	}
	
/*	else if ( input_field_format == "HDF5" ) 
			CV_Input_HDF5(N, 6, V1, V2, V3, W.V1, W.V2, W.V3);
	  */
	
	*V1 = *V1r;
	*V2 = *V2r;
	*V3 = *V3r;
	*W.V1 = *W.V1r;
	*W.V2 = *W.V2r;
	*W.V3 = *W.V3r;
	
	CV_Forward_transform(*VF_temp_r);
	W.CV_Forward_transform(*VF_temp_r);
	// *Vi contains the transformed data now
	
	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W);
	
	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;
	
	//	No need for dealiasing.. Dealiasing is done before product (ui*uj) computation.
}

//
//*********************************************************************************************
void  IncFluid::Init_cond_real_field(IncVF& W, IncSF& T)
{
	int kx, ky, kz;
	complx vz;
	
	(*V1) = 0.0;  (*V2) = 0.0;  (*V3) = 0.0;
	(*W.V1) = 0.0;  (*W.V2) = 0.0;  (*W.V3) = 0.0;
	(*T.F) = 0.0;
	
	if ( input_field_format == "ASCII" ) {
		Input_prefix(field_in_file);
		RV_input(field_in_file, *VF_temp, input_field_format);
		W.RV_input(field_in_file, *VF_temp, input_field_format);
		T.RS_input(field_in_file, *VF_temp, input_field_format);
	}
	
/*	else if ( input_field_format == "HDF5" ) 
			CV_Input_HDF5(N, 7, V1, V2, V3, W.V1, W.V2, W.V3, T.F);  */
		
	*V1 = *V1r;
	*V2 = *V2r;
	*V3 = *V3r;
	*W.V1 = *W.V1r;
	*W.V2 = *W.V2r;
	*W.V3 = *W.V3r;
	*T.F = *T.Fr;
	
	CV_Forward_transform(*VF_temp_r);
	W.CV_Forward_transform(*VF_temp_r);
	T.CS_Forward_transform(*VF_temp_r);
	// *Vi contains the transformed data now
	
	if (apply_realitycond_IC_switch == 1)
		Satisfy_reality_condition_field(W, T);
	
	if (my_id == master_id)
		cout  << "Reading of field configurations ended successfully" << endl;
	
	//	No need for dealiasing.. Dealiasing is done before product (ui*uj) computation.
}



//******************************** End of Init_cond_field.cc **********************************

