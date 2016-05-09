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

/*! \file IMHD_main.cc 
 * 
 * @brief Main program executing MHD turbulence
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */


#include "Keplerian_main.h"

// Global variables 
extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2, local_N2_start;
extern MPI_Status status;


extern Uniform<DP> SPECrand;					// for random vars--  declared in main.cc

//********************************************************************************************* 
 
int Keplerian_main(string data_dir_name)
{

	cout << "ENTERING IMHD MAIN my_id " << my_id << endl;

	//**************** Variable declarations *******************
	
	// Variables connected to the field variables and integration
	
	int			N[4];
	int			no_of_fields = 2;				// U, B fields
	DP			diss_coefficients[2];
	DP			hyper_diss_coefficients[2];
	int			hyper_diss_exponents[2];
	
	ifstream	para_file;					// field and input-output parameters stored here
	
	string		string_switches[MAXSIZE_STRING_SWITCHES_ARRAY];
	
	Array<int,1>	switches( MAXSIZE_SWITCHES_ARRAY);
	
	// solver para
	Array<int,1>	solver_meta_para(MAXSIZE_SOLVER_META_PARA);
	Array<int,1>	solver_int_para(MAXSIZE_SOLVER_INT_PARA);
	Array<DP,1>		solver_double_para(MAXSIZE_SOLVER_DOUBLE_PARA);
	string			solver_string_para[MAXSIZE_SOLVER_STRING_PARA];
	
	
	// Variables connected to the output functions
	
	//! Tinit, Tfinal, Tdt, Tdiagnostic_int in index 1..4
	Array<DP, 1>time_para(MAXSIZE_TIME_PARA_ARRAY);					
	time_para = 0.0;

	Array<DP,1>		time_save(MAXSIZE_TIME_SAVE_ARRAY);
	time_save		= 0.0;
	
	Array<int,1>	misc_output_para(MAXSIZE_MISC_OUTPUT_PARA);
	
	//!  contains no_output_k, no_output_r
	Array<int,1>	no_output_k_r(3);							
	Array<int, 2>	out_k_r_array(MAXSIZE_OUT_K_R_ARRAY,4);
	out_k_r_array		= 0;			// intialize
	
	
	
	// Variables connected to the energy transfer
	
	Array<int,1>	ET_parameters(MAXSIZE_ET_PARAMETERS);
	
	Array<DP,1>		ET_shell_radii_sector_array(MAXSIZE_ET_RADII_SECTOR_ARRAY);
	
	ET_shell_radii_sector_array = 0.0;



	// init-cond para
	Array<int,1>	init_cond_meta_para(MAXSIZE_INIT_COND_META_PARA);
	Array<int,1>	init_cond_int_para(MAXSIZE_INIT_COND_INT_PARA);
	Array<DP,1>		init_cond_double_para(MAXSIZE_INIT_COND_DOUBLE_PARA);
	string			init_cond_string_para[MAXSIZE_INIT_COND_STRING_PARA];
	
	// forcing para
	Array<int,1>	force_meta_para(MAXSIZE_FORCE_META_PARA);
	Array<int,1>	force_int_para(MAXSIZE_FORCE_INT_PARA);
	Array<DP,1>		force_double_para(MAXSIZE_FORCE_DOUBLE_PARA);
	string			force_string_para[MAXSIZE_FORCE_STRING_PARA];		
	
	
	//**************** Program starts here ****************
	
	// Read field para, input-output, force field parameters
	
	string		filename = "/in/para.d"; 
	filename = data_dir_name+ filename;  
	para_file.open(filename.c_str());		// filename = data_dir_name/in/para.d
		
	Read_para(para_file, 3, no_of_fields, N, string_switches, switches,
			  diss_coefficients, hyper_diss_coefficients, hyper_diss_exponents,
			  solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
			  time_para, time_save, 
			  misc_output_para, ET_parameters, ET_shell_radii_sector_array,
			  no_output_k_r, out_k_r_array,
			  init_cond_meta_para, init_cond_int_para, 
			  init_cond_double_para, init_cond_string_para,
			  force_meta_para, force_int_para, force_double_para, force_string_para);	

	string_switches[0] = data_dir_name;
	
	globalvar_anisotropy_switch = switches(14);
	globalvar_waveno_switch = switches(15);	
	globalvar_LES_switch = switches(3);	// LES switch 
	globalvar_fftw_original_switch = switches(26);	// 1 if original, 0 otherwise
	
	globalvar_exp_calc_done_switch = 0;		// exp(-nu k^2 dt) not computed yet.

	// Construct output prefix for all the output files
	string prefix_str,  nu_str, kappa_str;							
	ostringstream nu_buffer, kappa_buffer;
	
	nu_buffer << diss_coefficients[0];  
	nu_str = nu_buffer.str(); 
	
	kappa_buffer << diss_coefficients[1];  
	kappa_str = kappa_buffer.str();
	
	prefix_str = "%% nu = " + nu_str + " kappa = " + kappa_str;  
	// goes into all the output files
	
	// kfactor computations	
	DP	kfactor[4];
	
	if (solver_int_para(1) == 0)		// kfactor as given
	{
		kfactor[1] = solver_double_para(1);
		kfactor[2] = solver_double_para(2);
		kfactor[3] = solver_double_para(3);
			
	}
	
	else if (solver_int_para(1) == 1)	// kfactor multiplied by some predefined constant
	{
		kfactor[1] = 2*M_PI * solver_double_para(1);
		kfactor[2] = 2 * solver_double_para(2);
		kfactor[3] = 2*M_PI * solver_double_para(3);
	}
	
	if (my_id == master_id)
		cout << "kfactor: " << kfactor[1] << " "  << kfactor[2] << " "  << kfactor[3] 
			<< endl << endl;
    
    globalvar_Keplerian_omega = solver_double_para(4);
    globalvar_Keplerian_q_shear = solver_double_para(5);
	
    
    if (my_id == master_id) {	
		cout << " Omega =  " << globalvar_Keplerian_omega << endl;
		cout << " q_shear = " <<  globalvar_Keplerian_q_shear << endl;
	}
	
	string basis_type  = string_switches[1];
	
			
	// Local_N1, local_N2, local_N3 assignments 
	
	if (globalvar_fftw_original_switch == 1) {
		if (N[2] > 1) {
			int alloc_local;											  
			alloc_local = FFTW_MPI_LOCAL_SIZE_3D_TRANSPOSED_DP(N[1], N[2], N[3], MPI_COMM_WORLD,
									&local_N1, &local_N1_start, &local_N2, &local_N2_start);
			
			local_N3 = N[3]/2+1;
			local_N3_start = 0;
		} 
		
		else if (N[2] == 1) {
			int alloc_local;											  
			alloc_local = FFTW_MPI_LOCAL_SIZE_2D_TRANSPOSED_DP(N[1], N[3], MPI_COMM_WORLD,
									&local_N1, &local_N1_start, &local_N3, &local_N3_start);
			
			local_N2 = 1;
			local_N2_start = 0;
		}
	}
	
	else {
		if (N[2] > 1) {
			local_N1 = N[1]/numprocs;			
			local_N2 = N[2]/numprocs;
			local_N1_start = my_id * local_N1;
			local_N2_start = my_id * local_N2;
			
			local_N3 = N[3]/2+1;
			local_N3_start = 0;
		} 
		
		else if (N[2] == 1) {
			local_N1 = N[1]/numprocs;
			local_N2 = 1;
			local_N3 = (N[3]/2)/numprocs;
			
			local_N1_start = my_id * local_N1;
			local_N2_start = 0;
			local_N3_start = my_id * local_N3;	
		}	
	}
	
	
	if (my_id == master_id) 
		cout << "No or processors = " << numprocs << endl << endl;
	
	cout << "MY ID, local_N1, local_N1_start, local_N2, local_N2_start, local_N3, " 
		<< "local_N3_start = " << my_id << "  "
		<< local_N1 << "  " << local_N1_start << "  " << local_N2 << "  " << local_N2_start
		<< "  " << local_N3 << "  " << local_N3_start << endl << endl;

	int force_U_switch = switches(21);
	int force_W_switch = switches(22);
	
	// Constructors of Vector fields
	IncFluid  U(N, string_switches, switches, force_U_switch, kfactor, 
				diss_coefficients[0], hyper_diss_coefficients[0], hyper_diss_exponents[0],
				solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
				time_para, time_save,
				misc_output_para, no_output_k_r, out_k_r_array,
				ET_parameters, ET_shell_radii_sector_array,
				init_cond_meta_para, init_cond_int_para, 
				init_cond_double_para, init_cond_string_para,
				force_meta_para, force_int_para, force_double_para, force_string_para
				);		
	
	
	IncVF     W(N, string_switches, switches, force_W_switch, kfactor,
					diss_coefficients[1], hyper_diss_coefficients[1], hyper_diss_exponents[1],
					misc_output_para
				);					
	
	// Create FFTW plans
	
	U.Init_fftw_plan();	
	
	// Initialize vector fields;
	
	if (my_id == master_id)		
		U.Open_field_input_files();	
		
	U.Read_init_cond(W);
	U.Compute_divergence_field();
	
	if (my_id == master_id) {
		real((*W.V1)(0,0,0)) = solver_double_para(6);
		real((*W.V2)(0,0,0)) = solver_double_para(7);
		real((*W.V3)(0,0,0)) = solver_double_para(8);
	}
	
	cout << "mean-mag-field = " << (*W.V1)(0,0,0) << (*W.V2)(0,0,0) << (*W.V3)(0,0,0) << endl;
	
	if (my_id == master_id)	 
		U.Close_field_input_files();

//
//	 

	if (my_id == master_id)	
	{	
		U.Open_output_files(prefix_str);	
		cout << endl << "STARTING THE SIMULATION NOW" << endl;
	}	
	
	
	int  iter=0;  // iterations 
	U.Tnow = U.Tinit;
    globalvar_Tnow = U.Tnow;
	
	// Output at new time..
	U.Output_cout(W);
	
    //	U.Output_all_inloop(W);  // for initial cond

	do 
	{
		U.Compute_force(W);	
					
		U.Compute_nlin(W);						
		//		U.Output_field_k_inloop(W);					// T(k) in the output computation needs nlin
															// Output at the present time.		
		U.Add_force(W);	
		
	
	
		U.Tdt = U.Get_dt(W);
		U.Tnow = U.Tnow + U.Tdt;
		globalvar_Tnow = U.Tnow;
		iter++;
	
		U.Compute_pressure();  
		//	U.Output_pressure_spectrum_inloop();		
		// Output pressure spectrum at the present time
		//	cout << "Pressure " << (*U.F)(1,1,0) << endl;
	
		U.Time_advance(W);
		// FIELD AT new time.
		
		// Output at new time..
		U.CV_Compute_totalenergy();
		W.CV_Compute_totalenergy();
		if ((my_id == 0) && (isnan(U.CV_total_energy) || isnan(W.CV_total_energy)))
		{ 
			cout << "ERROR: Numerical Overflow " << endl;  break; 
		}
		
		U.Compute_divergence_field();
		//	W.Compute_divergence_field();
		
		if ((U.free_slip_verticalwall_switch == 1) && (U.basis_type == "SCFT"))
			U.free_slip_verticalwall_field(W);
		
		if (U.apply_realitycond_alltime_switch == 1)
			U.Satisfy_reality_condition_field(W);
		
		U.Output_all_inloop(W);
			
	}
	while (U.Tnow < U.Tfinal);
	

	if (my_id == master_id)		
		U.Close_output_files();
  
	return(1);
} 

//********************************** End of IMHD_main.cc **************************************	




