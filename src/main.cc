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


/*! \file  main.cc
 *
 * @brief  Main file for Tarang (sequential). Calls various modules.
 * @sa main.h
 *
 * @note Reads which turbulence to start (Fluid, MHD, Passive scalcar, RBslip),
 *					the dimensionality of simulation (2D or 3D), and
 *					which directory to use for input and output.
 *
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 */

#include "main.h"

// declarations of global MPI,  FFTW_PLAN_DP and random variables

int my_id;								// My process id
int numprocs;							// No of processors
const int master_id = 0;				// Id of master proc
ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
ptrdiff_t local_N2, local_N2_start;
ptrdiff_t local_N3, local_N3_start;
MPI_Status status;

// for fftw_original
FFTW_PLAN_DP r2c_plan_FOUR, c2r_plan_FOUR;

// for split fftw
FFTW_PLAN_DP r2c_2d_plan_FOUR, c2r_2d_plan_FOUR;
FFTW_PLAN_DP r2c_1d_plan_FOUR, c2r_1d_plan_FOUR;
FFTW_PLAN_DP c2c_1d_forward_plan_FOUR, c2c_1d_inverse_plan_FOUR;

// for SCFT
FFTW_PLAN_DP r2c_plan_SCFT, c2r_plan_SCFT, sintr_plan_SCFT, costr_plan_SCFT,
					isintr_plan_SCFT, icostr_plan_SCFT;  // i for inverse

FFTW_PLAN_DP r2c_1d_plan_SCFT, c2r_1d_plan_SCFT;  // for 2D

int		globalvar_fftw_original_switch;			
// assigned in the solvers,e.g. src/IFLUID/Ifluid_main.cc

int		globalvar_anisotropy_switch;			// 1,2,3 for x,y,z directions
int		globalvar_waveno_switch;				// 0 for actual (default), 1 for grid
int		globalvar_exp_calc_done_switch;
int		globalvar_LES_switch;

string	globalvar_prog_kind;
int		globalvar_T_exists_switch = 1;
int		globalvar_W_exists_switch = 1;

// for RB
string	globalvar_Pr_switch;					// Prandtl number switch (PRLARGNE..) for RB
string	globalvar_RB_Uscaling;					// UBscaling (ULARGE... ) for RB
DP		globalvar_Ra;							// Rayleigh number
DP		globalvar_r;							// normalized Rayleigh number
DP		globalvar_Pr;							// Prandtl number
DP		globalvar_temperature_grad;				// +1 for convection; -1 for stratification;
												// factor for u3 in temperature eqn
DP		globalvar_alpha_DT;						// nondim number alpha*DT

DP      globalvar_Keplerian_omega;
DP      globalvar_Keplerian_q_shear;
DP      globalvar_Tnow;


Uniform<DP> SPECrand;					// Global variable for random no generation


//*********************************************************************************************


int main(int argc, char** argv)
{

  	// MPI INIT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	string data_dir_name;					// data is stored here
	bool display_help=false;


	if (argc == 1){
	  ifstream prog_para_file;				// prog_para_file defined only in the master node
	  prog_para_file.open("prog_para.d");
	  Read_prog_para(prog_para_file, globalvar_prog_kind, data_dir_name);
	}
	else if (argc == 2){
	  if (string(argv[1]) == "--help" || string(argv[1]) == "-h")
	    display_help = true;
	}
	else if (argc == 3)
	{
	  globalvar_prog_kind = argv[1];
	  data_dir_name = argv[2];
	}
	else
	  display_help=true;

        time_t          start;
        time_t          end;

        time (&start);

		// All procs read program parameters

	if (globalvar_prog_kind == "INC_FLUID")
		Ifluid_main(data_dir_name);

	else if (globalvar_prog_kind == "INC_FLUID_DIAG")
		Ifluid_diag_main(data_dir_name);

	else if (globalvar_prog_kind == "INC_SCALAR")
		Iscalar_main(data_dir_name);

	else if (globalvar_prog_kind == "INC_SCALAR_DIAG")
		Iscalar_diag_main(data_dir_name);

	else if (globalvar_prog_kind == "INC_MHD")
		IMHD_main(data_dir_name);

	else if (globalvar_prog_kind == "INC_MHD_DIAG")
		IMHD_diag_main(data_dir_name);

	else if (globalvar_prog_kind == "RB_SLIP")
		RB_slip_main(data_dir_name);

	else if (globalvar_prog_kind == "RB_SLIP_DIAG")
		RB_slip_diag_main(data_dir_name);

	else if (globalvar_prog_kind == "NonBoussinesq")
		NonBoussinesq_main(data_dir_name);
	
	else if (globalvar_prog_kind == "KEPLERIAN")
		Keplerian_main(data_dir_name);

	else
		display_help=true;


        if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG")
		 || (globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG")
		 || (globalvar_prog_kind == "NonBoussinesq"))

		globalvar_T_exists_switch = 1;

	if ((globalvar_prog_kind == "INC_MHD") || (globalvar_prog_kind == "INC_MHD_DIAG") || (globalvar_prog_kind == "KEPLERIAN"))

		globalvar_W_exists_switch = 1;




	// start main program

	SPECrand.seed((unsigned int)time(0));				// Initialize random number seed



	if (display_help == true && my_id == master_id){
	  cout<<my_id<<endl;
	  cout<<endl<<"Zero or two arguments can be given to "<<argv[0]<<endl<<endl;

	  cout<<"Zero arguments:"<<endl;
	  cout<<"\tSolver and data directory are taken from prog_para.d"<<endl<<endl;

	  cout<<"Two arguments:"<<endl;
	  cout<<"\tSolver is the first argument"<<endl;
	  cout<<"\tData directory is second argument"<<endl;
	  cout<<"e.g. mpirun -np 4 "<<argv[0]<<" INC_FLUID `pwd`"<<endl<<endl;

	  cout<<"Solver is one of: INC_FLUID, INC_SCALAR, INC_MHD, RB_SLIP, NonBoussinesq."<<endl<<endl;
	  MPI_Finalize();
	  exit(1);
	}


	time(&end);

	DP vm, rss;
	process_mem_usage(vm, rss);

	DP dif = difftime(end, start);
	if (my_id == master_id){
		cout << endl << "My precision = " << MY_PRECISION << endl;

		cout << endl << "Memory usage by Process " << my_id
		<< "   VM: " << vm << "KB  RSS: " << rss << "KB" << endl;

		cout << endl << "Approximate total memory usage,  VM:"
		<< vm*numprocs << "KB   RSS:" << rss*numprocs << endl;

		cout << endl << "PROGRAM TERMINATING HERE" << endl << endl
			<< "TOTAL TIME ELAPSED: " << dif << " sec" << endl;
	}

	MPI_Finalize();
}



//********************************** End of main.cc *******************************************


