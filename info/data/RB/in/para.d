#para.d
#Only change the numerical fields -- Don't delete or add any line
# N[i]	
16	1	16


(1)solver-meta-para(int,float,string)		2	6	2
(2)Solver-para-int		1	1
(3)Solver-para-float		1.0	1.0	1.0	1e5	6.8	1 
(4)Solver-para-string		PRLARGE	ULARGE


(1)basis_type		SCFT
(2)alias_switch		DEALIAS
(3)Integration_scheme	RK4
(4)Input_number_mode	ASCII
(5)Output_no_mode	ASCII


(1)hyperdissipation 	0
(2)Low-dim		0
(3)LES			0
(4)Anisotropic-ring	0
(5)Anisotropic-cylinder-ring 	0
(6)structure-fn		1
(7)planar-structure-fn	0
(8)skpq-switch		0 
(9)Output-real-space-r 	0
(10)sincos-horizintal2d   0
(11)free-slip-3d	0
(12)helicity-flux	1
(13)helicity-shell_ET   0
(14)anisotropy-switch(1,2,3)		3
(15)waveno(Actual(0)_or_grid(1))	0
(16)fixed_dt_switch		0
(17)apply_reality_cond_IC_switch	1
(18)apply_reality_cond_alltime_switch	1
(19)output_vx_vy_switch		0
(20)input_vx_vy_switch		0
(21)force_U_switch		1
(22)force_W_switch		0
(23)force_T_switch		1
(24)structure_fn_approx_switch  0
(25)EnergyTr_switch		1
(26)globalvar_fftw_original_switch  0
 

#diss_coeff + diffusion_coefficient (If RB -- ignore it..)
0.1	0.1

#hyper dissipation coefficient + exp+ hyper diffusion coefficient +exp
0	0	0	0



U.Tinit 	0
U.Tfinal 	0.01
U.Tdt		0.01	
U.Tdiag_basic_init	0.0
U.Tdiag_advanced_init	0.0
U.Tstructure_fn_start	0.0


(1)U.Tglobal_save  		0.1
(2)U.Tfield_save		5000.0
(3)U.Tfield_frequent_save 	50.0
(4)U.Tfield_reduced_save	5000.0
(5)U.Trealfield_save  		1000.0
(6)U.Tfield_k_save		100.0
(7)U.Tfield_r_save		100.0
(8)U.Tspectrum_save 		0.1
(9)U.Tpressure_spectrum_save 	2000
(10)U.Tflux_save		0.1
(11)U.Tshell_to_shell_save 	5000.0
(12)U.Tring_to_ring_Ek_save 	5000
(13)U.Tring_to_ring_Ek_save 	5000
(14)U.TCylring_to_ring_Ek_save 	5000
(15)U.TCylring_to_ring_ET_save 	5000
(16)U.Tstruct_fn_save		5000
(17)U.Tmoment			5000
(18)U.Tskpq			5000
(19)U.Tcout_save 		0.01


(1)Ring-spectrum-input-scheme 	0
(2)Number-sectors-ring-spectrum	5
(3)Number-slabs-cyl 		5
(4)qmin 			1
(5)qmax 			2
(6)r_arraysize			6
(7)rmax/str_fn_r_max		2
(8)points_to_skip		1
(9-to-8+D)N_out_reduced		4	4	4


ET:(1)real-imag-switch	0
(2)no-spheres 		8
(3)no-shells 		8	
(4)no-ring-shells 	9
(5)no-sectors 		5
(6)cylinder-no-shells	6
(7)cylinder-no-slabs	5
(8)Shell-input-scheme		0
(9)Ring-shell-input-scheme	0
(10)Ring-ector-input-scheme	0
(11)Cylinder-shell-input-scheme	0
(12)Cylinder-slab-input-scheme	0


(1)no-output-k-modes	7
(2)no-output-r-points	1

1       1       0
2       0       0
2       1       0
1       2	0
3	1	0
1	3	0
4       0       0

1	1	0
	

(1)input-meta-para(int,float,string)		1	2	1
(2)input-para-int:field_input_proc,N_in_reduceced	5	
(3)input-para-float		0.1  0.1
(4)input-para-string		TEST


(1)force-meta-para(int,float,string)		1	2	1
(2)force-para-int:field_input_para,N_in_reduceced	51
(3)force-para-float		5.0	5.0
(4)force-para-string		TEST





