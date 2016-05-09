#para.d
#Only change the numerical fields -- Don't delete or add any line
# N[i]	
16	16	16


(1)solver-meta-para(int,float,string)		1	3	1
(2)Solver-para-int		0
(3)Solver-para-float		1.0	1.0	1.0
(4)Solver-para-string		TEST


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
(6)structure-fn		0
(7)planar-structure-fn	0
(8)skpq-switch		0 
(9)Output-real-space-r 	0
(10)sincos-horizintal2d   0
(11)free-slip-3d	0
(12)helicity-flux	0
(13)helicity-shell_ET   0
(14)anisotropy-switch(1,2,3)		2
(15)waveno(Actual(0)_or_grid(1))	0
(16)fixed_dt_switch		0
(17)apply_reality_cond_IC_switch	0
(18)apply_reality_cond_alltime_switch	0
(19)output_vx_vy_switch		0
(20)input_vx_vy_switch		0
(21)force_U_switch		1
(22)force_W_switch		1
(23)force_T_switch		1
(24)structure_fn_approx_switch  0
(25)EnergyTr_switch		1
(26)globalvar_fftw_original_switch  0

 
#diss_coeff + diffusion_coefficient (If RB -- ignore it..)
1.0

#Hyper_dissipation_switch + hyper dissipation coefficient
0	0


U.Tinit 	0.0
U.Tfinal 	0.001
U.Tdt		0.001	
U.Tdiag_basic_init	0.0
U.Tdiag_advanced_init	100.0
U.Tstructure_fn_start	0.0


(1)U.Tglobal_save  		0.01
(2)U.Tfield_save		0.001
(3)U.Tfield_frequent_save 	0.5
(4)U.Tfield_reduced_save	0.5
(5)U.Trealfield_save  		0.001
(6)U.Tfield_k_save		10
(7)U.Tfield_r_save		10
(8)U.Tspectrum_save 		0.5
(9)U.Tpressure_spectrum_save 	100
(10)U.Tflux_save		0.001
(11)U.Tshell_to_shell_save 	0.5
(12)U.Tring_to_ring_Ek_save 	100
(13)U.Tring_to_ring_Ek_save 	100
(14)U.TCylring_to_ring_Ek_save 	100
(15)U.TCylring_to_ring_ET_save 	100
(16)U.Tstruct_fn_save		0.02
(17)U.Tmoment			100
(18)U.Tskpq			100
(19)U.Tcout_save 		0.0001


(1)Ring-spectrum-input-scheme 	0
(2)Number-sectors-ring-spectrum	5
(3)Number-slabs-cyl 		5
(4)qmin 			1
(5)qmax 			2
(6)r_arraysize			8
(7)rmax/str_fn_r_max		1
(8)neighbour_index_1for nearest 1 
(9-to-8+D)N_out_reduced		4	4	4


ET:(1)real-imag-switch	0
(2)no-spheres 		10
(3)no-shells 		10	
(4)no-ring-shells 	9
(5)no-sectors 		5
(6)cylinder-no-shells	6
(7)cylinder-no-slabs	5
(8)Shell-input-scheme		0
(9)Ring-shell-input-scheme	0
(10)Ring-sector-input-scheme	0
(11)Cylinder-shell-input-scheme	0
(12)Cylinder-slab-input-scheme	0


(1)no-output-k-modes	1
(2)no-output-r-points	1

1	0	1

1	0	1
	

(1)input-meta-para(int,float,string)		1	5	1
(2)input-para-int:field_input_proc,N_in_reduceced	5	
(3)input-para-float		0.1  
(4)input-para-string		TEST


(1)force-meta-para(int,float,string)		1	0	1
(2)force-para-int:field_input_para,N_in_reduceced	0
(3)force-para-float		2 3 0.1 0
(4)force-para-string		TEST





