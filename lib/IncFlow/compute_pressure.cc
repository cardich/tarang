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

/*! \file compute_pressure.h 
 * 
 * @brief Functions for computing pressure and pressure spectrum.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date 30 August 2008
 * 
 * @bug  No known bugs
 */ 

//*********************************************************************************************

#include "IncVF.h"


/** @brief Compute pressure and puts it in the CSF of NLIN.
 *
 *	Before entering this function:  
 *		nlin \f$ \vec{N} \f$ contains the nonlinear of udot equation (in Fourier space).
 *		Fluid & Passive scalar: nlin = FT[U.grad U] 
 *		MHD: nlin = T[U.grad U - B.grad B]
 *
 *	Procedure:  \f$ p(K) = \mathcal{F}(\nabla \cdot \vec{N}) \f$, then \f$ p(K) = p(K)/K^2 \f$.
 *
 *  @return  pressure that is stored in F of CSF of NLIN class.   
 */
void IncVF::Compute_pressure()
{
    DP q_omega_t = globalvar_Keplerian_q_shear* globalvar_Keplerian_omega* globalvar_Tnow;
	DP q_omega = globalvar_Keplerian_q_shear* globalvar_Keplerian_omega;
	
    int Kx, Ky, Kz;
    DP kksqr, mu_sqr;
	
	if (globalvar_prog_kind == "KEPLERIAN") {
		Compute_divergence_nlin(); 
        
       // Add i q*omega*t*ky*nlin_x to div(nlin)
        Yderiv_FOUR(N, *nlin1, *VF_temp, kfactor);
        *F = *F + complex<DP>(q_omega_t,0)*(*VF_temp);
		
        // Add -i q*omega*t*ky*u_x to div(nlin)
		Yderiv_FOUR(N, *V1, *VF_temp, kfactor);
        *F = *F - complex<DP>(q_omega,0)*(*VF_temp);

        for (int l1=0; l1<local_N1; l1++) 
            for (int l2=0; l2<N[2]; l2++) 
                for (int l3=0; l3<N[3]/2; l3++) {
                    Kx = Get_kx_FOUR(l1, N)*kfactor[1];
                    Ky = Get_ky3D_FOUR(l2, N)*kfactor[2];
                    Kz = l3*kfactor[3];
                    kksqr = pow2(Kx) + pow2(Ky) + pow2(Kz);
                    
                    mu_sqr = kksqr + pow2(q_omega_t*Ky) + 2*q_omega_t*Kx*Ky;
                    (*F)(l1,l2,l3) /= mu_sqr;
                } 
		
		// To avoid division by zero
		if (my_id == master_id)   
			(*F)(0,0,0) = 0.0;
		
		if (alias_switch == "DEALIAS")   
			Dealias_array(basis_type, N, *F);		// Dealiase pressure
	}
	
	else {
		IncVF::Compute_divergence_nlin();   
		NLIN::CS_divide_ksqr();	
	
		if (alias_switch == "DEALIAS")   
			Dealias_array(basis_type, N, *F);		// Dealiase pressure
	}
}


//*********************************************************************************************

/// Computers pressure shell pectrum \f$ E^p (K) = 1/2 |p(K)|^2 \f$.
/// pressure is CSF inherited by NLIN
void IncVF::Compute_shell_spectrum_pressure()			// Note IncSF function
{
	Compute_shell_spectrum(basis_type, alias_switch, N, *F, 0, *CS_shell_ek, kfactor);			
	// Result in (*CS_iso_ek)
}



/// Computers pressure ring pectrum \f$ E^p (K) = 1/2 |p(K)|^2 \f$.
/// pressure is CSF inherited by NLIN
void IncVF::Compute_ring_spectrum_pressure()				// Note IncSF function
{
	Compute_ring_spectrum(basis_type, alias_switch, N, *F, *CS_sector_angle_array_spectrum, 0, 
								*CS_ring_ek, kfactor);		
}


/// Computers pressure cylindrical ring pectrum \f$ E^p (K) = 1/2 |p(K)|^2 \f$.
/// pressure is CSF inherited by NLIN
void IncVF::Compute_cylinder_ring_spectrum_pressure()				
{

	Compute_cylinder_ring_spectrum(basis_type, alias_switch, N, *F, 
										*CS_cylinder_kpll_array_spectrum, 0, 
										*CS_cylinder_ring_ek, kfactor);
			
}


//*******************************  End of compute_pressure.cc  ********************************

