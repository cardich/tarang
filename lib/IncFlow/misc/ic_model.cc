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

/*! \file  ic_model.cc
 * 
 * @brief Initial condition where low k modes are forced: TG, ABC flows etc.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug	 No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"

extern Uniform<DP> SPECrand;



//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x) sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_Taylor_Green_field(int k0, DP amp)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;				// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0)  = -I * (amp/8);
			(*V1)(lx_k0, ly_minus_k0, k0) = -I * (amp/8);
			
			(*V2)(lx_k0, ly_k0, k0)  = I * (amp/8);
			(*V2)(lx_k0, ly_minus_k0, k0) = -I * (amp/8);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, ly_k0, k0) = I * (amp/8);
			(*V1)(lx_minus_k0, ly_minus_k0, k0)= I * (amp/8);
			
			(*V2)(lx_minus_k0, ly_k0, k0) = I * (amp/8);
			(*V2)(lx_minus_k0, ly_minus_k0, k0)= -I * (amp/8);
		}
	}
	
	else if (basis_type == "SCFT")
	{
		int lx_k0 = Get_lx("SCFT", k0, N);
		int ly_k0 = Get_ly3D("SCFT", k0, N);
		int ly_minus_k0 = Get_ly3D("SCFT", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0) = I *(amp/8);
			(*V2)(lx_k0, ly_minus_k0, k0) = -I *(amp/8);
		}
	}

}


//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note  Vx = amp*sin(k0 x) sin(k0 y) cos(k0 z)
 *		Vy = amp*cos(k0 x) cos(k0 y) cos(k0 z)
 *		Vz = 0
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_Taylor_Green_field_mag_field(int k0, DP amp)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;				// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0)  = complx(-amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0)  = complx(amp/8, 0.0);
			(*V2)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_minus_k0, ly_minus_k0, k0)= complx(-amp/8, 0.0);
			
			(*V2)(lx_minus_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V2)(lx_minus_k0, ly_minus_k0, k0)= complx(amp/8, 0.0);
		}
	}
	
	else if (basis_type == "SCFT")
	{
		int lx_k0 = Get_lx("SCFT", k0, N);
		int ly_k0 = Get_ly3D("SCFT", k0, N);
		int ly_minus_k0 = Get_ly3D("SCFT", -k0, N);
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V1)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V1)(lx_k0, ly_minus_k0, k0) = complx(-amp/8, 0.0);
			
			(*V2)(lx_k0, ly_k0, k0) = complx(amp/8, 0.0);
			(*V2)(lx_k0, ly_minus_k0, k0) = complx(amp/8, 0.0);
		}
	}

}



//*********************************************************************************************

/** @brief Init condition for Taylor Green flow 
 * 
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @param k0
 *	@param amp  
 *
 */
void IncVF::Setup_ABC_field(int k0, DP amp, DP A, DP B, DP C)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0;			// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		if (my_id == master_id)  
		{
			(*V1)(0, ly_k0, 0)  =	complx(B * amp/2, 0.0);
			(*V1)(0, ly_minus_k0, 0) =  complx(B * amp/2, 0.0);
			(*V1)(0, 0, k0)  =  -I * (C * amp/2);
			
			(*V2)(0, 0, k0)  =   complx(C * amp/2, 0.0);
			
			(*V3)(0, ly_k0, 0)  =   -I * (B * amp/2);
			(*V3)(0, ly_minus_k0, 0) =    I * (B * amp/2);
		}
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) ) 
		{
			(*V2)(lx_k0, 0, 0)  =  -I * (A * amp/2);
			(*V3)(lx_k0, 0, 0)  =   complx(A * amp/2, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{	
			(*V2)(lx_minus_k0, 0, 0) =   I * (A * amp/2);
			(*V3)(lx_minus_k0, 0, 0) =   complx(A * amp/2, 0.0);
		}	
	}
	
	else if ((basis_type == "SCFT") && (my_id == master_id))
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(0);
	}
}



//*********************************************************************************************


void IncVF::Setup_SIX_MODE_field(int k0, DP amp101, DP amp011, DP amp112, DP h)
{

	*V1 = 0.0; 
	*V2 = 0.0;  
	*V3 = 0.0; 
	// initialize

	if (basis_type == "FOUR")
	{
		int lx_k0 = Get_lx("FOUR", k0, N);
		int lx_minus_k0 = Get_lx("FOUR", -k0, N);
		int ly_k0 = Get_ly3D("FOUR", k0, N);
		int ly_minus_k0 = Get_ly3D("FOUR", -k0, N);
		
		DP factor1 = 2.0/sqrt(2+h*h) * amp101;
		DP factor2 = 2.0/sqrt(2+h*h) * amp011;
		DP factor3 = 4.0/sqrt(6+10*h*h) * amp112;
		
		if (my_id == master_id)  
		{
			(*V1)(0, ly_k0, k0) = complx(-h*(factor2/4), 0.0);
			(*V1)(0, ly_minus_k0, k0) = complx(h*(factor2/4), 0.0);
			
			(*V2)(0, ly_k0, k0)  = (I) * (factor2/4);
			(*V2)(0, ly_minus_k0, k0) = (-I) * (factor2/4);
			
			(*V3)(0, ly_k0, k0)  = (-I) * (factor2/4);
			(*V3)(0, ly_minus_k0, k0) = (-I) * (factor2/4); 
		}
		
		if ( (lx_k0 >= 0) && (lx_k0 < local_N1) )
		{
			(*V1)(lx_k0, 0, k0)  = I * (factor1/4);
			(*V2)(lx_k0, 0, k0)  = complx(-h*(factor1/4), 0.0);
			(*V3)(lx_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(*V1)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8);
			(*V1)(lx_k0, ly_minus_k0, 2*k0) = I * (factor3/8);
			
			(*V2)(lx_k0, ly_k0, 2*k0)  = I * (factor3/8) + complx(h*factor3/4, 0.0);
			(*V2)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + complx(h*factor3/4, 0.0);
			
			(*V3)(lx_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) - complx(h*factor3/8, 0.0);
			(*V3)(lx_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) + complx(h*factor3/8, 0.0);
		}
		
		if ( (lx_minus_k0 >= 0) && (lx_minus_k0 < local_N1) ) 
		{
			(*V1)(lx_minus_k0, 0, k0)  = -I * (factor1/4);
			(*V2)(lx_minus_k0, 0, k0)  = complx(h*(factor1/4), 0.0);
			(*V3)(lx_minus_k0, 0, k0)  =  (-I) * (factor1/4);
			
			(*V1)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8);
			(*V1)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8);	
					
			(*V2)(lx_minus_k0, ly_k0, 2*k0)  = I * (factor3/8) - complx(h*factor3/4, 0.0);
			(*V2)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- complx(h*factor3/4, 0.0);
			
			(*V3)(lx_minus_k0, ly_k0, 2*k0)  = (-I) * (factor3/8) + complx(h*factor3/8, 0.0);
			(*V3)(lx_minus_k0, ly_minus_k0, 2*k0) = (-I) * (factor3/8) 
													- complx(h*factor3/8, 0.0);
		}	
		
	}
	
	else if ((basis_type == "SCFT") && (my_id == master_id))
	{
		cout << "Setup_SIX_MODE_field Not allowed in SCFT basis_type " << endl;
		exit(0);
	}

}


/**********************************************************************************************
 
 Based on Pope's book Turbulent flow, page 232.
 
 Model energy spectrum for using Pope's formula
 Sk(k) = C epsilon^2/3 k^-5/3 f_L(kL) f_eta(k eta)
 f_L = [(kL)/sqrt{(kL)^2 + c_L}]^{5/3+p_0}  
 f_eta = exp[ -beta {[(k eta)^4+c_eta^4]^1/4 - c_eta} ]
 c_L = 6.78, c_eta = 0.40, beta = 5.2, p0=2
 C = 1.5; Kolmogorov's constant
 
 Implemented on 9 Oct 2011.
 
 ***********************************************************************************************/


void IncVF::Model_initial_shell_spectrum_Pope(Array<DP,1> Sk, DP epsilon)
{
	
	DP k;	// radius
	DP L;	// size of the box
	DP fL, feta, fL_arg, feta_arg;
	DP p0=2;
	DP cL = 6.78;
	DP ceta = 0.40;
	DP beta = 5.2;
	DP Kolmogrov_const = 1.5;
	
	Sk = 0.0;
	
	if (N[2] > 1) 
		L =  pow(RV_L[1]*RV_L[2]*RV_L[3], (DP) 1./3);  
	else if (N[2] == 1)
		L =  sqrt(RV_L[1]*RV_L[2]); 
		
	DP eta = 0.42*pow(my_pow(dissipation_coefficient,3)/epsilon, (DP) 1./4);
	// Kolmogrov's length.. The factor 0.42 to make Pi(k)-> 0 as k-> infty.
	// Verma's on LM
	
	Sk(0) = 0.0;
	for (int i=1; i < (Sk.length())(0); i++) {
		k = 1.0*i;
		fL_arg = k*L/sqrt(my_pow(k*L,2)+ceta);
		fL = pow(fL_arg, (DP) 5.0/3+p0);
		feta_arg = pow( my_pow(k*eta,4) + my_pow(ceta,4), (DP) 1.0/4) - ceta;
		feta = exp(-beta *feta_arg);
		
		Sk(i) = Kolmogrov_const*pow(epsilon, (DP) 1.0/3)*pow(k, (DP) -5.0/3)*fL*feta;
	}
}

// for scalar field

void IncSF::Model_initial_shell_spectrum_Pope(Array<DP,1> Sk, DP epsilon)
{
	
	DP k;	// radius
	DP L;	// size of the box
	DP fL, feta, fL_arg, feta_arg;
	DP p0=2;
	DP cL = 6.78;
	DP ceta = 0.40;
	DP beta = 5.2;
	DP Kolmogrov_const = 1.5;
	
	Sk = 0.0;
	
	if (NIs[2] > 1) 
		L =  pow(RS_L[1]*RS_L[2]*RS_L[3], (DP) 1./3);  
	else if (NIs[2] == 1)
		L =  sqrt(RS_L[1]*RS_L[3]);  
	
	DP eta = 0.42*pow(my_pow(diffusion_coefficient,3)/epsilon, (DP) 1./4);
	// Kolmogrov's length.. The factor 0.42 to make Pi(k)-> 0 as k-> infty.
	// Verma's on LM
	
	Sk(0) = 0.0;
	for (int i=1; i < (Sk.length())(0); i++) {
		k = 1.0*i;
		fL_arg = k*L/sqrt(my_pow(k*L,2)+ceta);
		fL = pow(fL_arg, (DP) 5.0/3+p0);
		feta_arg = pow( my_pow(k*eta,4) + my_pow(ceta,4), (DP) 1.0/4) - ceta;
		feta = exp(-beta *feta_arg);
		
		Sk(i) = Kolmogrov_const*pow(epsilon, (DP) 1.0/3)*pow(k, (DP) -5.0/3)*fL*feta;
	}
	
}


//*******************************  End of ic_model.cc *****************************************




	
