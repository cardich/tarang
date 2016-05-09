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

/*! \file basis_basicfn_inline.h 
 * 
 * @brief basic inline functions and constant definitions that are common to all basis functions (MPI).
 * 
 * @version 4.0 Parallel version
 * @author  M. K. Verma
 * @date    August 2008
 * @bug		No known bug
 */ 



#ifndef _BASIC_BASIC_FN_INLINE_H
#define _BASIC_BASIC_FN_INLINE_H
  
#include <mpi.h>
#include <fftw3-mpi.h>

#include <iostream> 
#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <iomanip>

// external mpi and fftw vars like local_N1, r2c_plan_FOUR, c2r_plan_FOUR etc.
#include "extern_def_vars.h"


// using namespace blitz ;


//*********************************************************************************************

#ifndef _MYMIN
#define _MYMIN

	inline int min(int a, int b) { return a < b ? a : b; }		
	inline DP  min(DP a, DP b)   { return a < b ? a : b; }
	
	inline int minimum(ptrdiff_t a, int b) { return a < b ? a : b; }
	inline int maximum(ptrdiff_t a, int b) { return (b < a) ? a : b; }
	
#endif 

inline int newdivision(ptrdiff_t a, ptrdiff_t b) { return ((a%b) == 0) ? (a/b) : (a/b+1);}

/* inline DP my_pow(DP x, int n) {		// n=0,2,4 only
	if (x<MYEPS) return 0.0;
	else if (n==0) return 1.0;
	else if (n==1) return (x);
	else if (n==2) return (x*x);
	else if (n==3) return (x*x*x);
	else if (n==4) return (my_pow(x,2)*my_pow(x,2));
	else if (n==5) return (my_pow(x,4)*x);
	else if (n==6) return (my_pow(x,2)*my_pow(x,4));
	else if (n==7) return (my_pow(x,6)*x);
	else if (n==8) return (my_pow(x,4)*my_pow(x,4));
	else if (n==9) return (my_pow(x,8)*x);
	else if (n==10) return (my_pow(x,8)*my_pow(x,2));
	else return pow(x,n); // for -Wall. 
} */

inline DP my_pow(DP x, int n) {	
	double temp;
	switch (n) {
		case 0 : return 1.0; break;
		case 1 : return x; break;
		case 2 : return (x*x); break;
		case 3 : return (x*x*x); break;
		case 4 : temp = my_pow(x,2); return (temp*temp); break;	
		case 5 : return (my_pow(x,4)*x); break;
		case 6 : temp = my_pow(x,2); return (temp*temp*temp); break;	
		case 7 : return (my_pow(x,6)*x); break;
		case 8 : return (my_pow(x,4)*my_pow(x,4)); break;
		case 9 : return (my_pow(x,8)*x); break;
		case 10 : return (my_pow(x,8)*my_pow(x,2)); break;
		default : return pow(x,n); 	
	}
}	

//
//
/*! \brief Returns the polar angle, the angle K vector makes with the anisotropic axis
 * 
 * The range of angle is \f$ \Theta = 0:\Theta_{max}\f$.
 *
 * \return \f$ \tan^{-1}(K_{\perp}/K_{||}) \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	

inline DP Get_polar_angle(DP kkperp, DP kkpll)
{
	DP temp;
	
	if (abs(kkperp) > MYEPS)
	{
		if (abs(kkpll) > MYEPS)
			temp = atan(kkperp/kkpll);
		
		else
			temp = M_PI/2;
	}
	
	else		// along pll axis
	{	
		if (kkpll >= 0)
			temp = 0.0;
		
		else
			temp = M_PI;
	}
		
	return (temp >= 0) ? temp : temp + M_PI; 
}


//
//
/*! \brief Returns the azimuthal angle, the angle K_perp vector makes with the h1 axis.
 * 
 * The range of angle is \f$ [0: 2\pi] \f$.
 *
 * \return \f$ \cos^{-1}(K_{||1}/K_{\rho}) \f$.
 */	

inline DP Get_azimuthal_angle(DP kkh1, DP kkh2)
{
	DP angle;
	
	DP kkperp = sqrt(pow2(kkh1) + pow2(kkh2));
	
	if (abs(kkperp) > MYEPS)
	{
		angle = acos(kkh1/kkperp);

		if (kkh2 >= 0)
			return angle;
			
		// when kkh2 < 0, angle is negative.  Convert it positive by adding 2 pi.
		else
			return (2*M_PI - angle);
	}
	
	else // for origin-- to make sure program does not blow up. Avoid it.
		return  0.0;
}


#endif

//***************************** end of basis_basicfn_inline.h *********************************






