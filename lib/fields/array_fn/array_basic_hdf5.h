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

/*! \file field_basic.h
 * 
 * @brief Some basic Array operations: Array_real_mult, Output_asreal, Model_initial_shell_spectrum
 *
 * @author  M. K. Verma
 * @date Sept 2008
 *
 * @bugs Out of date functions Read_data_split_arrays_MPI, Write_data_split_arrays_MPI contain bugs.
 */ 
 

#ifndef _FIELD_BASIC
#define _FIELD_BASIC

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


#include "array_basic_inline.h"

using namespace blitz;



#endif


//*********************************   End of array_basic.h  ***********************************

