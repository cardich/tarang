# Tarang-4.0
# 
# Copyright (C) 2008, 2009  Mahendra K. Verma
#
# Mahendra K. Verma
# Indian Institute of Technology, Kanpur-208016
# UP, India
#
# mkv@iitk.ac.in
#
# This file is part of Tarang-4.0 .
#
# Tarang-4.0 is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# Tarang-4.0 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Tarang-4.0; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
#

#ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
#ADD_DEFINITIONS(-DD3 -DWAVENOACTUAL -DANISDIRN1)
#INCLUDE_DIRECTORIES(/Users/mkv/local/include )
#LINK_DIRECTORIES(/Users/mkv/local/lib )

add_library(fields 
		array_fn/array_basic_inline.h
		array_fn/array_basic.h array_fn/array_basic.cc
		struct_fn.h	struct_fn.cc
		planar_struct_fn.cc
		universal/universal_inline.h
    universal/universal_basic.h universal/universal_basic.cc
		universal/universal_energy.h universal/universal_energy.cc
		universal/universal_ET.h universal/universal_ET.cc
		universal/universal_misc.h universal/universal_misc.cc
		universal/universal_tr.h universal/universal_tr.cc
		tarang_hdf5.h tarang_hdf5.cc
		Cvf.h Cvf.cc
		Csf.h Csf.cc	 
		Rvf.h Rvf.cc
		Rsf.h Rsf.cc
		fields.h) 

#TARGET_LINK_LIBRARIES(fields fftw3_mpi fftw3  mpi m)

