# tarang
Spectral numerical methods for the solution of the Navier-Stokes equations and the simulation of turbulent flows.

TARANG is an object-oriented pseudospectral code for the solution of the Navier-Stokes equations and the simulation of turbulent flows, magnetohydrodynamics (MHD), convection, passive scalars, etc. Various basis functions (Chebyshev, Fourier etc) can be used for the computation.

TARANG uses C++ to deliver a reusable, object-oriented library, which can be used as a basis for the development of solvers
for incompressible fluid flows. Solvers containing complex features and boundary conditions, many field variables etc, are
easily constructed with this library. Ready-to-use executables for certain cases are offered.

External libraries:
blitz++ is used for the handling of arrays and mathematical functions in C++.
FFTW is used for FFT (Fast Fourier Transform) calculations.

For more information see documentation.
