# The code system ELSEPA

Original code by: Francesc Salvat, Aleksander Jablonski and Cedric J. Powell (September 27, 2004)

This code system computes scattering amplitudes and cross sections
for ELastic Scattering of Electrons and Positrons by neutral Atoms (Z=1
to 103), molecules and positive ions, and projectiles with kinetic
energies larger than about 10 eV. The calculations are usually performed
by means of relativistic (Dirac) partial wave analysis for an effective
local central-interaction potential. When the convergence of the partial
wave series is too slow, alternative approximate calculation methods are
applied. Scattering by molecules is described by using an independent-
atom approximation in which the scattering amplitudes are obtained by
adding coherently the waves scattered from all the atoms in the
molecule, and performing an average over random orientations of the
molecule.

The distribution package consists of FORTRAN 77 programs (source code
files) and numerical data files. To ensure portability, all files are in
text (`ASCII`) format. The programs conform to ISO Standard FORTRAN 77,
except for the fact that they use `COMPLEX*16` variables and intrinsic
functions (such as `DCMPLX`, `CDABS`, etc.) and `END DO` statements that are
available in all compilers which we have used. 


## FILES IN THE DISTRIBUTION PACKAGE

- `docs/readme.txt` ... the original readme file.
- `src/elsepa.f` ...  calculation of elastic scattering of electrons
and positrons by neutral atoms and positive ions.
Dirac partial wave analysis for real and complex
central potentials; high-energy factorizations.
- `src/elscata.f` ... main program for scattering by atoms and ions.
- `src/elscatm.f`... main program for scattering by molecules.
- `examples` ... examples of input data files for `elscata`.
- `examples/h2o.in` ... example of input data file for `elscatm`.
- `data/z_nnn.den` ... DF electron densities of neutral atoms (Z=1-103).
- `data/z_nnn.dfs` ... high-energy DF screening functions (Z=1-103).
`nnn` (three digits) is the atomic number Z.

The programs `elscata` and `elscatm` read data from input files,
whose structure is described in the heading comments of the FORTRAN 77
source files (`src/*.f`). Examples of input data files (`*.in`) are included
in the distribution package. The calculated differential cross section
for the energy `x.yyyezz` (in eV, E format) is written to an output file
named `dcs_xpyyyezz.dat`, in a format ready for visualization with a
plotting program. Several subroutines write partial information (phase
shifts, imaginary potential, ...) on the standard output unit (=6) as
the calculation progresses; they may also generate files with the
extension '.dat' with self-explanatory contents. This does not impair
the calculation speed and helps the user to get a feeling of the time
that a planned calculation may take and to track partial results of the
calculation. Typical running times on a Pentium 4, 2.8 GHz, are between
a few seconds and about two minutes, depending on the atomic number of
the target, the adopted potential model and the kinetic energy of the
projectile.


## OPERATION INSTRUCTIONS
- Inspect the `Makefile` for possible settings (compiler choice, flags).
- Build `ELSEPA` by running `make`.
- Optionally, install `ELSEPA` by running `make install`
- Make the files in the `data` directory findable by setting the `ELSEPA_DATA`
environment variable. For example, run the H<sub>2</sub>O example by doing:

    ELSEPA_DATA=data ./elscatm < examples/h2o.in


## DOCKER IMAGE

To build the [Docker](http://docker.com) image for ELSEPA, run:

    docker build -t elsepa .

If you want to be sure that the container works, start an interactive session and run the H<sub>2</sub>O example:

    docker run -i -t elsepa
    elscatm < /usr/share/elsepa/examples/h2o.in


## PYTHON SUPPORT

Use ELSEPA from Python via the [PyELSEPA](http://github.com/eScatter/pyelsepa) module.
This module supports parallel execution of ELSEPA through the use of the Docker container,
as well as caching of previous results.


## CITATION

If you use ELSEPA in your research, please cite this paper:
- Salvat, Jablonski and Powell, Computer Physics Communications, Volume 165, Issue 2, 15 January 2005, Pages 157–190, [sciencedirect](http://www.sciencedirect.com/science/article/pii/S0010465504004795)

