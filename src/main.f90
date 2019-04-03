
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! main routine for the Elk code
program main
use modmain
use mod_hdf5
use mod_timer
implicit none
! local variables
integer itask
call timer_start(t_runtime,.true.)
#ifdef _MAD_
call madness_init
call mpi_world_initialize(.false.)
#else
call mpi_world_initialize(.true.)
#endif
#ifdef _SIRIUS_
call sirius_platform_initialize(0)
#endif
if (iproc.eq.0) call timestamp(6,"[main] done mpi_world_initialize")
call hdf5_initialize
! read input files
call readinput
if (iproc.eq.0) call timestamp(6,"[main] done readinput")
call papi_initialize(npapievents,papievent)
! perform the appropriate task
do itask=1,ntasks
  task=tasks(itask)
  select case(task)
  case(-1)
    write(*,'("Elk version ",I1.1,".",I1.1,".",I2.2)') version
    stop
  case(0,1,2,3)
    call gndstate
  case(5,6)
    call hartfock
  case(10)
    call dos
  case(14)
    call writesf
  case(15,16)
    call writelsj
  case(20,21)
    call bandstr
  case(25)
    call effmass
  case(31,32,33)
    call rhoplot
  case(41,42,43)
    call potplot
  case(51,52,53)
    call elfplot
  case(61,62,63,162)
    call wfplot
  case(72,73,82,83,142,143,152,153)
    call vecplot
  case(91,92,93)
    call dbxcplot
  case(100,101)
    call fermisurf
  case(102)
    call fermisurfbxsf
  case(110)
    call mossbauer
  case(115)
    call writeefg
  case(120)
    call writepmat
  case(121)
    call dielectric
  case(122)
    call moke
  case(130)
    call writeexpiqr
  case(140)
    call elnes
  case(170)
    call writewfpw
  case(172)
    call dielectricq
  case(175)
    call yambo
  case(190)
    call geomplot
  case(200,201)
    call phonon
  case(210)
    call phdos
  case(220)
    call phdisp
  case(230)
    call writephn
  case(240)
    call epcouple
  case(245)
    call phlwidth
  case(250)
    call alpha2f
  case(300)
    call rdmft
  case(400)
    call writetensmom
  case(500)
    call testcheck
  case(700)
    call sic_gndstate
  case(701)
    call sic_main
  !case(702)
  !  call sic_test_vme
  !case(703)
  !  call sic_test_localize
  case(800)
    call response
  case(801)
    call crpa
  case(802)
    call gwmain
  case(804)
    call genscell
  case(805)
    call genwfdrc
  case(806)
    call writebz
  case(807,808,809)
    call writewann
  case(822, 829)
    call bandrlm
  case(811)
    call dosrlm  
  case(861,862)
    call wann_plot
  case(863)
    call wann_plot_3d
  case(880)
    call test_xc
  case(881)
    call test_bloch_wf
  case(882)
    call writewf
  !case(883)
  !  call test_potcoul
!#ifdef _MAD_
!  case(890)
!    call test_madness
!#endif

#ifdef _SIRIUS_
  case(2000)
    call sirius_dft_ground_state
  case(2001)
    call sirius_crpa
#endif

!  case(2001)
!    call test_sirius_band
!  case(2002)
!    call test_sirius_potential
  case default
    write(*,*)
    write(*,'("Error(main): task not defined : ",I8)') task
    write(*,*)
    stop
  end select
  call mpi_world_barrier
end do
#ifdef _SIRIUS_
call sirius_clear
#endif
call papi_finalize
call hdf5_finalize
call mpi_grid_finalize
call mpi_world_finalize
call timer_stop(t_runtime)
if (iproc.eq.0) write(*,'("Total execution time : ",F12.4," seconds")')timer_get_value(t_runtime)
stop
end program

!BOI
! !TITLE: {\huge{\sc The Elk Code Manual}}\\ \Large{\sc Version 1.0.17}\\ \vskip 0.5cm \includegraphics[height=1cm]{elk_silhouette.pdf}
! !AUTHORS: {\sc J. K. Dewhurst, S. Sharma} \\ {\sc L. Nordstr\"{o}m, F. Cricchio, F. Bultmark} \\ {\sc E. K. U. Gross}
! !AFFILIATION:
! !INTRODUCTION: Introduction
!   Welcome to the {\sf Elk} Code! {\sf Elk} is an all-electron full-potential
!   linearised augmented-plane-wave (FP-LAPW) code for determining the
!   properties of crystalline solids. It was developed originally at the
!   Karl-Franzens-Universit\"{a}t Graz as part of the {\sf EXCITING} EU Research
!   and Training Network project\footnote{{\sf EXCITING} code developed under
!   the Research and Training Network {\sf EXCITING} funded by the EU, contract
!   No. HPRN-CT-2002-00317}. The guiding philosophy during
!   the implementation of the code was to keep it as simple as possible for both
!   users and developers without compromising on its capabilities. All the
!   routines are released under either the GNU General Public License (GPL) or
!   the GNU Lesser General Public License (LGPL) in the hope that they may
!   inspire other scientists to implement new developments in the field of
!   density functional theory and beyond.
!
!   \section{Acknowledgements}
!   Lots of people contributed to the {\sf Elk} code with ideas, checking and
!   testing, writing code or documentation and general encouragement. They
!   include Claudia Ambrosch-Draxl, Clas Persson, Christian Brouder, Rickard
!   Armiento, Andrew Chizmeshya, Per Anderson, Igor Nekrasov, Sushil Auluck,
!   Frank Wagner, Fateh Kalarasse, J\"{u}rgen Spitaler, Stefano Pittalis,
!   Nektarios Lathiotakis, Tobias Burnus, Stephan Sagmeister, Christian
!   Meisenbichler, S\'{e}bastien Leb\`{e}gue, Yigang Zhang, Fritz K\"{o}rmann,
!   Alexey Baranov, Anton Kozhevnikov, Shigeru Suehara, Frank Essenberger,
!   Antonio Sanna, Tyrel McQueen, Tim Baldsiefen, Marty Blaber, Anton Filanovich
!   and Torbj\"{o}rn Bj\"{o}rkman. Special mention of David Singh's very useful
!   book on the LAPW method\footnote{D. J. Singh, {\it Planewaves,
!   Pseudopotentials and the LAPW Method} (Kluwer Academic Publishers, Boston,
!   1994).} must also be made. Finally we would like to acknowledge the generous
!   support of Karl-Franzens-Universit\"{a}t Graz, as well as the EU Marie-Curie
!   Research Training Networks initiative.
!
!   \vspace{24pt}
!   Kay Dewhurst\newline
!   Sangeeta Sharma\newline
!   Lars Nordstr\"{o}m\newline
!   Francesco Cricchio\newline
!   Fredrik Bultmark\newline
!   Hardy Gross
!
!   \vspace{12pt}
!   Halle and Uppsala, April 2010
!   \newpage
!
!   \section{Units}
!   Unless explicitly stated otherwise, {\sf Elk} uses atomic units. In
!   this system $\hbar=1$, the electron mass $m=1$, the Bohr radius $a_0=1$ and
!   the electron charge $e=1$ (note that the electron charge is positive, so
!   that the atomic numbers $Z$ are negative). Thus, the atomic unit of length
!   is 0.52917720859(36) \AA, and the atomic unit of energy is the Hartree which
!   equals 27.21138386(68) eV. The unit of the external magnetic fields is
!   defined such that one unit of magnetic field in {\tt elk.in} equals
!   1715.255397557 Tesla.
!
!   \section{Compiling and running {\sf Elk}}
!   \subsection{Compiling the code}
!   Unpack the code from the archive file. Run the command
!   \begin{verbatim}
!     setup
!   \end{verbatim}
!   in the {\tt elk} directory and select the appropriate system and compiler.
!   We highly recommend that you edit the file {\tt make.inc} and tune the
!   compiler options for your particular system. In particular, make use of
!   machine-optimised {\tt BLAS/LAPACK} libraries if they are available, but
!   make sure they are version $3.x$. Setting the {\sf OpenMP} options of your
!   compiler will enable {\sf Elk} to run in parallel mode on multiprocessor
!   systems. Following this, run
!   \begin{verbatim}
!     make
!   \end{verbatim}
!   This will hopefully compile the entire code and all the libraries into one
!   executable, {\tt elk}, located in the {\tt elk/src} directory. It will also
!   compile two useful auxiliary programs, namely {\tt spacegroup} for producing
!   crystal geometries from spacegroup data and {\tt eos} for fitting equations
!   of state to energy-volume data. If you want to compile everything all over
!   again, then run {\tt make clean} from the {\tt elk} directory, followed by
!   {\tt make}.
!
!   \subsection{Linking with the {\sf Libxc} functional library}
!   {\sf Libxc} is the ETSF library of exchange-correlation functionals. Elk can
!   use the complete set of LDA and GGA functionals available in {\sf Libxc}.
!   In order to do this, first download and compile {\sf Libxc}. This should
!   have produced the file {\tt libxc.a}. Copy this file to the {\tt elk/src}
!   directory and then uncomment the lines indicated in the file
!   {\tt elk/src/libxc.inc}. Once this is done, run {\tt make clean} followed by
!   {\tt make}. To select a particular functional of {\sf Libxc}, use the block
!   \begin{verbatim}
!     xctype
!      100 nx nc
!   \end{verbatim}
!   where {\tt nx} and {\tt nc} are, respectively, the numbers of the exchange
!   and correlation functionals in the {\sf Libxc} library. See the
!   documentation of the library for the functionals and their associated
!   numbers.
!
!   \subsection{Running the code}
!   As a rule, all input files for the code are in lower case and end with the
!   extension {\tt .in}. All output files are uppercase and have the extension
!   {\tt .OUT}. For most cases, the user will only need to modify the file
!   {\tt elk.in}. In this file input parameters are arranged in blocks.
!   Each block consists of a block name on one line and the block variables on
!   subsequent lines. Almost all blocks are optional: the code uses reasonable
!   default values in cases where they are absent. Blocks can appear in any
!   order, if a block is repeated then the second instance is used. Comment
!   lines can be included in the input file and begin with the {\tt !}
!   character.
!
!   The only other input files are those describing the atomic species which go
!   into the crystal. These files are found in the {\tt species} directory and
!   are named with the element symbol and the extension {\tt .in}, for example
!   {\tt Sb.in}. They contain parameters like the atomic charge, mass,
!   muffin-tin radius, occupied atomic states and the type of linearisation
!   required. Users should not have to modify these files in the majority of
!   cases.
!
!   The best way to learn to use {\sf Elk} is to run the examples included
!   with the package. These can be found in the {\tt examples} directory and use
!   many of the code's capabilities. The following section which describes all
!   the input parameters will be of invaluable assistance.
!
!   \section{Input blocks}
!   This section lists all the input blocks available. It is arranged with the
!   name of the block followed by a table which lists each parameter name, what
!   the parameter does, its type and default value. A horizontal line in the
!   table indicates a new line in {\tt elk.in}. Below the table is a brief
!   overview of the block's function.
!
!   \subsection{{\tt atoms}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nspecies} & number of species & integer & 0 \\
!   \hline
!   {\tt spfname(i)} & species filename for species {\tt i} & string & - \\
!   \hline
!   {\tt natoms(i)} & number of atoms for species {\tt i} & integer & - \\
!   \hline
!   {\tt atposl(j,i)} & atomic position in lattice coordinates for atom {\tt j}
!    & real(3) & - \\
!   {\tt bfcmt(j,i)} & muffin-tin external magnetic field in Cartesian
!    coordinates for atom {\tt j} & real(3) & - \\
!   \hline
!   \end{tabularx}\newline\newline
!   Defines the atomic species as well as their positions in the unit cell and
!   the external magnetic field applied throughout the muffin-tin. These fields
!   are used to break spin symmetry and should be considered infinitesimal as
!   they do not contribute directly to the total energy. Collinear calculations
!   are more efficient if the field is applied in the $z$-direction. One could,
!   for example, set up an anti-ferromagnetic crystal by pointing the field on
!   one atom in the positive $z$-direction and in the opposite direction on
!   another atom. If {\tt molecule} is {\tt .true.} then the atomic positions
!   are assumed to be in Cartesian coordinates. See also {\tt sppath},
!   {\tt bfieldc} and {\tt molecule}.
!
!   \subsection{{\tt autokpt}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt autokpt} & {\tt .true.} if the {\bf k}-point set is to be determined
!    automatically & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt radkpt} for details.
!
!   \subsection{{\tt autormt}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt autormt} & {\tt .true.} if muffin-tin radii should be determined
!    automatically & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt rmtapm} for details.
!
!   \subsection{{\tt autoswidth}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt autoswidth} & {\tt .true.} if the smearing parameter {\tt swidth}
!    should be determined automatically & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Calculates the smearing width from the $k$-point density,
!   $V_{\rm BZ}/n_k$; the valence band width, $W$; and an effective mass
!   parameter, $m^{*}$; according to
!   $$ \sigma=\frac{\sqrt{2W}}{m^{*}}\left(\frac{3}{4\pi}
!    \frac{V_{\rm BZ}}{n_k}\right)^{1/3}. $$
!   The variable {\tt mstar} then replaces {\tt swidth} as the control parameter
!   of the smearing width. A large value of $m^{*}$ gives a narrower smearing
!   function. Since {\tt swidth} is adjusted according to the fineness of the
!   $\mathbf{k}$-mesh, the smearing parameter can then be eliminated. It is not
!   recommended that {\tt autoswidth} be used in conjunction with the
!   Fermi-Dirac smearing function, since the electronic temperature will then be
!   a function of the $\mathbf{k}$-point mesh. See T. Bj\"orkman and
!   O. Gr\aa n\"as, {\it Int. J. Quant. Chem.} DOI: 10.1002/qua.22476 (2010) for
!   details. See also {\tt stype} and {\tt swidth}.
!
!   \subsection{{\tt avec}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt avec(1)} & first lattice vector & real(3) & $(1.0,0.0,0.0)$ \\
!   \hline
!   {\tt avec(2)} & second lattice vector & real(3) & $(0.0,1.0,0.0)$ \\
!   \hline
!   {\tt avec(3)} & third lattice vector & real(3) & $(0.0,0.0,1.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Lattice vectors of the crystal in atomic units (Bohr). If {\tt molecule} is
!   {\tt .true.} then these vectors are not used.
!
!   \subsection{{\tt beta0}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt beta0} & adaptive mixing parameter & real & $0.05$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This determines how much of the potential from the previous iteration is
!   mixed with the current potential during the self-consistent loop. It should
!   be made smaller if the calculation is unstable. See {\tt betamax} and also
!   the routine {\tt mixadapt}.
!
!   \subsection{{\tt betamax}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt betamax} & maximum adaptive mixing parameter & real & $0.5$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Maximum allowed mixing parameter used in routine {\tt mixadapt}.
!
!   \subsection{{\tt bfieldc}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt bfieldc} & global external magnetic field in Cartesian coordinates &
!    real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is a constant magnetic field applied throughout the entire unit cell
!   and enters the second-variational Hamiltonian as
!   $$ \frac{g_e}{4c}\,\vec{\sigma}\cdot{\bf B}_{\rm ext}, $$
!   where $g_e$ is the electron $g$-factor. This field is normally used to break
!   spin symmetry for spin-polarised calculations and considered to be
!   infinitesimal with no direct contribution to the total energy. In cases
!   where the magnetic field is finite (for example when computing magnetic
!   response) the external ${\bf B}$-field energy reported in {\tt INFO.OUT}
!   should be added to the total by hand. This field is applied throughout the
!   entire unit cell. To apply magnetic fields in particular muffin-tins use the
!   {\tt bfcmt} vectors in the {\tt atoms} block. Collinear calculations are
!   more efficient if the field is applied in the $z$-direction.
!
!   \subsection{{\tt chgexs}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt chgexs } & excess electronic charge & real & $0.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This controls the amount of charge in the unit cell beyond that required to
!   maintain neutrality. It can be set positive or negative depending on whether
!   electron or hole doping is required.
!
!   \subsection{{\tt deband}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt deband} & initial band energy step size & real & $0.0025$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   The initial step length used when searching for the band energy, which is
!   used as the APW and local-orbital linearisation energies. This is done by
!   first searching upwards in energy until the radial wavefunction at the
!   muffin-tin radius is zero. This is the energy at the top of the band,
!   denoted $E_{\rm t}$. A downward search is now performed from $E_{\rm t}$
!   until the slope of the radial wavefunction at the muffin-tin radius is zero.
!   This energy, $E_{\rm b}$, is at the bottom of the band. The band energy is
!   taken as $(E_{\rm t}+E_{\rm b})/2$. If either $E_{\rm t}$ or $E_{\rm b}$
!   cannot be found then the band energy is set to the default value.
!
!   \subsection{{\tt deltaem}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt deltaem} & the size of the ${\bf k}$-vector displacement used when
!    calculating numerical derivatives for the effective mass tensor & real &
!    $0.025$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt ndspem} and {\tt vklem}.
!
!   \subsection{{\tt deltaph}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt deltaph} & the size of the atomic displacement used for calculating
!    dynamical matrices & real & $0.03$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Phonon calculations are performed by constructing a supercell corresponding
!   to a particular ${\bf q}$-vector and making a small periodic displacement of
!   the atoms. The magnitude of this displacement is given by {\tt deltaph}.
!   This should not be made too large, as anharmonic terms could then become
!   significant, neither should it be too small as this can introduce numerical
!   error.
!
!   \subsection{{\tt dos}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nwdos} & number of frequency/energy points in the DOS or optics plot &
!    integer & $500$ \\
!   {\tt ngrdos} & effective {\bf k}-point mesh size to be used for Brillouin
!    zone integration & integer & $100$ \\
!   {\tt nsmdos} & level of smoothing applied to DOS/optics output & integer &
!    $0$ \\
!   \hline
!   {\tt wintdos} & frequency/energy window for the DOS or optics plot &
!    real(2) & $(-0.5,0.5)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   DOS and optics plots require integrals of the kind
!   $$ g(\omega_i)=\frac{\Omega}{(2\pi)^3}\int_{\rm BZ} f({\bf k})
!    \delta(\omega_i-e({\bf k}))d{\bf k}. $$
!   These are calculated by first interpolating the functions $e({\bf k})$ and
!   $f({\bf k})$ with the trilinear method on a much finer mesh whose size is
!   determined by {\tt ngrdos}. Then the $\omega$-dependent histogram of the
!   integrand is accumulated over the fine mesh. If the output function is noisy
!   then either {\tt ngrdos} should be increased or {\tt nwdos} decreased.
!   Alternatively, the output function can be artificially smoothed up to a
!   level given by {\tt nsmdos}. This is the number of successive 3-point
!   averages to be applied to the function $g$.
!
!   \subsection{{\tt dosmsum}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt dosmsum} & {\tt .true.} if the partial DOS is to be summed over $m$ &
!    logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   By default, the partial density of states is resolved over $(l,m)$ quantum
!   numbers. If {\tt dosmsum} is set to {\tt .true.} then the partial DOS is
!   summed over $m$, and thus depends only on $l$.
!
!   \subsection{{\tt epsband}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epsband} & convergence tolerance for determining band energies & real &
!    $1\times 10^{-6}$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   APW and local-orbital linearisation energies are determined from the band
!   energies. Good convergence of these energies, especially for low-lying,
!   states is required for accurate total energies. See also {\tt deband}.
!
!   \subsection{{\tt epschg}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epschg} & maximum allowed error in the calculated total charge beyond
!    which a warning message will be issued & real & $1\times 10^{-3}$ \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt epsforce}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epsforce} & convergence tolerance for the forces during a structural
!    optimisation run & real & $5\times 10^{-4}$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   If the mean absolute value of the atomic forces is less than {\tt epsforce}
!   then the structural optimisation run is ended. See {\tt tasks}.
!
!   \subsection{{\tt epslat}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epslat } & vectors with lengths less than this are considered zero &
!    real & $10^{-6}$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Sets the tolerance for determining if a vector or its components are zero.
!   This is to account for any numerical error in real or reciprocal space
!   vectors.
!
!   \subsection{{\tt epsocc}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epsocc} & smallest occupancy for which a state will contribute to the
!    density & real & $1\times 10^{-8}$ \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt epspot}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt epspot} & convergence criterion for the effective potential and field &
!    real & $1\times 10^{-6}$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   If the RMS change in the effective potential and magnetic field is smaller
!   than {\tt epspot}, then the self-consistent loop is considered converged
!   and exited. For structural optimisation runs this results in the forces
!   being calculated, the atomic positions updated and the loop restarted. See
!   also {\tt maxscl}.
!
!   \subsection{{\tt emaxelnes}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt emaxelnes} & maximum allowed initial-state eigenvalue for ELNES
!    calculations & real & $-1.2$ \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt fixspin}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt fixspin} & 0 for no fixed spin moment (FSM), 1 for total FSM, 2 for
!    local muffin-tin FSM, and 3 for both total and local FSM & integer & 0 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Set to 1, 2 or 3 for fixed spin moment calculations. To fix only the
!   direction and not the magnitude set to $-1$, $-2$ or $-3$. See also
!   {\tt momfix}, {\tt mommtfix}, {\tt taufsm} and {\tt spinpol}.
!
!   \subsection{{\tt fracinr}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt fracinr} & fraction of the muffin-tin radius up to which {\tt lmaxinr}
!    is used as the angular momentum cut-off & real & $0.25$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt lmaxinr}.
!
!   \subsection{{\tt gmaxvr}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt gmaxvr} & maximum length of $|{\bf G}|$ for expanding the interstitial
!    density and potential & real & $12.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This variable has a lower bound which is enforced by the code as follows:
!   $$ {\rm gmaxvr}\rightarrow\max\,({\rm gmaxvr},2\times{\rm gkmax}
!    +{\rm epslat}) $$
!   See {\tt rgkmax}.
!
!   \subsection{{\tt intraband}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt intraband} & {\tt .true.} if the intraband (Drude-like) contribution is
!    to be added to the dieletric tensor & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt isgkmax}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt isgkmax} & species for which the muffin-tin radius will be used for
!    calculating {\tt gkmax} & integer & $-1$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   By default ($-1$) the average muffin-tin radius is used for determining
!   {\tt gkmax} from {\tt rgkmax}. This can be changed by setting {\tt isgkmax}
!   either to the desired species number, or to $-2$ in which case the smallest
!   radius is used.
!
!   \subsection{{\tt kstlist}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt kstlist(i)} & $i$th ${\bf k}$-point and state pair & integer(2) &
!    $(1,1)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is a user-defined list of ${\bf k}$-point and state index pairs which
!   are those used for plotting wavefunctions and writing ${\bf L}$, ${\bf S}$
!   and ${\bf J}$ expectation values. Only the first pair is used by the
!   aforementioned tasks. The list should be terminated by a blank line.
!
!   \subsection{{\tt lda+u}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt ldapu} & type of LDA+$U$ calculation & integer & 0 \\
!   {\tt inptypelu} & type of input for LDA+U calculation & integer & 1 \\
!   \hline
!   {\tt is} & species number & integer & - \\
!   {\tt l} & angular momentum value & integer & -1 \\
!   {\tt u} & the desired $U$ value & real & $0.0$ \\
!   {\tt j} & the desired $J$ value & real & $0.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This block contains the parameters required for an LDA+$U$ calculation, with
!   the list of parameters for each species terminated with a blank line. The
!   type of double counting required is set with the parameter {\tt ldapu}.
!   Currently implemented are:\newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   0 & No LDA+$U$ calculation \\
!   1 & Fully localised limit (FLL) \\
!   2 & Around mean field (AFM) \\
!   3 & An interpolation between FLL and AFM \\
!   \end{tabularx}\newline\newline
!   The type of input parameters is set with the parameter {\tt inptypelu}.
!   The current possibilities are:\newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   1 & U and J \\
!   2 & Slater parameters \\
!   3 & Racah parameters \\
!   4 & Yukawa screening length \\
!   5 & U and determination of corresponding Yukawa screening length
!   \end{tabularx}\newline\newline
!   See (amongst others) {\it Phys. Rev. B} {\bf 67}, 153106 (2003),
!   {\it Phys. Rev. B} {\bf 52}, R5467 (1995), {\it Phys. Rev. B} {\bf 60},
!   10673 (1999), and {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
!   \subsection{{\tt lmaxapw}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lmaxapw} & angular momentum cut-off for the APW functions & integer &
!    $8$ \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt lmaxinr}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lmaxinr} & angular momentum cut-off for the muffin-tin density and
!    potential on the inner part of the muffin-tin & integer & 2 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Close to the nucleus, the density and potential is almost spherical and
!   therefore the spherical harmonic expansion can be truncated a low angular
!   momentum. See also {\tt fracinr}.
!
!   \subsection{{\tt lmaxmat}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lmaxmat} & angular momentum cut-off for the outer-most loop in the
!    hamiltonian and overlap matrix setup & integer & 5 \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt lmaxvr}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lmaxvr} & angular momentum cut-off for the muffin-tin density and
!    potential & integer & 7 \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt lmirep}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lmirep} & {\tt .true.} if the $Y_{lm}$ basis is to be transformed
!    into the basis of irreducible representations of the site symmetries for
!    DOS plotting & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   When lmirep is set to .true., the spherical harmonic basis is transformed
!   into one in which the site symmetries are block diagonal. Band characters
!   determined from the density matrix expressed in this basis correspond to
!   irreducible representations, and allow the partial DOS to be resolved into
!   physically relevant contributions, for example $e_g$ and $t_{2g}$.
!
!   \subsection{{\tt lradstp}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt lradstp} & radial step length for determining coarse radial mesh &
!    integer & 4 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Some muffin-tin functions (such as the density) are calculated on a coarse
!   radial mesh and then interpolated onto a fine mesh. This is done for the
!   sake of efficiency. {\tt lradstp} defines the step size in going from the
!   fine to the coarse radial mesh. If it is too large, loss of precision may
!   occur.
!
!   \subsection{{\tt maxitoep}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt maxitoep} & maximum number of iterations when solving the exact
!   exchange integral equations & integer & 120 \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt tau0oep}.
!
!   \subsection{{\tt maxscl}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt maxscl } & maximum number of self-consistent loops allowed & integer &
!    200 \\
!   \hline
!   \end{tabularx}\newline\newline
!   This determines after how many loops the self-consistent cycle will
!   terminate if the convergence criterion is not met. If {\tt maxscl} is $1$
!   then the density and potential file, {\tt STATE.OUT}, will {\bf not} be
!   written to disk at the end of the loop. See {\tt epspot}.
!
!   \subsection{{\tt mixtype}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt mixtype } & type of mixing required for the potential & integer & 1 \\
!   \hline
!   \end{tabularx}\newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   1 & Adaptive linear mixing \\
!   2 & Pulay mixing, {\it Chem. Phys. Lett.} {\bf 73}, 393 (1980) \\
!   \end{tabularx}
!
!   \subsection{{\tt molecule}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt molecule} & {\tt .true.} if the system is an isolated molecule &
!    logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   If {\tt molecule} is {\tt .true.}, then the atomic positions, ${\bf a}$,
!   given in the {\tt atoms} block are assumed to be in Cartesian coordinates.
!
!   \subsection{{\tt momfix}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt momfix} & the desired total moment for a FSM calculation &
!    real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Note that all three components must be specified (even for collinear
!   calculations). See {\tt fixspin}, {\tt taufsm} and {\tt spinpol}.
!
!   \subsection{{\tt mommtfix}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt is} & species number & integer & 0 \\
!   {\tt ia} & atom number & integer & 0 \\
!   {\tt mommtfix} & the desired muffin-tin moment for a FSM calculation &
!    real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   The local muffin-tin moments are specified for a subset of atoms, with the
!   list terminated with a blank line. Note that all three components must be
!   specified (even for collinear calculations). See {\tt fixspin}, {\tt taufsm}
!   and {\tt spinpol}.
!
!   \subsection{{\tt mstar}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt mstar} & Value of the effective mass parameter used for adaptive
!    determination of {\tt swidth} & real & $10.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt autoswidth}.
!
!   \subsection{{\tt mustar}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt mustar} & Coulomb pseudopotential, $\mu^*$, used in the
!    McMillan-Allen-Dynes equation & real & $0.15$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is used when calculating the superconducting critical temperature with
!   the formula [Phys. Rev. B 12, 905 (1975)]
!   $$ T_c=\frac{\omega_{\rm log}}{1.2 k_B}\exp\left[\frac{-1.04(1+\lambda)}
!    {\lambda-\mu^*(1+0.62\lambda)}\right], $$
!   where $\omega_{\rm log}$ is the logarithmic average frequency and $\lambda$
!   is the electron-phonon coupling constant.
!
!   \subsection{{\tt ndspem}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt ndspem} & the number of {\bf k}-vector displacements in each direction
!    around {\tt vklem} when computing the numerical derivatives for the
!    effective mass tensor & integer & 1 \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt deltaem} and {\tt vklem}.
!
!   \subsection{{\tt nempty}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nempty} & the number of empty states & integer & 5 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Defines the number of eigenstates beyond that required for charge
!   neutrality. When running metals it is not known {\it a priori} how many
!   states will be below the Fermi energy for each {\bf k}-point. Setting
!   {\tt nempty} greater than zero allows the additional states to act as a
!   buffer in such cases. Furthermore, magnetic calculations use the
!   first-variational eigenstates as a basis for setting up the
!   second-variational Hamiltonian, and thus {\tt nempty} will determine the
!   size of this basis set. Convergence with respect to this quantity should be
!   checked.
!
!   \subsection{{\tt ngridk}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt ngridk } & the {\bf k}-point mesh sizes & integer(3) & $(1,1,1)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   \vspace{5pt}
!   The ${\bf k}$-vectors are generated using
!   $$ {\bf k}=(\frac{i_1}{n_1},\frac{i_2}{n_2},\frac{i_3}{n_3})
!    +{\bf v}_{\rm off}, $$
!   where $i_j$ runs from 0 to $n_j-1$ and $0\le{\bf v}_{{\rm off};j}<1$ for
!   $j=1,2,3$. See also {\tt reducek} and {\tt vkloff}.
!
!   \subsection{{\tt ngridq}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt ngridq } & the phonon {\bf q}-point mesh sizes & integer(3) &
!    $(1,1,1)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Same as {\tt ngridk}, except that this mesh is for the phonon
!   {\bf q}-points. See also {\tt reduceq}.
!
!   \subsection{{\tt nosource}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nosource} & when set to {\tt .true.}, source fields are projected out
!    of the exchange-correlation magnetic field & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Experimental feature.
!
!   \subsection{{\tt nosym}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nosym} & when set to {\tt .true.} no symmetries, apart from the
!    identity, are used anywhere in the code & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt notes}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt notes(i)} & the $i$th line of the notes & string & - \\
!   \hline
!   \end{tabularx}\newline\newline
!   This block allows users to add their own notes to the file {\tt INFO.OUT}.
!   The block should be terminated with a blank line, and no line should exceed
!   80 characters.
!
!   \subsection{{\tt nprad}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nprad} & radial polynomial order & integer & 4 \\
!   \hline
!   \end{tabularx}\newline\newline
!   This sets the polynomial order for the predictor-corrector method when
!   solving the radial Dirac and Schr\"odinger equations, as well as for
!   performing radial interpolation in the plotting routines.
!
!   \subsection{{\tt nseqit}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nseqit} & number of iterations per self-consistent loop using the
!    iterative first-variational secular equation solver & integer & 6 \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt tseqit}.
!
!   \subsection{{\tt nstfsp}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nstfsp} & number of states to be included in the Fermi surface plot
!    file & integer & 6 \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt nwrite}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nwrite} & number of iterations after which {\tt STATE.OUT} is to be
!    written & integer & 0 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Normally, the density and potentials are written to the file {\tt STATE.OUT}
!   only after completion of the self-consistent loop. By setting {\tt nwrite}
!   to a positive integer the file will be written during the loop every
!   {\tt nwrite} iterations.
!
!   \subsection{{\tt optcomp}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt optcomp} & the components of the first- or second-order optical tensor
!   to be calculated & integer(3) & $(1,1,1)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This selects which components of the optical tensor you would like to plot.
!   Only the first two are used for the first-order tensor.
!
!   \subsection{{\tt phwrite}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nphwrt} & number of {\bf q}-points for which phonon modes are to be
!    found & integer & 1 \\
!   \hline
!   {\tt vqlwrt(i)} & the $i$th {\bf q}-point in lattice coordinates & real(3) &
!    $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is used in conjunction with {\tt task}=230. The code will write the
!   phonon frequencies and eigenvectors to the file {\tt PHONON.OUT} for all the
!   {\bf q}-points in the list. The {\bf q}-points can be anywhere in the
!   Brillouin zone and do not have to lie on the mesh defined by {\tt ngridq}.
!   Obviously, all the dynamical matrices have to be computed first using
!   {\tt task}=200.
!
!   \subsection{{\tt plot1d}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt nvp1d} & number of vertices & integer & 2 \\
!   {\tt npp1d} & number of plotting points & integer & 200 \\
!   \hline
!   {\tt vvlp1d(i)} & lattice coordinates for vertex {\tt i} & real(3) &
!    $(0.0,0.0,0.0)\rightarrow(1.0,1.0,1.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Defines the path in either real or reciprocal space along which the 1D plot
!   is to be produced. The user should provide {\tt nvp1d} vertices in lattice
!   coordinates.
!
!   \subsection{{\tt plot2d}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt vclp2d(1)} & first corner (origin) & real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   {\tt vclp2d(2)} & second corner & real(3) & $(1.0,0.0,0.0)$ \\
!   \hline
!   {\tt vclp2d(3)} & third corner & real(3) & $(0.0,1.0,0.0)$ \\
!   \hline
!   {\tt np2d} & number of plotting points in both directions & integer(2) &
!    $(40,40)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Defines the corners of a parallelogram and the grid size used for producing
!   2D plots.
!
!   \subsection{{\tt plot3d}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt vclp3d(1)} & first corner (origin) & real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   {\tt vclp3d(2)} & second corner & real(3) & $(1.0,0.0,0.0)$ \\
!   \hline
!   {\tt vclp3d(3)} & third corner & real(3) & $(0.0,1.0,0.0)$ \\
!   \hline
!   {\tt vclp3d(4)} & fourth corner & real(3) & $(0.0,0.0,1.0)$ \\
!   \hline
!   {\tt np3d} & number of plotting points each direction & integer(3) &
!    $(20,20,20)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Defines the corners of a box and the grid size used for producing 3D plots.
!
!   \subsection{{\tt primcell}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt primcell} & {\tt .true.} if the primitive unit cell should be found
!    & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Allows the primitive unit cell to be determined automatically from the
!   conventional cell. This is done by searching for lattice vectors among all
!   those which connect atomic sites, and using the three shortest which produce
!   a unit cell with non-zero volume.
!
!   \subsection{{\tt radkpt}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt radkpt } & radius of sphere used to determine k-point density & real &
!    $40.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Used for the automatic determination of the {\bf k}-point mesh. If
!   {\tt autokpt} is set to {\tt .true.} then the mesh sizes will be determined
!   by $n_i=\lambda/|{\bf A}_i|+1$.
!
!   \subsection{{\tt readalu}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt readalu} & set to {\tt .true.} if the interpolation constant for
!    LDA+$U$ should be read from file rather than calculated & logical &
!    {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   When {\tt ldapu}=3, the LDA+$U$ energy and potential are interpolated
!   between FLL and AFM. The interpolation constant, $\alpha$, is normally
!   calculated from the density matrix, but can also be read in from the file
!   {\tt ALPHALU.OUT}. This allows the user to fix $\alpha$, but is also
!   necessary when calculating forces, since the contribution of the potential
!   of the variation of $\alpha$ with respect to the density matrix is not
!   computed. See {\tt lda+u}.
!
!   \subsection{{\tt reducebf}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt reducebf} & reduction factor for the external magnetic fields & real &
!    $1.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   After each iteration the external magnetic fields are multiplied with
!   {\tt reducebf}. This allows for a large external magnetic field at the start
!   of the self-consistent loop to break spin symmetry, while at the end of the
!   loop the field will be effectively zero, i.e. infinitesimal. See
!   {\tt bfieldc} and {\tt atoms}.
!
!   \subsection{{\tt reducek}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt reducek} & type of reduction of the $k$-point set & integer & 1 \\
!   \hline
!   \end{tabularx}\newline\newline
!   Types of reduction are defined by the symmetry group used:
!   \newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   0 & no reduction \\
!   1 & reduce with full crystal symmetry group (including non-symmorphic
!    symmetries) \\
!   2 & reduce with symmorphic symmetries only, and any spin rotation has to
!    be the same as the spatial rotation
!   \end{tabularx}
!   See also {\tt ngridk} and {\tt vkloff}.
!
!   \subsection{{\tt reduceq}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt reduceq} & type of reduction of the $q$-point set & integer & 1 \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt reducek} and {\tt ngridq}.
!
!   \subsection{{\tt rgkmax}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt rgkmax} & $R^{\rm MT}_{\rm min}\times\max\{|{\bf G}+{\bf k}|\}$ &
!    real & $7.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This sets the maximum length for the ${\bf G}+{\bf k}$ vectors, defined as
!   {\tt rgkmax} divided by the smallest muffin-tin radius.
!
!   \subsection{{\tt rmtapm}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt rmtapm } & parameters governing the automatic generation of the
!   muffin-tin radii & real(2) & $(0.25, 0.95)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   When {\tt autormt} is set to true, the muffin-tin radii are found
!   automatically from the formula
!   $$ R_i\propto 1+\zeta|Z_i|^{1/3}, $$
!   where $Z_i$ is the atomic number of the $i$th species, $\zeta$ is stored in
!   {\tt rmtapm(1)} and the value which governs the distance between the
!   muffin-tins is stored in {\tt rmtapm(2)}. When {\tt rmtapm(2)} $=1$, the
!   closest muffin-tins will touch.
!
!   \subsection{{\tt scale}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt scale } & lattice vector scaling factor & real & $1.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Scaling factor for all three lattice vectors. Applied in conjunction with
!   {\tt scale1}, {\tt scale2} and {\tt scale3}.
!
!   \subsection{{\tt scale1/2/3}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt scale1/2/3 } & separate scaling factors for each lattice vector &
!    real & $1.0$ \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt scissor}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt scissor} & the scissors correction & real & $0.0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is the scissors shift applied to states above the Fermi energy. Affects
!   DOS, optics and band structure plots.
!
!   \subsection{{\tt scrpath}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt scrpath} & scratch space path & string & {\tt ./} \\
!   \hline
!   \end{tabularx}\newline\newline
!   This is the path to scratch space where the eigenvector file
!   {\tt EIGVEC.OUT} will be written. If the local directory is accessed via a
!   network then {\tt scrpath} can be set to a directory on the local disk, for
!   example {\tt /tmp/}. Note that the forward slash {\tt /} at the end of the
!   string must be included.
!
!   \subsection{{\tt spincore}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt spincore} & set to {\tt .true.} if the core should be spin-polarised
!    & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt spinorb}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt spinorb} & set to {\tt .true.} if a spin-orbit coupling is required
!    & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   If {\tt spinorb} is {\tt .true.}, then a $\boldsymbol\sigma\cdot{\bf L}$
!   term is added to the second-variational Hamiltonian. See {\tt spinpol}.
!
!   \subsection{{\tt spinpol}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt spinpol} & set to {\tt .true.} if a spin-polarised calculation is
!    required & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   If {\tt spinpol} is {\tt .true.}, then the spin-polarised Hamiltonian is
!   solved as a second-variational step using two-component spinors in the
!   effective magnetic field. The first variational scalar wavefunctions are
!   used as a basis for setting this Hamiltonian.
!
!   \subsection{{\tt spinsprl}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt spinsprl} & set to {\tt .true.} if a spin-spiral calculation is
!   required & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Experimental feature for the calculation of spin-spiral states. See
!   {\tt vqlss} for details.
!
!   \subsection{{\tt sppath}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt sppath} & path where the species files can be found & string &
!    {\tt ./} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Note that the forward slash {\tt /} at the end of the string must be
!   included.
!
!   \subsection{{\tt stype}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt stype} & integer defining the type of smearing to be used & integer &
!    $0$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   A smooth approximation to the Dirac delta function is needed to compute the
!   occupancies of the Kohn-Sham states. The variable {\tt swidth} determines
!   the width of the approximate delta function. Currently implemented are
!   \newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   0 & Gaussian \\
!   1 & Methfessel-Paxton order 1, Phys. Rev. B {\bf 40}, 3616 (1989) \\
!   2 & Methfessel-Paxton order 2 \\
!   3 & Fermi-Dirac
!   \end{tabularx}
!   See also {\tt autoswidth} and {\tt swidth}.
!
!   \subsection{{\tt swidth}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt swidth} & width of the smooth approximation to the Dirac delta
!    function & real & $0.01$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt stype} for details.
!
!   \subsection{{\tt tasks}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt task(i) } & the $i$th task & integer & $-1$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   A list of tasks for the code to perform sequentially. The list should be
!   terminated with a blank line. Each task has an associated integer as
!   follows:\newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   -1 & Write out the version number of the code. \\
!   0 & Ground state run starting from the atomic densities. \\
!   1 & Resumption of ground state run using density in {\tt STATE.OUT}. \\
!   2 & Structural optimisation run starting from the atomic densities, with
!    atomic positions written to {\tt GEOMETRY.OUT}. \\
!   3 & Resumption of structural optimisation run using density in
!    {\tt STATE.OUT} but with positions from {\tt elk.in}. \\
!   5 & Ground state Hartree-Fock run. \\
!   10 & Total, partial and interstitial density of states (DOS). \\
!   14 & Plots the smooth Dirac delta and Heaviside step functions used by the
!        code to calculate occupancies. \\
!   15 & Output ${\bf L}$, ${\bf S}$ and ${\bf J}$ total expectation values. \\
!   16 & Output ${\bf L}$, ${\bf S}$ and ${\bf J}$ expectation values for each
!        {\bf k}-point and state in {\tt kstlist}. \\
!   20 & Band structure plot. \\
!   21 & Band structure plot which includes angular momentum characters for
!    every atom. \\
!   25 & Compute the effective mass tensor at the {\bf k}-point given by
!    {\tt vklem}. \\
!   31, 32, 33 & 1/2/3D charge density plot. \\
!   41, 42, 43 & 1/2/3D exchange-correlation and Coulomb potential plots. \\
!   51, 52, 53 & 1/2/3D electron localisation function (ELF) plot. \\
!   61, 62, 63 & 1/2/3D wavefunction plot:
!    $\left|\Psi_{i{\bf k}}({\bf r})\right|^2$. \\
!   72, 73 & 2/3D plot of magnetisation vector field, ${\bf m}({\bf r})$. \\
!   82, 83 & 2/3D plot of exchange-correlation magnetic vector field,
!    ${\bf B}_{\rm xc}({\bf r})$. \\
!   91, 92, 93 & 1/2/3D plot of $\nabla\cdot{\bf B}_{\rm xc}({\bf r})$. \\
!   100 & 3D Fermi surface plot using the scalar product
!    $p({\bf k})=\Pi_i(\epsilon_{i{\bf k}}-\epsilon_{\rm F})$. \\
!   101 & 3D Fermi surface plot using separate bands (minus the Fermi
!    energy). \\
!   110 & Calculation of M\"{o}ssbauer contact charge densities and magnetic
!    fields at the nuclear sites. \\
!   115 & Calculation of the electric field gradient (EFG) at the nuclear
!    sites. \\
!   120 & Output of the momentum matrix elements
!    $\langle\Psi_{i{\bf k}}|-i\nabla|\Psi_{j{\bf k}}\rangle$. \\
!   121 & Linear optical response tensor. \\
!   122 & Magneto optical Kerr effect (MOKE) angle. \\
!   130 & Output matrix elements of the type
!    $\langle\Psi_{i{\bf k+q}}|\exp[i({\bf G+q})\cdot{\bf r}]|
!    \Psi_{j{\bf k}}\rangle$. \\
!   140 & Energy loss near edge structure (ELNES). \\
!   142, 143 & 2/3D plot of the electric field
!    ${\bf E}({\bf r})\equiv\nabla V_{\rm C}({\bf r})$. \\
!   152, 153 & 2/3D plot of
!    ${\bf m}({\bf r})\times{\bf B}_{\rm xc}({\bf r})$. \\
!   162 & Scanning-tunneling microscopy (STM) image. \\
!   190 & Write the atomic geometry to file for plotting with {\sf XCrySDen}
!    and {\sf V\_Sim}. \\
!   200 & Calculation of dynamical matrices on a {\bf q}-point set defined by
!    {\tt ngridq}. \\
!   210 & Phonon density of states. \\
!   220 & Phonon dispersion plot. \\
!   230 & Phonon frequencies and eigenvectors for an arbitrary
!    ${\bf q}$-point. \\
!   240 & Generate the ${\bf q}$-dependent phonon linewidths and electron-phonon
!    coupling constants and write them to file. \\
!   245 & Phonon linewidths plot. \\
!   250 & Eliashberg function $\alpha^2F(\omega)$, electron-phonon coupling
!    constant $\lambda$, and the McMillan-Allen-Dynes critical temperature
!    $T_c$. \\
!   300 & Reduced density matrix functional theory (RDMFT) calculation. \\
!   400 & Calculation of tensor moments and corresponding LDA+U Hartree-Fock
!    energy contributions.
!   \end{tabularx}
!
!   \subsection{{\tt tau0atm}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tau0atm} & the step size to be used for structural optimisation &
!    real & $0.2$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   The position of atom $\alpha$ is updated on step $m$ of a structural
!   optimisation run using
!   $$ {\bf r}_{\alpha}^{m+1}={\bf r}_{\alpha}^m+\tau_{\alpha}^m
!    \left({\bf F}_{\alpha}^m+{\bf F}_{\alpha}^{m-1}\right), $$
!   where $\tau_{\alpha}$ is set to {\tt tau0atm} for $m=0$, and incremented by
!   the same amount if the atom is moving in the same direction between steps.
!   If the direction changes then $\tau_{\alpha}$ is reset to {\tt tau0atm}.
!
!   \subsection{{\tt tauoep}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tauoep} & step length initial value and scaling factors for the OEP
!    iterative solver & real(3) & $(1.0, 0.2, 1.5)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   The optimised effective potential is determined using an interative method.
!   [Phys. Rev. Lett. 98, 196405 (2007)]. At the first iteration the step length
!   is set to {\tt tauoep(1)}. During subsequent iterations, the step length is
!   scaled by {\tt tauoep(2)} or {\tt tauoep(3)}, when the residual is
!   increasing or decreasing, respectively. See also {\tt maxitoep}.
!
!   \subsection{{\tt taufsm}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt taufsm} & the step size to be used when finding the effective magnetic
!   field in fixed spin moment calculations & real & $0.01$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   An effective magnetic field, ${\bf B}_{\rm FSM}$, is required for fixing the
!   spin moment to a given value, $\boldsymbol{\mu}_{\rm FSM}$. This is found by
!   adding a vector to the field which is proportional to the difference between
!   the moment calculated in the $i$th self-consistent loop and the required
!   moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left(
!    \boldsymbol{\mu}^i-\boldsymbol{\mu}_{\rm FSM}\right), $$
!   where $\lambda$ is proportional to {\tt taufsm}. See also {\tt fixspin},
!   {\tt momfix} and {\tt spinpol}.
!
!   \subsection{{\tt tfibs}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tfibs} & set to {\tt .true.} if the IBS correction to the force should
!    be calculated & logical & {\tt .true.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   Because calculation of the incomplete basis set (IBS) correction to the
!   force is fairly time-consuming, it can be switched off by setting
!   {\tt tfibs} to {\tt .false.} This correction can then be included only when
!   necessary, i.e. when the atoms are close to equilibrium in a structural
!   relaxation run.
!
!   \subsection{{\tt tforce}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tforce} & set to {\tt .true.} if the force should be calculated at the
!    end of the self-consistent cycle & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   This variable is automatically set to {\tt .true.} when performing
!   structural optimisation.
!
!   \subsection{{\tt tmomlu}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tmomlu} & set to {\tt .true.} if the tensor moments and the
!   corresponding decomposition of LDA+U energy should be calculated
!   at every iteration of the self-consistent cycle & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   This variable is useful to check the convergence of the tensor moments in
!   LDA+U caculations. Alternatively, with {\tt task} equal to 400, one can
!   calculate the tensor moments and corresponding LDA+U energy contributions
!   from a given density matrix and set of Slater parameters at the end of the
!   self-consistent cycle.
!
!   \subsection{{\tt tseqit}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tseqit} & set to {\tt .true.} if the first-variational secular equation
!    should be solved iteratively & logical & {\tt .false.} \\
!   \hline
!   \end{tabularx}\newline\newline
!   See also {\tt nseqit}.
!
!   \subsection{{\tt tshift}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt tshift} & set to {\tt .true.} if the crystal can be shifted so that the
!    atom closest to the origin is exactly at the origin &
!    logical & {\tt .true.} \\
!   \hline
!   \end{tabularx}
!
!   \subsection{{\tt vklem}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt vklem} & the ${\bf k}$-point in lattice coordinates at which to compute
!    the effective mass tensors & real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt deltaem} and {\tt ndspem}.
!
!   \subsection{{\tt vqlss}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt vqlss} & the ${\bf q}$-vector of the spin-spiral state in lattice
!    coordinates & real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Spin-spirals arise from spinor states assumed to be of the form
!   $$ \Psi^{\bf q}_{\bf k}({\bf r})=
!    \left( \begin{array}{c}
!    U^{{\bf q}\uparrow}_{\bf k}({\bf r})e^{i({\bf k+q/2})\cdot{\bf r}} \\
!    U^{{\bf q}\downarrow}_{\bf k}({\bf r})e^{i({\bf k-q/2})\cdot{\bf r}} \\
!    \end{array} \right). $$
!   These are determined using a second-variational approach, and give rise to a
!   magnetisation density of the form
!   $$ {\bf m}^{\bf q}({\bf r})=(m_x({\bf r})\cos({\bf q \cdot r}),
!    m_y({\bf r})\sin({\bf q \cdot r}),m_z({\bf r})), $$
!    where $m_x$, $m_y$ and $m_z$ are lattice periodic. See also {\tt spinprl}.
!
!   \subsection{{\tt vkloff}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt vkloff } & the ${\bf k}$-point offset vector in lattice coordinates &
!    real(3) & $(0.0,0.0,0.0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   See {\tt ngridk}.
!
!   \subsection{{\tt xctype}}
!   \begin{tabularx}{\textwidth}[h]{|l|X|c|c|}
!   \hline
!   {\tt xctype} & integers defining the type of exchange-correlation functional
!    to be used & integer(3) & $(3,0,0)$ \\
!   \hline
!   \end{tabularx}\newline\newline
!   Normally only the first value is used to define the functional type. The
!   other value may be used for external libraries. Currently implemented are:
!   \newline\newline
!   \begin{tabularx}{\textwidth}[h]{lX}
!   1 & No exchange-correlation funtional ($E_{\rm xc}\equiv 0$) \\
!   2 & LDA, Perdew-Zunger/Ceperley-Alder, {\it Phys. Rev. B} {\bf 23}, 5048
!    (1981) \\
!   3 & LSDA, Perdew-Wang/Ceperley-Alder, {\it Phys. Rev. B} {\bf 45}, 13244
!    (1992) \\
!   4 & LDA, X-alpha approximation, J. C. Slater, {\it Phys. Rev.} {\bf 81}, 385
!    (1951) \\
!   5 & LSDA, von Barth-Hedin, {\it J. Phys. C} {\bf 5}, 1629 (1972) \\
!   20 & GGA, Perdew-Burke-Ernzerhof, {\it Phys. Rev. Lett.} {\bf 77}, 3865
!    (1996) \\
!   21 & GGA, Revised PBE, Zhang-Yang, {\it Phys. Rev. Lett.} {\bf 80}, 890
!    (1998) \\
!   22 & GGA, PBEsol, Phys. Rev. Lett. 100, 136406 (2008) \\
!   26 & GGA, Wu-Cohen exchange (WC06) with PBE correlation, {\it Phys. Rev. B}
!    {\bf 73}, 235116 (2006) \\
!   30 & GGA, Armiento-Mattsson (AM05) spin-unpolarised functional,
!    {\it Phys. Rev. B} {\bf 72}, 085108 (2005) \\
!   100 & {\tt libxc} functionals, the second and third values of {\tt xctype}
!    define the exchange and correlation functionals in the {\tt libxc} library,
!    respectively \\
!   \end{tabularx}
!
!   \section{Contributing to {\sf Elk}}
!   Please bear in mind when writing code for the {\sf Elk} project that it
!   should be an exercise in physics and not software engineering. All code
!   should therefore be kept as simple and concise as possible, and above all it
!   should be easy for anyone to locate and follow the {\sf Fortran}
!   representation of the original mathematics. We would also appreciate the
!   following conventions being adhered to:
!   \begin{itemize}
!   \item Strict {\sf Fortran} 90/95 should be used. Features which are marked
!    as obsolescent in F90/95 should be avoided. These include assigned format
!    specifiers, labeled do-loops, computed goto statements and statement
!    functions.
!   \item Modules should be used in place of common blocks for declaring
!    global variables. Use the existing modules to declare new global variables.
!   \item Any code should be written in lower-case free form style, starting
!    from column one. Try and keep the length of each line to fewer than 80
!    characters using the \& character for line continuation.
!   \item Every function or subroutine, no matter how small, should be in its
!    own file named {\tt routine.f90}, where {\tt routine} is the function or
!    subroutine name. It is recommended that the routines are named so as to
!    make their purpose apparent from the name alone.
!   \item Use of {\tt implicit none} is mandatory. Remember also to define the
!    {\tt intent} of any passed arguments.
!   \item Local allocatable arrays must be deallocated on exit of the routine to
!    prevent memory leakage. Use of automatic arrays should be limited to arrays
!    of small size.
!   \item Every function or subroutine must be documented with the {\sf Protex}
!    source code documentation system. This should include a short \LaTeX\
!    description of the algorithms and methods involved. Equations which need to
!    be referenced should be labeled with {\tt routine\_1}, {\tt routine\_2}
!    etc. The authorship of each new piece of code or modification should be
!    indicated in the {\tt REVISION HISTORY} part of the header. See the
!    {\sf Protex} documentation for details.
!   \item Ensure as much as possible that a routine will terminate the program
!    when given improper input instead of continuing with erroneous results.
!    Specifically, functions should have a well-defined domain for which they
!    return accurate results. Input outside that domain should result in an
!    error message and termination.
!   \item Report errors prior to termination with a short description, for
!    example:
!    \begin{verbatim}
!     write(*,*)
!     write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
!     write(*,'(" for species ",I4)') is
!     write(*,*)
!    \end{verbatim}
!   \item Wherever possible, real numbers outputted as ASCII data should be
!    formatted with the {\tt G18.10} specifier.
!   \item Avoid redundant or repeated code: check to see if the routine you need
!    already exists, before writing a new one.
!   \item All reading in of ASCII data should be done in the subroutine
!    {\tt readinput}. For binary data, separate routines for reading and writing
!    should be used (for example {\tt writestate} and {\tt readstate}).
!   \item Input filenames should be in lowercase and have the extension
!    {\tt .in} . All output filenames should be in uppercase with the extension
!    {\tt .OUT} .
!   \item All internal units should be atomic. Input and output units should be
!    atomic by default and clearly stated otherwise. Rydbergs should not be used
!    under any circumstances.
!   \end{itemize}
!   \subsection{Licensing}
!    Routines which constitute the main part of the code are released under the
!    GNU General Public License (GPL). Library routines are released under the
!    less restrictive GNU Lesser General Public License (LGPL). Both licenses
!    are contained in the file {\tt COPYING}. Any contribution to the code must
!    be licensed at the authors' discretion under either the GPL or LGPL.
!    Author(s) of the code retain the copyrights. Copyright and (L)GPL
!    information must be included at the beginning of every file, and no code
!    will be accepted without this.
!
!EOI
