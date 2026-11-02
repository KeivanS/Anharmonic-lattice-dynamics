!==============================================================
 module constants
 implicit none
!--! specific precisions, usually same as real and double precision
! integer, parameter :: r6 = selected_real_kind(6)  ! for single precision
integer, parameter,public :: r15 = kind(1.0d0)  ! for double precision
! integer, parameter,public :: r15 = selected_real_kind(15) ! should work on any compiler
integer, parameter,public :: c15 = 8
!integer, parameter :: c6 = selected_complex_kind(6)
! integer, parameter :: dp=8 !(for nag compilers use kd=2 for double precision; 8 for gfortran)
! complex(kind=8) :: ci=dcmplx(0d0,1d0)
 complex(r15) :: ci=cmplx(0d0,1d0)
 real(r15) :: pi=3.14159265358979323846d0
 real(r15) :: h_plank= 6.62606896d-34
 real(r15) :: n_avog= 6.023d23
 real(r15) :: k_b = 1.3806504d-23       ! J/K
 real(r15) :: c_light = 2.99792458d+8   ! in (m/s)
 real(r15) :: hbar = 1.054571628d-34  !h_plank/2/pi
 real(r15) :: ee = 1.60217653d-19
 real(r15) :: eps0 = 8.854187817d-12
 real(r15) :: me = 9.1093826d-31
 real(r15) :: uma = 1.66053878d-27 ! kg
 real(r15) :: cnst= 521.1098918 ! sqrt(ee*1d20/uma)/pi/200/c_light ! converts sqrt(eV/A/A/uma) to cm^-1
 real(r15) :: ryd= 27.2116  !me*ee**3/hbar**2/16/pi/pi/eps0/eps0 (in eV)
 real(r15) :: ab = 0.529177  ! hbar*hbar/me/ee/ee*4*pi*eps0
 real(r15) :: kbe=8.617343e-05  ! = value of 1 Kelvin in eV
 real(r15) :: eps0scale =6.94461587186d-2 ! (4*pi*eps0)/ (ee*1d10)

 end module constants
!==============================================
