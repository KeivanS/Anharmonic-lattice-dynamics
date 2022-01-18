module constants
 implicit none
 real(8) :: pi=3.14159265358979323846d0
 complex(8) :: ci=cmplx(0d0,1d0)
 real(8) :: h_plank= 6.62606896d-34
 real(8) :: n_avog= 6.023d23
 real(8) :: k_b = 1.3806504d-23       ! J/K:
 real(8) :: c_light = 2.99792458d+8   ! in (m/s)
 real(8) :: hbar = 1.054571628d-34  !h_plank/2/pi  !
 real(8) :: ee = 1.60217653d-19
 real(8) :: eps0 = 8.854187817d-12
 real(8) :: me = 9.1093826d-31
 real(8) :: uma = 1.66053878d-27 ! kg
 !real(8) :: cnst = sqrt(ee*1d20/uma)/pi/200/c_light ! converts sqrt(eV/A/A/uma) to cm^-1
  real(8) :: cnst= 521.1098918
! real(8) :: ryd= me*ee**4/2/pi/pi/hbar**2/16/pi/pi/eps0/eps0
  real(8) :: ryd= 27.2116
! real(8) :: ab = hbar*hbar/me/ee/ee*4*pi*eps0
  real(8) :: ab = 0.529177
! kbe=8.617343e-05
!****************in a.u.*****************
 !real(8) :: k_b = 1d0
 !real(8) :: hbar = 1d0
!********************************************

 end module constants
 !===========================================================
 ! Module for saved data
module force_constants_module
! maximum number of shells of nearest neighbors
      integer maxneighbors
!     parameter(maxneighbors=18 )
! maximum number of atoms out to maxneighbors
      integer maxatoms,imaxat
!     parameter(maxatoms=2800 )
! op_matrix(k,j,i), matrix for the ith point operator
      double precision op_matrix(3,3,48)
      double precision op_kmatrix(3,3,48)
      integer lattpgcount
! isgopcount, number of operators in space group
      integer isgopcount
! isgop(i), point operation in ith space group operator
      integer isgop(48)
! sgfract(j,i), jth cartesian coordinate of fractional in ith space group
! operator
      double precision sgfract(3,48)
! iatomop(j,i), point operator that takes jth atom into ith atom
      integer, allocatable:: iatomop(:,:)
! atomopfract(k,j,i), kth cartesian coordinate of fractional to accompany
! iatomop
      double precision, allocatable :: atomopfract(:,:,:)
! natoms0, number of atoms in the primitive unit cell
      integer natoms0
! natoms, number of atoms out to maxneighbors
      integer natoms
! atompos(j,i), jth cartesian coordinate of ith atom
   !  double precision atompos(3,maxatoms)
      real(8), allocatable :: atompos(:,:)
! iatomcell(j,i), linear combination of basis vectors of the primitive lattice
! that takes us to the unit cell containing the ith atom
   !  integer iatomcell(3,maxatoms)
      integer, allocatable ::  iatomcell(:,:),iatomcell0(:),iatomneighbor(:,:)
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith
! atom
   !  integer iatomcell0(maxatoms)
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
end module force_constants_module
!*************************************************************************


module test
    use constants
    implicit none
    real(8) :: pi_test
    contains

    subroutine assign_pi
        pi_test=pi
    end subroutine assign_pi

end module test
