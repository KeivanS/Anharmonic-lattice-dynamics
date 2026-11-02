!! module born: dielectric tensor, Born charges, non-analytical term (shared by FOCEX and THERMACOND)
  module born
  use constants , only : pi,r15,ci
  implicit none
!  private
!  public :: non_anal
!  public :: non_analq0
!  public :: set_dynamical_NA_3N3N
!  public :: dyn_na5D_on_grid
!  public :: dyn_na5D_on_fine_grid
!  public :: dyn_coulomb_pure_step
!  public :: compute_phi_na_realspace

  real(r15) epsil(3,3),epsinv(3,3),bref(3,3)
  real(r15), allocatable:: dyn_naq0(:,:,:)
  integer born_flag,np
  complex(r15), allocatable :: dyn_na(:,:,:,:,:)  
  complex(r15), allocatable :: dyn_g(:,:,:,:,:)
  complex(r15), allocatable :: phi_bare(:,:,:,:,:)  ! bare one from fitting 
  complex(r15), allocatable :: phi_sr(:,:,:,:,:) 
  complex(r15), allocatable :: phi_periodic(:,:,:,:,:)  ! from FT of dyn_G phi_periodic(t,R+tp)==sum_L phi_bare(t,R+tp+L)

 contains

  subroutine allocate_fc_dyn(n,nr,ng)
    integer n,nr,ng
    if(allocated(phi_sr)) deallocate (phi_sr)
    if(allocated(phi_bare)) deallocate (phi_bare)
    if(allocated(dyn_g)) deallocate (dyn_g)
    if(allocated(dyn_na)) deallocate (dyn_na)
    allocate(phi_sr(n,n,3,3,nr),phi_bare(n,n,3,3,nr),dyn_g(n,n,3,3,ng),dyn_na(n,n,3,3,ng))
  end subroutine allocate_fc_dyn

end module born
