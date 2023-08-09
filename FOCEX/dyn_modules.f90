!==========================================================
 module om_dos
 integer wmesh,ndyn2
 real(8) wmax,width,etaz
 real(8), allocatable :: dos(:,:),om(:)

   contains

   subroutine set_omdos(la,mesh)
   implicit none
   integer la,mesh,i
   allocate(om(mesh),dos(0:la,mesh))
   ndyn2 = la  ! this is a copy of ndyn
   do i=1,mesh
      om(i) = wmax *(0.0001 + (i-1)/(mesh-1.))
   enddo
   end subroutine set_omdos

   subroutine write_dos(udos)
   implicit none
   integer i,la,udos
   real(8) sumdos

   write(udos,*)'# om(i),i,integrated_dos,total_dos,(dos(la,i),la=1,n))'
   sumdos=0
   do i=1,wmesh
      sumdos=sumdos+dos(0,i)*wmax/wmesh
      write(udos,3)om(i),i,sumdos,(dos(la,i),la=0,min(10,ndyn2))
   enddo
   sumdos=sumdos-dos(0,1)*wmax/wmesh/2d0

3  format(g11.5,2x,i5,99(1x,g10.4))
   end subroutine write_dos

 end module om_dos
!==========================================================
 module eigen
 integer ndyn,nkc2  ! this is a copy of nkc
! eigenval for the coarse mesh (nkc,kpc),

 real(8), allocatable :: velocibz(:,:,:),veloc(:,:,:)
 real(8), allocatable :: eigenval_bs(:,:),eigenval(:,:),eivalibz(:,:)
! sqr of phonon freqs in eV/A/A/m0
 complex(8), allocatable :: eigenvec_bs(:,:,:),eigenvec(:,:,:),grun(:,:),grun_bs(:,:),eivecibz(:,:,:)

    contains

    subroutine allocate_eig(nb,ni,nk) ! nb for band, nk for coarse mesh in FBZ
    integer nb,nk,ni
      allocate( eigenval(nb,nk),eigenvec(nb,nb,nk),veloc   (3,nb,nk), grun(nb,nk) )
      allocate( eivalibz(nb,ni),eivecibz(nb,nb,ni),velocibz(3,nb,ni) )
    end subroutine allocate_eig
!---------------------------------
    subroutine allocate_eig_bs(nb,nk,nv) ! nb for band, nk for band structure mesh
    integer nb,nk,nv
     if (nv.ne.0) then
      allocate( eigenval_bs(nb,nk),eigenvec_bs(nb,nv,nk), grun_bs(nb,nk) ,veloc(3,nb,nk) )
     else
      allocate( eigenval_bs(nb,nk),eigenvec_bs(1,1,1), grun_bs(nb,nk))
     endif
    end subroutine allocate_eig_bs
!---------------------------------
    subroutine deallocate_eig_bs ! nb for band, nk for band structure mesh
      if(allocated(veloc)) deallocate(veloc)
      deallocate( eigenval_bs, eigenvec_bs, grun_bs )
    end subroutine deallocate_eig_bs
!---------------------------------
    subroutine deallocate_eig ! nb for band, nk for coarse mesh in FBZ
      deallocate( eigenval,eigenvec,grun,veloc,eivalibz,eivecibz,velocibz )
    end subroutine deallocate_eig

!---------------------------------
 function onedim(nb,nk)
! calculates the running index onedim
! calculates the running index onedim
! use eigen
 implicit none
 integer nk,nb,onedim
 ! onedim=(nb-1)*nkc2+nk !if outermost loop is j=1,ndyn and innermost loop is i=1,nkc2
 onedim=(nk-1)*ndyn+nb   !if innermost loop is j=1,ndyn and outermost loop is i=1,nkc2
 end function onedim
!---------------------------------
 function nkpt(onedm,nb)
! assumes the outerloop is over kpoints and the innerloop over bands
! use eigen
 implicit none
 integer nkpt,onedm,nb
 nkpt = int((onedm-1)/nb)+1
 end function nkpt
!---------------------------------
 function nband(onedm,nb)
! assumes the outerloop is over kpoints and the innerloop over bands
! use eigen
 implicit none
 integer nband,onedm,nb
 nband = mod(onedm-1,nb)+1
 end function nband
!---------------------------------
 function mysqrt(x)
 implicit none
 real(8) x,mysqrt
 if(x.ge.0) then
    mysqrt=sqrt(x)
 else
    mysqrt=-sqrt(-x)
 endif
 end function mysqrt

 end module eigen
