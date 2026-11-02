!! dielectric.params: born_flag, the dielectric tensor and the Born charges; the inverse
!! dielectric tensor and the acoustic sum rule on the charges (shared by FOCEX and THERMACOND;
!! moved from FOCEX/ios.f90).  Needs natom_prim_cell and atom0 allocated.
 subroutine read_dielectric
!! reads born effective charge and dielectric constant from dielectric.params
 use born
 use atoms_force_constants
 use ios
 use constants, only : r15
 implicit none
 integer i,j,k
 real(r15) asr,aux(3,3)

 open(uborn,file='dielectric.params',status='old')
! allocate(zeu(3,3,natom_prim_cell))

 read(uborn,*) born_flag   ! if 0 use default
 if (born_flag.le.0) then
    write(ulog,*)' Born_flag  =',born_flag,' Born charges are not subtracted'
    write(6   ,*)' Born_flag  =',born_flag,' Born charges are not subtracted'
 elseif(born_flag.eq.1) then
    write(ulog,*)' Born_flag  =',born_flag,' Ewald force will be subtracted from input forces'
 else
    write(ulog,*)' Born_flag  =',born_flag,' non-analytical force will be subtracted from input forces'
    write(6   ,*)' Born_flag  =',born_flag,' non-analytical force will be subtracted from input forces'
 endif
! user's rule (2026-11-01): use 0 (non-polar) or 2 (polar); 1 only for testing; > 2 not recommended
 if (born_flag.gt.2) then
    write(ulog,*)'WARNING: born_flag=',born_flag,' > 2 is not recommended: use 0 or 2 (1 only for testing)'
    write(6   ,*)'WARNING: born_flag=',born_flag,' > 2 is not recommended: use 0 or 2 (1 only for testing)'
 endif

 do i=1,3
    read(uborn,*)(epsil(i,j),j=1,3)
 end do

 aux=epsil
 call inverse_real(aux,epsinv,3)

 do k=1,natom_prim_cell
    do i=1,3
       read(uborn,*)(atom0(k)%charge(i,j),j=1,3)  !(zeu(i,j,k),j=1,3)
    end do
 end do
 close(uborn)

! impose ASR on Born charges: extra-charge is divided equally between atoms
 do i=1,3
 do j=1,3
      asr=0
      do k=1,natom_prim_cell
         asr=asr+atom0(k)%charge(i,j)  !zeu(i,j,k)
      end do
      asr=asr/natom_prim_cell
      atom0(:)%charge(i,j)=atom0(:)%charge(i,j) - asr
 end do
 end do

 write(ulog,*)' READ dielectric  constant='
 do i=1,3
    write(ulog,3)epsil(i,:)
 enddo
 write(ulog,*)' INVERSE dielectric  constant='
 do i=1,3
    write(ulog,3)epsinv(i,:)
 enddo

 do k=1,natom_prim_cell
 write(ulog,*)' ASR-enforced Born charges for atom ',k,atom0(k)%name
 do i=1,3
    write(ulog,3)atom0(k)%charge(i,:)
 enddo
 enddo

3 format(9(1x,f11.6))

 end subroutine read_dielectric
