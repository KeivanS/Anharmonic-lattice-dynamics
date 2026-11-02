!! Small routines shared by FOCEX and THERMACOND (one copy; see ../CLAUDE.md, consolidation).
!! sort: THERMACOND's loop over n passes; FOCEX's n-1 passes left the n-th entry unsorted when
!! maxat > n (the padded arrays of make_grid_shell_ewald).
 subroutine sort(n,r,mcor,maxat)
! sorts the first n elements of array r in ascending order r(mcor(i)) is the ordered array
! make sure to initialize the array r to a huge number if n\=maxat
! r(mcor(1)) < r(mcor(2)) < r(mcor(3)) < ... < r(mcor(n))
 use constants, only : r15
  implicit none
  integer,intent(in):: maxat,n
  integer i,j,temp
  real(r15), intent(in) :: r(maxat)
  integer, intent(out) :: mcor(maxat)

  do i=1,maxat ! n
     mcor(i)=i
  enddo
  do i=1,n    ! n passes: with maxat > n the n-th entry is also the n-th smallest
     do j=i+1,maxat  ! was i+1,n
        if(r(mcor(i)).gt.r(mcor(j))) then
           temp=mcor(i)
           mcor(i)=mcor(j)
           mcor(j)=temp
        endif
     enddo
  enddo
!  r(mcor) is now the ordered array!
 end subroutine sort

 subroutine make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(r15) om

 om = (r1 .dot. (r2 .cross. r3))
 if (abs(om)/length(r1).lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_2PI: volume is zero; check your translation vectors'
    stop
 endif
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=2*pi*(r2 .cross. r3)/om
 g2=2*pi*(r3 .cross. r1)/om
 g3=2*pi*(r1 .cross. r2)/om
 om=abs(om)

 end subroutine make_reciprocal_lattice_2pi

 subroutine inverse_real(a,b,n)
!! inverts a real matrix a(nxn)
 use constants, only : r15
      implicit none
      integer, intent(in) :: n
      integer imax,k,j,i,ii
      real(r15), intent(inout) :: a(n,n)
      real(r15), intent(out) :: b(n,n)
      real(r15) amax,atmp,btmp,amult,div


      b=0d0
      do i=1,n
         b(i,i)=1d0
      enddo

  MAIN_LOOP: do k=1,n

         amax=-1.d-31
         imax=1
         if (k .ne. n) then

           do i=k,n
             if (abs(a(i,k)) .ge. amax) then
               amax=abs(a(i,k))
               imax=i
             endif
           enddo
           if (imax.ne.k) then
             do j=1,n
               atmp=a(imax,j)
               a(imax,j)=a(k,j)
               a(k,j)=atmp
               btmp=b(imax,j)
               b(imax,j)=b(k,j)
               b(k,j)=btmp
             enddo
           endif

         endif

         div=a(k,k)

         if(abs(div).gt.1.d-14) then

           do ii=1,n
             a(k,ii)=a(k,ii)/div
             b(k,ii)=b(k,ii)/div
           enddo
!          a(k,:)=a(k,:)/div
!          b(k,:)=b(k,:)/div

         endif

         do i=1,n
            amult=a(i,k)
            if (i.eq.k) cycle
            do j=1,n
               a(i,j)=a(i,j)-amult*a(k,j)
               b(i,j)=b(i,j)-amult*b(k,j)
            enddo
         enddo

   enddo MAIN_LOOP

 end subroutine inverse_real

 subroutine get_upper_bounds(x1,x2,x3,rcut,maxatm,nsmx)
!! finds the number of grid points in the sphere of radius rcut and m, the needed mesh
!! outputs are nsmx and maxatm: upper bounds in 1D and in 3D (used for 3 loops and allocation)
 use geometry
 use atoms_force_constants, only : natom_prim_cell
 use ios, only : ulog
 use constants, only : r15
 implicit none
 type(vector),intent(in) :: x1,x2,x3
 integer,intent(out) :: nsmx,maxatm ! upper boundary in 1D and in 3D (used for 3 loops and allocation)
 real(r15),intent(in) :: rcut
 integer m1,m2
 real(r15) om0,ls,lmax

 call calculate_volume(x1,x2,x3,om0)
 if(om0.myeq.0d0) then
    write(*,*)'get_upper_bounds: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',x1,x2,x3
    stop
 endif
 ls=om0**0.333333 ! typical unitcell size
 maxatm=(nint(12.56/3d0*(1.3*(rcut+ls))**3/om0)+1)*natom_prim_cell  ! upperbound to the number of atoms in the sphere

 lmax= max(length(x3),length(x1),length(x2))
! find shortest translation vector length
 ls=om0/lmax/lmax
 ls=min(ls,length(x3),length(x1),length(x2))
! write(ulog,3)'GET_UPPER_BOUNDS: lower bound to translation vectors length=',ls

 m1=nint(rcut/ls)+2
 m2=(nint((3/12.56*(maxatm/natom_prim_cell))**0.33333)+1)
 nsmx=max(m1,m2)
! maxp=(2*m+1)**3

 write(ulog,4)'GET_UPPER_BOUNDS:rcut,(rcut+ls)+2,nsmx,maxatom=',rcut,m1,nsmx,maxatm

3 format(a,9(1x,g14.7))
4 format(a,f9.4,1x,3i8)

 end subroutine get_upper_bounds

 function delta_k(i,j) result(k)  ! kronecker delta for integers
 use constants
 implicit none
 integer i,j,k

 if (i.eq.j) then 
   k=1
 else
   k=0
 endif

 end function delta_k

     FUNCTION my_erfc(x) result(S15ADF)
     use constants, only : r15
!     MARK 5A REVISED - NAG COPYRIGHT 1976
!     MARK 5C REVISED
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!     COMPLEMENT OF ERROR FUNCTION ERFC(X)
!
!     .. Scalar Arguments ..
      real(r15), intent(in) :: X
!     .. Local Scalars ..
      real(r15)   T, XHI, XLO, Y, s15adf
!     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP
!     .. Data statements ..
!     PRECISION DEPENDENT CONSTANTS
! 08   DATA XLO/-4.5D0/
! 12   DATA XLO/-5.25D0/
! 14   DATA XLO/-5.75D0/
      DATA XLO/-6.25D0/
! 18   DATA XLO/-6.5D0/
!
!     RANGE DEPENDENT CONSTANTS
      DATA XHI/ 2.66D+1 /
!     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
! R1   DATA XHI/13.0D0/
! R2   DATA XHI/9.5D0/
! R3   DATA XHI/13.0D0/
! R4   DATA XHI/25.0D0/
! R5   DATA XHI/26.0D0/
!     .. Executable Statements ..
!
!     TEST EXTREME EXITS
      IF (X.GE.XHI) then
           s15adf=0
           return
      endif
      IF (X.LE.XLO) then
           s15adf=2
           return
      endif
!
!     EXPANSION ARGUMENT
      T = 1.0D0 - 7.5D0/(ABS(X)+3.75D0)
!
!      * EXPANSION (0021) *
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 08E
! 08   Y = (((((((((((+3.1475326D-5)*T-1.3874589D-4)*T-6.4127909D-6)
! 08  *    *T+1.7866301D-3)*T-8.2316935D-3)*T+2.4151896D-2)
! 08  *    *T-5.4799165D-2)*T+1.0260225D-1)*T-1.6357229D-1)
! 08  *    *T+2.2600824D-1)*T-2.7342192D-1)*T + 1.4558972D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 12E
! 12   Y = ((((((((((((((((-4.21661579602D-8*T-8.63384346353D-8)
! 12  *    *T+6.06038693567D-7)*T+5.90655413508D-7)
! 12  *    *T-6.12872971594D-6)*T+3.73223486059D-6)
! 12  *    *T+4.78645837248D-5)*T-1.52546487034D-4)
! 12  *    *T-2.55222360474D-5)*T+1.80299061562D-3)
! 12  *    *T-8.22062412199D-3)*T+2.41432185990D-2)
! 12  *    *T-5.48023263289D-2)*T+1.02604312548D-1)
! 12  *    *T-1.63571895545D-1)*T+2.26008066898D-1)
! 12  *    *T-2.73421931495D-1)*T + 1.45589721275D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 14E
! 14   Y = (((((((((((((((-2.2356173494379D-9
! 14  *    *T+4.5302502889845D-9)*T+2.5918103316137D-8)
! 14  *    *T-6.3684846832869D-8)*T-1.7642194353331D-7)
! 14  *    *T+6.4907607131235D-7)*T+7.4296952017617D-7)
! 14  *    *T-6.1758018478516D-6)*T+3.5866167916231D-6)
! 14  *    *T+4.7895180610590D-5)*T-1.5246364229106D-4)
! 14  *    *T-2.5534256252531D-5)*T+1.8029626230333D-3)
! 14  *    *T-8.2206213481002D-3)*T+2.4143223946968D-2)
! 14  *    *T-5.4802326675661D-2)*T+1.0260431203382D-1
! 14   Y = (((Y*T-1.6357189552481D-1)*T+2.2600806691658D-1)
! 14  *    *T-2.7342193149541D-1)*T + 1.4558972127504D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+3.328130055126039D-10        &
     &    *T-5.718639670776992D-10)*T-4.066088879757269D-9)       &
     &    *T+7.532536116142436D-9)*T+3.026547320064576D-8)       &
     &    *T-7.043998994397452D-8)*T-1.822565715362025D-7)       &
     &    *T+6.575825478226343D-7)*T+7.478317101785790D-7)       &
     &    *T-6.182369348098529D-6)*T+3.584014089915968D-6)       &
     &    *T+4.789838226695987D-5)*T-1.524627476123466D-4)       &
     &    *T-2.553523453642242D-5)*T+1.802962431316418D-3)       &
     &    *T-8.220621168415435D-3)*T+2.414322397093253D-2
      Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)       &
     &    *T-1.635718955239687D-1)*T+2.260080669166197D-1)       &
     &    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 18E
! 18   Y = (((((((((((((((-1.58023488119651697D-11
! 18  *    *T-4.94972069009392927D-11)*T+1.86424953544623784D-10)
! 18  *    *T+6.29796246918239617D-10)*T-1.34751340973493898D-9)
! 18  *    *T-4.84566988844706300D-9)*T+9.22474802259858004D-9)
! 18  *    *T+3.14410318645430670D-8)*T-7.26754673242913196D-8)
! 18  *    *T-1.83380699508554268D-7)*T+6.59488268069175234D-7)
! 18  *    *T+7.48541685740064308D-7)*T-6.18344429012694168D-6)
! 18  *    *T+3.58371497984145357D-6)*T+4.78987832434182054D-5)
! 18  *    *T-1.52462664665855354D-4)*T-2.55353311432760448D-5
! 18   Y = ((((((((Y*T+1.80296241673597993D-3)
! 18  *    *T-8.22062115413991215D-3)
! 18  *    *T+2.41432239724445769D-2)*T-5.48023266949776152D-2)
! 18  *    *T+1.02604312032198239D-1)*T-1.63571895523923969D-1)
! 18  *    *T+2.26008066916621431D-1)*T-2.73421931495426482D-1)*T +
! 18  *     1.45589721275038539D-1
!
      S15ADF = EXP(-X*X)*Y
      IF (X.LT.0.0D0) then
            S15ADF = 2.0D0 - S15ADF
      endif

      END function my_erfc !s15adf
