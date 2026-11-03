 subroutine set_band_structure 
 use svd_stuff, only : map
 use kpoints
 use atoms_force_constants, only : natom_prim_cell,atom0 ,op_kmatrix
 use ios, only : ulog , uibs, uband, ugrun, write_out
 use eigen 
 use constants, only : r15,ee,pi,eps0
 use lattice, only : volume_r0,g0ws26,r01,r02,r03,primitivelattice
 use mech
 use fourier, only : nrgrid,rgrid,rws_weights ,ggrid,nggrid ,gws_weights,nreg,rgridreg
 use om_dos
 use params, only : tmin,tmax,include_fc,coef
 use born
 use tetrahedron
 implicit none
 integer i,j,k,uio,tau,taup,l,ll,narms,ier
 real(r15) vkms(3),vavg(3),qr(3),q2(3),sf,etot,free,cv,cp,entropy,eq_vol,tempk,cvte
 real(r15) rotc(3,3),atldc(3,3,3,3),elasticc(6,6),atld_sym(3,3,3,3),asym_c
 character(len=12) csys
 character(len=8) vtest
 real(r15) deteps3,grune,boft,boft0,epsilon0(3,3),velocfbz(3),mat(3,3)
 complex(r15) dyn_coul(3,3),ddn(6,6,3)
 real(r15), allocatable :: aux(:),array(:) ,eivaux(:),velaux(:,:)
 complex(r15), allocatable :: evcaux(:,:) 
 real(r15) kvecstar(3,48),tempi(3,3),v3(3),v4(3),velred(3)
 integer kvecop(48),f,total_arms
 integer, allocatable :: visit_count(:)
 real tim
 logical found

  rho_SI=sum(atom0(:)%mass)*uma/volume_r0*1d30   ! density in kg/m^3
  debug_kpoints=.true.

!-------------------------------------------------------------------------------------------
! do a band structure calculation along the symmetry directions defined in kpbs.params

  call make_kp_bs

  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen

  call get_environment_variable('FOCEX_VELTEST',vtest)
  if (vtest(1:1).eq.'1') open(330,file='veltest.dat')   ! finite-difference velocity check (get_frequencies)

  write(*,*)' calling get_freqs'
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc_bs,1)

 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' DYNAMICAL after get_frequencies (BS)   TIME  IS ',tim

! Band structure projection sorting, added by Bolin
  write(*,*) "Projection Band Sorting for band structures!"
  call band_sort_bs(nkp_bs,ndyn,kp_bs,eigenval_bs,eigenvec_bs,veloc_bs)

 
  open(uband,file='bs_freq.dat')
  call write_all_eigenvalues(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,veloc_bs,uband)
  close(uband)

  open(uband,file='bands.dat')
  write(uband,*)'# i,dk_bs(i),kp_bs(:,i),eival_bs(:,i)),length(veloc_bs(:,l,i))/1000'
  do i=1,nkp_bs
     write(uband,6)i,dk_bs(i),kp_bs(:,i),eigenval_bs(:,i),(length(veloc_bs(:,l,i))/1000,l=1,ndyn)
  enddo
  close(uband)

 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' DYNAMICAL after bandsort,              TIME  IS ',tim

! ------------ Mechanical properties ------------------------------------------------------
  uio=321
  open(uio,file='mech.dat')
  call allocate_mech(ndyn)
  if(ndyn.gt.3) then ! only makes sense for non-bravais lattices
    allocate(gama(ndyn-3,ndyn-3))
    call get_phi_zeta_Xi(uio) 
    call residuals (uio) 
    call mechanical(uio) 
  else
    call mechanical0(atld0,uio)
  endif

! the elastic tensor in the crystal's conventional frame (symmetry axes along x,y,z: the fewest
! independent constants); the frame-dependent quantities below (Zener ratio, sound speeds along
! 100, 110, 111) are taken there. The Hill moduli do not depend on the frame.
  call crystal_frame(rotc,csys)
  call symmetrize_tensor4(atld0,atldc,asym_c)     ! average over the crystal's point group
  write(ulog,*)'CRYSTAL FRAME: largest symmetry-forbidden part of the elastic tensor removed (GPa)=',asym_c
  write(uio ,*)'CRYSTAL FRAME: largest symmetry-forbidden part of the elastic tensor removed (GPa)=',asym_c
  atld_sym=atldc
  call rotate_tensor4(rotc,atld_sym,atldc)
  call convert_to_voigt(atldc,elasticc)
  call report_crystal_frame(uio,rotc,csys,elasticc)

! bulk mod(cubic)=B=(C11+2*C12+P)/3 ; shear modulus: G=(C11-C12+3C44)/5  
! Young Modulus E=9BG/(3B+G) ; Poisson=(3B-2G)/(6B+2G)

  write( uio,8)'B0P=dB/dP from cubic_elastic (will be overwritten) =',b0p

  vbulkmod = sum(elastic(1:3,1:3))/9
  write(ulog,8)'Voigt Bulk modulus(GPa) from elastic tensor=',vbulkmod
  write( uio,8)'Voigt Bulk modulus(GPa) from elastic tensor=',vbulkmod
  rbulkmod = 1/sum(compliance(1:3,1:3))
  write(ulog,8)'Reuss Bulk modulus(GPa) from complnc tensor=',rbulkmod
  write(uio ,8)'Reuss Bulk modulus(GPa) from complnc tensor=',rbulkmod
  hbulkmod=(vbulkmod+rbulkmod)/2
  write(ulog,8)'Hill  Bulk modulus(GPa) from        average=',hbulkmod
  write(uio ,8)'Hill  Bulk modulus(GPa) from        average=',hbulkmod

! B0' = dB/dP at T=0 from FC2 and FC3 (FOCEX/notes/generalized_huang.tex, Eq. Bp); overwrites cubic_elastic's
  if(include_fc(3).ne.0) call pressure_derivative_b(uio)

  vshearmod =(elastic(1,1)+elastic(2,2)+elastic(3,3)  &
&            -elastic(1,2)-elastic(2,3)-elastic(3,1)  &
&         +3*(elastic(4,4)+elastic(5,5)+elastic(6,6)) )/15d0
  write(ulog,8)'Voigt shear modulus(GPa) from elastic tensor=',vshearmod
  write( uio,8)'Voigt shear modulus(GPa) from elastic tensor=',vshearmod
  rshearmod = 15/(4*(compliance(1,1)+compliance(2,2)+compliance(3,3))  &
&                -4*(compliance(1,2)+compliance(2,3)+compliance(3,1))  &
&                +3*(compliance(4,4)+compliance(5,5)+compliance(6,6)) )
  write(ulog,8)'Reuss shear modulus(GPa) from complnc tensor=',rshearmod
  write(uio ,8)'Reuss shear modulus(GPa) from complnc tensor=',rshearmod
  hshearmod=(vshearmod+rshearmod)/2
  write(ulog,8)'Isotropic Hill  shear modulus(GPa) from    average=',hshearmod
  write(uio ,8)'Isotropic Hill  shear modulus(GPa) from    average=',hshearmod
  write(ulog,8)'Shear modulus(GPa) 1/S44,55,66  compliance tensor =',  &
 &        1/(compliance(4,4)),  1/(compliance(5,5)),  1/(compliance(6,6))  
  write(uio ,8)'Shear modulus(GPa) 1/S44,55,66  compliance tensor =',  &
 &        1/(compliance(4,4)),  1/(compliance(5,5)),  1/(compliance(6,6))  


  Youngmod=9*hbulkmod*hshearmod/(3*hbulkmod+hshearmod)
  write(ulog,8)'Average (from Hill) Young modulus(GPa) =9BG/(3B+G)=',youngmod
  write(uio ,8)'Average (from Hill) Young modulus(GPa) =9BG/(3B+G)=',youngmod
  write(ulog,8)'single crystal Young modulus(GPa) from 1/S11,22,33=',1/compliance(1,1),  &
&     1/compliance(2,2),1/compliance(3,3)
  write(uio ,8)'single crystal Young modulus(GPa) from 1/S11,22,33=',1/compliance(1,1),  &
&     1/compliance(2,2),1/compliance(3,3)


  poisson=(1.5*hbulkmod-hshearmod)/(3*hbulkmod+hshearmod)
  write(ulog,8)'Hill Poisson ratio=(1.5B-G)/(3B+G) for isotropic  =',poisson
  write(uio ,8)'Hill Poisson ratio=(1.5B-G)/(3B+G) for isotropic  =',poisson
  poisson = -compliance(1,2)/compliance(1,1)
  write(ulog,8)'Poisson ratio from -S12/S11 for single crystals   =',poisson
  write(uio ,8)'Poisson ratio from -S12/S11 for single crystals   =',poisson

  anisotropy_ratio=2*elasticc(4,4)/(elasticc(1,1)-elasticc(1,2))   ! crystal frame
  write(ulog,8)'Anisotropy ratio from 2C44/(C11-C12)              =',anisotropy_ratio
  write(uio ,8)'Anisotropy ratio from 2C44/(C11-C12)              =',anisotropy_ratio

  vl_iso=sqrt((hbulkmod+4/3d0*hshearmod)*1d9/rho_si) / 1000d0
  if (hshearmod.gt.0) then
     vt_iso=sqrt(1d9* hshearmod /rho_SI) /1000d0
  else
     vt_iso=1d-1  
  endif
  v_mean=1d0/((1d0/vl_iso**3+2d0/vt_iso**3)/3d0)**0.333333
  t_debye=hbar/k_b*(6*pi*pi*natom_prim_cell/volume_r0)**0.33333 *1d10 * v_mean*1000
  write(ulog,8)'isotropic V_l, V_t, V_mean (km/s) = ',vl_iso, vt_iso,v_mean
  write( uio,8)'isotropic V_l, V_t, V_mean (km/s) = ',vl_iso, vt_iso,v_mean
  xagne=vt_iso/vl_iso
  write(ulog,8)' xagne, T_D(K) =',xagne,t_debye
  write(uio ,8)' xagne, T_D(K) =',xagne,t_debye


  call elastic_longwave(uio,rotc)   ! C_ij from the sound speeds, E = 0 and with the NA term
  close(uio)

! now calculating average sound speeds and Debye temperature along 100, 110, and 111
 vavg = 0
 do f=0,2  ! loop for 100,110,111 for f=0,1,2 respectively

    q2(1)=1
    q2(2)=nint((f+0.6)/2) 
    q2(3)=nint((f-0.1)/2) 

    call sound_speeds(q2,atldc,vkms)     ! directions of the crystal frame
    vavg = vavg + vkms
 enddo
 vavg = vavg/3

 write(ulog,8)'vavg avged over 100,110,111 = ',vavg

 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' After mechanical and sound speeds,      TIME IS ',tim

  if(include_fc(3).ne.0) then
     write(*,*)' calling gruneisen'
     open(ugrun,file='bs_grun.dat')
     call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,grun_bs,ugrun,1)
     close(ugrun)
  endif

  call deallocate_eig_bs
  call deallocate_kp_bs  ! also deallocates dk_bs


   
! ---------------- Now defined the IBZ and weights, and eigenvalues within the IBZ -----------
! now take kpoints within the IBZ for lattice summations

  call set_omdos(ndyn,wmesh)

  allocate(kpc(3,nkc),wk(nkc),mapinv(nkc),mapibz(nkc))

! allocates eig and afunc(n1*n2*n3), dos_tet(wmesh)
  call allocate_tetra(nc(1),nc(2),nc(3),wmesh)

  call make_kp_reg_tet

  call get_weights(nkc,kpc,mapinv,mapibz)   ! generates kibz and weights
  call allocate_eig_ibz(ndyn,nibz) ;  allocate(dk_bs(nibz)) ; grunibz=0; dk_bs=0
  write(*,*)' calling get_freq in IBZ with nibz=', nibz 
  call get_frequencies(nibz,kibz,dk_bs,ndyn,eivalibz,ndyn,eivecibz,velocibz,0)

  open(uband,file='ibz_bands.dat')
  call write_all_eigenvalues(nibz,kibz,dk_bs,ndyn,eivalibz,velocibz,uband)
  close(uband)
  close(330)

  if(include_fc(3).ne.0) then
     open(ugrun,file='ibz_grun.dat')
     call gruneisen(nibz,kibz,dk_bs,ndyn,eivalibz,eivecibz,grunibz,ugrun,0)
     close(ugrun)
  endif

  allocate(aux(nibz),array(wmesh)); deallocate(dk_bs)
  dos = 0
  do j=1,ndyn
     aux = mysqrt(eivalibz(j,:)) * cnst
     call calculate_dos(nibz,aux ,wibz,wmesh,om,array)  ! Gaussian-smeared DOS
     dos(j,:)=array
     dos(0,:)=dos(0,:)+array
  enddo
  open(udos,file='ibz_dos.dat')
  call write_dos(udos)
  close(udos)

  deallocate(aux,array)  
  allocate(eivaux(ndyn),velaux(3,ndyn),evcaux(ndyn,ndyn))
  q2=normal*1d-6 ! almost zero!
  call get_freq(q2,ndyn,velaux,eivaux,evcaux)
  call dielectric(ndyn,eivaux,evcaux,wmesh,om,epsilon0) ! only need Gamma point spectrum
  call gamma_mode_symmetry(ndyn)                         ! IR / Raman activity at Gamma (gamma_modes.dat)
  call write_out(ulog,' Ion-clamped dielectric constant ',epsil)
  call write_out(ulog,' Total dielectric constant along normal ',epsilon0)
  deallocate(eivaux,evcaux,velaux)

 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' DYNAMICAL after IBZ DOS,               TIME  IS ',tim

! -------------- now DOS using the tetrahedron method: first use symmetry(in mapibz) to get eival(FBZ)
  dos = 0; afunc=1d0 ! weight in front of the delta function
  do i=1,wmesh
  do j=1,ndyn
     do k=1,nkc
        eig(k) = mysqrt(eivalibz(j,mapibz(k))) * cnst
     enddo
     call tet_sum(om(i),nkc,eig,afunc,dos(j,i))  ! tetrahedron DOS
     dos(0,i)=dos(0,i)+dos(j,i)
  enddo
  enddo
  open(udos,file='fbz_dos.dat')
  call write_dos(udos)
  close(udos)

! -------------- THERMAL CONDUCTANCE : count number of modes along "normal" direction  
! 
! first normalize the normal vector
  normal=normal/(length(normal)+1d-12)
  allocate(visit_count(nkc))
  dos = 0
  do j=1,-1 !ndyn  ! sum over bands
     eig=0d0; afunc=0d0; visit_count=0; total_arms=0
     do k=1,nibz   ! sum over k vectors in the irreducible wedge
        call getstar(kibz(:,k),primitivelattice,narms,kvecstar,kvecop)
        if(debug_kpoints) then
          total_arms=total_arms+narms
          if(abs(narms/dble(nkc)-wibz(k)) .gt. 1d-6) then
             write(*,5) 'ERROR: weight mismatch for IBZ point ', k,narms/real(nkc),wibz(k),narms/real(nkc)-wibz(k)
          endif
        endif
        do ll=1,narms  ! sum over the stars of kibz obtained by symmetry
           vkms = matmul( op_kmatrix(:,:,kvecop(ll)) , velocibz(:,j,k) )
           q2   = matmul( op_kmatrix(:,:,kvecop(ll)) , kibz(:,k) )
! id this kibz and arm with actual k
           found=.false.
           do l=1,nkc
              if(length(kpc(:,l)-q2).lt.1d-6) then
                 found=.true.
                 exit
              endif
           enddo
           if(.not. found) then
              write(*,5)'could not find k corresponding to kibz,arm ',ll,kibz(:,k)
           !  stop
           else
              visit_count(l)=visit_count(l)+1
           endif
           if (l.le.nkc) then
              eig(l) = mysqrt(eivalibz(j,k)) * cnst
              afunc(l)=dot_product(normal,vkms) ! weight in front of the delta function
           else
              write(*,*) 'Could not find l; stopping!'
           !  stop
           endif
        enddo 
     enddo

! Now check to see if any l point is visited once and only once
   if(debug_kpoints) then
     do l = 1, nkc
        if (visit_count(l) .eq. 0) then
           write(*,*) 'ERROR: k-point ', l, ' never visited: ', kpc(:,l)
        elseif (visit_count(l) .gt. 1) then
           write(*,*) 'ERROR: k-point ', l, ' visited ', visit_count(l), ' times: ', kpc(:,l)
        endif
     enddo
! Summary check
     if (total_arms .ne. nkc) then
        write(*,*) 'ERROR: sum of star sizes = ', total_arms, ' != nkc = ', nkc
        stop
     endif
     if (sum(visit_count) .ne. nkc) then
        write(*,*) 'ERROR: total visits = ', sum(visit_count), ' but nkc = ', nkc
        stop
     endif
     if (any(visit_count .ne. 1)) then
        write(*,*) 'ERROR: visit_count is not uniformly 1'
        stop
     endif
     write(*,*) 'OK: all ', nkc, ' k-points visited exactly once'
   endif

     do i=1,wmesh
        call tet_sum(om(i),nkc,eig,afunc,dos(j,i))  ! tetrahedron DOS
        dos(0,i)=dos(0,i)+dos(j,i)
     enddo
  enddo

  open(udos,file='conductance.dat')
  call write_dos(udos)
  close(udos)

2 format(a,2i5,9(1x,f10.5))
  call deallocate_tetra


 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' DYNAMICAL after FBZ DOS_TETRA ,        TIME  IS ',tim

! ------------- Thermodynamic properties within the QHA
  open(utherm,file='thermo_QHA.dat')
  open(uio,file='temp_free.dat')
  write(utherm,'(a)')"# T(K) , E_eq(kJ/mol) , F_eq(kJ/mol),V_eq/V0,strain , Cv, "// &
 &      "Cp(J/K.mol), S(J/K.mol) , lin_CTE=dlnL/dT(1/K),Grun , B_T(T)/B_T(T->0), B_T(GPa), B0(Reuss,static,GPa)"
  write(uio,*)'# temperature, strain , E_free(kJ/mol) , Pressure (GPa) , i(P=0),B0,B_T(strain),gruneisen,(dP/dT)_V,'// &
 &   ' F0,P0 (static), F_agne,P_agne (harmonic-only term), F_ph,P_ph (phonons); F in kJ/mol, P in GPa'
  tempk=tmin 
  temploop: do i=1,500 ! this allows to from 5K to Tmax 
     if(tempk.lt.40) then ! at t<150 put more T values differing by 5 K
        tempk=tempk+2
     elseif(tempk.lt.150) then ! at t<150 put more T values differing by 5 K
        tempk=tempk+5
     else
        tempk=tempk+50  
     endif
     if (tempk .gt. Tmax) exit temploop

     boft=hbulkmod

     call QHA_free_energy_vs_strain(nibz,wibz,ndyn,eivalibz,real(grunibz),tempk, &
&              free,etot,cv,cp,eq_vol,entropy,cvte,grune,boft,uio)  ! eq_vol is really eq_vol/volum_r0

! B_T normalized like V_eq (to the lowest temperatures, T < 8 K), so column 11 goes to 1 at T -> 0 (user,
! 2026-11-02; was B_T/B0(Reuss) at V0, /= 1 at T=0 because of zero-point and residual pressure)
     if (i.eq.1 .or. tempk.lt.8) boft0=boft
     write(utherm,7) tempk,etot,free,eq_vol,(eq_vol)**0.33333-1,cv,cp,entropy,cvte/3,grune,boft/boft0, &
 &                   boft,rbulkmod
  enddo temploop

 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' DYNAMICAL after QHA THERMO,            TIME  IS ',tim

  close(uio)
  close(utherm)
5 format(a,i5,99(1x,g11.4))
6 format(i5,9999(1x,g11.4))   ! bands.dat: one line per k-point (99 wrapped beyond 47 bands)
7 format(99(1x,g11.5))
8 format(a,9(1x,g11.4))

  call deallocate_eig_ibz 
  deallocate(kibz,mapinv,mapibz,kpc,wk,wibz,dos,om)  

 end subroutine set_band_structure
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg,ipath)
! also outputs the group velocities from HF theorem (exact) in km/s
! Where each point is evaluated (the k-point itself; Gamma as a directional limit, the band-path
! segment for ipath=1, else uinc) is set by gamma_limit and freq_vel (module phonon_rules, common/,
! shared with THERMACOND). FOCEX_VELTEST=1 also writes the finite-difference check of the
! velocities (veltest.dat).
 use params
 use constants
 use born
 use geometry
 use lattice, only : volume_g0
 use atoms_force_constants, only : natom_prim_cell, atom0
 use ios, only : ulog
 use phonon_rules, only : gamma_limit, freq_vel, freq_routine
 implicit none
 integer, intent(in) :: nkp,ndn,nv,ipath ! no of wanted eivecs
 real(r15), intent(in) :: kp(3,nkp),dk(nkp)
 real(r15), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(r15), intent(out) :: eivec(ndn,ndn,nkp)
 procedure(freq_routine) :: get_freq
 integer i
 real(r15) eivl(ndn),vel(3,ndn),q(3),dir(3),gsize,mass(natom_prim_cell)
 complex(r15) eivc(ndn,ndn)
 character(len=8) veltest

 call get_environment_variable('FOCEX_VELTEST',veltest)
 gsize=volume_g0**0.33333333
 mass=atom0(1:natom_prim_cell)%mass

 kloop: do i=1,nkp
    call gamma_limit(i,nkp,kp,ipath,q,dir)
    call freq_vel(get_freq,q,dir,gsize,natom_prim_cell,mass,ndn, &
 &                vg(:,:,i),eival(:,i),eivec(:,1:nv,i),ulog)
    if (veltest(1:1).eq.'1') call finitedif_vel(kp(:,i),ndn,vel,eivl,eivc)
 enddo kloop

 end subroutine get_frequencies
!=====================================================
  subroutine finitedif_vel(qq,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
! it would give zero near band crossings, thus HF is a better way to do it
  use constants
  use lattice, only : cart2red_g,volume_g0
  use geometry
  implicit none
  integer, intent(in) :: ndn
  real(r15),intent(in) :: qq(3)
  real(r15),intent(out) :: vgr(3,ndn),evl0(ndn)
  complex(r15),intent(out) :: evc0(ndn,ndn)
  integer i,l
  real(r15) q1(3),dq,dq0,v2(3,ndn),evlp(ndn),evlm(ndn) ,shft(3),qr(3),v1(3,ndn)
  complex(r15) evct(ndn,ndn)
! a fixed direction whose components are mutually incommensurate: sqrt(2)-1,
! sqrt(3)-1, sqrt(5)-2. Shifting along it leaves the point on no symmetry plane
! or axis of any lattice, which is what lifts the band degeneracies.
  real(r15), parameter :: uinc(3)= (/ 0.41421356237309505d0,  &
  &                                   0.73205080756887729d0,  &
  &                                   0.23606797749978970d0 /)

! shift slightly off the high-symmetry directions so that degenerate bands are
! separated before the finite differences are taken. This used to be a random
! shift (call random_number), which made the frequencies, the DOS and the QHA
! properties differ from one run to the next; the acoustic branches near Gamma
! are small enough that the shift was their whole value. The direction is now
! fixed, so runs are reproducible.
 dq0=3d-4 *volume_g0**0.33333333  ! to break the degeneracies
 shft=dq0*uinc/sqrt(dot_product(uinc,uinc))
 qr=qq + shft
 dq=3d-5 *volume_g0**0.33333333

      write(*,4)' get_freq k= ',qr

 call get_freq(qr,ndn,v1,evl0,evc0) 
! eivals being non-degenerate no rediagonalization of dD/dk is needed

! return
!! CAN USE BAND CONNECTION FOR CORRECT EVALUATION OF VELOCITIES

  q1=qr
  do i=1,3
     q1(i)=qr(i)+dq
     call get_freq(q1,ndn,v2,evlp,evct)  ! v2 used here as dummy
     q1(i)=qr(i)-dq
     call get_freq(q1,ndn,v2,evlm,evct)
     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi *c_light
     q1(i)=qr(i) ! back to qr
  enddo

  write(330,4)'#la , v_FD ; v_PERT *** FOR q_red= ',cart2red_g(qr) ,' *** qcart=',qr 
  do l=1,ndn
     write(330,5)' ',l,evl0(l),vgr(:,l), length(vgr(:,l)) , v1(:,l),length(v1(:,l))
  enddo
! return the Hellmann-Feynman velocities; the finite differences above are only a check
! (written to veltest.dat), and they underestimate v near band crossings
  vgr=v1


 4 format(a,3(1x,f8.4),a,3(1x,g11.4))
 5 format(a,i5,g11.4,2(5x,f11.3,3(1x,f11.3)))
 6 format(a,3(1x,f9.3),1x,g12.5,9(1x,f9.4))

  end subroutine finitedif_vel
!===========================================================
 subroutine get_freq(kp,ndn,vg,eival,eivec)
!! given a kpoint kp, it sets up the dynamical matrix and diagonalizes it
!! eival is in eV/uma/Ang^2 ; group velocity calculated from Hellman-Feynmann theorem
 use params
 use ios
 use eigen, only: mysqrt
 use constants
 use geometry, only : is_integer,v2a,length
 use born
 use lattice, only : cart2red,r01,r02,r03
 use atoms_force_constants, only : natom_prim_cell
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: kp(3)
 real(r15), intent(out):: eival(ndn),vg(3,ndn)
 complex(r15), intent(out) :: eivec(ndn,ndn)
 integer j,l,ier,nd2,mp(ndn),al,be,tau,ia,ib 
 integer ,save :: nksave
 real(r15) absvec,eivl(ndn),kdotr,vg2(3,ndn),vtmp(ndn)  ! automatic arrays
 complex(r15) dynmat(ndn,ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn),dynm0(ndn,ndn),sumd
 complex(r15) dn(3,3),ddn(3,3,3)

 nd2 = min(ndn,12)

 nksave=nksave+1
! write(*,3)'get_freq k=',nksave,kp
!   write(*,*) ' calling set_dynamical_matrix'

! both of these should work the same since we are working on the extended grid
 if (born_flag.gt.10) then
   call set_dynamical_matrix_hard(kp,dynmat,ndn,ddyn)  ! uses phi_sr=phi_bare-phi_na (5D arrays)
 else  
   call set_dynamical_matrix_std (kp,dynmat,ndn,ddyn)  ! uses fc2 which has NA forces subtracted
 endif   

!  call set_dynamical_matrix_convolution(kp,dynmat,ndn,ddyn)
!  call step2smooth(kp,ndn,dynmat,ddyn)

   call make_hermitian(ndn,dynmat,ddyn)

   if(length(kp).lt.1d-12)  then
      call check_asr_dynmat(ndn,dynmat,sumd,1,ier)
      if(ier.ne.0) write(ulog,*)' GET_FREQ: asr broken!! ', sumd
   endif

    if (verbose) then
       write(ulog,3)' ===================================================='
       write(ulog,3)' THE DYNAMICAL MATRIX for KP# (reduced,cart)=',nksave,cart2red(kp,'g'),kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ===================================================='
       do al=1,3  ! ddyn should be antisymmetric wr swithing of l,j
          write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,al ),j=1,nd2)
          enddo
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,temp,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = eivl(mp(j))  ! so that all frequencies are positive
       eivec(:,j) = eivc(:,mp(j))
    enddo

  if (verbose) then

    do l=1,ndn
        write(ulog,3)' GET_FREQ:-----------------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo

! now calculate the square of frequencies based on dynmat
     temp = matmul(transpose(conjg(eivec)),matmul(dynmat,eivec))
     absvec = sum(cdabs(temp))
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e=',j, temp(j,j)
     enddo

!return

     do j=1,ndn
       do l=1+j,ndn-1
          if (abs(temp(l,j)).gt.1d-5*max(1d0,absvec)) then
             call warn2(ulog,'GET_FREQ: e.D.e has off-diagonal elements!!')
             write(ulog,7)'l,j,[e.D.e](l,j)=',l,j,temp(l,j)
             write(*   ,7)'l,j,[e.D.e](l,j)=',l,j,temp(l,j)
     !       stop
          endif
       enddo
     enddo

  endif

! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
! V_pert(k,l)= < e_l(k) | dD/dk | e_l(k) > / 2 om(k,l)
!   write(*,*) ' calculating group velocities'
  do al=1,3
      temp = matmul(transpose(conjg(eivec)),matmul(ddyn(:,:,al),eivec))
! eival are non-degenerate since kpoints have been randomly shifted from high symmetry positions=> no diagonalization needed
!     call diagonalize(ndn,temp,vtmp,ndn,eivc,ier)
      do l=1,ndn
         vg(al,l) = real(temp(l,l)) 
         vg(al,l) = vg(al,l)/2/mysqrt(abs(eival(l))+1d-24)*cnst*1d-10*100*2*pi*c_light
      enddo
  enddo

 if (verbose) then
    write(ulog,*)' three components of the velocity of all bands in m/s are '
    do al=1,3
      write(ulog,5)'alpha,v_alpha(lambda) (pert)=',al,vg(al,:)
    enddo
 endif


 3 format(a,i5,9(1x,f11.4))
 4 format(99(1x,2(1x,g11.4),1x))
 5 format(a,i5,99(1x,2(1x,g11.4),1x))
 6 format(a,99(1x,i5))
 7 format(a,2i6,2x,99(1x,f8.4))
 8 format(99(1x,2(1x,f9.4),1x))
 9 format(a,i5,99(1x,2(1x,f9.4),1x))

 end subroutine get_freq
!==========================================================
 subroutine get_freq_fromfit(kp,ndn,vg,eival,eivec)
!! given a kpoint kp, it sets up the dynamical matrix and diagonalizes it
!! eival is in eV/uma/Ang^2 ; group velocity calculated from Hellman-Feynmann theorem
 use params
 use ios
 use eigen, only: mysqrt
 use constants
 use geometry, only : is_integer,v2a,length
 use born
 use lattice, only : cart2red,r01,r02,r03
 use atoms_force_constants, only : natom_prim_cell,atom0
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: kp(3)
 real(r15), intent(out):: eival(ndn),vg(3,ndn)
 complex(r15), intent(out) :: eivec(ndn,ndn)
 integer j,l,ier,nd2,mp(ndn),al,be,ga,tau,taup,ia,ib ,srow,scol
 integer ,save :: nksave
 real(r15) absvec,eivl(ndn),kdotr,vg2(3,ndn),vtmp(ndn),mi,mj
 complex(r15) dynmat(ndn,ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn),dynm0(ndn,ndn),sumd
 complex(r15) dn(3,3),ddn(3,3,3),dyn4(natom_prim_cell,natom_prim_cell,3,3),ddyn4(natom_prim_cell,natom_prim_cell,3,3,3)

 nd2 = min(ndn,12)

 nksave=nksave+1
! write(*,3)'get_freq k=',nksave,kp
!   write(*,*) ' calling set_dynamical_matrix'

   call fourier_fc2fit_4D(kp,dyn4,ddyn4) 
   
 do tau=1,natom_prim_cell
    srow=3*(tau-1)+1
    mi = atom0(tau )%mass
 do taup=1,natom_prim_cell
    scol=3*(taup-1)+1
    mj = atom0(taup)%mass
 
    dynmat(srow:srow+2,scol:scol+2)=  dyn4(tau,taup,:,:) /sqrt(mi*mj)
    do ga=1,3
       ddyn  (srow:srow+2,scol:scol+2,ga)=  ddyn4(tau,taup,:,:,ga) /sqrt(mi*mj)
    enddo
 enddo
 enddo

   call make_hermitian(ndn,dynmat,ddyn)

   if(length(kp).lt.1d-12)  then
      call check_asr_dynmat(ndn,dynmat,sumd,1,ier)
      if(ier.ne.0) write(ulog,*)' GET_FREQ_fromfit: asr broken!! ', sumd
   endif

    if (verbose) then
       write(ulog,3)' ===================================================='
       write(ulog,3)' THE DYNAMICAL MATRIX for KP# (reduced,cart)=',nksave,cart2red(kp,'g'),kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ===================================================='
       do al=1,3  ! ddyn should be antisymmetric wr swithing of l,j
          write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,al ),j=1,nd2)
          enddo
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,temp,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = eivl(mp(j))  ! so that all frequencies are positive
       eivec(:,j) = eivc(:,mp(j))
    enddo

  if (verbose) then

    do l=1,ndn
        write(ulog,3)' GET_FREQ_fromfit:-----------------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo

! now calculate the square of frequencies based on dynmat
     temp = matmul(transpose(conjg(eivec)),matmul(dynmat,eivec))
     absvec = sum(cdabs(temp))
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e=',j, temp(j,j)
     enddo

!return

     do j=1,ndn
       do l=1+j,ndn-1
          if (abs(temp(l,j)).gt.1d-5*max(1d0,absvec)) then
             call warn2(ulog,'GET_FREQ_fromfit: e.D.e has off-diagonal elements!!')
             write(ulog,7)'l,j,[e.D.e](l,j)=',l,j,temp(l,j)
             write(*   ,7)'l,j,[e.D.e](l,j)=',l,j,temp(l,j)
     !       stop
          endif
       enddo
     enddo

  endif

! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
! V_pert(k,l)= < e_l(k) | dD/dk | e_l(k) > / 2 om(k,l)
!   write(*,*) ' calculating group velocities'
  do al=1,3
!     temp = matmul(transpose(conjg(eivec)),matmul(ddyn(:,:,al),eivec))
!     call diagonalize(ndn,temp,vtmp,ndn,eivc,ier)
! assume ddyn is non-degenerate since kpoints have been randomly shifted from high symmetry positions
      do l=1,ndn
         vg(al,l) = real( dot_product(conjg(eivec(:,l)),matmul(ddyn(:,:,al),eivec(:,l))) )
         vg(al,l) = vg(al,l)/2/sqrt(abs(eival(l))+1d-14)*cnst*1d-10*100*2*pi*c_light
!        vg(al,l) =  vtmp(l)/2/mysqrt(abs(eival(l))+1d-24)*cnst*1d-10*100*2*pi*c_light
      enddo
  enddo

 if (verbose) then
    write(ulog,*)' three components of the velocity of all bands in m/s are '
    do al=1,3
      write(ulog,5)'alpha,v_alpha(lambda) (pert)=',al,vg(al,:)
    enddo
 endif


 3 format(a,i5,9(1x,f11.4))
 4 format(99(1x,2(1x,g11.4),1x))
 5 format(a,i5,99(1x,2(1x,g11.4),1x))
 6 format(a,99(1x,i5))
 7 format(a,2i6,2x,99(1x,f8.4))
 8 format(99(1x,2(1x,f9.4),1x))
 9 format(a,i5,99(1x,2(1x,f9.4),1x))

 end subroutine get_freq_fromfit
!==========================================================
 subroutine set_dynamical_matrix_hard(kpt,dynmat,ndim,ddyn)
!! uses the hard or step phase convention and adds the short-range part of FCs to the NA part
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use fourier
 use ewald
 use born, only : phi_sr,born_flag
 use constants, only : r15,ci
 use ios, only : ulog,write_out
 implicit none
 integer, intent(in) :: ndim
 real(r15), intent(in) :: kpt(3)
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3) 
 complex(r15) dn(3,3),ddn(3,3,3),sr(3,3),dsr(3,3,3),zz,sumd
 complex(r15) , allocatable :: auxr(:)
 integer tau,taup,al,be,ga,ir,srow,scol
 real(r15) mi,mj,rr(3),srmax,anmax,asrs(3,3),asrc(3,3)

 ! Safety checks
 if (.not. allocated(phi_sr)) then
    write(*,*) 'ERROR in set_dynamical_matrix_hard: phi_sr not allocated!'
    stop
 endif
 if (.not. allocated(rgrid_xtnd)) then
    write(*,*) 'ERROR in set_dynamical_matrix_hard: rgrid_xtnd not allocated!'
    stop
 endif
 if (.not. allocated(rws_weights_xtnd)) then
    write(*,*) 'ERROR in set_dynamical_matrix_hard: rws_weights_xtnd not allocated!'
    stop
 endif

 ! Allocate auxr with current grid size
 if (allocated(auxr)) deallocate(auxr)
 allocate(auxr(nrgrid_xtnd))

 if( born_flag.gt.0) then
    call set_dynamical_NA_3N3N(kpt,dynmat,ndim,ddyn)
 else 
    dynmat =0; ddyn=0
 endif
 
! write(*,*) 'Debug: nrgrid_xtnd=', nrgrid_xtnd, ' size(phi_sr,5)=', size(phi_sr,5)
  write(*,'(a,999(f6.2))') 'kp_set_dyn_mat_hard=', kpt,cart2red(kpt,'g')

 if(length(kpt).lt.1d-12)  then
    call check_asr_dynmat(ndim,dynmat,sumd,1,ir)  ! 1 because it has mass denominator
    if(ir.ne.0) write(ulog,*)' SETDYN_HARD: asr broken in dyn_na!! ', sumd
 endif

 srmax=0;anmax=0
 do tau=1,natom_prim_cell
    srow=3*(tau-1)+1
    mi = atom0(tau )%mass
    if(length(kpt).lt.1d-6) asrs=0;asrc=0
 do taup=1,natom_prim_cell
    scol=3*(taup-1)+1
    mj = atom0(taup)%mass
 
! the Short-Range part is contained in sr
    sr=cmplx(0.0_r15,0.0_r15); dsr=cmplx(0.0_r15,0.0_r15)
    do al=1,3
    do be=1,3
       ! Bounds check
       if (size(phi_sr,1) < tau .or. size(phi_sr,2) < taup) then
          write(*,*) 'ERROR: phi_sr indices out of bounds!'
          write(*,*) 'tau=', tau, 'taup=', taup
          write(*,*) 'size(phi_sr,1)=', size(phi_sr,1), 'size(phi_sr,2)=', size(phi_sr,2)
          stop
       endif
       if (size(phi_sr,5) /= nrgrid_xtnd) then
          write(*,*) 'ERROR: phi_sr size mismatch!'
          write(*,*) 'size(phi_sr,5)=', size(phi_sr,5), 'nrgrid_xtnd=', nrgrid_xtnd
          stop
       endif
       
       call fourier_r2k_c_single(kpt,phi_sr(tau,taup,al,be,:),nrgrid_xtnd,rgrid_xtnd,rws_weights_xtnd,sr(al,be))
       do ga=1,3
          auxr = phi_sr(tau,taup,al,be,:) * ci*rgrid_xtnd(ga,:)
          call fourier_r2k_c_single(kpt,auxr,nrgrid_xtnd,rgrid_xtnd,rws_weights_xtnd,dsr(al,be,ga))
       enddo
    enddo
    enddo

    if(length(kpt).lt.1d-6) asrs = asrs + sr

    anmax = anmax + sum(abs(dynmat))/natom_prim_cell**2

    dynmat(srow:srow+2,scol:scol+2)=  dynmat(srow:srow+2,scol:scol+2)  +  sr  /sqrt(mi*mj)
    ddyn  (srow:srow+2,scol:scol+2,:)=  ddyn(srow:srow+2,scol:scol+2,:)+ dsr  /sqrt(mi*mj)
 
    srmax = srmax + sum(abs(sr))

 enddo
    if(length(kpt).lt.1d-6) write(ulog,6)'SDYNHARD tau,ASR_Coulomb=',tau,asrc
    if(length(kpt).lt.1d-6) write(ulog,6)'SDYNHARD tau,ASR_Short  =',tau,asrs
 enddo

 if (abs(srmax).gt.1d-6) write(*,7) '==== kp,SR,NA,NA/SR:',kpt,srmax,anmax,anmax/(srmax+1d-12)

 if (allocated(auxr)) deallocate(auxr)

6 format(a,i4,99(1x,f9.5))
7 format(a,99(1x,g11.4))

 end subroutine set_dynamical_matrix_hard
!==========================================================
 subroutine set_dynamical_matrix_smooth(kpt,dynmat,ndim,ddyn)
! calculates dynmat with fitted FCs and using the smooth phase
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15,ci
 use ios, only : ulog
 use born, only: dyn_naq0,dyn_na,born_flag
 use ewald, only : eta,ng_ewald,g_ewald ,nr_ewald,r_ewald ,ewald_2nd_deriv_Gnew,ewald_2nd_deriv_hat
 implicit none
 integer, intent(in) :: ndim
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3)
 real(r15), intent(in) :: kpt(3)
 complex(r15) junk
 real(r15) mi,mj,rr(3),rfold(3),dta(3)
 integer tau,j,taup,al,be,ga,i3,j3,t,cnt,ti,ired,g,nb,srow,scol
 logical insid
 complex(r15) phase ,d2ew(3,3),d3ew(3,3,3)

4 format(a,4(1x,i5),9(2x,f9.4))
8 format(a,99(1x,f9.3))
9 format(i6,99(1x,f9.3))

 ddyn   = dcmplx(0d0,0d0)
 dynmat = dcmplx(0d0,0d0)
 do tau = 1,natom_prim_cell
    srow=3*(tau-1)+1
    mi = atom0(tau)%mass
    do al=1,3
       i3=3*(tau-1)+al
! write(ulog,*) 'i,al,mass=',tau,al,i3,mi
    gloop: do g=1,map(2)%ngr
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          if( map(2)%keep(counteri(2,g,ti) ) .eq. 0 ) cycle 
          ired = map(1)%nkeptind + current(2,g,ti)
          if( sum(map(2)%keep(1:counteri(2,g,ti) )) .gt. map(2)%nkeptind) then
             write(*,*)'SET_DYNMAT: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
             be = map(2)%gr(g)%ixyz(2,t)
             j  = map(2)%gr(g)%iat(2,t)
             taup = iatomcell0(j)  
             mj = atom0(taup)%mass
             j3 = be+3*(taup-1)
! k1 6-22-24  this is with the smooth phase
             rr = atompos(:,j) - atompos(:,tau)   ! this is R+taup-tau  ! new phase convention
             call check_inside_ws(rr,rws26,insid,nb) 
             if(insid) then
                rfold = rr
             else
                rfold = fold_ws(rr,rws26) 
                write(ulog,8)'#*#*#*#*#*#* SMOOTH : rr outside WS, need to fold:',rr,rfold
                write(*   ,8)'#*#*#*#*#*#* SMOOTH : rr outside WS, need to fold:',rr,rfold
             endif

! k1 6-22-24  this is ith the smooth phase
             junk = fcs(ired) * map(2)%gr(g)%mat(t,ti) * exp(ci*(kpt .dot. rfold))/sqrt(mi*mj)  ! &
             dynmat(i3,j3) = dynmat(i3,j3) + junk
             ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rfold(:)
          enddo tloop
       enddo tiloop
    enddo gloop
    enddo

    if(born_flag.ne.0) then

       do taup=1,natom_prim_cell
          scol=3*(taup-1)+1
          mj = atom0(taup)%mass
          call ewald_2nd_deriv_hat (kpt,tau,taup,nr_ewald,ng_ewald,r_ewald,g_ewald,eta,d2ew,d3ew)
          d2ew = matmul(matmul(atom0(tau)%charge,d2ew),transpose(atom0(taup)%charge)) 
          do ga=1,3
             d3ew(:,:,ga)=matmul(matmul(atom0(tau)%charge,d3ew(:,:,ga)),transpose(atom0(taup)%charge)) 
          enddo
          if(tau.eq.taup) d2ew=d2ew-dyn_naq0(tau,:,:)

          dta= atompos(:,taup)-atompos(:,tau) 
          phase= exp(ci*dot_product(kpt,dta)) /sqrt(mi*mj)  
          dynmat(srow:srow+2,scol:scol+2) = dynmat(srow:srow+2,scol:scol+2) + d2ew(:,:)*phase
          do ga = 1,3
          ddyn(srow:srow+2,scol:scol+2,ga) = ddyn(srow:srow+2,scol:scol+2,ga) + (d3ew(:,:,ga)+ci*dta(ga)*d2ew(:,:))*phase
          enddo
       enddo

    endif

 enddo


3 format(a,i5,9(f8.3))
5 format(a,5i5,9(f8.3))

 end subroutine set_dynamical_matrix_smooth
!==========================================================
 subroutine fourier_fc2fit_4D(kpt,dynmat,ddyn)
!! calculates from the fitted FC2s, phi_bare, the dynamical matrix without NA, mass term 
!! or soft phase; it is just the fourier transform of the FC2s at kpt; it is G0-periodic
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15,ci
 use ios, only : ulog
 use born, only : phi_bare
 use fourier, only : rgrid_xtnd,nrgrid_xtnd,nrgrid,rgrid,rws_weights
 implicit none
 complex(r15), intent(out) :: dynmat(natom_prim_cell,natom_prim_cell,3,3),ddyn(natom_prim_cell,natom_prim_cell,3,3,3) 
 real(r15), intent(in) :: kpt(3)
 complex(r15) junk
 real(r15) mi,mj,rr(3),rfold(3),rs(3) !,dr(3)
 integer tau,j,taup,al,be,i3,j3,t,ti,ired,g,nb,igrid,ir
 logical insid,found_in_grid

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 ddyn   = dcmplx(0d0,0d0)
 dynmat = dcmplx(0d0,0d0); phi_bare =0
 do tau = 1,natom_prim_cell
 do al  = 1,3
    i3 = al+3*(tau-1)
    gloop: do g=1,map(2)%ngr
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          if( map(2)%keep(counteri(2,g,ti) ) .eq. 0 ) cycle 
          ired = map(1)%nkeptind + current(2,g,ti)
          if( sum(map(2)%keep(1:counteri(2,g,ti) )) .gt. map(2)%nkeptind) then
             write(*,*)'SET_DYNMAT4D: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
             be = map(2)%gr(g)%ixyz(2,t)
             j  = map(2)%gr(g)%iat(2,t)
             taup = iatomcell0(j)    
             j3 = be+3*(taup-1)
             rr = atompos(:,j) - atompos(:,taup) ! this is just R  used for standard phase; makes a G-periodic dynmat
             rs = atompos(:,j) - atompos(:,tau) ! this is R+taup-tau  used for soft phase

! a member of a kept group may lie outside the supercell's WS cell when the supercell lacks some of
! the crystal's rotations; it is summed at its own R (counted once in the log, setup_FC2_in_supercell)
!             call check_inside_ws(rr,rws26,insid,nb) 
!             if(.not.insid) then 
!                write(ulog,4) 'WARNING: R not found in grid:t,tp,j,ired,R,R_red=', &
! &                  tau, taup, j, ired, rr, cart2red(rr,'r'),length(rr)
!                write(ulog,*) '  FC value would be:', fcs(ired) * map(2)%gr(g)%mat(t,ti)
!                write(*,*)' fourier_fc2fit_4D: R not inside!',cart2red(rr,'r')
!       !        stop
!             endif

! Add this check (only do it once, not in every loop):
!  found_in_grid = .false.
!  do ir=1,nrgrid
!     if (length(rr - rgrid(:,ir)) < 1d-6) then
!        found_in_grid = .true.
!        exit
!     endif
!  enddo
   
!  if (.not. found_in_grid) then
!     write(*   ,*) 'FC pair R NOT in rgrid!!:', tau,taup,j,rr
!     write(ulog,*) 'FC pair R NOT in rgrid!!:', tau,taup,j,rr
!     stop
!  endif


             junk = fcs(ired) * map(2)%gr(g)%mat(t,ti) * exp(ci*(kpt .dot. rr))  ! makes dynmat G0-periodic
             dynmat(tau,taup,al,be) = dynmat(tau,taup,al,be) + junk
             ddyn(tau,taup,al,be,:) = ddyn(tau,taup,al,be,:) + junk*ci*rr(:)

          enddo tloop

       enddo tiloop
    enddo gloop
 enddo
 enddo

3 format(a,2i4,9(1x,f9.4))
6 format(9(1x,f9.4))

 end subroutine fourier_fc2fit_4D
!==========================================================
 subroutine get_phi_bare_5D(nat,nr,grd,phi) 
!! calculates from the fitted FC2s, the 5-dimensional phi_bare(tau,taup,al,be,rgrid) on a grid defined by grd 
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15,ci
 use ios, only : ulog
 use born, only : phi_bare
 use fourier, only : rgrid_xtnd,nrgrid_xtnd
 implicit none
 integer, intent(in) :: nr,nat
 real(r15), intent(in) :: grd(3,nr) 
 complex(r15), intent(out) :: phi(nat,nat,3,3,nr) 
 real(r15) rr(3)
 integer tau,i0,j,taup,al,be,g,t,ti,ired,ir,nfound,nmissing,nb
 logical found,insid

   if(nat.ne.natom_prim_cell ) then ! .or. nr.ne.nrgrid_xtnd) then
      write(*,*)' nat.ne.natom_prim_cell; check your arguments:',nat,natom_prim_cell
 !    write(*,*)' nr .ne.nrgrid_xtnd ; check your arguments:',nr,nrgrid_xtnd
      stop
    endif

    nfound=0; nmissing=0
   
    phi = dcmplx(0d0,0d0);
 do tau = 1,natom_prim_cell
 do al  = 1,3
    gloop: do g=1,map(2)%ngr
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          if( map(2)%keep(counteri(2,g,ti) ) .eq. 0 ) cycle 
          ired = map(1)%nkeptind + current(2,g,ti)
          if( sum(map(2)%keep(1:counteri(2,g,ti) )) .gt. map(2)%nkeptind) then
             write(*,*)'get_PHI_BARE_5D: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
             be  = map(2)%gr(g)%ixyz(2,t) 
             j   = map(2)%gr(g)%iat(2,t)
             taup= iatomcell0(j)    
!        write(*,*)'get_phi_bare_5d: t,tp,al,be,j=',tau,taup,al,be,j 
             rr = atompos(:,j) - atompos(:,taup) ! this is just R  

      !      call check_inside_ws(rr,rws26,insid,nb) 
      !      if(.not.insid) then 
      !         write(ulog,4) 'WARNING: R not inside WS :t,tp,j,ired,R,R_red=', &
 &    !             tau, taup, j, ired, rr, cart2red(rr,'r'),length(rr)
      !         write(ulog,*) '  FC value would be:', fcs(ired) * map(2)%gr(g)%mat(t,ti)
      !         write(*   ,4) 'WARNING: R not inside WS :t,tp,j,ired,R,R_red=', &
 &    !             tau, taup, j, ired, rr, cart2red(rr,'r'),length(rr)
      !         write(*   ,*) '  FC value would be:', fcs(ired) * map(2)%gr(g)%mat(t,ti)
      !      endif
             found=.false.
             do ir=1,nr
                if (length(rr - grd(:,ir)) < 1d-3) then
                   phi(tau,taup,al,be,ir) = phi(tau,taup,al,be,ir) + &
 &                            fcs(ired) * map(2)%gr(g)%mat(t,ti)
                   found=.true.
                   nfound=nfound+1
                   exit
                endif
             enddo

             if (.not. found) then
               nmissing = nmissing + 1
               write(ulog,4) 'WARNING: R not found in grid:t,tp,j,ired,R,R_red=', &
 &                tau, taup, j, ired, rr, cart2red(rr,'r'),length(rr)
               write(ulog,*) '  FC value would be:', fcs(ired) * map(2)%gr(g)%mat(t,ti)
               write(*,4) 'WARNING: R not found in grid:t,tp,j,ired,R,R_red=', &
 &                tau, taup, j, ired, rr, cart2red(rr,'r'),length(rr)
               write(*,*) 'Check the grid; it seems to have been exceeded!! '
    !          stop
             endif
   
          enddo tloop
       enddo tiloop
    enddo gloop
 enddo
 enddo

write(*,*) 'FC terms found in grid:', nfound
write(*,*) 'FC terms MISSING from grid:', nmissing

write(ulog,*) 'FC terms found in grid:', nfound
write(ulog,*) 'FC terms MISSING from grid:', nmissing

3 format(a,2i4,9(1x,f9.4))
4 format(a,4i4,9(1x,f9.4))
6 format(9(1x,f9.4))

 end subroutine get_phi_bare_5D
!==========================================================
 subroutine set_dynamical_matrix_convolution(kp,dynmat,ndn,ddyn)
!! calculates the dynamical matrix at arbitrary q from its knowledge on a grid associated 
!! with the supercell by doing a convolution: D(q)=sum_G* D(G*) S(q-G*) w(G*) 
!! where w is the weight of the wavevector G* and S is the structure factor: S(G*)=0; S(0)=1=S(G)
 use geometry, only : v2a,length
 use atoms_force_constants, only : natom_prim_cell,atom0
 use born, only : phi_sr,born_flag,dyn_g,dyn_na
 use fourier, only : nggrid,ggrid ,gws_weights
 use constants, only : r15,ci
 implicit none
 real(r15), intent(in):: kp(3)
 integer, intent(in) :: ndn
 complex(r15), intent(out) :: dynmat(ndn,ndn),ddyn(ndn,ndn,3)
 integer al,be,tau,taup,g,ia,ib
 complex(r15) dn(3,3),ddn(3,3,3),sr(3,3),dsr(3,3,3),zz,aux(nggrid)
 real(r15) mi,mj

 do tau=1,ndn/3
    mi = atom0(tau )%mass
 do taup=1,ndn/3
    mj = atom0(taup)%mass
!   call dyn_coulomb_structfactor(kp,tau,taup,dn,ddn)
 do al=1,3
 do be=1,3
    ia=al+3*(tau-1)
    ib=be+3*(taup-1)

!   aux = (dyn_g(tau,taup,al,be,:) +  dyn_na(tau,taup,al,be,:)) ! not divided by mass
    aux =  dyn_na(tau,taup,al,be,:)   ! must be an array defined on ggrid
    call d_sf_conv(kp,aux,dn(al,be),ddn(al,be,:),nggrid,ggrid,gws_weights)

! the SR part is contained in sr
!   sr=0; dsr=0
!   do igrid=1,nrgrid
!      rr=rgrid(:,igrid)
!      zz = exp(ci*dot_product(kpt,rr))
!      sr = sr + phi_sr(tau,taup,:,:,igrid) * zz
!      do ga=1,3
!         dsr(:,:,ga)=dsr(:,:,ga) + phi_sr(tau,taup,:,:,igrid) * zz*ci*rr(ga)
!      enddo
!   enddo

    dynmat(ia,ib)=dn(al,be) / sqrt(mi*mj)
    ddyn(ia,ib,:)=ddn(al,be,:) / sqrt(mi*mj)
 enddo
 enddo
 enddo
 enddo

 end subroutine set_dynamical_matrix_convolution
!==========================================================
 subroutine d_sf_conv(kpt,dystar,dyn,ddyn,ng,grd,weig)
! calculates dynmat using the convoulution of dyn(G*) with structure factor(q-G*)
 use params
 use lattice
 use fourier
 implicit none

 integer, intent(in) :: ng
 real(r15), intent(in) :: kpt(3),grd(3,ng),weig(ng)
 complex(r15), intent(in) :: dystar(ng) 
 complex(r15), intent(out) ::   dyn ,ddyn(3)
 integer g
 real(r15) sf,dsf(3)
! complex(r15) sf,dsf(3)

! write(*,3)'d_sf_conv k,k_red=',ng,kpt,cart2red(kpt,'g')
 dyn=0;ddyn=0
 do g=1,ng
    call structure_factor_recip(kpt-grd(:,g),nrgrid,rgrid,rws_weights,sf,dsf)
    dyn= dyn+dystar(g)*sf*weig(g)
!   write(*,3)'g,dyn_na(g),sf,w(g)=',g, dystar(g),sf,weig(g)
    ddyn=ddyn+dystar(g)*dsf*weig(g)
 enddo
 
 3 format(a,i5,9(1x,f11.4))

 end subroutine d_sf_conv
!==========================================================
 subroutine set_dynamical_matrix_conv_arr(kpt,dystar,dyn,ng,grd)
! calculates dynmat with fitted FCs and using the standard phase
 use params
 use lattice
 use fourier
 use atoms_force_constants, only : natom_prim_cell
 implicit none

 integer, intent(in) :: ng
 real(r15), intent(in) :: kpt(3),grd(3,ng)
 complex(r15), intent(in) :: dystar(natom_prim_cell,natom_prim_cell,3,3,ng) 
 complex(r15), intent(out) ::   dyn(natom_prim_cell,natom_prim_cell,3,3,ng) 

 complex(r15) junk
 real(r15) mi,mj,rr(3),rfold(3),dta(3)
 integer g
 complex(r15) st,dst(3)

 dyn=0
 do g=1,ng
    call structure_factor_c(kpt-grd(:,g),nrgrid,rgrid,st,dst)
    dyn(:,:,:,:,g)=dyn(:,:,:,:,g)+dystar(:,:,:,:,g)*st
 enddo
 
 end subroutine set_dynamical_matrix_conv_arr
!==========================================================
 subroutine set_dynamical_matrix_std(kpt,dynmat,ndim,ddyn)
! calculates dynmat with fitted FCs and using the standard phase
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15,ci
 use ios, only : ulog
 use born, only: dyn_naq0,dyn_na,born_flag
 use ewald, only : eta,ng_ewald,g_ewald ,nr_ewald,r_ewald ,ewald_2nd_deriv_Gnew,ewald_2nd_deriv_hat
 use fourier, only: nrgrid,rgrid,rws_weights,nggrid,ggrid,gws_weights
 implicit none
 integer, intent(in) :: ndim
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3)
 real(r15), intent(in) :: kpt(3)
 complex(r15) junk
 real(r15) mi,mj,rr(3),rfold(3),dta(3),rs(3)
 integer tau,j,taup,al,be,i3,j3,t,cnt,ti,ired,g,nb,ir
 logical insid,found_in_grid 
 complex(r15) dna,ddna(3),auxg(nggrid) !phase,,dn(3,3),ddn(3,3,3)

4 format(a,4(1x,i5),9(2x,f9.4))
8 format(a,99(1x,f9.3))
9 format(i6,99(1x,f9.3))

 ddyn   = dcmplx(0d0,0d0)
 dynmat = dcmplx(0d0,0d0)
 call set_dynamical_NA_3N3N(kpt,dynmat,ndim,ddyn)  ! with mass term

 do tau = 1,natom_prim_cell
    mi = atom0(tau)%mass
 do al  = 1,3
    i3 = al+3*(tau-1)
! write(ulog,*) 'i,al,mass=',tau,al,i3,mi
    gloop: do g=1,map(2)%ngr
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          if( map(2)%keep(counteri(2,g,ti) ) .eq. 0 ) cycle 
          ired = map(1)%nkeptind + current(2,g,ti)
          if( sum(map(2)%keep(1:counteri(2,g,ti) )) .gt. map(2)%nkeptind) then
             write(*,*)'SET_DYNMAT_STD: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
             be = map(2)%gr(g)%ixyz(2,t)
             j  = map(2)%gr(g)%iat(2,t)
             taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
             mj = atom0(taup)%mass
             j3 = be+3*(taup-1)
             rr = atompos(:,j) - atompos(:,taup) ! this is just R 
             rs = atompos(:,j) - atompos(:,tau ) ! this is R+taup-tau 

! a member of a kept group may lie outside the supercell's WS cell when the supercell lacks some of
! the crystal's rotations; it is summed at its own R (counted once in the log, setup_FC2_in_supercell)
         !   call check_inside_ws(rr,rws26,insid,nb) 
! by constriction it should be inside 
         !   if(.not.insid) then 
         !      write(*,3)' SET_DYNMAT_STD: R not inside!,t,tp,j,R',tau,taup,j,cart2red(rr,'r')
!        !      stop
         !   endif

! Add this check (only do it once, not in every loop):
!            found_in_grid = .false.
!            do ir=1,nrgrid
!               if (length(rr - rgrid(:,ir)) < 1d-6) then
!                 found_in_grid = .true.
!                 exit
!               endif
!            enddo
   
             junk = fcs(ired) * map(2)%gr(g)%mat(t,ti) * exp(ci*(kpt .dot. rr))/sqrt(mi*mj) 
             dynmat(i3,j3) = dynmat(i3,j3) + junk !+ dna/sqrt(mi*mj) 
             ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:) !+ ddna/sqrt(mi*mj) 

          enddo tloop
       enddo tiloop
    enddo gloop
 enddo   ! al loop
 enddo   ! tau loop


3 format(a,3i5,9(f8.3))
5 format(a,5i5,9(f8.3))

 end subroutine set_dynamical_matrix_std
!===========================================================
 subroutine set_dynamical_matrix_new(kpt,dynmat,ndim,ddyn)
!! uses fourier interpolation from phi_periodic and then adds back the NA correction, followed by smooth phase at the end
 use params
 use lattice
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15,ci
 use ios, only : ulog,write_out
 use born, only : dyn_naq0,born_flag,phi_periodic
 use fourier, only : nrgrid,rgrid,rws_weights
 use ewald, only : eta,ng_ewald,g_ewald ,nr_ewald,r_ewald ,ewald_2nd_deriv_Gnew,ewald_2nd_deriv_hat
 implicit none
 integer, intent(in) :: ndim
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3)
 real(r15), intent(in) :: kpt(3)
 real(r15) mi,mj,rr(3),dta(3),qq(3,1),error,rg(3),rfold(3)
 integer tau,i,taup,al,be,ga,i3,j3,nb
 complex(r15) phase ,d2ew(3,3),d3ew(3,3,3),dynp(1,1,3,3,1)
 logical insid

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 ddyn   = dcmplx(0d0,0d0)
 dynmat = dcmplx(0d0,0d0)
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
    mi = atom0(tau )%mass
    mj = atom0(taup)%mass
    dta= atompos(:,taup)-atompos(:,tau)
    phase= exp(ci*dot_product(kpt,dta))  ! this is used for correct group velocities

    if(born_flag.ne.0) then
   !   call ewald_2nd_deriv_Gnew(kpt,tau,taup,ng_ewald,g_ewald,eta,d2ew,d3ew)
       call ewald_2nd_deriv_hat (kpt,tau,taup,nr_ewald,ng_ewald,r_ewald,g_ewald,eta,d2ew,d3ew)
       d2ew = matmul(matmul(atom0(tau)%charge,d2ew),transpose(atom0(taup)%charge)) 
       do ga=1,3
          d3ew(:,:,ga)=matmul(matmul(atom0(tau)%charge,d3ew(:,:,ga)),transpose(atom0(taup)%charge)) 
       enddo

       if(tau.eq.taup) d2ew=d2ew-dyn_naq0(tau,:,:)

    else
       d2ew=0; d3ew=0
    endif

!   call write_out(6,'d2ew ',d2ew)

    do al=1,3
    do be=1,3
       i3 = al+3*(tau -1)
       j3 = be+3*(taup-1)

       do i=1,nrgrid
          rr=rgrid(:,i)+dta 
          call check_inside_ws(rr,rws26,insid,nb) 
          if(.not. insid) then
             rfold = fold_ws(rr,rws26) 
          else
             rfold=rr
          endif
! k1 6-22-24  this is ith the smooth phase
  
          dynmat(i3,j3) = dynmat(i3,j3)+phi_periodic(tau,taup,al,be,i) &   ! phi_periodic has the proper weights at grid i
&         *exp(ci*dot_product(kpt,rfold))/sqrt(mi*mj)  
!&         *phase*exp(ci*dot_product(kpt,rgrid(:,i)))/sqrt(mi*mj)  

          ddyn(i3,j3,:) = ddyn(i3,j3,:)+phi_periodic(tau,taup,al,be,i) & 
&         *exp(ci*dot_product(kpt,rfold))/sqrt(mi*mj) * ci*rfold(:)
!&         *phase*exp(ci*dot_product(kpt,rgrid(:,i)))/sqrt(mi*mj) * ci*(rgrid(:,i)+dta)
       enddo

! add NA correction from ewald sums
       dynmat(i3,j3) = dynmat(i3,j3) + d2ew(al,be)*phase/sqrt(mi*mj)
       ddyn(i3,j3,:) = ddyn(i3,j3,:) + (d3ew(al,be,:)+ci*dta*d2ew(al,be)) &
&                         *phase/sqrt(mi*mj)

  !    if(tau.eq.taup .and. born_flag.ne.0) then
  !       dynmat(i3,j3)=dynmat(i3,j3)-dyn_naq0(tau,al,be)*phase/sqrt(mi*mj)
  !    endif

    enddo
    enddo
 enddo
 enddo

3 format(a,i5,9(f8.3))
5 format(a,5i5,9(f8.3))

 end subroutine set_dynamical_matrix_new
!===========================================================
 subroutine make_hermitian(n,dyn,ddyn)
!! measures how far dyn is from hermitian, warns if max|D-D^+| > tol*max|D|, and
!! replaces dyn by (D+D^+)/2 (loop over the upper triangle). A large deviation points to an error elsewhere (e.g. in the
!! phase factors), not to round-off. ddyn is not touched.
 use constants
 use ios, only : ulog
 implicit none
 integer, intent(in):: n
 complex(r15), intent(inout) :: dyn(n,n),ddyn(n,n,3)
 real(r15), parameter :: tol=1d-6  ! relative; about the accuracy of the atomic coordinates
 real(r15) nrm1,dev,devmax,mi,mj
 complex(r15) dt,djt
 integer t,j,tmax,jmax
 integer, save :: nwarn=0

 nrm1 = max(maxval(cdabs(dyn)),tiny(1d0))
 devmax=0; tmax=1; jmax=1; dt=0; djt=0
! upper triangle (diagonal included): measure the deviation, then set D(t,j) and D(j,t) from
! the average of the real parts and half the difference of the imaginary parts
 do t=1,n
 do j=t,n
    dev=cdabs(dyn(t,j)-conjg(dyn(j,t)))
    if(dev.gt.devmax) then
       devmax=dev; tmax=t; jmax=j; dt=dyn(t,j); djt=dyn(j,t)
    endif
    mi=(aimag(dyn(t,j))-aimag(dyn(j,t)))/2
    mj=(real (dyn(t,j))+real (dyn(j,t)))/2
    dyn(t,j)=dcmplx(mj, mi)
    dyn(j,t)=dcmplx(mj,-mi)
 enddo
 enddo
 if (devmax .gt. tol*nrm1) then
    nwarn=nwarn+1
    write(ulog,'(a,es10.3,a,2i4,a,2(2es12.4,2x))')' WARNING MAKE_HERMITIAN: max|D-D^+|/max|D|=', &
 &        devmax/nrm1,' at i,j=',tmax,jmax,' D(i,j),D(j,i)=',dt,djt
    if(nwarn.eq.1) write(*,*)'WARNING MAKE_HERMITIAN: dynamical matrix not hermitian; see log'
 endif

 end subroutine make_hermitian
!===========================================================
 subroutine dyn_coulomb_structfactor(q,tau,taup,dn,ddn)
!! calculates the q component of the Coulomb dynamical matrix to add if born_flag=7 
! dyn(tau,taup) = [z(tau)^al,ga.q^ga][z(taup)^be.de.q^de] * sf(q)/eps0/(q.epsr.q)/volume_r0
! DOES NOT INCLUDE THE MASS DENOMINATOR, and uses the step phase
! use atoms_force_constants
 use geometry, only : v2a
 use lattice, only :  volume_r0 ,g0ws26 !, volume_r, cart2red_g ,r0ws26
 use fourier, only : nrgrid,rgrid,rws_weights ,ggrid,nggrid,gws_weights  !,nreg,rgridreg
 use constants , only : pi,r15,ci
 use ios, only : ulog,write_out
 use atoms_force_constants, only : atom0
! use ewald, only : naq
 implicit none
 integer, intent(in) :: tau,taup
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: dn(3,3),ddn(3,3,3) 
 integer ga,ig
 real(r15) sf,dsf(3),dta(3),damp
 complex(r15) d3(3,3) ,dd3(3,3,3)

 damp=10000000 ! make the gaussian=1

 do ig=1,nggrid
    call naq(ggrid(:,ig),tau,taup,damp,dn,ddn) ! this is in step phase
    call structure_factor_recip(q-ggrid(:,ig),nrgrid,rgrid,rws_weights,sf,dsf)
    dn = dn*sf*gws_weights(ig)
    do ga=1,3
       ddn(:,:,ga) = ddn(:,:,ga)*dsf(ga)*gws_weights(ig)
    enddo
 enddo

 end subroutine dyn_coulomb_structfactor
!============================================================
 subroutine gruneisen(nkp,kp,dk,ndn,eivl,eivc,grn,ugr2,ipath)
!! mode Gruneisen parameters at the k-points kp (eivl, eivc from get_frequencies with the same ipath),
!! computed by the routine shared with THERMACOND (common/gruneisen.f90, module gruneisen_rules):
!! gamma = -<e|W(k)|e>/(6 omega^2), W from the cubic FCs times the positions y = x + xi (internal
!! relaxation per unit isotropic strain), with the phase exp(ik.R) of set_dynamical_matrix_std; the 3 rigid
!! translations at Gamma get the long-wave limit along the direction of approach.
!! Until 2026-11-01 W used exp(iq.(R+taup-tau)) while the eigenvectors come from exp(ik.R): gamma was
!! wrong away from Gamma and for the acoustic branches near it (exact-model test, HfO2).
 use ios
 use lattice, only : volume_g0
 use geometry
 use atoms_force_constants
 use svd_stuff
 use eigen
 use mech, only : xi
 use constants
 use gruneisen_rules, only : dynmat_routine, w_terms_set, internal_strain_trace, mode_gruneisen
 implicit none
 integer, intent(in) :: nkp,ndn,ugr2,ipath
 real(r15), intent(in) :: kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(r15), intent(in) :: eivc(ndn,ndn,nkp)
 complex(r15), intent(out) :: grn(ndn,nkp)
 procedure(dynmat_routine) :: get_dynmat
 integer ik,tau,la,al,be,ga,j,k,taup,k0,t,g,ti,ired,n2,n3,pass,nat
 real(r15) mi,mj,y,xir(ndn),mass(natom_prim_cell),dev
 integer, allocatable :: row(:),col(:)
 real(r15), allocatable :: val(:),ph(:,:)

 nat=natom_prim_cell
 mass=atom0(1:nat)%mass
! FC2 terms -> internal relaxation of the sublattices per unit isotropic strain
 do pass=1,2
    n2=0
    do g=1,map(2)%ngr
    do ti=1,map(2)%ntind(g)
       if(map(2)%keep(counteri(2,g,ti)).ne.1) cycle
       ired = map(1)%nkeptind + current(2,g,ti)
       do t=1,map(2)%nt(g)
          if (map(2)%gr(g)%iat(1,t).gt.nat) cycle   ! first atom in the primitive cell only
          if (map(2)%gr(g)%mat(t,ti).eq.0) cycle
          n2=n2+1
          if (pass.eq.2) then
             tau=map(2)%gr(g)%iat(1,t); al=map(2)%gr(g)%ixyz(1,t)
             j=map(2)%gr(g)%iat(2,t);   ga=map(2)%gr(g)%ixyz(2,t)
             row(n2)=al+3*(tau-1); col(n2)=ga+3*(iatomcell0(j)-1)
             val(n2)=fcs(ired)*map(2)%gr(g)%mat(t,ti)
             ph(:,n2)=atompos(:,j)-atompos(:,tau)
          endif
       enddo
    enddo
    enddo
    if (pass.eq.1) allocate(row(n2),col(n2),val(n2),ph(3,n2))
 enddo
 call internal_strain_trace(nat,n2,row,col,val,ph,xir)
 if (ndn.gt.3) then   ! mech's xi is set by get_phi_zeta_Xi only for more than one atom
    dev=0
    do t=1,ndn
       dev=max(dev,abs(xir(t)-(xi(t,1,1)+xi(t,2,2)+xi(t,3,3))))
    enddo
    write(ulog,'(a,es10.3,a)')' GRUNEISEN: internal relaxation xi (per unit isotropic strain), max |xi - Tr xi(mech)|=',dev,' A'
 endif
 deallocate(row,col,val,ph)

! FC3 terms of W: Psi(tau al, j be, k ga) (x_k + xi_k0)_ga / sqrt(m m'), phase vector R_j
 do pass=1,2
    n3=0
    do g=1,map(3)%ngr
    do t=1,map(3)%nt(g)
       tau=map(3)%gr(g)%iat(1,t)
       if (tau.gt.nat) cycle
       al=map(3)%gr(g)%ixyz(1,t); be=map(3)%gr(g)%ixyz(2,t); ga=map(3)%gr(g)%ixyz(3,t)
       j=map(3)%gr(g)%iat(2,t); k=map(3)%gr(g)%iat(3,t)
       taup=iatomcell0(j); k0=iatomcell0(k)
       mi=atom0(tau)%mass; mj=atom0(taup)%mass
       y=atompos(ga,k)+xir(ga+3*(k0-1))
       do ti=1,map(3)%ntind(g)
          if (map(3)%gr(g)%mat(t,ti).eq.0) cycle
          n3=n3+1
          if (pass.eq.2) then
             ired=counteri(3,g,ti) + map(1)%ntotind + map(2)%nkeptind
             row(n3)=al+3*(tau-1); col(n3)=be+3*(taup-1)
             val(n3)=fcs(ired)*map(3)%gr(g)%mat(t,ti)*y/sqrt(mi*mj)
             ph(:,n3)=atompos(:,j)-atompos(:,taup)
          endif
       enddo
    enddo
    enddo
    if (pass.eq.1) allocate(row(n3),col(n3),val(n3),ph(3,n3))
 enddo
 call w_terms_set(n3,row,col,val,ph)
 deallocate(row,col,val,ph)

 call mode_gruneisen(get_dynmat,nkp,kp,ipath,ndn,nat,mass,volume_g0**0.33333333d0,eivl,eivc,grn,ulog)

 write(ugr2,*)'# nk , dk(nk) ,     kp(1:3,nk)     , om(la,nk),Re(gruneisen(la,nk)) '
 do ik=1,nkp
    do la=1,ndn
       if (abs(aimag(grn(la,ik))) .gt. 5d-4) &
 &        write(ulog,6)' GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)
    enddo
    write(ugr2,8)ik,dk(ik),kp(:,ik),( mysqrt(eivl(la,ik))*cnst , la=1,ndn), &
 &              (real(grn(la,ik)),la=1,ndn)
 enddo

6 format(a,2i5,99(1x,g11.4))
8 format(i8,9999(1x,f9.4))   ! one line per k-point (66+6 wrapped beyond 33 bands)
 end subroutine gruneisen
!============================================================
 subroutine get_dynmat(kp,ndn,dynmat)
!! the dynamical matrix of get_freq at kp (non-analytic term included, made Hermitian)
 use born, only : born_flag
 use constants, only : r15
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: kp(3)
 complex(r15), intent(out) :: dynmat(ndn,ndn)
 complex(r15) ddyn(ndn,ndn,3)
 if (born_flag.gt.10) then
    call set_dynamical_matrix_hard(kp,dynmat,ndn,ddyn)
 else
    call set_dynamical_matrix_std (kp,dynmat,ndn,ddyn)
 endif
 call make_hermitian(ndn,dynmat,ddyn)
 end subroutine get_dynmat
!============================================================
 subroutine get_phi_zeta_Xi(uio) !ndn,atld0,gama,phi,zeta,teta,xi,qiu,uio)
!! calculates some matrices useful for later operations, in addition to the  "standard" elastic constants
 use ios , only: ulog, write_out
 use lattice, only : volume_r0 
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use params, only : verbose
 use mech
 use eigen, only : ndyn
 use linalgb
 implicit none
 integer, intent(in) :: uio !,ndn
! real(r15), intent(out) :: gama(max(1,ndn-3),max(1,ndn-3)),phi(ndn,ndn),zeta(ndn,ndn,3), &
!&            xi(ndn,3,3),teta(ndn,ndn,3),atld0(3,3,3,3),qiu(ndn,3,3)
 integer tau,taup,al,be,ga,de,j,t,g,ti,s,cnt2,ired,la,nl,nc,nq !,cnt
 real(r15)  rij(3),junk,sumt,sump,tm(max(1,ndyn-3),max(1,ndyn-3))
 real(r15)  matr(3,3),c1(6,6),c0(6,6),aux(ndyn,ndyn),constr(3,3,3) 

 write(  * ,*)' ********** ENTERING get_phi_zeta_x *************'
 write(ulog,*)' ********** ENTERING get_phi_zeta_x *************'

 zeta=0; phi=0; teta=0; atld0=0; bulkmod=0;
 do al=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
 gloop: do g=1,map(2)%ngr
    do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
       if(map(2)%keep(counteri(2,g,ti)).ne.1) cycle
       ired = map(1)%nkeptind + current(2,g,ti)  
       if( sum(map(2)%keep(1:counteri(2,g,ti) )) .gt. map(2)%nkeptind) then
          write(*,*)'GET_PHI_ZETA: size of fc2 exceeded, ired=',ired
          stop
       endif
       tloop: do t=1,map(2)%nt(g)
          if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          if                                     ( ga .ne. map(2)%gr(g)%ixyz(2,t) ) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          rij = atompos(:,j) - atompos(:,tau)
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti) 
          nl=al+3*(tau-1) ; nc=ga+3*(taup-1)    ! dynmat dimensions
          phi(nl,nc) =phi(nl,nc)+junk ! = sum_R fc2(0,tau;R,Taup)
          do la=1,3
             teta(nl,nc,la)=teta(nl,nc,la)+junk * atompos(la,j)  ! = sum_R fc2(0,tau;R,Taup) (R+taup)_la
             zeta(nl,nc,la)=zeta(nl,nc,la)+junk * rij(la)  ! = sum_R fc2(0,tau;R,Taup) (R+taup-tau)_la
          enddo
          bulkmod=bulkmod-junk*rij(al)*rij(ga)/18d0
          do be=1,3
          do de=1,3
             atld0(al,be,ga,de)=atld0(al,be,ga,de) - junk*rij(be)*rij(de)/2
          enddo
          enddo
       enddo tloop
    enddo
 enddo gloop
 enddo
 enddo
 enddo


! qiu(tau al;ga,la) = Q_{tau al;ga la} = sum_{R taup} fc2(0 tau al;R taup ga) (R+taup-tau)_la
 qiu=0
 do la=1,3
    do al=1,3
    do ga=1,3
    do tau=1,natom_prim_cell
       do taup=1,natom_prim_cell
          nl=al+3*(tau-1) ; nc=ga+3*(taup-1)    ! dynmat dimensions
          qiu(nl,ga,la)=qiu(nl,ga,la)+zeta(nl,nc,la)
       enddo
    enddo
    enddo
    enddo
 enddo

! ASR: sum_taup phi(tau al,taup ga) = sum_{R taup} fc2(0 tau al;R taup ga) = 0
 sump=0
 do al=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
    junk=0
    do taup=1,natom_prim_cell
       junk=junk+phi(al+3*(tau-1),ga+3*(taup-1))
    enddo
    if(abs(junk).gt.1d-5) write(ulog,6)'ASR violated in phi: tau,al,ga,sum_taup phi (eV/Ang^2)=',tau,al,ga,junk
    sump=max(sump,abs(junk))
 enddo
 enddo
 enddo
 write(ulog,*)'GET_PHI_ZETA: max |sum_taup phi(tau,taup)| (ASR, eV/Ang^2)=',sump
 if(sump.gt.1d-5) then
    call warn3(   6,'ASR violated: max |sum_taup phi(tau,taup)| (eV/Ang^2)=',sump)
    call warn3(ulog,'ASR violated: max |sum_taup phi(tau,taup)| (eV/Ang^2)=',sump)
 endif

! check if phi is symmetric
 junk=maxval((phi-transpose(phi))*(phi-transpose(phi)))
 if (junk.gt.1d-12) then
    call write_out(   6,' PHI NOT SYMMETRIC ',phi)
    call write_out(ulog,' PHI NOT SYMMETRIC ',phi)
    stop
 endif
 call write_out(uio,' PHI=d^2F/d u0_tau d u0_taup=sum_R phi(tau,R+taup)  (eV/Ang^2)',phi)

 do la =1,3 ! for each la, zeta(:,:,la) is ANTIsymmetric wrt first 2 indices
    call write_out(uio,' ZETA=sum_R phi(tau,R+taup)(R+taup-tau) (eV/Ang)',zeta(:,:,la))
    aux=zeta(:,:,la)+transpose(zeta(:,:,la))
    junk=maxval(aux*aux)
    if (junk.gt.1d-12) then
       call warn5(   6,' Zeta NOT ANTISYMMETRIC:z+z^T ',aux)
       call warn5(ulog,' Zeta NOT ANTISYMMETRIC:z+z^T ',aux)
       stop
    endif
 enddo

! the antisymmetric part of qiu in its last 2 indices is fixed by the residual forces (rotational
! invariance); it is checked in residuals (from teta), where pi0 is known. Only the symmetric part
! couples to the (Cauchy) strain: qiu is symmetrized here, before xi, the residual stress and QGQ
 do s=1,ndyn
    matr=qiu(s,:,:)
    write(ulog,*)'GET_PHI_ZETA: s, max |Qiu(s,:,:)-Qiu(s,:,:)^T|=',s,maxval(abs(matr-transpose(matr)))
    qiu(s,:,:)=(matr+transpose(matr))/2
 enddo

2 format(a,2i5,9(1x,f10.4))
3 format(a,3i5,9(1x,f10.4))
4 format(a,4i5,9(1x,f10.4))
6 format(a,3i5,9(1x,g11.4))

! to get gama, invert phi: gam.phi=1a ---------------------------------------
  gam=phi  ! force constants and the images sum_R K(tau,R+taup)
  gama=0 ! gama.phi=1
  if(ndyn.gt.3) then
     gama=gam(4:ndyn,4:ndyn)
     call inverse_real(gama,tm,ndyn-3) 
     gam(1:3,:)=0 ; gam(:,1:3)=0
     gam(4:ndyn,4:ndyn)=tm  
     gama=tm  
     junk=maxval((gam-transpose(gam))*(gam-transpose(gam)))
     if (junk.gt.1d-12) then
        call write_out(   6,' |G-G^T|^2 ',junk)
        call write_out(ulog,' |G-G^T|^2 ',junk)
        call write_out(   6,' GAMA NOT SYMMETRIC ',gam)
        call write_out(ulog,' GAMA NOT SYMMETRIC ',gam)
        stop
     endif
!    call write_out(uio ,' Gamma=1/PHI (A^2/eV) ',gam)
!    call write_out(ulog,' Gamma (should be symmetric) ',gam)
  endif

! this part calculates xi(tau,ga;al,be) = du(tau,ga)/d eta_al,be = - Gam.qiu  
! should be symmetric under al <-> be since Qiu already symmetrized
! sum_tau qiu(tau ga;al,be) = 0 follows from translational and rotational invariance (checked on the
! symmetric part)
 xi=0; constr=0
 do al=1,3
 do be=1,3
 do ga=1,3
    do tau=1,natom_prim_cell
       s=3*(tau-1)+ga
       xi(s,al,be)= -dot_product(gam(s,:),qiu(:,al,be))
       constr(al,be,ga)=constr(al,be,ga)+qiu(s,al,be)
    enddo
    if(abs(constr(al,be,ga)).gt.1d-5) write(ulog,6)'sum_tau Qiu(tau ga;al,be) ne 0: al,be,ga,sum (eV/Ang)=', &
 &                                     al,be,ga,constr(al,be,ga)
 enddo
 enddo
 enddo
 write(ulog,*)'GET_PHI_ZETA: max |sum_tau Qiu(tau,:;:,:)| (eV/Ang)=',maxval(abs(constr))
 if(maxval(abs(constr)).gt.1d-5) then
    call warn3(   6,'ASR/rotational invariance violated: max |sum_tau Qiu| (eV/Ang)=',maxval(abs(constr)))
    call warn3(ulog,'ASR/rotational invariance violated: max |sum_tau Qiu| (eV/Ang)=',maxval(abs(constr)))
 endif

 atld0   = atld0   /volume_r0*ee*1d30*1d-9
 atld_raw=atld0
 write(*,*)' GET_PHI_ZETA: symmetrizing atld0'
 call symmetrize4(3,atld0) 
 bulkmod = bulkmod /volume_r0*ee*1d30*1d-9
 call convert_to_voigt(atld0,c0)

 call write_out (ulog,' GET_PHI_ZETA_XI: Bulk modulus(GPa)  ',bulkmod)
 call write_out (ulog,' Elastic Tensor (standard term; old formula) in GPa, in voigt ',c0)

 call cubic_elastic

1 format(a,99(1x,f10.4))
5 format(a,i5, 99(1x,f10.4))
 end subroutine get_phi_zeta_Xi
!============================================================
 subroutine read_stress_params(uio,s,ex)
!! stress of the reference structure from stress.params (3 lines of 3, GPa, s = dE/(V d eps),
!! tensile > 0; the file of the Huang rows), symmetrized; s=0 and ex=.false. without the file
 use ios, only : ulog, write_out
 use constants, only : r15
 use linalgb, only : symmetrize2
 use params, only : ref_relaxed
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(out) :: s(3,3)
 logical, intent(out) :: ex
 integer i
 s=0
 if(ref_relaxed) then  ! relax_reference relaxed the structure: the stress of the new reference is 0
    ex=.true.
    write(uio,*)'stress.params not used: the reference was relaxed by RELAX_REFERENCE (stress 0)'
    return
 endif
 inquire(file='stress.params',exist=ex)
 if(.not.ex) return
 open(357,file='stress.params',status='old')
 do i=1,3
    read(357,*) s(i,:)
 enddo
 close(357)
 call write_out(uio ,' stress.params: stress of the reference structure (GPa) ',s)
 call write_out(ulog,' stress.params: stress of the reference structure (GPa) ',s)
 call symmetrize2(3,s)
 end subroutine read_stress_params
!============================================================
 subroutine check_huang_stress(uio,s)
!! Huang conditions, generalized to a stressed reference (FOCEX/notes/generalized_huang.tex; the same
!! relation as the Huang rows of set_huang_inv_constraints): [ab,gl]-[gl,ab] = d_ab s_gl - d_gl s_ab,
!! [ab,gl] = -1/(2V) sum phi_ab x_g x_l = atld_raw(a,g,b,l) (GPa); s the reference stress (GPa).
!! born_flag 2: the fitted FC2 is the total minus the Born-dipole FCs, whose brackets
!! [ab,gl]_B (coulomb_lr) carry their own asymmetry; the relation is then that of the sum (as the Huang rows).
 use ios, only : ulog
 use constants, only : r15
 use mech, only : atld_raw
 use born, only : born_flag
 use lattice, only : volume_r0
 use coulomb_lr, only : lr_brackets
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(in) :: s(3,3)
 integer al,be,ga,de,delta_k
 real(r15) junk,hmax,hdev,amax,sb(3,3,3,3)
 hmax=0; hdev=0; amax=maxval(abs(atld_raw))
 sb=0
 if (born_flag.eq.2) then
    call lr_brackets(sb)
    sb=-sb/(2*volume_r0)*160.21766208d0       ! S_B (eV) -> brackets [ab,gd]_B = sb(a,b,g,d) (GPa)
    write(uio,*)'Huang check: born_flag 2, the Born-dipole brackets are added to the fitted ones'
 endif
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
    junk=atld_raw(al,ga,be,de)-atld_raw(ga,al,de,be) + sb(al,be,ga,de)-sb(ga,de,al,be)
    hmax=max(hmax,abs(junk))
    junk=junk-(delta_k(al,be)*s(ga,de)-delta_k(ga,de)*s(al,be))
    hdev=max(hdev,abs(junk))
 enddo
 enddo
 enddo
 enddo
 write(ulog,4)'Huang: max|[ab,gl]-[gl,ab]|, max|[ab,gl]-[gl,ab]-(d_ab s_gl-d_gl s_ab)|, max|[ab,gl]| (GPa)=', &
 &            hmax,hdev,amax
 write(uio ,4)'Huang: max|[ab,gl]-[gl,ab]|, max|[ab,gl]-[gl,ab]-(d_ab s_gl-d_gl s_ab)|, max|[ab,gl]| (GPa)=', &
 &            hmax,hdev,amax
 if(hdev.gt.1d-3*amax) then
    call warn3(   6,'Huang conditions (with the reference stress) violated by (GPa)',hdev)
    call warn3(ulog,'Huang conditions (with the reference stress) violated by (GPa)',hdev)
 endif
4 format(a,99(1x,g11.4))
 end subroutine check_huang_stress
!============================================================
 subroutine mechanical0(atld0,uio)
!! Born-Huang brackets -0.5*phi_ik R_jl R_jl and the elastic constants of a Bravais lattice
!! (natom_prim_cell=1); the sound-speed route for any lattice is elastic_longwave

 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use mech, only : elastic,compliance,rho_SI,atld_raw
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(out) :: atld0(3,3,3,3) 
 integer tau,taup,al,be,ga,de,j,t,g,ti,ired,la,nl,nc,ier
 logical ex
 real(r15)  rij(3),junk,constr(3,3,3),sref(3,3)

 write(  * ,*)' ********** ENTERING MECHANICAL0 *************'
 write(uio,*)' ********** ENTERING MECHANICAL0 *************'
 write(uio,*)'density in kg/m^3 is=',rho_SI

! calculation of the standard Bravais formula in atld0
 atld0=0; constr=0  ! atld0 = -0.25 sum_r,tau,taup phi(tau;R+taup)*(R+taup-tau)  ------------
 do al=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
! cnt2=0  ! cumulative number of independent terms up to group (g-1)
 gloop: do g=1,map(2)%ngr
    do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
! (this test was before the ti loop, with ti undefined: one-atom cells stopped here)
       if(map(2)%keep(counteri(2,g,ti)).ne.1) cycle
       ired = map(1)%nkeptind + current(2,g,ti)   ! this is the corresponding index of i in ared
       if( sum(map(2)%keep(1:counteri(2,g,ti))) .gt. map(2)%nkeptind) then
          write(*,*)'MECHANICAL0: size of fc2 exceeded, ired=',ired
          stop
       endif
       tloop: do t=1,map(2)%nt(g)
          if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          if( map(2)%gr(g)%ixyz(2,t) .ne. ga) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          rij = atompos(:,j) - atompos(:,tau)
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti)
          nl=al+3*(tau-1) ; nc=ga+3*(taup-1)  ! dynmat dimensions
          do la=1,3
             constr(al,ga,la)=constr(al,ga,la) + junk * rij(la)  ! = sum_R,tau,taup fc2(0,tau;R,Taup) (R+taup-tau)_la
          enddo
          do be=1,3
          do de=1,3
             atld0(al,be,ga,de)=atld0(al,be,ga,de) - junk*(rij(be)*rij(de))/2
          enddo
          enddo
       enddo tloop
    enddo
 enddo gloop
 enddo
 enddo
 enddo
 atld0=atld0/volume_r0*ee*1d30*1d-9
 atld_raw=atld0
 write(*,*)' MECHANICAL0: symmetrizing atld0'
 call symmetrize4(3,atld0) 
! stressed reference (stress.params; pi0=0 in a Bravais lattice, so sigma_res = sigma_ref):
! C from the Born-Huang brackets with the stress terms (elastic_born_huang), and the Huang check
 call read_stress_params(uio,sref,ex)
 call elastic_born_huang(atld_raw,sref,atld0)   ! replaces Sym(A) - Sym(d_ag s_bd) (wrong, 2026-11-01)
 call check_huang_stress(uio,sref)

 do la =1,3
    call write_out(uio,' CONSTR=\sum_R,tau,taup phi_ij * R_ij (eV/Ang) = 0',constr(:,:,la))
 enddo

 call convert_to_voigt(atld0,elastic)
 call write_out (uio,' Elastic Tensor (standard term; old formula) in GPa, in voigt ',elastic)

 call xmatinv(6,elastic,compliance,ier)

 if(ier.eq.0) call write_out (uio,' Compliance Tensor(1/GPa) ',compliance)

 end subroutine mechanical0
!============================================================
 subroutine residuals(uio) !ndn,xi,zeta,phi,gam,sigma0,y0,pi0,uio)
!! given phi,gam=1/phi,and xi=d_u/d_eta, calculates residual forces pi0, stresses sigma0 and displacements y0 
!! sigma0=sum pi*tau ; y0=-gama*pi ; y_eq=y0-gama*qiu*u_eq ; a'*u_eq=(qiu*gama*pi-sigma0)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use linalgb, only : symmetrize2,symmetrize_res
 use mech
 use eigen, only : ndyn
 use born, only : born_flag
 use coulomb_lr, only : lr_moments
 implicit none
 integer, intent(in) :: uio !,ndn 
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,ired,voigt,la,nl,nc,nu,s,nq,delta_k,cnt2
 real(r15)  rij(3),junk,res(3,3) ,mat2(3,3),fsum(3),rotmax,sdft(3,3)
 real(r15), allocatable :: mo(:,:,:,:)
 logical ex

 write(  * ,*)' ********** ENTERING RESIDUALS *************, ndyn=',ndyn
 write(uio,*)' ********** ENTERING RESIDUALS *************, ndyn=',ndyn

! pi0 is the residual force ------------------------------------
 pi0=0
 do tau=1,natom_prim_cell
 do al=1,3
    i=al+3*(tau-1)
    cnt2=0
    do g=1,map(1)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(1)%ntind(g-1)
    do t=1,map(1)%nt(g)
       if ( map(1)%gr(g)%ixyz(1,t) .ne. al ) cycle
       if ( map(1)%gr(g)%iat(1,t) .ne. tau ) cycle
       do ti=1,map(1)%ntind(g)
          ired = cnt2+ti    ! this is the corresponding index of i in ared
     !    ired = map(1)%ntotind + current(2,g,ti)
          pi0(i)=pi0(i)+fcs(ired) * map(1)%gr(g)%mat(t,ti)
       enddo
    enddo
    enddo
 enddo
    write(uio,3)'Residual -force(eV/Ang): tau,pi0(tau,:)=',tau,(pi0(3*(tau-1)+al),al=1,3)
 enddo

! ASR of rank 1: total force on the primitive cell sum_tau pi0(tau,al) = 0
 do al=1,3
    fsum(al)=sum(pi0(al:ndyn:3))
 enddo
 write(ulog,4)'RESIDUALS: sum_tau pi0(tau,:) (ASR, eV/Ang)=',fsum
 if(maxval(abs(fsum)).gt.1d-5) then
    call warn4(   6,'ASR violated: sum_tau pi0(tau,:) (eV/Ang)=',fsum)
    call warn4(ulog,'ASR violated: sum_tau pi0(tau,:) (eV/Ang)=',fsum)
 endif

3 format(a,i4,99(1x,f11.5))
4 format(a,99(1x,g11.4))

 y0=0 ! y0 = -Gama*pi correction to equilibrium positions at eta=0----------------------
 do tau=2,natom_prim_cell
    do al=1,3
    i=3*(tau-1)+al
    do j=4,ndyn
       y0(i)=y0(i)-gama(i-3,j-3)*pi0(j)
    enddo
    enddo
    write(uio,3)' Residual displacement, y0(tau,:)=-gama*Pi(tau)=',tau,(y0(3*(tau-1)+al),al=1,3)
 enddo

! residual Stress tensor (under no strain) ------------------------------------
 sigma0=0
 do al=1,3
 do be=1,3
    do tau=1,natom_prim_cell
         sigma0(al,be)=sigma0(al,be)+atompos(be,tau)*pi0(3*(tau-1)+al) !  &
    enddo
 enddo
 enddo
 call write_out(ulog,' sigma(eta=0) before symmetrization (eV) ',sigma0)
! should be symmetric according to rotational invariance
 call symmetrize2(3,sigma0)   ! equivalent to enforcing rotational invariance on pi0
! sum_tau pi0 tau depends on which image of each atom is in the primitive cell when pi0 /= 0, and
! is zero for a Bravais lattice under any stress. The stress of the reference structure from the
! DFT calculation, if given in stress.params (3x3, GPa, sigma=dE/(V d eps): tensile > 0), replaces it.
 call read_stress_params(uio,sdft,ex)
 if(ex) then
    call write_out(uio ,' FC1 estimate sum_tau pi0 tau/V (GPa), replaced by stress.params ', &
 &                 sigma0/volume_r0*1d30*ee*1d-9)
    sigma0 = sdft*volume_r0/(1d30*ee*1d-9)   ! in eV
 else
    write(uio,*)'RESIDUALS: no stress.params; reference stress = sum_tau pi0 tau/V (FC1 only)'
 endif

 qiuv=0
 do al=1,3
 do be=al,3
    sigmav(voigt(al,be))=sigma0(al,be)
    qiuv(:,voigt(al,be))=qiu(:,al,be)
 enddo
 enddo
 do la=1,6
    sigmav(la)=sigmav(la)+dot_product(qiuv(:,la),y0)
 enddo
 call write_out(ulog,'Residual stress sigma(eta) (eV) ',sigmav)
 sigma0 = sigma0/volume_r0*1d30*ee*1d-9 
 sigmav = sigmav/volume_r0*1d30*ee*1d-9 
 call write_out(uio,' sigma(eta=0) after symmetrization (GPa) ',sigma0)
 call write_out(uio,' sigma(eta  ) = sigma0-Pi Gamma Qiu(GPA) ',sigmav)

! Huang conditions with the reference stress of stress.params (0 without it, as in the fit)
 call check_huang_stress(uio,sdft)
 call write_out(uio,' Qiu(1:ndyn,1:6)=d^2F/d u0_tau d eta (GPa) ',qiuv)

! rotational invariance of rank 2 (theory.rst): sum_{R taup} fc2(0 tau al;R taup be)(R+taup)_ga eps(be,ga,nu)
! + pi0(tau,be) eps(be,al,nu) = 0, i.e. Qiu(tau al;be,ga)-Qiu(tau al;ga,be) + pi0(tau,be) delta(al,ga)
! - pi0(tau,ga) delta(al,be) = 0  (teta and qiu sums differ by tau*sum_taup phi = 0 by ASR)
! enforces mat2(a,b)-mat2(b,a)=res(a,b)-res(b,a); consider mat=mat2-res ; symmetrize mat then mat2=sym(mat)+res
! born_flag 2: the fitted FC2 is the total minus the Born-dipole FCs; the relation holds
! for the sum (as the rotational rows, invariance_ranks), so their moments mo(tau,al,be,ga) = sum_j Phi_B(tau al,j be)
! (x_j-x_tau)_ga (coulomb_lr) are added (x_tau sum_j Phi_B = 0 by ASR, as for teta)
 allocate(mo(natom_prim_cell,3,3,3)); mo=0
 if (born_flag.eq.2) then
    call lr_moments(mo)
    write(uio,*)'Rotational check: born_flag 2, the Born-dipole moments are added to the fitted ones'
 endif
 s=0; rotmax=0
 do al=1,3
    do tau=1,natom_prim_cell
       do be=1,3
       do ga=1,3
          junk=mo(tau,al,be,ga)-mo(tau,al,ga,be)
          do taup=1,natom_prim_cell
            junk=junk+teta(3*(tau-1)+al,3*(taup-1)+be,ga)-teta(3*(tau-1)+al,3*(taup-1)+ga,be)
          enddo
          junk=junk + (pi0(3*(tau-1)+be)*delta_k(al,ga)-pi0(3*(tau-1)+ga)*delta_k(al,be))
          if(abs(junk).gt.1d-5) then
             write(ulog,'(a,4i4,g12.4)')'RESIDUALS: rot invce violation in theta:tau,al,be,ga,rot (eV/Ang)=',tau,al,be,ga,junk
             s=1
          endif
          rotmax=max(rotmax,abs(junk))
       enddo
       enddo
    enddo
 enddo
 deallocate(mo)
 write(ulog,4)'RESIDUALS: max rotational invariance violation of Qiu and pi0 (eV/Ang)=',rotmax
 if(s.ne.0) then
    call warn3(   6,'rotational invariance violated: max |Q(ab)-Q(ba)+pi0 terms| (eV/Ang)=',rotmax)
    call warn3(ulog,'rotational invariance violated: max |Q(ab)-Q(ba)+pi0 terms| (eV/Ang)=',rotmax)
 endif

! if(s.ne.0) then 
!   write(*,*)' symmetrizing now teta '
!! Impose rotational invariance on teta: sum_taup teta(tau,al;taup,be;ga) + pi(taup,be) delta(al,ga) symm in be<->ga
!    do al=1,3
!    do tau=1,natom_prim_cell
!       la=al+3*(tau-1)
!       do be=1,3
!       do ga=1,3
!          res(be,ga)= -( pi0(3*(tau-1)+be)*delta_k(al,ga)-pi0(3*(tau-1)+ga)*delta_k(al,be) )
!          mat2(be,ga)=0
!          do taup=1,natom_prim_cell
!            mat2(be,ga)=mat2(be,ga)+teta(3*(tau-1)+al,3*(taup-1)+be,ga)-teta(3*(tau-1)+al,3*(taup-1)+ga,be)
!          enddo
!       enddo
!       enddo
!       call symmetrize_res(mat2,res)
!    enddo
!    enddo
!
! endif

 write(uio,*)'Suggested relaxations: eta=1/(A-QGQ) (QGPi-sigma0); u0=-G(Pi+Q eta)'

 end subroutine residuals
!===================================================
 subroutine mechanical(uio) !ndn,atld1,sigma0,phi,zeta,xi,qiu,gama,elastic,uio)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use mech
 use eigen, only : ndyn
 implicit none
 integer, intent(in) :: uio !,ndn 
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,ired,la,nl,nc,nu,s,nq,delta_k,ier
 real(r15) c1(6,6),c2(6,6),c3(6,6),cq(6,6) ,atld2(3,3,3,3),atld3(3,3,3,3),qgq(3,3,3,3),sres(3,3)


 call convert_to_voigt(atld0,c1)
 call write_out (uio,' Elastic Tensor (first term in GPa) in voigt ',c1)

 qgq=0 

 atld2=0    ! cross terms : phi*(R+tau)*Xi = zeta*delta (Xi) needs symmtrization --------------------
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
    qgq(al,be,ga,de)=dot_product(qiu(4:ndyn,al,be),matmul(gama,qiu(4:ndyn,ga,de)))
 enddo
 enddo
 enddo
 enddo
 qgq  =qgq  /volume_r0*ee*1d30*1d-9
 write(*,*)' symmetrizing QGQ'
 call symmetrize4(3,qgq) 

  write(uio,*)'checking symmetry for QGQ'
  do al=1,3
  do be=1,3
  do ga=1,3
  do de=1,3
    if (abs(qgq(al,be,ga,de)-qgq(ga,de,al,be)).gt.1d-4)  &
 &       write(uio,*)'ab<->gd ',al,be,ga,de,qgq(al,be,ga,de),qgq(ga,de,al,be)
    if (abs(qgq(al,be,ga,de)-qgq(be,al,ga,de)).gt.1d-4)  &
 &       write(uio,*)' a<->b ',al,be,ga,de,qgq(al,be,ga,de),qgq(be,al,ga,de)
    if (abs(qgq(al,be,ga,de)-qgq(al,be,de,ga)) .gt.1d-4)  &
 &       write(uio,*)' g<->d ',al,be,ga,de,qgq(al,be,ga,de),qgq(al,be,de,ga)
  enddo
  enddo
  enddo
  enddo

 call convert_to_voigt(qgq,cq)
 call write_out (uio,' Elastic Tensor correction QGQ in voigt ',cq)

! elastic constants (Lagrangian strain) of a possibly stressed reference: Born-Huang brackets of the
! unsymmetrized A plus the stress terms (elastic_born_huang), minus QGQ; s = relaxed residual stress sigmav
 sres(1,1)=sigmav(1); sres(2,2)=sigmav(2); sres(3,3)=sigmav(3)
 sres(2,3)=sigmav(4); sres(1,3)=sigmav(5); sres(1,2)=sigmav(6)
 sres(3,2)=sres(2,3); sres(3,1)=sres(1,3); sres(2,1)=sres(1,2)
 call elastic_born_huang(atld_raw,sres,atld0)   ! replaces Sym(A) - Sym(d_ag s_bd) (wrong, 2026-11-01)
 
  write(uio,*)'checking symmetry for atld0'
  do al=1,3
  do be=1,3
  do ga=1,3
  do de=1,3
    if (abs(atld0(al,be,ga,de)-atld0(ga,de,al,be)).gt.1d-4)  &
 &       write(uio,*)'ab<->gd ',al,be,ga,de,atld0(al,be,ga,de),atld0(ga,de,al,be)
    if (abs(atld0(al,be,ga,de)-atld0(be,al,ga,de)).gt.1d-4)  &
 &       write(uio,*)' a<->b ',al,be,ga,de,atld0(al,be,ga,de),atld0(be,al,ga,de)
    if (abs(atld0(al,be,ga,de)-atld0(al,be,de,ga)) .gt.1d-4)  &
 &       write(uio,*)' g<->d ',al,be,ga,de,atld0(al,be,ga,de),atld0(al,be,de,ga)
  enddo
  enddo
  enddo
  enddo

 atld0=atld0-qgq
! call symmetrize4(3,atld0)
! these should be symmetric wr (al,be <-> ga,de)
! can(?) symmetrize wr (al<->be) and (ga<->de) 
!
 call convert_to_voigt(atld0,elastic)
 call write_out (uio,' Elastic Tensor A0-QGQ (GPa) ',elastic)

! to avoid signular matrix for 2D (in xy plane)
 if (elastic(3,3) .myeq. 0d0) elastic (3,3)=1d-2
 if (elastic(4,4) .myeq. 0d0) elastic (4,4)=1d-2
 if (elastic(5,5) .myeq. 0d0) elastic (5,5)=1d-2
 if (elastic(6,6) .myeq. 0d0) elastic (6,6)=1d-2

 call xmatinv(6,elastic,compliance,ier)

 call write_out (uio,' Compliance Tensor(1/GPa) ',compliance)

 u0v=-matmul(compliance,sigmav)

 call write_out (uio,' Total residual strain',u0v)

 y0(4:ndyn) =y0(4:ndyn)-matmul(gama,matmul(qiuv(4:ndyn,:),u0v))

 call write_out (uio,' Total residual displacements (Ang)',y0)

 end subroutine mechanical
!============================================================
 subroutine sound_speeds(q,calbe,vkms)
!! uses the elastic constant tensor to find sound speeds along q
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use eigen, only : mysqrt
 use constants, only : r15,ee,uma,hbar,k_b,pi
 use mech
 implicit none
 real(r15), intent(in) :: calbe(3,3,3,3),q(3) !elastic tensor in GPa
 real(r15), intent(out) :: vkms(3)
 integer i,j,k,l,ier,ndn
 real(r15)  rj(3),junk,qhat(3),vl,vt,va,td 
 real(r15) vbar, eivl(3) ,lame_la,lame_mu
 complex(r15)  eivc(3,3)
 complex(r15) dynr2(3,3)

 write(ulog,*)' ********** ENTERING SOUND SPEEDS *************'
 qhat= q/length(q) !(g01+g02)/length(g01+g02)
 ndn=3

  dynr2=0
  do i=1,3
  do j=1,3
     do k=1,3
     do l=1,3
        dynr2(i,j)=dynr2(i,j)+calbe(i,k,l,j)*qhat(l)*qhat(k)*1d9   ! to convert to Pa
     enddo
     enddo
  enddo
  enddo

! diagonalize dynr2 to get speeds of sound
  call diagonalize(ndn,dynr2,eivl,ndn,eivc,ier)

! eivals are in units of phi*R^2, i.e. in eV
  vkms(:)=mysqrt(eivl(:))/sqrt(rho_SI)/1000 ! to convert to (km/s)
  write(ulog,3)' sound speeds along ',q,' in (km/s)=',vkms

  vbar = (3d0/sum(1d0/vkms**3))**0.3333333
  td=hbar/k_b*(6*pi*pi*natom_prim_cell/volume_r0)**0.3333333 *1d10 * vbar*1000
  write(ulog,8)' vbar(km/s), Debye Temperature    =',vbar,td

3 format(2(a,3(1x,g11.4)))
8 format(a,9(1x,g11.4))

 end subroutine sound_speeds
!============================================================
 subroutine elastic_longwave(uio,rotc)
!! elastic constants and the deviatoric stress from the long-wave limit of the code's own dynamical matrix
!! (the one of the band structure: fitted FC2 + Ewald/NA term), 2026-11-02 (user's request).
!! Along 13 directions qh (axes, face and body diagonals) the acoustic matrix A(qh) (omega^2 = q^2 eig A,
!! optical modes folded out exactly; gruneisen_rules acoustic_tensor) gives
!!   Gamma_ik = rho A_ik = sum_jl B[(ik),(jl)] qh_j qh_l = (C_ijkl + d_ik s_jl) qh_j qh_l,  rho = M/V,
!! B = the 36 relaxed Born-Huang brackets, fitted by least squares (78 equations, rank 36). Then
!!   B[ab,cd] - B[cd,ab] = d_ab s_cd - d_cd s_ab   (stressed Huang relation, force-free reference)
!! gives the deviatoric stress s_dev (15 equations, 5 unknowns; the pressure is not seen by long waves),
!! and C (21) follows from B - d_ik s_jl with s = s_dev + the pressure of FOCEX's own residual stress.
!! Checked on the exact model (REFERENCE/exact_model/longwave): s_dev to 1e-13 GPa, C to 5e-10 GPa
!! (ZnO, HfO2 stressed, CaF2), and against the strain Hessian at the relaxed ZnO model (1e-3 GPa).
!! With residual forces (FC1 /= 0) the asymmetry also contains a force term (antisymmetric part of Q).
!! Two results when Born charges are on (born_flag 2,3,6,7,8, ewald_2nd_deriv_hat):
!!   E = 0 : the non-analytic G = 0 (macroscopic-field) term left out (ewald_no_macro): the relaxed-ion
!!           elastic constants at zero field, as VASP reports them;
!!   with NA: D as in the dispersion, i.e. C^E stiffened by the ionic piezoelectric term only (no
!!           electronic clamped-ion part in FOCEX), direction dependent: for reference.
 use ios, only : ulog, write_out
 use lattice, only : volume_r0,volume_g0
 use atoms_force_constants, only : natom_prim_cell,atom0
 use born, only : born_flag
 use ewald, only : ewald_no_macro
 use mech, only : sigmav
 use constants, only : r15,ee
 use gruneisen_rules, only : dynmat_routine, acoustic_tensor
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(in) :: rotc(3,3)
 procedure(dynmat_routine) :: get_dynmat
 integer, parameter :: nd=13
 real(r15) dirs(3,nd),qh(3),mass(natom_prim_cell),rho,conv,gsize,gam(3,3,nd),eg(3), &
 &   ata(36,36),atai(36,36),atb(36),row(36),x(36),bb(6,6),b4(3,3,3,3),dl(3,3), &
 &   sm(6,6),smi(6,6),sb(6),sr(6),sl(3,3),sf(3,3),sfd(3,3),pr,res,nrm, &
 &   ca(21,21),cai(21,21),cb(21),cr(21),cx(21),cv(6,6),c4(3,3,3,3),cs(3,3,3,3),cro(3,3,3,3),cc(6,6),dev,dif
 integer i,j,k,l,a,b,c,d,n,ns,mode,nmode,ip(6,6),np,ii,jj,voigt,ier,vij,vil
 complex(r15) am(3,3),ev(3,3)
 logical e0,ex
 character(len=40) lab

 dirs=reshape((/ 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0, 1d0,1d0,0d0, 1d0,-1d0,0d0, 1d0,0d0,1d0, &
 &   1d0,0d0,-1d0, 0d0,1d0,1d0, 0d0,1d0,-1d0, 1d0,1d0,1d0, 1d0,1d0,-1d0, 1d0,-1d0,1d0, -1d0,1d0,1d0 /),(/3,nd/))
 do i=1,natom_prim_cell
    mass(i)=atom0(i)%mass
 enddo
 rho=sum(mass)/volume_r0                  ! amu/A^3; rho A in eV/A^3
 conv=ee*1d30*1d-9                         ! eV/A^3 -> GPa
 gsize=volume_g0**0.33333333d0
 dl=0; dl(1,1)=1; dl(2,2)=1; dl(3,3)=1
 np=0                                       ! independent pairs (I<=J) of the Voigt-symmetric C
 do ii=1,6
 do jj=ii,6
    np=np+1; ip(ii,jj)=np; ip(jj,ii)=np
 enddo
 enddo
! FOCEX's own stress (GPa): the relaxed residual stress of residuals, or stress.params for one atom per cell
 sf=0
 if (3*natom_prim_cell.gt.3) then
    sf(1,1)=sigmav(1); sf(2,2)=sigmav(2); sf(3,3)=sigmav(3)
    sf(2,3)=sigmav(4); sf(1,3)=sigmav(5); sf(1,2)=sigmav(6)
    sf(3,2)=sf(2,3); sf(3,1)=sf(1,3); sf(2,1)=sf(1,2)
 else
    call read_stress_params(ulog,sf,ex)
 endif
 pr=(sf(1,1)+sf(2,2)+sf(3,3))/3
 sfd=sf-pr*dl
 e0 = born_flag.eq.2 .or. born_flag.eq.3 .or. born_flag.eq.6 .or. born_flag.eq.7 .or. born_flag.eq.8
 nmode=1; if (e0) nmode=2

 write(uio,*)'ELASTIC_LONGWAVE: C_ij and deviatoric stress from the sound speeds (long-wave limit of D, 13 directions)'
 if (born_flag.ne.0 .and. .not.e0) write(uio,*)'ELASTIC_LONGWAVE: E = 0 switch not available for born_flag=', &
 &   born_flag,'; the result includes the non-analytic term'
 do mode=1,nmode
    ewald_no_macro = (mode.eq.1 .and. e0)
    if (born_flag.eq.0) then
       lab='(no Born charges)'
    elseif (ewald_no_macro) then
       lab='E = 0 (G=0 macroscopic term left out)'
    else
       lab='with the non-analytic term (as bands)'
    endif
! Christoffel tensors and the fit of the 36 brackets B[(ik),(jl)]
    ata=0; atb=0; nrm=0
    do n=1,nd
       qh=dirs(:,n)/sqrt(dot_product(dirs(:,n),dirs(:,n)))
       call acoustic_tensor(get_dynmat,qh,3*natom_prim_cell,natom_prim_cell,mass,gsize,am,ns)
       gam(:,:,n)=real(am)*rho*conv
       if (ns.ne.3) write(ulog,*)'ELASTIC_LONGWAVE: WARNING: near-zero optical modes at Gamma, ns=',ns
       am=gam(:,:,n)
       call diagonalize(3,am,eg,0,ev,ier)
       write(ulog,'(a,3f8.4,a,3f12.4)')' ELASTIC_LONGWAVE: qh=',qh,'  rho v^2 (GPa)=',eg
       if (minval(eg).lt.0) write(uio,'(a,3f8.4,a,3f11.3)')' ELASTIC_LONGWAVE: WARNING: rho v^2 < 0 along ', &
 &        qh,' (GPa): ',eg
       do i=1,3
       do k=i,3
          row=0
          do j=1,3
          do l=1,3
             row(6*(voigt(i,k)-1)+voigt(j,l))=row(6*(voigt(i,k)-1)+voigt(j,l))+qh(j)*qh(l)
          enddo
          enddo
          do ii=1,36
             ata(ii,:)=ata(ii,:)+row(ii)*row
          enddo
          atb=atb+row*gam(i,k,n)
          nrm=nrm+gam(i,k,n)**2
       enddo
       enddo
    enddo
    call inverse_real(ata,atai,36)
    x=matmul(atai,atb)
    do ii=1,6
    do jj=1,6
       bb(ii,jj)=x(6*(ii-1)+jj)
    enddo
    enddo
    do a=1,3; do b=1,3; do c=1,3; do d=1,3
       b4(a,b,c,d)=bb(voigt(a,b),voigt(c,d))
    enddo; enddo; enddo; enddo
! deviatoric stress from B[ab,cd]-B[cd,ab] = d_ab s_cd - d_cd s_ab, with Tr s = 0 (least squares, 6 unknowns)
    sm=0; sb=0; res=0
    do a=1,3; do b=1,3; do c=1,3; do d=1,3
       sr=0
       do ii=1,6
          sr(ii)=dl(a,b)*e6(ii,c,d)-dl(c,d)*e6(ii,a,b)
       enddo
       do ii=1,6
          sm(ii,:)=sm(ii,:)+sr(ii)*sr
       enddo
       sb=sb+sr*(b4(a,b,c,d)-b4(c,d,a,b))
    enddo; enddo; enddo; enddo
    sr=(/1d0,1d0,1d0,0d0,0d0,0d0/)
    do ii=1,6
       sm(ii,:)=sm(ii,:)+sr(ii)*sr
    enddo
    call inverse_real(sm,smi,6)
    sr=matmul(smi,sb)
    sl(1,1)=sr(1); sl(2,2)=sr(2); sl(3,3)=sr(3); sl(2,3)=sr(4); sl(1,3)=sr(5); sl(1,2)=sr(6)
    sl(3,2)=sl(2,3); sl(3,1)=sl(1,3); sl(2,1)=sl(1,2)
    do a=1,3; do b=1,3; do c=1,3; do d=1,3
       res=res+(b4(a,b,c,d)-b4(c,d,a,b)-dl(a,b)*sl(c,d)+dl(c,d)*sl(a,b))**2
    enddo; enddo; enddo; enddo
! C (21) from B - d_ik s_jl = (C_ijkl + C_ilkj)/2, s = s_dev(long wave) + the pressure of FOCEX's stress
    ca=0; cb=0
    do i=1,3; do k=1,3; do j=1,3; do l=1,3
       cr=0
       vij=ip(voigt(i,j),voigt(k,l)); vil=ip(voigt(i,l),voigt(k,j))
       cr(vij)=cr(vij)+0.5d0; cr(vil)=cr(vil)+0.5d0
       do ii=1,21
          ca(ii,:)=ca(ii,:)+cr(ii)*cr
       enddo
       cb=cb+cr*(b4(i,k,j,l)-dl(i,k)*(sl(j,l)+pr*dl(j,l)))
    enddo; enddo; enddo; enddo
    call inverse_real(ca,cai,21)
    cx=matmul(cai,cb)
    do ii=1,6
    do jj=1,6
       cv(ii,jj)=cx(ip(ii,jj))
    enddo
    enddo
    do i=1,3; do j=1,3; do k=1,3; do l=1,3
       c4(i,j,k,l)=cv(voigt(i,j),voigt(k,l))
    enddo; enddo; enddo; enddo
    dif=0
    do i=1,3; do k=1,3; do j=1,3; do l=1,3
       dif=dif+(b4(i,k,j,l)-dl(i,k)*(sl(j,l)+pr*dl(j,l))-(c4(i,j,k,l)+c4(i,l,k,j))/2)**2
    enddo; enddo; enddo; enddo
    call symmetrize_tensor4(c4,cs,dev)
    call rotate_tensor4(rotc,cs,cro)
    call convert_to_voigt(cro,cc)
    write(uio,*)'ELASTIC_LONGWAVE: ',trim(lab)
    call write_out(uio,' ELASTIC_LONGWAVE: brackets B[(ik),(jl)] (GPa), rows ik, columns jl, Voigt order',bb)
    write(uio,'(a,es10.2)')' ELASTIC_LONGWAVE: max |B - B^T| (GPa) =',maxval(abs(bb-transpose(bb)))
    call write_out(uio,' ELASTIC_LONGWAVE: deviatoric stress from the asymmetry of B (GPa)',sl)
    call write_out(uio,' ELASTIC_LONGWAVE: deviatoric part of FOCEX residual stress (GPa)',sfd)
    write(uio,'(a,es10.2,a,f10.4,a)')' ELASTIC_LONGWAVE: rms misfit of the asymmetry relation (GPa) =', &
 &     sqrt(res/81),';  pressure part used for C: Tr(s)/3 =',pr,' GPa (FOCEX residual stress)'
    write(uio,'(a,es10.2,a,f8.3)')' ELASTIC_LONGWAVE: rms misfit of the 21-C fit to B - d s (GPa) =', &
 &     sqrt(dif/81),'; symmetry-forbidden part of C removed (GPa)=',dev
    write(uio,*)'ELASTIC_LONGWAVE: CRYSTAL FRAME C_ij (GPa)'
    do i=1,6
       write(uio,'(6(1x,f10.3))') (merge(cc(i,j),0d0,abs(cc(i,j)).ge.1d-3*maxval(abs(cc))),j=1,6)
    enddo
 enddo
 ewald_no_macro=.false.

 contains
 real(r15) function e6(m,p,q)
!! (p,q) component of the symmetric unit tensor of Voigt index m
 integer, intent(in) :: m,p,q
 e6=0
 if (voigt(p,q).eq.m) e6=1
 end function e6
 end subroutine elastic_longwave
!============================================================
 subroutine QHA_free_energy_vs_strain(nk,wk,ndyn,eival,grin,temprk,eq_free,eq_etot,eq_cv,eq_cp,eq_vol,entropy,cvte,eq_grun,boft,uio)
!! Quasi-harmonic thermodynamics at temperature temprk (rewritten 2026-11-01 with the user):
!!   F(V) = E_static(V) + F_phonon(V,T) [+ F_agne] on 61 isotropic strains -5%..+10%; omega(V) from gamma
!!   E_static: Murnaghan EOS, B0 = Reuss modulus (exact for a hydrostatic load), B0' from FC2+FC3
!!     (pressure_derivative_b) or 1 + 2 gamma (Dugdale-MacDonald, Leont'ev gamma) without FC3; P0 = residual pressure
!!   P(V) = -dF/dV is summed analytically at each strain; V(T) is the root of P = 0 (local cubic fit),
!!   B_T = -V dP/dV at V(T) (the isothermal bulk modulus), alpha_V = (dP/dT)_V / B_T = d ln V/dT
!!   (thermodynamic identity; (dP/dT)_V summed analytically), linear CTE = alpha_V/3 = d ln L/dT,
!!   Cp = Cv + V T B_T alpha_V^2. Cv, <gamma> (Cv-weighted), E, F, S interpolated to V(T).
!!   Harmonic-only fits (all gamma = 0): every mode gets Leont'ev's gamma = 1.5(3-4x^2)/(1+2x^2), x = vt/vl,
!!   plus P_agne = NkT/V (F_agne = -NkT ln(V/V0)).
!!   Until 2026-11-01: V(T) as a Boltzmann-weighted average over the strains with an artificial kT, the CTE
!!   as the secant eq_strain/T, B(T) = B0(1 - strain (4+3<gamma>)), static energy quadratic in dV (B0' = -1).
!! Outputs: eq_vol = V(T)/V0, cvte = alpha_V (1/K), boft = B_T (GPa), eq_* per mole of atoms.
 use ios
 use params
 use lattice
 use constants
 use atoms_force_constants, only : natom_prim_cell
 use mech, only : sigma0,rbulkmod,b0p,xagne
 use eigen, only : mysqrt
 implicit none
 integer, intent(in) :: nk,ndyn,uio
 real(r15), intent(in) :: wk(nk),eival(ndyn,nk),grin(ndyn,nk), temprk
 real(r15), intent(out):: cvte,eq_free,eq_vol,eq_cv,eq_cp,eq_etot,entropy,eq_grun
 real(r15), intent(inout):: boft
 integer, parameter :: imax=20, ng=3*imax+1
 integer b,k,nat,i,i0,j,it
 real(r15) x0,x1,cv_nk,hw,nbe,pres0,dvbyv,ev2kjpermol,geff,gr2,bp,b0,v0,vr,fstat,pstat,cvsum,grsum
 real(r15) strain(ng),free(ng),pres(ng),etot(ng),cvv(ng),grv(ng),dpdt(ng),bloc(ng), &
 &         fst(ng),pst(ng),fag(ng),pag(ng),grn(ndyn,nk)
 real(r15) xs(4),c(4),seq,dpds,veq,bt,alv,w(4),f_eq,e_eq,cv_eq,dpdt_eq,gsum_eq
 logical agne
 real(r15), save :: vol0=0
 logical, save :: first=.true.

 nat = ndyn/3
 ev2kjpermol = n_avog/(nat*1000d0)   ! J per cell -> kJ per mole of atoms
 if (temprk.le.0) then
    write(ulog,*)'temperature not in the proper range!!',temprk
    stop
 endif
 v0=volume_r0*1d-30                   ! m^3
 pres0=-(sigma0(1,1)+sigma0(2,2)+sigma0(3,3))/3   ! GPa, residual pressure of the reference
 b0=rbulkmod
 agne=.false.
 if (maxval(abs(grin)).lt.1d-6) then  ! harmonic theory only
    agne=.true.
    gr2=1.5*(3-4*xagne*xagne)/(1+2*xagne*xagne)
    grn=gr2
 else
    grn=grin
 endif
! B0': from FC2+FC3 (pressure_derivative_b); without FC3 the Dugdale-MacDonald relation
! gamma = B0'/2 - 1/2 with Leont'ev's gamma (user, 2026-11-02; was 4)
 bp=b0p
 if (include_fc(3).eq.0) bp=1+2*gr2
 if (abs(bp).lt.1d-6) bp=1d-6
 if (abs(bp-1).lt.1d-6) bp=1+1d-6
 if (first) then
    write(ulog,'(a,f12.4,a,f10.5,a,f10.4,a,l2)')' QHA: B0 (Reuss, GPa)=',b0,'  B0p=',bp,'  P0 (GPa)=',pres0, &
 &        '  harmonic-only (Leont''ev gamma + P_agne):',agne
    if (include_fc(3).eq.0) write(ulog,*)' QHA: no FC3: B0'' = 1 + 2 gamma (Dugdale-MacDonald, Leont''ev gamma)'
    first=.false.
 endif

 do i=1,ng
    strain(i)= (i-1-imax)/(20d0*imax)
    vr=(1+strain(i))**3                   ! V/V0
    dvbyv=vr-1
! static Murnaghan: P = P0 + (B0/B')[(V0/V)^B' - 1], E = -int P dV (E(V0)=0)
    pstat=pres0 + b0/bp*(vr**(-bp)-1)
    fstat=( b0/bp*( (vr**(1-bp)-1)/(bp-1) + dvbyv ) - pres0*dvbyv )*v0*1d9 * ev2kjpermol
    free(i)=fstat; pres(i)=pstat; etot(i)=fstat
    fst(i)=fstat; pst(i)=pstat; fag(i)=0; pag(i)=0; dpdt(i)=0
    if (agne) then
       pag(i)=k_b*natom_prim_cell*temprk/(v0*vr)*1d-9
       fag(i)=-k_b*natom_prim_cell*temprk*log(vr)*ev2kjpermol
       pres(i)=pres(i)+pag(i); free(i)=free(i)+fag(i)
       dpdt(i)=k_b*natom_prim_cell/(v0*vr)*1d-9
    endif
    cvsum=0; grsum=0
    do k=1,nk
    do b=1,ndyn
       if (eival(b,k).eq.0) cycle   ! rigid translation at k=0 (freq_vel): not a vibration
       x0=h_plank*mysqrt(abs(eival(b,k)))*cnst*100*c_light/(k_b*temprk)
       if (x0.gt.0) then
          if (strain(i)*grn(b,k).lt.0) then
             x1=x0*(1-3*grn(b,k)*strain(i))
             geff=grn(b,k)*(1+strain(i))/(1-3*grn(b,k)*strain(i))
          else
             x1=x0/(1+3*grn(b,k)*strain(i))
             geff=grn(b,k)*(1+strain(i))/(1+3*grn(b,k)*strain(i))
          endif
          if (x1.le.x0/100) then
             x1=x0/100d0; geff=0
          endif
       else
          x1=abs(x0)/100; geff=0
       endif
       if (x1.gt.60) then
          cv_nk=0; hw=x1*k_b*temprk
          pres(i)=pres(i)+hw*0.5d0*wk(k)*geff/(v0*vr)*1d-9
          etot(i)=etot(i)+hw*0.5d0*wk(k)*ev2kjpermol
          free(i)=free(i)+k_b*temprk*(x1/2)*wk(k)*ev2kjpermol
       else
          cv_nk=x1*x1*nbe(x1,1d0,classical)*(1+nbe(x1,1d0,classical))
          hw=x1*k_b*temprk
          pres(i)=pres(i)+hw*(nbe(x1,1d0,classical)+0.5d0)*wk(k)*geff/(v0*vr)*1d-9
          etot(i)=etot(i)+hw*(nbe(x1,1d0,classical)+0.5d0)*wk(k)*ev2kjpermol
          free(i)=free(i)+k_b*temprk*(x1/2+log(1-exp(-x1)))*wk(k)*ev2kjpermol
       endif
       cvsum=cvsum+cv_nk*wk(k)
       grsum=grsum+cv_nk*wk(k)*grn(b,k)
       dpdt(i)=dpdt(i)+k_b*cv_nk*wk(k)*geff/(v0*vr)*1d-9      ! (dP/dT)_V, GPa/K
    enddo
    enddo
    cvv(i)=cvsum
    grv(i)=grsum/max(cvsum,1d-300)
 enddo
! B_T at each strain: -V dP/dV = -(1+s)/3 dP/ds (central differences)
 do i=1,ng
    j=min(max(i,2),ng-1)
    bloc(i)=-(1+strain(i))/3*(pres(j+1)-pres(j-1))/(strain(j+1)-strain(j-1))
 enddo
! V(T): root of P(s) = 0 by a cubic through the 4 points around the sign change
 i0=0
 do i=1,ng-1
    if (pres(i).ge.0 .and. pres(i+1).lt.0) then
       i0=i; exit
    endif
 enddo
 if (i0.eq.0) then
    write(ulog,*)' QHA: WARNING: P(V) has no root in the strain range at T=',temprk,'; endpoint used'
    i0=merge(1,ng-1,pres(1).lt.0)
    seq=strain(i0)
 else
    j=min(max(i0-1,1),ng-3)
    xs=strain(j:j+3)
    seq=strain(i0)+(strain(i0+1)-strain(i0))*pres(i0)/(pres(i0)-pres(i0+1))
    do it=1,30                        ! Newton on the Lagrange cubic of P
       call lagrange4(xs,seq,w,c)
       dpds=dot_product(c,pres(j:j+3))
       seq=seq-dot_product(w,pres(j:j+3))/dpds
    enddo
 endif
 j=min(max(i0-1,1),ng-3)
 xs=strain(j:j+3)
 call lagrange4(xs,seq,w,c)
 dpds=dot_product(c,pres(j:j+3))
 bt=-(1+seq)/3*dpds                  ! GPa
 dpdt_eq=dot_product(w,dpdt(j:j+3))
 alv=dpdt_eq/bt                      ! 1/K
 f_eq=dot_product(w,free(j:j+3)); e_eq=dot_product(w,etot(j:j+3))
 cv_eq=dot_product(w,cvv(j:j+3)); eq_grun=dot_product(w,grv(j:j+3))
 veq=(1+seq)**3
 if (temprk.lt.8 .or. vol0.eq.0) vol0=veq
 eq_vol=veq/vol0
 eq_free=f_eq; eq_etot=e_eq
 eq_cv=cv_eq*k_b*ev2kjpermol*1000                    ! J/K/mol
 eq_cp=eq_cv+v0*veq*temprk*bt*1d9*alv*alv*ev2kjpermol*1000
 entropy=(eq_etot-eq_free)/temprk*1000               ! J/K/mol
 cvte=alv
 boft=bt
 do i=1,ng
    write(uio,3)temprk,strain(i),free(i),pres(i),dble(i0),b0,bloc(i),grv(i),dpdt(i), &
 &              fst(i),pst(i),fag(i),pag(i),free(i)-fst(i)-fag(i),pres(i)-pst(i)-pag(i)
 enddo
 write(ulog,'(a,f9.2,a,f12.8,a,f12.4,a,es12.4,a,f10.5)')' QHA: T=',temprk,' strain=',seq,' B_T=',bt, &
 &     ' alpha_V=',alv,' <gamma>=',eq_grun
3 format(99(2x,g11.4))
 end subroutine QHA_free_energy_vs_strain
!------------------------------------------------------------
 subroutine lagrange4(xs,x,w,dw)
!! weights w (value) and dw (first derivative) at x of the cubic through the 4 points xs
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: xs(4),x
 real(r15), intent(out) :: w(4),dw(4)
 integer i,j,k
 real(r15) p,s
 do i=1,4
    p=1; s=0
    do j=1,4
       if (j.eq.i) cycle
       p=p*(x-xs(j))/(xs(i)-xs(j))
    enddo
    w(i)=p
    do k=1,4                          ! derivative: sum over the left-out factor
       if (k.eq.i) cycle
       p=1/(xs(i)-xs(k))
       do j=1,4
          if (j.eq.i .or. j.eq.k) cycle
          p=p*(x-xs(j))/(xs(i)-xs(j))
       enddo
       s=s+p
    enddo
    dw(i)=s
 enddo
 end subroutine lagrange4
!=======================================================
 subroutine check_mdyn(ndyn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 use kpoints
 use constants, only : r15
 implicit none
 integer ndyn,nk,j,k,l,i,ier,i1,j1,k1,mi,inside
 real(r15) kp(3,nk),eival(ndyn,nk),eivl(ndyn)
 complex(r15) eivec(ndyn,ndyn,nk),dynm(ndyn,ndyn),eivc(ndyn,ndyn),d2(ndyn,ndyn)

 do i=1,nk
    dynm=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        dynm(j,k)=dynm(j,k)+eival(l,i)*eivec(j,l,i)*conjg(eivec(k,l,i))
      enddo
    enddo
    enddo
    call get_k_info(-kp(:,i),NC,mi,i1,j1,k1,inside)
    d2=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        d2(j,k)=d2(j,k)+eival(l,mi)*eivec(j,l,mi)*conjg(eivec(k,l,mi))
      enddo
    enddo
    enddo

    do j=1,ndyn
    do k=1,ndyn
       if(abs( d2(j,k)-conjg(dynm(j,k))) .gt.1d-5) then
         write(*,4)'CHECK_MDYN: j,k,d(-q),d(q)=',j,k,d2(j,k),dynm(j,k)
       endif
    enddo
    enddo
!   return

    call diagonalize(ndyn,dynm,eivl,ndyn,eivc,ier)
    do j=1,ndyn
       if(abs(eivl(j)-eival(j,i)).gt.1d-4) then
          write(*,3)'CHECK_MDYN:j,eivl,eival=',j,eivl(j),eival(j,i)
       endif
       do k=1,ndyn
         if(abs(eivc(k,j)-eivec(k,j,i)).gt.1d-4) then
         if(abs(eivc(k,j)+eivec(k,j,i)).gt.1d-4) then
            write(*,4)'CHECK_MDYN:j,k,eiv(j),eivecs(k,j)=',j,k,eivl(j),eivc(k,j),eivec(k,j,i)
         endif
         endif
       enddo
    enddo
 enddo
3 format(a,i6,9(1x,f15.6))
4 format(a,2i6,f15.7,9(1x,f11.4))

 end subroutine check_mdyn
!===========================================================
 subroutine write_eivecs(ndyn,nk,kp,eigenvec)
 use lattice
 use ios
 use geometry
 use constants
 use kpoints ! for shift
 implicit none
 integer nk,ndyn,i,la,j,i2,j2,k2,inside,l,i3,j3,k3
 real(r15) kp(3,nk),w(3)
 complex(r15) eigenvec(ndyn,ndyn,nk)

 write(ulog,*)'WRITE_EIGVECS: eivc(q)-conjg(eivc(-q)),eivc(q),eivc(-q)'
do i=1,nk
! bring kp(j) into the FBZ

   call get_k_info(kp(:,i)-shift,NC,l,i3,j3,k3,inside)
   if (i.ne.l) then
      write(ulog,*) 'i .ne. l ',i,l
      stop
   endif
!  call bring_to_cell_c(-kp(:,i),g1,g2,g3,r1,r2,r3,w)
!  w =w /2/pi
   w = -kp(:,i)+shift
   call get_k_info(w,NC,j,i2,j2,k2,inside)
!  write(ulog,5)'k,-k=',i3,j3,k3,' ** ',i2,j2,k2,kp(:,i),w

  if (j.lt.nk) then
   do la=1,ndyn
   do l=1,ndyn
      if(abs(eigenvec(l,la,i)-conjg(eigenvec(l,la,j))).gt.1d-4) then
      if(abs(eigenvec(l,la,i)+conjg(eigenvec(l,la,j))).gt.1d-4) then
          write(ulog,8)'la,l,k,-k=', la,l,i,j,eigenvec(l,la,i)-conjg(eigenvec(l,la,j)) &
      &     ,eigenvec(l,la,i),eigenvec(l,la,j)
      endif
      endif
   enddo
   enddo
  else
    write(ulog,4)'kpoint i in FBZ=',i,kp(:,i)
    write(ulog,4)'was taken to -k=',j,w
  endif
enddo
4 format(a,i5,2x,99(1x,f8.3))
5 format(a,3i5,a,3i5,2x,99(1x,f8.3))
8 format(a,4i6,99(1x,f8.3))

 end subroutine write_eivecs
!==========================================================
 subroutine enforce_asr_on_phi
! enforces asr on force constants instead of the dynamical matrix
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15
 implicit none
 real(r15) mi,mj,nrm1,rr(3),asr,junk,oldfc
 integer i0,j,j0,al,be,i3,i4,j3,t,ti,ired,g,iredsave

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 do i0=1,natom_prim_cell
!   write(*,3)'i0,ri0=',i0,atom0(i0)%equilibrium_pos ,atompos(:,i0)
 do al=1,3
 do be=1,3
    i3 = al+3*(i0-1)
    i4 = be+3*(i0-1) ! used for asr
    mi = atom0(i0)%mass
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    asr = 0 
    gloop: do g=1,map(2)%ngr
    tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
       if(map(2)%keep(counteri(2,g,ti)).ne.1) cycle
       ired = map(1)%nkeptind + current(2,g,ti)    ! this is the corresponding index of i in ared
       tloop: do t=1,map(2)%nt(g)
          if ( i0 .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t)  &
&             .or. be .ne. map(2)%gr(g)%ixyz(2,t) ) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          j0 = iatomcell0(j)
! this produces the diagonal term fc(i0,al;i0,be)
          if ( i0 .eq. j ) iredsave=ired
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)
          rr = atompos(:,j) - atompos(:,j0)
          junk = fcs(ired) * map(2)%gr(g)%mat(t,ti) 
          asr=asr+junk
       enddo tloop
    enddo tiloop
    enddo gloop
    oldfc=fcs(iredsave)
    fcs(iredsave)=fcs(iredsave)-asr
    if(abs(asr).gt.1d-6) then
       write(*   ,5)'i0,al,be,asr_corr,oldfc,newfc=',i0,al,be,asr,oldfc,fcs(iredsave)
       write(ulog,5)'i0,al,be,asr_corr,oldfc,newfc=',i0,al,be,asr,oldfc,fcs(iredsave)
    endif
 enddo
 enddo
 enddo
5 format(a,3i4,9(1x,g13.6))
 end subroutine enforce_asr_on_phi
!===========================================================
 subroutine enforce_asr_simple(n,dyn)
! this dyn(q=0) has the mass denominator
 use atoms_force_constants , only :atom0
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 integer nat,al,be,io,jo,lo,tau,taup
 complex(r15), intent(inout) :: dyn(n,n)
 complex(r15) sumd

 nat=n/3
 do tau=1,nat
 do al=1,3
    io=al+3*(tau-1)
    do be=1,3
       lo=be+3*(tau-1)
       sumd=0
       do taup=1,nat
       !  if (iat.eq.jat) cycle
          jo=be+3*(taup-1)
          sumd=sumd+dyn(io,jo) * sqrt(atom0(tau)%mass*atom0(taup)%mass)
       enddo
       write(*,4)'enforce_asr: old, new d(i,i)=',io,lo,dyn(io,lo),sumd
       dyn(io,lo)=dyn(io,lo) - sumd / atom0(tau)%mass
    enddo
 enddo
 enddo
4 format(a,2i4,9(1x,g11.4))

 end subroutine enforce_asr_simple
!===========================================================
 subroutine check_asr_dynmat(n,dyn,sumd,imass,ier)
 use constants, only : r15
 use atoms_force_constants , only :atom0
 use ios, only : ulog
 implicit none
 integer, intent(in) :: n,imass
 integer, intent(out) :: ier
 complex(r15), intent(inout) :: dyn(n,n),sumd
 integer nat,al,be,io,jo,tau,taup
 real(r15) maxdyn,sqm

 maxdyn=maxval(abs(dyn))
 nat=n/3
 ier=0
 do tau=1,nat
 do al=1,3
    io=al+3*(tau-1)
    do be=1,3
!      write(*,*)'i,al,be=',tau,al,be
       sumd=cmplx(0d0,0d0)
       do taup=1,nat
          jo=be+3*(taup-1)
          if(imass.eq.0) then  ! dyn was not divided by mass term
             sqm = 1
          else
             sqm = sqrt(atom0(tau)%mass*atom0(taup)%mass)
          endif
          sumd=sumd+dyn(io,jo) * sqm
       enddo
       if(abs(sumd).gt.1d-4*maxdyn*sqm) then
          ier=1
          call warn2(ulog,"check_asr_dynmat: asr violation")
          call warn2(6   ,"check_asr_dynmat: asr violation")
          write(   *,5)'ASR_CHECK: tau,al,be,sum,norm=',tau,al,be,sumd,maxdyn*sqm
          write(ulog,5)'ASR_CHECK: tau,al,be,sum,norm=',tau,al,be,sumd,maxdyn*sqm
!         call enforce_asr_simple(n,dyn)
       else
          write(*,5)'ASR_CHECK PASSED: i0,ial,be,sum=',tau,al,be,sumd
       endif
    enddo
 enddo
 enddo

5 format(a,3i3,9(1x,f10.6))

 end subroutine check_asr_dynmat
!===========================================================
subroutine band_sort_bs(nkp,ndyn,kp,eival,eivec,vel)
!! orders the bands along the path by eigenvector continuity (2026-10-22; replaces the greedy
!! three-step version, which could give one eigenvalue to two bands and drop another):
!! - inside each degenerate set at k (|w_i-w_j| < degtol) the eigenvectors are first rotated to
!!   match the previous point's vectors best (polar decomposition of the overlap block), so the
!!   arbitrary basis of a degenerate subspace does not decide the connection;
!! - bands at k-1 and k are then matched by the assignment with the largest total overlap
!!   |<e_b(k-1)|e_m(k)>|^2 (Hungarian algorithm): always a permutation;
!! - a connection across more than 5 x the average level spacing (bandwidth/ndyn, as in the old
!!   version) is penalized; a small energy term breaks ties between equal overlaps. Without the
!!   window, BZO connected 295 to 607 1/cm on the first step away from Gamma (LO-TO);
!! - at a jump in the path (|dk| > 4 x the median step) the order restarts from the energies.
!! Only the order changes: the eigenvalues, eigenvectors and velocities written are the ones of
!! get_frequencies; the rotated vectors are used only as the reference for the next point.
 use constants, only : r15,cnst
 use ios, only : ulog
implicit none
integer, intent(in) :: nkp,ndyn
real(r15), intent(in) :: kp(3,nkp)
real(r15), intent(inout) :: eival(ndyn,nkp) ,vel(3,ndyn,nkp)
complex(r15), intent(inout) :: eivec(ndyn,ndyn,nkp)
real(r15), parameter :: degtol=1d-3   ! 1/cm, same as THERMACOND
integer b,k,m,g,i0,i1,ier,ng,nweak
integer, allocatable :: emap(:,:),srt(:),top(:),perm(:)
real(r15), allocatable :: w(:),cost(:,:),dk(:),wt(:),lam(:),dks(:),wref(:)
complex(r15), allocatable :: ref(:,:),e(:,:),a(:,:),mm(:,:),h(:,:),q(:,:),u(:,:)
real(r15) dkmed,c,window

 allocate(emap(ndyn,nkp),srt(ndyn),perm(ndyn),w(ndyn),cost(ndyn,ndyn),dk(nkp),dks(nkp),wref(ndyn))
 allocate(ref(ndyn,ndyn),e(ndyn,ndyn))
 dk=0
 do k=2,nkp
    dk(k)=sqrt(sum((kp(:,k)-kp(:,k-1))**2))
 enddo
 dks=dk; call sort_real(nkp-1,dks(2:nkp))
 dkmed=dks(1+nkp/2)

 do b=1,ndyn
    emap(b,1)=b
 enddo
 ref=eivec(:,:,1)
 wref=sign(sqrt(abs(eival(:,1))),eival(:,1))*cnst
 window=5*(maxval(sign(sqrt(abs(eival)),eival))-minval(sign(sqrt(abs(eival)),eival)))*cnst/ndyn
 write(ulog,*)'BAND_SORT_BS: connections across more than ',window,' 1/cm are penalized'
 nweak=0

 kloop: do k=2,nkp
    if (dk(k).gt.4*dkmed) then   ! a jump in the path: start again from the energy order
       do b=1,ndyn
          emap(b,k)=b
       enddo
       ref=eivec(:,:,k)
       wref=sign(sqrt(abs(eival(:,k))),eival(:,k))*cnst
       cycle kloop
    endif
    e=eivec(:,:,k)
    w=sign(sqrt(abs(eival(:,k))),eival(:,k))*cnst
    call sort_index(ndyn,w,srt)
! rotate inside each degenerate set
    i0=1
    do while (i0.le.ndyn)
       i1=i0
       do while (i1.lt.ndyn)
          if (abs(w(srt(i1+1))-w(srt(i0))).ge.degtol) exit
          i1=i1+1
       enddo
       ng=i1-i0+1
       if (ng.gt.1) call align_block(srt(i0:i1),ng)
       i0=i1+1
    enddo
! optimal assignment of the reference bands b to the bands m at k
    do b=1,ndyn
    do m=1,ndyn
       cost(b,m)=1-abs(dot_product(ref(:,b),e(:,m)))**2 + 1d-3*abs(wref(b)-w(m))/window
       if (abs(wref(b)-w(m)).gt.window) cost(b,m)=cost(b,m)+10
    enddo
    enddo
    call hungarian_min(ndyn,cost,perm)
    do b=1,ndyn
       c=abs(dot_product(ref(:,b),e(:,perm(b))))**2        ! overlap of the chosen connection
       emap(b,k)=perm(b)
       ref(:,b)=e(:,perm(b))
       wref(b)=w(perm(b))
       if (c.lt.0.3) then
          nweak=nweak+1
          write(ulog,'(a,i6,3(1x,f9.4),a,i4,a,f6.3)')' BAND_SORT_BS: weak connection at k',k,kp(:,k), &
 &              ' band',b,' overlap',c
       endif
    enddo
 enddo kloop
 if (nweak.gt.0) write(*,*)'BAND_SORT_BS: ',nweak,' connections with overlap < 0.3; see log'

 do k=1,nkp
    eival(:,k)=eival(emap(:,k),k)
    eivec(:,:,k)=eivec(:,emap(:,k),k)
    vel(:,:,k)=vel(:,emap(:,k),k)
 enddo
 deallocate(emap,srt,perm,w,cost,dk,dks,ref,e,wref)

 contains

 subroutine align_block(idx,n)
!! rotates the degenerate eigenvectors e(:,idx) by U (unitary) to maximize Re Tr(M U), M = the
!! overlap of the n reference vectors closest to the subspace with it: U = M^+ (M M^+)^(-1/2)
 integer, intent(in) :: n,idx(n)
 integer j,l,mx(1)
 allocate(a(ndyn,n),wt(ndyn),top(n),mm(n,n),h(n,n),q(n,n),u(n,n),lam(n))
 a=matmul(conjg(transpose(ref)),e(:,idx))       ! a(b,j) = <ref_b|e_j>
 do j=1,ndyn
    wt(j)=sum(abs(a(j,:))**2)
 enddo
 do j=1,n
    mx=maxloc(wt); top(j)=mx(1); wt(mx(1))=-1
 enddo
 mm=a(top,:)
 h=matmul(mm,conjg(transpose(mm)))
 call diagonalize(n,h,lam,n,q,ier)
 if (ier.eq.0 .and. minval(lam).gt.1d-8) then
    do l=1,n
       q(:,l)=q(:,l)/sqrt(sqrt(lam(l)))
    enddo
    u=matmul(conjg(transpose(mm)),matmul(q,conjg(transpose(q))))
    e(:,idx)=matmul(e(:,idx),u)
 endif                                           ! else (subspace nearly orthogonal): keep the basis
 deallocate(a,wt,top,mm,h,q,u,lam)
 end subroutine align_block

end subroutine band_sort_bs
!===========================================================
 subroutine hungarian_min(n,a,assign)
!! assignment minimizing sum_i a(i,assign(i)) (Kuhn-Munkres with potentials, O(n^3))
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: a(n,n)
 integer, intent(out) :: assign(n)
 real(r15) u(0:n),v(0:n),minv(0:n),cur,delta
 integer p(0:n),way(0:n),i,j,i0,j0,j1
 logical used(0:n)
 u=0; v=0; p=0; way=0
 do i=1,n
    p(0)=i; j0=0; minv=huge(1d0); used=.false.
    do
       used(j0)=.true.; i0=p(j0); delta=huge(1d0); j1=0
       do j=1,n
          if (.not.used(j)) then
             cur=a(i0,j)-u(i0)-v(j)
             if (cur.lt.minv(j)) then
                minv(j)=cur; way(j)=j0
             endif
             if (minv(j).lt.delta) then
                delta=minv(j); j1=j
             endif
          endif
       enddo
       do j=0,n
          if (used(j)) then
             u(p(j))=u(p(j))+delta; v(j)=v(j)-delta
          else
             minv(j)=minv(j)-delta
          endif
       enddo
       j0=j1
       if (p(j0).eq.0) exit
    enddo
    do
       j1=way(j0); p(j0)=p(j1); j0=j1
       if (j0.eq.0) exit
    enddo
 enddo
 do j=1,n
    assign(p(j))=j
 enddo
 end subroutine hungarian_min
!===========================================================
 subroutine sort_index(n,x,idx)
!! idx = permutation that sorts x in ascending order (insertion sort; n is small)
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: x(n)
 integer, intent(out) :: idx(n)
 integer i,j,t
 do i=1,n
    idx(i)=i
 enddo
 do i=2,n
    t=idx(i); j=i-1
    do while (j.ge.1)
       if (x(idx(j)).le.x(t)) exit
       idx(j+1)=idx(j); j=j-1
    enddo
    idx(j+1)=t
 enddo
 end subroutine sort_index
!===========================================================
 subroutine sort_real(n,x)
!! sorts x in ascending order in place (insertion sort)
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(inout) :: x(n)
 integer i,j
 real(r15) t
 do i=2,n
    t=x(i); j=i-1
    do while (j.ge.1)
       if (x(j).le.t) exit
       x(j+1)=x(j); j=j-1
    enddo
    x(j+1)=t
 enddo
 end subroutine sort_real
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 use constants, only : r15
 implicit none
 integer, intent(in) :: n,nv
 integer, intent(out) :: ier
 complex(r15), intent(in) :: mat(n,n)
 complex(r15), intent(out) :: eivec(n,n)
 real(r15), intent(out) :: eival(n)
! This is used by eigch
! real(r15), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(r15), allocatable :: rwork(:)
  complex(r15), allocatable :: work(:)
  integer lwork  !, zero


! n = size(mat(:,1))
 if (n .ne. size(eival) ) then
    write(   *,*)' EIGCH, size inconsistency:mat,eival=',n, size(eival)
    stop
 endif

! tol = -1.d0
! zero= 0  ! for now do not compute any eigenvectors
! ier=0
! allocate(w(n,7),lw(n))
! call eigch(mat,n,n,n,nv,tol,w,lw,eival,eivec,ier)
! deallocate(w,lw)
! return

  lwork = max(1,2*n-1)
  allocate(work(lwork),rwork(max(1,3*n-2)))
  if(nv.eq.0) then
    call zheev('N','U',n,mat,n,eival,work,lwork,rwork,ier)
  else
    call zheev('V','U',n,mat,n,eival,work,lwork,rwork,ier)
    eivec = mat
  endif
  deallocate(work,rwork)
4 format(i5,1x,3(f6.3,1x,f6.3,4x))

 end subroutine diagonalize
!==========================================================
 subroutine calculate_dos(mx,eival,wkp,mesh,ene,ds)
! use gaussian broadening
 use ios
 use params
 use om_dos, only : width
 use constants
 implicit none
 integer, intent(in):: mx,mesh
 real(r15),intent(in) :: wkp(mx),ene(mesh),eival(mx)
 real(r15),intent(inout):: ds(mesh)
 integer i,j
 real(r15) x,delta

! wkp=1d0/(nkx*nky*nkz)
! write(udos,2)'# wkp,width =',wkp,width
  if( width.lt.(ene(2)-ene(1))/2 )  width=(ene(2)-ene(1))/2

    do i=1,mesh
       ds(i) = 0
       do j=1,mx
          x = (eival(j) - ene(i))/width   ! these are all in cm^-1
!         x = (eival(j) - ene*ene)/width/width
          if ( abs(x) .gt. 5 ) cycle
          ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx
       enddo
    enddo

2 format(a,999(1x,g11.4))
3 format(i5,9(3x,g11.4))

 end subroutine calculate_dos
!============================================================
 subroutine dielectric(n,eival,eivec,mesh,om,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natom_prim_cell, atom0
 use lattice, only: volume_r0
 use params, only : coef
 use born, only : epsil
 use om_dos, only : width
 use constants, only : r15,pi,cnst,ci
 use kpoints, only : normal
 use eigen , only : mysqrt
 implicit none
 integer, intent(in):: n,mesh
 real(r15),intent(in) :: om(mesh),eival(n)
 complex(r15),intent(in) :: eivec(n,n)
 real(r15),intent(out) :: epsilon0(3,3)
 complex(r15) chi(3,3),epsnormal !,intent(out) :: chi(3,3,mesh)
 real(r15) d1,d2,z(n,3),reflectivity
 integer i,j,k,b,al,be,ga,de,tau,taup

 if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts!! '
    call write_out(ulog,' imag(eivec) ',aimag(eivec))
    call write_out(   6,'DIELECTRIC:q=0 imag(eivec) ',aimag(eivec))
!   stop
 endif

! calculate Born chargs in the eivec basis 
 write(ulog,*)'#================== IR intensities versus  frequency(1/cm) =================='
 z=0
 do b=4,n
    do de=1,3
       do ga=1,3
       do tau =1,natom_prim_cell
          j=ga+3*(tau-1)
          z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
       enddo
       enddo
    enddo
    write(ulog,3) mysqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)),z(b,:)
 enddo
 write(ulog,*)'#================ w_la , I_la , Z_la ======================================='

! calculate and write susceptibility on a frequency mesh
 open(345,file='chi_real.dat')
 open(346,file='chi_imag.dat')
 write(345,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , eps_nn, R_n' 
 write(346,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz, n(w), k(w)' 
 do i=1,mesh
    chi=0
    do al=1,3
    do be=1,3
    do b=4,n
       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &
 &      (-om(i)*om(i)+eival(b)*cnst*cnst-ci*om(i)*width)
    enddo
    enddo
    enddo
    chi=chi* cnst*cnst * coef
    if(i.eq.1) epsilon0=epsil+real(chi)
    epsnormal=dot_product(normal,matmul((chi+epsil),normal))
    reflectivity= abs((sqrt(epsnormal)-1)/(sqrt(epsnormal)+1))**2
    write(345,2)om(i),real(chi),epsnormal,reflectivity 
    write(346,2)om(i),aimag(chi),real(sqrt(epsnormal)) ,aimag(sqrt(epsnormal)) 
 enddo


3 format(f10.4,1x,g11.4,3x,3(1x,g11.4))
2 format(99(1x,g11.4))

 end subroutine dielectric
!===============
 function parlinski_sf(q)  result(sf)
!! with the soft phase convention
 use constants, only : r15
 use lattice, only :  volume_r0,g0ws26
 use geometry, only: length,svector
 implicit none
 integer i
 real(r15), intent(in):: q(3)  
 real(r15) damp,s,dsf(3),expo
 type(svector) sf

   damp=0.1  *volume_r0**0.666666
   sf%s= exp(-2*length(q)**2*damp)
   sf%v= -4*q(:)*damp*sf%s  ! uses spherically symmetry to preserve symmetry
! these added terms donot break symmetry and have G-vector periodicity!
   do i=1,26  
      expo=exp(-2*length(q+g0ws26(:,i))**2*damp)
      sf%s=sf%s + expo
      sf%v=sf%v - expo*4*(q+g0ws26(:,i))*damp
   enddo


 end function parlinski_sf
!===============
 subroutine remove_eiktau(n,q,dyn,dyn0)
! divides the dynamical matrix by the phase term e^ik.(tau-taup)
 use constants, only : r15,ci
 use params, only : tolerance
 use atoms_force_constants, only : natom_prim_cell, atompos
 implicit none
 integer , intent(in) :: n
 real(r15), intent(in):: q(3)  
 complex(r15), intent(in):: dyn(n,n)   
 complex(r15), intent(out):: dyn0(n,n)   
 real(r15) rr(3),dmax 
 integer al,be,tau,taup,i,j,err

 dmax=maxval(abs(dyn))
 err=0
 do al=1,3
 do tau =1,natom_prim_cell
    i=al+3*(tau-1)
 do be=1,3
 do taup=1,natom_prim_cell
    j=be+3*(taup-1)
    
    rr = atompos(:,taup) - atompos(:,tau)
    dyn0(i,j)=dyn(i,j)*exp(-ci*dot_product(q,rr)) !* sqrt(atom0(tau)%mass*atom0(taup)%mass)
 enddo
 enddo
 enddo
 enddo

5 format(a,2i5,4(1x,f10.5))

 end subroutine remove_eiktau
!===============
 subroutine putback_eiktau(n,q,dyn,dyn0)
 use constants, only : r15,ci
 use params, only : tolerance
 use atoms_force_constants, only : natom_prim_cell, atompos
 implicit none
 integer , intent(in) :: n
 real(r15), intent(in):: q(3)  
 complex(r15), intent(out):: dyn(n,n)   
 complex(r15), intent(in):: dyn0(n,n)   
 real(r15) rr(3),dmax 
 integer al,be,tau,taup,i,j

 dmax=maxval(abs(dyn))
 do al=1,3
 do tau =1,natom_prim_cell
    i=al+3*(tau-1)
 do be=1,3
 do taup=1,natom_prim_cell
    j=be+3*(taup-1)
    
    rr = atompos(:,taup) - atompos(:,tau)
    dyn(i,j)=dyn0(i,j)*exp(ci*dot_product(q,rr)) !/ sqrt(atom0(tau)%mass*atom0(taup)%mass)
 enddo
 enddo
 enddo
 enddo

 end subroutine putback_eiktau
!============================================================
 subroutine elastic_born_huang(araw,s,c)
!! elastic constants c (GPa) for the Lagrangian strain from the unsymmetrized
!! araw(a,b,g,d) = A_{ab,gd} = -1/(2V) sum Phi^{ag} d_b d_d (a,g displacements), i.e. the brackets
!! [ab,gd] = A_{ag,bd} (Born and Huang), and the reference stress s (GPa):
!!   c_abgd = Sym{ [ag,bd] + [bg,ad] - [ab,gd] + (d_ab s_gd + d_gd s_ab)/2 - d_ag s_bd - d_bd s_ag }
!! Sym = average over a<->b, g<->d and (ab)<->(gd). Checked on the exact model (REFERENCE/exact_model/bprime/
!! elastic_check.py): exact (5e-10 GPa) for relaxed ZnO, CaF2 and stressed CaF2; with residual forces
!! (stressed ZnO, HfO2 references) 0.2 GPa off, the part FOCEX handles through the relaxation.
!! Until 2026-11-01 FOCEX used Sym(A) - Sym(d_ag s_bd): off-diagonal and shear constants wrong
!! (relaxed ZnO model: C12 171 instead of 221 GPa, C44 108 instead of 78).
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: araw(3,3,3,3),s(3,3)
 real(r15), intent(out) :: c(3,3,3,3)
 real(r15) t(3,3,3,3),dl(3,3)
 integer a,b,g,d
 dl=0; dl(1,1)=1; dl(2,2)=1; dl(3,3)=1
 do a=1,3
 do b=1,3
 do g=1,3
 do d=1,3
    t(a,b,g,d)=araw(a,b,g,d)+araw(b,a,g,d)-araw(a,g,b,d) &
 &            +(dl(a,b)*s(g,d)+dl(g,d)*s(a,b))/2 - dl(a,g)*s(b,d) - dl(b,d)*s(a,g)
 enddo
 enddo
 enddo
 enddo
 do a=1,3
 do b=1,3
 do g=1,3
 do d=1,3
    c(a,b,g,d)=(t(a,b,g,d)+t(b,a,g,d)+t(a,b,d,g)+t(b,a,d,g))/4
 enddo
 enddo
 enddo
 enddo
! and over (ab)<->(gd): with a non-hydrostatic stress the expression above is not symmetric under it, a
! second derivative is (2026-11-02, exact model ZnO at 14.97/14.97/14.87 GPa, internal coordinates relaxed:
! H13 51.452, H31 51.420, exact 51.436 eV per cell = the average; hydrostatic stress: no change)
 c=(c+reshape(c,[3,3,3,3],order=[3,4,1,2]))/2
 end subroutine elastic_born_huang
!============================================================
 subroutine pressure_derivative_b(uio)
!! B and B' = dB/dP of the reference at T=0 from FC2 and FC3 (FOCEX/notes/generalized_huang.tex Part II;
!! method checked on the exact model, REFERENCE/exact_model/bprime/w3_check.py and mkrelaxed.py):
!!   eta1 = -Ht^-1 v, u1 = xi.eta1: strain and internal relaxation per unit pressure (path p1)
!!   the crystal deformed by +-lam p1 has FC2 Phi +- lam dPhi, dPhi_{tau j} = sum_k Psi_{tau j k}(S_k-S_tau),
!!   S_j = eta1.x_j + u1(tau_j), relative vectors d +- lam (eta1.d + u1_j - u1_tau), pressure P0 +- lam;
!!   B(+-lam) = V/(v^T Ht^-1 v), Ht = V (c + S2) + P V V2, with c its relaxed elastic constants (Born-Huang
!!   brackets + stress terms - QGQ, as elastic_born_huang), S2 = d^2(sigma:eta^2/2), V2 = (Tr e)^2 - Tr e^2;
!!   B' = (B(+lam) - B(-lam))/(2 lam)  (FC4 does not enter B').
!! Evaluating p1^T H[dPhi] p1 with the plain sum -1/2 sum dPhi d d instead gives a wrong W3 (exact ZnO model:
!! B' -0.064 instead of -0.157), the same reason as Sym(A) for the elastic constants.
!! Assumes a hydrostatic reference stress -P0 (sigma0) and no residual forces; deviations are reported.
 use ios, only : ulog
 use lattice, only : volume_r0
 use atoms_force_constants
 use svd_stuff
 use constants, only : r15
 use eigen, only : ndyn
 use mech
 implicit none
 integer, intent(in) :: uio
 real(r15), parameter :: gpa=160.21766208d0, lam=1d-4     ! 1 eV/A^3 in GPa; pressure step (eV/A^3)
 real(r15) cs(3,3,3,3),e1(6),eta(3,3),p0,v0,b0,bpm(2),sdev,psi,dk(3),sg
 real(r15), allocatable :: u1(:,:),ph(:,:,:),dph(:,:,:),dv(:,:)
 integer, allocatable :: slot(:,:),pt(:),pj(:)
 integer i,j,k,g,t,ti,tau,al,be,ga,ired,nat,np,is

 nat=natom_prim_cell
 v0=volume_r0
 p0=-(sigma0(1,1)+sigma0(2,2)+sigma0(3,3))/3/gpa          ! eV/A^3
 sdev=0
 do i=1,3
    do j=1,3
       sdev=max(sdev,abs(sigma0(i,j)+merge(p0*gpa,0d0,i.eq.j)))
    enddo
 enddo
! pairs (tau, j) of FC2 and FC3 terms, tau in the primitive cell
 allocate(slot(nat,size(atompos,2)))
 slot=0; np=0
 call collect_pairs(2)
 call collect_pairs(3)
 allocate(pt(np),pj(np),ph(3,3,np),dph(3,3,np),dv(3,np),u1(3,nat))
 slot=0; np=0
 call collect_pairs(2)
 call collect_pairs(3)
 ph=0; dph=0
 do i=1,np
    dv(:,i)=atompos(:,pj(i))-atompos(:,pt(i))
 enddo
 do g=1,map(2)%ngr
 do ti=1,map(2)%ntind(g)
    if(map(2)%keep(counteri(2,g,ti)).ne.1) cycle
    ired=map(1)%nkeptind+current(2,g,ti)
    do t=1,map(2)%nt(g)
       tau=map(2)%gr(g)%iat(1,t); if (tau.gt.nat) cycle
       is=slot(tau,map(2)%gr(g)%iat(2,t))
       al=map(2)%gr(g)%ixyz(1,t); be=map(2)%gr(g)%ixyz(2,t)
       ph(al,be,is)=ph(al,be,is)+fcs(ired)*map(2)%gr(g)%mat(t,ti)
    enddo
 enddo
 enddo
! B and eta1 of the reference
 call b_of_state(0d0,b0,e1)
 eta(1,1)=e1(1); eta(2,2)=e1(2); eta(3,3)=e1(3)
 eta(2,3)=e1(4); eta(3,2)=e1(4); eta(1,3)=e1(5); eta(3,1)=e1(5); eta(1,2)=e1(6); eta(2,1)=e1(6)
 u1=0
 if (ndyn.gt.3) then
    do tau=1,nat
    do al=1,3
       u1(al,tau)=sum(xi(al+3*(tau-1),:,:)*eta)
    enddo
    enddo
 endif
! dPhi from the FC3 terms
 do g=1,map(3)%ngr
 do t=1,map(3)%nt(g)
    tau=map(3)%gr(g)%iat(1,t); if (tau.gt.nat) cycle
    j=map(3)%gr(g)%iat(2,t); k=map(3)%gr(g)%iat(3,t)
    al=map(3)%gr(g)%ixyz(1,t); be=map(3)%gr(g)%ixyz(2,t); ga=map(3)%gr(g)%ixyz(3,t)
    dk=matmul(eta,atompos(:,k)-atompos(:,tau))+u1(:,iatomcell0(k))-u1(:,tau)
    is=slot(tau,j)
    do ti=1,map(3)%ntind(g)
       if (map(3)%gr(g)%mat(t,ti).eq.0) cycle
       ired=counteri(3,g,ti) + map(1)%ntotind + map(2)%nkeptind
       psi=fcs(ired)*map(3)%gr(g)%mat(t,ti)
       dph(al,be,is)=dph(al,be,is)+psi*dk(ga)
    enddo
 enddo
 enddo
 do i=1,2
    sg=merge(1d0,-1d0,i.eq.1)
    call b_of_state(sg*lam,bpm(i),e1)
 enddo
 b0p=(bpm(1)-bpm(2))/(2*lam)
 write(ulog,'(a,f14.6,a,f12.6)')' PRESSURE_DERIVATIVE_B: B (GPa) =',b0*gpa,'   B0''=dB/dP =',b0p
 write(ulog,'(a,3es12.4,a,es12.4)')'   eta1 diagonal per GPa =',eta(1,1)/gpa,eta(2,2)/gpa,eta(3,3)/gpa, &
 &   '  max|u1| per GPa (A)=',maxval(abs(u1))/gpa
 write(uio ,'(a,f14.6,a,f12.6)')' B (GPa) from the relaxed elastic constants =',b0*gpa,'   B0''=dB/dP from FC2+FC3 =',b0p
 if (sdev.gt.1d-3) write(ulog,*)' PRESSURE_DERIVATIVE_B: WARNING: reference stress not hydrostatic, max dev (GPa)=',sdev
 if (ndyn.gt.3) then
    if (maxval(abs(pi0)).gt.1d-4) write(ulog,*)' PRESSURE_DERIVATIVE_B: WARNING: residual forces (eV/A)=', &
 &        maxval(abs(pi0)),'; the formula assumes none'
 endif
 deallocate(slot,pt,pj,ph,dph,dv,u1)

 contains

 subroutine collect_pairs(rank)
 integer, intent(in) :: rank
 integer gg,tt,ta,jj
 do gg=1,map(rank)%ngr
 do tt=1,map(rank)%nt(gg)
    ta=map(rank)%gr(gg)%iat(1,tt); if (ta.gt.nat) cycle
    jj=map(rank)%gr(gg)%iat(2,tt)
    if (slot(ta,jj).eq.0) then
       np=np+1; slot(ta,jj)=np
       if (allocated(pt)) then
          pt(np)=ta; pj(np)=jj
       endif
    endif
 enddo
 enddo
 end subroutine collect_pairs

 subroutine b_of_state(dp,bb,ee1)
!! B and eta1 of the crystal deformed by dp*p1 (dp = 0: the reference)
 real(r15), intent(in) :: dp
 real(r15), intent(out) :: bb,ee1(6)
 real(r15) a4(3,3,3,3),c4(3,3,3,3),q(3*nat,3,3),phi0(3*nat,3*nat),gm(3*nat-3,3*nat-3),gi(3*nat-3,3*nat-3), &
 &         sig(3,3),vol,dd(3),mc(6,6),mv2(6,6),ht(6,6),hti(6,6),vv(6),fe(3,3),ek(3,3),el(3,3),qs(3*nat,3,3)
 integer ip,a,bq,gq,d,s1,s2,kk,ll,j0
 integer, parameter :: vo(2,6)=reshape((/1,1, 2,2, 3,3, 2,3, 1,3, 1,2/),(/2,6/))
 fe=0; fe(1,1)=1; fe(2,2)=1; fe(3,3)=1
 vol=v0
 if (dp.ne.0) then
    fe=fe+dp*eta
    vol=v0*(fe(1,1)*(fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))-fe(1,2)*(fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1)) &
 &          +fe(1,3)*(fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1)))
 endif
 a4=0; q=0; phi0=0
 do ip=1,np
    dd=dv(:,ip)
    if (dp.ne.0) dd=dd+dp*(matmul(eta,dv(:,ip))+u1(:,iatomcell0(pj(ip)))-u1(:,pt(ip)))
    j0=iatomcell0(pj(ip))
    do a=1,3
    do gq=1,3
       psi=ph(a,gq,ip)+dp*dph(a,gq,ip)
       if (psi.eq.0) cycle
       s1=a+3*(pt(ip)-1); s2=gq+3*(j0-1)
       phi0(s1,s2)=phi0(s1,s2)+psi
       do bq=1,3
          q(s1,gq,bq)=q(s1,gq,bq)+psi*dd(bq)
          do d=1,3
             a4(a,bq,gq,d)=a4(a,bq,gq,d)-0.5d0*psi*dd(bq)*dd(d)
          enddo
       enddo
    enddo
    enddo
 enddo
 a4=a4/vol*gpa
 sig=0; sig(1,1)=-(p0+dp)*gpa; sig(2,2)=sig(1,1); sig(3,3)=sig(1,1)
 call elastic_born_huang(a4,sig,c4)                      ! GPa
 if (nat.gt.1) then                                      ! relaxation of the internal coordinates
    gm=phi0(4:,4:)
    call inverse_real(gm,gi,3*nat-3)
    do s1=1,3*nat
       qs(s1,:,:)=(q(s1,:,:)+transpose(q(s1,:,:)))/2
    enddo
    do a=1,3
    do bq=1,3
    do gq=1,3
    do d=1,3
       c4(a,bq,gq,d)=c4(a,bq,gq,d)-dot_product(qs(4:,a,bq),matmul(gi,qs(4:,gq,d)))/vol*gpa
    enddo
    enddo
    enddo
    enddo
 endif
! Ht in the tensor strain components: V (c + S2) + P V V2, engineering factor 2 on shears for c
 mv2=0
 mv2(1,2)=1; mv2(1,3)=1; mv2(2,3)=1; mv2(2,1)=1; mv2(3,1)=1; mv2(3,2)=1
 mv2(4,4)=-2; mv2(5,5)=-2; mv2(6,6)=-2
 do kk=1,6
 do ll=1,6
    mc(kk,ll)=c4(vo(1,kk),vo(2,kk),vo(1,ll),vo(2,ll))/gpa*merge(2d0,1d0,kk.gt.3)*merge(2d0,1d0,ll.gt.3)
    ek=0; ek(vo(1,kk),vo(2,kk))=1; ek(vo(2,kk),vo(1,kk))=1
    el=0; el(vo(1,ll),vo(2,ll))=1; el(vo(2,ll),vo(1,ll))=1
    mc(kk,ll)=mc(kk,ll)+sum(sig/gpa*(matmul(ek,el)+matmul(el,ek)))/2
 enddo
 enddo
 ht=vol*mc+(p0+dp)*vol*mv2
 call inverse_real(ht,hti,6)
 vv=0; vv(1:3)=vol
 ee1=-matmul(hti,vv)
 bb=vol/dot_product(vv,matmul(hti,vv))
 end subroutine b_of_state

 end subroutine pressure_derivative_b
!===========================================================
 subroutine cubic_elastic 
!! calculates b0p=dB/dP=-d ln B/d ln V=(-1/3B) dB/deta =-(1/27BV0) \sum psi R1 R2 R3 
 use ios , only: ulog, write_out
 use lattice, only : volume_r0 
 use constants, only : ee
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params, only : verbose
 use mech
 use eigen, only : ndyn
 use linalgb
 use params, only : include_fc
 implicit none
 integer tau,taup,taus,al,be,ga,t,g,ti,ired,j,k
 real(r15)  rij(3),rik(3),junk,b0s,rr(3)

 if(include_fc(3).ne.1) then
    write(ulog,*)'CUBIC_ELASTIC: include_fc(3)=1 not selected: ',include_fc(3) 
    write(ulog,*)'CUBIC_ELASTIC:  returning to main...'
    return
 endif
   
 write(  * ,*)' ********** ENTERING cubic elastic *************'
 write(ulog,*)' ********** ENTERING cubic elastic hbulkmod(GPa)=',hbulkmod

 b0p=0;b0s=0
 do al=1,3
 do be=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
 gloop: do g=1,map(3)%ngr
    tloop: do t=1,map(3)%nt(g)
       if ( tau .ne. map(3)%gr(g)%iat(1,t) ) cycle tloop
       if ( al .ne. map(3)%gr(g)%ixyz(1,t) ) cycle tloop
       if ( be .ne. map(3)%gr(g)%ixyz(2,t) ) cycle tloop
       if ( ga .ne. map(3)%gr(g)%ixyz(3,t) ) cycle tloop
       j  = map(3)%gr(g)%iat(2,t)
       k  = map(3)%gr(g)%iat(3,t)
       taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
       rr  = atompos(:,tau)- atompos(:,1)
       rij = atompos(:,j)  - atompos(:,1)
       rik = atompos(:,k)  - atompos(:,1)
    do ti=1,map(3)%ntind(g)  ! index of independent terms in that group g
       ired = map(1)%nkeptind + map(2)%nkeptind + ti
       junk = fcs(ired)* map(3)%gr(g)%mat(t,ti)
       b0s = b0s+junk*atompos(al,tau)*atompos(be,j)*atompos(ga,k)
       b0p = b0p+junk*rr(al)*rij(be)*rik(ga)
    enddo
    enddo tloop
 enddo gloop
 enddo
 enddo
 enddo
 enddo 

 write(ulog,*) 'CUBIC_ELASTIC: before rescaling b0p,b0s=',b0p,b0s

 b0p=b0p *ee/volume_r0*1d30*1d-9 / (-27* bulkmod)
 b0s=b0s *ee/volume_r0*1d30*1d-9 / (-27* bulkmod)
 
 write(ulog,*) 'CUBIC_ELASTIC: b0p,b0s=',b0p,b0s

 end subroutine cubic_elastic
!--------------------------------------------------
 subroutine energy_strain(strain,ene,deriv) 
!! computes the zero-temperature expansion of the potential energy in powers of strain (only if u^0=0)
!! units are the same as E i.e. in eV for both ene and deriv, which is the coefficient of the powers of strain in the Taylor expansion
!! attention: a residual displacement u0 will affect the expansion and is not included 
 use ios , only: ulog, write_out
 use lattice, only : volume_r0 
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use params, only : verbose
 use mech
 use eigen, only : ndyn
 use linalgb
 use params, only : include_fc
 implicit none
 real(r15), intent(in) :: strain 
 real(r15), intent(out) :: ene,deriv(maxrank)
 integer be,ga,t,g,ti,ired,j,ka,cnt2
 integer al(maxrank),tau(maxrank),j_sc(maxrank),res,rnk,rk
 real(r15)  rij(3),rik(3),junk,grad ,ux

  ene =0 ; deriv=0

  deriv(1)=trace(sigma0) *volume_r0/1d30/ee/1d-9  ! to convert back to eV

  res=map(1)%nkeptind ! accumulated # of terms from pervious ranks
  do rnk=2,maxrank   !************* each rank contributes a power of strain
  if(include_fc(rnk).ne.0) then

     grad=0
     do ga=1,natom_prim_cell
     do be=1,3
        tau(1)=ga
        al(1)=be

        cnt2= 0       ! accumulated # of terms from previous groups within the same rank: rnk
        do g=1,map(rnk)%ngr  
           if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
           do t=1,map(rnk)%nt(g) ! sum over all terms in that group
              if ( tau(1).ne. map(rnk)%gr(g)%iat (1,t)   .or.  &
          &        al (1).ne. map(rnk)%gr(g)%ixyz(1,t) ) cycle
              do ti=1,map(rnk)%ntind(g)
! ired is the full index of the indep FC coming in the A*FC=b matrix product
             !   if(rnk.eq.2) then
                    ired = res + current(rnk,g,ti)
             !   else
             !      ired = res + cnt2 + ti  ! same as res+sum(map(rnk)%ntind(1:g-1))+ti 
             !   endif 
!                ux=(atompos(al(1),tau(1))-atompos(al(1),1))*strain
                 ux=(atompos(al(1),tau(1)))*strain
                 do rk=2,rnk
                    j       =  map(rnk)%gr(g)%iat (rk,t)
                    al (rk) =  map(rnk)%gr(g)%ixyz(rk,t)
!                   ux=ux*(atompos(al(rk),j)-atompos(al(rk),1))*strain/dble(rk)
                    ux=ux*(atompos(al(rk),j))*strain/dble(rk)
                 enddo
                 junk=fcs(ired)*map(rnk)%gr(g)%mat(t,ti)
                 ene =ene + junk*ux
                 if (strain.ne.0) grad=grad+junk*ux/strain*rnk
              enddo
           enddo
        enddo

     enddo
     enddo
     if(rnk.ge.2) deriv(rnk-1)=grad

!    if(rnk.eq.2) then
        res = res + map(rnk)%nkeptind
!    else
!       res = res + map(rnk)%ntotind 
!    endif   
!     if(g.gt.1) res=res+cnt2+ map(rnk)%ntind(g-1)
!     write(726,3)'ENERGY_STRAIN: rank, res,energy,grad=',rnk,tau,al,res,ene,grad
  endif
  enddo

3 format(a,4i4,9(1x,g12.5))
 end subroutine energy_strain
!==================================================
 subroutine fix_gauge
 use atoms_force_constants, only : natom_prim_cell
 use constants, only : r15
 use born, only : dyn_na, dyn_naq0,born_flag,bref
 use ios, only : ulog,write_out
 use params, only: verbose
 implicit none
 real(r15) q0(3),b2(3,3)
 complex(r15) dn0(3,3),ddn0(3,3,3)

   q0=(/1.0e-4_r15,0.0_r15,0.0_r15/) 
   call non_anal(q0,1,1,dn0,ddn0,1,0)  ! Parlinski
   bref=real(dn0)
   call write_out(ulog,"DNA 1 ",bref)
   call non_anal(q0,1,1,dn0,ddn0,born_flag,0)
   call write_out(ulog,"DNA 2 ",real(dn0))
   bref=bref-real(dn0)
   call write_out(ulog,"difference 1-2 Bref ",bref)

   call non_anal(q0,1,2,dn0,ddn0,born_flag,0)
   call write_out(ulog,"DNA 2 OD ",real(dn0))

   q0=(/3.0e-5_r15,7.0e-5_r15,5.0e-5_r15/) 
   call non_anal(q0,1,1,dn0,ddn0,1,0)
   b2=real(dn0)
   call non_anal(q0,1,1,dn0,ddn0,born_flag,0)
   b2=b2-real(dn0)
   call write_out(ulog,"B2   ",b2)

3 format(a,3i4,9(1x,f9.4))

 end subroutine fix_gauge
!==================================
 subroutine smooth2step(q,ndn,dyn,ddyn)
!! changes the phase of the dynamical matrix from smooth to step 
 use atoms_force_constants
 use constants, only : r15,ci
 use geometry
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: q(3)
 complex(r15), intent(inout) :: dyn(ndn,ndn), ddyn(ndn,ndn,3)
 integer tau,taup,al,be,nat,ia,ib
 real(r15) dtau(3)

 nat = ndn/3

 do tau=1,nat
 do taup=1,nat
 do al=1,3
 do be=1,3
    ia=al+3*(tau-1)
    ib=be+3*(taup-1)
    dtau = v2a(atom0(taup)%equilibrium_pos-atom0(tau)%equilibrium_pos)
    dyn(ia,ib)=dyn(ia,ib) * exp(-ci*dot_product(q,dtau))
    ddyn(ia,ib,:)=ddyn(ia,ib,:) - ci*dtau(:) * dyn(ia,ib)
 enddo
 enddo
 enddo
 enddo

 end subroutine smooth2step
!==================================
 subroutine step2smooth(q,ndn,dyn,ddyn)
!! changes the phase of the dynamical matrix from step to smooth
 use atoms_force_constants
 use constants, only : r15,ci
 use geometry
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: q(3)
 complex(r15), intent(inout) :: dyn(ndn,ndn), ddyn(ndn,ndn,3)
 integer tau,taup,al,be,nat,ia,ib
 real(r15) dtau(3)

 nat = ndn/3

 do tau=1,nat
 do taup=1,nat
 do al=1,3
 do be=1,3
    ia=al+3*(tau-1)
    ib=be+3*(taup-1)
    dtau = v2a(atom0(taup)%equilibrium_pos-atom0(tau)%equilibrium_pos)
    dyn(ia,ib)=dyn(ia,ib) * exp(ci*dot_product(q,dtau))
    ddyn(ia,ib,:)=ddyn(ia,ib,:) + ci*dtau(:) * dyn(ia,ib)
 enddo
 enddo
 enddo
 enddo

 end subroutine step2smooth
!===========================================================
 subroutine crystal_frame(rot,csys)
!! rot (rows: the new x, y, z in the present Cartesian frame) puts the crystal's symmetry axes
!! along x, y, z -- the conventional crystallographic setting in which tensors (the elastic
!! constants) have the fewest independent components:
!!   cubic: 4-fold (or, for T/Th, 2-fold) axes on x,y,z; hexagonal, trigonal, tetragonal: the
!!   6-, 3- or 4-fold axis on z, a 2-fold axis (else the shortest lattice vector) on x;
!!   orthorhombic: the three 2-fold axes; monoclinic: the 2-fold axis on y; triclinic: unchanged.
!! Improper operations count with their proper part -S (a mirror gives a 2-fold axis).
 use atoms_force_constants, only : op_matrix,isgop,isgopcount
 use lattice, only : r01,r02,r03
 use geometry
 use constants, only : r15,pi
 implicit none
 real(r15), intent(out) :: rot(3,3)
 character(len=*), intent(out) :: csys
 real(r15) s(3,3),ax(3,48),a(3),c1(3),c2(3),c3(3),x(3),y(3),z(3),tr,cs,dt,best
 integer ordr(48),na,i,j,n,n2,n3,n4,n6,i1,i2,i3,iz,ix
 logical found
 na=0
 do i=1,isgopcount
    s=op_matrix(:,:,isgop(i))
    dt=s(1,1)*(s(2,2)*s(3,3)-s(2,3)*s(3,2))-s(1,2)*(s(2,1)*s(3,3)-s(2,3)*s(3,1)) &
 &    +s(1,3)*(s(2,1)*s(3,2)-s(2,2)*s(3,1))
    if (dt.lt.0) s=-s
    tr=s(1,1)+s(2,2)+s(3,3)
    cs=max(-1d0,min(1d0,(tr-1)/2))
    if (cs.gt.0.999d0) cycle                            ! identity
    n=nint(2*pi/acos(cs))
    do j=1,3 ; s(j,j)=s(j,j)-1 ; enddo                ! axis: null vector of S-1
    c1=cross3(s(1,:),s(2,:)) ; c2=cross3(s(2,:),s(3,:)) ; c3=cross3(s(3,:),s(1,:))
    a=c1
    if (sum(c2**2).gt.sum(a**2)) a=c2
    if (sum(c3**2).gt.sum(a**2)) a=c3
    a=a/sqrt(sum(a**2))
    found=.false.
    do j=1,na                                         ! same axis already found
       if (abs(dot_product(a,ax(:,j))).gt.0.999d0) then
          ordr(j)=max(ordr(j),n) ; found=.true. ; exit
       endif
    enddo
    if (.not.found) then
       na=na+1 ; ax(:,na)=a ; ordr(na)=n
    endif
 enddo
 n2=count(ordr(1:na).eq.2) ; n3=count(ordr(1:na).eq.3)
 n4=count(ordr(1:na).eq.4) ; n6=count(ordr(1:na).eq.6)
 iz=0 ; ix=0
 if (n3.ge.4) then
    csys='cubic'
    n=4 ; if (n4.lt.3) n=2
    do i=1,na
       if (ordr(i).ne.n) cycle
       if (iz.eq.0) then
          iz=i
       elseif (ix.eq.0 .and. abs(dot_product(ax(:,i),ax(:,iz))).lt.1d-3) then
          ix=i
       endif
    enddo
 elseif (n6.ge.1 .or. n3.eq.1 .or. n4.eq.1) then
    csys='hexagonal'
    if (n6.eq.0 .and. n3.eq.1) csys='trigonal'
    if (n4.eq.1) csys='tetragonal'
    do i=1,na
       if (ordr(i).ge.3) iz=i
    enddo
    do i=1,na
       if (ordr(i).eq.2 .and. abs(dot_product(ax(:,i),ax(:,iz))).lt.1d-3) then
          ix=i ; exit
       endif
    enddo
 elseif (n2.ge.3) then
    csys='orthorhombic'
    iz=1
    do i=2,na
       if (abs(dot_product(ax(:,i),ax(:,iz))).lt.1d-3) then
          ix=i ; exit
       endif
    enddo
 elseif (n2.eq.1) then
    csys='monoclinic'
 else
    csys='triclinic'
 endif
 rot=0 ; rot(1,1)=1 ; rot(2,2)=1 ; rot(3,3)=1
 if (csys.eq.'triclinic') return
 if (csys.eq.'monoclinic') then                     ! unique 2-fold axis on y
    y=ax(:,1)
    call shortest_perp(y,z)
    x=cross3(y,z)
 else
    z=ax(:,iz)
    if (ix.gt.0) then
       x=ax(:,ix)
    else
       call shortest_perp(z,x)
    endif
    x=x-dot_product(x,z)*z ; x=x/sqrt(sum(x**2))
    y=cross3(z,x)
 endif
 rot(1,:)=x ; rot(2,:)=y ; rot(3,:)=z
 contains
 function cross3(u,v) result(w)
 real(r15), intent(in) :: u(3),v(3)
 real(r15) w(3)
 w=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)]
 end function cross3
 subroutine shortest_perp(u,w)
!! w: the shortest lattice vector (combinations -2..2 of r0i) perpendicular to u, normalized
 real(r15), intent(in) :: u(3)
 real(r15), intent(out) :: w(3)
 real(r15) v(3)
 best=huge(1d0) ; w=0
 do i1=-2,2 ; do i2=-2,2 ; do i3=-2,2
    v=i1*v2a(r01)+i2*v2a(r02)+i3*v2a(r03)
    if (sum(v**2).lt.1d-8) cycle
    if (abs(dot_product(v,u)).gt.1d-4*sqrt(sum(v**2))) cycle
    if (sum(v**2).lt.best-1d-8) then
       best=sum(v**2) ; w=v
    endif
 enddo ; enddo ; enddo
 if (best.lt.huge(1d0)) then
    w=w/sqrt(sum(w**2))
 else                                              ! no lattice vector found: any perpendicular
    w=cross3(u,[1d0,0d0,0d0])
    if (sum(w**2).lt.1d-6) w=cross3(u,[0d0,1d0,0d0])
    w=w/sqrt(sum(w**2))
 endif
 end subroutine shortest_perp
 end subroutine crystal_frame
!===========================================================
 subroutine rotate_tensor4(r,a,b)
!! b_ijkl = r_ia r_jb r_kc r_ld a_abcd (r rows: new axes in the old frame)
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: r(3,3),a(3,3,3,3)
 real(r15), intent(out) :: b(3,3,3,3)
 real(r15) t1(3,3,3,3),t2(3,3,3,3),t3(3,3,3,3)
 integer i,j,k,l,m
 t1=0 ; t2=0 ; t3=0 ; b=0
 do i=1,3 ; do j=1,3 ; do k=1,3 ; do l=1,3 ; do m=1,3
    t1(i,j,k,l)=t1(i,j,k,l)+r(i,m)*a(m,j,k,l)
 enddo ; enddo ; enddo ; enddo ; enddo
 do i=1,3 ; do j=1,3 ; do k=1,3 ; do l=1,3 ; do m=1,3
    t2(i,j,k,l)=t2(i,j,k,l)+r(j,m)*t1(i,m,k,l)
 enddo ; enddo ; enddo ; enddo ; enddo
 do i=1,3 ; do j=1,3 ; do k=1,3 ; do l=1,3 ; do m=1,3
    t3(i,j,k,l)=t3(i,j,k,l)+r(k,m)*t2(i,j,m,l)
 enddo ; enddo ; enddo ; enddo ; enddo
 do i=1,3 ; do j=1,3 ; do k=1,3 ; do l=1,3 ; do m=1,3
    b(i,j,k,l)=b(i,j,k,l)+r(l,m)*t3(i,j,k,m)
 enddo ; enddo ; enddo ; enddo ; enddo
 end subroutine rotate_tensor4
!===========================================================
 subroutine report_crystal_frame(uio,rot,csys,cc)
!! writes the rotation into the crystal frame, the rotated cell (reduced coordinates do not
!! change) and the elastic constants there, with the number of distinct nonzero ones
 use ios, only : ulog
 use lattice, only : r01,r02,r03,r1conv,r2conv,r3conv
 use atoms_force_constants, only : natom_prim_cell,atom0
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(in) :: rot(3,3),cc(6,6)
 character(len=*), intent(in) :: csys
 real(r15) cmax,vals(21),c
 integer i,j,k,nv,u,iu
 logical new
 cmax=maxval(abs(cc)) ; nv=0
 do i=1,6 ; do j=i,6
    c=cc(i,j)
    if (abs(c).lt.1d-3*cmax) cycle
    new=.true.
    do k=1,nv
       if (abs(abs(c)-abs(vals(k))).lt.1d-3*cmax) new=.false.
    enddo
    if (new) then
       nv=nv+1 ; vals(nv)=c
    endif
 enddo ; enddo
 do iu=1,2
    u=ulog ; if (iu.eq.2) u=uio
    write(u,*)'CRYSTAL FRAME: ',trim(csys),'; rotation (rows: new x,y,z in the frame of structure.params)'
    write(u,'(3(1x,f12.8))') (rot(i,:),i=1,3)
    write(u,*)'CRYSTAL FRAME: primitive vectors (A) ; conventional vectors (A), both in this frame'
    write(u,'(3(1x,f13.7),4x,3(1x,f13.7))') matmul(rot,v2a(r01)),matmul(rot,v2a(r1conv))
    write(u,'(3(1x,f13.7),4x,3(1x,f13.7))') matmul(rot,v2a(r02)),matmul(rot,v2a(r2conv))
    write(u,'(3(1x,f13.7),4x,3(1x,f13.7))') matmul(rot,v2a(r03)),matmul(rot,v2a(r3conv))
    write(u,*)'CRYSTAL FRAME: Cartesian atom positions in it (reduced coordinates are unchanged)'
    do i=1,natom_prim_cell
       write(u,'(i4,1x,a,3(1x,f13.7))') i,atom0(i)%name,matmul(rot,v2a(atom0(i)%equilibrium_pos))
    enddo
    write(u,*)'CRYSTAL FRAME: elastic constants C_ij (GPa), entries below 1e-3 of the largest set to 0'
    do i=1,6
       write(u,'(6(1x,f10.3))') (merge(cc(i,j),0d0,abs(cc(i,j)).ge.1d-3*cmax),j=1,6)
    enddo
    write(u,*)'CRYSTAL FRAME: number of distinct nonzero elastic constants =',nv
 enddo
 end subroutine report_crystal_frame
!===========================================================
 subroutine symmetrize_tensor4(a,b,dev)
!! b = average of S(a) over the crystal's point-group operations (Cartesian op_matrix of the
!! space-group elements): the part of a allowed by symmetry; dev = max |a - b| removed
 use atoms_force_constants, only : op_matrix,isgop,isgopcount
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: a(3,3,3,3)
 real(r15), intent(out) :: b(3,3,3,3),dev
 real(r15) t(3,3,3,3)
 integer i
 b=0
 do i=1,isgopcount
    call rotate_tensor4(op_matrix(:,:,isgop(i)),a,t)
    b=b+t
 enddo
 b=b/isgopcount
 dev=maxval(abs(a-b))
 end subroutine symmetrize_tensor4
!============================================================
 subroutine gamma_mode_symmetry(ndn)
!! IR and Raman activity of the phonon modes at Gamma from symmetry alone (2026-11-02, user).
!! At q = 0 the modes transform under the space group modulo lattice translations (isomorphic to the
!! point group); operation g takes atom k to img(k,g) (r' = S r + sgfract within tol_sym, as
!! check_born_symmetry) and
!! acts on a displacement pattern by Gam(g)(3(img-1)+a, 3(k-1)+b) = S(a,b). For each set of degenerate
!! modes (eigenvectors E of D at exactly q = 0, where the non-analytic term is direction-averaged and
!! keeps the full symmetry) the characters chi(g) = Tr(E^+ Gam(g) E) give
!!   irreducibility  (1/|G|) sum |chi|^2            (1: one irrep; 2: a pair of complex-conjugate irreps
!!                                                    degenerate by time reversal, or an accidental one;
!!                                                    the acoustic set: sum of squared multiplicities,
!!                                                    e.g. C6v A1+E1 -> 2, C2h Au+2Bu -> 5)
!!   IR active       (1/|G|) sum chi^* Tr S          > 0  (contains the vector representation)
!!   Raman active    (1/|G|) sum chi^* ((Tr S)^2 + Tr S^2)/2 > 0  (symmetric second-rank tensor)
!!   parity          chi(inversion)/dim, when the inversion is an operation (g = +1, u = -1)
!! With Born charges, the IR strength sum over the set of |Z_m|^2, Z_m = sum Z*_tau e_tau(m)/sqrt(m_tau)
!! (as in dielectric). The irreducible representation of each set and the character table (computed by
!! point_group_rules, common/point_group.f90). Written to gamma_modes.dat and the log.
 use ios, only : ulog
 use born, only : born_flag
 use lattice, only : cart_to_prim, prim_to_cart
 use atoms_force_constants, only : atom0, natom_prim_cell, op_matrix, isgop, isgopcount, sgfract
 use geometry, only : v2a
 use constants, only : r15, cnst
 use params, only : tol_sym
 use eigen, only : mysqrt
 use point_group_rules, only : pg_setup, pg_reduce, pg_write_table, pg_nirr, pg_label, pg_partner
 use phonon_rules, only : translation_weight, translation_min
 implicit none
 integer, intent(in) :: ndn
 integer isg, k, kp, n, i, j, a, b, n1, n2, ng, de, ga, ninv, unit
 integer, allocatable :: img(:,:)
 real(r15) s(3,3), d(3), kz(3), vg(3,ndn), eival(ndn), om(ndn), trs, ts2, dt, tol, nirr, nir, nram, &
 &         zz, par, omax
 real(r15), allocatable :: gam(:,:,:)
 complex(r15) eivec(ndn,ndn), chi(48), zm(3)
 character(len=8) :: sir, sram, spar
 character(len=40) :: sirr
 character(len=400) :: gtot
 real(r15), allocatable :: ops(:,:,:), mult(:)
 complex(r15), allocatable :: chitot(:)
 integer ier
 real(r15) mass(natom_prim_cell)
 logical acou

 n = natom_prim_cell
 allocate(img(n,isgopcount), gam(ndn,ndn,isgopcount))
 gam = 0 ; ninv = 0
 do isg = 1, isgopcount
    s = op_matrix(:,:,isgop(isg))
    dt = s(1,1)*(s(2,2)*s(3,3)-s(2,3)*s(3,2))-s(1,2)*(s(2,1)*s(3,3)-s(2,3)*s(3,1)) &
 &     + s(1,3)*(s(2,1)*s(3,2)-s(2,2)*s(3,1))
    if (dt.lt.0 .and. abs(s(1,1)+s(2,2)+s(3,3)+3).lt.1d-6) ninv = isg     ! S = -1
    do k = 1, n
       img(k,isg) = 0
       do kp = 1, n
          if (atom0(kp)%at_type .ne. atom0(k)%at_type) cycle
          d = matmul(cart_to_prim, matmul(s, v2a(atom0(k)%equilibrium_pos)) + sgfract(:,isg) &
 &                                 - v2a(atom0(kp)%equilibrium_pos))
          if (sqrt(sum(matmul(prim_to_cart, d - nint(d))**2)) .le. tol_sym) then
             img(k,isg) = kp ; exit
          endif
       enddo
       if (img(k,isg) .eq. 0) then
          write(ulog,*)' GAMMA_MODE_SYMMETRY: no image of atom ',k,' under op ',isg,'; not done'
          return
       endif
       gam(3*(img(k,isg)-1)+1:3*img(k,isg), 3*(k-1)+1:3*k, isg) = s
    enddo
 enddo

 allocate(ops(3,3,isgopcount), chitot(isgopcount))
 do isg = 1, isgopcount
    ops(:,:,isg) = op_matrix(:,:,isgop(isg))
    chitot(isg) = gam(1,1,isg)*0
    do i = 1, ndn
       chitot(isg) = chitot(isg) + gam(i,i,isg)
    enddo
 enddo
 call pg_setup(isgopcount, ops, ulog, ier)
 if (ier .ne. 0) write(ulog,*)' GAMMA_MODE_SYMMETRY: point group setup failed (ier=',ier,'); no irreps'
 allocate(mult(max(1,pg_nirr)))

 kz = 0
 call get_freq(kz, ndn, vg, eival, eivec)        ! exactly q = 0
 do k = 1, n
    mass(k) = atom0(k)%mass
 enddo
 om = mysqrt(eival)*cnst                          ! 1/cm (negative = unstable)
 omax = maxval(abs(om))
 tol = max(1d-3, 1d-6*omax)                       ! degeneracy tolerance (1/cm)

 open(newunit=unit, file='gamma_modes.dat')
 write(unit,'(a)')'# phonon modes at Gamma (q = 0; polar modes: non-analytic term direction-averaged, between TO and LO)'
 write(unit,'(a,i3,a)')'# symmetry from the ',isgopcount,' space-group operations modulo lattice translations'
 if (ier .eq. 0) then
    write(unit,'(a)')'# character table (computed; labels: A/B/E/T, g/u, prime, serial number, not Mulliken numbering)'
    call pg_write_table(unit, '#  ')
    call pg_write_table(ulog, ' GAMMA_MODE_SYMMETRY: ')
    call pg_reduce(chitot, mult)
    gtot = ''
    do j = 1, pg_nirr
       if (nint(mult(j)) .gt. 0) then
          write(sirr,'(i0,a)') nint(mult(j)), trim(pg_label(j))
          if (nint(mult(j)) .eq. 1) sirr = pg_label(j)
          if (len_trim(gtot) .gt. 0) gtot = trim(gtot)//' + '
          gtot = trim(gtot)//trim(sirr)
       endif
    enddo
    write(unit,'(a)')'# all modes (incl. acoustic): Gamma = '//trim(gtot)
    write(ulog,*)' GAMMA_MODE_SYMMETRY: Gamma = '//trim(gtot)
 endif
 write(unit,'(a)')'# modes   freq(1/cm)  dim  irred   IR     Raman  parity  IR strength sum|Z_m|^2 '// &
 &                 '(e^2/amu; with Born charges)  irrep'
 write(ulog,*)' GAMMA_MODE_SYMMETRY: modes, freq(1/cm), dim, irreducibility (1 = one irrep), IR, Raman, parity, IR strength'
 n1 = 1
 do while (n1 .le. ndn)
    n2 = n1
    do while (n2 .lt. ndn)
       if (abs(om(n2+1) - om(n1)) .gt. tol) exit
       n2 = n2 + 1
    enddo
    ng = n2 - n1 + 1
    nirr = 0 ; nir = 0 ; nram = 0
    do isg = 1, isgopcount
       chi(isg) = 0
       do i = n1, n2
          chi(isg) = chi(isg) + dot_product(eivec(:,i), matmul(gam(:,:,isg), eivec(:,i)))
       enddo
       s = op_matrix(:,:,isgop(isg))
       trs = s(1,1)+s(2,2)+s(3,3)
       ts2 = 0
       do a = 1, 3
       do b = 1, 3
          ts2 = ts2 + s(a,b)*s(b,a)
       enddo
       enddo
       nirr = nirr + abs(chi(isg))**2
       nir  = nir  + real(conjg(chi(isg))*trs)
       nram = nram + real(conjg(chi(isg))*(trs**2 + ts2)/2)
    enddo
    nirr = nirr/isgopcount ; nir = nir/isgopcount ; nram = nram/isgopcount
    sir = 'no' ; sram = 'no' ; spar = '-'
    if (nir .gt. 0.5d0)  sir = 'yes'
    if (nram .gt. 0.5d0) sram = 'yes'
    acou = .true.             ! rigid translations: mass-weighted translation weight > translation_min
    do i = n1, n2             ! (as freq_vel; their omega at q = 0 need not be 0, e.g. with Born charges)
       if (translation_weight(n, mass, eivec(:,i)) .le. translation_min) acou = .false.
    enddo
    if (acou) then
       sir = 'acoustic' ; sram = '-'
    endif
    if (ninv .gt. 0) then
       par = real(chi(ninv))/ng
       if (par .gt. 0.5d0) spar = 'g'
       if (par .lt. -0.5d0) spar = 'u'
    endif
    zz = 0
    if (born_flag .ne. 0) then
       do i = n1, n2
          zm = 0
          do k = 1, n
          do de = 1, 3
          do ga = 1, 3
             zm(de) = zm(de) + atom0(k)%charge(de,ga)*eivec(ga+3*(k-1),i)/sqrt(atom0(k)%mass)
          enddo
          enddo
          enddo
          zz = zz + real(dot_product(zm,zm))
       enddo
    endif
    sirr = ''
    if (ier .eq. 0) then
       call pg_reduce(chi(1:isgopcount), mult)
       do j = 1, pg_nirr
          if (nint(mult(j)) .gt. 0) then
             if (len_trim(sirr) .gt. 0) sirr = trim(sirr)//'+'
             if (nint(mult(j)) .gt. 1) write(sirr,'(a,i0)') trim(sirr), nint(mult(j))
             sirr = trim(sirr)//trim(pg_label(j))
          endif
       enddo
    endif
    write(unit,'(i3,a,i3,f12.3,i5,f7.3,2x,a8,1x,a8,1x,a6,es14.5,2x,a)') n1,'-',n2, sum(om(n1:n2))/ng, ng, nirr, &
 &        sir, sram, spar, zz, trim(sirr)
    write(ulog,'(a,i3,a,i3,f12.3,i5,f7.3,2x,a8,1x,a8,1x,a6,es14.5,2x,a)')' GAMMA_MODE_SYMMETRY:',n1,'-',n2, &
 &        sum(om(n1:n2))/ng, ng, nirr, sir, sram, spar, zz, trim(sirr)
    if (abs(nirr - nint(nirr)).gt.1d-3 .or. abs(nir - nint(nir)).gt.1d-3 .or. abs(nram - nint(nram)).gt.1d-3) &
 &     write(ulog,*)' GAMMA_MODE_SYMMETRY: WARNING: non-integer multiplicities (degenerate set split?)', &
 &     nirr, nir, nram
    n1 = n2 + 1
 enddo
 close(unit)
 deallocate(img, gam, ops, chitot, mult)
 end subroutine gamma_mode_symmetry
