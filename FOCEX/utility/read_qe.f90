! from the output file of PWSCF
 program read_force_position_data

  implicit none
  integer utraj,upos,maxl,i,t,j,k,natom
  real(8),allocatable :: pos(:,:)
  real(8) forc(3),toten,a0,ab,ryd
  character line*99,fin*66
  logical found,convert

 ab=0.5291772108
 ryd=13.6056923
 convert=.FALSE.
!     print*,' input the name of the PWSCF output file '
      read(*,*) fin
      fin=fin(1:len_trim(fin))
 !    print*, 'reading as input filename ',fin,' of length=',len_trim(fin)
      utraj = 10; upos = 20
      maxl = 999999999
  open(utraj,file=fin,status='old')
  open(upos,file='pos-forc.dat',status='unknown')

! first get the number of atoms
   fstloop: do j=1,maxl
       read(utraj,'(a)')line
       call findword('number of atoms/cell      =',line,found,k)
       if (found) then
          read(line(33:98),*)natom
          exit fstloop
       endif
   enddo fstloop
   write(*,*)' read the number of atoms, now rewinding '
   rewind(utraj)
   print*,' Detected ',natom,' atoms'
   allocate ( pos(natom,3))

! then find the position-force data
   t=0
   do j=1,maxl

! find and read positions
       read(utraj,'(a)',end=88)line
       call findword('lattice parameter',line,found,k)
       if (found) then
           read(line(33:98),*)a0
           write(*,*)' lattice parameter is=',a0,' a.u.'
           a0 = a0*ab   ! convert to angstroem
           write(*,*)' lattice parameter is=',a0,' Ang'
           cycle
       endif
       call findword(' site n.     atom ',line,found,k)
       if (found) then
           do i=1,natom
              read(utraj,'(a)')line
              read(line(39:98),*) (pos(i,k),k=1,3)
           enddo
           pos = pos*a0   ! converts them from reduced to cartesian
       endif
       call findword('!    total energy',line,found,k)
       if (found) then
           read(line(33:98),*)toten
           if (line(51:52).eq.'Ry') toten = toten*ryd  ! convert to eV if necessary
           write(*,*)' j,TOTEN(eV)=',j,toten
           cycle
       endif
       call findword('Forces acting on atoms',line,found,k)
       if (found) then
        call findword('Ry/au',line,found,k)
!       if (line(30:34) .eq. 'Ry/au') convert=.TRUE.  ! convert to eV/Ang if necessary
        if (found) convert=.TRUE.  ! convert to eV/Ang if necessary
        t = t+1
        write(*,*) ' # POSITION (Ang)     TOTAL FORCE (eV/Ang) '
        write(*,*) t,toten,'= t,Etot(eV)'
        write(upos,*) ' # POSITION (Ang)     TOTAL FORCE (eV/Ang) '
        write(upos,*) t,toten,'= t,Etot(eV)'
        kloop: do k=1,5000
          read(utraj,'(a)',end=88)line
          if(line(6:19) .eq. 'atom    1 type') then
             i=1
             read(line(33:98),*) forc(1:3) 
             if ( convert ) forc=forc*ryd/ab  ! convert to eV/Ang if necessary
             write(*,6) pos(i,1:3),forc(1:3) 
             write(upos,6) pos(i,1:3),forc(1:3) 
             exit kloop
          endif
        enddo kloop
        do i=2,natom
           read(utraj,'(a)',end=88)line
           read(line(33:98),*) forc(1:3)
           if ( convert ) forc=forc*ryd/ab  ! convert to eV/Ang if necessary
           write(*,6) pos(i,1:3),forc(1:3) 
           write(upos,6) pos(i,1:3),forc(1:3) 
        enddo
       endif
   enddo
88 write(*,*)' reached the end of OUTCAR file after steps= ',t
 6 format(3(1x,g12.6),3x,3(1x,g14.8))
   write(*,*)' last line read was '
   write(*,'(a)') line
   close(utraj)
   close(upos)

  end program read_force_position_data
!============================================================
 subroutine findword(word,line,found,i)
 implicit none
 logical found
 character(*) line, word
 integer i,l,k
 l = len_trim(line); k=len_trim(word)
 found = .False.
 do i=1,l-k+1
    if(line(i:i+k-1) .eq. word(1:k)) then
       found = .True.
       exit
    endif
 enddo
 end subroutine findword
!============================================================

