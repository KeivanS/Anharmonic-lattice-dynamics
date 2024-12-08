 program read_force_position_data
! this program reads the OUTCAR file (vasp output) and writes 
! the position-force data in a file called pos-forc.dat
  implicit none
   integer utraj,upos,maxl,i,t,j,k,natom
    real(8) pos(3),forc(3),toten
     character line*99
      logical found

       utraj = 10; upos = 20
        maxl = 999999999
 open(utraj,file='OUTCAR',status='old')
  open(upos,file='pos-forc.dat',status='unknown')

   do j=1,maxl
       read(utraj,'(a)',end=98)line
           call findword('NIONS =',line,found,k)
       if (found) then
              write(*,*)' found is true ; at line=',j
             read(line(k+8:98),*)natom
            write(*,*)' NATOMS=',natom
           exit
       endif
        enddo
98 write(*,*)' read the number of atoms, now rewinding '
 rewind(utraj)

 ! now get the FORCES from OUTCAR file
  t=0
   do j=1,maxl
       read(utraj,'(a)',end=88)line
           call findword('free  energy   TOTEN  =',line,found,k)
       if (found) then
              read(line(k+25:98),*)toten
             write(*,*)' j,TOTEN=',j,toten
            cycle
        endif
    call findword('POSITION',line,found,i)
        if (found) then
       read(utraj,*)
              t = t+1
             write(upos,*) ' # POSITION      TOTAL FORCE '
            write(upos,*) t,toten
           do i=1,natom
              read(utraj,*) pos(1:3),forc(1:3)
                 write(upos,6) pos(1:3),forc(1:3)
        enddo
    endif
     enddo
     88 write(*,*)' reached the end of OUTCAR file after steps= ',t
     6 format(2(2x,3(2x,g16.10)))
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

