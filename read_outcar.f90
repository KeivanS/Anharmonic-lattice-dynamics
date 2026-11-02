 program read_force_position_data
  implicit none
   integer utraj,upos,maxl,i,t,j,k,natom
    real(8) toten
    real(8), allocatable :: pos(:,:),forc(:,:)
     character line*99
      logical found,pending

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
 allocate(pos(3,natom),forc(3,natom))

  t=0; toten=0.d0; pending=.false.
   do j=1,maxl
       read(utraj,'(a)',end=88)line
           call findword('free  energy   TOTEN  =',line,found,k)
       if (found) then
              read(line(k+23:98),*)toten
             write(*,*)' j,TOTEN=',j,toten
             if (pending) then
                t = t+1
                write(upos,*) ' # POSITION      TOTAL FORCE '
                write(upos,*) t,toten
                do i=1,natom
                   write(upos,6) pos(1:3,i),forc(1:3,i)
                enddo
                pending=.false.
             endif
            cycle
        endif
    call findword('TOTAL-FORCE',line,found,i)
        if (found) then
           read(utraj,*)
           do i=1,natom
              read(utraj,*) pos(1:3,i),forc(1:3,i)
           enddo
           pending=.true.
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

