 program read_force_position_data
  implicit none
   integer utraj,upos,maxl,i,t,j,k,natom,nf,ne
    real(8) toten
    real(8), allocatable :: pos(:,:),forc(:,:)
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
 allocate(pos(3,natom),forc(3,natom))

! ionic step t = the t-th TOTAL-FORCE block and the t-th TOTEN; each is found by its own search from
! the top of the file (find_nth rewinds), so their order in the OUTCAR does not matter
  nf=count_word(utraj,'TOTAL-FORCE')
  ne=count_word(utraj,'free  energy   TOTEN  =')
  write(*,*)' TOTAL-FORCE blocks, TOTEN lines =',nf,ne
  if (nf.ne.ne) write(*,*)' WARNING: their numbers differ; the first min of the two steps are written'
  do t=1,min(nf,ne)
     call find_nth(utraj,'TOTAL-FORCE',t,line,k)
     read(utraj,*)
     do i=1,natom
        read(utraj,*) pos(1:3,i),forc(1:3,i)
     enddo
     call find_nth(utraj,'free  energy   TOTEN  =',t,line,k)
     read(line(k+23:98),*)toten
     write(*,*)' step,TOTEN=',t,toten
     write(upos,*) ' # POSITION      TOTAL FORCE '
     write(upos,*) t,toten
     do i=1,natom
        write(upos,6) pos(1:3,i),forc(1:3,i)
     enddo
  enddo
  write(*,*)' steps written to pos-forc.dat: ',min(nf,ne)
6 format(2(2x,3(2x,g16.10)))
  close(utraj)
  close(upos)

 contains
!------------------------------------------------------------
 integer function count_word(u,word)
!! number of lines of the file (from its top) that contain word
 integer, intent(in) :: u
 character(*), intent(in) :: word
 character line*99
 logical found
 integer k
 rewind(u)
 count_word=0
 do
    read(u,'(a)',end=1)line
    call findword(word,line,found,k)
    if (found) count_word=count_word+1
 enddo
1 rewind(u)
 end function count_word
!------------------------------------------------------------
 subroutine find_nth(u,word,n,line,k)
!! rewinds and reads to the n-th line containing word; the file is left right after it, word at line(k:)
 integer, intent(in) :: u,n
 character(*), intent(in) :: word
 character(*), intent(out) :: line
 integer, intent(out) :: k
 logical found
 integer m
 rewind(u)
 m=0
 do
    read(u,'(a)',end=1)line
    call findword(word,line,found,k)
    if (found) then
       m=m+1
       if (m.eq.n) return
    endif
 enddo
1 write(*,*)' FIND_NTH: occurrence',n,' of "',word,'" not found'
 stop
 end subroutine find_nth

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

