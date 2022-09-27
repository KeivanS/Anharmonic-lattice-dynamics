program offdisplacement
integer(8) invnumber,i,j,k,p,fcid,total_lines,reason,lines,header,numFC,indexp, &
lend,fid,numxyz,N_snapshots,size_of_supercell,fdconsval,startfd,fidfc,amat_col,nsnapshot,fider
integer(8) allsnap,newsnaper
character :: fname*20, each_line*199, fdconstraint*199, word*4
real(8), allocatable :: amatrix(:,:), forceconstant(:), bmat(:), axm(:), & 
prodmbmat(:),snap_err(:),fctemp(:),force_temp(:),snaper(:)
logical examat,exfcval,found, exfc
real(8) aval, bval, cval, nsnapval, snaperror, avgerr,avg_erroir,var,sd,sdval
inquire(file='amatrx.dat',exist=examat)
inquire(file='force-value.dat',exist=exfcval)
inquire(file='fc-temp-all.dat',exist=exfc)
if(examat .and. exfcval .and. exfc) then
write(*,*) "amatrx.dat, force-value.dat and fc-temp-all.dat present in the current working directory"
else
write(*,*) "amatrx.dat or force-value.dat or fc-temp-all.dat not present in the current working directory, rerun fc234-13 again"
stop
endif
fid=123
fidfc=456
fider=789
fcid=369
allsnap=8639
newsnaper=333
total_lines=1234567899
numxyz=3
header=7
lines=0
!---------------------------------------------
open(fidfc,file='force-value.dat')
read(fidfc,*,iostat=reason) size_of_supercell,numFC
allocate(forceconstant(numFC),fctemp(numFC))
do i=1, numFC
read(fidfc,*,iostat=reason) forceconstant(i)
enddo
amat_col=numFC+1
!---------------------------------------------
! This is to read number of force displacement constraint in the Ax=b system
open(fid,file='amatrx.dat')
do i=1,2
read(fid,'(A)') fdconstraint
enddo
close(fid)
! ----------------------------------------------------------------------------------
open(fid, file='amatrx.dat')
do i=1,total_lines
read(fid,'(A)',iostat=reason) each_line
if (reason>0) then
write(*,*) "Something is wrong with amatrx.dat file, please check this file"
exit
elseif (reason<0) then
write(*,*) "...END OF FILE..."
exit
else
lines=lines+1
endif
enddo
write(*,*) lines
lend=len(fdconstraint)-len(trim(fdconstraint))
word='size'
call findword(word,fdconstraint,found,indexp)
if (found) then
read(fdconstraint(indexp+15:len(trim(fdconstraint))),*) fdconsval
endif
write(*,*) fdconsval
startfd=lines-fdconsval
write(*,*) startfd
close(fid)
invnumber=startfd-invnumber
! ------------------------------------------------------------------------------------------------------------
! READ TEMPERATURE DEPENDENT FORCE CONSTANT FROM THE fc-temp-all.dat FILE
open(fcid,file='fc-temp-all.dat')
do i=1, numFC
read(fcid,*) fctemp(i)
enddo
close(fcid)
! ------------------------------------------------------------------------------------------------------------
! Below this line is for separating each of the force constant
allocate(amatrix(fdconsval,amat_col))
j=1
open(fid,file='amatrx.dat')
do i=1,startfd
read(fid,'(A)',iostat=reason) each_line
enddo
do i=1, fdconsval
read(fid,*,iostat=reason) amatrix(i,:)
if (reason > 0) then
write(*,*) "Something is wrong with the file"
exit
elseif (reason < 0) then
write(*,*) "...END OF FILE..."
exit
else
j=j+1
endif
enddo
close(fid)
allocate(force_temp(fdconsval))
! -------------------------------------------------------------------------------------------------------------------------------
! EVALUATE THE A*x A--> displacement matrix and x--> force constant AT THAT TEMPERATURE TO GET FORCE VALUE AT THAT TEMPERATURE
do i=1,fdconsval
force_temp(i)=dot_product(amatrix(i,1:numFC-1),fctemp(1:numFC))
enddo
! --------------------------------------------------------------------------------------------------------------------------------
allocate(bmat(fdconsval))
allocate(prodmbmat(fdconsval))
allocate(axm(fdconsval))
!------------------------------------------------------------------------------------------------------------------------------------------------------------
! Get the value of b(DFT force), prod-b (Ax-b) [displacement matrix]*FCs-b
! since there are 3*num_snapshots*super_cell_size+invariance_relation_line
do i=1,invnumber
read(fidfc,*) aval, bval, cval
enddo 
do i=1,fdconsval
read(fidfc,*,iostat=reason) bmat(i), axm(i), prodmbmat(i)
if (reason < 0) then
write(*,*) "...END OF FILE... for force-value.dat"
endif
enddo
close(fidfc)
!------------------------------------------------------------------------------------------------------------------------------------------------------------
nsnapval=fdconsval/(3*size_of_supercell)
nsnapshot=int(nsnapval,8)
allocate(snap_err(nsnapshot))
allocate(snaper(nsnapshot))
p=0
avgerr=0.0d0
avg_error=0.0d0
do i=1,nsnapshot
snaperror=0.0d0
do j=1,3*size_of_supercell
p=p+1
!snaperror=snaperror+prodmbmat(p)*prodmbmat(p) ! if in case for 0 K use this one 
snaperror=snaperror+(bmat(i)-force_temp(i))*(bmat(i)-force_temp(i)) ! this is for temperature dependent force constant case
enddo
snaper(i)=sqrt(snaperror/(3*size_of_supercell))
snap_err(i)=sqrt(snaperror)
avgerr=avgerr+snap_err(i)
avg_error=avg_error+snaper(i)
enddo
avg_error=avg_error/nsnapshot
!----------Calculate the variance in error from all the snapshot--------------------------------
var=0.0d0
do i=1,nsnapshot
var=var+(snaper(i)-avg_error)*(snaper(i)-avg_error)
enddo
sd=sqrt(var/(nsnapshot-1))
!------------------------------------------------------------------------------------------------------------------------------------------------------------
open(allsnap,file='allsnapshots.txt')
do i=1,nsnapshot
write(allsnap,*) i,",",snap_err(i)
enddo
close(allsnap)
!------------------------------------------------------------------------------------------------------------------------------------------------------------
avgerr=avgerr/nsnapshot
j=0
write(*,*) "*********************SNAPSHOTS GREATER THAN AVERAGE ERROR**********************************"
open(fider,file='snap.txt')
write(fider,*) "Snapshot, error"
do i=1,nsnapshot
if (snap_err(i) .lt. avgerr) then
write(fider,*) i, snap_err(i)
j=j+1
endif
if (snap_err(i) .gt. sd) then
write(*,*) i
endif
enddo
write(*,*) "***************************************STATS*********************************************"
write(*,*) "The average error is: ",avgerr
write(*,*) "The number of snapshots below average error is: ",j
close(fider)
! ---- WRITE IN ANOTHER FILE ABOUT EACH SNAPSHOT ERROR -----
open(newsnaper,file='newsnap.txt')
do i=1, nsnapshot
write(newsnaper,*) i,",",snaper(i)
enddo
close(newsnaper)
!------------------------------------------------------------------------------------------------------------------------------------------------------------
deallocate(amatrix,snaper,forceconstant,bmat,prodmbmat,axm,snap_err,fctemp,force_temp)
end program offdisplacement

subroutine findword(word,line,found,indexp)
! says if the string "word" exists in string "line"
 implicit none
 logical found
 character(*) line, word
 integer i,l,k
 integer(8) indexp
 l = len_trim(line); k=len_trim(word)
 found = .False.
 do i=1,l-k+1
    if(line(i:i+k-1) .eq. word(1:k)) then
       found = .True.
       indexp=i
       exit
    endif
 enddo
end subroutine findword
