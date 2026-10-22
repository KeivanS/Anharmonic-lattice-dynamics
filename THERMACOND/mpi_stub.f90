!! Single-process stand-in for the MPI library, used by "make serial" so that THERMACOND
!! builds with a plain Fortran compiler.  It provides only what kappa.f90 calls; every
!! collective operation reduces to a copy, since there is one process (rank 0 of 1).
 module mpi
 implicit none
 integer, parameter :: MPI_COMM_WORLD=0, MPI_INTEGER=1, MPI_DOUBLE_PRECISION=2
 integer, parameter :: MPI_STATUS_SIZE=1

 contains

 subroutine MPI_INIT(ierr)
 integer ierr
 ierr=0
 end subroutine MPI_INIT

 subroutine MPI_FINALIZE(ierr)
 integer ierr
 ierr=0
 end subroutine MPI_FINALIZE

 subroutine MPI_COMM_SIZE(comm,n,ierr)
 integer comm,n,ierr
 n=1 ; ierr=0
 end subroutine MPI_COMM_SIZE

 subroutine MPI_COMM_RANK(comm,r,ierr)
 integer comm,r,ierr
 r=0 ; ierr=0
 end subroutine MPI_COMM_RANK

 subroutine MPI_BARRIER(comm,ierr)
 integer comm,ierr
 ierr=0
 end subroutine MPI_BARRIER

 subroutine MPI_ABORT(comm,code,ierr)
 integer comm,code,ierr
 write(*,*)'MPI_ABORT (serial build), code ',code
 error stop 1
 end subroutine MPI_ABORT

 subroutine MPI_ALLGATHER(sendval,ns,stype,recv,nr,rtype,comm,ierr)
 integer sendval,ns,stype,recv(*),nr,rtype,comm,ierr
 recv(1)=sendval ; ierr=0
 end subroutine MPI_ALLGATHER

 subroutine MPI_ALLGATHERV(send,ns,stype,recv,counts,displs,rtype,comm,ierr)
 real(8) send(*),recv(*)
 integer ns,stype,counts(*),displs(*),rtype,comm,ierr
 recv(displs(1)+1:displs(1)+ns)=send(1:ns) ; ierr=0
 end subroutine MPI_ALLGATHERV

 subroutine MPI_GATHERV(send,ns,stype,recv,counts,displs,rtype,root,comm,ierr)
 real(8) send(*),recv(*)
 integer ns,stype,counts(*),displs(*),rtype,root,comm,ierr
 recv(displs(1)+1:displs(1)+ns)=send(1:ns) ; ierr=0
 end subroutine MPI_GATHERV

 end module mpi
