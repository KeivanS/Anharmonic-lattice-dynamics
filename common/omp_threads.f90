!======================================================
 subroutine omp_threads_choose(envname,nthr,src)
!! Number of OpenMP threads (FOCEX, THERMACOND; 2026-11-04), chosen here, not from OMP_NUM_THREADS
!! (often set to 1 globally for MPI codes), and set with omp_set_num_threads:
!!   the variable envname (FOCEX_NUM_THREADS, TCOND_NUM_THREADS)   if set: explicit choice
!!   1                       if started by an MPI launcher (OMPI_COMM_WORLD_SIZE, PMI_SIZE or
!!                           MPI_LOCALNRANKS set): the cores belong to the other ranks
!!   SLURM_CPUS_PER_TASK     if set: the cores Slurm gave this task
!!   n - max(1,n/5)          otherwise, n = omp_get_num_procs: about 20 % (at least 1) left for
!!                           other work (user: 8 of the 10 cores of the Mac)
!! Without -fopenmp: 1, nothing set.
!$ use omp_lib
 implicit none
 character(*), intent(in) :: envname
 integer, intent(out) :: nthr
 character(*), intent(out) :: src
 integer ios,l
 character(64) v
 nthr=1 ; src='built without OpenMP'
!$ nthr=omp_get_num_procs() ; nthr=max(1,nthr-max(1,nthr/5))
!$ write(src,'(a,i0,a)')'default: ',omp_get_num_procs(),' cores minus about 20 %'
!$ call get_environment_variable('SLURM_CPUS_PER_TASK',v,l)
!$ if (l.gt.0) then
!$    read(v,*,iostat=ios) nthr ; if (ios.ne.0 .or. nthr.lt.1) nthr=1
!$    src='SLURM_CPUS_PER_TASK'
!$ endif
!$ call get_environment_variable('OMPI_COMM_WORLD_SIZE',v,l)
!$ if (l.eq.0) call get_environment_variable('PMI_SIZE',v,l)
!$ if (l.eq.0) call get_environment_variable('MPI_LOCALNRANKS',v,l)
!$ if (l.gt.0) then
!$    nthr=1 ; src='started by an MPI launcher'
!$ endif
!$ call get_environment_variable(envname,v,l)
!$ if (l.gt.0) then
!$    read(v,*,iostat=ios) nthr ; if (ios.ne.0 .or. nthr.lt.1) nthr=1
!$    src=envname
!$ endif
!$ call omp_set_num_threads(nthr)
 end subroutine omp_threads_choose
