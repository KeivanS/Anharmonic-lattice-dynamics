module mpi_params
    implicit none

    integer :: nprocs, mpi_rank, mpi_err
    integer :: kpart, kremain, mynumk, mystartk, myendk, myeivals, myeivecs

    integer, allocatable, dimension(:) :: offsets, vals_displs, vecs_displs
    integer, allocatable, dimension(:) :: vals_count, vecs_count
end module mpi_params
