!----------------------------------------------------
      module kp
      use Structure_info
      implicit none
      integer nkpoint
      real(8),DIMENSION(d) :: kp_gamma
      type(vector),allocatable :: kvector(:) !each element of this array itself is a d-dimensional vector
      !real(8),allocatable :: kvector(:)
      !real(8) , allocatable :: wk(:)

      contains

        subroutine allocatek(n)
            implicit none
            integer :: i
            integer,intent(in) :: n

            IF(ALLOCATED(kvector)) DEALLOCATE(kvector)

            allocate ( kvector(n))!,wk(n) )
            DO i=1,n
                kvector(i)%component=0
                !wk(i)=0
            END DO
            return
        end subroutine allocatek
!-------------------------------------------------------------------------------------------
        subroutine generate_irk(n)
        ! generate the irreducible n kpoints in (0,pi/R]
            implicit none
            integer :: i
            integer, intent(in) :: n
            real(8) :: R,interval
            R=SQRT(trans_vec(:,1).dot.trans_vec(:,1)) !WARNING: since the R here is updated, thus k is updated using the new R
            interval=pi/n/R                           !in the case where k is dot product with R, the R should also use updated R
            do i=1,n
                kvector(i)%component(1)=pi/R-(i-1)*interval
            end do
        end subroutine
!-------------------------------------------------------------------------------------------
        subroutine fill_k(n)
            implicit none
            integer :: i,j
            integer,intent(in) :: n
            do i=2,int((n+1)/2)
                kvector(n+2-i)%component=-kvector(i)%component
            end do

        end subroutine fill_k

        subroutine shift_k(n)
            implicit none
            integer,intent(in)::n
            integer::i
            real(8)::shift
            shift=(trans_vec(:,1).dot.trans_vec(:,1))/(n+1)
            do i=1,n
                kvector(i)%component(1)=kvector(i)%component(1)-shift
            end do
        end subroutine shift_k

        subroutine eliminate_gamma(n)
            implicit none
            integer,intent(in) :: n
            integer :: i,j
            real(8) :: shift(d),temp(d)
            shift=kvector(2)%component
            kvector(1)%component=kvector(int(n/2)+2)%component-shift
        !###########################re-order k###########################
            temp=kvector(1)%component;j=0
            do i=1,n-1
                if(kvector(i+1)%component(1).gt.0) then
                    kvector(i)=kvector(i+1)
                    j=i
                end if
            end do
            kvector(j+1)%component=temp
        end subroutine eliminate_gamma
!-------------------------------------------------------------------------------------
      end module kp
	  !==========================

