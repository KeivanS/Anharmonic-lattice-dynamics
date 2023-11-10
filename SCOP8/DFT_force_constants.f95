MODULE DFT_force_constants
    !!This module is for read in 'real' force constants from DFT calculation(the fc#.dat files)
    USE Structure_info !We use 'atom_number, unit_cell_number, d' of this module

    IMPLICIT NONE

    INTEGER :: fc_terms(4),rank,ifc2_terms
    INTEGER :: eff_fc2_terms !number of fc2 that has 1st atom within primitive cell
    INTEGER :: ifc2_atompairs
    REAL(8) :: max_fc2=0d0
    REAL(8) :: max_fc3=0d0
    REAL(8),DIMENSION(6,6) :: elastic,compliance
    !!Below are FORTRAN TYPE definition for force constants
!---------------------------------------------------------------------------------------------------
    TYPE fc2_index  !index info directly reads from file, will be used to map fc2_value
        INTEGER :: junk,group
        INTEGER :: iatom_number,jatom_number,iatom_xyz,jatom_xyz
        REAL(8) :: phi_temp
    END TYPE

    TYPE fc2_value  !with index info, quickly retrieve corresponding fc value
        REAL(8),DIMENSION(d,d) :: phi !if 3d case, then phi(1,1) means phi_xx, et cetera
    END TYPE

    TYPE fc2_group
        INTEGER :: group(2)
        REAL(8) :: mat(2)
    END TYPE
!----------------------------------------------------------------------------------------------------
    TYPE fc3_index
        INTEGER :: junk,group
        INTEGER :: iatom_number,jatom_number,katom_number,iatom_xyz,jatom_xyz,katom_xyz
        REAL(8) :: psi_temp
    END TYPE

    TYPE fc3_value
        REAL(8),DIMENSION(d,d,d) :: psi
    END TYPE
!-----------------------------------------------------------------------------------------------------
    TYPE fc4_index
        INTEGER :: junk,group
        INTEGER :: iatom_number,jatom_number,katom_number,latom_number,iatom_xyz,jatom_xyz,katom_xyz,latom_xyz
        REAL(8) :: chi_temp
    END TYPE

    TYPE fc4_value
        REAL(8),DIMENSION(d,d,d,d) :: chi
    END TYPE
!------------------------------------------------------------------------------------------------------
    !! Below are variable declarations
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: myfc2_index, oldfc2_index
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: myfc2_value,myefc2_value, prev_fc2, inv_phiTau
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: eff_fc2_index !to store all fc2 that atom1 within p cell
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: indiefc2_index !to store all the independent fc2 index
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: old_indiefc2_index
    TYPE(fc2_group),DIMENSION(:,:,:,:),ALLOCATABLE :: myfc2_group !to categorize all fc2 by indie_index
    TYPE(fc2_group),DIMENSION(:,:,:,:),ALLOCATABLE :: old_myfc2_group

    TYPE(fc3_index),DIMENSION(:),ALLOCATABLE :: myfc3_index, oldfc3_index
    TYPE(fc3_value),DIMENSION(:,:,:),ALLOCATABLE :: myfc3_value,myefc3_value, prev_fc3
    INTEGER, DIMENSION(:),ALLOCATABLE :: fc3_raw_idx, fc3_unique_idx
    INTEGER, DIMENSION(:),ALLOCATABLE :: old_fc3_unique_idx !to store all unique atom index that gives non-zero fc3

    TYPE(fc4_index),DIMENSION(:),ALLOCATABLE :: myfc4_index, oldfc4_index
    TYPE(fc4_value),DIMENSION(:,:,:,:),ALLOCATABLE :: myfc4_value,myefc4_value, prev_fc4
    INTEGER, DIMENSION(:),ALLOCATABLE :: fc4_raw_idx, fc4_unique_idx
    INTEGER, DIMENSION(:),ALLOCATABLE :: old_fc4_unique_idx !to store all unique atom index that gives non-zero fc4


CONTAINS
!******************************************************************************************************
    FUNCTION get_letter(direction) RESULT(letter)
        !!Utility function: for convert direction from integer(123) to character(xyz)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: direction
        CHARACTER :: letter
        SELECTCASE(direction)
            CASE(1)
                letter = 'x'
            CASE(2)
                letter = 'y'
            CASE(3)
                letter = 'z'
            CASE DEFAULT
                WRITE(*,*)'ERROR when convert character'
                STOP
        ENDSELECT
    END FUNCTION get_letter
!-------------------------------------------------------------------------------------
    SUBROUTINE unique_sort(arrayIn, arrayOut)
    !!Utility subroutine: unique_sort
!!read in an unsorted array with duplicated elements and return array with only unique values, sorted
        IMPLICIT NONE
        INTEGER :: i=0, min_val, max_val
        INTEGER, INTENT(in) :: arrayIn(:)
        INTEGER, DIMENSION(:), ALLOCATABLE :: unique
        INTEGER, INTENT(out), DIMENSION(:), ALLOCATABLE :: arrayOut

        ALLOCATE(unique(SIZE(arrayIn)))

        min_val = minval(arrayIn) - 1
        max_val = maxval(arrayIn)

        DO WHILE (min_val < max_val)
            i = i + 1
            min_val = minval(arrayIn, mask=arrayIn>min_val)
            unique(i) = min_val
        END DO

        ALLOCATE(arrayOut(i),source=unique(1:i))
        i=0 ! To reset i is necessary, typical Fortran bs again...
        DEALLOCATE(unique)
    END SUBROUTINE unique_sort


    SUBROUTINE unique_sort2(arrayIn, arrayOut)
    !!unique_sort for double value
        IMPLICIT NONE
        INTEGER :: i=0
        REAL(8) :: min_val, max_val
        REAL(8), INTENT(in) :: arrayIn(:)
        REAL(8), DIMENSION(:), ALLOCATABLE :: unique
        REAL(8), INTENT(out), DIMENSION(:), ALLOCATABLE :: arrayOut

        ALLOCATE(unique(SIZE(arrayIn)))

        min_val = minval(arrayIn) - 1d0
        max_val = maxval(arrayIn)

        DO WHILE ((min_val - max_val) .lt. 0)
            i = i + 1
            min_val = minval(arrayIn, mask=((arrayIn - min_val) .gt. 0))
            unique(i) = min_val
        END DO

        ALLOCATE(arrayOut(i),source=unique(1:i))
        i=0 ! To reset i is necessary, typical Fortran bs again...
        DEALLOCATE(unique)
    END SUBROUTINE unique_sort2
!--------------------------------------------------------------------------
    FUNCTION find_loc(array, element) RESULT(idx)
    !!Utility subroutine: find_loc
    !!find the index of element in array
    !!NOTICE: intrinsic function FINDLOC() is not supported until F2008,
    !!This function may not be used in current build
        IMPLICIT NONE
        INTEGER, INTENT(in) :: array(:), element
        INTEGER :: idx, i=0
        DO i=1, SIZE(array)
            IF (array(i)==element) THEN
                idx = i
                EXIT
            END IF
        END DO
        IF(.NOT.ANY(array==element)) THEN
            WRITE(*,*)'Could not find corresponding atom!'
        END IF
    END FUNCTION find_loc


    FUNCTION find_loc2(array, element) RESULT(idx)
    !!find_loc for double value
        IMPLICIT NONE
        REAL(8), INTENT(in) :: array(:), element
        INTEGER :: idx, i=0
        DO i=1, SIZE(array)
            IF (array(i)==element) THEN
                idx = i
                EXIT
            END IF
        END DO
        IF(.NOT.ANY(array==element)) THEN
            WRITE(*,*)'Could not find corresponding atom!'
            idx = 0
        END IF
    END FUNCTION find_loc2
!-----------------------------------------------------------------------------------------------

    FUNCTION include_arrays(child, parent) RESULT(inc)
    !!check if child array is included in the parent array
        IMPLICIT NONE
        INTEGER,DIMENSION(:),INTENT(in) :: child, parent
        LOGICAL :: inc
        INTEGER :: i
        IF(SIZE(parent).lt.SIZE(child)) inc = .false.

        DO i=1,SIZE(child)
            IF(.NOT.ANY(parent==child(i))) THEN
                inc = .false.
                RETURN
            END IF
        END DO
        inc = .true.
    END FUNCTION include_arrays
!-----------------------------------------------------------------------------------------------
    SUBROUTINE get_atoms_fcs
    !!get all the fc-related atoms index given initial nshells(:,:)
        IMPLICIT NONE
        INTEGER :: i,j, idx, rnk
        INTEGER,ALLOCATABLE, DIMENSION(:) :: atoms

        DO rnk = 2, 4
            idx = 0
            DO j = 1, map(rnk)%ngr
                DO i = 1, map(rnk)%nt(j)
                    idx = idx + 1
                END DO
            END DO

            ALLOCATE(atoms(idx))
            idx = 0
            DO j = 1, map(rnk)%ngr
                DO i = 1, map(rnk)%nt(j)
                    idx = idx + 1
                    atoms(idx) = map(rnk)%gr(j)%iat(2,i)
                END DO
            END DO

            SELECTCASE(rnk)
                CASE(2)
                    CALL unique_sort(atoms,atoms_fc2)
                CASE(3)
                    CALL unique_sort(atoms,atoms_fc3)
                CASE(4)
                    CALL unique_sort(atoms,atoms_fc4)
            ENDSELECT

            DEALLOCATE(atoms)
        END DO

    END SUBROUTINE get_atoms_fcs
!-----------------------------------------------------------------------------
    SUBROUTINE get_atoms_otf(rnk,atoms_otf)
    !!get all the fc-related atoms index given nshells(:,:) on the fly, from map(:)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: rnk
        INTEGER :: i,j, idx
        INTEGER,ALLOCATABLE, DIMENSION(:) :: atoms
        INTEGER,INTENT(out),DIMENSION(:),ALLOCATABLE :: atoms_otf

        idx = 0
        DO j = 1, map(rnk)%ngr
            DO i = 1, map(rnk)%nt(j)
                idx = idx + 1
            END DO
        END DO

        ALLOCATE(atoms(idx))
        idx = 0
        DO j = 1, map(rnk)%ngr
            DO i = 1, map(rnk)%nt(j)
                idx = idx + 1
                atoms(idx) = map(rnk)%gr(j)%iat(2,i) !record the second atom idx
            END DO
        END DO

        CALL unique_sort(atoms, atoms_otf)

        DEALLOCATE(atoms)

    END SUBROUTINE get_atoms_otf
!===============================================================================================
    SUBROUTINE read_force_constants
    !!The major subroutine to read FCs from fc#.dat files(rank#=2,3,4)
        IMPLICIT NONE

        INTEGER :: i,j,k,l,ufc1,ufc2,ufc3,ufc4,mx
        INTEGER :: atom1,atom2,atom3,atom4  !for temporary use
        INTEGER :: direction1,direction2,direction3,direction4 !for temporary use
        REAL(8),DIMENSION(:),ALLOCATABLE :: temporary
        REAL(8) :: check
        CHARACTER(LEN=99) line
        INTEGER,ALLOCATABLE,DIMENSION(:) :: atoms, directions

        check = 0d0

         ufc1=11
         ufc2=12
         ufc3=13
         ufc4=14
         CALL read_input_fit !this is original read_input_fit with some modifications
!WRITE(*,*) 'read_input_fit Ends Normally'
         CALL read_structure !it will only use one number: tot_atom_number
!WRITE(*,*) 'read_structure Ends Normally'
        CALL transform_input_structure !newly added 09/15/2023
         CALL make_reciprocal_lattice(r01,r02,r03,g01,g02,g03)
         CALL set_neighbor_list
         ! inputs: atompos,
         ! outpus: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n
         CALL write_neighbors

         WRITE(*,*) 'dimension=',d
         WRITE(*,*) 'total atom number=',tot_atom_number

!------------------------------------------------------------------------------------------------
         rank=2
         OPEN(ufc2,file='fc2.dat',status='old',action='read')
         READ(ufc2,'(a)') line
         mx =10000000 ! nterms(rank) has to be initialized
         DO j=1,mx
            READ(ufc2,*,END=92) temporary
         END DO
92       mx=j-1
         REWIND(ufc2)

         fc_terms(rank)=mx

         ALLOCATE(myfc2_index(mx),myfc2_value(atom_number,tot_atom_number)) !save memory space
         !obviously, SIZE(myfc2_value) < mx
         DO i=1,atom_number
            DO j=1,tot_atom_number
                myfc2_value(i,j)%phi=0
            END DO
         END DO
         !myfc2_value(:,:)%phi(1:,1:)=0
         i=1;j=1

         READ(ufc2,'(a)') line
         WRITE(*,'(a)') line
         WRITE(*,*) '********** FCs for rank=2:' ,mx,'  ************'

         eff_fc2_terms = 0
         DO i=1,mx
            READ(ufc2,*) myfc2_index(i)%junk,myfc2_index(i)%group,myfc2_index(i)%iatom_number,myfc2_index(i)%iatom_xyz,&
            &            myfc2_index(i)%jatom_number,myfc2_index(i)%jatom_xyz,myfc2_index(i)%phi_temp

            atom1=myfc2_index(i)%iatom_number !first atom index
            IF(atom1.gt.atom_number) CYCLE !if the first atom is not within central unit-cell, not select
            atom2=myfc2_index(i)%jatom_number !second atom index
            direction1=myfc2_index(i)%iatom_xyz !first direction index
            direction2=myfc2_index(i)%jatom_xyz !second direction index

            eff_fc2_terms = eff_fc2_terms + 1
            myfc2_value(atom1,atom2)%phi(direction1,direction2)=myfc2_index(i)%phi_temp !give value to phi(:,:)

            IF(ABS(myfc2_value(atom1,atom2)%phi(direction1,direction2)).gt.max_fc2) THEN
                max_fc2 = ABS(myfc2_value(atom1,atom2)%phi(direction1,direction2))
            END IF

         END DO

         !record all the fc2 that has atom 1 within primitive cell
         ALLOCATE(eff_fc2_index(eff_fc2_terms))
         j = 0
         DO i=1, mx
            atom1=myfc2_index(i)%iatom_number !first atom index
            atom2=myfc2_index(i)%jatom_number !second atom index
            direction1=myfc2_index(i)%iatom_xyz !first direction index
            direction2=myfc2_index(i)%jatom_xyz !second direction index

            IF(atom1.le.atom_number) THEN
                j = j + 1
                eff_fc2_index(j)%group = myfc2_index(i)%group
                eff_fc2_index(j)%iatom_number = atom1
                eff_fc2_index(j)%jatom_number = atom2
                eff_fc2_index(j)%iatom_xyz = direction1
                eff_fc2_index(j)%jatom_xyz = direction2
            END IF

         END DO
!----------------------------------------------------------------------------------------------------
         rank=3
         OPEN(ufc3,file='fc3.dat',status='old',action='read')
         READ(ufc3,'(a)') line
         mx =10000000 ! nterms(rank) has to be initialized
         DO j=1,mx
            READ(ufc3,*,END=93) temporary
         END DO
93       mx=j-1
         REWIND(ufc3)

         fc_terms(rank)=mx

         ALLOCATE(myfc3_index(mx))
!         ALLOCATE(myfc3_value(atom_number,tot_atom_number,tot_atom_number))
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                DO k=1,tot_atom_number
!                    myfc3_value(i,j,k)%psi=0
!                END DO
!            END DO
!         END DO
         !myfc3_value%psi=0
         ALLOCATE(fc3_raw_idx(mx))

         READ(ufc3,'(a)') line
         WRITE(*,'(a)') line
         WRITE(*,*) '********** FCs for rank=3:' ,mx,'***************'

         i=1;j=1;k=1
         DO i=1,mx
            READ(ufc3,*) myfc3_index(i)%junk,myfc3_index(i)%group,myfc3_index(i)%iatom_number,myfc3_index(i)%iatom_xyz,&
            &            myfc3_index(i)%jatom_number,myfc3_index(i)%jatom_xyz,&
            &            myfc3_index(i)%katom_number,myfc3_index(i)%katom_xyz,myfc3_index(i)%psi_temp

            atom1=myfc3_index(i)%iatom_number !first atom index
            fc3_raw_idx(i) = atom1
!            IF(atom1.gt.atom_number) CYCLE
!            atom2=myfc3_index(i)%jatom_number !second atom index
!            atom3=myfc3_index(i)%katom_number
!            direction1=myfc3_index(i)%iatom_xyz !first direction index
!            direction2=myfc3_index(i)%jatom_xyz !second direction index
!            direction3=myfc3_index(i)%katom_xyz

            IF(ABS(myfc3_index(i)%psi_temp).gt.max_fc3) THEN
                max_fc3 = ABS(myfc3_index(i)%psi_temp)
            END IF

         END DO

         CALL unique_sort(fc3_raw_idx, fc3_unique_idx)

         ALLOCATE(myfc3_value(atom_number, SIZE(fc3_unique_idx), SIZE(fc3_unique_idx)))
         i=1;j=1;k=1
         DO i=1,atom_number
         DO j=1,SIZE(fc3_unique_idx)
         DO k=1,SIZE(fc3_unique_idx)
            myfc3_value(i,j,k)%psi = 0d0
         END DO
         END DO
         END DO
         DO i=1, mx
            atom1=myfc3_index(i)%iatom_number !first atom index
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc3_index(i)%jatom_number !second atom index
            atom3=myfc3_index(i)%katom_number
            direction1=myfc3_index(i)%iatom_xyz !first direction index
            direction2=myfc3_index(i)%jatom_xyz !second direction index
            direction3=myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)
                myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)=myfc3_index(i)%psi_temp
            END IF
         END DO
!---------------------------------------------------------------------------------------------------------
!similar here for fc4

         rank=4
         OPEN(ufc4,file='fc4.dat',status='old',action='read')
         READ(ufc4,'(a)') line
         mx =10000000 ! nterms(rank) has to be initialized
         DO j=1,mx
            READ(ufc4,*,END=94) temporary
         END DO
94       mx=j-1
         REWIND(ufc4)

         fc_terms(rank)=mx
!WRITE(*,*) 'Check if anything wrong','mx=',mx,&
!&'tot_atom_number=',tot_atom_number

         ALLOCATE(myfc4_index(mx))
!         ALLOCATE(myfc4_value(atom_number,tot_atom_number,tot_atom_number,tot_atom_number))
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                DO k=1,tot_atom_number
!                    DO l=1,tot_atom_number
!                        myfc4_value(i,j,k,l)%chi=0
!                    END DO
!                END DO
!            END DO
!         END DO
         ALLOCATE(fc4_raw_idx(mx))

         READ(ufc4,'(a)') line
         WRITE(*,'(a)') line
         WRITE(*,*) '********** FCs for rank=4:' ,mx,'***************'

         i=1;j=1;l=1;k=1
         DO i=1,mx
            READ(ufc4,*) myfc4_index(i)%junk,myfc4_index(i)%group,myfc4_index(i)%iatom_number,myfc4_index(i)%iatom_xyz,&
            &            myfc4_index(i)%jatom_number,myfc4_index(i)%jatom_xyz,&
            &            myfc4_index(i)%katom_number,myfc4_index(i)%katom_xyz,&
            &            myfc4_index(i)%latom_number,myfc4_index(i)%latom_xyz,myfc4_index(i)%chi_temp

            atom1=myfc4_index(i)%iatom_number !first atom index
            fc4_raw_idx(i) = atom1

         END DO
         CALL unique_sort(fc4_raw_idx,fc4_unique_idx)
         ALLOCATE(myfc4_value(atom_number,SIZE(fc4_unique_idx),SIZE(fc4_unique_idx),SIZE(fc4_unique_idx)))
         i=1;j=1;l=1;k=1
         DO i=1,atom_number
         DO j=1,SIZE(fc4_unique_idx)
         DO k=1,SIZE(fc4_unique_idx)
         DO l=1,SIZE(fc4_unique_idx)
            myfc4_value(i,j,k,l)%chi = 0d0
         END DO
         END DO
         END DO
         END DO
         DO i=1,mx
            atom1=myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc4_index(i)%jatom_number !second atom index
            atom3=myfc4_index(i)%katom_number
            atom4=myfc4_index(i)%latom_number
            direction1=myfc4_index(i)%iatom_xyz !first direction index
            direction2=myfc4_index(i)%jatom_xyz !second direction index
            direction3=myfc4_index(i)%katom_xyz
            direction4=myfc4_index(i)%latom_xyz
            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(fc4_unique_idx,atom2)
                atom3 = find_loc(fc4_unique_idx,atom3)
                atom4 = find_loc(fc4_unique_idx,atom4)
                myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)=myfc4_index(i)%chi_temp
            END IF
         END DO
!------------------------------------------------------------------------------------------------
         WRITE(*,*) '**********End of check read_force_constants**************'
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)


         CLOSE(ufc2)
         CLOSE(ufc3)
         CLOSE(ufc4)
         !variational_parameters_size(3)=atom_number*d*tot_atom_number*d !necessary for the 1 by 1 correspondence for Broyden


   END SUBROUTINE read_force_constants
!******************************************************************************************************
    SUBROUTINE read_map
    !! The major subroutine to read maps.dat file
        IMPLICIT NONE

        INTEGER :: indie_number,indie_flag,term_flag,i,weight_number
        INTEGER :: linenumber,mapcheck
        CHARACTER(LEN=99) :: line,temporary
        REAL(8) :: weight(2)
        INTEGER :: indie(2)
        INTEGER :: atom1,atom2,direction1,direction2

        IF(.NOT.ALLOCATED(myfc2_group)) ALLOCATE(myfc2_group(d,tot_atom_number,d,tot_atom_number))
        !initialize
        DO direction1=1,d
        DO atom1=1,tot_atom_number
        DO direction2=1,d
        DO atom2=1,tot_atom_number
            myfc2_group(direction1,atom1,direction2,atom2)%group=0
            myfc2_group(direction1,atom1,direction2,atom2)%mat=0d0
        END DO
        END DO
        END DO
        END DO

        atom1=0;atom2=0;direction1=0;direction2=0
        weight=0

        mapcheck=82
        linenumber=0

        OPEN(umap,FILE='maps.dat',STATUS='unknown') !already opened in the main
!        OPEN(mapcheck,FILE='maps_readCheck.dat',STATUS='unknown',ACTION='write')
        !Get the total line number
        DO WHILE(.True.)
            READ(umap,*,END=82) temporary
            linenumber=linenumber+1
        END DO
82 WRITE(*,*) 'READ MAP FILE SUCCESSFULLY: maps.dat total line:',linenumber
        REWIND(umap)

!read the first few lines until reach the '====RANK, GROUP=====' part
        READ(umap,'(A99)') line
        linenumber=linenumber-1
        DO WHILE(.True.)
            IF(line(1:1).eq.'*') THEN

                IF(line(19:19).eq.'2') THEN
                    BACKSPACE(umap)
                    READ(umap,'(A96,I8)') line,indie_number !refer to [other3_nshells.f90] line 1091
!                    WRITE(mapcheck,*) 'RANK 2 GROUP NUMBER: ', indie_number
                    IF(.not.ALLOCATED(indiefc2_index)) ALLOCATE(indiefc2_index(indie_number)) !at most
                END IF

            READ(umap,'(A99)') line !continue read next line
            linenumber=linenumber-1

            ELSEIF(line(1:1).eq.'r') THEN
                READ(umap,'(A99)') line !continue read next line
                linenumber=linenumber-1

            ELSEIF(line(1:2).eq.'  ') THEN !This is for the strange behavior of maps.dat
                READ(umap,'(A99)') line !continue read next line
                linenumber=linenumber-1

            ELSE
                EXIT

            END IF
        END DO


!----------------------------the '====RANK,Group=====' part-----------------------------
!------------------------------check if starts from rank1-------------------------------

        BACKSPACE(umap) !that '====RANK,Group====='line is read by list-directed way, need to be re-read
        linenumber=linenumber+1
!        WRITE(*,*) 'left lines:', linenumber
        DO WHILE(.True.)
            READ(umap,'(A99)') line
            linenumber=linenumber-1
            IF(line(47:47).eq.'1') THEN !'==== line'
                !WRITE(mapcheck,*) 'RANK  ',line(47:47), '  Group  ', line(51:51)
                CYCLE
            ELSEIF(line(2:2).eq.'-') THEN
                IF(line(12:12).eq.'I') THEN
                    !WRITE(mapcheck,*) '# of independent terms:  ', line(40:40)
                    CYCLE
                ELSE
                    !WRITE(mapcheck,*) '# of all terms, MAT: ', line(40:40)
                    CYCLE
                END IF
            ELSEIF(line(2:2).eq.' ') THEN
!                WRITE(mapcheck,*) line
                CYCLE
            ELSE
                IF(line(47:47).eq.'1') THEN
                    !WRITE(mapcheck,*) 'RANK  ',line(47:47), '  Group  ', line(50:50)
                    CYCLE
                ELSE !have read the '=====RANK,GROUP====' line for RANK 2, which basically is line(47:47).eq.'2'
                    EXIT
                END IF
            END IF
        END DO

!---------------------------------the rank2 part-----------------------------------
        indie_number=0
        BACKSPACE(umap) !that '====RANK,Group====='line is read by list-directed way, need to be re-read
        linenumber=linenumber+1
!        WRITE(mapcheck,*) 'left lines:',linenumber
        DO WHILE(linenumber.gt.0)
            READ(umap,'(A99)') line
!            WRITE(mapcheck,*) line
            linenumber=linenumber-1
            IF(line(2:2).eq.'=' .AND. line(47:47).eq.'2') THEN
!            WRITE(mapcheck,*) 'RANK  ',line(47:47), '  Group  ', line(51:51)
            ELSEIF(line(2:2).eq.'-') THEN
                IF(line(12:12).eq.'I') THEN !take record
                    indie_flag=ICHAR(line(40:40))-48 !indie_flag may not exceed 9
                    weight_number=indie_flag
!                    WRITE(mapcheck,*) '# of independent terms:  ', indie_flag !it could be 1 or 2
                    DO WHILE(indie_flag.gt.0)
                        READ(umap,'(A99)') line
                        indie_flag=indie_flag-1
                        linenumber=linenumber-1
                        indie_number=indie_number+1
!                        WRITE(mapcheck,*) line
                        CALL Get_FC2index(line,direction1,atom1,direction2,atom2)

                        indiefc2_index(indie_number)%iatom_xyz=direction1
                        indiefc2_index(indie_number)%iatom_number=atom1
                        indiefc2_index(indie_number)%jatom_xyz=direction2
                        indiefc2_index(indie_number)%jatom_number=atom2
                        indiefc2_index(indie_number)%group=indie_number
                    END DO

                ELSE !line(12:12).eq.'A', categorize
                    BACKSPACE(umap)
                    READ(umap,'(A36,I4)') line,term_flag
!                    WRITE(mapcheck,*) '# of all terms, MAT: ', term_flag !it could be a large integer with multiple digits
                    DO WHILE(term_flag.gt.0)
                        IF(weight_number.eq.1) THEN
                            READ(umap,'(A29,F7.3)') line,weight(1)
                            weight(2)=0d0
                            indie(1)=indie_number
                            indie(2)=0
                        ELSE
                            READ(umap,'(A29,2F7.3)')line,weight(:) !could be 2 different number sqs to 1
                            indie(1)=indie_number-1
                            indie(2)=indie_number
                        END IF
                        term_flag=term_flag-1
                        linenumber=linenumber-1
!                        WRITE(mapcheck,*) term_flag,line!,mat
                        CALL Get_FC2index(line,direction1,atom1,direction2,atom2)
                        myfc2_group(direction1,atom1,direction2,atom2)%group=indie
                        myfc2_group(direction1,atom1,direction2,atom2)%mat=weight
                    END DO
                END IF

            ELSE
                EXIT
            END IF
        END DO

        !***modification, newly added
        DO i=1,SIZE(indiefc2_index)
            atom1 = indiefc2_index(i)%iatom_number
            atom2 = indiefc2_index(i)%jatom_number
            direction1 = indiefc2_index(i)%iatom_xyz
            direction2 = indiefc2_index(i)%jatom_xyz
            indiefc2_index(i)%phi_temp = myfc2_value(atom1,atom2)%phi(direction1,direction2)
        END DO
        !***


!            WRITE(mapcheck,*) 'Total Independent FC2 indexes: ', indie_number
            ifc2_terms=indie_number
            variational_parameters_size(3)=ifc2_terms
!            WRITE(mapcheck,*)'Simple Check',myfc2_group(1,1,2,2)%group
            CLOSE(umap)
!            CLOSE(mapcheck)

        END SUBROUTINE read_map
!******************************************************************************************************
SUBROUTINE Get_FC2index(line,direction1,atom1,direction2,atom2)
!! Convert FC2 indexes back to string format, not used
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: line
    INTEGER :: i,j,temp_digit(5)
    INTEGER,INTENT(OUT) :: atom1,atom2,direction1,direction2

    atom1=0;atom2=0;direction1=0;direction2=0
    temp_digit=0

    i=10 !the first direction 'x/y/z'
    SELECTCASE(line(i:i))
        CASE('x')
            direction1=1
        CASE('y')
            direction1=2
        CASE('z')
            direction1=3
        CASE DEFAULT
            WRITE(*,*)'ERROR when convert string'
            STOP
    ENDSELECT

    j=0
    DO WHILE(line(i+1:i+1).ne.'d')
        i=i+1 !go to next digit
        j=j+1
        temp_digit(j)=ICHAR(line(i:i))-48
        !notice the atom1 index might be multi-digits, but the last digit always followed by a 'd'
        atom1=atom1*10+temp_digit(j) !read one more digit means the original value*10+this digit
    END DO

    i=i+2 !skip the 'd', to the next 'x/y/z'
    SELECTCASE(line(i:i))
        CASE('x')
            direction2=1
        CASE('y')
            direction2=2
        CASE('z')
            direction2=3
        CASE DEFAULT
            WRITE(*,*) 'ERROR when convert string'
            STOP
    ENDSELECT

    j=0
    DO WHILE(line(i+1:i+1).ne.' ')
        i=i+1 !go to next digit
        j=j+1
        temp_digit(j)=ICHAR(line(i:i))-48
        !notice the atom1 index might be multi-digits, but the last digit always followed by a 'd'
        atom2=atom2*10+temp_digit(j) !read one more digit means the original value*10+this digit
    END DO

    !WRITE(*,'(4(A,I6))')'direction1: ',direction1, ', atom1: ', atom1,', direction2: ',direction2,', atom2: ',atom2
END SUBROUTINE Get_FC2index
!------------------------------------------------------------------------------------------------
SUBROUTINE Get_FC2pairs
!!This subroutine is for given initial indiefc2_index(:) calculate the different atom pairs
!!then multiply by 6 to get max possible variational_parameters_size(3)
    IMPLICIT NONE

    TYPE tuple
        INTEGER :: member1,member2
    END TYPE

    INTEGER :: i,atom1,atom2,j,counter
    TYPE(tuple),ALLOCATABLE,DIMENSION(:) :: atom_pairs
    INTEGER :: member1, member2
    LOGICAL :: add
    !max possible different atom pairs <= # of indie fc2
    ALLOCATE(atom_pairs(SIZE(indiefc2_index)))
    !initialize to be 0, 'cause no atom label = 0
    DO j=1,SIZE(atom_pairs)
        atom_pairs(j)%member1 = 0
        atom_pairs(j)%member2 = 0
    END DO
    !iterate through all indie fc2
    counter = 1
    add = .TRUE.
    DO i=1,SIZE(indiefc2_index)
        atom1 = indiefc2_index(i)%iatom_number
        atom2 = indiefc2_index(i)%jatom_number
        DO j=1,SIZE(atom_pairs)
            member1 = atom_pairs(j)%member1
            member2 = atom_pairs(j)%member2
            IF(member1.eq.atom1 .AND. member2.eq.atom2) THEN
                add = .FALSE.
                EXIT
            END IF
        END DO
        IF(add) THEN
            atom_pairs(counter)%member1 = atom1
            atom_pairs(counter)%member2 = atom2
            counter = counter + 1
        ELSE
            add = .TRUE.
            CYCLE
        END IF
    END DO
    ifc2_atompairs = counter-1
!    DO j=1,SIZE(atom_pairs)
!        WRITE(*,*) 'atom1=',atom_pairs(j)%member1,'atom2=',atom_pairs(j)%member2
!    END DO
END SUBROUTINE Get_FC2pairs
!------------------------------------------------------------------------------------------------
SUBROUTINE Extend_indiefc2(ex_indiefc2_index)
!!This subroutine is for at any time, extend the current indiefc2_index to accomodate/fill in
!!all possible indiefc2 directional components, even if they are 0
    IMPLICIT NONE
    TYPE tuple
        INTEGER :: member1,member2
    END TYPE

    INTEGER :: i,atom1,atom2,j,counter,temp,k,direction1,direction2
    TYPE(tuple),ALLOCATABLE,DIMENSION(:) :: atom_pairs
    INTEGER :: member1, member2
    LOGICAL :: add
    !expanded independent fc2 index
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ex_indiefc2_index

    !stage 1: firstly extract all the different atom-pairs from indiefc2_index
    ALLOCATE(atom_pairs(SIZE(indiefc2_index))) !max possible different atom pairs <= # of indie fc2
    !initialize to be 0, 'cause no atom label = 0
    DO j=1,SIZE(atom_pairs)
        atom_pairs(j)%member1 = 0
        atom_pairs(j)%member2 = 0
    END DO
    !iterate through all indie fc2
    counter = 1 !position label
    add = .TRUE. !new pair?  label
    DO i=1,SIZE(indiefc2_index)
        atom1 = indiefc2_index(i)%iatom_number
        atom2 = indiefc2_index(i)%jatom_number
        DO j=1,SIZE(atom_pairs)
            member1 = atom_pairs(j)%member1
            member2 = atom_pairs(j)%member2
            IF(member1.eq.atom1 .AND. member2.eq.atom2) THEN
                add = .FALSE.
                EXIT
            END IF
        END DO
        IF(add) THEN
            atom_pairs(counter)%member1 = atom1
            atom_pairs(counter)%member2 = atom2
            counter = counter + 1
        ELSE
            add = .TRUE.
            CYCLE
        END IF
    END DO

    !stage 2: secondly use all those different atom-pairs in atom_pairs to
    !generate ex_indiefc2_index(:)
    temp = 6*(counter-1)
    ALLOCATE(ex_indiefc2_index(temp))
    DO i=1,SIZE(atom_pairs)
        atom1 = atom_pairs(i)%member1
        atom2 = atom_pairs(i)%member2
        DO j=1,6
            temp = 6*(i-1)+j !current ex_indiefc2_index label
            ex_indiefc2_index(temp)%iatom_number = atom1
            ex_indiefc2_index(temp)%jatom_number = atom2
            SELECTCASE(j)
                CASE(1)
                    ex_indiefc2_index(temp)%iatom_xyz = 1
                    ex_indiefc2_index(temp)%jatom_xyz = 1
                CASE(2)
                    ex_indiefc2_index(temp)%iatom_xyz = 1
                    ex_indiefc2_index(temp)%jatom_xyz = 2
                CASE(3)
                    ex_indiefc2_index(temp)%iatom_xyz = 1
                    ex_indiefc2_index(temp)%jatom_xyz = 3
                CASE(4)
                    ex_indiefc2_index(temp)%iatom_xyz = 2
                    ex_indiefc2_index(temp)%jatom_xyz = 2
                CASE(5)
                    ex_indiefc2_index(temp)%iatom_xyz = 2
                    ex_indiefc2_index(temp)%jatom_xyz = 3
                CASE(6)
                    ex_indiefc2_index(temp)%iatom_xyz = 3
                    ex_indiefc2_index(temp)%jatom_xyz = 3
            ENDSELECT
            !don't forget to assign the group value and most importantly,
            !the fc2 value as well
            !firstly initialize them to be 0
            ex_indiefc2_index(temp)%group = 0
            ex_indiefc2_index(temp)%phi_temp = 0d0 !important
            !then try to find the match in original indiefc2_index
            direction1 = ex_indiefc2_index(temp)%iatom_xyz
            direction2 = ex_indiefc2_index(temp)%jatom_xyz
            DO k=1,SIZE(indiefc2_index)
                IF(atom1.ne.indiefc2_index(k)%iatom_number) CYCLE
                IF(atom2.ne.indiefc2_index(k)%jatom_number) CYCLE
                IF(direction1.ne.indiefc2_index(k)%iatom_xyz) CYCLE
                IF(direction2.ne.indiefc2_index(k)%jatom_xyz) CYCLE
                !found match
                ex_indiefc2_index(temp)%group = indiefc2_index(k)%group
                ex_indiefc2_index(temp)%phi_temp = indiefc2_index(k)%phi_temp !important
            END DO
        END DO
    END DO

    !for check
!    WRITE(34,*)'=====Extend indiefc2_index check====='
!    DO i=1,SIZE(ex_indiefc2_index)
!        WRITE(34,*)'atom1',ex_indiefc2_index(i)%iatom_number,'atom2',ex_indiefc2_index(i)%jatom_number,&
!        &'     ',get_letter(ex_indiefc2_index(i)%iatom_xyz),get_letter(ex_indiefc2_index(i)%jatom_xyz)
!        WRITE(34,*)'fc2 group=',ex_indiefc2_index(i)%group
!        WRITE(34,*)'fc2 value=',ex_indiefc2_index(i)%phi_temp
!        WRITE(34,*)'----------------------------------------------'
!    END DO
!    WRITE(34,*)
!    WRITE(34,*)'=====Original indiefc2_index check====='
!    DO i=1,SIZE(indiefc2_index)
!        WRITE(34,*)'atom1',indiefc2_index(i)%iatom_number,'atom2',indiefc2_index(i)%jatom_number,&
!        &'     ',get_letter(indiefc2_index(i)%iatom_xyz),get_letter(indiefc2_index(i)%jatom_xyz)
!        WRITE(34,*)'fc2 group=',indiefc2_index(i)%group
!        WRITE(34,*)'fc2 value=',indiefc2_index(i)%phi_temp
!        WRITE(34,*)'----------------------------------------------'
!    END DO
END SUBROUTINE Extend_indiefc2
!------------------------------------------------------------------------------------------------
FUNCTION find_indiefc2(atoms,xyzs) RESULT(found)
!!Utility function, given atoms(2),xyzs(2) label,
!!search in current indiefc2_index(:), then output found or not
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: atoms(:), xyzs(:)
    LOGICAL :: found

    INTEGER :: i,atom1,atom2,direction1,direction2


    DO i=1,SIZE(indiefc2_index)
        atom1 = indiefc2_index(i)%iatom_number
        IF(atom1.ne.atoms(1)) CYCLE
        atom2 = indiefc2_index(i)%jatom_number
        IF(atom2.ne.atoms(2)) CYCLE
        direction1 = indiefc2_index(i)%iatom_xyz
        IF(direction1.ne.xyzs(1)) CYCLE
        direction2 = indiefc2_index(i)%jatom_xyz
        IF(direction2.ne.xyzs(2)) CYCLE

        found = .true.
        RETURN
    END DO
    found = .false.
END FUNCTION find_indiefc2
!------------------------------------------------------------------------------------------------
!****************************************Below are not used***************************************
!   SUBROUTINE Force_ASR
!         IMPLICIT NONE
!         TYPE(fc2_value),DIMENSION(:),ALLOCATABLE :: check_fc2
!         TYPE(fc3_value),DIMENSION(:,:),ALLOCATABLE :: check_fc3
!         TYPE(fc4_value),DIMENSION(:,:,:),ALLOCATABLE :: check_fc4
!         INTEGER :: i,j,k,l,m,n,o
!         INTEGER :: new_i, new_j, new_k, new_l, new_m, new_n, new_o
!         INTEGER :: atom1,atom2,atom3,atom4
!         INTEGER :: direction1,direction2,direction3,direction4
!
!!-------------------------------------------------------------------------------------------------------
!!check fc2 list and force ASR
!         ALLOCATE(check_fc2(atom_number))
!         DO i=1,atom_number
!            check_fc2(i)%phi=0
!         END DO
!         DO i=1,fc_terms(2)
!            atom1=myfc2_index(i)%iatom_number
!            IF(atom1.gt.atom_number) cycle
!            atom2=myfc2_index(i)%jatom_number
!            direction1=myfc2_index(i)%iatom_xyz
!            direction2=myfc2_index(i)%jatom_xyz
!            check_fc2(atom1)%phi(direction1,direction2)=check_fc2(atom1)%phi(direction1,direction2)+&
!            &                                           myfc2_value(atom1,atom2)%phi(direction1,direction2)
!            !for a fixed atom1 and direction 1&2, sum up atom2, get all the sum first
!         END DO
!         DO i=1,atom_number
!            DO j=1,d
!                DO k=1,d
!                    IF(check_fc2(i)%phi(j,k).ne.0d0) THEN
!                        myfc2_value(i,i)%phi(j,k)=myfc2_value(i,i)%phi(j,k)-check_fc2(i)%phi(j,k)
!                        check_fc2(i)%phi(j,k)=0
!                    END IF
!                END DO
!            END DO
!         END DO
!         DEALLOCATE(check_fc2)
!!-------------------------------------------------------------------------------------------------------
!!check fc3 list and force ASR
!         ALLOCATE(check_fc3(atom_number,tot_atom_number))
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                check_fc3(i,j)%psi=0
!            END DO
!         END DO
!         DO i=1,fc_terms(3)
!            atom1=myfc3_index(i)%iatom_number
!            IF(atom1.gt.atom_number) cycle
!            atom2=myfc3_index(i)%jatom_number
!            atom3=myfc3_index(i)%katom_number
!            direction1=myfc3_index(i)%iatom_xyz
!            direction2=myfc3_index(i)%jatom_xyz
!            direction3=myfc3_index(i)%katom_xyz
!            atom2=find_loc(fc3_unique_idx,atom2)
!            atom3=find_loc(fc3_unique_idx,atom3)
!            check_fc3(atom1,atom2)%psi(direction1,direction2,direction3)= &
!            &                       check_fc3(atom1,atom2)%psi(direction1,direction2,direction3)+ &
!            &                       myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
!         END DO
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                IF(.NOT.ANY(fc3_unique_idx==j)) THEN
!                    CYCLE
!                ELSE
!                    new_j=find_loc(fc3_unique_idx,j)
!                END IF
!                DO k=1,d
!                    DO l=1,d
!                        DO m=1,d
!                            IF(check_fc3(i,j)%psi(k,l,m).ne.0) THEN
!                                myfc3_value(i,new_j,new_j)%psi(k,l,m)=myfc3_value(i,new_j,new_j)%psi(k,l,m)-&
!                                &check_fc3(i,new_j)%psi(k,l,m)
!                                check_fc3(i,new_j)%psi(k,l,m)=0
!                            END IF
!                        END DO
!                    END DO
!                END DO
!            END DO
!         END DO
!         DEALLOCATE(check_fc3)
!!-------------------------------------------------------------------------------------------------------
!!check fc4 list and force ASR
!         ALLOCATE(check_fc4(atom_number,tot_atom_number,tot_atom_number))
!         !***************
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                DO k=1,tot_atom_number
!                    check_fc4(i,j,k)%chi=0
!                END DO
!            END DO
!         END DO
!         DO i=1,fc_terms(4)
!            atom1=myfc4_index(i)%iatom_number
!            IF(atom1.gt.atom_number) cycle
!            atom2=myfc4_index(i)%jatom_number
!            atom3=myfc4_index(i)%katom_number
!            atom4=myfc4_index(i)%latom_number
!            direction1=myfc4_index(i)%iatom_xyz
!            direction2=myfc4_index(i)%jatom_xyz
!            direction3=myfc4_index(i)%katom_xyz
!            direction4=myfc4_index(i)%latom_xyz
!            atom2=find_loc(fc4_unique_idx,atom2)
!            atom3=find_loc(fc4_unique_idx,atom3)
!            atom4=find_loc(fc4_unique_idx,atom4)
!            check_fc4(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4)=&
!            & check_fc4(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4)+ &
!            & myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
!         END DO
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                IF (.NOT.ANY(fc4_unique_idx==j)) THEN
!                    CYCLE
!                ELSE
!                    new_j=find_loc(fc4_unique_idx,j)
!                END IF
!                DO k=1,tot_atom_number
!                    IF(.NOT.ANY(fc4_unique_idx==k)) THEN
!                        CYCLE
!                    ELSE
!                        new_k=find_loc(fc4_unique_idx,k)
!                    END IF
!                    DO l=1,d
!                        DO m=1,d
!                            DO n=1,d
!                                DO o=1,d
!                                    IF(check_fc4(i,new_j,new_k)%chi(l,m,n,o).ne.0) THEN
!                                        myfc4_value(i,new_j,new_k,new_k)%chi(l,m,n,o)=&
!                                        &myfc4_value(i,new_j,new_k,new_k)%chi(l,m,n,o)-&
!                                        &check_fc4(i,new_j,new_k)%chi(l,m,n,o)
!                                        check_fc4(i,new_j,new_k)%chi(l,m,n,o)=0
!                                    END IF
!                                END DO
!                            END DO
!                        END DO
!                    END DO
!                END DO
!            END DO
!         END DO
!         !***************
!         DO i=1,fc_terms(4)
!            atom1=myfc4_index(i)%iatom_number
!            IF(atom1.gt.atom_number) cycle
!            atom2=myfc4_index(i)%jatom_number
!            atom3=myfc4_index(i)%katom_number
!            atom4=myfc4_index(i)%latom_number
!            direction1=myfc4_index(i)%iatom_xyz
!            direction2=myfc4_index(i)%jatom_xyz
!            direction3=myfc4_index(i)%katom_xyz
!            direction4=myfc4_index(i)%latom_xyz
!            atom2=find_loc(fc4_unique_idx,atom2)
!            atom3=find_loc(fc4_unique_idx,atom3)
!            atom4=find_loc(fc4_unique_idx,atom4)
!            check_fc4(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4)=&
!            & check_fc4(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4)+ &
!            & myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
!         END DO
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                IF (.NOT.ANY(fc4_unique_idx==j)) THEN
!                    CYCLE
!                ELSE
!                    new_j=find_loc(fc4_unique_idx,j)
!                END IF
!                DO k=1,tot_atom_number
!                    IF(.NOT.ANY(fc4_unique_idx==k)) THEN
!                        CYCLE
!                    ELSE
!                        new_k=find_loc(fc4_unique_idx,k)
!                    END IF
!                    DO l=1,d
!                        DO m=1,d
!                            DO n=1,d
!                                DO o=1,d
!                                    IF(check_fc4(i,new_j,new_k)%chi(l,m,n,o).ne.0) THEN
!                                        myfc4_value(i,new_j,new_j,new_k)%chi(l,m,n,o)=&
!                                        &myfc4_value(i,new_j,new_j,new_k)%chi(l,m,n,o)-&
!                                        &check_fc4(i,new_j,new_k)%chi(l,m,n,o)
!                                        check_fc4(i,new_j,new_k)%chi(l,m,n,o)=0
!                                    END IF
!                                END DO
!                            END DO
!                        END DO
!                    END DO
!                END DO
!            END DO
!         END DO
!         !****************
!         DO i=1,fc_terms(4)
!            atom1=myfc4_index(i)%iatom_number
!            IF(atom1.gt.atom_number) cycle
!            atom2=myfc4_index(i)%jatom_number
!            atom3=myfc4_index(i)%katom_number
!            atom4=myfc4_index(i)%latom_number
!            direction1=myfc4_index(i)%iatom_xyz
!            direction2=myfc4_index(i)%jatom_xyz
!            direction3=myfc4_index(i)%katom_xyz
!            direction4=myfc4_index(i)%latom_xyz
!            atom2=find_loc(fc4_unique_idx,atom2)
!            atom3=find_loc(fc4_unique_idx,atom3)
!            atom4=find_loc(fc4_unique_idx,atom4)
!            check_fc4(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4)=&
!            & check_fc4(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4)+ &
!            & myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
!         END DO
!         DO i=1,atom_number
!            DO j=1,tot_atom_number
!                IF (.NOT.ANY(fc4_unique_idx==j)) THEN
!                    CYCLE
!                ELSE
!                    new_j=find_loc(fc4_unique_idx,j)
!                END IF
!                DO k=1,tot_atom_number
!                    IF(.NOT.ANY(fc4_unique_idx==k)) THEN
!                        CYCLE
!                    ELSE
!                        new_k=find_loc(fc4_unique_idx,k)
!                    END IF
!                    DO l=1,d
!                        DO m=1,d
!                            DO n=1,d
!                                DO o=1,d
!                                    IF(check_fc4(i,new_j,new_k)%chi(l,m,n,o).ne.0) THEN
!                                        myfc4_value(i,new_k,new_j,new_k)%chi(l,m,n,o)=&
!                                        &myfc4_value(i,new_k,new_j,new_k)%chi(l,m,n,o)-&
!                                        &check_fc4(i,new_j,new_k)%chi(l,m,n,o)
!                                        check_fc4(i,new_j,new_k)%chi(l,m,n,o)=0
!                                    END IF
!                                END DO
!                            END DO
!                        END DO
!                    END DO
!                END DO
!            END DO
!         END DO
!         DEALLOCATE(check_fc4)
!   END SUBROUTINE Force_ASR
!******************************************************************************************************************
   SUBROUTINE Symmetrize_FCs
    !!Symmetrize FC2s based on indexes F_ij^xy = F_ji^yx

        IMPLICIT NONE
        INTEGER :: x,y,z,w

        ALLOCATE(myefc2_value(atom_number,tot_atom_number))
        ALLOCATE(myefc3_value(atom_number,tot_atom_number,tot_atom_number))
        ALLOCATE(myefc4_value(atom_number,tot_atom_number,tot_atom_number,tot_atom_number))


        DO x=1,d
        DO y=1,d
            myefc2_value%phi(x,y)=0.5*(myfc2_value%phi(x,y)+myfc2_value%phi(y,x))
        DO z=1,d
            myefc3_value%psi(x,y,z)=0.5*(myfc3_value%psi(x,y,z)+myfc3_value%psi(y,x,z))
        DO w=1,d
            myefc4_value%chi(x,y,z,w)=0.5*(myfc4_value%chi(x,y,z,w)+myfc4_value%chi(y,x,z,w))
        END DO !w
        END DO !z
        END DO !y
        END DO !x
        !deallocate myfc#_value arrays

    END SUBROUTINE Symmetrize_FCs
!****************************************************************************************************************
 subroutine read_fc23
  !! legacy subroutinie that read fc2.dat and fc3.dat, then store the info in
 !! legacy variables
 use io2
 implicit none
 character line*90
 integer t,rank,res,j,mx

  open(ufc2,file='fc2.dat'    ,status='old')
  open(ufc3,file='fc3.dat'    ,status='old')

! read(ufc2,'(a)') line
! read(ufc,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
! read(ufc,'(a)') line
! read(ufc,*)nterms(1),nterms(2),nterms(3),nterms(4)
! read(ufc,'(a)') line
! read(ufc,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)
 res = 0
!----------------------------------------
! rank=1
! if ( include_fc(rank) .eq. 1 ) then
! mx = nterms(rank)
! allocate(iatomterm_1(rank,mx),ixyzterm_1(rank,mx),igroup_1(mx), &
! & ampterm_1(mx),fcs_1(ngroups(rank)))
! read(ufc,*)line
! do j=1,nterms(rank)
!       read(ufc,*)t,igroup_1(t), &
!& iatomterm_1(1,t),ixyzterm_1(1,t),  &
!& fcs_1(igroup_1(t)),ampterm_1(t)
! enddo
! res = res + igroup_1(nterms(rank))
! endif
!----------------------------------------
 rank=2
! if ( include_fc(rank) .eq. 1 ) then
 mx = nterms(rank)
 read(ufc2,'(a)') line
 do j=1,mx
    read(ufc2,*,end=92)t
 enddo
92 mx=j-1
 rewind(ufc2)
 nterms(rank)=mx
 write(ulog,*)'Rank=2, reading nterms(2)=',mx,' terms'

 IF(ALLOCATED(iatomterm_2)) DEALLOCATE(iatomterm_2)
 IF(ALLOCATED(ixyzterm_2)) DEALLOCATE(ixyzterm_2)
 IF(ALLOCATED(igroup_2)) DEALLOCATE(igroup_2)
 allocate(iatomterm_2(rank,mx),ixyzterm_2(rank,mx),igroup_2(mx)) !,ampterm_2(mx))

 IF(ALLOCATED(fcs_2)) DEALLOCATE(fcs_2)
 IF(ALLOCATED(grun_fc)) DEALLOCATE(grun_fc)
 allocate(fcs_2(mx),grun_fc(mx))
 read(ufc2,'(a)') line
! write(*  ,'(a)') line
! write(*,*) '********** FCs for rank=2',mx

 do j=1,mx
       read(ufc2,*)t,igroup_2(j), &
& iatomterm_2(1,j),ixyzterm_2(1,j),  &
& iatomterm_2(2,j),ixyzterm_2(2,j),  &
!& fcs_2(igroup_2(j)),ampterm_2(j)
& fcs_2(j)
 enddo
 res = res + igroup_2(nterms(rank))
! endif
!----------------------------------------
 rank=3
 if ( include_fc(rank) .eq. 1 ) then
 mx = nterms(rank)
! write(*,*) '********** FCs for rank=3',mx
 read(ufc3,'(a)') line
 do j=1,nterms(rank)
       read(ufc3,*,end=93)t
 enddo
93 mx=j-1
 rewind(ufc3)
 nterms(rank)=mx
 write(ulog,*)'Rank=3, reading nterms(3)=',mx,' terms'

 IF(ALLOCATED(iatomterm_3)) DEALLOCATE(iatomterm_3)
 IF(ALLOCATED(ixyzterm_3)) DEALLOCATE(ixyzterm_3)
 IF(ALLOCATED(igroup_3)) DEALLOCATE(igroup_3)
 allocate(iatomterm_3(rank,mx),ixyzterm_3(rank,mx),igroup_3(mx) ) !,ampterm_3(mx),fcs_3(ngroups(rank)))

IF(ALLOCATED(fcs_3)) DEALLOCATE(fcs_3)
 allocate(fcs_3(mx))
 read(ufc3,'(a)') line
!write(*  ,'(a)') line
 do j=1,mx
!       read(ufc3,*)t,igroup_3(t), &
!& iatomterm_3(1,t),ixyzterm_3(1,t),  &
!& iatomterm_3(2,t),ixyzterm_3(2,t),  &
!& iatomterm_3(3,t),ixyzterm_3(3,t),  &
!& fcs_3(igroup_3(t)),ampterm_3(t)
       read(ufc3,*)t,igroup_3(j), &
& iatomterm_3(1,j),ixyzterm_3(1,j),  &
& iatomterm_3(2,j),ixyzterm_3(2,j),  &
& iatomterm_3(3,j),ixyzterm_3(3,j),  &
& fcs_3(j)
!& fcs_3(igroup_3(j)),ampterm_3(j)
! write(*,*) j,fcs_3(j)
 enddo
 res = res + igroup_3(nterms(rank))
 endif
 write(ulog,*)'READ_FC23: done!, number of groups read is=',res
 write(ulog,*)'READ_FC23:',fcs_2

 close(ufc2)
 close(ufc3)
 end subroutine read_fc23
!----------------------------------------------------------------------------------------
 SUBROUTINE check_read_fcs(rnk)
  !!mainly to check the effects of <fix_asr_fc#> subroutines have on the fcs
 !!show the difference between 'original value' and 'after fix value' of fcs
    IMPLICIT NONE
    INTEGER,INTENT(in) :: rnk
    INTEGER :: atom1,atom2,atom3,atom4
    INTEGER :: xyz1,xyz2,xyz3,xyz4
    INTEGER :: i
    INTEGER :: new_atom2,new_atom3,new_atom4
    OPEN(71,FILE='check_readfcs.dat',STATUS='unknown',POSITION='append',ACTION='write')
    !check myfc3_value
    IF(rnk.eq.3) THEN
        DO i=1,SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number
            xyz1 = myfc3_index(i)%iatom_xyz
            xyz2 = myfc3_index(i)%jatom_xyz
            xyz3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                new_atom2 = find_loc(fc3_unique_idx,atom2)
                new_atom3 = find_loc(fc3_unique_idx,atom3)
                IF(ABS(myfc3_value(atom1,new_atom2,new_atom3)%psi(xyz1,xyz2,xyz3)&
                &-myfc3_index(i)%psi_temp).ge.1d-8) THEN
                    WRITE(71,*)'found one! ',i
                    WRITE(71,*)'read: ',myfc3_value(atom1,new_atom2,new_atom3)%psi(xyz1,xyz2,xyz3)
                    WRITE(71,*)'original: ',myfc3_index(i)%psi_temp
                END IF
            END IF
        END DO
    END IF

    !check myfc4_value
    IF(rnk.eq.4) THEN
        DO i=1,SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number

            xyz1 = myfc4_index(i)%iatom_xyz
            xyz2 = myfc4_index(i)%jatom_xyz
            xyz3 = myfc4_index(i)%katom_xyz
            xyz4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                new_atom2 = find_loc(fc4_unique_idx,atom2)
                new_atom3 = find_loc(fc4_unique_idx,atom3)
                new_atom4 = find_loc(fc4_unique_idx,atom4)
                IF(ABS(myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)&
                &-myfc4_index(i)%chi_temp).ge.1d-8) THEN
                    WRITE(71,*)'found one! ',i
                    WRITE(71,*)'read: ',myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)
                    WRITE(71,*)'original: ',myfc4_index(i)%chi_temp
                END IF
            END IF
        END DO
    END IF

    CLOSE(71)
 END SUBROUTINE check_read_fcs
!----------------------------------------------------------------------------------------
SUBROUTINE check_read_fc2
!!a subroutine to check <fc2.dat> read correctly
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: atom1,atom2
    INTEGER :: direction1, direction2
    OPEN(24,FILE='check_read_fc2.dat')
    WRITE(24,*) 'term,group,atom1,xyz1,atom2,xyz2,fc2_value'
    DO i=1,SIZE(myfc2_index)
        atom1 = myfc2_index(i)%iatom_number
        atom2 = myfc2_index(i)%jatom_number
        direction1 = myfc2_index(i)%iatom_xyz
        direction2 = myfc2_index(i)%jatom_xyz
        IF(atom1.le.atom_number) THEN
            WRITE(24,7) i,myfc2_index(i)%group,atom1,direction1,atom2,direction2,&
            & myfc2_value(atom1,atom2)%phi(direction1,direction2)
        ELSE
            WRITE(24,8) i,myfc2_index(i)%group,atom1,direction1,atom2,direction2,&
            & 'not used'
        ENDIF
    END DO
7 format(2(i6),2(i7,i2),3x,e14.8)
8 format(2(i6),2(i7,i2),3x,a8)
    CLOSE(24)
END SUBROUTINE check_read_fc2
!----------------------------------------------------------------------------------------
SUBROUTINE check_trans_fcs(rnk)
!!check if translational symmetry is still satisfied or not
    IMPLICIT NONE
    INTEGER,INTENT(in) :: rnk
    INTEGER :: atom1,atom2,atom3,atom4
    INTEGER :: xyz1,xyz2,xyz3,xyz4
    INTEGER :: i
    INTEGER :: new_atom2,new_atom3,new_atom4
    OPEN(72,FILE='check_transfcs.dat',STATUS='unknown',POSITION='append',ACTION='write')
    !check FC3
    IF(rnk.eq.3) THEN
        DO i=1,SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number
            xyz1 = myfc3_index(i)%iatom_xyz
            xyz2 = myfc3_index(i)%jatom_xyz
            xyz3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                new_atom2 = find_loc(fc3_unique_idx,atom2)
                new_atom3 = find_loc(fc3_unique_idx,atom3)
                IF(ABS(myfc3_value(atom1,new_atom2,new_atom3)%psi(xyz1,xyz2,xyz3)&
                &-myfc3_value(atom1,new_atom3,new_atom2)%psi(xyz1,xyz3,xyz2)).ge.1d-8) THEN
                    WRITE(71,*)'FC3 found one! ',i
                    WRITE(71,*)'read: ',myfc3_value(atom1,new_atom2,new_atom3)%psi(xyz1,xyz2,xyz3)
                    WRITE(71,*)'swapped 23: ',myfc3_value(atom1,new_atom3,new_atom2)%psi(xyz1,xyz3,xyz2)
                END IF
            END IF
        END DO
    END IF

    !check FC4
    IF(rnk.eq.4) THEN
        DO i=1,SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number

            xyz1 = myfc4_index(i)%iatom_xyz
            xyz2 = myfc4_index(i)%jatom_xyz
            xyz3 = myfc4_index(i)%katom_xyz
            xyz4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                new_atom2 = find_loc(fc4_unique_idx,atom2)
                new_atom3 = find_loc(fc4_unique_idx,atom3)
                new_atom4 = find_loc(fc4_unique_idx,atom4)
                IF(ABS(myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)&
                &-myfc4_value(atom1,new_atom3,new_atom2,new_atom4)%chi(xyz1,xyz3,xyz2,xyz4)).ge.1d-8) THEN
                    WRITE(71,*)'FC4 found one! ',i
                    WRITE(71,*)'read: ',myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)
                    WRITE(71,*)'swapped 23: ',myfc4_value(atom1,new_atom3,new_atom2,new_atom4)%chi(xyz1,xyz3,xyz2,xyz4)
                ELSEIF (ABS(myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)&
                &-myfc4_value(atom1,new_atom2,new_atom4,new_atom3)%chi(xyz1,xyz2,xyz4,xyz3)).ge.1d-8) THEN
                    WRITE(71,*)'FC4 found one! ',i
                    WRITE(71,*)'read: ',myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)
                    WRITE(71,*)'swapped 34: ',myfc4_value(atom1,new_atom2,new_atom4,new_atom3)%chi(xyz1,xyz2,xyz4,xyz3)
                ELSEIF (ABS(myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)&
                &-myfc4_value(atom1,new_atom4,new_atom3,new_atom2)%chi(xyz1,xyz4,xyz3,xyz2)).ge.1d-8) THEN
                    WRITE(71,*)'FC4 found one! ',i
                    WRITE(71,*)'read: ',myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(xyz1,xyz2,xyz3,xyz4)
                    WRITE(71,*)'swapped 24: ',myfc4_value(atom1,new_atom4,new_atom3,new_atom2)%chi(xyz1,xyz4,xyz3,xyz2)
                END IF
            END IF
        END DO
    END IF
    WRITE(72,*)'==============================================================='

    CLOSE(72)
END SUBROUTINE check_trans_fcs
!-----------------------------------------------------------------------------------
 SUBROUTINE GetInversePhi
  !! Get inverse matrix inv_phiiTau for FC2, used for guess start calculation
    IMPLICIT NONE
    INTEGER :: i,j
    INTEGER :: atom1,atom2,xyz1,xyz2
    INTEGER :: n !matrix rank,n=3*number of tau
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: a,b,check !matrix for phiTau and inv_phiTau
    n = d*atom_number
    ALLOCATE(a(n,n),b(n,n))
    IF(.not.allocated(inv_phiTau)) ALLOCATE(inv_phiTau(atom_number,atom_number))
    !make phiTau into matrix a
    DO i=1,n
        DO j=1,n
            atom1 = INT((i+d-1)/d)
            atom2 = INT((j+d-1)/d)
            xyz1 = MOD((i+d-1),d)+1
            xyz2 = MOD((j+d-1),d)+1
            a(i,j) = myfc2_value(atom1,atom2)%phi(xyz1,xyz2)
        END DO
    END DO
    CALL invers_r(a,b,n)

    !save matrix b into inv_phiTau
    DO i=1,n
        DO j=1,n
            atom1 = INT((i+d-1)/d)
            atom2 = INT((j+d-1)/d)
            xyz1 = MOD((i+d-1),d)+1
            xyz2 = MOD((j+d-1),d)+1
            inv_phiTau(atom1,atom2)%phi(xyz1,xyz2) = b(i,j)
        END DO
    END DO

    !check
!    ALLOCATE(check(n,n))
!    DO i=1,n
!        DO j=1,n
!            atom1 = INT((i+d-1)/d)
!            atom2 = INT((j+d-1)/d)
!            xyz1 = MOD((i+d-1),d)+1
!            xyz2 = MOD((j+d-1),d)+1
!            a(i,j) = myfc2_value(atom1,atom2)%phi(xyz1,xyz2)
!            b(i,j) = inv_phiTau(atom1,atom2)%phi(xyz1,xyz2)
!        END DO
!    END DO
!
!    DO i=1,n
!        DO j=1,n
!            check(i,j) = a(i,:).dot.b(:,j)
!        END DO
!        WRITE(*,*)'check(',i,',:)=',check(i,:)
!    END DO
!    DEALLOCATE(check)

    DEALLOCATE(a,b)
 END SUBROUTINE GetInversePhi
 !-----------------------------------------------------------------------------------
 SUBROUTINE GetElastic_simple
  !!approximate elastic constant from simple formula, using initial fc2
    IMPLICIT NONE
    INTEGER :: i,j,al,be,ga,de
    INTEGER :: v1, v2, voigt
    REAL(8) :: Rij_ga,Rij_de,cell_volume
    REAL(8),DIMENSION(6,6) :: temp !call inverse matrix will destroy original matrix, so

    CALL calculate_volume(r1,r2,r3,cell_volume)

    elastic = 0d0;compliance = 0d0!don't forget to initialize as 0

    !4 nested xyz loop
    DO al=1,3
    DO be=1,3
    DO ga=1,3
    DO de=1,3

    v1 = voigt(al,be)
    v2 = voigt(ga,de)

    !atom ij sum
    DO i=1,atom_number
    DO j=1,tot_atom_number
        Rij_ga = every_atom(j)%R(ga)+every_atom(j)%tau(ga)-every_atom(i)%R(ga)-every_atom(i)%tau(ga)
        Rij_de = every_atom(j)%R(de)+every_atom(j)%tau(de)-every_atom(i)%R(de)-every_atom(i)%tau(de)
        elastic(v1,v2) = elastic(v1,v2) + myfc2_value(i,j)%phi(al,be)*Rij_ga*Rij_de/cell_volume
    END DO !atom i loop
    END DO !atom j loop

    END DO !de loop
    END DO !ga loop
    END DO !be loop
    END DO !al loop

    !get inverse elastic
    temp = elastic
    CALL invers_r(temp, compliance,6)

 END SUBROUTINE GetElastic_simple
!========================================================================================
    SUBROUTINE print_indieFC2
    !! print all the irreducible FC2
        IMPLICIT NONE
        INTEGER :: i

        OPEN(117,FILE='print_indieFC2.txt',STATUS='unknown')
        DO i=1, SIZE(indiefc2_index)
            WRITE(117,*) i,',',indiefc2_index(i)%phi_temp
        END DO
        CLOSE(117)
    END SUBROUTINE print_indieFC2
!========================================================================================
SUBROUTINE find_largest_fc2atomidx
    !this subroutine is used for check the largest(farthest) atom label of fc2
    IMPLICIT NONE
    INTEGER :: i,largest

    largest = 1
    DO i=1, SIZE(myfc2_index)
        largest = MAX(largest, myfc2_index(i)%iatom_number, myfc2_index(i)%jatom_number)
    END DO
    WRITE(*,*) 'largest atom index for fc2 data:', largest
END SUBROUTINE find_largest_fc2atomidx
!========================================================================================
END MODULE DFT_force_constants
