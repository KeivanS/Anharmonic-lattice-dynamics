!unit_number2=34

MODULE force_update
    USE VA_math
    REAL(8) R_0
    REAL(8) :: resolution=1e-4

    CONTAINS
!==========================================================================================================
    SUBROUTINE all_fc_update
        IMPLICIT NONE
        INTEGER :: rnk, i, j
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER, ALLOCATABLE,DIMENSION(:) :: atoms, xyzs

INTEGER :: new_atom1, new_atom2, new_atom3
INTEGER,ALLOCATABLE,DIMENSION(:) :: directions
REAL(8) :: temp

        ![prepare] alloc & re-alloc, remap indexes
        CALL prepare_fc4
        WRITE(34,*)"fc4 successfully prepared!"
        CALL prepare_fc3
        WRITE(34,*)"fc3 successfully prepared!"
        CALL prepare_fc2
        WRITE(34,*)"fc2 successfully prepared!"
        ![update] fc value one by one, which is low efficient
        !input args should be NEW atom_idx
        !so the best way is to loop through map()

        !NOTICE: remember to update myfc#_index(i) after value update
        !-------------------fc4-------------------------
        rnk = 4
        ALLOCATE(atoms(rnk),xyzs(rnk))
        DO j=1, map(rnk)%ngr
            DO i=1, map(rnk)%nt(j)
                atoms(1) = map(rnk)%gr(j)%iat(1,i)
                IF(atoms(1).gt.atom_number) CYCLE
                atoms(2) = map(rnk)%gr(j)%iat(2,i)
                atoms(3) = map(rnk)%gr(j)%iat(3,i)
                atoms(4) = map(rnk)%gr(j)%iat(4,i)

                xyzs(1) = map(rnk)%gr(j)%ixyz(1,i)
                xyzs(2) = map(rnk)%gr(j)%ixyz(2,i)
                xyzs(3) = map(rnk)%gr(j)%ixyz(3,i)
                xyzs(4) = map(rnk)%gr(j)%ixyz(4,i)
                CALL uni_fc_update(rnk, atoms, xyzs)
            END DO
        END DO
        DEALLOCATE(atoms,xyzs)

        DO i=1, SIZE(myfc4_index)
            atom1=myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc4_index(i)%jatom_number
            atom3=myfc4_index(i)%katom_number
            atom4=myfc4_index(i)%latom_number
            direction1=myfc4_index(i)%iatom_xyz
            direction2=myfc4_index(i)%jatom_xyz
            direction3=myfc4_index(i)%katom_xyz
            direction4=myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(fc4_unique_idx,atom2)
                atom3 = find_loc(fc4_unique_idx,atom3)
                atom4 = find_loc(fc4_unique_idx,atom4)
                myfc4_index(i)%chi_temp = &
                &myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
            ELSE
                myfc4_index(i)%chi_temp = 0d0
            END IF
        END DO
WRITE(34,*)"fc4 update finished!"
        !------------------fc3-------------------------
        rnk = 3
        ALLOCATE(atoms(rnk),xyzs(rnk))
        DO j=1, map(rnk)%ngr
            DO i=1, map(rnk)%nt(j)
                atoms(1) = map(rnk)%gr(j)%iat(1,i)
                IF(atoms(1).gt.atom_number) CYCLE
                atoms(2) = map(rnk)%gr(j)%iat(2,i)
                atoms(3) = map(rnk)%gr(j)%iat(3,i)

                xyzs(1) = map(rnk)%gr(j)%ixyz(1,i)
                xyzs(2) = map(rnk)%gr(j)%ixyz(2,i)
                xyzs(3) = map(rnk)%gr(j)%ixyz(3,i)

                CALL uni_fc_update(rnk, atoms, xyzs)
            END DO
        END DO
        DEALLOCATE(atoms,xyzs)

!----------------------------additional check----------------------------------
!IF(.NOT.ALLOCATED(atoms)) ALLOCATE(atoms(3))
!IF(.NOT.ALLOCATED(directions)) ALLOCATE(directions(3))
!OPEN(373, FILE='additional_checkfc3.txt',STATUS='unknown')
!WRITE(373,*)'=====term,(atom,xyz)*3,fc3_value difference====='
!DO i=1,SIZE(myfc3_index)
!    atoms(1) = myfc3_index(i)%iatom_number
!    IF(atoms(1).gt.atom_number) CYCLE
!    atoms(2) = myfc3_index(i)%jatom_number
!    atoms(3) = myfc3_index(i)%katom_number
!    directions(1) = myfc3_index(i)%iatom_xyz
!    directions(2) = myfc3_index(i)%jatom_xyz
!    directions(3) = myfc3_index(i)%katom_xyz
!
!    !the 'if' condition clause below should always be satisfied
!    IF(ANY(fc3_unique_idx==atoms(2)).AND.ANY(fc3_unique_idx==atoms(3))) THEN
!        atom1=find_loc(fc3_unique_idx,atoms(1))
!        atom2=find_loc(fc3_unique_idx,atoms(2))
!        atom3=find_loc(fc3_unique_idx,atoms(3))
!        temp = myfc3_value(atom1,atom2,atom3)%psi(directions(1),directions(2),directions(3))
!        temp = temp - myfc3_index(i)%psi_temp
!        IF(ABS(temp).gt.0d0) THEN
!        WRITE(373,*)i,atoms(1),get_letter(directions(1)),atoms(2),get_letter(directions(2)),&
!        &atoms(3),get_letter(directions(3)),temp
!        END IF
!    END IF
!END DO
!DEALLOCATE(atoms,directions)
!CLOSE(373)
!STOP
!------------------------------------------------------------------------------
WRITE(34,*)"fc3 update finished!"
!------------------------------------------------------------------------------
        DO i=1, SIZE(myfc3_index)
            atom1=myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc3_index(i)%jatom_number
            atom3=myfc3_index(i)%katom_number
            direction1=myfc3_index(i)%iatom_xyz
            direction2=myfc3_index(i)%jatom_xyz
            direction3=myfc3_index(i)%katom_xyz

            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)
                myfc3_index(i)%psi_temp = &
                &myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
            ELSE
                myfc3_index(i)%psi_temp = 0d0
            END IF
        END DO

        !------------------fc2-------------------------
        rnk = 2
        ALLOCATE(atoms(rnk),xyzs(rnk))
        DO j=1, map(rnk)%ngr
            DO i=1, map(rnk)%nt(j)
                atoms(1) = map(rnk)%gr(j)%iat(1,i)
                IF(atoms(1).gt.atom_number) CYCLE
                atoms(2) = map(rnk)%gr(j)%iat(2,i)

                xyzs(1) = map(rnk)%gr(j)%ixyz(1,i)
                xyzs(2) = map(rnk)%gr(j)%ixyz(2,i)

                CALL uni_fc_update(rnk, atoms, xyzs)
            END DO
        END DO
        DEALLOCATE(atoms,xyzs)

!----------------------------additional check----------------------------------
WRITE(34,*)"Now perform additional check for fc2"
OPEN(302, FILE='additional_checkfc2.txt')
WRITE(302,*)'=====term,(atom,xyz)*2,fc2_value difference====='
DO i=1,SIZE(myfc2_index)
    atom1 = myfc2_index(i)%iatom_number
    IF(atom1.gt.atom_number) CYCLE
    atom2 = myfc2_index(i)%jatom_number
    direction1 = myfc2_index(i)%iatom_xyz
    direction2 = myfc2_index(i)%jatom_xyz
    temp = myfc2_value(atom1,atom2)%phi(direction1,direction2)
    temp = temp - myfc2_index(i)%phi_temp
    WRITE(302,*)i,atom1,get_letter(direction1),atom2,get_letter(direction2),&
    &ABS(temp)
END DO
CLOSE(302)
!------------------------------------------------------------------------------
        DO i=1, SIZE(myfc2_index)
            atom1=myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc2_index(i)%jatom_number
            direction1=myfc2_index(i)%iatom_xyz
            direction2=myfc2_index(i)%jatom_xyz
            myfc2_index(i)%phi_temp = myfc2_value(atom1,atom2)%phi(direction1,direction2)
        END DO
        WRITE(34,*)"All FCs Updated"

        WRITE(*,*)"=======fc2 update checked======="
        CALL fc3_update_check
        CALL fc4_update_check
    END SUBROUTINE all_fc_update
!---------------------------------------------------
    SUBROUTINE prepare_fc2
        IMPLICIT NONE
        INTEGER :: i, j, idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2
        INTEGER :: region1, region2, region3
        INTEGER :: tot_indie_fc2, ntindp, tot_fc2
        INTEGER,DIMENSION(:), ALLOCATABLE :: atoms, directions

        OPEN(362,FILE='fc2_update_check.txt',STATUS='unknown',ACTION = 'write')
        WRITE(362,*)'===== # RANK 2 tensors :term,group,(iatom,ixyz)_2 ====='

        ![1] Store old fc2_group and old indiefc2_index info
        IF(ALLOCATED(old_indiefc2_index)) DEALLOCATE(old_indiefc2_index)
        ALLOCATE(old_indiefc2_index,source=indiefc2_index)
        IF(ALLOCATED(old_myfc2_group)) DEALLOCATE(old_myfc2_group)
        ALLOCATE(old_myfc2_group, source=myfc2_group)

        ![2] Remap 'indiefc2_index' info from <setup_maps>
        tot_indie_fc2 = 0
        DO i = 1, map(2)%ngr
            tot_indie_fc2 = tot_indie_fc2 + map(2)%ntind(i)
        END DO
        IF(ALLOCATED(indiefc2_index)) DEALLOCATE(indiefc2_index)
        ALLOCATE(indiefc2_index(tot_indie_fc2))
        !'tot_indie_fc2' is the 'indie_number' in <read_map>
        ifc2_terms = tot_indie_fc2
        variational_parameters_size(3) = ifc2_terms
        IF(ALLOCATED(myfc2_group)) DEALLOCATE(myfc2_group)
        ALLOCATE(myfc2_group(d,tot_atom_number,d,tot_atom_number))
        WRITE(362,*)'================new indiefc2 got from <setup_maps>========================'

        idx = 0
        DO j = 1, map(2)%ngr
            DO i = 1, map(2)%ntind(j)
                idx = idx + 1
                indiefc2_index(idx)%group = idx

                atom1 = map(2)%gr(j)%iatind(1,i)
                atom2 = map(2)%gr(j)%iatind(2,i)
                direction1 = map(2)%gr(j)%ixyzind(1,i)
                direction2 = map(2)%gr(j)%ixyzind(2,i)

                indiefc2_index(idx)%iatom_number = atom1
                indiefc2_index(idx)%iatom_xyz = direction1
                indiefc2_index(idx)%jatom_number = atom2
                indiefc2_index(idx)%jatom_xyz = direction2

                myfc2_group(direction1,atom1,direction2,atom2)%group = (/j,0/)
                myfc2_group(direction1,atom1,direction2,atom2)%mat = (/1d0,0d0/)
                !------------------------------------------------------------
                WRITE(362,*)idx
                WRITE(362,*)'atom1,direction1,atom2,direction2'
                WRITE(362,*)atom1,get_letter(direction1),atom2,get_letter(direction2)
                !------------------------------------------------------------

            END DO
        END DO

        ![3] Remap 'myfc2_group'
        region2 = 0
        region3 = 0
        DO j = 1, map(2)%ngr
            ntindp = map(2)%ntind(j)
            DO i = 1, map(2)%nt(j)
                atom1 = map(2)%gr(j)%iat(1,i)
                direction1 = map(2)%gr(j)%ixyz(1,i)
                atom2 = map(2)%gr(j)%iat(2,i)
                direction2 = map(2)%gr(j)%ixyz(2,i)
                IF(ntindp.eq.1) THEN
                    myfc2_group(direction1,atom1,direction2,atom2)%group = (/j,0/)
                    myfc2_group(direction1,atom1,direction2,atom2)%mat = (/map(2)%gr(j)%mat(i,1),0d0/)
                ELSE
                    myfc2_group(direction1,atom1,direction2,atom2)%group = (/j-1,j/)
                    myfc2_group(direction1,atom1,direction2,atom2)%mat&
                    & = (/map(2)%gr(j)%mat(i,1),map(2)%gr(j)%mat(i,2)/)
                END IF

                WRITE(362,*)'-----fc2 generated by <setup_maps>-----'
                WRITE(362,*)'atom1, xyz1, atom2, xyz2'
                WRITE(362,*)atom1,get_letter(direction1),atom2,get_letter(direction2)
                atoms = (/atom1,atom2/)
                directions = (/direction1,direction2/)
                IF(findAtom_inRegion12(2,atoms,directions)) THEN
                    region2 = region2 + 1
                    WRITE(362,*) 'Found by remapping to old ones: this fc2 existed in the last iteration'
                ELSE
                    region3 = region3 + 1
                    WRITE(362,*) 'Not found by remapping to old ones'
                END IF

            END DO
        END DO

        ![4] Store previous 'myfc2_value' and 'myfc2_index'
        !also reallocate them based on the new size
        !because tot#of fc2 may increase, tot#of atoms may increase
        IF(ALLOCATED(oldfc2_index)) DEALLOCATE(oldfc2_index)
        ALLOCATE(oldfc2_index, source=myfc2_index)
        IF(ALLOCATED(prev_fc2)) DEALLOCATE(prev_fc2)
        ALLOCATE(prev_fc2, source=myfc2_value)

        tot_fc2 = 0
        DO j=1, map(2)%ngr
            tot_fc2 = tot_fc2 + map(2)%nt(j)
        END DO
        IF(ALLOCATED(myfc2_index)) DEALLOCATE(myfc2_index)
        ALLOCATE(myfc2_index(tot_fc2))

        IF(ALLOCATED(myfc2_value)) DEALLOCATE(myfc2_value)
        ALLOCATE(myfc2_value(atom_number,tot_atom_number))
        DO i=1,atom_number
            DO j=1,tot_atom_number
                myfc2_value(i,j)%phi=0
            END DO
         END DO
        ![5] Update 'myfc2_index' by remapping from 'oldfc2_index'
        !also initialize by corresponding old value
        idx = 0
        DO j = 1, map(2)%ngr
            DO i = 1, map(2)%nt(j)
                idx = idx + 1

                new_atom1 = map(2)%gr(j)%iat(1,i)
                new_atom2 = map(2)%gr(j)%iat(2,i)
                direction1 = map(2)%gr(j)%ixyz(1,i)
                direction2 = map(2)%gr(j)%ixyz(2,i)

                myfc2_index(idx)%iatom_number = new_atom1
                myfc2_index(idx)%jatom_number = new_atom2
                myfc2_index(idx)%iatom_xyz = direction1
                myfc2_index(idx)%jatom_xyz = direction2

                atom1 = findAtom_inOld(new_atom1)
                atom2 = findAtom_inOld(new_atom2)

                !case1: if can map to old ones
                IF(atom1.ne.0 .AND. atom2.ne.0) THEN
                    myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) =&
                    & prev_fc2(atom1,atom2)%phi(direction1,direction2)
                    myfc2_index(idx)%phi_temp = prev_fc2(atom1,atom2)%phi(direction1,direction2)
                ELSE
                    myfc2_index(idx)%phi_temp = 0d0
                END IF
            END DO
        END DO
!------------------------------old method--------------------------------
!        DO idx = 1, SIZE(oldfc2_index)
!            atom1 = oldfc2_index(idx)%iatom_number
!            new_atom1 = findAtom_inNew(atom1)
!            IF(new_atom1.gt.atom_number) CYCLE
!            myfc2_index(idx)%iatom_number = new_atom1
!
!            atom2 = oldfc2_index(idx)%jatom_number
!            new_atom2 = findAtom_inNew(atom2)
!            myfc2_index(idx)%jatom_number = new_atom2
!
!            !directional index will be the same
!            direction1 = oldfc2_index(idx)%iatom_xyz
!            direction2 = oldfc2_index(idx)%jatom_xyz
!
!            myfc2_index(idx)%iatom_xyz = direction1
!            myfc2_index(idx)%jatom_xyz = direction2
!
!            ! if can map both old atoms to the new atoms
!            IF(new_atom1.ne.0 .and. new_atom2.ne.0) THEN
!                !initialize myfc2_value with old one
!                myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) = oldfc2_index(idx)%phi_temp
!                !record phi_temp for the next iteration
!                myfc2_index(idx)%phi_temp = prev_fc2(atom1,atom2)%phi(direction1,direction2)
!            ELSE
!                !cant map, then set it to 0
!                myfc2_index(idx)%phi_temp = 0d0
!            END IF
!        END DO
!--------------------------------------------------------------------------
        ![6] Check reversely
        region1 = 0
        DO i=1,SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz

            atoms = (/atom1,atom2/)
            directions = (/direction1,direction2/)
            IF(findAtom_inMap(2,atoms,directions)) THEN
                WRITE(362,*)'-----matched, continue------'
            ELSE
                region1 = region1 + 1
                WRITE(362,*)'-----found a fc2 that belongs to remapped but is not generated by the <setup_maps>-----'
            END IF
        END DO

        WRITE(362,*)'total number of indie fc2 = ', ifc2_terms
        WRITE(362,*)'========================================================================'
        WRITE(362,*)'region1: ',region1
        WRITE(362,*)'region2: ',region2
        WRITE(362,*)'region3: ',region3

        DEALLOCATE(atoms,directions)
        CLOSE(362)
    END SUBROUTINE prepare_fc2
!---------------------------------------------------
    SUBROUTINE prepare_fc3
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: old_idx, new_idx
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2, new_atom3
        INTEGER :: idx1,idx2,idx3,idx4
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp,atoms,directions


!        OPEN(363,FILE='fc3_update_info.txt',STATUS='unknown',ACTION = 'write')
!        WRITE(363,*)'=====Unique atomic index for fc3 in fc3.dat file, sorted====='
!        DO i=1,SIZE(fc3_unique_idx)
!            WRITE(363, *)'#',i,'  ',fc3_unique_idx(i)
!        END DO

        ![1] Store old unique fc3 atomic index
        IF(ALLOCATED(old_fc3_unique_idx)) DEALLOCATE(old_fc3_unique_idx)
        ALLOCATE(old_fc3_unique_idx(SIZE(fc3_unique_idx)))
        old_fc3_unique_idx = fc3_unique_idx
        ![2]&[3] Get new unique fc3 atomic index, reallocate myfc3_value
        !NOTICE: in this way new unique fc3 might be more than old ones
        !So need to reallocate myfc3_value later
        IF(ALLOCATED(fc3_unique_idx)) DEALLOCATE(fc3_unique_idx)
        CALL get_atoms_otf(3,fc3_unique_idx)

!------------------------------Old method, not good-------------------------------------------
!        ![2] Try to mapping these unique atomic index to new ones
!        !NOTICE: <atompos_Update> has to be called before, to use prev_atom(:)
!        ALLOCATE(temp(SIZE(fc3_unique_idx)))
!        WRITE(363,*)'=====Mapped unique atomic index, unsorted====='
!        DO i=1, SIZE(fc3_unique_idx)
!            old_idx = fc3_unique_idx(i)
!            new_idx = findAtom_inNew(old_idx)
!            temp(i) = new_idx ! record mapped index
!            WRITE(363,*)'old index:', old_idx, 'new index:', new_idx
!        END DO
!
!        ![3] Update fc3_unique_idx(:)
!        DEALLOCATE(fc3_unique_idx)
!        CALL unique_sort(temp, fc3_unique_idx)
!        WRITE(363,*)'=====Unique atomic index for fc3 in the new iteration, sorted====='
!        DO i=1,SIZE(fc3_unique_idx)
!            WRITE(363, *)'#',i,'  ',fc3_unique_idx(i)
!        END DO
!---------------------------------------------------------------------------------------------

        ![4] Update the group info of all fc3 for the new set of atomic index
        !...using myfc3_index(:)

        ![4.1]Store previous myfc3_value and reallocate it
        IF(ALLOCATED(prev_fc3)) DEALLOCATE(prev_fc3)
        ALLOCATE(prev_fc3,source=myfc3_value)
        IF(ALLOCATED(myfc3_value)) DEALLOCATE(myfc3_value)
        ALLOCATE(myfc3_value(atom_number, SIZE(fc3_unique_idx), SIZE(fc3_unique_idx)))

        ![4.2]Store previous myfc3_index into oldfc3_index
        IF(ALLOCATED(oldfc3_index)) DEALLOCATE(oldfc3_index)
        ALLOCATE(oldfc3_index,source=myfc3_index)

        ![4.3]Update myfc3_index and re-initialize myfc3_value, hard
        DO i=1,SIZE(oldfc3_index)
            atom1 = oldfc3_index(i)%iatom_number !first old atom index
            atom1 = findAtom_inNew(atom1) !first new atom index
            myfc3_index(i)%iatom_number = atom1 !store new one in myfc3_index(:)
            IF(atom1.gt.atom_number) CYCLE

            atom2 = oldfc3_index(i)%jatom_number !second old atom index
            atom2 = findAtom_inNew(atom2) !second new atom index
            myfc3_index(i)%jatom_number = atom2 !store new one in myfc3_index(:)

            atom3 = oldfc3_index(i)%katom_number !third old atom index
            atom3 = findAtom_inNew(atom3) !third new atom index
            myfc3_index(i)%katom_number = atom3 !store new one in myfc3_index(:)

            direction1=oldfc3_index(i)%iatom_xyz !first direction index
            direction2=oldfc3_index(i)%jatom_xyz !second direction index
            direction3=oldfc3_index(i)%katom_xyz !third direction index

            myfc3_index(i)%iatom_xyz = direction1
            myfc3_index(i)%jatom_xyz = direction2
            myfc3_index(i)%katom_xyz = direction3

            !NOTICE:fc3_unique_idx is already updated at this step
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)
                !re-initialize new fc3 value to be the same with the mapped old one
                myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3) = oldfc3_index(i)%psi_temp

            END IF

            !And don't forget to record psi_temp in this iteration for 'myfc3_index' as well
            !Here is what 'old_fc3_unique_idx' and 'prev_fc3' really for
            atom2 = oldfc3_index(i)%jatom_number !re-get second old atom index
            atom3 = oldfc3_index(i)%katom_number !re-get third old atom index

            IF(ANY(old_fc3_unique_idx==atom2).AND.ANY(old_fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(old_fc3_unique_idx,atom2)
                atom3 = find_loc(old_fc3_unique_idx,atom3)

                myfc3_index(i)%psi_temp = prev_fc3(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
            ELSE
                myfc3_index(i)%psi_temp = 0d0
            END IF

        END DO

!        CLOSE(363)

    END SUBROUTINE prepare_fc3
!---------------------------------------------------
    SUBROUTINE prepare_fc4
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: old_idx, new_idx
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER,ALLOCATABLE,DIMENSION(:) :: atoms, directions

!        OPEN(364,FILE='fc4_update_info.txt',STATUS='unknown',ACTION = 'write')
!        WRITE(364,*)'=====Unique atomic index for fc4 in fc4.dat file, sorted====='
!        DO i=1,SIZE(fc4_unique_idx)
!            WRITE(364, *)'#',i,'  ',fc4_unique_idx(i)
!        END DO

        ![1]Store old unique fc4 atomic index
        IF(ALLOCATED(old_fc4_unique_idx)) DEALLOCATE(old_fc4_unique_idx)
        ALLOCATE(old_fc4_unique_idx, source=fc4_unique_idx)

        ![2]&[3] Get new unique fc4 atomic index
        !NOTICE: in this way new unique fc4 might be more than old ones
        !So need to reallocate myfc4_value later
        IF(ALLOCATED(fc4_unique_idx)) DEALLOCATE(fc4_unique_idx)
        CALL get_atoms_otf(4,fc4_unique_idx)

!------------------------------Old method, not good-------------------------------------------
!        ![2]Try to mapping these unique atomic index to new ones
!        !NOTICE: <atompos_Update> has to be called before, to use prev_atom(:)
!        ALLOCATE(temp(SIZE(fc4_unique_idx)))
!        WRITE(364,*)'=====Mapped unique atomic index, unsorted====='
!        DO i=1, SIZE(fc4_unique_idx)
!            old_idx = fc4_unique_idx(i)
!            new_idx = findAtom_inNew(old_idx)
!            temp(i) = new_idx ! record mapped index
!            WRITE(364,*)'old index:', old_idx, 'new index:', new_idx
!        END DO
!
!        ![3]Update fc4_unique_idx(:)
!        DEALLOCATE(fc4_unique_idx)
!        CALL unique_sort(temp, fc4_unique_idx)
!        WRITE(364,*)'=====Unique atomic index for fc4 in the new iteration, sorted====='
!        DO i=1,SIZE(fc4_unique_idx)
!            WRITE(364, *)'#',i,'  ',fc4_unique_idx(i)
!        END DO
!--------------------------------------------------------------------------------------------

        ![4]Update the group info of all fc4 for the new set of atomic index
        !...using myfc4_index(:)

        ![4.1]Store previous myfc4_value and reallocate it
        IF(ALLOCATED(prev_fc4)) DEALLOCATE(prev_fc4)
        ALLOCATE(prev_fc4, source=myfc4_value)
        IF(ALLOCATED(myfc4_value)) DEALLOCATE(myfc4_value)
        ALLOCATE(myfc4_value(atom_number,SIZE(fc4_unique_idx),SIZE(fc4_unique_idx),SIZE(fc4_unique_idx)))
        ![4.2]Store previous myfc4_index into oldfc4_index
        IF(ALLOCATED(oldfc4_index)) DEALLOCATE(oldfc4_index)
        ALLOCATE(oldfc4_index, source=myfc4_index)

        ![4.3]Update myfc4_index and myfc4_value, hard
        DO i=1,SIZE(oldfc4_index)
            atom1 = oldfc4_index(i)%iatom_number !first old atom index
            atom1 = findAtom_inNew(atom1) !first new atom index
            myfc4_index(i)%iatom_number = atom1 !store new one in myfc4_index(:)
            IF(atom1.gt.atom_number) CYCLE

            atom2 = oldfc4_index(i)%jatom_number !second old atom index
            atom2 = findAtom_inNew(atom2) !second new atom index
            myfc4_index(i)%jatom_number = atom2 !store new one in myfc4_index(:)

            atom3 = oldfc4_index(i)%katom_number !third old atom index
            atom3 = findAtom_inNew(atom3) !third new atom index
            myfc4_index(i)%katom_number = atom3 !store new one in myfc4_index(:)

            atom4 = oldfc4_index(i)%latom_number !fourth old atom index
            atom4 = findAtom_inNew(atom4) !fourth new atom index
            myfc4_index(i)%latom_number = atom4 !store new one in myfc4_index(:)

            ! directional index will always be the same
            direction1=oldfc4_index(i)%iatom_xyz !first direction index
            direction2=oldfc4_index(i)%jatom_xyz !second direction index
            direction3=oldfc4_index(i)%katom_xyz !third direction index
            direction4=oldfc4_index(i)%latom_xyz !fourth direction index

            myfc4_index(i)%iatom_xyz = direction1
            myfc4_index(i)%jatom_xyz = direction2
            myfc4_index(i)%katom_xyz = direction3
            myfc4_index(i)%latom_xyz = direction4

            !NOTICE:fc4_unique_idx is already updated at this step
            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(fc4_unique_idx,atom2)
                atom3 = find_loc(fc4_unique_idx,atom3)
                atom4 = find_loc(fc4_unique_idx,atom4)
                !update new fc4 value to be the same with the mapped old one
                myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)=oldfc4_index(i)%chi_temp
            END IF

            !And don't forget to record psi_temp in this iteration for 'myfc4_index' as well
            !Here is what 'old_fc4_unique_idx' and 'prev_fc4' for
            atom2 = oldfc4_index(i)%jatom_number !re-get second old atom index
            atom3 = oldfc4_index(i)%katom_number !re-get third old atom index
            atom4 = oldfc4_index(i)%latom_number !re-get fourth old atom index

            IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
            &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(old_fc4_unique_idx,atom2)
                atom3 = find_loc(old_fc4_unique_idx,atom3)
                atom4 = find_loc(old_fc4_unique_idx,atom4)
                myfc4_index(i)%chi_temp = prev_fc4(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
            ELSE
                myfc4_index(i)%chi_temp = 0d0
            END IF

!            WRITE(364,*)i,myfc4_index(i)%iatom_number,get_letter(myfc4_index(i)%iatom_xyz),&
!            &myfc4_index(i)%jatom_number,get_letter(myfc4_index(i)%jatom_xyz),&
!            &myfc4_index(i)%katom_number,get_letter(myfc4_index(i)%katom_xyz),&
!            &myfc4_index(i)%latom_number,get_letter(myfc4_index(i)%latom_xyz),&
!            &myfc4_index(i)%chi_temp

        END DO

!        CLOSE(364)

        !----------------------------additional check----------------------------------
!        IF(.NOT.ALLOCATED(atoms)) ALLOCATE(atoms(4))
!        IF(.NOT.ALLOCATED(directions)) ALLOCATE(directions(4))
!        OPEN(374, FILE='additional_checkfc4.txt',STATUS='unknown')
!        WRITE(374,*)'=====term,(atom,xyz)*4,fc4_value====='
!        DO i=1,SIZE(myfc4_index)
!            atoms(1) = myfc4_index(i)%iatom_number
!            IF(atoms(1).gt.atom_number) CYCLE
!            atoms(2) = myfc4_index(i)%jatom_number
!            atoms(3) = myfc4_index(i)%katom_number
!            atoms(4) = myfc4_index(i)%latom_number
!            directions(1) = myfc4_index(i)%iatom_xyz
!            directions(2) = myfc4_index(i)%jatom_xyz
!            directions(3) = myfc4_index(i)%katom_xyz
!            directions(4) = myfc4_index(i)%latom_xyz
!            !the 'if' condition clause below should always be satisfied
!            IF(ANY(fc4_unique_idx==atoms(2)).AND.ANY(fc4_unique_idx==atoms(3))&
!            &.AND.ANY(fc4_unique_idx==atoms(4))) THEN
!                atom1 = find_loc(fc4_unique_idx,atoms(1))
!                atom2 = find_loc(fc4_unique_idx,atoms(2))
!                atom3 = find_loc(fc4_unique_idx,atoms(3))
!                atom4 = find_loc(fc4_unique_idx,atoms(4))
!
!!                WRITE(374,*)i,atoms(1),get_letter(direction1),atoms(2),get_letter(direction2),&
!!                &atoms(3),get_letter(direction3),atoms(4),get_letter(direction4),&
!!                &myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
!
!                WRITE(374,*)i,myfc4_index(i)%iatom_number,get_letter(myfc4_index(i)%iatom_xyz),&
!                &myfc4_index(i)%jatom_number,get_letter(myfc4_index(i)%jatom_xyz),&
!                &myfc4_index(i)%katom_number,get_letter(myfc4_index(i)%katom_xyz),&
!                &myfc4_index(i)%latom_number,get_letter(myfc4_index(i)%latom_xyz),&
!                &myfc4_index(i)%chi_temp
!            END IF
!        END DO
!        CLOSE(374)

    END SUBROUTINE prepare_fc4
!---------------------------------------------------
    SUBROUTINE uni_fc_update(rnk, atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: rnk, atoms(:), xyzs(:)
        SELECTCASE(rnk)
            CASE(2)
                CALL renew_fc2(atoms, xyzs)
            CASE(3)
                CALL renew_fc3(atoms, xyzs)
            CASE(4)
                CALL renew_fc4(atoms, xyzs)
        ENDSELECT

    END SUBROUTINE uni_fc_update
!---------------------------------------------------
    SUBROUTINE renew_fc2(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        INTEGER :: i, j, idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2

        REAL(8) :: temp
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau

!WRITE(34,*)"atom1,xyz1,atom2,xyz2",atoms(1),get_letter(xyzs(1)),atoms(2),get_letter(xyzs(2))

        ![1]get all the deviation u first
        ALLOCATE(deltau(d,SIZE(prev_atom)))
        DO atom1=1,SIZE(prev_atom)
            deltau(:,atom1)=strain(:,:).dot.prev_atom(atom1)%R+prev_atom(atom1)%tau+&
            &atomic_deviation(:,prev_atom(atom1)%type_tau)
        END DO
!WRITE(34,*)"check mark 1"
        ![2]fc4 term sum in the formula: 0.5*u*Chi*u
        DO idx = 1, SIZE(oldfc4_index)
            atom1=oldfc4_index(idx)%iatom_number
            new_atom1=findAtom_inNew(atom1)
            IF(new_atom1.ne.atoms(1)) CYCLE !Low efficient

            atom2=oldfc4_index(idx)%jatom_number
            new_atom2=findAtom_inNew(atom2)
            IF(new_atom2.ne.atoms(2)) CYCLE !Low efficient

            direction1 = oldfc4_index(idx)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficient
            direction2 = oldfc4_index(idx)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficient

            atom3=oldfc4_index(idx)%katom_number
            atom4=oldfc4_index(idx)%latom_number
            direction3 = oldfc4_index(idx)%katom_xyz
            direction4 = oldfc4_index(idx)%latom_xyz

            temp = 0d0
            !if corresponding fc4 term exists, add into
            IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
                &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                !get corresponding idx pos
                atom2 = find_loc(old_fc4_unique_idx,atom2)
                atom3 = find_loc(old_fc4_unique_idx,atom3)
                atom4 = find_loc(old_fc4_unique_idx,atom4)

                temp = 0.5*prev_fc4(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)&
                &*deltau(direction3,atom3)*deltau(direction4,atom4)

                !remember to recover the atom labels
                atom2=oldfc4_index(idx)%jatom_number
                atom3=oldfc4_index(idx)%katom_number
                atom4=oldfc4_index(idx)%latom_number
!WRITE(34,*)"atom2,atom3,atom4",atom2,atom3,atom4
!WRITE(34,*)"check term:",temp
            END IF

            myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) = &
            &myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) + temp
        END DO
!WRITE(34,*)"check mark 2"
        ![3]fc3 term sum in the formula: Psi*u
        DO idx = 1, SIZE(oldfc3_index)
            atom1 = oldfc3_index(idx)%iatom_number
            new_atom1 = findAtom_inNew(atom1)
            IF(new_atom1.ne.atoms(1)) CYCLE !Low efficiency

            atom2 = oldfc3_index(idx)%jatom_number
            new_atom2 = findAtom_inNew(atom2)
            IF(new_atom2.ne.atoms(2)) CYCLE !Low efficiency

            direction1 = oldfc3_index(idx)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficiency
            direction2 = oldfc3_index(idx)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficiency

            atom3 = oldfc3_index(idx)%katom_number
            direction3 = oldfc3_index(idx)%katom_xyz

            temp = 0d0
            !if corresponding fc3 term exists, add into
            IF(ANY(old_fc3_unique_idx==atom2).AND.ANY(old_fc3_unique_idx==atom3)) THEN
                !get corresponding idx pos
                atom2 = find_loc(old_fc3_unique_idx,atom2)
                atom3 = find_loc(old_fc3_unique_idx,atom3)

                temp = prev_fc3(atom1,atom2,atom3)%psi(direction1,direction2,direction3)*deltau(direction3,atom3)
                !remember to recover the atom labels
                atom2 = oldfc3_index(idx)%jatom_number
                atom3 = oldfc3_index(idx)%katom_number
            END IF

            myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) = &
            &myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) + temp
        END DO

        DEALLOCATE(deltau)
!WRITE(34,*)'Successfully updated this fc2!'
    END SUBROUTINE renew_fc2
!---------------------------------------------------
    SUBROUTINE renew_fc3(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        INTEGER :: i, j, idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2, new_atom3

        REAL(8) :: temp
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau

        ![1]get all the deviation u first
        ALLOCATE(deltau(d,SIZE(prev_atom)))
        DO atom1=1,SIZE(prev_atom)
            deltau(:,atom1)=strain(:,:).dot.prev_atom(atom1)%R+prev_atom(atom1)%tau+&
            &atomic_deviation(:,prev_atom(atom1)%type_tau)
        END DO
        ![2]update myfc3_value using formula
        !FORMULA: Psi_new = Psi_old + SUM(Chi_old*deltaU)
        ! Phi_new is initialized to be Phi_old in step 3.
        ! Following is Phi_new = Phi_new + SUM(Chi_old*deltaU)
        ! new fc4 info is stored in myfc4_index(:), old fc4 info is stored in oldfc4_index(:)
        ! new fc4 value is stored in myfc4_value(:), old fc4 value is stored in prev_fc4(:)
        DO i=1, SIZE(oldfc4_index)
            atom1=oldfc4_index(i)%iatom_number
            new_atom1=findAtom_inNew(atom1)
            IF(new_atom1.ne.atoms(1)) CYCLE !Low efficiency

            atom2=oldfc4_index(i)%jatom_number
            new_atom2=findAtom_inNew(atom2)
            IF(new_atom2.ne.atoms(2)) CYCLE !Low efficiency

            atom3=oldfc4_index(i)%katom_number
            new_atom3=findAtom_inNew(atom3)
            IF(new_atom3.ne.atoms(3)) CYCLE !Low efficiency

            direction1 = oldfc4_index(i)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficiency
            direction2 = oldfc4_index(i)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficiency
            direction3 = oldfc4_index(i)%katom_xyz
            IF(direction3.ne.xyzs(3)) CYCLE !Low efficiency

            atom4=oldfc4_index(i)%latom_number
            direction4 = oldfc4_index(i)%latom_xyz

            IF(ANY(fc3_unique_idx==new_atom2).AND.ANY(fc3_unique_idx==new_atom3)) THEN
                new_atom2 = find_loc(fc3_unique_idx,new_atom2)
                new_atom3 = find_loc(fc3_unique_idx,new_atom3)

                !need to make sure atomic index for fc4 are also valid
                !use 'old_fc4_unique_idx' to check
                IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
                &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                    atom2 = find_loc(old_fc4_unique_idx,atom2)
                    atom3 = find_loc(old_fc4_unique_idx,atom3)
                    atom4 = find_loc(old_fc4_unique_idx,atom4)

                    myfc3_value(new_atom1,new_atom2,new_atom3)%psi(xyzs(1),xyzs(2),xyzs(3)) = &
                    & myfc3_value(new_atom1,new_atom2,new_atom3)%psi(xyzs(1),xyzs(2),xyzs(3)) + &
                    & prev_fc4(atom1, atom2, atom3, atom4)%chi(direction1, direction2, direction3, direction4) * &
                    & deltau(direction4,atom4)
                END IF
            END IF
        END DO

        DEALLOCATE(deltau)

    END SUBROUTINE renew_fc3
!---------------------------------------------------
    SUBROUTINE renew_fc4(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        !Do nothing, fc4 is already updated in the prepare_fc4

    END SUBROUTINE renew_fc4
!================================================================================================================
! Below are subroutines related with update structure and fcs
!================================================================================================================
     SUBROUTINE trans_Update
        IMPLICIT NONE
        INTEGER :: i=0
        REAL(8) :: a, b, c
        !restore original
        cell_vec=cell_save
        trans_vec=trans_save

        !******not sure*****
        !* since original primitive matrix is not identity matrix
        !* but we'll be keeping using identity matrix for it since now
        !* we have to firstly change the trans_vec to match (a,b,c)
    !    a=latticeparameters(1)
    !    b=latticeparameters(2)
    !    c=latticeparameters(3)
    !    trans_vec(:,1)=(/a, 0d0, 0d0/)
    !    trans_vec(:,2)=(/0d0, b, 0d0/)
    !    trans_vec(:,3)=(/0d0, 0d0, c/)
        !*******************

        !update trans_vec and r1, r2, r3
        DO i=1,d
            trans_vec(:,i)=trans_vec(:,i)+(strain(:,:).dot.trans_vec(:,i))
        END DO

        r1%component=trans_vec(:,1)
        r2%component=trans_vec(:,2)
        r3%component=trans_vec(:,3)

        WRITE(*,*) "r1?",r1%component(:)
        WRITE(34,*) "r1?",r1%component(:)
        WRITE(*,*) "r2?",r2%component(:)
        WRITE(34,*) "r2?",r2%component(:)
        WRITE(*,*) "r3?",r3%component(:)
        WRITE(34,*) "r3?",r3%component(:)

        ! also update reciprocal vector and volume for
        ! other subroutines may will use them
        call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
        call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
        call calculate_volume(r1,r2,r3,volume_r)
        call calculate_volume(g1,g2,g3,volume_g)
        WRITE(*,*) "Check Product", volume_r*volume_g/(2*pi)**d !1
        WRITE(34,*) "Check Product", volume_r*volume_g/(2*pi)**d !1
    END SUBROUTINE trans_Update
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE convs_Update(conv_to_cart)! update conventional lattice vector
        IMPLICIT NONE
        INTEGER :: i
        REAL(8), INTENT(OUT) :: conv_to_cart(3,3)!translation vectors of conventional lattice
        REAL(8) :: latp(6), d2r

        d2r = 4d0*atan(1d0)/180d0 !convert degree to radian
        latp = latticeparameters
        latp(4:6) = latp(4:6)*d2r

        !firstly calculate conv_to_cart using old lattice params
        conv_to_cart = 0d0
        conv_to_cart(1,1)=latp(1)
        conv_to_cart(1,2)=latp(2)*dcos(latp(6))
        conv_to_cart(2,2)=latp(2)*dsin(latp(6))
        conv_to_cart(1,3)=latp(3)*dcos(latp(5))
        conv_to_cart(2,3)=latp(3)*(dcos(latp(4))-dcos(latp(6))*dcos(latp(5)))/dsin(latp(6))
        conv_to_cart(3,3)=sqrt(latp(3)**2-conv_to_cart(1,3)**2-conv_to_cart(2,3)**2)

        !secondly update conv_to_cart by directly multiply strain
        DO i=1, d
            conv_to_cart(:,i) = conv_to_cart(:,i) + (strain(:,:).dot.conv_to_cart(:,i))
        END DO
    END SUBROUTINE convs_Update
!-----------------------------------------------------------------------------------------------------
    SUBROUTINE latticeparameters_Update(conv_to_cart) ! update a, b, c, alpha, beta, gamma
        IMPLICIT NONE
        REAL(8) :: a, b, c, alpha, beta, garma
        REAL(8),INTENT(IN) :: conv_to_cart(3,3)
        REAL(8) :: r2d

        r2d =180d0/4d0/atan(1d0) !convert radian to degree

        ! reverse functions in <convs_Update>
        a = conv_to_cart(1,1)
        b = sqrt(conv_to_cart(1,2)**2 + conv_to_cart(2,2)**2)
        c = sqrt(conv_to_cart(3,3)**2 + conv_to_cart(1,3)**2 + conv_to_cart(2,3)**2)
        garma = atan(conv_to_cart(2,2)/conv_to_cart(1,2))
        beta = acos(conv_to_cart(1,3)/c)
        alpha = acos(conv_to_cart(2,3)/c*dsin(garma) + dcos(garma)*dcos(beta))

        alpha = alpha * r2d
        beta = beta * r2d
        garma = garma * r2d

        WRITE(*,*) "a,b,c=:",a,b,c
        WRITE(34,*) "a,b,c=:",a,b,c
        WRITE(*,*) "alpha, beta, gamma=:",alpha, beta, garma
        WRITE(34,*) "alpha, beta, gamma=:",alpha, beta, garma
        latticeparameters = (/a,b,c,alpha,beta,garma/)
    END SUBROUTINE latticeparameters_Update
!----------------------------------------------------------------------------------------------------
    SUBROUTINE atompos0_Update(conv_to_cart) ! update mipritive cell atoms position
        IMPLICIT NONE
        INTEGER :: i
        REAL(8), INTENT(IN) :: conv_to_cart(3,3)! this is conventional trans vec
        TYPE(vector) :: cR1, cR2, cR3, cG1, cG2, cG3 ! c means "conventional"
        !e.g. conv_to_cart(:,1) is R1

        ! firstly get conventional trans vectors from conv_to_cart
        cR1%component = conv_to_cart(:,1)
        cR2%component = conv_to_cart(:,2)
        cR3%component = conv_to_cart(:,3)
        WRITE(*,*) "R1?",cR1%component(:)
        WRITE(34,*) "R1?",cR1%component(:)
        WRITE(*,*) "R2?",cR2%component(:)
        WRITE(34,*) "R2?",cR2%component(:)
        WRITE(*,*) "R3?",cR3%component(:)
        WRITE(34,*) "R3?",cR3%component(:)
        ! secondly get reciprocal ones using <make_reciprocal_lattice>
        CALL make_reciprocal_lattice(cR1, cR2, cR3, cG1, cG2, cG3)
        WRITE(*,*) "G1?",cG1%component(:)
        WRITE(34,*) "G1?",cG1%component(:)
        WRITE(*,*) "G2?",cG2%component(:)
        WRITE(34,*) "G2?",cG2%component(:)
        WRITE(*,*) "G3?",cG3%component(:)
        WRITE(34,*) "G3?",cG3%component(:)
        ! thirdly update atom0(in cartesian) and atompos0(in reduced unit) accordingly
        DO i=1, natoms0
            !update my variable with the same meaning of atom0, but in cartesian
            iatom(i)%pos_tau = iatom(i)%pos_tau + atomic_deviation(:,i)

            WRITE(*,*)"old atompos0(",i,")=",atompos0(:,i)
            WRITE(34,*)"old atompos0(",i,")=",atompos0(:,i)
            !atompos0(i) uses reduced unit which are coefficients in front of conventional trans vectors
            !dot product with cartesian coords to get reduced units
            !e.g. atompos0(:,i) = atom0(i)%equilibrium_pos%component(:).dot.(/g1, g2, g3/)
            atompos0(1,i) = iatom(i)%pos_tau(:).dot.cG1%component(:)!1st coeff for ith atom
            atompos0(2,i) = iatom(i)%pos_tau(:).dot.cG2%component(:)!2nd ...
            atompos0(3,i) = iatom(i)%pos_tau(:).dot.cG3%component(:)!3rd ...

            !atom0 should also be in reduced unit
            atom0(i)%equilibrium_pos%component(:) = atompos0(:,i)
            WRITE(*,*)"new atompos0(",i,")=",atompos0(:,i)
            WRITE(34,*)"new atompos0(",i,")=",atompos0(:,i)
        END DO
    END SUBROUTINE atompos0_Update
!----------------------------------------------------------------------------------------------------
    ! primitivelattices can be left untouched
    SUBROUTINE primitivelattices_Update(conv_to_cart) ! keep the primitive lattice matrix identity
        IMPLICIT NONE
        INTEGER :: i
        REAL(8), INTENT(IN) :: conv_to_cart(3,3)! this is conventional trans vec
        TYPE(vector) :: cR1, cR2, cR3, cG1, cG2, cG3 ! c means "conventional"
        REAL(8) :: G_mat(3,3)
        !conventional trans vectors
        cR1%component = conv_to_cart(:,1)
        cR2%component = conv_to_cart(:,2)
        cR3%component = conv_to_cart(:,3)
        CALL make_reciprocal_lattice(cR1, cR2, cR3, cG1, cG2, cG3)
        G_mat(:,1)=cG1%component(:)
        G_mat(:,2)=cG2%component(:)
        G_mat(:,3)=cG3%component(:)

        primitivelattice(1,:) = trans_vec(1,:).dot.G_mat(:,:)
        primitivelattice(2,:) = trans_vec(2,:).dot.G_mat(:,:)
        primitivelattice(3,:) = trans_vec(3,:).dot.G_mat(:,:)

        WRITE(*,*) "primitivelattice(1,:)=",primitivelattice(1,:)
        WRITE(34,*) "primitivelattice(1,:)=",primitivelattice(1,:)
        WRITE(*,*) "primitivelattice(2,:)=",primitivelattice(2,:)
        WRITE(34,*) "primitivelattice(2,:)=",primitivelattice(2,:)
        WRITE(*,*) "primitivelattice(3,:)=",primitivelattice(3,:)
        WRITE(34,*) "primitivelattice(3,:)=",primitivelattice(3,:)
    END SUBROUTINE primitivelattices_Update
!------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
    !re-allocate iatomcello, iatomcell, atompos
    !get new atompos
    SUBROUTINE fcinit_Update
        IMPLICIT NONE
!        INTEGER :: counter
!        counter = 0
        maxatoms=2800; imaxat=1

        DO WHILE (imaxat.ne.0) !imaxat is a flag number, it equals 0 when some conditions are satisfied
!            counter = counter + 1
            maxatoms=maxatoms+300 !???
!            WRITE(34,*)' counter=',counter
            WRITE(6,*)' maxatoms=',maxatoms
!            WRITE(34,*)' maxatoms=',maxatoms
            IF (allocated(iatomcell0)) DEALLOCATE(iatomcell0)
            IF (allocated(iatomcell))  DEALLOCATE(iatomcell)
            IF (allocated(atompos))    DEALLOCATE(atompos)
            ALLOCATE(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
            ! run previous subroutines
            CALL force_constants_init(latticeparameters,primitivelattice,natoms0,  &
            &     atom_type,atompos0)  !(n,natom_type,atompos0)
!            WRITE(34,*)' imaxat=',imaxat
        END DO
        WRITE(34,*)'fcinit_Update finished!'
    END SUBROUTINE fcinit_Update
!----------------------------------------------------------------------------------------------------
    SUBROUTINE get_maxneighbors
        IMPLICIT NONE
        INTEGER :: i
        REAL(8) :: smallest

        !update maxneighbors
        smallest = MIN(lengtha(trans_vec(:,1)),lengtha(trans_vec(:,2)),lengtha(trans_vec(:,3)))
        maxneighbors = INT(R_0/smallest) + 2

        !temporary check: hardcoded maxneighbors=48
        maxneighbors = 84
    END SUBROUTINE
!----------------------------------------------------------------------------------------------------
    SUBROUTINE atompos_Update(conv_to_cart) ! update every atom position
        IMPLICIT NONE
        INTEGER :: i
        REAL(8), INTENT(IN) :: conv_to_cart(3,3)! this is conventional trans vec

        !Firstly, store all the previous 'every_atom(:)' for the force update mapping part
        IF(ALLOCATED(prev_atom)) DEALLOCATE(prev_atom)
        ALLOCATE(prev_atom(tot_atom_number))
        DO i=1, tot_atom_number
            prev_atom(i) = every_atom(i)
        END DO

        !Secondly, the total atom number will change after structure update
        tot_atom_number = natoms

        !Thirdly, reallocate my 'every_atom(:)' type
        IF(ALLOCATED(every_atom)) DEALLOCATE(every_atom)
        ALLOCATE(every_atom(tot_atom_number))

        !Finally, store all information in newly allocated 'every_atom(:)'
        DO i=1, tot_atom_number
            every_atom(i)%label_number = i
            every_atom(i)%x =atompos(1,i)!cartesian x coordinate
            every_atom(i)%y =atompos(2,i)!cartesian y coordinate
            every_atom(i)%z =atompos(3,i)!cartesian z coordinate
            every_atom(i)%n1 = iatomcell(1,i) !n1
            every_atom(i)%n2 = iatomcell(2,i) !n2
            every_atom(i)%n3 = iatomcell(3,i) !n3
            every_atom(i)%type_tau = iatomcell0(i) !atom type tau
            !below the trans_vec have to be updated before
            every_atom(i)%R = every_atom(i)%n1*trans_vec(:,1) + &
                &every_atom(i)%n2*trans_vec(:,2) + every_atom(i)%n3*trans_vec(:,3)
            !below the tau position has to be updated before
            every_atom(i)%tau = iatom(every_atom(i)%type_tau)%pos_tau
            !every_atom(i)%type_R will be updated in subroutine <cellvec_Update>
        END DO
        !for checking
        OPEN(222, file="new_lat_fc.dat",status='unknown',action='write')
        DO i=1, tot_atom_number
            WRITE(222,*) i,every_atom(i)%x, every_atom(i)%y, every_atom(i)%z
        END DO
        CLOSE(222)

    END SUBROUTINE atompos_Update
!---------------------------------------------------------------------------------------------------
    FUNCTION findAtom_inOld(new_idx) RESULT(atomIndex)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: new_idx
        INTEGER :: tau, n1, n2, n3
        INTEGER :: atomIndex

        tau = every_atom(new_idx)%type_tau
        n1 = every_atom(new_idx)%n1
        n2 = every_atom(new_idx)%n2
        n3 = every_atom(new_idx)%n3

        DO i=1, SIZE(prev_atom)
            IF(prev_atom(i)%type_tau .eq.tau) THEN
                IF(prev_atom(i)%n1 .eq. n1 .AND. &
                &prev_atom(i)%n2 .eq. n2 .AND. &
                &prev_atom(i)%n3 .eq. n3) THEN
                    atomIndex = prev_atom(i)%label_number
                    RETURN
                END IF
            END IF
        END DO

        WRITE(*,*) 'Could not find corresponding old atom'
        atomIndex = 0
    END FUNCTION findAtom_inOld

    FUNCTION findAtom_inNew(old_idx) RESULT(atomIndex)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(IN) :: old_idx
        INTEGER :: tau, n1, n2, n3
        INTEGER :: atomIndex

        tau = prev_atom(old_idx)%type_tau
        n1 = prev_atom(old_idx)%n1
        n2 = prev_atom(old_idx)%n2
        n3 = prev_atom(old_idx)%n3

        DO i=1, SIZE(every_atom)
            IF(every_atom(i)%type_tau .eq.tau) THEN
                IF(every_atom(i)%n1 .eq. n1 .AND. &
                &every_atom(i)%n2 .eq. n2 .AND. &
                &every_atom(i)%n3 .eq. n3) THEN
                    atomIndex = every_atom(i)%label_number
                    RETURN
                END IF
            END IF
        END DO

        WRITE(*,*) 'Could not find corresponding new atom'
        atomIndex = 0
    END FUNCTION findAtom_inNew

    FUNCTION findAtom_inMap(rnk,atoms,directions) RESULT(found)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: rnk
        INTEGER,INTENT(in) :: atoms(:),directions(:)
        INTEGER :: i,j
        LOGICAL :: found

        DO j=1,map(rnk)%ngr
            DO i=1,map(rnk)%nt(j)
                IF(ALL(map(rnk)%gr(j)%iat(:,i).eq.atoms).AND.&
                &ALL(map(rnk)%gr(j)%ixyz(:,i).eq.directions)) THEN
                    found = .TRUE.
                    RETURN
                ELSE
                    found = .FALSE.
                END IF
            END DO
        END DO

    END FUNCTION findAtom_inMap

    FUNCTION findAtom_inRegion12(rnk,atoms,directions) RESULT(found)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: rnk
        INTEGER,INTENT(in) :: atoms(:), directions(:)
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER,DIMENSION(:),ALLOCATABLE :: myatoms, mydirections

        LOGICAL :: found
        INTEGER :: i

        ALLOCATE(myatoms(rnk),mydirections(rnk))
        !for fc2
        IF(rnk.eq.2) THEN
            DO i=1,SIZE(myfc2_index)
                atom1 = myfc2_index(i)%iatom_number
                atom2 = myfc2_index(i)%jatom_number
                direction1 = myfc2_index(i)%iatom_xyz
                direction2 = myfc2_index(i)%jatom_xyz

                myatoms = (/atom1,atom2/)
                mydirections = (/direction1,direction2/)

                IF(ALL(myatoms.eq.atoms).AND.ALL(mydirections.eq.directions)) THEN
                    !I take all fc2, no need to check if its non-zero or not.
                    found = .TRUE.
                    RETURN
                ELSE
                    found = .FALSE.
                END IF
            END DO
        END IF
        !for fc3
        IF(rnk.eq.3) THEN
            DO i=1,SIZE(myfc3_index)
                atom1 = myfc3_index(i)%iatom_number
                atom2 = myfc3_index(i)%jatom_number
                atom3 = myfc3_index(i)%katom_number
                direction1 = myfc3_index(i)%iatom_xyz
                direction2 = myfc3_index(i)%jatom_xyz
                direction3 = myfc3_index(i)%katom_xyz

                myatoms = (/atom1,atom2,atom3/)
                mydirections = (/direction1,direction2,direction3/)

                IF(ALL(myatoms.eq.atoms).AND.ALL(mydirections.eq.directions)) THEN
                    !make sure that found fc3 has non-zero value,
                    !which means that corresponding myfc3_value(...) exists
                    IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                        found = .TRUE.
                        RETURN
                    ELSE
                        found = .FALSE.
                    END IF
                ELSE
                    found = .FALSE.
                END IF
            END DO
        END IF

        DEALLOCATE(myatoms, mydirections)
    END FUNCTION findAtom_inRegion12
!--------------------------------------------------------------------------------------------------
    !renew an array of atom index by remapping from old to new, using (n, tau)
    FUNCTION include_atoms_fc(rnk) RESULT(inc)
        IMPLICIT NONE
        INTEGER :: i,temp
        INTEGER, INTENT(in) :: rnk
        INTEGER,ALLOCATABLE,DIMENSION(:) :: array
        LOGICAL :: inc

        OPEN(624,FILE='inc_atoms_fc.txt',STATUS='unknown',POSITION='append', ACTION='write')
        SELECTCASE(rnk)
            CASE(2)
                ALLOCATE(array, source=atoms_fc2)
            CASE(3)
                ALLOCATE(array, source=atoms_fc3)
            CASE(4)
                ALLOCATE(array, source=atoms_fc4)
        ENDSELECT
        WRITE(624,*)'Current max shells:',maxneighbors
        WRITE(624,*) '=======RANK========',rnk
        DO i=1, SIZE(array)
            WRITE(624,*)'included atom# fc',array(i)
            temp = findAtom_inNew(array(i))
            IF(temp.eq.0) THEN
                inc = .false.
                WRITE(624,*) 'This atom is not found in the new atom list'
                RETURN
            END IF
            WRITE(624,*)'match new atom# in the new atom list',temp
            WRITE(624,*)
        END DO
        inc = .true.
        DEALLOCATE(array)
        CLOSE(624)
    END FUNCTION include_atoms_fc
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE get_nshells
        IMPLICIT NONE
        INTEGER :: rnk,i
        INTEGER,DIMENSION(:),ALLOCATABLE :: atoms_otf, array
        LOGICAL :: inc
        OPEN(625,FILE='log_nshell.txt',STATUS='unknown',ACTION='write',POSITION='append')

        !***firstly check before running <setup_maps> again, to see if mentioned atoms are
        !all included in the new atom list***
        !shouldn't comment out the following lines
!        WRITE(625,*) '=====Check pre-<setup_maps> stage======'
        DO rnk = 2, 4
            WRITE(625,*)'Current nshells value for rank',rnk,'is', nshells(rnk,1)
            DO WHILE(.not.include_atoms_fc(rnk))
                DO i=1,atom_number
                    nshells(rnk,i) = nshells(rnk,i) + 1
                END DO
                WRITE(625,*)'rank=',rnk,'nshells=',nshells(rnk,1)
                CALL fcinit_Update
                CALL atompos_Update(conv_to_cart)
            END DO
        END DO
WRITE(34,*)'Finish <get_nshells> stage1'

        CALL make_reciprocal_lattice(r01,r02,r03,g01,g02,g03)
        CALL set_neighbor_list
        CALL setup_maps
        IF(ALLOCATED(atoms_otf)) DEALLOCATE(atoms_otf)
WRITE(34,*)'Finish <get_nshells> post-stage1'

        !***secondly check after everytime running <setup_maps> (and those subroutines before it),
        !to see if mentioned atoms are all covered by the new nshells region
        WRITE(625,*) '=====Check loop-<setup_maps> stage======'
        DO rnk = 2, 4
            WRITE(625,*)'Current nshells value for rank',rnk,'is', nshells(rnk,1)

            !step 1, get new atoms list from updated map(:) object
            !which depends on the current nshells(:) value
            CALL get_atoms_otf(rnk,atoms_otf)
WRITE(34,*)'Finish <get_nshells> stage2-step1'

            !step 2, copy atoms_fc#, keep it untouched
            IF(ALLOCATED(array)) DEALLOCATE(array)
            SELECTCASE(rnk)
                CASE(2)
                    ALLOCATE(array, source=atoms_fc2)
                CASE(3)
                    ALLOCATE(array, source=atoms_fc3)
                CASE(4)
                    ALLOCATE(array, source=atoms_fc4)
            ENDSELECT
WRITE(34,*)'Finish <get_nshells> stage2-step2'


            !step 3, remap array
            DO i=1,SIZE(array)
               array(i) = findAtom_inNew(array(i))
            END DO
!!-----------------temporary check------------------
OPEN(818,FILE='checkFCatoms.dat',POSITION='append')
WRITE(818,*)'*********************************************************************'
WRITE(818,*) 'current rank:',rnk
WRITE(818,*) 'current nshells value for this rank:',nshells(rnk,1), nshells(rnk,2)
WRITE(818,*)'====== fc2 atoms that are included in initial input(remapped) ======'
DO i=1,SIZE(array)
   WRITE(818,*) '#',i,'atom number:', array(i)
END DO
WRITE(818,*)'====== fc2 atoms that are included in this nshells setting ======'
DO i=1,SIZE(atoms_otf)
    WRITE(818,*)'#',i,'atom number:', atoms_otf(i)
END DO
!-----------------------------------
            !loop check-increase
            !what to check: array is included in the atoms_otf?
            !what to increase: if not, increase nshells(rnk,:)
            DO WHILE(.not.include_arrays(array,atoms_otf))
                !all the atom in the unit cell
                DO i=1,atom_number
                    nshells(rnk,i) = nshells(rnk,i) + 1
                END DO
                !***modification
!                IF(nshells(rnk,i).ge.(maxneighbors*0.5)) THEN
!                    maxneighbors = maxneighbors + 2
!                END IF
                !***
                WRITE(625,*)'rank=',rnk,'nshells=',nshells(rnk,1)
                !pipeline up until <setup_maps>
                CALL fcinit_Update
                CALL atompos_Update(conv_to_cart)
WRITE(34,*)'*intermediate check:finish <atompos_Update>'
                CALL make_reciprocal_lattice(r01,r02,r03,g01,g02,g03)
WRITE(34,*)'*intermediate check:finish <make_reciprocal_lattice>'
                CALL set_neighbor_list
WRITE(34,*)'*intermediate check:finish <set_neighbor_list>'
                CALL setup_maps
WRITE(34,*)'*intermediate check:finish <setup_maps>'
                !update atoms_otf after increase nshells(rnk,:) and running <setup_maps>
                IF(ALLOCATED(atoms_otf)) DEALLOCATE(atoms_otf)
                CALL get_atoms_otf(rnk,atoms_otf)
WRITE(34,*)'*intermediate check:finish <get_atoms_otf>'
                !recopy array and remap after running <fcinit_Update> & <atompos_Update>
                IF(ALLOCATED(array)) DEALLOCATE(array)
                SELECTCASE(rnk)
                    CASE(2)
                        ALLOCATE(array, source=atoms_fc2)
                    CASE(3)
                        ALLOCATE(array, source=atoms_fc3)
                    CASE(4)
                        ALLOCATE(array, source=atoms_fc4)
                ENDSELECT
                DO i=1,SIZE(array)
                    array(i) = findAtom_inNew(array(i))
                END DO
!!-----------------temporary check------------------
WRITE(818,*)'*********************************************************************'
WRITE(818,*) 'current rank:',rnk
WRITE(818,*) 'current nshells value for this rank:',nshells(rnk,1), nshells(rnk,2)
WRITE(818,*)'====== fc2 atoms that are included in initial input(remapped) ======'
DO i=1,SIZE(array)
   WRITE(818,*) '#',i,'atom number:', array(i)
END DO
WRITE(818,*)'====== fc2 atoms that are included in this nshells setting ======'
DO i=1,SIZE(atoms_otf)
    WRITE(818,*)'#',i,'atom number:', atoms_otf(i)
END DO
!!-----------------------------------
            END DO
WRITE(34,*)'Finish <get_nshells> stage2-step3'
            !thirdly, update atoms_fc# by rank
            !notice here the array(:) already corresponds to the new atom list
            SELECTCASE(rnk)
                CASE(2)
                    IF(ALLOCATED(atoms_fc2)) DEALLOCATE(atoms_fc2)
                    ALLOCATE(atoms_fc2, source=array)
                CASE(3)
                    IF(ALLOCATED(atoms_fc3)) DEALLOCATE(atoms_fc3)
                    ALLOCATE(atoms_fc3, source=array)
                CASE(4)
                    IF(ALLOCATED(atoms_fc4)) DEALLOCATE(atoms_fc4)
                    ALLOCATE(atoms_fc4, source=array)
            ENDSELECT
        END DO
WRITE(34,*)'Finish <get_nshells> stage3'
        IF(ALLOCATED(array)) DEALLOCATE(array)
        CLOSE(818)
        CLOSE(625)
    END SUBROUTINE get_nshells
!---------------------------------------------------------------------------------------------------
    SUBROUTINE cellvec_Update ! very crucial
        IMPLICIT NONE
        INTEGER :: i, j, direction
        REAL(8) :: check(d)

        IF(ALLOCATED(cell_vec)) DEALLOCATE(cell_vec)
        IF(ALLOCATED(cell_save)) DEALLOCATE(cell_save)
        !the most number of cell_vec(3,:) is when every atom is in different unit cell
        !which is basically monoatomic cell
        ALLOCATE(cell_vec(d,tot_atom_number))
        ALLOCATE(cell_save(d,tot_atom_number))

        cell_vec(:,1)=every_atom(1)%R
        every_atom(1)%type_R=1 !center cell starts at R=0
        j=2
        Do i=2,tot_atom_number
            check = every_atom(i)%R-every_atom(i-1)%R ! you have to update R before
    inner:  DO direction=1,d
                IF(check(direction).ne.0) THEN
                    cell_vec(:,j)=every_atom(i)%R
                    every_atom(i)%type_R=j
                    j=j+1
                    exit inner
                ELSE
                    every_atom(i)%type_R=j-1
                END IF
            END DO inner
        END DO

        cell_save=cell_vec
        cell_number=j-2   !j=j+1, starts from j=1, so -2
    END SUBROUTINE cellvec_Update

    SUBROUTINE CheckCellvec(i)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: i
        INTEGER :: j
        CHARACTER :: name*2
        WRITE(name,'(i2)') i
        OPEN(715,FILE='CellvecCheck'//name//'.txt')

        DO j=1, tot_atom_number
            WRITE(715,*)'------------------------------------'
            WRITE(715,*)'Atom idx = ',j
            WRITE(715,*)'R = ',every_atom(j)%R
            WRITE(715,*)'type_R = ', every_atom(j)%type_R
            WRITE(715,*)'cell_vec = ', cell_vec(:,every_atom(j)%type_R)
            IF(ANY(cell_vec(:,every_atom(j)%type_R).ne.every_atom(j)%R)) THEN
                WRITE(715,*)'FOUND MISTAKE'
            END IF
        END DO
        CLOSE(715)

    END SUBROUTINE CheckCellvec
!--------------------------------------------------------------------------------------------------
    SUBROUTINE fc4_update
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: old_idx, new_idx
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp

        !Firstly store old unique fc4 atomic index
        OPEN(374,FILE='fc4_update_info.txt',STATUS='unknown',ACTION = 'write')
        WRITE(374,*)'=====Unique atomic index for fc4 in fc4.dat file, sorted====='
        DO i=1,SIZE(fc4_unique_idx)
            WRITE(374, *)'#',i,'  ',fc4_unique_idx(i)
        END DO

        IF(ALLOCATED(old_fc4_unique_idx)) DEALLOCATE(old_fc4_unique_idx)
        ALLOCATE(old_fc4_unique_idx(SIZE(fc4_unique_idx)))
        old_fc4_unique_idx = fc4_unique_idx

        !Secondly try to mapping these unique atomic index to new ones
        !NOTICE: <atompos_Update> has to be called before, to use prev_atom(:)
        ALLOCATE(temp(SIZE(fc4_unique_idx)))
        WRITE(374,*)'=====Mapped unique atomic index, unsorted====='
        DO i=1, SIZE(fc4_unique_idx)
            old_idx = fc4_unique_idx(i)
            new_idx = findAtom_inNew(old_idx)
            temp(i) = new_idx ! record mapped index
            WRITE(374,*)'old index:', old_idx, 'new index:', new_idx
        END DO

        !Thirdly update fc4_unique_idx(:)
        DEALLOCATE(fc4_unique_idx)
        CALL unique_sort(temp, fc4_unique_idx)
        WRITE(374,*)'=====Unique atomic index for fc4 in the new iteration, sorted====='
        DO i=1,SIZE(fc4_unique_idx)
            WRITE(374, *)'#',i,'  ',fc4_unique_idx(i)
        END DO

        !Fourthly update the group info of all fc4 for the new set of atomic index
        !...using myfc4_index(:)

        !1. Store previous fc4_value
        IF(ALLOCATED(prev_fc4)) DEALLOCATE(prev_fc4)
        ALLOCATE(prev_fc4(SIZE(myfc4_value, DIM=1),SIZE(myfc4_value, DIM=2), &
        &SIZE(myfc4_value, DIM=3), SIZE(myfc4_value, DIM=4)))

        prev_fc4 = myfc4_value
        !2. Store previous myfc4_index into oldfc4_index
        IF(ALLOCATED(oldfc4_index)) DEALLOCATE(oldfc4_index)
        ALLOCATE(oldfc4_index(SIZE(myfc4_index)))

        oldfc4_index = myfc4_index
        !3. Update myfc4_index and myfc4_value, hard
        DO i=1,SIZE(oldfc4_index)
            atom1 = oldfc4_index(i)%iatom_number !first old atom index
            atom1 = findAtom_inNew(atom1) !first new atom index
            myfc4_index(i)%iatom_number = atom1 !store new one in myfc4_index(:)
            IF(atom1.gt.atom_number) CYCLE

            atom2 = oldfc4_index(i)%jatom_number !second old atom index
            atom2 = findAtom_inNew(atom2) !second new atom index
            myfc4_index(i)%jatom_number = atom2 !store new one in myfc4_index(:)

            atom3 = oldfc4_index(i)%katom_number !third old atom index
            atom3 = findAtom_inNew(atom3) !third new atom index
            myfc4_index(i)%katom_number = atom3 !store new one in myfc4_index(:)

            atom4 = oldfc4_index(i)%latom_number !fourth old atom index
            atom4 = findAtom_inNew(atom4) !fourth new atom index
            myfc4_index(i)%latom_number = atom4 !store new one in myfc4_index(:)

            ! directional index will always be the same
            direction1=oldfc4_index(i)%iatom_xyz !first direction index
            direction2=oldfc4_index(i)%jatom_xyz !second direction index
            direction3=oldfc4_index(i)%katom_xyz !third direction index
            direction4=oldfc4_index(i)%latom_xyz !fourth direction index

            myfc4_index(i)%iatom_xyz = direction1
            myfc4_index(i)%jatom_xyz = direction2
            myfc4_index(i)%katom_xyz = direction3
            myfc4_index(i)%latom_xyz = direction4

            !NOTICE:fc4_unique_idx is already updated at this step
            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3).AND.ANY(fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(fc4_unique_idx,atom2)
                atom3 = find_loc(fc4_unique_idx,atom3)
                atom4 = find_loc(fc4_unique_idx,atom4)
                !update new fc4 value to be the same with the mapped old one
                myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)=oldfc4_index(i)%chi_temp
            END IF

            !And don't forget to record psi_temp in this iteration for 'myfc4_index' as well
            !Here is what 'old_fc4_unique_idx' and 'prev_fc4' for
            atom2 = oldfc4_index(i)%jatom_number !re-get second old atom index
            atom3 = oldfc4_index(i)%katom_number !re-get third old atom index
            atom4 = oldfc4_index(i)%latom_number !re-get fourth old atom index

            IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
            &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                atom2 = find_loc(old_fc4_unique_idx,atom2)
                atom3 = find_loc(old_fc4_unique_idx,atom3)
                atom4 = find_loc(old_fc4_unique_idx,atom4)
                myfc4_index(i)%chi_temp = prev_fc4(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
            ELSE
                myfc4_index(i)%chi_temp = 0d0
            END IF

        END DO

        DEALLOCATE(temp)
        CLOSE(374)
        WRITE(*,*)"=======fc4 finished update======="
    END SUBROUTINE fc4_update

    SUBROUTINE fc4_update_check
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER :: atom1, atom2, atom3,atom4
        INTEGER :: direction1, direction2, direction3,direction4
        INTEGER :: region3,region2,region1
        INTEGER, DIMENSION(:),ALLOCATABLE :: atoms, directions

        ALLOCATE(atoms(4),directions(4))

        OPEN(384, FILE='fc4_update_check.txt',STATUS='unknown')
        WRITE(384,*)'=====This subroutine is for finding difference in fc4====='
        WRITE(384,*)'=====that are generated by <setup_map> and ====='
        WRITE(384,*)'=====remapping the atom indexes====='

        region2 = 0
        region3 = 0
        DO j=1,map(4)%ngr
            DO i=1,map(4)%nt(j)
                atom1 = map(4)%gr(j)%iat(1,i)
                direction1 = map(4)%gr(j)%ixyz(1,i)
                atom2 = map(4)%gr(j)%iat(2,i)
                direction2 = map(4)%gr(j)%ixyz(2,i)
                atom3 = map(4)%gr(j)%iat(3,i)
                direction3 = map(4)%gr(j)%ixyz(3,i)
                atom4 = map(4)%gr(j)%iat(4,i)
                direction4 = map(4)%gr(j)%ixyz(4,i)
                WRITE(384,*)'-----fc4 generated by <setup_maps>-----'
                WRITE(384,*)'atom1, xyz1, atom2, xyz2, atom3, xyz3, atom4, xyz4'
                WRITE(384,*) atom1,get_letter(direction1),atom2,get_letter(direction2),&
                &atom3,get_letter(direction3), atom4, get_letter(direction4)

                atoms = (/atom1,atom2,atom3,atom4/)
                directions = (/direction1,direction2,direction3,direction4/)
                IF(findAtom_inRegion12(4,atoms,directions)) THEN
                    region2 = region2 + 1
                    WRITE(384,*) 'Found by remapping to old ones: this fc3 existed in the last iteration'
                ELSE
                    region3 = region3 + 1
                    WRITE(384,*) 'Not found by remapping to old ones'
                END IF
            END DO
        END DO

        region1 = 0
        DO i=1,SIZE(myfc4_index)
            atoms(1) = myfc4_index(i)%iatom_number
            IF(atoms(1).gt.atom_number) CYCLE
            atoms(2) = myfc4_index(i)%jatom_number
            atoms(3) = myfc4_index(i)%katom_number
            atoms(4) = myfc4_index(i)%latom_number
            directions(1) = myfc4_index(i)%iatom_xyz
            directions(2) = myfc4_index(i)%jatom_xyz
            directions(3) = myfc4_index(i)%katom_xyz
            directions(4) = myfc4_index(i)%latom_xyz
            !the 'if' condition clause below should always be satisfied
            IF(ANY(fc4_unique_idx==atoms(2)).AND.ANY(fc4_unique_idx==atoms(3))&
            &.AND.ANY(fc4_unique_idx==atoms(4))) THEN
                IF(findAtom_inMap(4,atoms,directions)) THEN
                    WRITE(384,*)'-----matched, continue------'
                ELSE
                    region1 = region1 + 1
                    WRITE(384,*)'-----found a fc4 that belongs to remapped but is not generated by the <setup_maps>-----'
                    !set it to 0?
    !                myfc4_value(atoms(1),atoms(2),atoms(3),atoms(4))%chi(&
    !                &directions(1),directions(2),directions(3),directions(4)) = 0d0
    !                myfc4_index(i)%chi_temp = 0d0
                END IF
            END IF
        END DO
        WRITE(384,*)'========================================================================'
        WRITE(384,*)'region1: ',region1
        WRITE(384,*)'region2: ',region2
        WRITE(384,*)'region3: ',region3


        CLOSE(384)
        WRITE(*,*)"=======fc4 update checked======="

!        !----------------------------additional check----------------------------------
!        OPEN(374, FILE='additional_checkfc4.txt',STATUS='unknown')
!        WRITE(374,*)'=====term,(atom,xyz)*4,fc4_value====='
!        DO i=1,SIZE(myfc4_index)
!            atoms(1) = myfc4_index(i)%iatom_number
!            IF(atoms(1).gt.atom_number) CYCLE
!            atoms(2) = myfc4_index(i)%jatom_number
!            atoms(3) = myfc4_index(i)%katom_number
!            atoms(4) = myfc4_index(i)%latom_number
!            directions(1) = myfc4_index(i)%iatom_xyz
!            directions(2) = myfc4_index(i)%jatom_xyz
!            directions(3) = myfc4_index(i)%katom_xyz
!            directions(4) = myfc4_index(i)%latom_xyz
!            !the 'if' condition clause below should always be satisfied
!            IF(ANY(fc4_unique_idx==atoms(2)).AND.ANY(fc4_unique_idx==atoms(3))&
!            &.AND.ANY(fc4_unique_idx==atoms(4))) THEN
!                atom1 = find_loc(fc4_unique_idx,atoms(1))
!                atom2 = find_loc(fc4_unique_idx,atoms(2))
!                atom3 = find_loc(fc4_unique_idx,atoms(3))
!                atom4 = find_loc(fc4_unique_idx,atoms(4))
!
!!                WRITE(374,*)i,atoms(1),get_letter(direction1),atoms(2),get_letter(direction2),&
!!                &atoms(3),get_letter(direction3),atoms(4),get_letter(direction4),&
!!                &myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)
!
!                WRITE(374,*)i,myfc4_index(i)%iatom_number,get_letter(myfc4_index(i)%iatom_xyz),&
!                &myfc4_index(i)%jatom_number,get_letter(myfc4_index(i)%jatom_xyz),&
!                &myfc4_index(i)%katom_number,get_letter(myfc4_index(i)%katom_xyz),&
!                &myfc4_index(i)%latom_number,get_letter(myfc4_index(i)%latom_xyz),&
!                &myfc4_index(i)%chi_temp
!            END IF
!        END DO
!        CLOSE(374)

        DEALLOCATE(atoms,directions)
    END SUBROUTINE fc4_update_check
!--------------------------------------------------------------------------------------------------
    SUBROUTINE fc3_update
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: old_idx, new_idx
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2, new_atom3
        INTEGER :: idx1,idx2,idx3,idx4
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp,atoms,directions
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau

        REAL(8) :: sum_fc4

        !Firstly store old unique fc3 atomic index
        OPEN(373,FILE='fc3_update_info.txt',STATUS='unknown',ACTION = 'write')
        WRITE(373,*)'=====Unique atomic index for fc3 in fc3.dat file, sorted====='
        DO i=1,SIZE(fc3_unique_idx)
            WRITE(373, *)'#',i,'  ',fc3_unique_idx(i)
        END DO

        IF(ALLOCATED(old_fc3_unique_idx)) DEALLOCATE(old_fc3_unique_idx)
        ALLOCATE(old_fc3_unique_idx(SIZE(fc3_unique_idx)))
        old_fc3_unique_idx = fc3_unique_idx

        !Secondly try to mapping these unique atomic index to new ones
        !NOTICE: <atompos_Update> has to be called before, to use prev_atom(:)
        ALLOCATE(temp(SIZE(fc3_unique_idx)))
        WRITE(373,*)'=====Mapped unique atomic index, unsorted====='
        DO i=1, SIZE(fc3_unique_idx)
            old_idx = fc3_unique_idx(i)
            new_idx = findAtom_inNew(old_idx)
            temp(i) = new_idx ! record mapped index
            WRITE(373,*)'old index:', old_idx, 'new index:', new_idx
        END DO

        !Thirdly update fc3_unique_idx(:)
        DEALLOCATE(fc3_unique_idx)
        CALL unique_sort(temp, fc3_unique_idx)
        WRITE(373,*)'=====Unique atomic index for fc3 in the new iteration, sorted====='
        DO i=1,SIZE(fc3_unique_idx)
            WRITE(373, *)'#',i,'  ',fc3_unique_idx(i)
        END DO

        !Fourthly update the group info of all fc3 for the new set of atomic index
        !...using myfc3_index(:)

        !1. Store previous fc3_value
        IF(ALLOCATED(prev_fc3)) DEALLOCATE(prev_fc3)
        ALLOCATE(prev_fc3(SIZE(myfc3_value, DIM=1),SIZE(myfc3_value, DIM=2), &
        &SIZE(myfc3_value, DIM=3)))

        prev_fc3 = myfc3_value
        !2. Store previous myfc3_index into oldfc3_index
        IF(ALLOCATED(oldfc3_index)) DEALLOCATE(oldfc3_index)
        ALLOCATE(oldfc3_index(SIZE(myfc3_index)))

        oldfc3_index = myfc3_index
        !3. update myfc3_index and re-initialize myfc3_value, hard
    !WRITE(373,*)'=====Re-Initialize====='
        DO i=1,SIZE(oldfc3_index)
            atom1 = oldfc3_index(i)%iatom_number !first old atom index
            atom1 = findAtom_inNew(atom1) !first new atom index
            myfc3_index(i)%iatom_number = atom1 !store new one in myfc3_index(:)
            IF(atom1.gt.atom_number) CYCLE

            atom2 = oldfc3_index(i)%jatom_number !second old atom index
            atom2 = findAtom_inNew(atom2) !second new atom index
            myfc3_index(i)%jatom_number = atom2 !store new one in myfc3_index(:)

            atom3 = oldfc3_index(i)%katom_number !third old atom index
            atom3 = findAtom_inNew(atom3) !third new atom index
            myfc3_index(i)%katom_number = atom3 !store new one in myfc3_index(:)

            direction1=oldfc3_index(i)%iatom_xyz !first direction index
            direction2=oldfc3_index(i)%jatom_xyz !second direction index
            direction3=oldfc3_index(i)%katom_xyz !third direction index

            myfc3_index(i)%iatom_xyz = direction1
            myfc3_index(i)%jatom_xyz = direction2
            myfc3_index(i)%katom_xyz = direction3

            !NOTICE:fc3_unique_idx is already remapped at this step
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)
    !WRITE(373,*)'--------------------------------------'
    !WRITE(373,*)'corresponding non_zero new fc3 index'
    !WRITE(373,*) atom1, atom2, atom3
    !WRITE(373,*) oldfc3_index(i)%psi_temp
                !re-initialize new fc3 value to be the same with the mapped old one
                myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3) = oldfc3_index(i)%psi_temp

            END IF


            !And don't forget to record psi_temp in this iteration for 'myfc3_index' as well
            !Here is what 'old_fc3_unique_idx' and 'prev_fc3' for
            atom2 = oldfc3_index(i)%jatom_number !re-get second old atom index
            atom3 = oldfc3_index(i)%katom_number !re-get third old atom index

            IF(ANY(old_fc3_unique_idx==atom2).AND.ANY(old_fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(old_fc3_unique_idx,atom2)
                atom3 = find_loc(old_fc3_unique_idx,atom3)
                !this value will become oldfc3_index(i)%psi_temp to initialize for the 'myfc3_value' in next iteration
                myfc3_index(i)%psi_temp = prev_fc3(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
            ELSE
                myfc3_index(i)%psi_temp = 0d0
            END IF

        END DO

        !4. Get atom deformation, for all previous atom
        ALLOCATE(deltau(d,SIZE(prev_atom)))
        DO atom1=1,SIZE(prev_atom)
            deltau(:,atom1)=strain(:,:).dot.prev_atom(atom1)%R+prev_atom(atom1)%tau+&
            &atomic_deviation(:,prev_atom(atom1)%type_tau)
        END DO

        !5. update myfc3_value using formula, hard
        !FORMULA: Phi_new = Phi_old + SUM(Chi_old*deltaU)
        ! Phi_new is initialized to be Phi_old in step 3.
        ! Following is Phi_new = Phi_new + SUM(Chi_old*deltaU)
        ! new fc4 info is stored in myfc4_index(:), old fc4 info is stored in oldfc4_index(:)
        ! new fc4 value is stored in myfc4_value(:), old fc4 value is stored in prev_fc4(:)
    !WRITE(373,*)'=====Compare====='
        DO i=1, SIZE(oldfc4_index)
            atom1=oldfc4_index(i)%iatom_number
            new_atom1=findAtom_inNew(atom1)
            IF(new_atom1.gt.atom_number) CYCLE

            atom2=oldfc4_index(i)%jatom_number
            new_atom2=findAtom_inNew(atom2)

            atom3=oldfc4_index(i)%katom_number
            new_atom3=findAtom_inNew(atom3)

            atom4=oldfc4_index(i)%latom_number

            direction1 = oldfc4_index(i)%iatom_xyz
            direction2 = oldfc4_index(i)%jatom_xyz
            direction3 = oldfc4_index(i)%katom_xyz
            direction4 = oldfc4_index(i)%latom_xyz
    !WRITE(373,*) '-----------------------------'
    !WRITE(373,*) 'old atomic index for fc4 sum '
    !WRITE(373,*)'atom1, atom2, atom3, atom4'
    !WRITE(373,*) atom1, atom2, atom3, atom4
    !WRITE(373,*) 'new atomic index for fc3 value '
    !WRITE(373,*)'new_atom1, new_atom2, new_atom3'
    !WRITE(373,*) new_atom1,new_atom2,new_atom3

            !Case1: if existent fc3
            IF(ANY(fc3_unique_idx==new_atom2).AND.ANY(fc3_unique_idx==new_atom3)) THEN
                new_atom2 = find_loc(fc3_unique_idx,new_atom2)
                new_atom3 = find_loc(fc3_unique_idx,new_atom3)
    !WRITE(373,*)'corresponding non_zero new fc3 index'
    !WRITE(373,*) new_atom1, new_atom2, new_atom3
    !WRITE(373,*)'Phi',myfc3_value(new_atom1,new_atom2,new_atom3)%psi(direction1,direction2,direction3)
    !WRITE(373,*)'Chi',prev_fc4(atom1, atom2, atom3, atom4)%chi(direction1, direction2, direction3, direction4)
    !WRITE(373,*)'deltaU',deltau(direction4,atom4)
    !WRITE(373,*)
                !need to make sure atomic index for fc4 are also valid
                !use 'old_fc4_unique_idx' to check
                IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
                &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                    atom2 = find_loc(old_fc4_unique_idx,atom2)
                    atom3 = find_loc(old_fc4_unique_idx,atom3)
                    atom4 = find_loc(old_fc4_unique_idx,atom4)

                    myfc3_value(new_atom1,new_atom2,new_atom3)%psi(direction1,direction2,direction3) = &
                    & myfc3_value(new_atom1,new_atom2,new_atom3)%psi(direction1,direction2,direction3) + &
                    & prev_fc4(atom1, atom2, atom3, atom4)%chi(direction1, direction2, direction3, direction4) * &
                    & deltau(direction4,atom4)
                END IF
            END IF

        END DO

        !*----------------------------need to check--------------------------------
        ALLOCATE(atoms(3),directions(3))
        !Case2: fc3 previously doesn't exist, but the sum of fc4 is non-zero thus a new fc3 appears
    !    DO atom1=1, atom_number
    !    DO atom2=1,SIZE(prev_atom)
    !    DO atom3=1,SIZE(prev_atom)
    !    DO direction1=1,d
    !    DO direction2=1,d
    !    DO direction3=1,d
    !        !check the sum of fc4
    !        sum_fc4 = 0d0
    !        DO atom4=1,SIZE(prev_atom)
    !        DO direction4=1,d
    !            !firstly check if this fc4 exist
    !            IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3).AND.ANY(old_fc4_unique_idx==atom4)) THEN
    !                !!remember to remap that atom2, atom3, atom4 back using find_loc!
    !                idx1 = find_loc(old_fc4_unique_idx,atom1)
    !                idx2 = find_loc(old_fc4_unique_idx,atom2)
    !                idx3 = find_loc(old_fc4_unique_idx,atom3)
    !                idx4 = find_loc(old_fc4_unique_idx,atom4)
    !
    !                sum_fc4 = sum_fc4 + prev_fc4(idx1, idx2, idx3, idx4)%chi(direction1, direction2, direction3, direction4) * &
    !                & deltau(direction4,atom4)
    !            END IF
    !
    !        END DO !xyz4
    !        END DO !atom4
    !        !secondly check if the sum is non-zero, new corresponding fc3_unique_idx has to be recorded
    !        IF(ABS(sum_fc4).gt.1e-7) THEN
    !            new_atom1 = findAtom_inNew(atom1)
    !            new_atom2 = findAtom_inNew(atom2)
    !            new_atom3 = findAtom_inNew(atom3)
    !            CALL fc3_update_unique(new_atom2) !automatically check and add it to fc3_unique_idx
    !            CALL fc3_update_unique(new_atom3) !automatically check and add it to fc3_unique_idx
    !            atoms = (/new_atom1,new_atom2,new_atom3/)
    !            directions = (/direction1,direction2,direction3/)
    !            CALL fc3_add(atoms,directions,sum_fc4)
    !        END IF
    !
    !
    !    END DO !xyz3
    !    END DO !xyz2
    !    END DO !xyz1
    !    END DO !atom3
    !    END DO !atom2
    !    END DO !atom1

        DEALLOCATE(atoms,directions)
        DEALLOCATE(deltau)
        DEALLOCATE(temp)
        CLOSE(373)

        WRITE(*,*)"=======fc3 finished update======="
    END SUBROUTINE fc3_update

    SUBROUTINE fc3_update_check
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER :: atom1, atom2, atom3, direction1, direction2, direction3
        INTEGER :: region3,region2,region1
        INTEGER, DIMENSION(:),ALLOCATABLE :: atoms, directions

        INTEGER :: check

        ALLOCATE(atoms(3),directions(3))

        OPEN(383, FILE='fc3_update_check.txt',STATUS='unknown')
        WRITE(383,*)'=====This subroutine is for finding difference in fc3====='
        WRITE(383,*)'=====that are generated by <setup_map> and ====='
        WRITE(383,*)'=====remapping the atom indexes====='
        check = 0
        region2 = 0
        region3 = 0
        DO j=1,map(3)%ngr
            DO i=1,map(3)%nt(j)
                atom1 = map(3)%gr(j)%iat(1,i)
                direction1 = map(3)%gr(j)%ixyz(1,i)
                atom2 = map(3)%gr(j)%iat(2,i)
                direction2 = map(3)%gr(j)%ixyz(2,i)
                atom3 = map(3)%gr(j)%iat(3,i)
                direction3 = map(3)%gr(j)%ixyz(3,i)
                WRITE(383,*)'-----fc3 generated by <setup_maps>-----'
                WRITE(383,*)'atom1, xyz1, atom2, xyz2, atom3, xyz3'
                WRITE(383,*) atom1,get_letter(direction1),atom2,get_letter(direction2),&
                &atom3,get_letter(direction3)

                atoms = (/atom1,atom2,atom3/)
                directions = (/direction1,direction2,direction3/)
                IF(findAtom_inRegion12(3,atoms,directions)) THEN
                    region2 = region2 + 1
                    WRITE(383,*) 'Found by remapping to old ones: this fc3 existed in the last iteration'
                ELSE
                    region3 = region3 + 1
                    WRITE(383,*) 'Not found by remapping to old ones'
                END IF

                IF(findAtom_inMap(3,atoms,directions)) THEN
                    check = check + 1
                END IF

            END DO
        END DO

        region1 = 0

        DO i=1,SIZE(myfc3_index)
            atoms(1) = myfc3_index(i)%iatom_number
            IF(atoms(1).gt.atom_number) CYCLE
            atoms(2) = myfc3_index(i)%jatom_number
            atoms(3) = myfc3_index(i)%katom_number
            directions(1) = myfc3_index(i)%iatom_xyz
            directions(2) = myfc3_index(i)%jatom_xyz
            directions(3) = myfc3_index(i)%katom_xyz

            !the 'if' condition clause below should always be satisfied
            IF(ANY(fc3_unique_idx==atoms(2)).AND.ANY(fc3_unique_idx==atoms(3))) THEN
                IF(findAtom_inMap(3,atoms,directions)) THEN
                    WRITE(383,*)'-----matched, continue------'
                ELSE
                    region1 = region1 + 1
                    WRITE(383,*)'-----found a fc3 that belongs to remapped but is not generated by the <setup_maps>-----'
                    !set it to 0?
    !                myfc3_value(atoms(1),atoms(2),atoms(3))%psi(&
    !                &directions(1),directions(2),directions(3)) = 0d0
    !                myfc3_index(i)%psi_temp = 0d0
                END IF
            END IF
        END DO
        WRITE(383,*)'========================================================================'
        WRITE(383,*)'region1: ',region1
        WRITE(383,*)'region2: ',region2
        WRITE(383,*)'region3: ',region3
        WRITE(383,*)'check: ',check

        CLOSE(383)

        fc3_1 = region1
        fc3_2 = region2
        fc3_3 = region3

        WRITE(*,*)"=======fc3 update checked======="
!        !----------------------------additional check----------------------------------
!        OPEN(373, FILE='additional_checkfc3.txt',STATUS='unknown')
!        WRITE(373,*)'=====term,(atom,xyz)*3,fc3_value====='
!        DO i=1,SIZE(myfc3_index)
!            atoms(1) = myfc3_index(i)%iatom_number
!            IF(atoms(1).gt.atom_number) CYCLE
!            atoms(2) = myfc3_index(i)%jatom_number
!            atoms(3) = myfc3_index(i)%katom_number
!            directions(1) = myfc3_index(i)%iatom_xyz
!            directions(2) = myfc3_index(i)%jatom_xyz
!            directions(3) = myfc3_index(i)%katom_xyz
!
!            !the 'if' condition clause below should always be satisfied
!            IF(ANY(fc3_unique_idx==atoms(2)).AND.ANY(fc3_unique_idx==atoms(3))) THEN
!                atom1=find_loc(fc3_unique_idx,atoms(1))
!                atom2=find_loc(fc3_unique_idx,atoms(2))
!                atom3=find_loc(fc3_unique_idx,atoms(3))
!
!                WRITE(373,*)i,atoms(1),get_letter(direction1),atoms(2),get_letter(direction2),&
!                &atoms(3),get_letter(direction3),myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
!            END IF
!        END DO
!        CLOSE(373)

        DEALLOCATE(atoms,directions)
    END SUBROUTINE fc3_update_check
!------------------------------------------------------------------------------------------------
    !utility subroutines
    !not really add in terms of set operation. It append a new element to the end of arrayIn,
    !output as arrayOut
    SUBROUTINE fortran_add(arrayIn,element,arrayOut)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: arrayIn(:),element
        INTEGER,INTENT(out),ALLOCATABLE,DIMENSION(:) :: arrayOut
        IF(.NOT.ANY(arrayIn==element)) THEN
            ALLOCATE(arrayOut(SIZE(arrayIn)+1),source=(/arrayIn,element/))
        ELSE
            ALLOCATE(arrayOut(SIZE(arrayIn)),source=arrayIn)
        END IF

    END SUBROUTINE fortran_add
!----------------------------------------------------------------------------------------------
    !utility subroutines for fc3
    SUBROUTINE fc3_update_unique(new_atom)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: new_atom
        INTEGER,ALLOCATABLE,DIMENSION(:) :: array

        CALL fortran_add(fc3_unique_idx,new_atom,array)

        DEALLOCATE(fc3_unique_idx)
        ALLOCATE(fc3_unique_idx(SIZE(array)),source=array)
        DEALLOCATE(array)
    END SUBROUTINE fc3_update_unique

    !when use this, atoms have to be the member of fc3_unique_idx
    SUBROUTINE fc3_add(atoms,directions,fc3_val)
        IMPLICIT NONE
        INTEGER,INTENT(in),DIMENSION(:)::atoms,directions
        REAL(8),INTENT(in)::fc3_val
        INTEGER :: atom1, atom2, atom3
        INTEGER :: i, j, k
        INTEGER :: whatever=0

        TYPE(fc3_value),DIMENSION(:,:,:),ALLOCATABLE :: array
        TYPE(fc3_index),DIMENSION(:),ALLOCATABLE :: array2

        atom1 = find_loc(fc3_unique_idx,atoms(1))
        atom2 = find_loc(fc3_unique_idx,atoms(2))
        atom3 = find_loc(fc3_unique_idx,atoms(3))

        !add new term to array then copy it to myfc3_value
        ALLOCATE(array(atom_number, SIZE(fc3_unique_idx), SIZE(fc3_unique_idx)))
        DO i=1,SIZE(myfc3_value,DIM=1)
        DO j=1,SIZE(myfc3_value,DIM=2)
        DO k=1,SIZE(myfc3_value,DIM=3)
            array(i,j,k)%psi = 0d0 !initialize to 0 first
            array(i,j,k) = myfc3_value(i,j,k) !then copy myfc3_value to array
        END DO
        END DO
        END DO
        array(atom1,atom2,atom3)%psi = 0d0
        array(atom1,atom2,atom3)%psi(directions(1),directions(2),directions(3)) = fc3_val
        DEALLOCATE(myfc3_value)
        ALLOCATE(myfc3_value,source=array)
        DEALLOCATE(array)
        !problem imbedded below
        !add new term to array2 then copy it to myfc3_index
        ALLOCATE(array2(SIZE(myfc3_index)+1))
        DO i=1,SIZE(myfc3_index)
            array2(i)=myfc3_index(i)
        END DO
        !after the loop, i=SIZE(myfc3_index)+1
        array2(i)%group = whatever !this info is never used
        array2(i)%iatom_number = atoms(1)
        array2(i)%jatom_number = atoms(2)
        array2(i)%katom_number = atoms(3)
        array2(i)%iatom_xyz = directions(1)
        array2(i)%jatom_xyz = directions(2)
        array2(i)%katom_xyz = directions(3)
        array2(i)%junk = whatever !this infor is never used
        array2(i)%psi_temp = fc3_val
        DEALLOCATE(myfc3_index)
        ALLOCATE(myfc3_index,source=array2)
        DEALLOCATE(array2)

    END SUBROUTINE fc3_add
!-------------------------------------------------------------------------------------------------

    SUBROUTINE fc2_update
        IMPLICIT NONE
        INTEGER :: i, j, idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2
        INTEGER :: region1, region2, region3
        INTEGER :: tot_indie_fc2, ntindp

        REAL(8) :: temp
        REAL(8),DIMENSION(:),ALLOCATABLE :: arrayIn, arrayOut
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau

        INTEGER,DIMENSION(:), ALLOCATABLE :: atoms, directions


        !***Firstly store old fc2_group and old indiefc2_index info
        OPEN(372,FILE='fc2_update_check.txt',STATUS='unknown',ACTION = 'write')
        WRITE(372,*)'===== # RANK 2 tensors :term,group,(iatom,ixyz)_2 ====='

        IF(ALLOCATED(old_indiefc2_index)) DEALLOCATE(old_indiefc2_index)
        ALLOCATE(old_indiefc2_index,source=indiefc2_index)

        IF(ALLOCATED(old_myfc2_group)) DEALLOCATE(old_myfc2_group)
        ALLOCATE(old_myfc2_group, source=myfc2_group)


        !***Secondly, remap 'indiefc2_index' info from <setup_maps>
        !number of fc2 groups: map(2)%ngr --> j
        !number of independent fc2 in each group: map(2)%ntind(j)
        tot_indie_fc2 = 0
        DO i = 1, map(2)%ngr
            tot_indie_fc2 = tot_indie_fc2 + map(2)%ntind(i)
        END DO
        IF(ALLOCATED(indiefc2_index)) DEALLOCATE(indiefc2_index)
        ALLOCATE(indiefc2_index(tot_indie_fc2))
        !'tot_indie_fc2' is the 'indie_number' in <read_map>
        ifc2_terms = tot_indie_fc2
        variational_parameters_size(3) = ifc2_terms
        IF(ALLOCATED(myfc2_group)) DEALLOCATE(myfc2_group)
        ALLOCATE(myfc2_group(d,tot_atom_number,d,tot_atom_number))
        WRITE(372,*)'================indiefc2 got from <setup_maps>========================'

        idx = 0
        DO j = 1, map(2)%ngr
            DO i = 1, map(2)%ntind(j)
                idx = idx + 1
                indiefc2_index(idx)%group = idx

                atom1 = map(2)%gr(j)%iatind(1,i)
                atom2 = map(2)%gr(j)%iatind(2,i)
                direction1 = map(2)%gr(j)%ixyzind(1,i)
                direction2 = map(2)%gr(j)%ixyzind(2,i)

                indiefc2_index(idx)%iatom_number = atom1
                indiefc2_index(idx)%iatom_xyz = direction1
                indiefc2_index(idx)%jatom_number = atom2
                indiefc2_index(idx)%jatom_xyz = direction2

                myfc2_group(direction1,atom1,direction2,atom2)%group = (/j,0/)
                myfc2_group(direction1,atom1,direction2,atom2)%mat = (/1d0,0d0/)
                !------------------------------------------------------------
                WRITE(372,*)idx
                WRITE(372,*)'atom1,direction1,atom2,direction2'
                WRITE(372,*)atom1,get_letter(direction1),atom2,get_letter(direction2)
                !------------------------------------------------------------

            END DO
        END DO
        !***Thirdly, remap 'myfc2_group'
        region2 = 0
        region3 = 0
        ALLOCATE(atoms(2),directions(2))
        DO j = 1, map(2)%ngr
            ntindp = map(2)%ntind(j)
            DO i = 1, map(2)%nt(j)
                atom1 = map(2)%gr(j)%iat(1,i)
                direction1 = map(2)%gr(j)%ixyz(1,i)
                atom2 = map(2)%gr(j)%iat(2,i)
                direction2 = map(2)%gr(j)%ixyz(2,i)
                IF(ntindp.eq.1) THEN
                    myfc2_group(direction1,atom1,direction2,atom2)%group = (/j,0/)
                    myfc2_group(direction1,atom1,direction2,atom2)%mat = (/map(2)%gr(j)%mat(i,1),0d0/)
                ELSE
                    myfc2_group(direction1,atom1,direction2,atom2)%group = (/j-1,j/)
                    myfc2_group(direction1,atom1,direction2,atom2)%mat&
                    & = (/map(2)%gr(j)%mat(i,1),map(2)%gr(j)%mat(i,2)/)
                END IF

                WRITE(372,*)'-----fc2 generated by <setup_maps>-----'
                WRITE(372,*)'atom1, xyz1, atom2, xyz2'
                WRITE(372,*)atom1,get_letter(direction1),atom2,get_letter(direction2)
                atoms = (/atom1,atom2/)
                directions = (/direction1,direction2/)
                IF(findAtom_inRegion12(2,atoms,directions)) THEN
                    region2 = region2 + 1
                    WRITE(372,*) 'Found by remapping to old ones: this fc2 existed in the last iteration'
                ELSE
                    region3 = region3 + 1
                    WRITE(372,*) 'Not found by remapping to old ones'
                END IF

            END DO
        END DO

    !---------------------------------------------------------------------------------------
    !-- Or we can forget about <setup_maps>, only run it once at the beginning,           --
    !-- get all the indie fc2 then track them all the ways through the whole calculation. --
    !---------------------------------------------------------------------------------------
    !DO i = 1, SIZE(old_indiefc2_index)
    !    atom1 = old_indiefc2_index(i)%iatom_number
    !    indiefc2_index(i)%iatom_number = atom1
    !    atom2 = old_indiefc2_index(i)%jatom_number
    !    atom2 = findAtom_inNew(atom2)
    !    IF(.NOT.atom2.eq.0) THEN
    !        indiefc2_index(i)%jatom_number = atom2
    !    ELSE
    !        WRITE(*,*) "atom2 cannot be found in the new atom list"
    !        STOP
    !    END IF
    !END DO
    ! No need to do anything with myfc2_group, there is also no need to reallocate
    !---------------------------------------------------------------------------------------
        !***Fourthly, store previous 'myfc2_value' and 'myfc2_index'
        IF(ALLOCATED(oldfc2_index)) DEALLOCATE(oldfc2_index)
        ALLOCATE(oldfc2_index(SIZE(myfc2_index)))
        IF(ALLOCATED(prev_fc2)) DEALLOCATE(prev_fc2)
        ALLOCATE(prev_fc2(SIZE(myfc2_value, DIM=1), SIZE(myfc2_value, DIM=2)))
        oldfc2_index = myfc2_index
        prev_fc2 = myfc2_value
        !***Fifthly, update 'myfc2_index' by remapping
        DO idx = 1, SIZE(oldfc2_index)
            atom1 = oldfc2_index(idx)%iatom_number
            new_atom1 = findAtom_inNew(atom1)
            IF(new_atom1.gt.atom_number) CYCLE
            myfc2_index(idx)%iatom_number = new_atom1

            atom2 = oldfc2_index(idx)%jatom_number
            new_atom2 = findAtom_inNew(atom2)
            myfc2_index(idx)%jatom_number = new_atom2

            !directional index will be the same
            direction1 = oldfc2_index(idx)%iatom_xyz
            direction2 = oldfc2_index(idx)%jatom_xyz

            myfc2_index(idx)%iatom_xyz = direction1
            myfc2_index(idx)%jatom_xyz = direction2

            ! if can map both old atoms to the new atoms
            IF(new_atom1.ne.0 .and. new_atom2.ne.0) THEN
                !initialize myfc2_value with old one
                myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) = oldfc2_index(idx)%phi_temp
                !record phi_temp for the next iteration
                myfc2_index(idx)%phi_temp = prev_fc2(atom1,atom2)%phi(direction1,direction2)
            ELSE
                !cant map
                myfc2_index(idx)%phi_temp = 0d0
            END IF
        END DO
        !***Sixthly, update 'myfc2_value' using all the info above + the formula of 1st order correction

        !Get atom deformation, for all previous atom
        ALLOCATE(deltau(d,SIZE(prev_atom)))
        DO atom1=1,SIZE(prev_atom)
            deltau(:,atom1)=strain(:,:).dot.prev_atom(atom1)%R+prev_atom(atom1)%tau+&
            &atomic_deviation(:,prev_atom(atom1)%type_tau)
        END DO
        !fc4 term sum in the formula
        DO idx = 1, SIZE(oldfc4_index)
            atom1=oldfc4_index(i)%iatom_number
            new_atom1=findAtom_inNew(atom1)
            IF(new_atom1.gt.atom_number) CYCLE

            atom2=oldfc4_index(i)%jatom_number
            new_atom2=findAtom_inNew(atom2)

            atom3=oldfc4_index(i)%katom_number

            atom4=oldfc4_index(i)%latom_number

            direction1 = oldfc4_index(i)%iatom_xyz
            direction2 = oldfc4_index(i)%jatom_xyz
            direction3 = oldfc4_index(i)%katom_xyz
            direction4 = oldfc4_index(i)%latom_xyz

            temp = 0d0
            !if corresponding fc4 term exists, add into
            IF(ANY(old_fc4_unique_idx==atom2).AND.ANY(old_fc4_unique_idx==atom3)&
                &.AND.ANY(old_fc4_unique_idx==atom4)) THEN
                !get corresponding idx pos
                atom2 = find_loc(old_fc4_unique_idx,atom2)
                atom3 = find_loc(old_fc4_unique_idx,atom3)
                atom4 = find_loc(old_fc4_unique_idx,atom4)

                temp = 0.5*prev_fc4(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)&
                &*deltau(direction3,atom3)*deltau(direction4,atom4)

                !remember to recover the atom labels
                atom2=oldfc4_index(i)%jatom_number
                atom3=oldfc4_index(i)%katom_number
                atom4=oldfc4_index(i)%latom_number
            END IF

            myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) = &
            &myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) + temp
        END DO
        !fc3 term sum in the formula
        DO idx = 1, SIZE(oldfc3_index)
            atom1 = oldfc3_index(idx)%iatom_number
            new_atom1 = findAtom_inNew(atom1)
            IF(new_atom1.gt.atom_number) CYCLE

            atom2 = oldfc3_index(idx)%jatom_number
            new_atom2 = findAtom_inNew(atom2)

            atom3 = oldfc3_index(idx)%katom_number

            direction1 = oldfc3_index(idx)%iatom_xyz
            direction2 = oldfc3_index(idx)%jatom_xyz
            direction3 = oldfc3_index(idx)%katom_xyz

            temp = 0d0
            !if corresponding fc3 term exists, add into
            IF(ANY(old_fc3_unique_idx==atom2).AND.ANY(old_fc3_unique_idx==atom3)) THEN
                !get corresponding idx pos
                atom2 = find_loc(old_fc3_unique_idx,atom2)
                atom3 = find_loc(old_fc3_unique_idx,atom3)

                temp = prev_fc3(atom1,atom2,atom3)%psi(direction1,direction2,direction3)*deltau(direction3,atom3)
                !remember to recover the atom labels
                atom2 = oldfc3_index(idx)%jatom_number
                atom3 = oldfc3_index(idx)%katom_number
            END IF

            myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) = &
            &myfc2_value(new_atom1,new_atom2)%phi(direction1,direction2) + temp

        END DO

        ALLOCATE(arrayIn(SIZE(myfc2_index)))
        region1 = 0
        !check reversely
        DO i=1,SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz

            atoms = (/atom1,atom2/)
            directions = (/direction1,direction2/)
            IF(findAtom_inMap(2,atoms,directions)) THEN
                WRITE(372,*)'-----matched, continue------'
            ELSE
                region1 = region1 + 1
                WRITE(372,*)'-----found a fc2 that belongs to remapped but is not generated by the <setup_maps>-----'
                !set it to 0?
    !            myfc2_value(atom1,atom2)%phi(direction1,direction2) = 0d0
    !            myfc2_index(i)%phi_temp = 0d0
            END IF
            myfc2_index(i)%phi_temp = myfc2_value(atom1,atom2)%phi(direction1,direction2)
            arrayIn(i) = myfc2_index(i)%phi_temp
        END DO

        WRITE(372,*)'===========indie fc2 got after 1st order correction=========='

        !clear and reallocate 'myfc2_group'
        IF(ALLOCATED(myfc2_group)) DEALLOCATE(myfc2_group) !the 'tot_atom_number' has changed in <atompos_update>
        ALLOCATE(myfc2_group(d,tot_atom_number,d,tot_atom_number))
        !initialize 'myfc2_group'
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

        !***Seventhly, update 'indiefc2_index' again
        CALL unique_sort2(arrayIn,arrayOut)
        tot_indie_fc2 = SIZE(arrayOut)
        IF(ALLOCATED(indiefc2_index)) DEALLOCATE(indiefc2_index)
        ALLOCATE(indiefc2_index(tot_indie_fc2))
        !'tot_indie_fc2' is the 'indie_number' in <read_map>
        ifc2_terms = tot_indie_fc2
        variational_parameters_size(3) = ifc2_terms
        idx = 0
        DO i=1,SIZE(arrayOut)
            idx = find_loc2(arrayIn,arrayOut(i))

            atom1 = myfc2_index(idx)%iatom_number
            atom2 = myfc2_index(idx)%jatom_number
            direction1 = myfc2_index(idx)%iatom_xyz
            direction2 = myfc2_index(idx)%jatom_xyz

            indiefc2_index(i)%group = i
            indiefc2_index(i)%iatom_number = atom1
            indiefc2_index(i)%iatom_xyz = direction1
            indiefc2_index(i)%jatom_number = atom2
            indiefc2_index(i)%jatom_xyz = direction2
            indiefc2_index(i)%phi_temp = arrayOut(i)

            myfc2_group(direction1,atom1,direction2,atom2)%group = (/i,0/)
            myfc2_group(direction1,atom1,direction2,atom2)%mat = (/1d0,0d0/)

            WRITE(372,*)
            WRITE(372,*)'corresponding myfc2_index=',idx
            WRITE(372,*)'atom1,direction1,atom2,direction2'
            WRITE(372,*)atom1,get_letter(direction1),atom2,get_letter(direction2)
            WRITE(372,*)'corresponding myfc2_value=',arrayOut(i)

        END DO

        !***Eighthly, need to regroup the rest of fc2, the reducible ones, myfc2_group
        DO i=1,SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz
            idx = find_loc2(arrayOut,myfc2_index(i)%phi_temp)
            !no way to recover the mat in this method
            myfc2_group(direction1,atom1,direction2,atom2)%group = (/idx,0/)
            myfc2_group(direction1,atom1,direction2,atom2)%mat = (/1d0,0d0/)
            !---quick check---
            If(idx.eq.0) THEN
                WRITE(372,*)"something is wrong",idx
                WRITE(372,*)"fc2 index =",i
                WRITE(372,*)"atom1,atom2,direction1,direction2"
                WRITE(372,*) atom1,atom2,get_letter(direction1),get_letter(direction2)
                WRITE(372,*) "fc2 value =",myfc2_index(i)%phi_temp
            End If
    !        IF(i.eq.610) THEN
    !            WRITE(372,*)"quick check"
    !            WRITE(372,*)"group#", idx
    !            WRITE(372,*)"fc2 value found by indiefc2",indiefc2_index(idx)%phi_temp
    !            WRITE(372,*)"fc2 value found by myfc2_value",myfc2_value(atom1,atom2)%phi(direction1,direction2)
    !            WRITE(372,*)"fc2 value found by myfc2_index",myfc2_index(i)%phi_temp
    !        END IF
        END DO

        WRITE(372,*)'total number of indie fc2 = ', ifc2_terms
        WRITE(372,*)'========================================================================'
        WRITE(372,*)'region1: ',region1
        WRITE(372,*)'region2: ',region2
        WRITE(372,*)'region3: ',region3
        DEALLOCATE(arrayIn,arrayOut)
        DEALLOCATE(atoms, directions)
        DEALLOCATE(deltau)
        CLOSE(372)

        fc2_1 = region1
        fc2_2 = region2
        fc2_3 = region3

        WRITE(*,*)"=======fc2 finished update======="
    END SUBROUTINE fc2_update
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------

!========================================================================================================
!--------------------------------------------ASR related-------------------------------------------------
!========================================================================================================

    SUBROUTINE asr_fc2
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: atom1, atom2, direction1, direction2
        TYPE(fc2_value),ALLOCATABLE,DIMENSION(:) :: fc2_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR check for fc2================'

        ALLOCATE(fc2_sum(atom_number))

        DO i=1, atom_number
            fc2_sum(i)%phi = 0d0
        END DO

        DO i=1, SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz
            fc2_sum(atom1)%phi(direction1,direction2) = fc2_sum(atom1)%phi(direction1,direction2) +&
            &myfc2_value(atom1,atom2)%phi(direction1,direction2)
        END DO

        DO i=1, atom_number
        DO direction1=1,d
        DO direction2=1,d
            IF(ABS(fc2_sum(i)%phi(direction1,direction2)).gt.resolution) THEN
                WRITE(34,*)"asr for atom1: ",i,"    xyz1: ",get_letter(direction1),&
                &"    xyz2: ",get_letter(direction2)
                WRITE(34,*)fc2_sum(i)%phi(direction1,direction2)
                WRITE(34,*)'---------------------------------------'
            END IF
        END DO
        END DO
        END DO

        DEALLOCATE(fc2_sum)
!        CLOSE(34)
    END SUBROUTINE asr_fc2

    SUBROUTINE asr_fc3
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER :: atom1, atom2,atom3, direction1, direction2, direction3
        INTEGER :: atom1_idx, atom2_idx, atom3_idx
        TYPE(fc3_value),ALLOCATABLE,DIMENSION(:,:) :: fc3_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR check for fc3================'

        ALLOCATE(fc3_sum(atom_number,tot_atom_number))

        DO i=1, atom_number
            DO j=1,tot_atom_number
                fc3_sum(i,j)%psi = 0d0
            END DO
        END DO

        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc3_unique_idx,atom1)

            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number
            direction1 = myfc3_index(i)%iatom_xyz
            direction2 = myfc3_index(i)%jatom_xyz
            direction3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2) .AND. ANY(fc3_unique_idx==atom3)) THEN
                atom2_idx = find_loc(fc3_unique_idx,atom2)
                atom3_idx = find_loc(fc3_unique_idx,atom3)
                fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) =&
                & fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) +&
                & myfc3_value(atom1_idx,atom2_idx,atom3_idx)%psi(direction1,direction2,direction3)
            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
            IF(ABS(fc3_sum(i,j)%psi(direction1,direction2,direction3)).gt.resolution) THEN
                WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,&
                &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
                &"    xyz3: ",get_letter(direction3)
                WRITE(34,*)fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'---------------------------------------'
            END IF
        END DO
        END DO
        END DO
        END DO
        END DO

        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc3_unique_idx,atom1)

            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number
            direction1 = myfc3_index(i)%iatom_xyz
            direction2 = myfc3_index(i)%jatom_xyz
            direction3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2) .AND. ANY(fc3_unique_idx==atom3)) THEN
                atom2_idx = find_loc(fc3_unique_idx,atom2)
                atom3_idx = find_loc(fc3_unique_idx,atom3)
                fc3_sum(atom1,atom3)%psi(direction1,direction2,direction3) =&
                & fc3_sum(atom1,atom3)%psi(direction1,direction2,direction3) +&
                & myfc3_value(atom1_idx,atom2_idx,atom3_idx)%psi(direction1,direction2,direction3)
            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
            IF(ABS(fc3_sum(i,j)%psi(direction1,direction2,direction3)).gt.resolution) THEN
                WRITE(34,*)"asr for atom1: ",i,"    atom3: ",j,&
                &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
                &"    xyz3: ",get_letter(direction3)
                WRITE(34,*)fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'---------------------------------------'
            END IF
        END DO
        END DO
        END DO
        END DO
        END DO
        DEALLOCATE(fc3_sum)
!        CLOSE(34)
    END SUBROUTINE asr_fc3

    SUBROUTINE asr_fc4
        IMPLICIT NONE
        INTEGER :: i,j,k
        INTEGER :: atom1, atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: atom1_idx, atom2_idx, atom3_idx, atom4_idx
        TYPE(fc4_value),ALLOCATABLE,DIMENSION(:,:,:) :: fc4_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR check for fc4================'

        ALLOCATE(fc4_sum(atom_number,tot_atom_number,tot_atom_number))

        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO
        !sum the 4th atom index
        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO k=1, tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom3: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)'---------------------------------------'
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !sum the 3rd atom index
        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO k=1, tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom4: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)'---------------------------------------'
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !sum the 2nd atom index
        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO k=1, tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom3: ",j,"    atom4: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)'---------------------------------------'
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        DEALLOCATE(fc4_sum)
!        CLOSE(34)
    END SUBROUTINE asr_fc4
!========================================================================================================
    SUBROUTINE fix_asr_fc2
        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: atom1, atom2, direction1, direction2
        TYPE(fc2_value),ALLOCATABLE,DIMENSION(:) :: fc2_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR fix for fc2================'

        ALLOCATE(fc2_sum(atom_number))

        DO i=1, atom_number
            fc2_sum(i)%phi = 0d0
        END DO

        DO i=1, SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz
            fc2_sum(atom1)%phi(direction1,direction2) = fc2_sum(atom1)%phi(direction1,direction2) +&
            &myfc2_value(atom1,atom2)%phi(direction1,direction2)
        END DO

        DO i=1, atom_number
        DO direction1=1,d
        DO direction2=1,d
            IF(ABS(fc2_sum(i)%phi(direction1,direction2)).gt.resolution) THEN
                WRITE(34,*)"asr for atom1: ",i,"    xyz1: ",get_letter(direction1),&
                &"    xyz2: ",get_letter(direction2)
                WRITE(34,*)fc2_sum(i)%phi(direction1,direction2)
                WRITE(34,*)'...is broken, now fix it'
                myfc2_value(i,i)%phi(direction1,direction2) = myfc2_value(i,i)%phi(direction1,direction2)-&
                &fc2_sum(i)%phi(direction1,direction2)
            END IF
        END DO
        END DO
        END DO

        DEALLOCATE(fc2_sum)
!        CLOSE(34)
    END SUBROUTINE fix_asr_fc2

    SUBROUTINE fix_asr_fc3
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER :: atom1, atom2,atom3, direction1, direction2, direction3
        INTEGER :: atom1_idx, atom2_idx, atom3_idx
        TYPE(fc3_value),ALLOCATABLE,DIMENSION(:,:) :: fc3_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR fix for fc3================'

        ALLOCATE(fc3_sum(atom_number,tot_atom_number))

        !fix 3rd idx, i neq j
        DO i=1, atom_number
            DO j=1,tot_atom_number
                fc3_sum(i,j)%psi = 0d0
            END DO
        END DO

        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc3_unique_idx,atom1)

            atom2 = myfc3_index(i)%jatom_number
            IF(atom2.eq.atom1) CYCLE
            atom3 = myfc3_index(i)%katom_number
!            IF(atom3.eq.atom1) CYCLE
            direction1 = myfc3_index(i)%iatom_xyz
            direction2 = myfc3_index(i)%jatom_xyz
            direction3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2) .AND. ANY(fc3_unique_idx==atom3)) THEN
                atom2_idx = find_loc(fc3_unique_idx,atom2)
                atom3_idx = find_loc(fc3_unique_idx,atom3)
                fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) =&
                & fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) +&
                & myfc3_value(atom1_idx,atom2_idx,atom3_idx)%psi(direction1,direction2,direction3)
            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
            IF(ABS(fc3_sum(i,j)%psi(direction1,direction2,direction3)).gt.resolution) THEN

                WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,&
                &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
                &"    xyz3: ",get_letter(direction3)
                WRITE(34,*)fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'...is broken, now fix it'
                !here is the sum non-zero, it means i,j must belong to fc3_unique_idx
                atom1_idx = find_loc(fc3_unique_idx,i)
                atom2_idx = find_loc(fc3_unique_idx,j)
                myfc3_value(atom1_idx,atom2_idx,atom2_idx)%psi(direction1,direction2,direction3) = &
                & myfc3_value(atom1_idx,atom2_idx,atom2_idx)%psi(direction1,direction2,direction3) - &
                & fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'---------------------------------------'
            END IF
        END DO
        END DO
        END DO
        END DO
        END DO

        !fix 2nd idx, i neq j
        DO i=1, atom_number
            DO j=1,tot_atom_number
                fc3_sum(i,j)%psi = 0d0
            END DO
        END DO

        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc3_unique_idx,atom1)

            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number
            IF(atom3.eq.atom1) CYCLE
            direction1 = myfc3_index(i)%iatom_xyz
            direction2 = myfc3_index(i)%jatom_xyz
            direction3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2) .AND. ANY(fc3_unique_idx==atom3)) THEN
            atom2_idx = find_loc(fc3_unique_idx,atom2)
            atom3_idx = find_loc(fc3_unique_idx,atom3)
                fc3_sum(atom1,atom3)%psi(direction1,direction2,direction3) =&
                & fc3_sum(atom1,atom3)%psi(direction1,direction2,direction3) +&
                & myfc3_value(atom1_idx,atom2_idx,atom3_idx)%psi(direction1,direction2,direction3)
            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
            IF(ABS(fc3_sum(i,j)%psi(direction1,direction2,direction3)).gt.resolution) THEN

                WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,&
                &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
                &"    xyz3: ",get_letter(direction3)
                WRITE(34,*)fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'...is broken, now fix it'
                !here is the sum non-zero, it means i,j must belong to fc3_unique_idx
                atom1_idx = find_loc(fc3_unique_idx,i)
                atom2_idx = find_loc(fc3_unique_idx,j)
                myfc3_value(atom1_idx,atom1_idx,atom2_idx)%psi(direction1,direction2,direction3) = &
                & myfc3_value(atom1_idx,atom1_idx,atom2_idx)%psi(direction1,direction2,direction3) - &
                & fc3_sum(i,j)%psi(direction1,direction2,direction3)
                WRITE(34,*)'---------------------------------------'

            END IF
        END DO
        END DO
        END DO
        END DO
        END DO

        !fix the first idx, the Psi_iii term
        DO i=1, atom_number
            DO j=1,tot_atom_number
                fc3_sum(i,j)%psi = 0d0
            END DO
        END DO

        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc3_unique_idx,atom1)

            atom2 = myfc3_index(i)%jatom_number
            IF(atom2.ne.atom1) CYCLE
            atom3 = myfc3_index(i)%katom_number

            direction1 = myfc3_index(i)%iatom_xyz
            direction2 = myfc3_index(i)%jatom_xyz
            direction3 = myfc3_index(i)%katom_xyz
            IF(ANY(fc3_unique_idx==atom2) .AND. ANY(fc3_unique_idx==atom3)) THEN
                atom2_idx = find_loc(fc3_unique_idx,atom2)
                atom3_idx = find_loc(fc3_unique_idx,atom3)
                fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) =&
                & fc3_sum(atom1,atom2)%psi(direction1,direction2,direction3) +&
                & myfc3_value(atom1_idx,atom2_idx,atom3_idx)%psi(direction1,direction2,direction3)
            END IF
        END DO

        DO i=1, atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
            IF(ABS(fc3_sum(i,i)%psi(direction1,direction2,direction3)).gt.resolution) THEN

                WRITE(34,*)"asr for atom1: ",i,"    atom2: ",i,&
                &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
                &"    xyz3: ",get_letter(direction3)
                WRITE(34,*)fc3_sum(i,i)%psi(direction1,direction2,direction3)
                WRITE(34,*)'...is broken, now fix it'
                !here is the sum non-zero, it means i,j must belong to fc3_unique_idx
                atom1_idx = find_loc(fc3_unique_idx,i)
                myfc3_value(atom1_idx,atom1_idx,atom1_idx)%psi(direction1,direction2,direction3) = &
                & myfc3_value(atom1_idx,atom1_idx,atom1_idx)%psi(direction1,direction2,direction3) - &
                & fc3_sum(i,i)%psi(direction1,direction2,direction3)
                WRITE(34,*)'---------------------------------------'

            END IF
        END DO
        END DO
        END DO
        END DO

        DEALLOCATE(fc3_sum)
!        CLOSE(34)
    END SUBROUTINE fix_asr_fc3

    SUBROUTINE fix_asr_fc4
        IMPLICIT NONE
        INTEGER :: i,j,k
        INTEGER :: atom1, atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: atom1_idx, atom2_idx, atom3_idx, atom4_idx
        TYPE(fc4_value),ALLOCATABLE,DIMENSION(:,:,:) :: fc4_sum

!        OPEN(34,FILE='output.txt',STATUS='old',action='write',POSITION='APPEND')
        WRITE(34,*)'===========ASR fix for fc4================'

        ALLOCATE(fc4_sum(atom_number,tot_atom_number,tot_atom_number))

        !fix 4th idx, j,k .ne. i and j.ne.k
        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO

        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            IF(atom2.eq.atom1) CYCLE
            atom3 = myfc4_index(i)%katom_number
            IF(atom3.eq.atom2 .or. atom3.eq.atom1) CYCLE
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
            IF(j.eq.i) CYCLE
        DO k=1, tot_atom_number
            IF(k.eq.i .or. k.eq.j) CYCLE
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom3: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)' ...is broken, now fix it'
            WRITE(34,*)'---------------------------------------'
            !here is the sum non-zero, it means i,j,k must belong to fc4_unique_idx
            atom1_idx = find_loc(fc4_unique_idx,i)
            atom2_idx = find_loc(fc4_unique_idx,j)
            atom3_idx = find_loc(fc4_unique_idx,k)
            myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) =&
            & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) - &
            & fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !fix 3rd idx,similar
        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO

        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            IF(atom2.eq.atom1) CYCLE
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number
            IF(atom4.eq.atom1 .or. atom4.eq.atom2) CYCLE
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom4)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
            IF(j.eq.i) CYCLE
        DO k=1, tot_atom_number
            IF(k.eq.i .or. k.eq.j) CYCLE
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom3: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)' ...is broken, now fix it'
            WRITE(34,*)'---------------------------------------'
            !here is the sum non-zero, it means i,j,k must belong to fc4_unique_idx
            atom1_idx = find_loc(fc4_unique_idx,i)
            atom2_idx = find_loc(fc4_unique_idx,j)
            atom3_idx = find_loc(fc4_unique_idx,k)
            myfc4_value(atom1_idx,atom2_idx,atom2_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) =&
            & myfc4_value(atom1_idx,atom2_idx,atom2_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) - &
            & fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !fix intermediate terms Chi_ijjj
        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO

        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            IF(atom2.eq.atom1) CYCLE
            atom3 = myfc4_index(i)%katom_number
            IF(atom3.eq.atom1 .or. atom3.ne.atom2) CYCLE
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
            IF(j.eq.i) CYCLE
        DO k=1, tot_atom_number
            IF(k.eq.i .or. k.ne.j) CYCLE
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom3: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)' ...is broken, now fix it'
            WRITE(34,*)'---------------------------------------'
            !here is the sum non-zero, it means i,j,k must belong to fc4_unique_idx
            atom1_idx = find_loc(fc4_unique_idx,i)
            atom2_idx = find_loc(fc4_unique_idx,j)
!            atom3_idx = find_loc(fc4_unique_idx,k)
            myfc4_value(atom1_idx,atom2_idx,atom2_idx,atom2_idx)%chi(direction1,direction2,direction3,direction4) =&
            & myfc4_value(atom1_idx,atom2_idx,atom2_idx,atom2_idx)%chi(direction1,direction2,direction3,direction4) - &
            & fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !fix 2nd idx,similar
        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO

        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            IF(atom3.eq.atom1) CYCLE
            atom4 = myfc4_index(i)%latom_number
            IF(atom4.eq.atom1) CYCLE
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom3,atom4)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO j=1,tot_atom_number
            IF(j.eq.i) CYCLE
        DO K=1, tot_atom_number
            IF(k.eq.i) CYCLE
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",j,"    atom3: ",k,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)' ...is broken, now fix it'
            WRITE(34,*)'---------------------------------------'
            !here is the sum non-zero, it means i,j,k must belong to fc4_unique_idx
            atom1_idx = find_loc(fc4_unique_idx,i)
            atom2_idx = find_loc(fc4_unique_idx,j)
            atom3_idx = find_loc(fc4_unique_idx,k)
            myfc4_value(atom1_idx,atom1_idx,atom2_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) =&
            & myfc4_value(atom1_idx,atom1_idx,atom2_idx,atom3_idx)%chi(direction1,direction2,direction3,direction4) - &
            & fc4_sum(i,j,k)%chi(direction1,direction2,direction3,direction4)
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        !fix 1st index, the Chi_iiii terms
        DO i=1, atom_number
            DO j=1,tot_atom_number
                DO k=1,tot_atom_number
                    fc4_sum(i,j,k)%chi = 0d0
                END DO
            END DO
        END DO

        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom1_idx = find_loc(fc4_unique_idx,atom1)

            atom2 = myfc4_index(i)%jatom_number
            IF(atom2.ne.atom1) CYCLE
            atom3 = myfc4_index(i)%katom_number
            IF(atom3.ne.atom2) CYCLE
            atom4 = myfc4_index(i)%latom_number
            direction1 = myfc4_index(i)%iatom_xyz
            direction2 = myfc4_index(i)%jatom_xyz
            direction3 = myfc4_index(i)%katom_xyz
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc4_unique_idx==atom2) .AND. ANY(fc4_unique_idx==atom3) &
            &.AND. ANY(fc4_unique_idx==atom4)) THEN

            atom2_idx = find_loc(fc4_unique_idx,atom2)
            atom3_idx = find_loc(fc4_unique_idx,atom3)
            atom4_idx = find_loc(fc4_unique_idx,atom4)
            fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) =&
            & fc4_sum(atom1,atom2,atom3)%chi(direction1,direction2,direction3,direction4) +&
    & myfc4_value(atom1_idx,atom2_idx,atom3_idx,atom4_idx)%chi(direction1,direction2,direction3,direction4)

            END IF
        END DO

        DO i=1, atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO direction3=1,d
        DO direction4=1,d

        IF(ABS(fc4_sum(i,i,i)%chi(direction1,direction2,direction3,direction4)).gt.resolution) THEN
            WRITE(34,*)"asr for atom1: ",i,"    atom2: ",i,"    atom3: ",i,&
            &"    xyz1: ",get_letter(direction1),"    xyz2: ",get_letter(direction2),&
            &"    xyz3: ",get_letter(direction3),"    xyz4: ",get_letter(direction4)
            WRITE(34,*)fc4_sum(i,i,i)%chi(direction1,direction2,direction3,direction4)
            WRITE(34,*)' ...is broken, now fix it'
            WRITE(34,*)'---------------------------------------'
            !here is the sum non-zero, it means i,j,k must belong to fc4_unique_idx
            atom1_idx = find_loc(fc4_unique_idx,i)
            myfc4_value(atom1_idx,atom1_idx,atom1_idx,atom1_idx)%chi(direction1,direction2,direction3,direction4) =&
            & myfc4_value(atom1_idx,atom1_idx,atom1_idx,atom1_idx)%chi(direction1,direction2,direction3,direction4) - &
            & fc4_sum(i,i,i)%chi(direction1,direction2,direction3,direction4)
        END IF

        END DO
        END DO
        END DO
        END DO
        END DO !7 nested loops

        DEALLOCATE(fc4_sum)
!        CLOSE(34)
    END SUBROUTINE fix_asr_fc4
!========================================================================================================
    SUBROUTINE kvector_Update(nband,nk)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: nband,nk
        REAL(8),ALLOCATABLE:: eigen_temp(:,:),integrate_dos(:),total_dos(:),afunc_dos(:),junk(:)
        REAL(8),DIMENSION(3,3) :: prim2cart

        CALL allocate_tetra(nc1,nc2,nc3,wmesh)

        IF(ALLOCATED(kpc)) DEALLOCATE(kpc)
        IF(ALLOCATED(wk)) DEALLOCATE(wk)
        IF(ALLOCATED(eigen_temp)) DEALLOCATE(eigen_temp)
        IF(ALLOCATED(total_dos)) DEALLOCATE(total_dos)
        IF(ALLOCATED(integrate_dos)) DEALLOCATE(integrate_dos)
        IF(ALLOCATED(afunc_dos)) DEALLOCATE(afunc_dos)
        IF(ALLOCATED(junk)) DEALLOCATE(junk)

        ALLOCATE(kpc(3,nk),wk(nk),eigen_temp(nband,nk) ,total_dos(wmesh),integrate_dos(wmesh),afunc_dos(wmesh),junk(nk))
        CALL make_kp_reg_tet ! shft should be initialized (nx,ny,nz,sx,sy,sz,kpc,wkt)

        !k points are stored in kpc(direction,ith), map to kvector(ith)%component(direction)
        CALL allocatek(nk)
        CALL tet_map(kvector)
        CALL get_weights(nk,kpc)

        !Another get_weights method
!        prim2cart(:, 1) = r01%component(:)
!        prim2cart(:, 2) = r02%component(:)
!        prim2cart(:, 3) = r03%component(:)
!        CALL get_weights3(nk,kpc,prim2cart,nibz)

        kp_gamma = kpc(:,1)

        WRITE(34,*)"kvector_Update successfully called"
    END SUBROUTINE kvector_Update
 !---------------------------------------------------------------------------
    SUBROUTINE kvector_Shift !shift all kvector_components by a random shift
        IMPLICIT NONE
        INTEGER :: i, shift_seed
        REAL(8) :: shift_range

        shift_range = 0.01
        shift_seed = 2021
        CALL srand(shift_seed)

        DO i=1, SIZE(kpc,DIM=2)
            kpc(:,i) = kpc(:,i) + shift_range*2*(/rand(),rand(),rand()/)-shift_range*(/1d0,1d0,1d0/)
        END DO

        DO i=1,SIZE(kvector)
            kvector(i)%component(:) = kvector(i)%component(:) + &
            &shift_range*2*(/rand(),rand(),rand()/)-shift_range*(/1d0,1d0,1d0/)
        END DO

        kp_gamma = kpc(:,1)
    END SUBROUTINE kvector_Shift
  !---------------------------------------------------------------------------
    SUBROUTINE check_degeneracy !check two successive eivals for each k, if difference small
        IMPLICIT NONE
        INTEGER :: k,l
        REAL(8) :: diff,threshold

        threshold = 1d-8

        OPEN(58,FILE='degeneracy.dat',STATUS='unknown',ACTION='write')
        WRITE(58,*) '==== lambda, k, eigen value ===='
        DO k=1,SIZE(eivals, dim=2)
            DO l=1,SIZE(eivals,dim=1)-1
                diff = ABS(eivals(l,k)-eivals(l+1,k))
                IF(diff.lt.threshold) THEN
                    WRITE(58,*)'lambda=',l,'&',l+1,'  k=',k,'  eivals=',eivals(l,k),eivals(l+1,k)
                END IF
            END DO
        END DO

        CLOSE(58)
    END SUBROUTINE check_degeneracy
  !---------------------------------------------------------------------------
    SUBROUTINE special_check !manually check the inverse fourier transform, 1d
        IMPLICIT NONE
        INTEGER :: i,k
        REAL(8),DIMENSION(:),ALLOCATABLE :: fc_pseudo,R_mesh,q_mesh
        COMPLEX(8),DIMENSION(:),ALLOCATABLE :: dm_pseudo,fc_check
        ALLOCATE(R_mesh(10),q_mesh(17),fc_pseudo(10),dm_pseudo(17),fc_check(10))
        !manually set fc, read from fc2.dat
        fc_pseudo(1)=10d0;fc_pseudo(2)=2d0;fc_pseudo(3)=0.5d0;fc_pseudo(4:10)=0d0
        !place the R points, from 1 to 10
        DO i=1,SIZE(R_mesh)
            R_mesh(i) = i-1
        END DO
        !generate q points, accordingly
        DO k=1,SIZE(q_mesh)
            q_mesh(k) = (k-1)*2*pi/SIZE(q_mesh)
        END DO
        !----------------------------------
        dm_pseudo(:) = 0d0
        !calculate dynmat pseudo from above
        DO k=1,SIZE(q_mesh)
            DO i=1,SIZE(R_mesh)
                dm_pseudo(k) = dm_pseudo(k) + fc_pseudo(i)*EXP(ci*q_mesh(k)*R_mesh(i))
            END DO
        END DO

        fc_check(:) = 0d0
        !inverse fourier transform to get 'calculated fc2'
        DO i=1,SIZE(R_mesh)
            DO k=1,SIZE(q_mesh)
                fc_check(i) = fc_check(i) + dm_pseudo(k)*EXP(-ci*q_mesh(k)*R_mesh(i))/SIZE(q_mesh)
            END DO
        END DO

        !output
        OPEN(74,FILE='simple_check.dat')
        DO i=1,SIZE(R_mesh)
            WRITE(74,*)'fc_pseudo=',fc_pseudo(i),'fc_calculated=',fc_check(i)
        END DO
        DO k=1,SIZE(q_mesh)
            WRITE(74,*)'q point=',q_mesh(k),'dynmat=',dm_pseudo(k)
        END DO
        CLOSE(74)

    END SUBROUTINE special_check
!========================================================================================================

    SUBROUTINE all_Update
        IMPLICIT NONE
        WRITE(34,*)'====================ENTER UPDATE ROUTINES==================='
        CALL trans_Update
        CALL convs_Update(conv_to_cart)
        CALL latticeparameters_Update(conv_to_cart)
        CALL primitivelattices_Update(conv_to_cart)
        CALL atompos0_Update(conv_to_cart)
!-------------------------------------------------------------------------------------
        CALL fcinit_Update
    !=======increase maxneighbors to keep total atom number same========condition old
!        DO WHILE(natoms.lt.tot_atom_number)
!            maxneighbors = maxneighbors + 1
!            WRITE(*,*)"Now increase maxneighbors to: ", maxneighbors
!            WRITE(34,*)"Now increase maxneighbors to: ", maxneighbors
!            CALL fcinit_Update
!            WRITE(34,*)"the natoms is now=",natoms
!            WRITE(34,*)"the tot_atom_number=",tot_atom_number
!        END DO
!        WRITE(34,*)"natoms loop exit!"
    !=======increase maxneighbors to keep R_0 value same========condition new
        CALL get_maxneighbors
        WRITE(*,*) 'Now maxneighbors = ', maxneighbors
        CALL fcinit_Update

        CALL atompos_Update(conv_to_cart)
        WRITE(34,*)"atompos_Update successfully called"
    !=======increase nshell(:,:) to include initial selected atoms for fc2,fc3,fc4========
        CALL get_nshells
        WRITE(34,*)"get_nshells successfully called"
!--------------------------------------------------------------------------------------------
        CALL cellvec_Update
        WRITE(34,*)"cellvec_Update successfully called"
    !=======fc update=======
        !old method:
!        CALL fc4_update
!        CALL fc3_update
!        CALL fc2_update
!        CALL fc4_update_check
!        CALL fc3_update_check
        !new method:
        CALL all_fc_update

        WRITE(34,*)'====================FINISH UPDATE ROUTINES==================='

    END SUBROUTINE all_Update
!===============================================================================================
END MODULE force_update
