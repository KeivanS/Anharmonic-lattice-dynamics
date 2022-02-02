## Generate new FCs with expanded supercell

#### Some thoughts

1. Assume new atom lists are provided after expanding supercell, and this new atom list is stored in object `every_atom(:)`, while previous atom list(before expansion) is stored in another object `prev_atom(:)`, notice that these two better be seperated. So the whole atom list is `every_atom(:)` + `prev_atom(:)`
2. Assume after expansion, the original center unitcell should remain still the center unitcell(because in the computational scheme, the first atom should always be in the center unitcell, a.k.a only these force constants are used)
3. One major subroutines `Force_Match` 
4. `Force_Match` input arguments should include `rnk`,`beacon(:)`,where `beacon` is an object that has member `n1`,`n2`,`n3`,`atom_type` these four can fully locate an atom in the supercell. Size of array `beacon(:)` depends on `rnk`, i.e. `rnk=4` then `beacon(3)` 
5. Take FC2 as example, run a search through new atom list `every_atom(:)`, for every `ni`,`nj`,`nk`, do a customized match check, which is (a specific utility function needed here)

```fortran
{ni,nj,nk} == {n1,n2,n3} and atom_type == atom_type and xyz1 == xyz1 and xyz2 == xyz2
```

where n1, n2, n3 belongs to 2nd atom in FC2 indexes. For the newly generated FC2 to match, make sure the atom types match, the two xyz indexes also match(? not sure).

```fortran
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
```

