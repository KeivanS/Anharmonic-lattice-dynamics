MODULE trans
!!assumption: a very redundant large initial atom list that has more shells than needed
!!not used assumption: new atom list(appended)
!!info: r1,r2,r3,g1,g2,g3 are all vector type
!! r1, r2, r3 is my trans_vec(:,:)
!! r1.dot.g2 = 0; r1.dot.g3 = 0; r1.dot.g1 = 2*pi, etc.
    USE DFT_force_constants
    IMPLICIT NONE

    REAL(8),DIMENSION(d,d) :: trans_mat !transverse matrix acts on original cell
    TYPE(vector) :: new_r1, new_r2, new_r3, new_g1, new_g2, new_g3
    TYPE(PAtom), DIMENSION(:), ALLOCATABLE :: new_iatom !new primitive cell atom (species)
    INTEGER, DIMENSION(:), ALLOCATABLE :: ref_atoms !record the found ref atoms' old atom_idx

    TYPE(EAtom), DIMENSION(:),ALLOCATABLE :: added_atom !temporary list, dynamic
    TYPE(fc2_index),DIMENSION(:), ALLOCATABLE :: new_fc2_index
    TYPE(fc3_index),DIMENSION(:), ALLOCATABLE :: new_fc3_index
    TYPE(fc4_index),DIMENSION(:), ALLOCATABLE :: new_fc4_index

    interface operator(.match.)
        module procedure matchatom
    end interface

    interface operator(.merge.)
        module procedure addatom
    end interface

    interface operator(.minus.)
        module procedure minusatom
    end interface

    interface operator(.fcmatch.)
        module procedure matchfc2, matchfc3, matchfc4
    end interface

CONTAINS
!************************************************************************************
    function addatom(v, w) result(merged)
    !!add (n1,n2,n3) of atom w to atom v, atom type is of w
    !!v is one of the new primitive cell atom(ref atom)
    !!w is atom 2 in old fc2, for example
        type(EAtom),intent(in) :: v,w
        type(EAtom) :: merged



        merged%label_number = 0!whatever, will be modified manually
        merged%n1 = v%n1 + w%n1 !will adjust later
        merged%n2 = v%n2 + w%n2 !will adjust later
        merged%n3 = v%n3 + w%n3 !will adjust later
        merged%type_R = w%type_R!whatever, not used
        merged%type_tau = w%type_tau !will adjust later

        merged%R = merged%n1*r1%component + merged%n2*r2%component + merged%n3*r3%component !w.a.l
        merged%tau = w%tau !will adjust later, but notice R+tau is the absolute position
        merged%x = merged%R(1) + merged%tau(1) !exact
        merged%y = merged%R(2) + merged%tau(2) !exact
        merged%z = merged%R(3) + merged%tau(3) !exact

    end function addatom

    function minusatom(v, w) result(subtract)
        type(EAtom), intent(in) :: v, w
        type(EAtom) :: subtract
        real(8), dimension(d) :: r
        integer :: n1,n2,n3,tau

        subtract%n1 = v%n1 - w%n1
        subtract%n2 = v%n2 - w%n2
        subtract%n3 = v%n3 - w%n3
        subtract%type_tau = v%type_tau
        !below are whatever, not used
        subtract%label_number = v%label_number
        subtract%R = v%R
        subtract%type_R = v%type_R
        subtract%tau = v%tau
        subtract%x = v%x
        subtract%y = v%y
        subtract%z = v%z

    end function minusatom

    function matchatom(v, w) result(found)
        type(EAtom), intent(in) :: v, w
        logical :: found

        found = .true.
        if(v%n1.ne.w%n1) found = .false.
        if(v%n2.ne.w%n2) found = .false.
        if(v%n3.ne.w%n3) found = .false.
        if(v%type_tau.ne.w%type_tau) found = .false.
    end function matchatom

    function matchfc2(v, w) result(found)
        type(fc2_index), intent(in) :: v, w
        logical :: found

        found = .true.
        if(v%iatom_number .ne. w%iatom_number) found = .false.
        if(v%jatom_number .ne. w%jatom_number) found = .false.
        if(v%iatom_xyz .ne. w%iatom_xyz) found = .false.
        if(v%jatom_xyz .ne. w%jatom_xyz) found = .false.

    end function matchfc2

    function matchfc3(v, w) result(found)
        type(fc3_index), intent(in) :: v, w
        logical :: found

        found = .true.
        if(v%iatom_number .ne. w%iatom_number) found = .false.
        if(v%jatom_number .ne. w%jatom_number) found = .false.
        if(v%katom_number .ne. w%katom_number) found = .false.
        if(v%iatom_xyz .ne. w%iatom_xyz) found = .false.
        if(v%jatom_xyz .ne. w%jatom_xyz) found = .false.
        if(v%katom_xyz .ne. w%katom_xyz) found = .false.

    end function matchfc3

    function matchfc4(v, w) result(found)
        type(fc4_index), intent(in) :: v, w
        logical :: found

        found = .true.
        if(v%iatom_number .ne. w%iatom_number) found = .false.
        if(v%jatom_number .ne. w%jatom_number) found = .false.
        if(v%katom_number .ne. w%katom_number) found = .false.
        if(v%latom_number .ne. w%latom_number) found = .false.
        if(v%iatom_xyz .ne. w%iatom_xyz) found = .false.
        if(v%jatom_xyz .ne. w%jatom_xyz) found = .false.
        if(v%katom_xyz .ne. w%katom_xyz) found = .false.
        if(v%latom_xyz .ne. w%latom_xyz) found = .false.

    end function matchfc4
!===================================================================================
     SUBROUTINE get_reduced_coor(r,cx,cy,cz,r1,r2,r3,inside)
    !! for a given r-vector, it finds its reduced unit assuming it was
    !! created as: r=cx*r1+cy*r2+cz*r3
    !! if the integer variable inside=1 then r is inside the primcell, defined
    !! by [0:r1[,[0:r2[,[0:r3[ ; if inside=0 it's outside

         IMPLICIT NONE
         REAL(8), INTENT(in):: r(3)
         TYPE(vector),INTENT(in):: r1,r2,r3
         TYPE(vector)gg1,gg2,gg3
         INTEGER, INTENT(out):: inside
         REAL(8), INTENT(out):: cx,cy,cz
         REAL(8) epsl

         epsl=5d-10

         CALL make_reciprocal_lattice(r1,r2,r3,gg1,gg2,gg3)  ! no factor of 2pi
         inside = 1

         cx = (r.dot.gg1) + epsl
         if (cx.lt.0 .or. cx.ge.1) inside=0

         cy = (r.dot.gg2) + epsl
         if (cy.lt.0 .or. cy.ge.1) inside=0

         cz = (r.dot.gg3) + epsl
         if (cz.lt.0 .or. cz.ge.1) inside=0

     END SUBROUTINE get_reduced_coor
!============================================================
!------------------------------------------------------------------------------------
    SUBROUTINE find_ref_atom(mat)
    !! get new trans vec and find the new ref atom(center atom)
        IMPLICIT NONE
        INTEGER :: i, inside, j
        REAL(8),DIMENSION(d,d), INTENT(in) :: mat
        REAL(8) :: counter, r(3), cx, cy, cz
        TYPE(PAtom) :: temp_atom

        !get new R1, R2, R3, G1, G2, G3
        new_r1%component(:) = mat.dot.r1%component(:)
        new_r2%component(:) = mat.dot.r2%component(:)
        new_r3%component(:) = mat.dot.r3%component(:)

        CALL make_reciprocal_lattice_2pi(new_r1,new_r2,new_r3,new_g1,new_g2,new_g3)

        !if the Fortran2003 feature is not supported, comment off these
!        counter = 0
!        DO i=atom_number+1, SIZE(every_atom)
!            r(:) = every_atom(i)%R(:) + every_atom(i)%tau(:)
!            CALL get_reduced_coor(r,cx,cy,cz,new_r1,new_r2,new_r3,inside)
!            IF (inside.eq.1) counter = counter + 1
!        END DO
!        ALLOCATE(new_iatom(counter))
!        ALLOCATE(ref_atoms(counter))

        !loop over all atoms to find the ref atom
        counter = 0; cx = 0; cy = 0; cz = 0
        DO i=atom_number+1,SIZE(every_atom)
            r(:) = every_atom(i)%R(:) + every_atom(i)%tau(:)
            CALL get_reduced_coor(r,cx,cy,cz,new_r1,new_r2,new_r3,inside)
            IF (inside.eq.1) THEN
                counter = counter + 1
                temp_atom%atom_type = every_atom(i)%type_tau + counter !assign new tau label
                temp_atom%mass = iatom(every_atom(i)%type_tau)%mass
                temp_atom%name = iatom(every_atom(i)%type_tau)%name
                temp_atom%pos_tau = every_atom(i)%R + every_atom(i)%tau !new tau position
                temp_atom%ttyp = iatom(every_atom(i)%type_tau)%ttyp
                !if there is F03 support, use these
                IF(.not.ALLOCATED(new_iatom)) THEN
                    ALLOCATE(new_iatom(1),ref_atoms(1))
                    new_iatom(1) = temp_atom
                    ref_atoms(1) = i
                ELSE
                    !append object and index,
                    new_iatom = [new_iatom,temp_atom]
                    ref_atoms = [ref_atoms,i]
                END IF
                !if there is no F03 support
!                new_iatom(counter) = temp_atom
!                ref_atoms(counter) = i

            ELSE
                CYCLE
            END IF
        END DO

        !allocate and record new_iatom
        IF(counter.eq.0) THEN
            WRITE(*,*) 'cannot find new ref atom'
            RETURN
        END IF
    END SUBROUTINE find_ref_atom
!------------------------------------------------------------------------------------
    SUBROUTINE scanForFC2(ref,tar,fc_old)
    !!do a scan over all old atom list k
    !!try to find the atom tar = atom k - atom ref
    !!if not find, generate a new atom k' that satisfies atom tar = atom k' - atom ref
    !!and append this k' to the added_atom(:)
    !!duplicate object info from fc_idx to fc_dup
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ref, tar
        TYPE(fc2_index),INTENT(in) :: fc_old
        TYPE(fc2_index) :: fc_dup

        INTEGER :: k,idx,i
        LOGICAL :: found, duplicated
        TYPE(EAtom) :: temp_atom

        TYPE(fc2_index),ALLOCATABLE,DIMENSION(:) :: temp_fc2
        TYPE(EAtom),ALLOCATABLE,DIMENSION(:) :: temp_added

loop:   DO k=1, tot_atom_number
            temp_atom = every_atom(k).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar)

            IF(found) THEN
                !everything matches here, the 2nd atom label is k
                fc_dup%jatom_number = k
                EXIT loop
            ELSEIF(k.eq.tot_atom_number.AND.(.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar)
                !but before create, firstly check if it's already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)
                        IF(found) THEN
                            fc_dup%jatom_number = added_atom(idx)%label_number
                            EXIT loop
                        END IF
                    END DO
                END IF

!                temp_atom = every_atom(ref).merge.every_atom(tar)
!WRITE(ulog,*)'size of this not/allocated added_atom',SIZE(added_atom),ALLOCATED(added_atom)
                !intrinsic ERROR!!: in Fortran, even for an unallocated array, it still has SIZE 1
!                temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1 !can't put it here due to weird error


                IF(.not.ALLOCATED(added_atom)) THEN
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom

                ELSE
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 2nd atom label is this newly generated atom label
                fc_dup%jatom_number = temp_atom%label_number
WRITE(ulog,*)'added atom index',fc_dup%jatom_number
                EXIT loop !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and k is less than tot_atom_number, go to next k
            END IF
        END DO loop


        !duplicate other info to this fc_dup object
        fc_dup%group = fc_old%group
        fc_dup%iatom_number = ref
        fc_dup%iatom_xyz = fc_old%iatom_xyz
        fc_dup%jatom_xyz = fc_old%jatom_xyz
        fc_dup%phi_temp = fc_old%phi_temp

!************************************************************************************
        !there is one scenario, the fc_dup is already included in old fcs
        !check both in old fcs and newly added fcs to see duplicated or not
        duplicated = .false.
oldfcs: DO i=1,SIZE(myfc2_index)
            IF(fc_dup.fcmatch.myfc2_index(i)) THEN
                duplicated = .true.
                EXIT oldfcs
            END IF
        END DO oldfcs

newfcs: IF(ALLOCATED(new_fc2_index).AND. (.not.duplicated)) THEN
            DO i=1, SIZE(new_fc2_index)
                IF(fc_dup.fcmatch.new_fc2_index(i)) THEN
                    duplicated = .true.
                    EXIT newfcs
                END IF
            END DO
        END IF newfcs
!************************************************************************************
        IF(.not.duplicated) THEN

            !add new fc2
            IF(.not.ALLOCATED(new_fc2_index)) THEN
                ALLOCATE(new_fc2_index(1))
                new_fc2_index(1) = fc_dup
            ELSE
                !simple append
                new_fc2_index = [new_fc2_index, fc_dup]
                !if F03 supported, I don't have to do this
    !            ALLOCATE(temp_fc2, source=new_fc2_index)
    !            DEALLOCATE(new_fc2_index)
    !            ALLOCATE(new_fc2_index(SIZE(temp_fc2)+1))
    !            new_fc2_index = (/temp_fc2,fc_dup/)
    !            DEALLOCATE(temp_fc2)
            END IF
        ENDIF

    END SUBROUTINE scanForFC2
!------------------------------------------------------------------------------------
    SUBROUTINE scanForFC3(ref,tar1,tar2,fc_old)
    !!similar with fc2, but more layers of complexity we have two target here
    !!tar1 -> atom2, tar2 -> atom3
        IMPLICIT NONE
        INTEGER,INTENT(in) :: ref, tar1, tar2
        TYPE(fc3_index),INTENT(in) :: fc_old
        TYPE(fc3_index) :: fc_dup

        INTEGER :: k, l, idx, i
        LOGICAL :: found, duplicated
        TYPE(EAtom) :: temp_atom

        TYPE(fc3_index),ALLOCATABLE,DIMENSION(:) :: temp_fc3
        TYPE(EAtom),ALLOCATABLE,DIMENSION(:) :: temp_added

        !first, match atom2 (tar1) to fc_dup%jatom_number
loopk:   DO k=1, tot_atom_number
            temp_atom = every_atom(k).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar1)

            IF(found) THEN
                !everything matches here, the 2nd atom label is k
                fc_dup%jatom_number = k
                EXIT loopk
            ELSEIF(k.eq.tot_atom_number .AND. (.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar1)
                !but before create, firstly check if it is already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)

!WRITE(ulog,*)'tar fc3',fc_old%iatom_number,fc_old%jatom_number,fc_old%katom_number
!WRITE(ulog,*)'tar1 atom',every_atom(tar1)%n1,every_atom(tar1)%n2,every_atom(tar1)%n3,every_atom(tar1)%type_tau
!WRITE(ulog,*)'added atom',added_atom(idx)%n1,added_atom(idx)%n2,added_atom(idx)%n3,added_atom(idx)%type_tau
!WRITE(ulog,*) found
                        IF(found) THEN
                            fc_dup%jatom_number = added_atom(idx)%label_number
                            EXIT loopk
                        END IF
                    END DO
                END IF



                IF(.not.ALLOCATED(added_atom)) THEN
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom
                ELSE
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 2nd atom label is this newly generated atom label
                fc_dup%jatom_number = temp_atom%label_number
                EXIT loopk !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and k is less than tot_atom_number, go to next k
            END IF
        END DO loopk

        !second, match atom3(tar2) to fc_dup%katom_number
loopl:   DO l=1, tot_atom_number
            temp_atom = every_atom(l).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar2)

            IF(found) THEN
                !everything matches here, the 3rd atom label is l
                fc_dup%katom_number = l
                EXIT loopl
            ELSEIF(l.eq.tot_atom_number .AND. (.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar2)
                !but before create, firstly check if it is already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)
                        IF(found) THEN
                            fc_dup%katom_number = added_atom(idx)%label_number
                            EXIT loopl
                        END IF
                    END DO
                END IF

                IF(.not.ALLOCATED(added_atom)) THEN
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom
                ELSE
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 3rd atom label is this newly generated atom label
                fc_dup%katom_number = temp_atom%label_number
                EXIT loopl !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and l is less than tot_atom_number, go to next l
            END IF
        END DO loopl

        !duplicate other info to this fc_dup object
        fc_dup%group = fc_old%group
        fc_dup%iatom_number = ref
        fc_dup%iatom_xyz = fc_old%iatom_xyz
        fc_dup%jatom_xyz = fc_old%jatom_xyz
        fc_dup%katom_xyz = fc_old%katom_xyz
        fc_dup%psi_temp = fc_old%psi_temp

!************************************************************************************
        !there is one scenario, the fc_dup is already included in old fcs
        !check both in old fcs and newly added fcs to see duplicated or not
        duplicated = .false.
oldfcs: DO i=1,SIZE(myfc3_index)
            IF(fc_dup.fcmatch.myfc3_index(i)) THEN
                duplicated = .true.
                EXIT oldfcs
            END IF
        END DO oldfcs

newfcs: IF(ALLOCATED(new_fc3_index).AND. (.not.duplicated)) THEN
            DO i=1, SIZE(new_fc3_index)
                IF(fc_dup.fcmatch.new_fc3_index(i)) THEN
                    duplicated = .true.
                    EXIT newfcs
                END IF
            END DO
        END IF newfcs
!************************************************************************************
        IF(.not.duplicated) THEN

            IF(.not.ALLOCATED(new_fc3_index)) THEN
                ALLOCATE(new_fc3_index(1))
                new_fc3_index(1) = fc_dup
            ELSE
                !simple append
                new_fc3_index = [new_fc3_index, fc_dup]
                !if F03 supported, I don't have to do this
    !            ALLOCATE(temp_fc3, source=new_fc3_index)
    !            DEALLOCATE(new_fc3_index)
    !            ALLOCATE(new_fc3_index(SIZE(temp_fc3)+1))
    !            new_fc3_index = (/temp_fc3,fc_dup/)
    !            DEALLOCATE(temp_fc3)
            END IF

        END IF

    END SUBROUTINE scanForFC3
!------------------------------------------------------------------------------------
    SUBROUTINE scanForFC4(ref, tar1, tar2, tar3, fc_old)
        !!similar with fc3, but more layers of complexity we have two target here
        !!tar1 -> atom2, tar2 -> atom3, tar3 -> atom4
        IMPLICIT NONE
        INTEGER,INTENT(in) :: ref, tar1, tar2, tar3
        TYPE(fc4_index),INTENT(in) :: fc_old
        TYPE(fc4_index) :: fc_dup

        INTEGER :: k, l, m, idx, i
        LOGICAL :: found, duplicated
        TYPE(EAtom) :: temp_atom

        !if F03 supported, I don't have to do this
        TYPE(fc4_index),ALLOCATABLE,DIMENSION(:) :: temp_fc4
        TYPE(EAtom),ALLOCATABLE,DIMENSION(:) :: temp_added


        !first, match atom2 (tar1) to fc_dup%jatom_number
loopk:   DO k=1, tot_atom_number
            temp_atom = every_atom(k).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar1)

            IF(found) THEN
                !everything matches here, the 2nd atom label is k
                fc_dup%jatom_number = k
                EXIT loopk
            ELSEIF(K.eq.tot_atom_number.AND.(.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar1)
                !but before create, firstly check if it is already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)
                        IF(found) THEN
                            fc_dup%jatom_number = added_atom(idx)%label_number
                            EXIT loopk
                        END IF
                    END DO
                END IF

                IF(.not.ALLOCATED(added_atom)) THEN
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom
                ELSE
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 2nd atom label is this newly generated atom label
                fc_dup%jatom_number = temp_atom%label_number
                EXIT loopk !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and k is less than tot_atom_number, go to next k
            END IF
        END DO loopk

        !second, match atom3(tar2) to fc_dup%katom_number
loopl:   DO l=1, tot_atom_number
            temp_atom = every_atom(l).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar2)

            IF(found) THEN
                !everything matches here, the 3rd atom label is l
                fc_dup%katom_number = l
                EXIT loopl
            ELSEIF(l.eq.tot_atom_number .AND. (.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar2)
                !but before create, firstly check if it is already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)
                        IF(found) THEN
                            fc_dup%katom_number = added_atom(idx)%label_number
                            EXIT loopl
                        END IF
                    END DO
                END IF

                IF(.not.ALLOCATED(added_atom)) THEN
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom
                ELSE
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 3rd atom label is this newly generated atom label
                fc_dup%katom_number = temp_atom%label_number
                EXIT loopl !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and l is less than tot_atom_number, go to next l
            END IF
        END DO loopl

        !third, match atom4(tar3) to fc_dup%latom_number
loopm:  DO m=1, tot_atom_number
            temp_atom = every_atom(m).minus.every_atom(ref)
            found = temp_atom.match.every_atom(tar3)

            IF(found) THEN
                !everything matches here, the 3rd atom label is l
                fc_dup%latom_number = m
                EXIT loopm
            ELSEIF(m.eq.tot_atom_number .AND. (.not.found)) THEN
                !if it's already the last existed atom and still no found,manually create an atom
                temp_atom = every_atom(ref).merge.every_atom(tar3)
                !but before create, firstly check if it is already created before
                IF(ALLOCATED(added_atom)) THEN
                    DO idx=1,SIZE(added_atom)
                        found = temp_atom.match.added_atom(idx)
                        IF(found) THEN
                            fc_dup%latom_number = added_atom(idx)%label_number
                            EXIT loopm
                        END IF
                    END DO
                END IF

                IF(.not.ALLOCATED(added_atom)) THEN
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + 1
                    ALLOCATE(added_atom(1))
                    added_atom(1) = temp_atom
                ELSE
                    !current atom number
                    temp_atom%label_number = SIZE(every_atom) + SIZE(added_atom) + 1
                    !append new atom and append to a temporary list
                    added_atom = [added_atom,temp_atom]
                    !if F03 supported, I don't have to do this
!                    ALLOCATE(temp_added,source=added_atom)
!                    DEALLOCATE(added_atom)
!                    ALLOCATE(added_atom(SIZE(temp_added)+1))
!                    added_atom = (/temp_added,temp_atom/)
!                    DEALLOCATE(temp_added)
                END IF

                !the 3rd atom label is this newly generated atom label
                fc_dup%latom_number = temp_atom%label_number
                EXIT loopm !not really needed since it's already the last iteration
            ELSE
                CYCLE !only if not match and m is less than tot_atom_number, go to next m
            END IF
        END DO loopm

        !duplicate other info to this fc_dup object
        fc_dup%group = fc_old%group
        fc_dup%iatom_number = ref
        fc_dup%iatom_xyz = fc_old%iatom_xyz
        fc_dup%jatom_xyz = fc_old%jatom_xyz
        fc_dup%katom_xyz = fc_old%katom_xyz
        fc_dup%latom_xyz = fc_old%latom_xyz
        fc_dup%chi_temp = fc_old%chi_temp

!************************************************************************************
        !there is one scenario, the fc_dup is already included in old fcs
        !check both in old fcs and newly added fcs to see duplicated or not
        duplicated = .false.
oldfcs: DO i=1,SIZE(myfc4_index)
            IF(fc_dup.fcmatch.myfc4_index(i)) THEN
                duplicated = .true.
                EXIT oldfcs
            END IF
        END DO oldfcs

newfcs: IF(ALLOCATED(new_fc4_index).AND. (.not.duplicated)) THEN
            DO i=1, SIZE(new_fc4_index)
                IF(fc_dup.fcmatch.new_fc4_index(i)) THEN
                    duplicated = .true.
                    EXIT newfcs
                END IF
            END DO
        END IF newfcs
!************************************************************************************
        IF(.not.duplicated) THEN

            IF(.not.ALLOCATED(new_fc4_index)) THEN
                ALLOCATE(new_fc4_index(1))
                new_fc4_index(1) = fc_dup
            ELSE
                !simple append
                new_fc4_index = [new_fc4_index,fc_dup]
                !if F03 supported, I don't have to do this
    !            ALLOCATE(temp_fc4, source=new_fc4_index)
    !            DEALLOCATE(new_fc4_index)
    !            ALLOCATE(new_fc4_index(SIZE(temp_fc4)+1))
    !            new_fc4_index = (/temp_fc4,fc_dup/)
    !            DEALLOCATE(temp_fc4)
            END IF

        END IF
    END SUBROUTINE scanForFC4
!------------------------------------------------------------------------------------
    SUBROUTINE duplicate_FCs
    !!duplicate rank 2, 3, 4 FCs w.r.t new ref atoms
        IMPLICIT NONE
        INTEGER :: i,j,k,l,m
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: tau1,ntau1

        INTEGER :: atm_idx

        !rank 2
        DO i=1, SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            tau1 = every_atom(atom1)%type_tau
            atom2 = myfc2_index(i)%jatom_number

            DO j=1,SIZE(ref_atoms)
                atm_idx = ref_atoms(j)
                ntau1 = every_atom(atm_idx)%type_tau
                IF(ntau1.ne.tau1) CYCLE

                CALL scanForFC2(atm_idx,atom2,myfc2_index(i))

            END DO !ref atom set
        END DO !old fc2 loop

        !rank 3
        DO i=1, SIZE(myfc3_index)
            atom1 = myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            tau1 = every_atom(atom1)%type_tau
            atom2 = myfc3_index(i)%jatom_number
            atom3 = myfc3_index(i)%katom_number

            DO j=1,SIZE(ref_atoms)
                atm_idx = ref_atoms(j)
                ntau1 = every_atom(atm_idx)%type_tau
                IF(ntau1.ne.tau1) CYCLE

                CALL scanForFC3(atm_idx,atom2,atom3,myfc3_index(i))

            END DO !ref atom set
        END DO !old fc3 loop

        !rank 4
        DO i=1, SIZE(myfc4_index)
            atom1 = myfc4_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            tau1 = every_atom(atom1)%type_tau
            atom2 = myfc4_index(i)%jatom_number
            atom3 = myfc4_index(i)%katom_number
            atom4 = myfc4_index(i)%latom_number

            DO j=1,SIZE(ref_atoms)
                atm_idx = ref_atoms(j)
                ntau1 = every_atom(atm_idx)%type_tau
                IF(ntau1.ne.tau1) CYCLE

                CALL scanForFC4(atm_idx,atom2,atom3,atom4,myfc4_index(i))

            END DO !ref atom set
        END DO !old fc4 loop

    END SUBROUTINE duplicate_FCs
!------------------------------------------------------------------------------------
    SUBROUTINE combineEverything
    !!append does not work
    !!append added_atom(:), new_fc#_index
    !!also update number of atoms
        IMPLICIT NONE
        TYPE(PAtom),DIMENSION(:),ALLOCATABLE ::temp_iatom
        TYPE(EAtom),DIMENSION(:),ALLOCATABLE :: temp_atom
        TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: temp_fc2_index
        TYPE(fc3_index),DIMENSION(:),ALLOCATABLE :: temp_fc3_index
        TYPE(fc4_index),DIMENSION(:),ALLOCATABLE :: temp_fc4_index

        IF(ALLOCATED(new_iatom)) THEN
            ALLOCATE(temp_iatom, source=(/iatom, new_iatom/))
            DEALLOCATE(iatom)
            ALLOCATE(iatom(SIZE(temp_iatom)))
            iatom = temp_iatom
            atom_number = SIZE(iatom)
            DEALLOCATE(temp_iatom)
        END IF

        IF(ALLOCATED(added_atom)) THEN
            ALLOCATE(temp_atom, source=(/every_atom, added_atom/))
            DEALLOCATE(every_atom)
            ALLOCATE(every_atom(SIZE(temp_atom)))
            every_atom = temp_atom
            tot_atom_number = SIZE(every_atom)
            DEALLOCATE(temp_atom)
        END IF

        ALLOCATE(temp_fc2_index, source=(/myfc2_index, new_fc2_index/))
        DEALLOCATE(myfc2_index)
        ALLOCATE(myfc2_index(SIZE(temp_fc2_index)))
        myfc2_index = temp_fc2_index
        DEALLOCATE(temp_fc2_index)

        ALLOCATE(temp_fc3_index, source=(/myfc3_index, new_fc3_index/))
        DEALLOCATE(myfc3_index)
        ALLOCATE(myfc3_index(SIZE(temp_fc3_index)))
        myfc3_index = temp_fc3_index
        DEALLOCATE(temp_fc3_index)

        ALLOCATE(temp_fc4_index, source=(/myfc4_index, new_fc4_index/))
        DEALLOCATE(myfc4_index)
        ALLOCATE(myfc4_index(SIZE(temp_fc4_index)))
        myfc4_index = temp_fc4_index
        DEALLOCATE(temp_fc4_index)

    END SUBROUTINE combineEverything
!------------------------------------------------------------------------------------
    SUBROUTINE swapAtomLabel(i,j)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: i,j
        INTEGER :: temp_label

        temp_label = every_atom(i)%label_number
        every_atom(i)%label_number = every_atom(j)%label_number
        every_atom(j)%label_number = temp_label
    END SUBROUTINE swapAtomLabel

    SUBROUTINE swapAtom(i,j)!swap the object
        IMPLICIT NONE
        TYPE(EAtom) :: i,j
        TYPE(EAtom) :: temp_atom

        temp_atom = i
        i = j
        j = temp_atom

    END SUBROUTINE swapAtom

    SUBROUTINE swapForceLabel(i,j)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: i,j
        INTEGER :: idx
        INTEGER :: atom1, atom2, atom3, atom4

        DO idx=1,SIZE(myfc2_index)
            atom1 = myfc2_index(idx)%iatom_number
            atom2 = myfc2_index(idx)%jatom_number

            IF(atom1.eq.i) myfc2_index(idx)%iatom_number = j
            IF(atom1.eq.j) myfc2_index(idx)%iatom_number = i
            IF(atom2.eq.i) myfc2_index(idx)%jatom_number = j
            IF(atom2.eq.j) myfc2_index(idx)%jatom_number = i
        END DO

        DO idx=1,SIZE(myfc3_index)
            atom1 = myfc3_index(idx)%iatom_number
            atom2 = myfc3_index(idx)%jatom_number
            atom3 = myfc3_index(idx)%katom_number

            IF(atom1.eq.i) myfc3_index(idx)%iatom_number = j
            IF(atom1.eq.j) myfc3_index(idx)%iatom_number = i
            IF(atom2.eq.i) myfc3_index(idx)%jatom_number = j
            IF(atom2.eq.j) myfc3_index(idx)%jatom_number = i
            IF(atom3.eq.i) myfc3_index(idx)%katom_number = j
            IF(atom3.eq.j) myfc3_index(idx)%katom_number = i
        END DO

        DO idx=1,SIZE(myfc4_index)
            atom1 = myfc4_index(idx)%iatom_number
            atom2 = myfc4_index(idx)%jatom_number
            atom3 = myfc4_index(idx)%katom_number
            atom4 = myfc4_index(idx)%latom_number

            IF(atom1.eq.i) myfc4_index(idx)%iatom_number = j
            IF(atom1.eq.j) myfc4_index(idx)%iatom_number = i
            IF(atom2.eq.i) myfc4_index(idx)%jatom_number = j
            IF(atom2.eq.j) myfc4_index(idx)%jatom_number = i
            IF(atom3.eq.i) myfc4_index(idx)%katom_number = j
            IF(atom3.eq.j) myfc4_index(idx)%katom_number = i
            IF(atom4.eq.i) myfc4_index(idx)%latom_number = j
            IF(atom4.eq.j) myfc4_index(idx)%latom_number = i
        END DO
    END SUBROUTINE swapForceLabel
!------------------------------------------------------------------------------------
    FUNCTION ifAllInt(n1,n2,n3) RESULT(yes)
    !!see if n1, n2, n3 are all integer
        IMPLICIT NONE
        REAL(8), INTENT(in) :: n1, n2, n3
        LOGICAL :: yes
        REAL(8) :: epsl
        epsl = 1d-5 !?what value

        yes = .true.
        IF((n1-floor(n1)).gt.epsl) yes=.false.
        IF((n2-floor(n2)).gt.epsl) yes=.false.
        IF((n3-floor(n3)).gt.epsl) yes=.false.

    END FUNCTION ifAllInt

    SUBROUTINE get_typeAndmore(r,n1,n2,n3,tau)
    !!this should run after combineEverything
    !!use this after calling find_ref_atom(where new translational vectors calculated)
    !!for renew atom list
    !!for new translational vectors new_r1, new_r2, new_r3 and given atom position r(:)
    !!determine new_n1, new_n2, new_n3 and new_tau
        IMPLICIT NONE
        REAL(8), INTENT(in), DIMENSION(d) :: r
        INTEGER, INTENT(out) :: n1,n2,n3,tau
        REAL(8) :: temp_r(d)
        INTEGER :: i
        LOGICAL :: found
        REAL(8) :: cx, cy, cz
        INTEGER :: inside !no use here

        found = .false.
        DO i=1,atom_number
            temp_r = r - iatom(i)%pos_tau
            CALL get_reduced_coor(temp_r,cx,cy,cz,new_r1,new_r2,new_r3,inside)
            IF(.not.ifAllInt(cx,cy,cz)) THEN
                CYCLE !cx, cy, cz are not all integers
            ELSE
                found = .true. !found this atom inside old atom type
                tau = i
                n1 = floor(cx)
                n2 = floor(cy)
                n3 = floor(cz)
            END IF
        END DO

        IF(.not.found) WRITE(*,*) "something is wrong" !still can't find correct attype
    END SUBROUTINE get_typeAndmore
!------------------------------------------------------------------------------------
    SUBROUTINE relabelAtoms
    !!this should run after combineEverything
    !!with the new ref_atoms, atom types are expanded
    !!with the new translational vectors, every n1, n2, n3 should be updated
    !!swap atom label if needed for newly added primitive cell atoms

        IMPLICIT NONE
        INTEGER :: i
        INTEGER :: n1, n2, n3, tau
        REAL(8) :: r(d)

        DO i=1,SIZE(iatom)
            iatom(i)%atom_type = i
        END DO

        DO i=1, SIZE(every_atom)
            r = every_atom(i)%R + every_atom(i)%tau
            CALL get_typeAndmore(r,n1,n2,n3,tau)
            every_atom(i)%n1 = n1
            every_atom(i)%n2 = n2
            every_atom(i)%n3 = n3
            every_atom(i)%type_tau = tau

            every_atom(i)%R = n1*new_r1 + n2*new_r2 + n3*new_r3
            every_atom(i)%tau = iatom(tau)%pos_tau
            every_atom(i)%tau = r - every_atom(i)%R !these two should equal, can be used to check

            IF(n1.eq.0 .AND. n2.eq.0 .AND. n3.eq.0) THEN
                !this is a new primitive cell atom, switch its primitive atom number
                IF(every_atom(i)%label_number.ne.tau) THEN
                    !if it's equal then that's what we want(original primitive cell atoms auto satisfy this)
!                    CALL swapAtomLabel(i,tau)
                    CALL swapAtom(every_atom(i),every_atom(tau))
                    CALL swapForceLabel(i,tau)
                END IF
            END IF
        END DO

    END SUBROUTINE relabelAtoms
!======================================output subroutines========================================
    SUBROUTINE printParams
    !!modified based on <write_input_fit>
        IMPLICIT NONE
        integer i,counter,label,unit_params
        real(8) scal,junk
        character jjj*1

        real(8) a,b,c,alpha,beta,gama !calculate the new lattice parameters
        real(8) cx,cy,cz
        integer inside

        unit_params=49
        open(unit_params,file='structure_new.params')
        a = length(new_r1)
        b = length(new_r2)
        c = length(new_r3)
        alpha = 180/pi*ACOS((new_r2.dot.new_r3)/b/c)
        beta = 180/pi*ACOS((new_r1.dot.new_r3)/a/c)
        gama = 180/pi*ACOS((new_r1.dot.new_r2)/a/b)

        latticeparameters = (/a,b,c,alpha,beta,gama/)
        primitivelattice(1,:) = (/1,0,0/)
        primitivelattice(2,:) = (/0,1,0/)
        primitivelattice(3,:) = (/0,0,1/)
        lattice_parameter = 1d0

        write(unit_params,1) latticeparameters   ! (a,b,c,alpha,beta,gamma)
        write(unit_params,2) primitivelattice     ! prim latt vectors(3:3) in terms of conventional above
        write(unit_params,'(f2.0)') lattice_parameter !scale factor for lattparams (must be consistent with POSCAR data)
        scal = lattice_parameter
!        write(unit_params,'(i2)') maxneighbors  !used in lat_fc.dat
        write(unit_params,'(4(i1,1x))') include_fc    ! if=1 include this rank
    ! read(uparams,*) nshells       ! # of neigr-shells to use for each rank
        write(unit_params,'(4(i1,1x))')  itrans,irot,ihuang,enforce_inv !not used? flags for including translational and rotational invce constraints.
    ! read(uparams,*) junk
    ! read(uparams,*) junk
        write(unit_params,*) nouse1,nouse2 ! this really the flags...&
                                   !&tolerance for equating (two coordinates in general), margin for eliminating a FC
!        write(unit_params,'(es8.2)') svdcut   !svd cutoff for the smallest eigenvalue to be included
        write(unit_params,*) fdfiles,verbose     !number of force-displacement files
        write(unit_params,'(i1)') natom_type   ! # of different elements present in prim cell, I call it 'atom_number'
    !read(uparams,*) jjj
    ! allocate(natom(natom_type))

        write(unit_params,'(f6.2)') (mas(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]

        write(unit_params,'(a2)') (atname(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]

        !----------modified---------
        write(unit_params,'(i1)') atom_number        ! # of atoms in primitive cell,coordinates are in reduced units (of cubic cell)
     !nshells(4,20) are declared in module[params]
        WRITE(unit_params,*) nouse3
        do i=2, atom_number
            nshells(:,i) = nshells(:,1)
        end do

!        write(unit_params,'(20(i2,1x))') nshells(1,1:atom_number)       ! # of neigr-shells to use for each rank
        write(unit_params,'(20(i2,1x))') nshells(2,1:atom_number)       ! # of neigr-shells to use for each rank
        write(unit_params,'(20(i2,1x))') nshells(3,1:atom_number)       ! # of neigr-shells to use for each rank
        write(unit_params,'(20(i2,1x))') nshells(4,1:atom_number)       ! # of neigr-shells to use for each rank


        do i=1,atom_number
        ! positions here must be D format in conventional cell for use by fcs_init
            label = i
            CALL get_reduced_coor(iatom(i)%pos_tau,cx,cy,cz,new_r1,new_r2,new_r3,inside)
            write(unit_params,5) label,iatom(i)%ttyp,cx,cy,cz
            WRITE(*,*) 'writing  atom #',i, counter
        enddo

        close(unit_params)

1 FORMAT(3(f4.2,1x),3(f3.0,1x))
2 FORMAT(9(f2.0,1x))
5 FORMAT(2(i2,1x),3(f5.2,1x))
    END SUBROUTINE printParams
!--------------------------------------------------------------------------------------------------------------------
    SUBROUTINE printAtoms
    !!modified based on <write_lat_fc>
        IMPLICIT NONE
        real(8) rij
        integer i,j,ngrps(4),ntrms(4),j_sc,maxshell

        maxshell = 9 !hard coded for not knowing how is this different from maxhells = 27 (also hard coded)

        open(ufco,file='lat_fc_new.dat',status='unknown')

        write(ufco,*)' Crystal data: translation vectors of the primitive cell '
        write(ufco,9)new_r1
        write(ufco,9)new_r2
        write(ufco,9)new_r3
        write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
        write(ufco,*)atom_number
        do i=1,atom_number
            write(ufco,6) i,iatom(i)%name,iatom(i)%ttyp,iatom(i)%pos_tau,iatom(i)%mass
        enddo
        write(ufco,*)' Crystal data written ************************************'
        write(ufco,*)' Included ranks of FCs '
        write(ufco,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
        write(ufco,*)' Number of FCs for each rank '
        write(ufco,*) 0,SIZE(myfc2_index),SIZE(myfc3_index),SIZE(myfc4_index)
        write(ufco,*)' Number of independent FCs for each rank '
        write(ufco,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)
        write(ufco,*)' maxshells,Neighborshell atoms: i,x,y,z,type_tau,n1,n2,n3 '
        write(ufco,*)maxshell, tot_atom_number
        do i=1,tot_atom_number
            rij = length(every_atom(i)%R+every_atom(i)%tau-every_atom(1)%R-every_atom(1)%tau)
            write(ufco,7) i,(every_atom(i)%R(j)+every_atom(i)%tau(j),j=1,3),every_atom(i)%type_tau,&
            &every_atom(i)%n1,every_atom(i)%n2,every_atom(i)%n3,rij
!        write(ufco,7)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),rij !legacy
        enddo

        close(ufco)

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,3(2x,f19.10),i6)
9 format(9(2x,f19.10))
    END SUBROUTINE printAtoms
!--------------------------------------------------------------------------------------------------------------------
    SUBROUTINE printFCs
    !!modified based on <write_output_fc>
        IMPLICIT NONE

        integer rank,t,i,res,j,term2,g
        real(8) rij,bunit,one,dij,trace,fcd

        one = 1d0

        rank=2
        open(ufc2,file='fc2_new.dat',status='unknown',action='write')
        write(ufc2,*)'# RANK 2 tensors :term,group,(iatom,ixyz)_2 d2U/dx_{i,alpha}^2'
        do t=1,SIZE(myfc2_index)
            write(ufc2,2)t,myfc2_index(t)%group, &
            & myfc2_index(t)%iatom_number,myfc2_index(t)%iatom_xyz, &
            & myfc2_index(t)%jatom_number,myfc2_index(t)%jatom_xyz, &
            & myfc2_index(t)%phi_temp, one
        end do
        close(ufc2)

        rank=3
        open(ufc3,file='fc3_new.dat',status='unknown',action='write')
        write(ufc3,*)' RANK 3 tensors :term,group, (iatom,ixyz)_3 d3U/dx_{i,alpha}^3'
        do t=1,SIZE(myfc3_index)
            write(ufc3,3)t,myfc3_index(t)%group, &
            & myfc3_index(t)%iatom_number, myfc3_index(t)%iatom_xyz, &
            & myfc3_index(t)%jatom_number, myfc3_index(t)%jatom_xyz, &
            & myfc3_index(t)%katom_number, myfc3_index(t)%katom_xyz, &
            & myfc3_index(t)%psi_temp, one
        enddo
        close(ufc3)

        rank=4
        open(ufc4,file='fc4_new.dat',status='unknown',action='write')
        write(ufc4,*)' RANK 4 tensors :term, group, (iatom,ixyz)_4 d4U/dx_{i,alpha}^4'
        do t=1,SIZE(myfc4_index)
            write(ufc4,4)t,myfc4_index(t)%group, &
            & myfc4_index(t)%iatom_number, myfc4_index(t)%iatom_xyz, &
            & myfc4_index(t)%jatom_number, myfc4_index(t)%jatom_xyz, &
            & myfc4_index(t)%katom_number, myfc4_index(t)%katom_xyz, &
            & myfc4_index(t)%latom_number, myfc4_index(t)%latom_xyz, &
            & myfc4_index(t)%chi_temp, one
        enddo
        close(ufc4)

2 format(i6,1x,i5,2(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
3 format(i6,1x,i5,3(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
4 format(i6,1x,i5,4(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)

    END SUBROUTINE printFCs
!--------------------------------------------------------------------------------------------------------------------

END MODULE trans
