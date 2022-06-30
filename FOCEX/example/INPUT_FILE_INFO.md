###### Example for input parameter file to run FOrce Constants EXtraction(FOCEX)

###### params.inp

1. `1 1 1 90 90 90`  &rarr; crystal structure x,y,z, $ \alpha $, $ \beta $, $ \gamma $

2. `0 0.5 0.5 0.5 0 0.5 0.5 0.5  0` &rarr; primitive lattice vector for example here in Ge

3. `5.762862` &rarr; lattice scale factor, must be consistent with POSCAR file

4. `9` &rarr; number of maximum shell to include in the force constant calculation from the supercell

5. `1 1 1 0` &rarr;  1st, 2nd, 3rd and 4th order force constants, 1 is include and 0 is do not include in the fitting

6. `1 0 0 0` &rarr;  flags for including translational and rotational and Huang invariance constraints, 1 is to include and 0 to not include

7. `0.0001 1.e-5`  &rarr; tolerance for equating (two coordinates in general), margin for eliminating a FC

8. `1d-9` &rarr;     svd cutoff for the smallest eigenvalue to be included

9. `1  .TRUE.`

10. `1` &rarr; number of types of atom

11. `72.64` &rarr; Mass of each individual atom

12. `Ge` &rarr; Name of atoms

13. `2` &rarr;    number of atoms in the primitive cell, coordinates are in reduced units (of cubic cell)

14. `1 1`

15. `5 5` &rarr; number of neighbors for calculating force constants

16. `1 1`

17. `1 1`

18. `1 1 0.00 0.00 0.00` &rarr; position of first atom in the primitive lattice

19. `2 1 0.25 0.25 0.25` &rarr; position of second atom in the primitive lattice

20. `300.00` &rarr; Temperature to evaluate force constant at, force constant files won’t be generated at this point but independent force constant displayed at the terminal

    

    ###### POSCAR1

    This file is the POSCAR file for supercell

    ###### OUTCAR1

    This file consists of position of displaced atom along x, y and z-direction and $F_x$, $F_y$ and $F_z$ along with the energy of each of the displaced structure. For example:

    `\# POSITION   TOTAL FORCE`

       `1    -289.18629538 =t, Etot(eV)`

      `2.89296000    0.00000000    2.88142999   -0.11758299    -0.00000000    -0.00000000`   

      `2.88142999    0.00000000    8.64428999    4.96000000E-004 -0.0000000    -0.00000000`   

      `2.88142999    5.76285999    2.88142999    4.96000000E-004 -0.00000000    -0.00000000 `  

      `2.88142999    5.76285999    8.64428999    -4.564000003E-003 -0.00000000    -0.00000000`

    `…               …                     …                 …        …             …`

    Here, the first line is the header for position and force for each atom in x, y and z direction. The second line consists of the energy of structure and the lines after second are positions (first three columns, x,y and z) and forces (the last three columns Fx, Fy and Fz).

    

    



