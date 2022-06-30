The example consists of force constants calculation of Ge. Inside Ge is `VASP_FILES` for evaluating the force constants using DFT and inside folder `FOCEX_example`  are the run for FOCEX for different neighbors 5 (Ge-11105), 6 (Ge-11106), 7 (Ge-11107), 8 (Ge-11108), 9 (Ge-11119) respectively. To run the FOCEX code navigate to any directory inside `FOCEX_example\Single_Displacement\Ge-1105`. 

- copy `params.inp`, `POSCAR1`,`OUTCAR1` and `fc234-13` binary to the directory
- run `fc234-14`

After successful run following output files should be available

- <fc1.dat>,<fc2.dat>,<fc3.dat>,<fc4.dat> are full list of force constants obtained from the setting
- <fc1_fit.dat>,<fc2_fit.dat>,<fc3_fit.dat>,<fc4_fit.dat> are independent set of force constants obtained from the setting params.inp
- <lat_fc.dat> contains some structure info of crystal and a  complete atom list of the supercell, stored as objects in this code
- <amatrx.dat> matrix file for solving `Ax=b` where `x` is the displacement of atoms, `b` the force from DFT calculation and `A` the force constant obtained by fitting using SVD
- <primlatt.xyz> crystal structure of primite cell interms of xyz format
- <latfc.xyz> same as that of <lat_fc.dat> but in xyz format
- <corresp.dat> mapping of the primitive lattice cell with the supercell
- <maps.dat> mapping of sign for the forces for each rank with the supercell
- <poscar.xyz> poscar file of the supercell for visualization
- <svd-all.dat> solution of force constant obtained from SVD with error and variance information on force constants.
- <log.dat> all the output logs written in this file
- <times.dat> runtime information of the code

Now to evaluate the phonon dispersion and density of states one needs to use the `THERMACOND` code. After FOCEX run is done, copy `params.born` , `params.phon`, `kpbs.in` files and `thermal-split` binary from `THERMACOND` to the working directory

- run `thermal-split` to obtain the phonon dispersion (phononlog.dat) and density of states(dos.dat), gruneisen (all_grun.dat) files
- To plot the phonon dispersion, and density of states file use the gnuplot file in the example