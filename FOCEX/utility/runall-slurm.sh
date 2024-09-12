#!/bin/bash
# calculate all poscar_* in parallel via job array
# submit a final job to combine the force-positions

# input files
CELL="cell.inp"
SNAPS="snaps.inp"
SUPERCELL="supercell.inp"

# Slurm parameters
JOB1="job.slurm"
JOB2="cat.slurm"
ACCOUNT="myaccount"
PARTITION="standard"
CORES=8

# check input files
if [ ! -e $CELL ]; then
    echo "You must provide $CELL"
    exit 1
fi

if [ ! -e $SNAPS ]; then
    echo "$SNAPS not found. Writing default template."
    cat >$SNAPS <<EOF
500    # Avg desired frequency (1/cm)  (500/cm = 15 THz)
300    # temperature in K for canonical sampling of displacements and velocities
15     # desired number of snapshots needed (typically 20-50 should suffice)
EOF
fi

if [ ! -e $SUPERCELL ]; then
    echo "$SUPERCELL not found. Writing default template."
    cat >$SUPERCELL <<EOF
3 0 0
0 3 0
0 0 3
EOF
fi

sc_snaps.x

# check VASP input
for i in INCAR KPOINTS POTCAR; do
    if [ ! -e $i ]; then
        echo "$i not found. Aborting."
        exit
    fi
done

N=$(ls poscar_*|wc -l)
N=$((N-1))

cat >$SLURM <<EOF
#!/bin/bash
#SBATCH -A $ACCOUNT
#SBATCH -p $PARTITION
#SBATCH -N 1 
#SBATCH --ntasks-per-node=$CORES
#SBATCH -t 1:0:0
#SBATCH -a 0-$N

i=$(printf %03g \$SLURM_ARRAY_TASK_ID)
mkdir -p \$i
cp INCAR KPOINTS POTCAR \$i
cp poscar_\$i \$i/POSCAR
cd \$i

module purge
module load vasp

# read vasp output OUTCAR and extract positions and forces into pos-forc.dat
srun vasp_std && read_outcar.x || exit 1
EOF

# submit job array and capture job ID
JOBID=$(sbatch --parsable $JOB1)
echo Submitted batch job $JOBID

cat >$JOB1 <<EOF
#!/bin/bash
#SBATCH -A $ACCOUNT
#SBATCH -p $PARTITION
#SBATCH -c 1
#SBATCH -t 0:10:0

# collect the force-positions from every snapshot into a single FORCEDISP1 file to be read by FOCEX
cat */pos-forc.dat >> FORCEDISP1    
EOF

sbatch -d afterok:$JOBID $JOB2
