#!/bin/bash
export OMP_NUM_THREADS=${VASP_OMP_NUM_THREADS:-1}

cd "/Users/ke4c/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/AGENT/Si_Test"
"/Users/ke4c/myenv/bin/python3" "/Users/ke4c/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/AGENT/focex-agent.py" run "/Users/ke4c/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/AGENT/Si_Test"
