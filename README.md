![ALATDYN Logo](https://github.com/KeivanS/Anharmonic-lattice-dynamics/blob/main/docs/source/_static/img/aladyn-logo.png)

## For Users
Anharmonic-LATtice-DYNamics (ALATDYN) is a lattice dynamics code. It calculates thermodynamic and thermal transport properties of solid crystalline materials from data on their force and potential energy as a function of positions, using density functional theory or classical molecular dynamics as input data.

For installation and usage instructions, please read our [documentation](https://aladyn.readthedocs.io/en/latest/index.html). In short:

```bash
git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
cd Anharmonic-lattice-dynamics
make            # FOCEX + utilities, SC-SNAPS, THERMACOND, ANFOMOD -> ~/ALADYN/bin
                # or one part: make focex | make sc_snaps | make thermacond
```

If you have any queries and want to know more about our code, please visit [our forum](https://matsci.org/c/aladyn/57) hosted at the MATSCI Community Discourse.

### The FOCEX agent: from a crystal structure to phonons

[`FOCEX/AGENT/focex-agent.py`](FOCEX/AGENT) runs the whole FOCEX workflow from one short instruction file:

- **from_scratch**: from a POSCAR, it runs the VASP relaxation, SCF and Born charges; the SC-SNAPS snapshots; one VASP run per snapshot; the force collection; and FOCEX. It then plots the dispersion, density of states, Grüneisen parameters, elastic constants and thermal expansion.
- **extract_only**: from existing supercell forces (`POSCAR1`, `FORCEDISP1`), it rotates the data into the FOCEX frame and runs FOCEX and the plots.

Every step is checked and logged, and the workflow can be restarted at any point; it runs on a workstation or on a SLURM cluster. When a run stops, `FOCEX/AGENT/diagnose.py` finds where and why, and writes a script that repairs and resubmits it.

The agent uses [VASP-Flow](https://github.com/KeivanS/VASP-Flow) for the VASP inputs and batch scripts; no prior experience with it is needed. The manual is [`FOCEX/AGENT/focexagent_manual.pdf`](FOCEX/AGENT/focexagent_manual.pdf); the same text is on [Read the Docs](https://aladyn.readthedocs.io/en/latest/runfocexagent.html). It covers installation, the two modes, every keyword, and the diagnosis.

## For Developers

### Code contribution

The developers should follow these conventions when adding their contribution to the project:
- **Documentation**: at the beginning of every fortran subroutine add a comment which explains what the subroutine does and why. The comments at the top of the subroutine should start with `!!`

    As an example:
    ```
    Subroutine svd(n,m,a,b,x)
    !! this subroutine performs a singular value decomposition to solve a set of overdetermined linear equations defined by the matrix a(n,m) and array b(n)
    !! it is used to solve for the force constants given force-displacement data
    ```
  
    In the above example, the first line explains what it does and the second line explains why it is used.

- **Style**: Standard fortran programming style should be adopted: 

    All subroutines must start with `IMPLICIT NONE`

    Intent of every variable must be specified (`IN`, `OUT`, or `INOUT`)

    `END` statements should be followed with the name of the corresponding subroutine

### Local build
To build the docs webpage locally:

```bash
pip install -r docs/requirements.txt
cd docs
make html
```

Open build/html/index.html in your browser.

## Citations

- K. Esfarjani and H.T. Stokes, _Phys Rev B_ **77**, 144112 (2008). [doi:10.1103/PhysRevB.77.144112](https://doi.org/10.1103/PhysRevB.77.144112)
- K. Esfarjani, H. T. Stokes, and G. Chen, _Phys Rev B_ **84**, 085204 (2011). [doi:10.1103/PhysRevB.84.085204](https://doi.org/10.1103/PhysRevB.84.085204)
- K. Esfarjani, H. Stokes, S.N. Sadeghi, Y. Liang, B. Timalsina, H. Meng, J. Shiomi, B. Liao, and R. Sun, _Computer Physics Communications_ **312**, 109575 (2025). [doi:10.1016/j.cpc.2025.109575](https://doi.org/10.1016/j.cpc.2025.109575)
