"""Test structures for poscar2focex.x: one per Setyawan-Curtarolo lattice / sub-case (built with
pymatgen from a space group, lattice parameters and Wyckoff positions), as POSCARs in DIR/<name>/.

    python make_structures.py DIR
"""
import numpy as np, warnings, os, sys
warnings.filterwarnings('ignore')
from pymatgen.core import Structure, Lattice
L=Lattice
cases = {
 'Ge_Fd-3m':    (227, L.cubic(5.65), ['Ge'], [[0,0,0]]),
 'NaCl_Fm-3m':  (225, L.cubic(5.64), ['Na','Cl'], [[0,0,0],[.5,.5,.5]]),
 'Fe_Im-3m':    (229, L.cubic(2.87), ['Fe'], [[0,0,0]]),
 'BZO_Pm-3m':   (221, L.cubic(4.19), ['Ba','Zr','O'], [[0,0,0],[.5,.5,.5],[.5,.5,0]]),
 'FeS2_Pa-3':   (205, L.cubic(5.42), ['Fe','S'], [[0,0,0],[.385,.385,.385]]),
 'TiO2_P42mnm': (136, L.tetragonal(4.59,2.96), ['Ti','O'], [[0,0,0],[.305,.305,0]]),
 'BCT1_I4mmm':  (139, L.tetragonal(4.0,3.0), ['Ca','Ge'], [[0,0,0],[0,0,.38]]),
 'BCT2_I4mmm':  (139, L.tetragonal(3.25,4.95), ['In'], [[0,0,0]]),
 'GeSe_Pnma':   (62, L.orthorhombic(10.8,3.85,4.4), ['Ge','Se'], [[.12,.25,.11],[.85,.25,.48]]),
 'ORCF1_Fddd':  (70, L.orthorhombic(5,7,9), ['S'], [[.1,.2,.3]]),
 'ORCF2_Fddd':  (70, L.orthorhombic(6,7,8), ['S'], [[.1,.2,.3]]),
 'ORCI_Immm':   (71, L.orthorhombic(3,4,5), ['Cu','O'], [[0,0,0],[.1,.2,.3]]),
 'ORCC_Cmcm':   (63, L.orthorhombic(3,5,4), ['Zr','O'], [[0,.1,.25],[0,.4,.25]]),
 'wurtzite':    (186, L.hexagonal(3.25,5.21), ['Zn','O'], [[1/3,2/3,0],[1/3,2/3,.382]]),
 'MoTe2_P63mmc':(194, L.hexagonal(3.52,13.97), ['Mo','Te'], [[1/3,2/3,.25],[1/3,2/3,.62]]),
 'CdI2_P-3m1':  (164, L.hexagonal(4.24,6.84), ['Cd','I'], [[0,0,0],[1/3,2/3,.25]]),
 'Bi_RHL1':     (166, L.hexagonal(4.54,11.86), ['Bi'], [[0,0,.234]]),
 'RHL2_R-3m':   (166, L.hexagonal(5.0,3.0), ['Al'], [[0,0,.2]]),
 'MCL_P2m':     (10, L.monoclinic(3,4,5,100), ['K','Cl'], [[.1,.2,.3],[.4,.5,.15]]),
 'MCLC_C2m_a':  (12, L.monoclinic(10,3,4,105), ['Al','O'], [[.1,0,.3],[.4,0,.15]]),
 'MCLC_C2m_b':  (12, L.monoclinic(5,8,6,110), ['Al','O'], [[.1,0,.3],[.4,.2,.15]]),
 'MCLC_C2m_c':  (12, L.monoclinic(3,9,4,95), ['Al','O'], [[.1,0,.3],[.4,.2,.15]]),
 'TRI_a':       (2, L.from_parameters(3,4,5,80,85,95), ['Si','O'], [[.1,.2,.3],[.4,.1,.15]]),
 'TRI_b':       (2, L.from_parameters(4,5,6,100,110,105), ['Si','O'], [[.1,.2,.3],[.4,.1,.15]]),
 'MCLC5_C2m':   (12, L.monoclinic(3,3,3,95), ['Al','O'], [[0,0,0],[.1,0,.3]]),
 # ORCF3 is the boundary 1/a^2 = 1/b^2 + 1/c^2
 'ORCF3_Fmmm':  (69, L.orthorhombic(1/np.sqrt(1/25+1/36),5,6), ['Cu'], [[0,0,0]]),
}
out=sys.argv[1]
for name,(sg,lat,sp,co) in cases.items():
    s=Structure.from_spacegroup(sg,lat,sp,co)
    s=s.get_primitive_structure()
    d=os.path.join(out,name); os.makedirs(d,exist_ok=True)
    s.to(filename=os.path.join(d,'POSCAR'),fmt='poscar')
print(len(cases),'cases')
