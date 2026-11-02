"""Re-apply consolidation steps 2-5 to THERMACOND (run in the THERMACOND folder of a main checkout):
constants, params, geometry, ios from ../common; THERMACOND-only variables in tc_params, tc_ios."""
import re, glob
def module_cut(t, start, end, repl):
    s = t.index(start); e = t.index(end) + len(end); return t[:s] + repl + t[e:]
files = [f for f in sorted(glob.glob('*.f90')) if f not in ('kappa9.f90', 'mpi_stub.f90')]
# --- step 2: constants
p = 'modules.f90'; t = open(p).read()
t = module_cut(t, ' module constants', ' end module constants', '! module constants: ../common/constants.f90 (shared with FOCEX)')
# --- step 3: params
t = module_cut(t, 'module params', 'end module params', """module tc_params
!! THERMACOND inputs not in the shared params module (../common/params.f90)
 real(8) lmicron
 integer ntemp,calc_cross,nsmax
 real(8) rcut(4),tau0,wshift(3)
 real(8) q_cross(3)
 real(8) alphamix
end module tc_params""")
t = t.replace("""       if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then""",
              """       ! (not called by THERMACOND; nshells is the shared 2-D array nshells(rank,atom))
       if ( iatomneighbor(i0,j) .eq. nshells(2,i0) .and. dist(j) .gt. rmax ) then""")
# --- step 4: geometry
t = module_cut(t, ' module geometry', ' end module geometry', '! module geometry: ../common/geometry.f90 (shared with FOCEX)')
# --- step 5: ios
t = module_cut(t, 'module ios', 'end module ios', """module tc_ios
!! THERMACOND unit numbers not in the shared ios module (../common/ios_mod.f90); udos, uibz and
!! ugrun keep the THERMACOND values (the shared module has others)
 integer, parameter::ubs=51,udos=41,ueiv=61,ukap=71,uslf=3000, &
&         ufc0=28,ualph=27,uvel=31,umod=2000,  &
&         fbz=93,bzf=95,uibz=94,ibs=96,ucors=65,uv3=85,ugrun=81,urate=88, &
&         ucof=87,ucross=89, k_kappa=9000
 integer, parameter:: debug=9999, self_detail=9990, iself_bs=9980, ksubset_bs_inp=9970    ! sy.
end module tc_ios""")
open(p, 'w').write(t)
P = 'use params, only : tolerance,margin,nconfigs,classical,fdfiles,threemtrx,lamin,lamax,nshells,include_fc,tmin,tmax,verbose'
I = 'use ios, only : ulog,uposcar,uparams,utraj,ufco,umap,umatrx,utimes,ufc1,ufc2,ufc3,ufc4,ufc,ufit1,ufit2,ufit3,ufit4,uborn,write_out'
for f in files:
    L = open(f).read().split('\n'); out = []
    for l in L:
        m = re.match(r'^(\s*)use\s+params\s*,\s*only\s*:\s*nsmax\s*$', l, re.I)
        if m: out.append(m.group(1) + 'use tc_params, only : nsmax'); continue
        m = re.match(r'^(\s*)use\s+params\s*(!.*)?$', l, re.I)
        if m: out += [m.group(1) + P, m.group(1) + 'use tc_params']; continue
        m = re.match(r'^(\s*)use\s+ios\s*(!.*)?$', l, re.I)
        if m: out += [m.group(1) + I, m.group(1) + 'use tc_ios']; continue
        out.append(l)
    open(f, 'w').write('\n'.join(out))
t = open('rest.f90').read()
t = t.replace(" write(ulog,*)'include_fc  =',include_fc\n", " write(ulog,*)'include_fc  =',include_fc(1:4)   ! THERMACOND uses ranks 1-4\n")
open('rest.f90', 'w').write(t)
# the Makefile (COMMON list and dependency chain) is kept from the worktree
