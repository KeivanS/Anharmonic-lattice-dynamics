Theory
======

The force constants are defined as derivatives of the potential
energy with respect to atomic displacements about their equilibrium
position. We write the potential energy in the following form:

.. math::
    V = V_0 + \sum_i \Pi_i u_i +  \frac{1}{2!} \sum_{ij}\Phi_{ij} \, u_i u_j +  \frac{1}{3!} \sum_{ijk}\Psi_{ijk} \,
    u_i u_j u_k + \frac{1}{4!} \sum_{ijkl}\chi_{ijkl} \, u_i u_j u_k u_l

where the roman index :math:`i` labels the triplet :math:`(R,\tau,\alpha)` with :math:`R`
being a translation vector of the primitive lattice, :math:`\tau`  an atom within the
primitive unit cell, and :math:`\alpha` the cartesian coordinate of the atomic
displacement :math:`u`. In other words,  :math:`u_i` is the displacement of the
atom :math:`(R,\tau)` in the direction :math:`\alpha` from its equilibrium or
any reference position.
:math:`\Phi,\Psi` and :math:`\chi` are respectively the harmonic, cubic and
quartic FCs, whereas :math:`\Pi` is the negative of the residual force on atom :math:`i`,
and is zero if the potential energy :math:`V` is expanded around its minimum
or equilibrium configuration.
In clusters or molecules the formalism is the same, only the
translation vector :math:`R` needs to be dropped.

The resulting force on atom :math:`i` would then be:

.. math::
    F_i = - \frac{\partial V}{\partial u_i} = -\Pi_i -\sum_j  \Phi_{ij} \,u_j - \frac{1}{2} \sum_{jk}\Psi_{ijk} \,
    u_j u_k - \frac{1}{3!} \sum_{jkl} \chi_{ijkl} \,  u_j u_k u_l

General symmetries of the problem
---------------------------------
The symmetries of the FCs are deduced from
rotational and translational invariance of the system, in addition to the
symmetries of the crystal itself. These symmetries are as follows:

1. Invariance under the permutation of indices:

.. math::
    \Phi_{ij} &= \Phi_{ji} \\
    \Psi_{ijk} &= \Psi_{ikj} = \Psi_{jik} = \Psi_{kji} = ... \\
    \chi_{ijkl} &= \chi_{ikjl}=\chi_{ijlk}=\chi_{jikl} = ...

where :math:`i,j,k` and :math:`l`
refer to neighboring atoms.
This comes about because the force constants are
derivatives of the potential energy, and one can change the order of
differentiation for such analytic function.

2. Invariance under an arbitrary translation of the system

Translational invariance of the whole system (even if not an ordered crystal)
also implies :math:`V(u) = V(u+c)` and :math:`F(u) = F(u+c)` (which is easier to use)
where :math:`u(t)` are the dynamical variables, and :math:`c` is a constant arbitrary
shift vector (may be incommensurate with :math:`R`).

This implies the well-known "acoustic sum rule" (ASR) generalized to higher order FCs:

.. math::
    0 &= \sum_{\tau} \, \Pi_{0 \tau}^{ \alpha} \,\,\, \forall (\alpha) \,\, => {\rm Total \,~force \,~on\,~unit\,~cell = 0}  \\
    0 &= \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta} \,\,\, \forall (\alpha\beta,\tau) \\
    0 &= \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \,\,\, \forall(\alpha\beta\gamma,R_1,\tau_1) \\
    0 &= \sum_{R_3,\tau_3} \, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta} \,\,\,  \forall(\alpha\beta\gamma\delta,R_1R_2,\tau_1\tau_2)

so that diagonal terms in these tensors can be defined in terms of their
off-diagonal elements, for arbitrary cartesian components.

3. Invariance under an arbitrary rotation of the system

Likewise if the system is rotated arbitrarily, the total
energy and forces should not change. This leads to the following
relations [Ref]_:

.. math::
    0 &= \sum_{\tau} \,  \Pi_{0 \tau}^{ \alpha}\, \tau^{\beta}\, \epsilon^{\alpha\beta\nu}, \,\,\, \forall ( \nu) \,\, ({\rm Torque \,on\,unit\,cell}=0) \\
    0 &=  \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta}\, (R_1\tau_1)^{\gamma}\, \epsilon^{\beta\gamma\nu}
     + \Pi_{0 \tau}^{ \beta }\, \epsilon^{\beta \alpha \nu}  \,\,\, \forall (\alpha\nu,\tau) \\
    0 &= \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \, (R_2\tau_2)^{\delta} \, \epsilon^{\gamma\delta\nu} + \Phi_{0\tau,R_1\tau_1}^{\gamma\beta}\, \epsilon^{\gamma\alpha\nu}  + \Phi_{0\tau,R_1\tau_1}{\alpha\gamma}\, \epsilon^{\gamma\beta\nu} \,\,\, \forall(\alpha\nu,R_1,\tau\tau_1) \\
    0 &= \sum_{R_3\tau_3}\, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta}\,
    (R_3\tau_3)^{\mu}\, \epsilon^{\delta\mu\nu}
    + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\delta\beta\gamma} \,\epsilon^{\delta\alpha\nu} +
     \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\delta\gamma} \,\epsilon^{\delta\beta\nu}   \\
    & \quad + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\delta} \,\epsilon^{\delta\gamma\nu}
    \,\, \, \forall(\alpha\nu,R_1R_2,\tau\tau_1\tau_2)

Here :math:`\epsilon^{\alpha\beta\gamma}` is the anti-symmetric Levy-Civita symbol, and
:math:`(R\tau)^{\alpha}` refers to the cartesian component :math:`\alpha` of the reference position
vector of the atom :math:`\tau` in unit cell defined by :math:`R`. Moreover, an implicit
summation over repeated cartesian indices is implied.

As we see, rotational invariance relates the second to the
third order terms, and the third to the fourth order terms, implying that
if the expansion of the potential energy is truncated after the fourth order
terms, we need to start with the fourth order terms, and application of
rotational invariance rules will give us constraints on third and
second order FCs respectively.

Point and/or space group symmetries
-----------------------------------

The other constraints come from symmetry operations, such as lattice
translation, rotation, mirror
or any symmetry operation of the space and/or point group of the crystal and/or molecule
which leaves the latter invariant.
Invariance under a translation of the system by any translation
lattice vector :math:`R` implies the following relations:

.. math::
    \Pi_{R\tau}^{\alpha} &= \Pi_{0 \tau}^{ \alpha} \,\, \forall (R \tau \alpha)  \\
    \Phi_{R\tau,R_1\tau_1}^{\alpha\beta} &= \Phi_{0 \tau ,R_1-R \tau_1}^{\alpha\beta}  \\
    \Psi_{R\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} &= \Psi_{0 \tau ,R_1-R \tau_1,R_2-R \tau_2}^{\alpha\beta\gamma}  \\
    \chi_{R\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta} &= \chi_{0 \tau ,R_1-R \tau_1,R_2-R \tau_2, R_3-R \tau_3}^{\alpha\beta\gamma\delta}

So in an ideal crystal, this reduces the number of force constants considerably
(by the number of unit cells, to be exact), meaning that we will use for the atoms in all the other cells the same FCs
as those we have kept for the atoms in the "central" cell.

If a rotation or mirror symmetry operation is denoted by :math:`S`,
we must have:

.. math::
    \Pi_{S\tau}^{\alpha} &= \sum_{\alpha'} \,\Pi_{\tau}^{\alpha'} \, {\cal S}_{\alpha,\alpha'} \\
    \Phi_{ S\tau, S\tau_1}^{\alpha\beta} &= \sum_{\alpha'\beta'} \,\Phi_{\tau,\tau_1}^{\alpha'
    \beta'} \, {\cal S}_{\alpha,\alpha'}\, {\cal S}_{\beta,\beta' }  \\
    \Psi_{ S\tau, S\tau_1, S\tau_2}^{\alpha\beta\gamma} &= \sum_{\alpha'\beta'\gamma'} \,\Psi_{\tau,\tau_1\tau_2}^{\alpha' \beta'\gamma'} \, {\cal S}_{\alpha,\alpha'}\, {\cal S}_{\beta,\beta' } \, {\cal S}_{\gamma,\gamma' }   \\
    \chi_{ S\tau, S\tau_1, S\tau_2, S\tau_3}^{\alpha\beta\gamma\delta} &= \sum_{\alpha'\beta'\gamma'\delta'} \,\chi_{\tau,\tau_1,\tau_2,\tau_3}^{\alpha' \beta' \gamma' \delta'} \, {\cal S}_{\alpha,\alpha'}\, {\cal S}_{\beta,\beta' } \, {\cal S}_{\gamma,\gamma'}\, {\cal S}_{\delta,\delta' }

where :math:`{\cal S}_{\alpha,\alpha'}` are the :math:`3\times3` matrix elements of the
symmetry operation :math:`S`.
These symmetry relations impose a linear set of constraints on the force
constants.

Any physically correct model of force constants
must satisfy the invariance relations. On the other hand, we do
approximations by truncating the range of FCs and their order in the
Taylor expansion of the potential. Therefore imposing the constraints
will move their value away from the true value, but has the advantage
that they are physically correct, and will for instance reproduce the linear
dispersion of acoustic phonons near :math:`k=0`. So one should keep in mind
that an unrealistic truncation to a too short of a range will produce
results in disagreement with experiments.

.. [Ref] G. Leibfried and W. Ludwig, in Solid State Physics, edited by F. Seitz and D. Turnbull (Academic, New York, 1961), Vol. 12.
