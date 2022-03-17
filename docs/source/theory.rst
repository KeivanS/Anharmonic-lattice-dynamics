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
    \Phi_{ij} = \Phi_{ji} \\
    \Psi_{ijk} = \Psi_{ikj} = \Psi_{jik} = \Psi_{kji} = ... \\
    \chi_{ijkl} = \chi_{ikjl}=\chi_{ijlk}=\chi_{jikl} = ...

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
    0 = \sum_{\tau} \, \Pi_{0 \tau}^{ \alpha} \,\,\, \forall (\alpha) \,\, => {\rm Total \,~force \,~on\,~unit\,~cell = 0}  \\
    0 = \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta} \,\,\, \forall (\alpha\beta,\tau) \\
    0 = \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \,\,\, \forall(\alpha\beta\gamma,R_1,\tau_1) \\
    0 = \sum_{R_3,\tau_3} \, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta} \,\,\,  \forall(\alpha\beta\gamma\delta,R_1R_2,\tau_1\tau_2)

so that diagonal terms in these tensors can be defined in terms of their 
off-diagonal elements, for arbitrary cartesian components.

3. Invariance under an arbitrary rotation of the system

Likewise if the system is rotated arbitrarily, the total
energy and forces should not change. This leads to the following
relations [Ref]_:

.. math::
    0 = \sum_{\tau} \,  \Pi_{0 \tau}^{ \alpha}\, \tau^{\beta}\, \epsilon^{\alpha\beta\nu}, \,\,\, \forall ( \nu) \,\, ({\rm Torque \,on\,unit\,cell}=0)

.. math::
    0 =  \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta}\, (R_1\tau_1)^{\gamma}\, \epsilon^{\beta\gamma\nu}
     + \Pi_{0 \tau}^{ \beta }\, \epsilon^{\beta \alpha \nu}  \,\,\, \forall (\alpha\nu,\tau)

.. math::
    0 = \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \, (R_2\tau_2)^{\delta} \, \epsilon^{\gamma\delta\nu} + \Phi_{0\tau,R_1\tau_1}^{\gamma\beta}\, \epsilon^{\gamma\alpha\nu}  + \Phi_{0\tau,R_1\tau_1}{\alpha\gamma}\, \epsilon^{\gamma\beta\nu} \,\,\, \forall(\alpha\nu,R_1,\tau\tau_1)

.. math::
    0 = \sum_{R_3\tau_3}\, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta}\,
    (R_3\tau_3)^{\mu}\, \epsilon^{\delta\mu\nu}
    + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\delta\beta\gamma} \,\epsilon^{\delta\alpha\nu} +
     \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\delta\gamma} \,\epsilon^{\delta\beta\nu}   \\
    + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\delta} \,\epsilon^{\delta\gamma\nu} 
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

.. [Ref] G. Leibfried and W. Ludwig, in Solid State Physics, edited by F. Seitz and D. Turnbull (Academic, New York, 1961), Vol. 12.
