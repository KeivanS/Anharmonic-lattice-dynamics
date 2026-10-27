Theory
======

The force constants are defined as derivatives of the potential
energy with respect to atomic displacements about their equilibrium
position. We write the potential energy in the following form:

.. math::
    :label: actualPotential

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

This implies the well-known "acoustic sum rule" (ASR) generalized to higher order FCs
[Leibfried-Ludwig]_:

.. math::
    0 &= \sum_{\tau} \, \Pi_{0 \tau}^{ \alpha} \,\,\, \forall (\alpha) \,\, => {\rm Total \,~force \,~on\,~unit\,~cell = 0}  \\
    0 &= \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta} \,\,\, \forall (\alpha\beta,\tau) \\
    0 &= \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \,\,\, \forall(\alpha\beta\gamma,\tau,R_1\tau_1) \\
    0 &= \sum_{R_3,\tau_3} \, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta} \,\,\,  \forall(\alpha\beta\gamma\delta,\tau,R_1\tau_1,R_2\tau_2)

so that diagonal terms in these tensors can be defined in terms of their
off-diagonal elements, for arbitrary cartesian components.

3. Invariance under an arbitrary rotation of the system

Likewise if the system is rotated arbitrarily, the total
energy and forces should not change. This leads to the following
relations [Leibfried-Ludwig]_:

.. math::
    0 &= \sum_{\tau} \,  \Pi_{0 \tau}^{ \alpha}\, \tau^{\beta}\, \epsilon^{\alpha\beta\nu}, \,\,\, \forall ( \nu) \,\, ({\rm Torque \,on\,unit\,cell}=0) \\
    0 &=  \sum_{R_1,\tau_1} \, \Phi_{0\tau,R_1 \tau_1}^{\alpha\beta}\, (R_1\tau_1)^{\gamma}\, \epsilon^{\beta\gamma\nu}
     + \Pi_{0 \tau}^{ \beta }\, \epsilon^{\beta \alpha \nu}  \,\,\, \forall (\alpha\nu,\tau) \\
    0 &= \sum_{R_2,\tau_2} \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\gamma} \, (R_2\tau_2)^{\delta} \, \epsilon^{\gamma\delta\nu} + \Phi_{0\tau,R_1\tau_1}^{\gamma\beta}\, \epsilon^{\gamma\alpha\nu}  + \Phi_{0\tau,R_1\tau_1}^{\alpha\gamma}\, \epsilon^{\gamma\beta\nu} \,\,\, \forall(\alpha\beta\nu,\tau,R_1\tau_1) \\
    0 &= \sum_{R_3\tau_3}\, \chi_{0\tau,R_1\tau_1,R_2\tau_2,R_3\tau_3}^{\alpha\beta\gamma\delta}\,
    (R_3\tau_3)^{\mu}\, \epsilon^{\delta\mu\nu}
    + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\delta\beta\gamma} \,\epsilon^{\delta\alpha\nu} +
     \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\delta\gamma} \,\epsilon^{\delta\beta\nu}   \\
    & \quad + \, \Psi_{0\tau,R_1\tau_1,R_2\tau_2}^{\alpha\beta\delta} \,\epsilon^{\delta\gamma\nu}
    \,\, \, \forall(\alpha\beta\gamma\nu,\tau,R_1\tau_1,R_2\tau_2)

Here :math:`\epsilon^{\alpha\beta\gamma}` is the anti-symmetric Levi-Civita symbol, and
:math:`(R\tau)^{\alpha}` refers to the cartesian component :math:`\alpha` of the reference position
vector of the atom :math:`\tau` in unit cell defined by :math:`R`. Moreover, an implicit
summation over repeated cartesian indices is implied.

As we see, rotational invariance relates the first to the second order terms, the second to the
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

Variational approach using Bogoliubov's inequality
--------------------------------------------------
The self-consistent phonon theory implemented in SCOP8, and the equations of this section, are
described in [Esfarjani-Liang2022]_. It is based on Bogoliubov's inequality

.. math::
    F\le F_{trial}= F_0 +\langle V-V_0 \rangle_0

The sharp brackets in above equation mean thermal average with respect to the density matrix :math:`\rho_0` associated with the variational or trial harmonic Hamiltonian :math:`H_0`

.. math::
    \langle A \rangle_0=Tr \rho_0 A; \,\, \rho_0=e^{-H_0/kT} /Z ; \,\, Z={\rm Tr}\, e^{-H_0/kT} ;\, F_0=-kT {\rm ln} Z 

As can be seen in above equations, :math:`\rho_0` is the trial density matrix and Z is the corresponding partition function. :math:`F_{trial}` is the effective or trial free energy that is optimized during iterations. :math:`F_0` and :math:`V_0` are the trial harmonic free energy and potential energy. :math:`V` is the (approximated) real or actual potential energy in a Taylor expansion form to describe the real crystal with force constants :math:`\Phi`,:math:`\Psi` and :math:`X` and selected variational parameters (strain, internal coordinates).

In order to include structural phase transitions and thermal expansion, we introduce two additional sets of variational parameters: the strain tensor :math:`\eta` and :math:`u_i^0`, the internal relaxation of atoms beyond the strain effect (captured by :math:`\eta R_i` due to changes in the unit cell shape or volume, and so that the general displacement of atom :math:`i` from its reference position as a result of cell deformation can be written as: :math:`u_i(t) = \eta R_i + u_i^0 + (1+\eta) y_i(t)=  S_i+ y_i(t)` where :math:`y_i` means the dynamic displacement of atom i, so by definition :math:` \langle y_i \rangle_0 =0`. :math:`R_i` refers to the lattice translation vector of cell containing :math:`i`. :math:`S_i=  \eta R_i + u_i^0` is the static displacement of atom :math:`i` introduced for brevity of notations.
The trial harmonic potential contains the variational parameters :math:`K_{ij}` and is defined as:

.. math::
    V_0 =\sum_{ij} \,\frac{1}{2!} K_{ij}\, y_i\, y_j 

For the trial harmonic Hamiltonian, the free energy can be calculated analytically:

.. math::
    F_0= kT \sum_{k,\lambda} {\rm ln} \left[ 2 {\rm sinh} \frac{\beta \hbar \omega_{k\lambda} }{2}\right]

where :math:`\beta=1/k_BT` and :math:` \omega_{\lambda} ` is the harmonic frequency of mode :math:`\lambda` obtained from diagonalizing the dynamical matrix associated with the trial force constants:  :math:`D(k) = \sum_R \frac{K_{ij}}{\sqrt {m_i m_j}}\, \epsilon^{i k \cdot (R_i-R_j)}`. :math:`k` refers to the k vector on a selected k mesh in reciprocal space. The matrix :math:`D` being Hermitian, it has real eigenvalues denoted by :math:`\omega_{\lambda}^2` and eigenvectors :math:`epsilon_{i\alpha,\lambda}` where :math:`\alpha` is the cartesian coordinate and :math:`\lambda` refers to the vibrational mode (:math:`i\alpha` can be understood as line index and :math:`\lambda` a the column index of the unitary eigenvector matrix :math:`e`). We also need the thermal averages of the trial and anharmonic potentials :math:`V_0` and :math:`V`. In terms of the eigenvalues  and eigenvectors of the above dynamical matrix, we have:

.. math::
    :label: yy

     \langle y_i^\alpha y_j^\beta \rangle= \frac{\hbar}{2} \sum_{k,\lambda}  \frac{ (2 n_{k\lambda}+1) }{\omega_{k\lambda}} \,\frac{\epsilon_{i\alpha,k\lambda} \, \epsilon_{j\beta,k\lambda}^\dagger}{\sqrt{m_i m_j}}e^{i\vec{k}\cdot(\vec{R}_j-\vec{R}_i)}

Thus the thermal averaged harmonic or trial potential :math:`\langle V_0\rangle` can be simply written into :math:`\langle V_0 \rangle = \sum_{ij} \,\frac{1}{2!} K_{ij}\, \langle y_i y_j \rangle` . Similarly, the thermal averaged actual potential :math:`\langle V \rangle` can be presented in Eq.[\ref{Vaverage}] as:

.. math:: 
    :label: Vaverage

    \langle V \rangle &= \frac{1}{2}\sum_{ij} \,  \Phi_{ij}\, \Big(( S_i S_j+ \langle y_i\, y_j \rangle \Big) 
    + \frac{1}{6}\sum_{ijk}  \Psi_{ijk}\,\Big(S_i S_j + 3 \langle y_i \,y_j \rangle \Big) \,S_k \nonumber\\
    &+\frac{1}{24} \sum_{ijkl}  X_{ijkl}\, \Big( S_i S_j S_k S_l + 6 S_i S_j \langle y_k \,y_l \rangle +  \langle y_i y_j y_k  y_l \rangle \Big)

For the sake of simplicity, all the cartesian related notations have been  symmetrized and omitted. In fact, terms like :math:`\eta` are rank 2 tensors in cartesian coordinates which means all the matrix multiplications and gradients calculations have to been handled carefully in real code implementation.

Choice of variational parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
From here we can use variational approach in an iterative manner to
minimize the trial Free energy :math:`F_{trial}` by solving the
gradients of it with respect to different variational parameters. The
Broyden method [Broyden]_ is being used in the code for this non-linear root
finding process. As mentioned previously, the variational parameters in
this formalism are internal relaxation of atoms :math:`u_i^0`, strain
tensor :math:`\eta` and trial harmonic force constants :math:`K_{ij}`.
Thus we have the following gradients formulas resulting from the
minimization condition:

.. math::
   :label: dfdu

   \frac{\partial F_{trial} }{\partial u_i^0 } =\frac{\partial  \langle V \rangle }{\partial u_i^0 } =0

.. math::
   :label: dfdeta

   \frac{\partial F_{trial} }{\partial \eta } =\frac{\partial \langle V \rangle }{\partial \eta }=0

.. math::
   :label: dfdk

   \frac{\partial F_{trial} }{\partial K_{ij} }  =0 


Same as before, cartesian related labelling and notations are overly
symmetrized and omitted for better readability. Eq. :eq:`dfdu` 
and Eq. :eq:`dfdeta` are straightforward and self-explanatory
equilibrium conditions; these are implemented in codes following the
exact representations which are easy to understand and to maintain.
However, these exact representations are omitted here as they could be
easily formulated from Eq. :eq:`Vaverage`. What worth
mentioning is the third gradients Eq. :eq:`dfdk` since trial
force constant :math:`K_{ij}`, although being a variational parameter,
is implicit in the expression :math:`F_{trial}` because of the term
:math:`\langle y_i y_j\rangle`. Refer to Eq. :eq:`yy`, using the
identity
:math:`2 n_{\lambda}+1= {\rm coth} \frac{ \beta \hbar \omega_{\lambda} }{2}`,
one can see that :math:`\langle y_i y_j\rangle` is a function of
eigenvalues :math:`\omega_\lambda` and eigenvectors
:math:`e_{i\alpha,\lambda}`, they both originated from diagonalization
of dynamical matrix D(k) while D(k) itself depends solely on the trial
force constant :math:`K_{ij}`. So it’s not easy to write out
Eq. :eq:`dfdk` in an explicit way that is necessary for code
implementation. However, by incorporating chain rule and
the free energy equation, we find a one to one equivalence
between Eq. :eq:`dfdk` and the following equation:

.. math:: 
    :label: dfdk2

    \frac{\partial \langle V \rangle }{\partial \langle y_iy_j\rangle } - \frac{1}{2} K_{ij} =0 \label{dfdk2}

Apparently, Eq. :eq:`Vaverage` shows that
:math:`\langle V\rangle` is an explicit polynomial of
:math:`\langle yy\rangle`, so this new equation can be easily translated
into codes in a straightforward way. One should be aware that,
Eq. :eq:`dfdk2` is by no means a ’gradient’ formula of
anything. It’s just a one to one correspondence to
Eq. :eq:`dfdk`, in other words, if Eq. :eq:`dfdk` is
satisfied (after fully optimizing the free energy), then
Eq. :eq:`dfdk2` will simultaneously be satisfied, and vice
versa. There are no further physical significance behind
Eq. :eq:`dfdk2`. In fact, the trial force constants
:math:`K_{ij}` and the thermally averaged dynamic displacement
:math:`\langle y_i y_j \rangle` can be considered interchangeable from a
mathematical standpoint. This provides two options when selecting the
variational parameters for self-consistent computations. The advantage
of using :math:`\langle y_i y_j \rangle` becomes more evident as it
directly appears as explicit parameters in the potential energy
expression in Eq. :eq:`Vaverage`, unlike the trial force
constants :math:`K_{ij}`, which do not have a direct manifestation in
energy or gradient formulas.

Some may argue that a significant benefit of favoring
:math:`\langle y_i y_j \rangle` over :math:`K_{ij}` is that it
eliminates the need for diagonalizing the dynamic matrix from
:math:`K_{ij}` to obtain new sets of phonon frequencies, thus bypassing
the entire computation in reciprocal space, which by employing the
Broyden method, the new value of :math:`\langle y_i y_j \rangle` can be
obtained straightforwardly. However, there are three main reasons why it
is preferable to still use trial force constants :math:`K_{ij}` as
variational parameters instead of :math:`\langle y_i y_j \rangle`.
First, as discussed earlier, there is a possibility of divergent
:math:`\langle y_i y_j \rangle` values for certain k points, especially
the gamma point. This could make it challenging for the quasi-Newton
iterative method to reconcile these outlier
:math:`\langle y_i y_j \rangle` values with the remaining normal
:math:`\langle y_i y_j \rangle` values, leading to convergence
difficulties. On the other hand, by choosing :math:`K_{ij}` as the
variational parameter, it becomes possible to regulate these abnormal
:math:`\langle y_i y_j \rangle` values using the aforementioned
analytical approximation, resulting in a smoother free energy landscape
that facilitates the Broyden method in finding the global minimum.
Second, using :math:`\langle y_i y_j \rangle` as the variational
parameter means bypassing the examination of phonon frequencies within
iterations. This eliminates the ability to assess the stability of
phonons at each iteration, potentially leading to the computation
becoming trapped in an unstable structure. Third, from a mathematical
perspective, finding an appropriate way to update the
:math:`K_{ij}`\ term in Eq. :eq:`dfdk2` without falling into
circular reasoning poses a challenge. In fact, it is analytically more
feasible to calculate :math:`\langle y_i y_j \rangle` from trial force
constants :math:`K_{ij}` than the other way around.

Preserve acoustic sum rule
^^^^^^^^^^^^^^^^^^^^^^^^^^
The acoustic sum rule holds significant importance in this computational
approach, and it can pose considerable challenges when attempting to
enforce it using real space Taylor series methods. Nevertheless, during
the self-consistent calculation, it is unnecessary to enforce it if the
initial input data already adheres to the acoustic sum rule (for the
sake of simplicity, it will be referred as ASR henceforth), which can be
found in [Leibfried-Ludwig]_ since
it originates from translational invariance. The trial force constants
:math:`K_{ij}` can be expressed as [Esfarjani-Liang2022]_:

.. math::
    :label: trialFC

    K_{i,j}^{\alpha,\beta} \approx\;\; \Phi_{i,j}^{\alpha,\beta} + \sum_{k}\Psi_{i,j,k}^{\alpha,\beta,\gamma}(S_k^\gamma- S_i^\gamma)+
       \frac{1}{2}\sum_{k,l}X_{i,j,k,l}^{\alpha,\beta,\gamma,\delta} [(S_k-S_i)^\gamma(S_l-S_i)^\delta + \langle y_k^\gamma y_l^\delta \rangle]
       \label{eq:trialFC}

which needs to preserve the ASR across the whole iterative procedures of
Broyden method. Thus, each term in Eq. :eq:`trialFC`
needs to satisfy ASR independently. The first term is trivial as long as
the input real force constant :math:`\Phi_{i,j}^{\alpha,\beta}`
satisfies ASR. For the second term, we can switch the position of
indexes into
:math:`\sum_{j}\,\sum_{k}\Psi_{i,j,k}^{\alpha,\beta,\gamma}S_k^{\gamma} = \sum_{k}S_{k}^\gamma\,\sum_{j}\Psi_{i,j,k}^{\alpha,\beta,\gamma}`,
which means the only requirement is the same: input cubic force
constants :math:`\Psi` has to satisfy ASR. A similar line of reasoning
can be applied to the final term in Eq. :eq:`trialFC`
thus establishing that the trial force constant :math:`K_{i,j}` despite
its lack of symmetry and unrestricted variability, automatically
satisfies the acoustic sum rule (ASR) during the self-consistent
iterations. This raises the question of how we can ensure or enforce the
ASR in the input force constants data, namely :math:`\Phi, \Psi, X`.
Several subroutines are specifically implemented for this task, based on
the idea that is to assign specific values to those force constants that
involve identical atomic indexes. The rules that can also be found in
the source code file ’force_update.f90’, are written as below:

.. math::

    &\Phi_{i,i} = -\sum_{j,j\neq i}\Phi_{i,j},\; \forall i \\
    \\
    &\text{  } \begin{cases}
    \Psi_{i,j,j} = -\sum_{k,i\neq j, k\neq j}\Psi_{i,j,k}, \;\forall i,j \nonumber \\
    \Psi_{i,i,k} = -\sum_{j,i\neq j, k\neq i}\Psi_{i,j,k}, \; \forall i,k \nonumber \\
    \Psi_{i,i,i} = -\sum_{k,k\neq i}\Psi_{i,i,k}, \;\forall i \nonumber \\
    \end{cases}\\
    \\
    &\text{  } \begin{cases}
    X_{i,j,k,k} = -\sum_{l,k\neq j, k\neq i, j\neq i,l\neq k}X_{i,j,k,l}, \; \forall i,j,k \nonumber \\
    X_{i,j,j,l} = -\sum_{k,l\neq j, l\neq i, j\neq i,k\neq j}X_{i,j,k,l}, \; \forall i,j,l \nonumber \\
    X_{i,j,j,j} = -\sum_{l,j\neq i, l\neq j }X_{i,j,j,l}, \; \forall i,j \nonumber \\
    X_{i,i,k,l} = -\sum_{j,k\neq i, l\neq i, j\neq i}X_{i,j,k,l}, \; \forall i,k,l \nonumber \\
    X_{i,i,i,i} = -\sum_{l,l\neq i}X_{i,i,i,l}, \; \forall i 
    \end{cases}

With this in mind, a translational form modification should also be
exerted on the potential energy formula
Eq. :eq:`Vaverage` and the free energy
gradient formulas Eq. :eq:`dfdu`, Eq. :eq:`dfdeta`
and Eq. :eq:`dfdk`. This can be done by replacing the every
occurrence of statistic displacement product :math:`S_iS_j` with
:math:`-\frac{1}{2}(S_j-S_i)^2` and etc, similar to
Eq. :eq:`trialFC` which is already in this
translational invariant representation. This will inevitably lead to a
more complicated expression especially for the gradient formulas, but it
works the best with the input force constants data generated by
**FOCEX** and has the minimum numerical deficiencies. The explicit
formulas are skipped here for the purpose of this paper, but can be
found in the source code file ’VA_math.f95’.


Divergence in dynamic displacement correlations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the numerical computation of dynamic displacements correlation
:math:`\langle y_i^\alpha y_j^\beta\rangle`, certain divergences will
happen, perforce, due to the gamma point (k=0) term in the sum in
formula Eq. :eq:`yy`. The reason behind this is that
:math:`\omega_{k=0,\lambda=acoustic}` should be 0 determined by the ASR.
Also notice that the Bose-Einstein coefficient
:math:`n_\lambda = (e^{\beta\hbar\omega_{k\lambda)}}-1)^-1` will also
diverge for small phonon frequency at high temperature limit. To address
this matter in a more suitable manner, we attempted to resolve it by
excluding the gamma point term and substituting it with a spherical
integral. This integral can be analytically derived as a function of
temperature. It is well-established that at the gamma point, only the
three acoustic phonon bands, including one longitudinal and two
transverse modes, tend to approach zero. The 3 eigenvectors
:math:`\epsilon_\lambda^{\tau\alpha}` at gamma point should satisfy
:math:`\sum_{\tau\alpha}e_{\lambda}^{\tau\alpha}e_{\lambda'}^{\tau\alpha}=\delta_{\lambda\lambda'}`
where :math:`\lambda` can be (1,2,3) or :math:`(//,\perp_1,\perp_2)`
that corresponds to 1 longitudinal and 2 transverse modes at gamma
point. The sets of eigenvectors
:math:`\{\epsilon_{//},\epsilon_{\perp1},\epsilon_{\perp2}\}` are
arbitrarily based on this rule. For this code package, we choose the
following set:

.. math::

   \begin{cases}
           \epsilon_{//}=(sin\theta cos\phi,sin\theta sin\phi,cos\theta)\\
           \epsilon_{\perp1}=(cos\theta cos\phi,cos\theta sin\phi, -sin\theta)\\
           \epsilon_{\perp2}=(-sin\phi,cos\phi)
       \end{cases}

Taking that into consideration, now Eq. :eq:`yy` can be written
into:

.. math:: <y_{0\tau_i}^\alpha y_{R\tau_j}^\beta> =  \frac{1}{N_k\sqrt{m_{\tau_i}m_{\tau_j}}} \sum_{\lambda,k\ne0}[\frac{\hbar}{2\omega_{k\lambda}}(2n_{k\lambda}+1)\epsilon_{k\lambda}^{\dagger\alpha}(\tau_i)\epsilon_{k\lambda}^{\beta}(\tau_j)cos(\vec{k}\cdot\vec{R}_j)]+\text{$\gamma$\;point fix}

where the ’\ :math:`\gamma` point fix’ part can be approximated properly
in an isotropic spherical integral shown below (call it
:math:`I_{0\tau_i,R\tau_j}^{\alpha\beta}`):

.. math::

   \begin{aligned}
        I_{0\tau_i,R\tau_j}^{\alpha\beta} &\approx \frac{1}{N_0N_k\sqrt{m_{\tau_i}m_{\tau_j}}}\sum_{Acoustic}\iiint\frac{\hbar}{2\omega_{k\rightarrow 0,\lambda}}(2n_{k\rightarrow 0\lambda}+1)\epsilon_\lambda^\alpha \epsilon_\lambda^\beta e^{i\vec{k}\cdot(\vec{R}_2-\vec{R}_1)}_{(k\rightarrow 0)}dv+...\nonumber
   \end{aligned}

here :math:`N_0` means the number of atom in a unit cell, :math:`N_k` is
the size of k-mesh. It can be proved that for this matrix
:math:`I_{0\tau_i,R\tau_j}`, all the off-diagonal terms are 0 only the
xx, yy and zz component remain non-zero. A series of estimation can be
utilized given that small k point condition and high temperature limit:

-  Term
   :math:`\frac{\hbar}{2\omega_{k\lambda}}(2n_{k\lambda}+1)=\frac{\hbar}{2\omega_{k\lambda}}coth\frac{\beta\hbar\omega_{k\lambda}}{2}\approx\frac{k_B T}{\omega_{k\lambda}^2}=\frac{k_BT}{C_\lambda^2 k^2}`
   given small :math:`\hbar\omega/k_BT`

-  This spherical integral is over Brillouin zone(BZ) substitute
   :math:`I\simeq \frac{1}{\Delta^3}\int_0^{K_\Delta}d^3k\frac{...}{k^2}=\frac{4\pi K_\Delta}{\Delta^3}`

-  Notice that one BZ mesh volume
   :math:`\Delta^3=\frac{4\pi K_{\Delta}^3}{3}=\frac{\Omega_g}{N_k}` so
   above result is just :math:`\frac{3}{K_\Delta^2}`

-  The exponential term can be simplified to 1 given the fact that k is
   approaching 0. Notice that an analytic result is still easy to
   formulate even without this simplification.

As a result, we end up with the final gamma point approximation formula
stated below:

.. math::
    :label: gammaFix

    \begin{aligned}
       \text{xx or yy direction   }I_{0\tau_i,R\tau_j}^{xx\text{ or }yy} &=\frac{1}{N_k\sqrt{m_{\tau_i}m_{\tau_j}}}\frac{3k_BT}{N_0K_{\Delta}^2}(\frac{1}{3C_{//}^2}+\frac{1}{6C_{\perp1}^2}+\frac{1}{2C_{\perp2}^2}) \nonumber \\
       \text{zz direction    }I_{0\tau_i,R\tau_j}^{zz} &=\frac{1}{N_k\sqrt{m_{\tau_i}m_{\tau_j}}}\frac{3k_BT}{N_0K_{\Delta}^2}(\frac{1}{3C_{//}^2}+\frac{2}{3C_{\perp1}^2})\label{eq:gammaFix}
    \end{aligned}

where :math:`C_\lambda` is the group velocity that could be inferred
from
:math:`C^2=\frac{1}{2}\frac{d\omega^2}{d^2k}\approx\frac{\omega^2_{k=0}+\omega^2_{k=2K}-2\omega^2_{k=K}}{2K_{int}^2}`.
Notice that the ’k’ in :math:`C^2` should be calculate by interval
:math:`K_{int}`\ (euclidean distance between neighbor k points on this
k-mesh). This is different from the :math:`K_\Delta` calculate by
cubic-root the BZ volume. The latest modification based on
Eq. :eq:`gammaFix` has been employed in the source
code file ’VA_math.f95’ and the convergence has been verified by
different k-mesh size choices.

Elastic constants
-----------------


It is possible to extract the elastic constants from the knowledge of the force constants [Wallace1998]_. Below, we outline a simpler method and provide the final formulas [Esfarjani2025]_.

Under a uniform applied strain, denoted by :math:`\eta_{\alpha\beta} = \partial u_{\alpha} / \partial r_{\beta}`, a point :math:`x` of the medium is moved to :math:`x'_{\alpha} = x_{\alpha} + \eta_{\alpha\beta} x_{\beta}`, and the total energy of the harmonic crystal is increased, by definition of the elastic constants, by:

.. math::
   \frac{\Delta E}{V_0} = \frac{C_{\alpha\beta,\gamma\delta} \epsilon_{\alpha\beta} \epsilon_{\gamma\delta}}{2},

where :math:`V_0` is the unit cell volume and the Cauchy strain tensor :math:`\epsilon` is defined as :math:`\epsilon_{\alpha\beta} = (\eta_{\alpha\beta} + \eta_{\beta\alpha}) / 2`. The latter is used because the antisymmetric part of the strain tensor :math:`\eta` represents a pure rotation and does not contribute to the total energy change. So, using the Cauchy strain, the 9 strains are reduced to 6 independent ones.

Under such strain, we can use Equation :eq:`actualPotential` to calculate the potential energy increase by replacing :math:`u_i` by:

.. math::
   u_{R\tau\alpha}(t) = \eta_{\alpha\beta} (R+\tau)^{\beta} + u_{\tau\alpha}^0 + y_{R\tau\alpha}(t) = S_{R\tau\alpha} + y_{R\tau\alpha}(t),

where :math:`S` is the static displacement, containing an extra term :math:`u^0` which represents the relaxation of atom :math:`\tau` within the primitive cell, in addition to that due to the uniform deformation dictated by :math:`\eta`. This can occur in low-symmetry crystals, where atoms may not completely follow the uniform deformation and require an extra displacement :math:`u_{\tau\alpha}^0` to minimize the potential energy. The last term :math:`y_{R\tau\alpha}(t)` represents the dynamical motion of the atom :math:`R\tau`, or phonon degree of freedom, which has a time average of zero by construction. Plugging into Equation :eq:`actualPotential`, we can derive the effective harmonic potential at finite temperatures as the coefficient of the :math:`y^2` term, and the elastic tensor as the second derivative of the average potential energy with respect to the applied (Cauchy) strain:

.. math::
   C_{\alpha\beta,\gamma\delta} = \frac{1}{V_0} \frac{\partial^2 E}{\partial \epsilon_{\alpha\beta} \partial \epsilon_{\gamma\delta}}.

Note that as a second derivative, :math:`C` must be symmetric under swapping of the order of differentiation :math:`\alpha\beta \leftrightarrow \gamma\delta`. Due to the symmetry of the Cauchy strain tensor :math:`\epsilon` itself, the elastic tensor :math:`C` will also be invariant with respect to swaps of :math:`\alpha \leftrightarrow \beta` and :math:`\gamma \leftrightarrow \delta` separately.

...

Finally, the isotropic and uniaxial elastic constants are denoted respectively by the bulk and Young moduli :math:`B` and :math:`Y`.

**Bulk Modulus**

The isothermal bulk modulus is defined as:

.. math::
   B_T = -V \left(\frac{dP}{dV} \right)_{P=0,T}.

The isotropic compression condition translates into :math:`\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = -P`, :math:`\sigma_{\alpha\beta} = 0` for :math:`\alpha \ne \beta`. Since the volume change is given by :math:`\Delta V / V = \epsilon_{xx} + \epsilon_{yy} + \epsilon_{zz} = \mathrm{Tr}(\epsilon)`, one can express :math:`B` in terms of the elastic or compliance tensor elements as:

.. math::
   B = \frac{\mathrm{Tr}(\sigma)/3}{\mathrm{Tr}(\epsilon)} = \frac{1}{\sum_{i,j=1,2,3} S_{ij}} = B_{Reuss},

where :math:`S` is the inverse of :math:`C` when written as a :math:`6 \times 6` matrix in the Voigt notation, and is called the *compliance tensor*. This is the definition of the Reuss value; the Voigt value is :math:`B_{Voigt} = \frac{1}{9} \sum_{i,j=1,2,3} C_{ij}`. The two are not equal in general. For crystals of cubic symmetry they are, and simplify to:

.. math::
   \frac{C_{11} + 2C_{12}}{3}.

...

**Poisson Ratio**

For a uniaxial deformation along the x-axis, the Poisson ratio is defined as:

.. math::
   \nu = -\left(\frac{\partial \epsilon_{2}}{\partial \epsilon_{1}} \right)_{\sigma_{i \ne 1} = 0} = -\frac{S_{21}}{S_{11}}.

For isotropic materials, the value is typically reported as:

.. math::
   \nu_{\mathrm{iso}} = \frac{1.5B - \mu_{\mathrm{iso}}}{3B + \mu_{\mathrm{iso}}}.




Gruneisen parameter
-------------------

The mode Gruneisen parameters are obtained from the phonons and the cubic force constants. They represent the change in the phonon frequencies under applied strain. Even though the strain is usually considered to be isotropic compression or dilation, one can also define the response to an arbitrary strain. By definition, we have:

.. math::
   \gamma_{\lambda,\mu\nu} = - \frac{d \ln \,\omega_{\lambda} }{3 \, d \epsilon_{\mu\nu}}.

Usually, the volume derivative, which corresponds to an isotropic strain, is considered:

.. math::
   \gamma_{\lambda} = - \frac{d \ln \,\omega_{\lambda} }{d \ln \, V},

where :math:`V` is the unit cell volume. Applying a uniform hydrostatic strain:

.. math::
   \epsilon_{\alpha \beta} = \delta_{\alpha \beta} \, \frac{dV}{3V} = \delta_{\alpha \beta} \, \frac{d \ln V}{3}.

Similar to the group velocity, we start with the square of the frequency written as the product of the dynamical matrix by eigenvectors:

.. math::
   \omega^2 = e \cdot D \cdot e.

We then consider the change in :math:`D` under strain:

.. math::
   \frac{dD}{d\eta} = \sum_R \frac{d\Phi}{d\eta} \, e^{i k \cdot R},

and use:

.. math::
   \frac{d\Phi}{d\eta} = \Psi \left(R + \tau + \frac{d u^0}{d \eta}\right).

Note that the latter term, which is non-zero for low-symmetry structures having their own internal relaxation under strain :math:`u^0(\eta)`, is often neglected in the literature.

In the previous section, we derived this equation for :math:`u^0` in the harmonic approximation:

.. math::
   u^0(\eta) = -\Gamma (\Pi + Q \eta).

Therefore:

.. math::
   \frac{d u^0}{d \eta} = -\Gamma Q.

Using the Hellmann-Feynman theorem applied to the derivative of the dynamical matrix, and substituting the above relation for the strain derivative of the harmonic force constants, we finally find:

.. math::
   \gamma_{\lambda}(k) = -\frac{1}{6 \, \omega_{\lambda}^2(k)} \sum_{\tau, R_1 \tau_1, R_2 \tau_2}
   \frac{\Psi_{\tau, R_1 \tau_1, R_2 \tau_2}^{\alpha \alpha_1 \alpha_2}}{\sqrt{m_{\tau} m_{\tau_1}}} \,
   e_{\tau \alpha, \lambda}(k)^* \, e_{\tau_1 \alpha_1, \lambda}(k) \, e^{i k \cdot (R_1 + \tau_1 - \tau)} \,
   \left[(R_2 + \tau_2)^{\alpha_2} - \sum_{\mu=1,2,3} \Gamma_{\tau_2 \alpha_2, .} Q_{., \mu \mu} \right].

The temperature-dependent Gruneisen parameter is defined as the thermal average of the above modes weighted by their heat capacity:

.. math::
   \gamma(T) = \frac{\sum_{k \lambda} c_{k \lambda}(T) \gamma_{\lambda}(k)}{\sum_{k \lambda} c_{k \lambda}(T)},

where the heat capacity per mode is:

.. math::
   c_{k \lambda}(T) = k_B \frac{x_{k \lambda}^2}{\sinh^2 x_{k \lambda}}, \quad x_{k \lambda} = \frac{\hbar \omega_{\lambda}(k)}{2 k_B T}.


Lattice thermal conductivity (THERMACOND)
-----------------------------------------

THERMACOND computes the lattice thermal conductivity from the harmonic and cubic force
constants of FOCEX by solving the linearized phonon Boltzmann transport equation (BTE) over the
irreducible wedge of the Brillouin zone, with the energy-conserving delta functions integrated by
the tetrahedron method [Sadeghi2025]_. A phonon state :math:`k = (\lambda, \mathbf k)` is a
branch :math:`\lambda` at wave vector :math:`\mathbf k`, with frequency :math:`\omega_k`, group
velocity :math:`\mathbf v_k` and equilibrium (Bose-Einstein) occupation :math:`n^0_k`.

**The linearized BTE.** In a small temperature gradient the occupation is
:math:`n_k = n^0_k + n^1_k`, and in the steady state the drift balances the collisions:

.. math::
   :label: lbte

   -\mathbf v_k \cdot \nabla T \, \frac{\partial n^0_k}{\partial T}
   = -\frac{\partial n^1_k}{\partial t}\bigg|_{coll} = \sum_{k'} C_{kk'} \, n^1_{k'} .

Writing :math:`n^1_k = n^0_k (n^0_k + 1)\, \psi_k`, the collision term is

.. math::
   :label: collision

   \sum_{k'} C_{kk'} n^1_{k'} = \sum_{k'k''} \Big[ P^{k''}_{kk'} (\psi_k + \psi_{k'} - \psi_{k''})
   + \tfrac{1}{2} P^{k'k''}_k (\psi_k - \psi_{k'} - \psi_{k''}) \Big]
   + \sum_{k'} P^{k'}_k (\psi_k - \psi_{k'}) ,

with the rates of the coalescence and decay three-phonon processes and of the mass-disorder
(isotope) scattering given by Fermi's golden rule,

.. math::

   \begin{aligned}
   P^{k''}_{kk'} &= 2\pi\, n^0_k\, n^0_{k'} (n^0_{k''} + 1)\, |V_3(-k,-k',k'')|^2\,
                    \delta(\omega_k + \omega_{k'} - \omega_{k''}) ,\\
   P^{k'k''}_k   &= 2\pi\, n^0_k (n^0_{k'} + 1)(n^0_{k''} + 1)\, |V_3(-k,k',k'')|^2\,
                    \delta(\omega_k - \omega_{k'} - \omega_{k''}) ,\\
   P^{k'}_k      &= \frac{\pi}{2}\, \omega_k \omega_{k'}\, n^0_k (n^0_{k'} + 1)
                    \Big( \sum_\tau g_\tau |\mathbf e^*_{\tau k} \cdot \mathbf e_{\tau k'}|^2 \Big)\,
                    \delta(\omega_k - \omega_{k'}) ,
   \end{aligned}

where :math:`V_3` is the three-phonon matrix element built from the cubic force constants and
the eigenvectors, and :math:`g_\tau = \sum_i f_{\tau,i} (1 - M_{\tau,i}/M_{\tau,avg})^2` the mass
variance of the isotopes of atom :math:`\tau`. (THERMACOND 2.0 does not include the isotope term
yet: the conductivity jobs stop if it is requested on line 11 of ``latdyn.params``.)

**Relaxation-time approximation (RTA).** If every mode but :math:`k` is kept at equilibrium
(:math:`\psi_{k'} = \psi_{k''} = 0`), the collision matrix is diagonal,
:math:`\partial n^1_k / \partial t|_{coll} = Q_k \psi_k` with

.. math::

   Q_k = \sum_{k'k''} \Big( P^{k''}_{kk'} + \tfrac{1}{2} P^{k'k''}_k \Big) + \sum_{k'} P^{k'}_k ,
   \qquad \tau_k = \frac{n^0_k (n^0_k + 1)}{Q_k} ,

and the conductivity is

.. math::
   :label: kappa_rta

   \kappa^{\alpha\beta}_{RTA} = \frac{1}{\Omega N_k} \sum_k C_k \, v^\alpha_k v^\beta_k \, \tau_k ,
   \qquad C_k = \hbar\omega_k \frac{\partial n^0_k}{\partial T} ,

with :math:`\Omega` the volume of the primitive cell and :math:`N_k` the number of k-points.
The RTA treats the momentum-conserving normal processes as resistive and so usually
underestimates :math:`\kappa`; it is accurate when Umklapp and impurity scattering dominate.

**Full solution.** With :math:`\psi_k = \sum_\alpha F_{k\alpha}\, \partial T/\partial x_\alpha`,
Eq. :eq:`lbte` becomes a linear system for the distribution :math:`\mathbf F_k`. Only the
irreducible points :math:`k^*` are kept: a vector at :math:`k = \mathcal S\, k^*` is obtained by
the point-group operation :math:`\mathcal S`, :math:`\mathbf v_k = \mathcal S \mathbf v_{k^*}`,
:math:`\mathbf F_k = \mathcal S \mathbf F_{k^*}`, so that

.. math::
   :label: bte_ibz

   \begin{aligned}
   F^\alpha_{k^*} = F^{0\alpha}_{k^*} &+ \frac{1}{Q_{k^*}} \sum_{k'k''} \Big[
     P^{k''}_{k^*k'} \big( \mathcal S^{\alpha\beta}_{k''} F^\beta_{k''^*} - \mathcal S^{\alpha\beta}_{k'} F^\beta_{k'^*} \big)
     + \tfrac{1}{2} P^{k'k''}_{k^*} \big( \mathcal S^{\alpha\beta}_{k'} F^\beta_{k'^*} + \mathcal S^{\alpha\beta}_{k''} F^\beta_{k''^*} \big) \Big] \\
     &+ \frac{1}{Q_{k^*}} \sum_{k'} P^{k'}_{k^*}\, \mathcal S^{\alpha\beta}_{k'} F^\beta_{k'^*} ,
   \end{aligned}

with the RTA solution

.. math::

   F^{0\alpha}_{k^*} = \frac{\hbar\omega_{k^*} n^0_{k^*} (n^0_{k^*} + 1)\, v^\alpha_{k^*}}{k_B T^2 Q_{k^*}} ,

and

.. math::
   :label: kappa_full

   \kappa^{\alpha\beta} = \frac{1}{N_k \Omega} \sum_{k^*, \mathcal S} \hbar\omega_{k^*}\, n^0_{k^*} (n^0_{k^*} + 1)\,
   \mathcal S^{\alpha\alpha'} v^{\alpha'}_{k^*} \, \mathcal S^{\beta\beta'} F^{\beta'}_{k^*} .

:math:`\mathbf F^0` alone gives back the RTA. THERMACOND solves Eq. :eq:`bte_ibz` in two ways
(``latdyn.params`` line 8, :doc:`runthermacond`):

* **iterative** (job 3) [Omini1995]_ [Broido2007]_: starting from :math:`\mathbf F^0`, the
  right-hand side is evaluated and substituted again until the conductivity changes by less than
  a tolerance. It needs no stored collision matrix, but its convergence is not guaranteed in
  general [Fugallo2013]_ and it needs dense meshes;
* **direct** (job 2): the collision matrix projected onto the irreducible wedge,
  :math:`M_{k^*k'^*} = \sum_{k''} C_{k^*k''} \mathcal S_{k''k'^*}`, a square matrix of the number
  of irreducible points times the number of modes, gives
  :math:`\mathbf v_{k^*} \partial n^0_{k^*}/\partial T = \sum_{k'^*} M_{k^*k'^*} \mathbf F_{k'^*}`,
  solved by singular value decomposition or conjugate gradients.

**Degenerate modes.** Eigenvectors of degenerate modes are defined only up to a rotation within
their degenerate set, so a result that depends on one velocity per mode is not unique there.
THERMACOND averages :math:`|V_3|^2` and the scattering rates over each degenerate set, so that all
partners of a set have the same lifetime, and takes the conductivity of the set from the trace
over the set of the velocity matrix (the trace form of the RTA term of Eq. :eq:`kappa_wigner`
below), which does not depend on the choice of eigenvectors.

**Three contributions to** :math:`\kappa`. Simoncelli, Marzari and Mauri [Simoncelli2019]_ derived
from the Wigner phase-space distribution a conductivity that, besides the particle-like
propagation of phonons, contains the wave-like tunnelling between pairs of modes at the same
wave vector. In the RTA it reads

.. math::
   :label: kappa_wigner

   \kappa^{\alpha\beta} = \frac{1}{\Omega N_k} \sum_{\mathbf k} \sum_{s,s'}
   \frac{\omega_s + \omega_{s'}}{4} \left( \frac{C_s}{\omega_s} + \frac{C_{s'}}{\omega_{s'}} \right)
   \mathrm{Re}\big( V^\alpha_{ss'} V^\beta_{s's} \big)\,
   \frac{\tfrac{1}{2}(\Gamma_s + \Gamma_{s'})}{(\omega_s - \omega_{s'})^2 + \tfrac{1}{4}(\Gamma_s + \Gamma_{s'})^2} ,

with :math:`s, s'` the branches at :math:`\mathbf k`, :math:`\Gamma_s = 1/\tau_s` the linewidths
and :math:`V_{ss'} = \langle e_s | \partial D / \partial \mathbf k | e_{s'} \rangle / (\omega_s + \omega_{s'})`
the velocity operator, whose diagonal elements are the group velocities. The terms with
:math:`s = s'` (and pairs of degenerate modes) give the RTA conductivity :math:`\kappa_{RTA}` of
Eq. :eq:`kappa_rta` (Simoncelli et al. call it the populations conductivity, :math:`\kappa_P`);
the others give the coherences conductivity :math:`\kappa_C`. The conductivity THERMACOND reports
is thus made of three parts:

.. math::
   :label: kappa_three

   \kappa = \underbrace{\kappa_{RTA}}_{\text{RTA}}
          + \underbrace{(\kappa_{full} - \kappa_{RTA})}_{\text{beyond the RTA}}
          + \underbrace{\kappa_C}_{\text{coherences}} ,

with :math:`\kappa_{full}` the full solution of Eq. :eq:`kappa_full`:

* :math:`\kappa_{RTA}`, the particle-like conductivity in the RTA (column "RTA" of
  ``kappa_vs_T.dat`` and column ``kappa_P`` of ``kappa_wigner.dat``);
* :math:`\kappa_{full} - \kappa_{RTA}`, the correction of the full solution of the BTE, iterative
  or direct, over the RTA: the collective response to momentum-conserving normal processes, large
  when they dominate (``kappa_vs_T.dat`` gives :math:`\kappa_{full}`, i.e. :math:`\kappa_{RTA}`
  plus this correction);
* :math:`\kappa_C`, the wave-like part inherited from the Wigner formulation (``kappa_wigner.dat``),
  evaluated with the RTA linewidths. It is negligible in simple crystals with well-separated
  branches and becomes important in complex, strongly anharmonic or low-conductivity materials,
  where the linewidths are comparable to the spacing of the branches.

.. note::

   **Units of the linewidth.** :math:`\Gamma_s` in Eq. :eq:`kappa_wigner` is the full width at
   half maximum in angular frequency, :math:`\Gamma_s = 2\Sigma''_s = 1/\tau_s`, with
   :math:`\omega` in rad/s; the half width is :math:`\Sigma'' = 1/(2\tau)`. In cycles (THz) the
   same widths are :math:`\Gamma/2\pi = 1/(2\pi\tau)` and :math:`\Sigma''/2\pi = 1/(4\pi\tau)`.
   THERMACOND takes :math:`\tau_s` from the RTA, so that the :math:`s = s'` terms of
   Eq. :eq:`kappa_wigner`, :math:`C_s v_s^2/\Gamma_s`, are exactly the RTA conductivity.


.. [Leibfried-Ludwig] G. Leibfried and W. Ludwig, Theory of anharmonic effects in crystals, in Solid State Physics, edited by F. Seitz and D. Turnbull (Academic, New York, 1961), Vol. 12, p. 275. `doi:10.1016/S0081-1947(08)60656-6 <https://doi.org/10.1016/S0081-1947(08)60656-6>`__

.. [Broyden] Broyden, C.G., 1965. A class of methods for solving nonlinear simultaneous equations. Mathematics of Computation 19, 577. 

.. [Esfarjani-Liang2022] K. Esfarjani and Y. Liang, Equilibrium and non-equilibrium lattice dynamics of anharmonic systems, Entropy 24, 1585 (2022). `doi:10.3390/e24111585 <https://doi.org/10.3390/e24111585>`__

.. [Esfarjani2025] K. Esfarjani, H. Stokes, S. Nayeb Sadeghi, Y. Liang, B. Timalsina, H. Meng, J. Shiomi, B. Liao and R. Sun, ALATDYN: A set of Anharmonic LATtice DYNamics codes to compute thermodynamic and thermal transport properties of crystalline solids, Comput. Phys. Commun. 312, 109575 (2025). `doi:10.1016/j.cpc.2025.109575 <https://doi.org/10.1016/j.cpc.2025.109575>`__

.. [Wallace1998] Wallace, D.C., 1998. Thermodynamics of Crystals. Dover Publications, New York.

.. [Sadeghi2025] S. Nayeb Sadeghi, S. Lee and K. Esfarjani, THERMACOND, a code to compute lattice thermal conductivity from harmonic and anharmonic force constants, npj Comput. Mater. 11, 303 (2025). `doi:10.1038/s41524-025-01673-8 <https://doi.org/10.1038/s41524-025-01673-8>`__

.. [Simoncelli2019] M. Simoncelli, N. Marzari and F. Mauri, Unified theory of thermal transport in crystals and glasses, Nat. Phys. 15, 809 (2019).

.. [Omini1995] M. Omini and A. Sparavigna, An iterative approach to the phonon Boltzmann equation in the theory of thermal conductivity, Physica B 212, 101 (1995).

.. [Broido2007] D. A. Broido, M. Malorny, G. Birner, N. Mingo and D. A. Stewart, Intrinsic lattice thermal conductivity of semiconductors from first principles, Appl. Phys. Lett. 91, 231922 (2007).

.. [Fugallo2013] G. Fugallo, M. Lazzeri, L. Paulatto and F. Mauri, Ab initio variational approach for evaluating lattice thermal conductivity, Phys. Rev. B 88, 045430 (2013).
