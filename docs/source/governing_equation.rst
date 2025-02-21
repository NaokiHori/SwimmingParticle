
.. include:: /reference.txt

##################
Governing Equation
##################

**********************
Differential Equations
**********************

We consider a two-dimensional domain where a circular solid particle is positioned at the center.
The liquid is assumed to be incompressible, leading to:

.. math::

    \pder{}{u_i}{x_i}
    =
    0.

We also assume that the Reynolds number is sufficiently small, resulting in the Stokes equation as the momentum balance:

.. math::

    0_i
    =
    -
    \pder{}{p}{x_i}
    +
    \pder{}{}{x_j}
    \pder{}{u_i}{x_j}.

The particle emits a solute, which is transported into the liquid.
The concentration of the solute (e.g., the non-dimensional mass per unit area) follows the advection-diffusion equation:

.. math::

    \pder{}{c}{t}
    +
    u_j
    \pder{}{c}{x_j}
    =
    \frac{1}{Pe}
    \pder{}{}{x_j}
    \pder{}{c}{x_j}.

The azimuthal direction is periodic by definition, while the domain is bounded by walls in the radial direction, specifically at :math:`\vr = 1` and :math:`\vr = R`.
Their implications are separately discussed below for the velocity and the concentration fields.

**************
Velocity Field
**************

The walls satisfy the impermeability condition:

.. math::

    \vat{\ur}{\vr = 1}
    =
    \vat{\ur}{\vr = R}
    =
    0.

The outer wall follows the no-slip condition:

.. math::

    \vat{\ut}{\vr = R}
    =
    0,

while a slip velocity, proportional to the azimuthal concentration gradient, is enforced on the inner wall:

.. math::

    \vat{\ut}{\vr = 1}
    =
    \vat{
        \frac{1}{r}
        \pder{}{c}{\vt}
    }{\vr = 1}.

Since the Stokes equation is a pure boundary-value problem, the entire flow field is uniquely determined once the aforementioned boundary conditions are prescribed.
In particular, the solution can be conveniently expressed using the stream function :math:`\psi \left( \vr, \vt, t \right)`, given by:

.. math::

    \psi \left( \vr, \vt, t \right)
    =
    \sum_k
    \Psi_k \left( \vr, t \right)
    \expp{I k \vt},

where :math:`I` is the imaginary unit, and :math:`\Psi_k` is the stream function in the frequency space.
Because of :math:`\psi \in \mathbb{R}`, :math:`\Psi_k` satisfies the complex-conjugate property, and thus it is sufficient to consider :math:`k \ge 0`.

As derived in :ref:`the appendix <appendix_stream_function>`, this is given by

.. math::

    \Psi_k
    =
    \begin{cases}
        k = 0
        &
        \text{unused}, \\
        k = 1
        &
        A_1 \vr^{-1}
        +
        B_1 \vr
        +
        E_1 \vr \log \vr
        +
        D_1 \vr^3, \\
        \text{otherwise}
        &
        A_k \vr^{-k}
        +
        B_k \vr^k
        +
        C_k \vr^{2 - k}
        +
        D_k \vr^{2 + k},
    \end{cases}

where the coefficients :math:`A_k, B_k, C_k, D_k, E_1` are determined using the boundary conditions at :math:`\vr = 1` and :math:`\vr = R`.
Assuming that :math:`R` is infinity, we request :math:`B_k, E_1, D_k` to be zero so that the solution converges, leading to

.. math::

    \Psi_k
    =
    A_k \vr^{-k}
    +
    C_k \vr^{2 - k}.

By imposing the boundary conditions at :math:`\vr = 1`:

.. math::

    \begin{pmatrix}
        1 & 1 \\
        - k & 2 - k \\
    \end{pmatrix}
    \begin{pmatrix}
        A_k \\
        C_k
    \end{pmatrix}
    =
    \begin{pmatrix}
        0 \\
        - \left( U_{\vt}^s \right)_k,
    \end{pmatrix}

we obtain

.. math::

    \Psi_k \left( \vr, t \right)
    =
    \frac{1 - \vr^2}{2 \vr^{k}}
    \left( U_{\vt}^s \right)_k.

In the above equation, :math:`U_{\vt}^s` is the azimuthal velocity along the particle surface, expressed in the frequency domain:

.. math::

    \left( U_{\vt}^s \right)_k
    =
    I k C_k^s,

or in the physical space:

.. math::

    \vat{\ut}{r = 1}
    =
    \vat{
        \frac{1}{\vr}
        \pder{}{\ut}{\vt}
    }{
        r = 1
    }.

:math:`C_k^s` represents the concentration field on the particle surface (at :math:`\vr = 1`) in the frequency domain:

.. math::

    c \left( \vr = 1, \vt, t \right)
    =
    \sum_k
    C_k^s \expp{I k \vt}.

*******************
Concentration Field
*******************

We impose the following conditions on the inner and outer boundaries, respectively:

.. math::

    \vat{\pder{}{c}{\vr}}{\vr = 1}
    =
    -1,

and

.. math::

    \vat{c}{\vr = R}
    =
    0.

.. note::

    In practice, we set :math:`R` to a finite value, which contradicts the obtained stream function derived under the assumption that :math:`R \to \infty`.
    This is originated by the fact that no analytical solution is available due to the non-linear term in the advection-diffusion equation.
    We assume, nevertheless, that the approximation remains sufficiently accurate for large values of :math:`R`.

