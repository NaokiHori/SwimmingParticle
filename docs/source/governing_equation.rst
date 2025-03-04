
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

*******************
Boundary Conditions
*******************

The azimuthal direction is periodic by definition.
The radial boundary conditions and their implications are discussed below.

==============
Velocity Field
==============

The domain is bounded by walls in the radial direction, specifically at :math:`\vr = 1` and :math:`\vr = R`.
These walls satisfy the impermeability condition:

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

Since the Stokes equation is a pure boundary-value problem, the flow field is uniquely determined once the aforementioned boundary conditions are prescribed.
In particular, the solution can be conveniently expressed using the stream function :math:`\psi \left( \vr, \vt, t \right)`, given by:

.. math::

    \psi \left( \vr, \vt, t \right)
    =
    \sum_k
    \Psi_k \left( \vr, t \right)
    \expp{I k \vt},

where :math:`I` is the imaginary unit, and

.. math::

    \Psi_k \left( \vr, t \right)
    =
    \frac{1 - \vr^2}{2 \vr^{\left| k \right|}}
    I k M C_k^s,

following |HU2019|.

Here, :math:`C_k^s` represents the concentration field on the particle surface (at :math:`\vr = 1`) in the frequency domain:

.. math::

    c \left( \vr = 1, \vt, t \right)
    =
    \sum_k
    C_k^s \expp{I k \vt}.

Note that, although this condition is not strictly met due to the presence of the outer wall, we assume it as a good approximation when :math:`R` is sufficiently large.

===================
Concentration Field
===================

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

    Despite not used in this project for now, the velocity components are obtained by differentiating the stream function in the frequency space:

    .. math::

        \left( U_{\vr} \right)_k
        \equiv
        -
        k^2
        \frac{
            1 - {\vr}^2
        }{
            2 {\vr}^{\left| k \right| + 1}
        }
        M C_k^s,

    .. math::

        \left( U_{\vt} \right)_k
        \equiv
        \left[
            \left( 1 - \frac{\left| k \right|}{2} \right)
            \frac{1}{{\vr}^{\left| k \right| - 1}}
            +
            \frac{\left| k \right|}{2} \frac{1}{{\vr}^{\left| k \right| + 1}}
        \right]
        i k M C_k^s.

    It is readily apparent that they satisfy the desired boundary conditions on the particle surface.

