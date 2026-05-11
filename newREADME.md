# Aeroelastic Flutter Simulation (C++)

Reduced-order aeroelastic flutter simulation of a flexible aircraft wing using custom C++ code, tetrahedral mesh deformation, quasi-steady aerodynamics, and VTK visualization.

---

# Demo

## Flutter

![Flutter](docs/gif/flutter.gif)

---


# Features

- Custom C++ aeroelastic solver
- Reduced-order structural dynamics
- Bending + torsion wing modes
- Quasi-steady aerodynamic model
- Aeroelastic flutter instability
- Effective angle of attack model
- Aerodynamic phase lag
- Tetrahedral wing mesh
- VTK / ParaView visualization
- Stress computation inside tetrahedral elements
- Flutter-driven stress growth
- Self-excited oscillations
- Critical flutter velocity detection

---

# Physical Model

The project implements a simplified but physically meaningful aeroelastic flutter model.

The goal is to reproduce:

- damped oscillations at low flow velocity;
- self-excited oscillations near critical velocity;
- unstable flutter growth at high velocity.

---

# Coordinate System

- **X** — airflow direction
- **Y** — wing span direction
- **Z** — vertical direction

Freestream velocity:

$ \mathbf{U}_\infty = (U,0,0) $

The wing root is clamped at the root section.

---

# Structural Model

The wing is represented using a reduced-order modal model with two coupled degrees of freedom:

- bending
- torsion

---

## Bending Mode

Wing deflection:

$ w(y,t) = q_b(t)\phi_b(y) $

where:

- $q_b$ — generalized bending coordinate
- $\phi_b(y)$ — bending mode shape

Cantilever mode shape:

$ \phi_b(\xi) = \xi^2(3 - 2\xi) $

with:

$ \xi = \frac{-y}{L} $

---

## Torsion Mode

Wing twist:

$ \theta(y,t) = q_t(t)\phi_t(y) $

where:

- $q_t$ — generalized torsional coordinate
- $\phi_t(y)$ — torsion mode shape

Current implementation:

$ \phi_t(\xi) = \xi $

or optionally:

$ \phi_t(\xi) = \xi^2 $

---

## Small-Angle Kinematics

Wing torsion uses linearized rotation around the spanwise axis:

$ x' = x + z\theta $

$ z' = z - x\theta $

which is valid for small angular deformation.

---

# Structural Dynamics

The coupled structural system is:

$ M\ddot q + C\dot q + Kq = Q_{aero} $

with:

$ q =
\begin{bmatrix}
q_b \\
q_t
\end{bmatrix}
$

Diagonal reduced-order matrices are currently used:

$ M =
\begin{bmatrix}
m_b & 0 \\
0 & I_t
\end{bmatrix}
$

$ C =
\begin{bmatrix}
c_b & 0 \\
0 & c_t
\end{bmatrix}
$

$ K =
\begin{bmatrix}
k_b & 0 \\
0 & k_t
\end{bmatrix}
$

Time integration:

- semi-implicit Euler scheme

---

# Aerodynamic Model

The solver uses a quasi-steady aerodynamic model.

---

## Dynamic Pressure

$ q_\infty = \frac12 \rho U^2 $

---

## Effective Angle of Attack

The key aeroelastic mechanism is based on the effective angle of attack:


$$ \alpha_{eff}
=
\theta
+
\frac{-\dot w + x\dot\theta}{U}
$$

This includes:

- geometric torsion;
- vertical bending velocity;
- rotational velocity contribution.

This creates aerodynamic feedback between:

- bending;
- torsion;
- aerodynamic lift;
- aerodynamic moment.

---

## Aerodynamic Lag

To reproduce unsteady aeroelastic effects, a first-order aerodynamic lag model is used:

$$ \tau_a \dot\alpha_{lag}
+
\alpha_{lag}
=
\alpha_{eff}
$$

Lift is computed using the lagged angle:

$ C_L = C_{L\alpha}\alpha_{lag} $

This introduces:

- phase delay;
- aerodynamic memory;
- effective negative damping;
- flutter instability.

---

## Lift Force

Local lift force:

$ L = q_\infty C_L A $

where:

- $A$ — nodal aerodynamic area.

---

## Aerodynamic Torsional Moment

Aerodynamic moment around torsional axis:

$ M_y = -(x - x_{ac})L $

This couples aerodynamic lift to torsional motion.

---

# Flutter Mechanism

The flutter instability emerges from the feedback loop:

```text
Bending velocity
    ↓
Effective angle of attack
    ↓
Lift variation
    ↓
Aerodynamic torsional moment
    ↓
Wing twist
    ↓
Additional lift variation
    ↓
Growing oscillation
```

At low velocity:

$ C_{total} > 0 $

Oscillations decay.

Near critical velocity:

$ C_{total} \approx 0 $

Sustained oscillations appear.

Above flutter speed:

$ C_{total} < 0 $

Oscillation amplitude grows exponentially.

---

# Stress Computation

Stress is computed directly inside tetrahedral elements.

---

## Bending Stress

$ \sigma = Ez\kappa $

where:

- $E$ — Young modulus
- $\kappa$ — local curvature

---

## Torsional Shear Stress

$ \tau = Gr\frac{d\theta}{dy} $

where:

- $G$ — shear modulus
- $r$ — distance from torsional axis

---

## Equivalent Stress

Von Mises equivalent stress:

$$ \sigma_{eq}
=
\sqrt{\sigma^2 + 3\tau^2}
$$

Nodal stresses are obtained via tetrahedral averaging.

---

# Visualization

The project exports VTK files for ParaView visualization.

Available fields:

- displacement
- pressure
- effective angle of attack
- lift
- bending deformation
- torsional deformation
- equivalent stress

---

# Example Flutter Regimes

| Flow Velocity | Behavior |
|---|---|
| Low $U$ | Damped oscillation |
| Near $U_{crit}$ | Sustained oscillation |
| High $U$ | Divergent flutter |

---

# Numerical Stability Notes

The simulation uses:

- linear elasticity
- small-angle approximation
- reduced-order structural dynamics
- quasi-steady aerodynamics

The model is intentionally simplified to remain:

- computationally lightweight;
- numerically stable;
- easy to extend.

---

# Technologies

- C++
- VTK
- ParaView
- Gmsh
- Tetrahedral FEM mesh

---

# Future Improvements

Potential future extensions:

- nonlinear structural dynamics
- modal coupling matrices
- Theodorsen unsteady aerodynamics
- geometric nonlinearity
- adaptive timestep
- failure/destruction model
- multi-mode structural basis
- GPU acceleration

---

# Screenshots

## Mesh

<p align="center">
  <img src="docs/images/mesh.png" width="800"/>
</p>

---

## Stress Field

<p align="center">
  <img src="docs/images/stress.png" width="800"/>
</p>

---
