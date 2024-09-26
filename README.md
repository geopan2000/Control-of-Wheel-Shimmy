# Modelling and Control of Wheel Shimmy

## Overview
This repository contains the MATLAB simulations and analysis for the project on wheel shimmy modeling and control. The coursework involves deriving equations of motion for a multibody mechanical system, analyzing its stability under various conditions, and applying feedback control to stabilize the system. The work also examines the influence of stiffness, forward speed, and relaxation length on the system's stability.

## Coursework Tasks

### 1. Equations of Motion
- Derive the equations of motion for a multibody system, including the conversion between Cartesian and polar coordinates and the velocity and acceleration vectors.

### 2. Constant Forward Speed
- Apply a force \( F_x \) to maintain constant forward speed \( \dot{x} = v \). Derive the expression for \( F_x \) by setting \( \ddot{x} = 0 \).

### 3. Equations with Non-Holonomic Constraint
- Introduce a non-holonomic constraint to the system, setting the lateral velocity \( v_{\text{lat}} = 0 \). Derive the new set of equations under this constraint.

### 4. Model Examination
- Examine the properties of the system model, including the mass, inertia, and dimensions. Verify the model by analyzing its behavior in MATLAB.

### 5. Varying Stiffness
- Analyze how changing the stiffness affects the system's eigenvalues and stability. Simulate the system for varying stiffness values and identify the point where the system becomes marginally stable.

### 6. Validation of Stability
- Validate the stability of the system for stiffness values of \( k = 30,000 \) and \( k = 100,000 \) by examining the system's eigenvalues and simulating its response.

### 7. Forward Speed Analysis
- Investigate the effect of varying the forward speed on the system's stability. Simulate the system for a range of forward speeds and analyze the resulting eigenvalues.

### 8. Relaxation Length Analysis
- Examine how changes in the relaxation length affect the system's stability. Simulate the system for different relaxation lengths and identify the critical length beyond which the system becomes unstable.

### 9. Feedback Control
- Implement feedback control using the control law \( T = -\lambda \dot{\delta} \) to stabilize the system. Use Nyquist stability analysis to determine the required feedback gain \( \lambda \) to ensure marginal stability.

  

