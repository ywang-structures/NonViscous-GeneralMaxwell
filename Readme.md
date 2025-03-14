## Introduction
This repository provides MATLAB scripts for analyzing the relationship between a pair of structures with different damping models but providing equivalent dynamical behaviors. Within the pair, one structure adopts the exponential non-viscous damping model and the other adopts a generalized Maxwell damping model.  Accompanying the manuscript, one m-file (Example_3DOF.m) analyzes the relationship between a pair of 3-DOF examples using the two different damping models.  The other m-file (Example_5DOF.m) analyzes a pair of 5-DOF examples.

## Companion Monograph
Y. Otsuki and Y. Wang. (2025). "On the relationship between exponential non-viscous damping and generalized Maxwell model", ASCE Journal of Engineering Mechanics (in review).

## Contents
### `Example_3DOF.m`

Using a pair of 3-DOF systems, this script first demonstrates the procedure for obtaining the transformation matrix that relates the internal velocities of exponential non-viscous damping to the damper velocities of the generalized Maxwell model. Additionally, it analyzes the modal response of both oscillatory and non-oscillatory modes in exponential non-viscous damping. The results presented in Section 7 of the manuscript are generated using this script.

<img src="Figures/3DOF.png?raw=false" width="480" height="360" />

### `Example_5DOF.m`

This script demonstrates the procedure for obtaining the transformation matrix between the pair of 5-DOF systems. The configurations are devised to be more generic, including multiple connections with the ground. Dynamic response is simulated using the 1940 El Centro NS ground motion. The results of Section 8 of the manuscript are generated using this script.

<img src="Figures/5DOF.png?raw=false" width="550" height="360" />
