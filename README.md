# Physics Simulations

This repository contains a collection of numerical physics simulations focused on gravitational N-body systems and orbital dynamics. The code spans coursework projects and research work completed during my internship at the Center for Computational Astrophysics (CfCA).

All simulations implement custom numerical integrators in C or Python, primarily using leapfrog and Hermite integration schemes.

## Projects
### Collisionless Cold Collapse (C)

Simulation of the collapse of a collisionless isothermal particle cloud under gravity. Implements a leapfrog integrator for long-term stability in Hamiltonian systems.

### Galaxy Collision (C)

N-body simulation of colliding Plummer sphere galaxies using a leapfrog integration scheme. Includes utilities for initialization, force calculation, and time evolution of particle systems.

### Kepler Two-Body / Hermite Integrator Validation (C)

Two-body orbital simulation used to validate a predictor-corrector Hermite integrator (P(EC)^N scheme). Includes orbital element conversion utilities and error analysis for long-term orbital stability.

### N-Body Force Simulator (Python / Jupyter)

Capstone coursework project implementing a generalized N-body force simulation in Python using a 4th-order Runge-Kutta integrator. Includes a custom particle system and visualization of multi-body gravitational dynamics.

## Notes

These projects were developed independently as part of coursework and research. While implemented separately, they all explore numerical integration methods for gravitational systems and orbital dynamics.
