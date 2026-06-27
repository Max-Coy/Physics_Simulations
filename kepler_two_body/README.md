# Kepler Two-Body / Hermite Integrator Testbed

C implementation of a two-body gravitational system used to validate a 4th-order Hermite N-body integrator with adaptive timestepping.

The integrator implements a predictor-corrector (PEC)^n Hermite scheme with individual particle timesteps chosen as powers of two, enabling adaptive resolution across dynamical timescales.

In addition to the integrator, the project includes utilities for converting between orbital elements and Cartesian coordinates, supporting analysis of Keplerian orbits and energy conservation.

This system was used to test numerical stability and accuracy of long-term orbital integration.