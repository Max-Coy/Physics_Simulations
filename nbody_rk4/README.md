# N-Body Force Simulation Framework

Python-based N-body simulation framework implementing numerical integration of particle systems under arbitrary force laws.

The system uses a custom particle class and a 4th-order Runge-Kutta integrator to solve equations of motion for interacting bodies. While originally developed as a coursework project, the simulation is designed to support general force functions beyond classical gravity.

A key feature of the framework is its modular force architecture, allowing users to define custom interactions such as gravitational, electrostatic, and synthetic force models within the same simulation loop.

The implementation is written in a Jupyter notebook and includes tools for visualization and animation using FFmpeg.