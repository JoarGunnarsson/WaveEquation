# Wave Equation Numerical Solver

This project provides a MATLAB script to numerically solve the wave equation in a 2D domain. The solver is capable of handling both Dirichlet and Neumann boundary conditions, along with arbitrary initial conditions. Additionally, it supports optional features such as forcing a sinus wave motion in the center of the domain and introducing an optional energy friction term. This script uses the central finite differences scheme. This can all be found in the file `WaveEquationNumerical.m`. This repository also includes an example analytical solution of the wave equation. More information can be found in the file `WaveEquationAnalytical.m`.

## Usage

To use the wave equation solver, follow these steps:

1. Open the MATLAB script `WaveEquationNumerical.m` in a MATLAB environment.

2. Adjust the parameters in the script to fit your specific problem:
   - `L_x`: Length of the domain in the x-direction.
   - `L_y`: Length of the domain in the y-direction.
   - `T`: Total simulation time.
   - `N_x`: Number of spatial points in the x-direction.
   - `N_y`: Number of spatial points in the y-direction.
   - `N_t`: Number of time steps.
   - `c`: The wave speed.
   - Make sure that `stability_constant` is less than 1. This can be changed by modifying the number of spatial points, time steps, and the wave speed.
   - `mu`: Parameter determining the energy loss due to friction.
   - Any additional changes to the boundary- and initial- conditions.
3. Run the first section of the script to compute the numerical solution to the wave equation.
4. Run the second section in order to visualize the solution in 3D.

## Results

The first section of the script generates a numerical solution to the wave equation.

## Visualization

The second section of the script script generates a 3D plot of the wave evolution over time.

## Additional Notes

- The script provides an example of an initial velocity function (`v_0`) that you can modify to suit your specific needs.

- The solver supports both Dirichlet and Neumann boundary conditions. You can specify these conditions in the script as needed.

- Optional features, such as forcing a wave motion in the center of the domain and introducing an energy friction term, are demonstrated in the script and can be modified or disabled as required.
