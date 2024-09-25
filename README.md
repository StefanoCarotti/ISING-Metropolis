# Ising Model 3D Monte Carlo Simulation

## Overview
This project implements a Monte Carlo simulation for the 3D Ising model using the Metropolis algorithm. It calculates various thermodynamic properties of the system, including magnetization, energy, susceptibility, and heat capacity, at different temperatures (represented by β values). The code also computes the Binder cumulant to identify phase transitions in the system.

## Features
- **3D Ising Model Simulation**: Simulates a cubic lattice with spins that interact according to the Ising Hamiltonian.
- **Metropolis Algorithm**: Used for updating the spin configurations based on thermal fluctuations.
- **Observable Quantities**: Computes the following observables:
  - Magnetization ⟨M⟩
  - Energy ⟨E⟩
  - Susceptibility χ
  - Heat Capacity C
  - Binder Cumulant U
- **Phase Transition Detection**: Includes graphs of physical quantities to identify phase transitions at critical β values.

## Installation
To run the code, ensure you have MATLAB or Octave installed. Clone the repository and simply run the script in your MATLAB environment.

## Usage
1. Define the number of spins `N`, coupling constant `J`, and number of sweeps.
2. Adjust the temperature range by modifying the `beta` variable.
3. Run the script to perform the Monte Carlo sweeps and collect statistical measurements.
4. Visualize the results using the generated graphs for magnetization, energy, susceptibility, heat capacity, and Binder cumulant.

## Output
- **Magnetization**: Plots of average magnetization as a function of β, compared to theoretical predictions.
- **Energy**: Plots of average energy per spin, showing how the system’s internal energy changes with temperature.
- **Susceptibility**: Graph of magnetic susceptibility, identifying the critical temperature.
- **Heat Capacity**: Visualization of the system's specific heat, another indicator of phase transitions.
- **Binder Cumulant**: Plot of the Binder cumulant, used to locate the critical temperature for phase transitions.
