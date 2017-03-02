# Rb-vapor-cell-modelling
Modelling of electrical susceptibility (absorption/refraction) and relaxation processes in Rb vapor cells. Also provides optimisation analysis of (microwave) magnetometer sensitivity. 

There are four example 'front-end' scripts:

cell_data:  Computes various properties for a vapor cell, and prints these in the command window

OD_spectra: Plots the optical density (absorption) spectrum for a Rb vapor cell: D1 or D2 line; arbitrary cell temperature, cell length, buffer gas mixtures; arbitrary background dc magnetic field amplitude. 

OD_vs_temperature: Plots the optical density (absorption) for a Rb vapor cell as a function of temperature: D1 or D2 line; arbitrary laser frequency, cell length, buffer gas mixtures; arbitrary background dc magnetic field amplitude

mw_sensitivity_optimisation: Produces 2D plots of microwave magnetic field sensitivity in a vapor cell as a function of buffer gas filling and vapor cell temperature. Both the photon-shot-noise and atomic-projection-noise limited sensitivities are given. The script also provides overlays of the corresponding spatial resolution, and plots of the optical density and T1 time.

