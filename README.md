# Optiflow

Optiflow started as a tool for optimising wind turbine airfoils and hydrofoils. It mains applications include:

* Routines for rotor planform optimization applicable to conceptual design of wind (and tidal) turbines
* Minimization routines for data-driven calibration of viscous-inviscid interaction codes (stochastic and deterministic)

Most applications combine external viscous-inviscid codes (e.g. Xfoil or Rfoil) with robust internal solvers: 
* An implicit steady state blade-element-momentum (BEM) solver with explicit residuals
* A Hess-Smith panel solver used to estimated the effect of inflow curvature on airfoils
* An integral boundary layer solver used to assess the stall margin of optimized airfoils
* A hybrid integral+finite-difference boundary layer code for modelling vortex generators

# Requirement
