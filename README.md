# Optiflow

 10.5281/zenodo.3520280

Optiflow started as a tool for optimising wind turbine airfoils and hydrofoils. It mains applications include:

* Routines for rotor planform optimization applicable to conceptual design of wind (and tidal) turbines
* Minimization routines for data-driven calibration of viscous-inviscid interaction codes (stochastic and deterministic)

Most applications combine external viscous-inviscid codes (e.g. Xfoil or Rfoil) with robust internal solvers: 
* An implicit steady state blade-element-momentum (BEM) solver with explicit residuals
* A Hess-Smith panel solver used to estimated the effect of inflow curvature on airfoils
* An integral boundary layer solver used to assess the stall margin of optimized airfoils
* A hybrid integral+finite-difference boundary layer code for modelling vortex generators
* A two stage least-sqaures solver for robust acquisition of Bernstein parameters used Class-Shape-Transform representation of 

Complex cases and workflows can be constructed thanks to the object-oriented nature of the code. The following algorithms are employed, amongst other:

* Genetic multi-objective (gamultiobj, NSGA-II) and gradient single-objective (fmincon, hybrid line-shearch with BFGS hessian estimation) optimization algorithms provided the by Matlab Optimization Toolbox. 
* Stochastic gradient descent minimization relies on an inhouse implementation of Hinton's algorithms.
* Shape parametrization is mostly handled with the class-shape-transform (CST) technique, complemented by curvilinear coordinate transformations (Lamé 1859) for flap deformations.
* Paralelization is mostly coarse grained, and has been shown to exhibits good scaling on shared memory architectures with up to 128 cores in the Azure cloud. Scaling is better for code calibration, while optimization run most efficiently at n time 32 cores however.

All code is shared under the attached MIT license within the preparation of a publication about the hydrofoil and tidal turbine design, which will extend:

* Espenica F., Santos Pereira R.B., de Oliveira G.L., Baltazar J. and J. A. C. Falcão de Campos. "Design and Optimization of Hydrofoils Tailored for Marine Current Turbines" in Proceedings of th 13th European Wave and Tidal Energy Conference, Naples, 2019


### How do I get set up? ###

Ideally, your system should meet the following requirements:  
  1 Recent version of Mac Os X (eg. 10.12) or Linux (eg. Ubuntu 16.04LTS)  
  2 Recent version of Matlab (eg. R2018a) with all optimization toolboxes
  3 Access to Rfoil MK3 compiled with GNU Fortran 4.8+ 

To start, the best is to explaore the example case files which should work out of the box:
  1 Open Matlab and cd to the root of the repository (always do this on a local copy!)
  2 Open a case file, for example *A0_start_example.m*. The file header describes the case.
  3 Normally, you can just press run to initialize the case. It should work out of the box!
  4 Running the full case requires that you uncomment the last lines of most example cases (or insert them manually in the matlab prompt), because a full optimization case can take several hours or days (weeks for code calibration, if few cores are available)
  5 Don't forget to adapt N_cores (close to example file header) to the number of physical cores on your system (usually half the number of logical processors)

What if my system does not fulfill the ideal requirements. Troubleshooting:
  1 Problem: I am running windows. Solution: activate the ubuntu subsystem for windows and adapt simulation_worker class as needed.
  2 Problem: I don't have matlab or the global optimization toolbox. Solution: use Octave instead, most code will run out of the box, but some editions will be needed to the gamultiobj_manager class.
  3 Problem: I don't have access to Rfoil MK3. Solution: use Xfoil compiled with double precision using gfortran (Intel fortran is less reliable during parallel execution for this particular case). Minor editions will be required to the simulation_worker and simulation_protocol classes.

### What about examples? ###

Yes! The code is provided with several example cases *A0-P0* located in repository root to showcase some of the codes capabilities. A few examples (this list will be completed):

  A0. Baseline Parallel multiobjective optimization case  
  B0. Parallel multiobjective optimization case two multiple polars per operating point (explicit treatment of roughness sensitivity)  
  C0. Parallel multiobjective optimization case with many polars per operating point and dynamic shape deformation to optimize airfoils with **FLAPS** *(new)*   
  D0. Parallel multiobjective optimization case using a simulation_protocol_dynamizer object to optimize suction region lenght and location (Gael's MSc thesis cases)  
  E0. Parallel multiobjective optimization case using a shape_dynamizer object to restrict design space to **Prescribed Thickness** *(new)*  airfoils  
  F0. Baseline Serial multiobjective optimization case  

Deeper descriptions are provided in each case files header, comments and includes (for C).

### How should I cite this work? ###

Follow the MIT license, and use the Zenendo DOI identifier of this repository if want to cite this work. Also, feel welcome to cite one or more of the following publications if they are relevant to your work (they may also help understand what is going on here):

* G. de Oliveira, R. Pereira, W.A. Timmer and R.P.J.M. van Rooij. Improved airfoil polar predictions with data-driven boundary-layer closure relations. Journal of Physics Conference Series 1037(2):022009 , 2018
* R. Pereira, G. de Oliveira, W. A. Timmer and E. Quaeghebeur. Probabilistic Design of Airfoils for Horizontal Axis Wind Turbines. Journal of Physics Conference Series 1037(2):022042, 2018
* R. Pereira, W.A. Timmer, G. de Oliveira and G.J.W. van Bussel. Design of HAWT airfoils tailored for Active Flow Control. Wind Energy, 2017
* G. de Oliveira, W. A. Timmer and B.W. van Oudheusden. Integral Equations for Boundary Layers with Streamwise Vortices. Proceedings of teh 52nd 3AF International Conference on Applied Aerodynamics, Lyon, 2017.
* G. de Oliveira, M. Kotsonis and B.W. Van Oudheusden. Laminar Boundary Layer Flow with DBD Plasma Actuation: A Similarity Equation. Springer Lecture Notes in Computational Science (120), 2017.
* E.C. Battle, R. Pereira, M. Kotsonis and G. de Oliveira. Airfoil Optimisation for DBD Plasma Actuator in a Wind Energy Environment: Design and Experimental Study. Proceedings of 55th AIAA Aerospace Sciences Meeting, 2017.
* S. Bal, R. Pereira, G. de Oliveira and D. Ragni. Investigating the Influence of DBD Plasma Actuators on the Skin Friction in Integral Boundary Layer Formulation. Proceedings of 54th AIAA Aerospace Sciences Meeting, AIAA2016-2162, 2016.
* G. de Oliveira, R. Pereira, D. Ragni and M. Kotsonis. Modeling DBD Plasma Actuators in Integral Boundary Layer Formulation for Application in Panel Methods. Proceedings of the 46th AIAA Plasmadynamics and Lasers Conference, AIAA2015-3367, 2015.

### What about included data? ###

All data included here was digitized or generated by the present authors, except for the NACA-TR824 polar curves which where digitized by Gregory P.D. Siemens (University of Saskatchewan, 1994) and released into the the Public Domain. 

The authors believe reproducing data that was digitized from published literature is compliant with the intent and terms of original reports, and with the terms of the MIT license of the present repository. Including this data in input files is necessary to demonstrate the functionality of the code, and ensure reproducibility of the research mentioned in the  aforementioned publications. Reproducibiliy is not entirely complete, though, since some of the above publications employed data cannot be redistributed by the present authors (those files were removed prior to release, and replaced by similar open versions when possible).

Please contact us if you see any reason for concern regarding data sources. We encourage users to employ the readme_data_source.txt files at the root of relevant folders to cite original sources when performing academic work.

### Developper advice ###

The code is object oriented and resorts to templating and class overloading in so far as Matlab permits without hurting backward compatibility excessively. It is therefore suggested that new entrants spend some time trying to grasp the main coding conventions, which are explained in comments throughout the code. Extensive comments should make everything rather straightforward.


