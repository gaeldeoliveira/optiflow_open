# OptiFLOW #

OptiFLOW is a flexible Parallel MultiObjective Airfoil Optimization System.

### Where is the stuff related to papers? ###

 1. */probabilistic* - Torque 2018 article "Probabilistic design for wind energy applications"
 1. */dev* - Torque 2018 article "Improved airfoil polar predictions with data-driven closure relations"

### How do I get set up? ###

 1. Make sure your system meets the following requirements:  
  1.1 Recent version of Mac Os X (eg. 10.12) or Linux (eg. Ubuntu 16.04LTS)  
  1.2 Recent version of Matlab (eg. R2018a) with several toolboxes
  1.3 Access to Rfoil MK3 (preferably Nsira edition) with GNU fortran compiler ()
 2. Get the code from the repository:  
  2.1 By downloading the whole repository form the downloads section, [here](https://bitbucket.org/gaeldeoliveira/optiflow_release/downloads)  
  2.2 Better, using a mercurial versionning client, as explained [here](https://bitbucket.org/gaeldeoliveira/optiflow_release/wiki/Using%20the%20Repository)  
 3. Explore the example case files! All of them should work out of the box:  
  3.1 Open Matlab and cd to the root of the repository (always do this on a local copy!)
  3.2 Open some case file, for example *A0_start_example.m*  
  3.3 Run the file!

### Are there any examples? ###

Yes! The code is provided with six example cases *A0-F0* located in repository root to showcase some of the codes capabilities:

  A0. Baseline Parallel multiobjective optimization case  
  B0. Parallel multiobjective optimization case two multiple polars per operating point (explicit treatment of roughness sensitivity)  
  C0. Parallel multiobjective optimization case with many polars per operating point and dynamic shape deformation to optimize airfoils with **FLAPS** *(new)*   
  D0. Parallel multiobjective optimization case using a simulation_protocol_dynamizer object to optimize suction region lenght and location (Gael's MSc thesis cases)  
  E0. Parallel multiobjective optimization case using a shape_dynamizer object to restrict design space to **Prescribed Thickness** *(new)*  airfoils  
  F0. Baseline Serial multiobjective optimization case  

Deeper descriptions are provided in each case files header, comments and includes (for C).


### Contribution guidelines ###

You are encouraged to contribute to the long term success of this code, by providing feedback and contributing to new developments if you want:

* Bugs can be filed in the Issues section, [here](https://bitbucket.org/gaeldeoliveira/optiflow_release/issues?status=new&status=open)
* Requests for new features can be filed in Issues or Pull Requests section, [here](https://bitbucket.org/gaeldeoliveira/optiflow_release/pull-requests)

You have write access to this repository, and you can use it to:

 * Add new user functions to the *./user_src* folder
 * Add new case examples to the repository root *./* folder
 * It is recommended to create a new branch whenever you want to modify the main source base, located in the *./src* folder (a chat with Gael or Ricardo can also be a good idea when you go down that path!)!


### Who do I talk to? ###

The code authors/extenders/maintainers, Gael de Oliveira and Ricardo Pereira, preferably by email:
```
#!mailto:gaeldeoliveira@gmail.com

gaeldeoliveira@gmail.com
G.L.deoliveiraandrade@tudelft.nl
```
```
#!mailto:asimov1984@gmail.com
asimov1984@gmail.com
R.B.SantosPereira@tudelft.nl
```

### Licensing Terms ###

The code was initially developed by Gael de Oliveira during his TU-Delft MSc graduation at Actiflow and subsequently extended and maintained by Ricardo Pereira. The contribution of all parties to the code's success is acknowledged. 

As such, all publications involving work for which part or all of this code have been used must include Gael de Oliveira and Ricardo Pereira in the authors list (last authors is ok!) unless otherwise agreed in written form.

This code may not be shared to third parties, inside or outside of the TU-Delft, without prior consent of the authors.