Rotor-Stator Impact Code
========================
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-purple.svg)](https://github.com/SV3A/masterproject-contact-code/blob/master/LICENSE)

A rotordynamical object-oriented program for simulating impact between a rotor and a stator.
The code was created as part of the master's thesis "[Analytical and Numerical Modeling of Contact Forces in Rotordynamics with Experimental Verification](doc/svea-abma-thesis.pdf)" conducted at the Technical University of Denmark (DTU) in 2019.

The specific rotor dynamics are implemented via the equations of motion in `dydt.m`, the currently implemented system is seen in the figure below.

![rotor model](doc/rotor-mod.png)

### Usage
- An example of a program is given in `sim1.m`
- System parameters are set in `settings.toml`

--------------------------------------------------
### Simulation Example
Below a result obtained from the code is seen:
![trajectory simulation](doc/animation.gif)
