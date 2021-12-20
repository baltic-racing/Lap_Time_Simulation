# Lap_Time_Simulation

Laptime simulation which is capable of simulating the Acceleration, AutoX and Endurance events of the Formula Student on different tracks.
The simulation features different predefined setups and the option to simulate combustion as well as electric racecars.
The setups are fully customizable and all variables can be changed without coding using the GUIs.
A single run can be simulated as well as an one or two dimensional parameter sweep for sensitivity analysis.
The software also features a complete GUI in which all parameters of the simulation can be changed and where all data can be inspected.

The simulation features a GUI with a completely by the user configurable Color Scheme including Dark and Bright Mode.

The Results can be viewed via the built in predefined plots or with the track analyzer or the real time lap replay as well as with the Simulink Data Inspector.

Copyright © 2021 Baltic Racing developed by: 

Lead Developer:
* [**Eric Dornieden**](https://github.com/builder1one)

Team member:
* Lukas Deeken
* Danesh Umarani

# Feature Overview

The simulation features many GUIs for different tasks a small selection of the different GUI windows which are currently implemented can be found below.

## Main GUI
![Main GUI](images/GUIHowToPage.png)  

* Main Gui with all important links and a documentation how-to use the lap time simulation. 

## Setup GUI
![Simulation Setup GUI](images/GUISimulationSetupPage.png)

* GUI used to start the simulations with options to select the track, the car setup and define simulation parameters.

## Suspension GUI
![Suspension GUI](images/SuspensionGUI.png)

* The Suspension GUI allows to edit the suspension setup of the race car including the input of all kinematic points. With that data important parameters for the handling of the car like roll or pitch centers are calculated automically and the complete suspension is drawn as a dynmic 3D Plot.

## Drivatrain GUI
![Drivetrain GUI](images/DrivetrainGUI.png)

* The different setup GUIs which can be opened from the Car Setup tab GUI allow the change of the parameters used by the simulation.

## Result GUI
![Result GUI](images/GUIResultPage.png)

* Result GUI to view the results from a simulation. Allows to plot different graphs for a single run or to directly compare multiple runs. It also features buttons to open the other GUIs which can be furthermore used to analyze the data from the simulation.
