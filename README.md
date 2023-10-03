# TreadmillingEvo
Evolution algorithm and code used for Hadjivasiliou and Kruse 2023

Code was written on c++. 

The main evolution algorithm runs in main.cpp that tracks the filament size and position of all bound subunits. For a given set of parameters the algorithm computes a mean value for the filament length. the parameters are then mutated and a new mean length is computed. Then a selection step determines whether the new parameters are accepted or not. The process continues until the mean length is within a window of the target. 

Global.hpp contains definitions of contants and global parameters used in the program. 

Filaments.cpp and Filaments.hpp contain functions that update the filament state using the Gillespie scheme. 

Evolution.cpp and Evolution.hpp contains functions that impose mutation on paremeters and selection. 

Functions.cpp and Functions.hpp contain a library of functions that need to be called throughout the program. 

InitiationFunctions.cpp and InitiationFunctions.hpp contain functions that are used to initiate arrays.

Write.cpp and Write.hpp contan function that ouput the evolved parameters and filament dynamics over time. 

This program was used to perform the evolution experiments in Hadjivasiliou and Kruse 2023. Its main output is the evolution of the parameter values over evolution time. Those can be plugged into  a simpler version of this code without the evolution component to extract more infromation about each individual parameter set. 
