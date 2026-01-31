# DataDrivenEndoscopyDepartmentResourcesOptimization-MILP

MILP Simulation – Execution and Configuration Guide

This repository contains the MATLAB implementation of a MILP-based scheduling simulation.
To successfully run the simulation, the required input data and configuration files must be correctly set up as described below.

1. Required Files and Directory Structure

To execute the MILP simulation, all input data must be provided in .m format and placed in the same directory as the main script:

Main script:
MDPI_prova04_AppliedSciences11.m

Input data file(s):
Provided as MATLAB .m files and loaded directly by the main script.
When running reduced test cases, the file Dades_COL_GAS_little_OK.m is used.

2. Basic Configuration Parameters

The main script allows basic parameterisation through the following variables:

little = false;              % Enables execution on a reduced dataset
num_days_desired = 1;        % Length of the scheduling horizon (in days)
default_inP = 20;            % Initial number of procedures to be packed per MILP iteration
num_rooms = 2;               % Number of available procedure rooms


Setting little = true activates a reduced dataset intended for debugging or rapid testing.

num_days_desired defines the temporal scope of the scheduling block.

default_inP controls the initial problem size per iteration.

num_rooms specifies the number of rooms available for scheduling.

3. Simulation Parameters

Additional configuration options are available in the PARAMETRES SIMULACIÓ section of the main script.

In this section, the user can modify:

Standard procedure durations

Total available time per day

Slot granularity (eslot duration)

Service start time

Other operational and scheduling parameters

These parameters allow the simulation to be adapted to different organisational or operational scenarios.

4. Output and Results

Simulation results are automatically stored in a separate directory for each MILP iteration, ensuring full traceability of intermediate and final outputs.

Each iteration folder contains the following files and summaries:

descartats
Total number of procedures that could not be scheduled by the MILP model.

descartats_filtre_data
Same as descartats, but restricted to procedures falling strictly within the exact scheduling period defined for the iteration.

mylog.txt
Complete console output generated during execution, including solver messages and diagnostic information.

reserves_vs_demanda
Analysis comparing unscheduled demand against the amount of pre-planned reserved time, highlighting capacity deficits or surpluses.

Resultat_planificacio_total.xlsx
Final aggregated scheduling output produced by the MILP, containing the complete planned solution.

resultat_planificacioXX
Partial MILP scheduling outputs generated during intermediate iterations (indexed by XX).

tauladebug.xlsx
Detailed log of input data, decisions, and key metrics recorded at each MILP iteration loop.

tauladebugMILP
Comprehensive record of MILP execution parameters and solver results, exported in:

MATLAB format (.m)

JSON (.json)

Excel (.xlsx)

This structured output layout supports reproducibility, debugging, and post-hoc analysis of both model behaviour and solver performance.

5. Notes

Ensure that all required .m data files are present before execution.

Results depend on the chosen parameter configuration and dataset size.

The reduced dataset option is recommended for validation and debugging purposes before full-scale execution.
