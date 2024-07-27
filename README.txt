This folder contains the codes necessary to calculate the Polynomial Chaos Expansion (PCE) metamodel of the following publication:

da Silva AHA, Pohl DP, Stojadinović B. Displacement prediction equations for seismic design of single friction pendulum base-isolated structures. Earthquake Engng Struct Dyn. 2024;1-24. https://doi.org/10.1002/eqe.4201

The Python functions necessary to perform the calculations are contained in the file PCE_calculator.py.
The PCE coefficients and basis functions are contained in the subfolder CSVfiles.
There are three files with examples:
 (1) PCE_LR_COMPARISON: is a jupyter notebook that shows how to obtain the fragility curve coefficients from the PCE model for a specific set of input parameters. Next, it also calculates the results of the LR model for comparison. It plots the fragility curves obtained with both models. 
 (2) PAPER_EXAMPLE_1: shows the calculations performed in the paper mentioned above (Section 5).
 (2) DISSERTATION_EXAMPLE: shows the calculations performed in Silva's Ph.D. dissertation (Chapter 4.5).
