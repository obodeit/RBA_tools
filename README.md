# RBA_tools

This tool serves as an interface to resource-allocation modelling with the RBA-method (https://rba.inrae.fr).
It includes methods to solve, manipulate and analyse RBA-models and to export model information and simulation results into various formats.
A requirement for the usage of this tool, is to install the rba library from this repository:
https://github.com/SysBioInra/RBApy

## Repository structure
### rbatools
This is the library with all the necessary classes for the usage of RBA_tools.
### html_documentation
This is a directory with html-files, serving as documentation of the rbatools library.
It can be viewed by navigating into this folder and opening the index.html file in a browser (via double click).
### Example_workflow_rba_tools.ipynb
This is a Jupyter notebook with example applications and workflows of rbatools.

## Model availability
RBA models can be obtained from the repository https://github.com/SysBioInra/Bacterial-RBA-models
The Example_workflow_rba_tools notebook requires the previosly mentioned Bacterial-RBA-models repository to be placed in the same location, as this repository.
