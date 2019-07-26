# psi_transmission

This repository contains scripts, Jupyter notebooks, and other documents used for analysis of the UCN guide transmission experiment, carried out in December of 2017 at the [Paul Scherrer Institute](https://www.psi.ch/en). This file provides an overview of the repository contents.

## data 

Notably missing from the repository are the experimental data, which can be found on the main UCN cluster at the following path: **/ucn/orithyia_data/psi_transmission/name-of-data-folder**

There are two such data folders:

1. **data_ucn** - UCN count data as collected by a main detector and a monitor detector
2. **data_p_beam** - Proton beam current measurements collected throughout the experiment.
   
There is a separate README.md file within each of these data sub-directories, describing the respective contents. 

## elog

There are .html, .txt, and .pdf versions of the elog available in the **elog** sub-folder.

## notebooks and scripts

In the **notebooks** and **scripts** folders, there are tools and demonstrations of the pre-storage lifetime and transmission analyses performed. Refer to the self-contained documentation for details.

## documents and img

Some relevant articles, and reports are provided. 

## psi_transmission_report

A git submodule linking to the repository containing a report on the analysis performed. 

## to-do

- re-walk through all of the error analysis
- confirm solid works model, put necessary simulation stuff here, .gitignore the rest with path instructions
- begin simulations
- histogram of proton beam
- better error presentation of the sD2 normalization slopes ('bands? - RP')
- a commented version of the elog pdf, clarifying any ambiguous statements
- play with changing the t=0 time and the consequences for the sD2 normalization analysis
- proper commenting of a header docstring on functions.py. Should list all functions.
- monitor detector? Some description of what we found there must be put in a notebook.