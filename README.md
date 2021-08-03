# EventShapes_Higgs--gg_NLO
Source code for EERAD3 to compute 3-jets event shapes NLO QCD corrections to Higgs decays into gluons. Reference: link_to_thesis.

The documentation of the main program EERAD3 can be found at: https://arxiv.org/abs/1402.4140
The full program is available at: https://eerad3.hepforge.org/ (as well as the NNLO QCD correction to hadronic Z decays).

To run the present version, namely to obtain NLO QCD corrections to event shape for hadronic Higgs decays to two gluons, execute the following line from
a shell (Linux) where all the files in the git-repository are located:

### Run the make file to compile the code. You can change the files to be compiled. For the H->gg NLO QCD corrections these are sigHG.f and aversub0.f.
### The compiler can be changed as well but we recommend the latest version of gfortran with the option given in the make file.

$ make


### Use the eerad3.input or the other input files to select the precision (number of points) and the event shape to be produced.
### Run the program. You can change the statistical seed for the pseudo-random MonteCarlo generation. More information at https://arxiv.org/abs/1402.4140

$ ./eerad3 
