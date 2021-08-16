# EventShapes_HiggsTOgg_QCD_NLO
Source code for EERAD3 to compute 3-jets event shapes QCD corrections up to NLO to Higgs decays into gluons. Reference: Event_Shapes_in_hadronic_Higgs_decays.pdf (file in this repository).

The documentation of the main program EERAD3 can be found at: https://arxiv.org/abs/1402.4140. The full program is available at: https://eerad3.hepforge.org/ (as well as the NNLO QCD corrections to hadronic Z decays).

Mind that EERAD3 is meant for NNLO QCD corrections. Hence, in the full program above, more files are present with respect to the ones in this repository, giving that our present analysis extends up to NLO precision. These files have been removed in the Makefile of this repository. In the documentation of EERAD3, as far as this repository is concerned, only the general and NLO discussion are needed. In this repository, we included the program eerad3 as well as all (and only) the file needed to compute NLO corrections in QCD to 3-jets events in Higgs decays into two gluons.


To run the present version, namely to obtain NLO QCD corrections to event shape for hadronic Higgs decays to two gluons, please do the following (the lines which begin with $ are to be executed from a shell (Linux*) where all the files in the git-repository are located):

#### Download the zip repository with all the files. Go to the folder where it is stored in your laptop. Open a terminal. Extract the files with the following lines:

$ unzip EventShapes_HiggsTOgg_QCD_NLO-main.zip

#### Run the make file to compile the code. You can change the files to be compiled. For the H->gg NLO QCD corrections these are sigHG.f and aversub0.f. The compiler can be changed as well but we recommend the latest version of gfortran with the option given in the make file. In the folder where the zip repository has been extracted, open a terminal and execute:

$ make


#### Use the eerad3.input or the other input files to select the precision (number of points) and the event shape to be produced. You can change the statistical seed for the pseudo-random MonteCarlo generation (more information at https://arxiv.org/abs/1402.4140, section 4). Run the program:

$ ./eerad3 




P.S. * We are no expert of Windows. We recommend to install minGW, a Windows port of the GNU Compiler Collection (GCC).
