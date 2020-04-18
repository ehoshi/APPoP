#!/bin/bash

# this script only contains a line that looks exactly like typed on the commandline
# to run the simulation program.

# also, don't forget to put full path name to the executable unless it is located in the 
# current working directory

# Using the output from eq, setup the production input files

# Run equilibration 

sed -e s/A_PLACE/eq/ -e s/TEMP_PLACE/100/ -e s/TF_PLACE/F/ -e s/STEP_PLACE/90000000/ -e s/SUFFIX_PLACE/prod/ stdin.temp > stdin.prod
exec ./flucq_Linux < stdin.prod
