#!/bin/bash

# this script only contains a line that looks exactly like typed on the commandline
# to run the simulation program.

# also, don't forget to put full path name to the executable unless it is located in the 
# current working directory

# NOTE: make stdin.temp, then modify it to stdin.eq
# then setup the equilibration input files
# which is mostly done. Just change the 
# output suffix .test -> .eq

# Run equilibration 

sed -e s/A_PLACE/begin/ -e s/TF_PLACE/T/ -e s/STEP_PLACE/1000/ -e s/TEMP_PLACE/298/ -e s/SUFFIX_PLACE/eq/ stdin.temp > stdin.eq
exec ./flucq_Linux < stdin.eq
