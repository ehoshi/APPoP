#!/bin/bash

rm props.prod
rm var.txt
rm unc.txt
./props_Linux < props.in > props.out 2>&1
cat props.prod | awk '/^<(PE|P|D)>/{print $3}' > var.txt
cat props.prod | awk '/^<(PE|P|D)>/{print $5}' > unc.txt
