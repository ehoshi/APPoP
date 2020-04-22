#!/bin/bash

COORDS=coord.out
SEDIN=temp_models.dat
SEDOUT=models.dat


AWK=awk
BC=bc

SIG_PLACE() {
	$AWK "NR==1" $1
}

EPS_PLACE() {
	$AWK "NR==2" $1
}

qH_PLACE() {
	$AWK "NR==3" $1
}

qO_PLACE() {
	qH=$(qH_PLACE $1)
	echo "-2.0 * $qH" | $BC
}

substitutions=""
for var in SIG_PLACE EPS_PLACE qH_PLACE qO_PLACE; do
	replacement=$($var $COORDS)
	echo "Replacing $var with $replacement"
	substitutions="$substitutions -e s/$var/$replacement/g"
done
sed $substitutions "$SEDIN" > "$SEDOUT"
