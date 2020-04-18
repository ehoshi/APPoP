# output the file name at the start of each new file
FNR == 1 { print FILENAME }

# we want timestamps that appear after the "====..." line,
/^={39}$/ { capture = 1 }
# but before an empty line
/^$/ { capture = 0 }

# also want timestamps between pairs of "++++..." lines
/^\+{49}$/ { capture = !capture }

# when capture is enabled output any line starting with a number
# (i.e., a timestamp)
capture && /^[+-]?[0-9]+/

