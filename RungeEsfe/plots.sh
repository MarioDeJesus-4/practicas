#!/bin/bash

# Count the number of \.dat files
ntmp=$(ls *\.dat | wc -w)
n=${ntmp}
# echo "***HOLA***" "n=$n ntmp=$ntmp"

# Use seq and parallel to process each \.dat file
seq 1 $ntmp | parallel -j 12 '
    file=$(ls *\.dat | awk "NR=={}")
    out=$(echo $file | sed "s/\.dat/\.png/")
    echo $out
    echo {}
    echo "
    set terminal png
    set title \"$file\"
    set output \"$out\"
    # plot [0:2] [-0.1:1.1] \"$file\" using 1:4 with lines
    plot [0.1:10.0] [0.0:10] \
         \"$file\" using 1:2 with lines title \"u1\", \
         \"$file\" using 1:3 with lines title \"u2\", \
         \"$file\" using 1:4 with lines title \"u3\"
    " | gnuplot
'
