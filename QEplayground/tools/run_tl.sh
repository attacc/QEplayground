#!/bin/bash

export QE=/home/attacc/SOFTWARE/qe-6.1/bin/pw.x

export ic=0
for filename in *in_TL*; do
    ic=$[$ic+1]
    $QE -inp $filename > output_TL$ic
done
