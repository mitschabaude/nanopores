#!/bin/bash

# anzahl punkte
N=2

points=$(python print_test_points.py N $N)

for p in $points; do
    echo $p
    python ff_one.py x0 $p
done
