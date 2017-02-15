#!/bin/bash

# anzahl punkte
N=1

points=$(python print_test_points.py N $N)

for p in $points; do
    echo $p
    python ff_one.py cache False dim 2 h 2. Nmax 1e3 x0 $p
done
