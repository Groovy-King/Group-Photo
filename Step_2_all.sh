#!/bin/bash

for i in {1..9}; do
    for j in {1..3}; do
        sbatch Step_2.sl $i $j
    done
done
