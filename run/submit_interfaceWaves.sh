#!/bin/bash

KS=(0.15 1.50 15.00)

for K in "${KS[@]}" ; do
    sh ./run/run_activeFolding.sh -K ${K}
done
