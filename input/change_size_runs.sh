#!/bin/bash


while read run; do
    scancel --name=$run
    rm -f /scratch-new/formanek/HYBRIDSIM/runs/$run/log.out
    rm -f /scratch-new/formanek/HYBRIDSIM/runs/$run/err.out
    rm -f /scratch-new/formanek/HYBRIDSIM/runs/$run/stats
    rm -f /scratch-new/formanek/HYBRIDSIM/runs/$run/configs/*
    sed -i 's/BoxX = 50/BoxX = 100/' /scratch-new/formanek/HYBRIDSIM/runs/$run/parameters.dat
    echo /scratch-new/formanek/HYBRIDSIM/runs/$run >> runs_to_submit.dat
done < runs-to-restart
