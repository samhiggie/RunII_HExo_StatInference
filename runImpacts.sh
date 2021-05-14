#!/bin/bash
if [ -z $1 ]
then
    mainout='AMass_gauss'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi 
#text2workspace.py datacard_full_${mainout}_mmmt.txt
#
#combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -t -1 --expectSignal=1 --doInitialFit --robustFit 1 --X-rtd MINIMIZER_analytic
#
#combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -t -1 --expectSignal=1 --doFits --robustFit 1 --X-rtd MINIMIZER_analytic
#
#combineTool.py -M Impacts -n  mass_a40_mmmt -d datacard_full_${mainout}_mmmt.root -m  40 -o testimpacts
#
#plotImpacts.py -i testimpacts -o testimpacts_12May2021

text2workspace.py datacard_full_mass_a40_${mainout}_mmmt.txt

combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -t -1 --doInitialFit  --freezeParameters MH

combineTool.py -M Impacts -m 40 -n mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -t -1 --doFits --freezeParameters MH

combineTool.py -M Impacts -n  mass_a40_mmmt -d datacard_full_mass_a40_${mainout}_mmmt.root -m  40 -o testimpacts

plotImpacts.py -i testimpacts -o testimpacts_12May2021

