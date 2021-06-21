echo "attempting extrapolation"
datacardnum=20
if [ -z $1 ]
then
    mainout='test'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi 
for i in {18..62}
do
if [ $((${i} % 5)) -eq 0 ]
then
    datacardnum=${i}
fi 
echo "working on mass point "${i}
#combine -M AsymptoticLimits -m ${i} -n _mass_${i}_mmmt --run blind --freezeParameters MH datacard_full_mass_${datacardnum}_${mainout}_mmmt.txt 
#combine -M AsymptoticLimits -m ${i} -n _mass_${i}_mmmt --run blind datacard_full_mass_${datacardnum}_${mainout}_mmmt.txt 
combine -M AsymptoticLimits -m ${i} -n _mass_${i}_mmmt --X-rtd ADDNLL_RECURSIVE=0 --freezeParameters MH --run blind datacard_full_${mainout}.txt  
done 
    


echo "hadding the files"
rm higgsCombine_aa_mmmt_best.root
hadd higgsCombine_aa_mmmt_best.root higgsCombine_mass_*_mmmt.AsymptoticLimits.mH*.root
root -l -b -q 'plotLimit.C+("aa","mmmt",2)'
