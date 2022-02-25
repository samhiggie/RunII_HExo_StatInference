g++ limit.cc -o limit.exe `root-config --cflags --glibs` -lRooFit -lRooFitCore -lboost_program_options
##./limit.exe -o 2017_mmmt -d foo --input-file skimmed_2017_paperex_mmmt_combined.root 

#./limit.exe -o 2017_mmmt -d 2017_mmmt_plots  --irbkgo 2 --irbkgt bernstein --bkgo 2 --bkgt bernstein  --input-file skimmed_2017_paperex_mmmt_combined.root --channel mmmt
