#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <boost/program_options.hpp>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "time.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAxis.h"

///////////////////////////////////////////
//
//supporting header files 
///////////////////////////////////////////

#include "outTuple.h"
#include "inTuple.h"
#include "shape.h"



int main(int ac, char* av[]) {

    //namespaces
    using namespace std;
    using namespace ROOT; 
    using namespace RooFit;
    namespace po = boost::program_options;

    // command line options? 
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "help message")
        ("quiet,q", "suppress messages to stdout")
        ("output-file,o",  po::value< std::string >(), "output file path: /path/to/file")
        ("input-file", po::value< std::string >(), "Path to data file: /path/to/data/file")
        ("output-dir,d", po::value< std::string >(), "Path to data file: /path/to/data/file")
        ("irbkgo,irbkgo", po::value< int >(), "order poly")
        ("irbkgt,irbkgt", po::value< string >(), "type")
        ("bkgo,bkgo", po::value< int >(), "order bkg poly")
        ("bkgt,bkgt", po::value< string >(), " bkg type")
        ("channel,channel", po::value< string >(), "channel")
        ;

    po::positional_options_description p;
    p.add("input-file", 1);
    po::variables_map vm;
    //p.add("output-file", 1);

    po::store(po::command_line_parser(ac, av).
            options(desc).positional(p).run(), vm);
    po::notify(vm);
    std::string inputFile;

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("input-file"))
    {
        
        inputFile = vm["input-file"].as<std::string>();
        cout << "Input files are: " 
             << inputFile << "\n";
    }
    string channel = vm["channel"].as<string>();
    //pass by reference to function ... it expects a pointer * 
    //shape * ZZ = new shape(&vm,&desc,&p,"irBkg");
    //shape * ZZ = new shape(vm,"irBkg");
    //ZZ->setsystematic("Nominal");
    //ZZ->connectFile(inputFile);
    //ZZ->fillTree(channel+"_inclusive");
    //ZZ->fillPDFs(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
    //ZZ->fillDataset();
    //ZZ->fitToData();
    //ZZ->finalFitToData();
    //ZZ->printParamsNLLs();
    //ZZ->createPlots(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
    // 
    //shape * FF = new shape(vm,"Bkg");
    //FF->setsystematic("Nominal");
    //FF->connectFile(inputFile);
    //FF->fillTree(channel+"_inclusive");
    //FF->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
    //FF->fillDataset();
    //FF->fitToData();
    //FF->finalFitToData();
    //FF->printParamsNLLs();
    //FF->createPlots(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");

    //vector<shape*> signals;
    map<string,shape*> signals;
    vector<string> masses = {
            //"a15",
            "a20",
            "a25",
            "a30",
            "a35",
            "a40",
            "a45",
            "a50",
            "a55",
            "a60"
            };
    for(auto const& x: masses){
        cout<<"working on signal "<<x<<endl;
        //signals.push_back(new shape(vm,x));
        signals[x] = new shape(vm,x);
    } 
    int mass;
    //shape * signal = new shape(vm,"a20");
    //signal->setsystematic("Nominal");
    //signal->connectFile(inputFile);
    //signal->fillTree(channel+"_inclusive");
    //signal->fillPDFs("gaussian",20,"a20");
    //signal->fillDataset();
    //signal->fitToData();
    //signal->finalFitToData("gaussian");
    //signal->printParamsNLLs();
    //signal->createPlots("gaussian",20,"a20");

    for(auto const& sh: signals){
        auto s = sh.second;
        mass = stoi((s->shape_dist).substr(1,2)); 
        s->setsystematic("Nominal");
        s->connectFile(inputFile);
        s->fillTree(channel+"_inclusive");
        s->fillPDFs("gaussian",mass,s->shape_dist);
        s->fillDataset();
        s->fitToData();
        s->finalFitToData("gaussian");
        s->printParamsNLLs();
        s->createPlots("gaussian",mass,s->shape_dist);
    }
    office * off = new office(vm); 
    off->interpolateParameters(signals);

    return 0;
}
int oldy(int ac, char* av[]){
    string in = *(av+1);  // av starts at 1 remember :)
    TFile* fileIn = TFile::Open(in.c_str(),"READ");
    std::cout<<"file in "<<fileIn->GetName()<<std::endl;

    string out = *(av+2);
    TFile* fileOut = TFile::Open(out.c_str(),"RECREATE");

    //input tree
    TTree * inTree = (TTree*) fileIn->Get("Events");
    inTuple* intuple = new inTuple(inTree);

    //outtuple
    TTree *tree = new TTree("ntuple","data");
    outTuple* otuple = new outTuple(tree);
    fileOut->cd();

    //event loop
    for (int ent = 0;ent<inTree->GetEntries();ent++){

        float idkNumber=1.0;
        if(ent%1000==0){
            //fprintf(stdout,"\r Processed events: %8d ",ent);
            std::cout<<"processed events: "<<ent<<std::endl;
        }

        otuple->fill(intuple,ent);
    }

    otuple->write(fileOut);//check header for this

    fileOut->Close();
    return 0;
}
