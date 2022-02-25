#ifndef shape_h
#define shape_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <boost/program_options.hpp>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TROOT.h>

#include "time.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TString.h"
#include "TList.h"
#include "TAxis.h"
#include "TAxis.h"
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <map>
#include <utility>
#include <boost/program_options.hpp>

#include "inTuple.h"

namespace po = boost::program_options;
using namespace std;
using namespace ROOT;
using namespace RooFit;

/////////////////
//unpacking by vector 
//https://stackoverflow.com/questions/11044504/any-solution-to-unpack-a-vector-to-function-arguments-in-c
////////////////
#include <iostream>
#include <utility>
#include <vector>
#include <cassert>

namespace util {
template <typename ReturnType, typename... Args>
struct function_traits_defs {
  static constexpr size_t arity = sizeof...(Args);

  using result_type = ReturnType;

  template <size_t i>
  struct arg {
    using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
  };
};

template <typename T>
struct function_traits_impl;

template <typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(*)(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename T, typename V = void>
struct function_traits
    : function_traits_impl<T> {};

template <typename T>
struct function_traits<T, decltype((void)&T::operator())>
    : function_traits_impl<decltype(&T::operator())> {};

template <size_t... Indices>
struct indices {
  using next = indices<Indices..., sizeof...(Indices)>;
};
template <size_t N>
struct build_indices {
  using type = typename build_indices<N - 1>::type::next;
};
template <>
struct build_indices<0> {
  using type = indices<>;
};
template <size_t N>
using BuildIndices = typename build_indices<N>::type;

namespace details {
template <typename FuncType,
          typename VecType,
          size_t... I,
          typename Traits = function_traits<FuncType>,
          typename ReturnT = typename Traits::result_type>
ReturnT do_call(FuncType& func,
                VecType& args,
           indices<I...> ) {
  assert(args.size() >= Traits::arity);
  return func(args[I]...);
}
}  // namespace details

template <typename FuncType,
          typename VecType,
          typename Traits = function_traits<FuncType>,
          typename ReturnT = typename Traits::result_type>
ReturnT unpack_caller(FuncType& func,
                VecType& args) {
  return details::do_call(func, args, BuildIndices<Traits::arity>());
}
}  // namespace util
TPaveText add_lumi(int intyear,bool doRatio){
    string year;
    year = to_string(intyear);
    float lowX=0.65;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.60;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC");
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.SetTextSize(0.04);
    lumi.SetTextFont (   42 );
    if (year == "2016")
        lumi.AddText(((year)+" 35.9 fb^{-1} (13 TeV)").c_str());
    if (year == "2017")
        lumi.AddText(((year)+" 41.8 fb^{-1} (13 TeV)").c_str());
    if (year == "2018")
        lumi.AddText(((year)+" 59.7 fb^{-1} (13 TeV)").c_str());
    return lumi;
}

TPaveText add_CMS(bool doRatio){
    float lowX=0.17;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.17;
        lowY=0.835;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC");
    lumi.SetTextFont(61);
    lumi.SetTextSize(0.06);
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.AddText("CMS");
    return lumi;
}

TPaveText add_Preliminary(string channel, bool doRatio){
    float lowX=0.45;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.30;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC");
    lumi.SetTextFont(52);
    lumi.SetTextSize(0.04);
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.AddText(("Preliminary "+channel).c_str());
    return lumi;
}


class shape{
    public:
        po::variables_map shape_vm;
        string           output="default";//string for naming
        string           outputdir="default_folder";
        string           channel = "mmmt";
        //po::positional_options_description *shape_p;
        //po::options_description *shape_desc;
        Int_t            treeNum;
        TTree            *ttree;
        string           shape_name;
        string           shape_dist;
        string           type="shapetype";//RooFit Pdf class
        string           systematic="Nominal";//systematic for shapes

        TFile           *tfile;
        RooRealVar *poi;
        RooRealVar *norm;
        RooRealVar *eventweight;
        map<string,RooAbsPdf*> pdf;
        map<string,float> nll;
        map<string,float> deltanll; // dictionary for the change in nlls from the pdf fits
        map<string,TList*> coeffs; // dictionary for RooRealVars
        map<string,TList*> finalcoeffs; // dictionary for RooRealVars
        RooDataSet *data;
        map<string,RooFitResult*> fitresults; // dictionary to hold fit results from the PDFs;
        map<string,RooFitResult*> finalfitresults;// dictionary to hold fit results from the PDFs;
        map<string,RooPlot*> plotframe;
        map<string,TCanvas*> canvas;
        float rangeMin  = 0.0;
        float rangeMax  = 0.0;
        float shape_binning  = 16;

        //current variables - make this into a class later
        //datatype
        unsigned int     value;

        //branch
        TBranch     b_value;

        //new variables
        float  newvar;

        //new branch
        TBranch b_newvar;

        //////////////////////////////////////
        // function methods and inheritance
        //////////////////////////////////////

        shape(  po::variables_map vm,
                //po::options_description *desc,
                //po::positional_options_description *p,
                string dist);
        ~shape();

        void Init(TTree *ttree);
        void connectFile(string path);
        void setsystematic(string sys);
        void fillTree(string pathtotree);
        void fillPDFs(string type,int order,string dist);
        void fillDataset();
        void fitToData();
        void finalFitToData(string type);
        void printParamsNLLs();
        void createPlots(string type,int order,string dist);
        Int_t GetEntry(int entry);
        Long64_t GetEntries();
        TTree* GetTree();
        void write(TFile* tfile);
        void fill(inTuple* iTuple, int entry);
        //other containers
};

//constructors
shape::shape(po::variables_map vm,
    //po::options_description *desc,
    //po::positional_options_description *p,
    string dist){
    //setting the program options 

    shape_vm = vm;
    //shape_p = p;
    //shape_desc = desc;

    //arglist = vm;
    shape_dist = dist;
    cout<<"creating shape with distribution "<<dist<<endl;

    output = shape_vm["output-file"].as<std::string>();
    outputdir = shape_vm["output-dir"].as<std::string>();
    channel = shape_vm["channel"].as<std::string>();
    cout<<"creating directory "<<outputdir<<endl;
    //output = arglist["output-file"].as<std::string>();
    //creating a directory
    const int dir_err = mkdir(outputdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        printf("Error creating directory!\n");
        //exit(1);
    }
  
  
    system("dir");

    //setting the vars
    poi = new RooRealVar("mll","m_{#mu#mu}", 18, 62);
    norm = new RooRealVar("norm",   "Normalization",1.0,0.0,10);
    eventweight = new RooRealVar("weight",   "weight",0.0,5.0);

}
shape::~shape(){}

// Functions
void shape::connectFile(string path)
{

    std::cout<<"trying to open "<<path.c_str()<<std::endl;
    //unique::ptr<TFile> tfile(TFile::Open(path.c_str()));
    //tfile->TFile::Open(path.c_str(),"read");
    tfile = new TFile(path.c_str(),"READ");
    tfile->cd();
    tfile->ls();

    return;
}
void shape::fillTree(string pathtotree)
{
    cout<<"trying to get "<<(systematic+"_"+shape_dist)<<endl;
    ttree = (TTree*)tfile->Get((pathtotree+"/"+systematic+"_"+shape_dist).c_str());
    cout<<"ttree entries "<<ttree->GetEntries()<<endl;
    return;
}
void shape::setsystematic(string sys)
{
    systematic = sys;
    return;
}
void shape::fillPDFs(string type,int order,string dist)
{
     
    shape_name = systematic+"_"+shape_dist+"_"+output;
    cout<< "shapename "<<shape_name<<endl;
    //RooRealVar * v = new RooRealVar("var","var",1.0, -10000.0, 10000.0);
    TList * tempvector_coeffs = new TList(); 
    TList * tempvector_finalcoeffs = new TList(); 

    if(type.compare("bernstein")==0){
        rangeMin = 18.0;
        rangeMax = 62.0;
        //vector<RooRealVar> tempvector_coeffs; 
        //vector<RooFormulaVar> tempvector_finalcoeffs; 

        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_coeffs->Add(
                    new RooRealVar(
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                                1.0, -10000.0, 10000.0)
                    );
        }
        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_finalcoeffs->Add(
                new RooFormulaVar(
                (TString)("c"+to_string(cnum)+"_sq_"+shape_dist+"_"+systematic+"_"+output),
                "@0*@1",
                RooArgList(*((RooRealVar*)tempvector_coeffs->At(cnum)),*((RooRealVar*)tempvector_coeffs->At(cnum))))
                );
            } 

        //std::pair<std::map<char,int>::iterator,bool> ret;
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        finalcoeffs[type+to_string(order)] = tempvector_finalcoeffs;//empty vector of type RooRealVar 
        
        pdf[type+to_string(order)] =  new RooBernstein(
                                        (systematic+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        (systematic+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        *poi,
                                        //util::unpack_caller(
                                        //RooArgList,tempvector_finalcoeffs.at(type+to_string(order)))
                                        RooArgList(*tempvector_finalcoeffs)
                    );

        //pdf.insert(
        //    pair<string,RooBernstein> (
        //        type+to_string(order),
        //        new RooBernstein(
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            &poi,
        //            int an = util::unpack_caller(
        //                RooArgList,tempvector_finalcoeffs.at(type+to_string(order))
        //                )
        //            )
        //        )
        //    );
        }

    if(type.compare("gaussian")==0){
        cout<<"mass num? "<<shape_dist.substr(1,2)<<endl;
        float mass = stoi((shape_dist.substr(1,2)));
        rangeMin = mass - 2.0;
        rangeMax = mass + 2.0;
        shape_binning = 25;
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                            stof(shape_dist.substr(1,2)),
                            stof(shape_dist.substr(1,2))-2.0,
                            stof(shape_dist.substr(1,2))+2.0)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                            0.5,
                            0.0,
                            2.0)
                            //0.3,
                            //0.1,
                            //0.5)
                );
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        pdf[type+to_string(order)] = new RooGaussian(
                                            (TString)("sigfit"),
                                            (TString)("sigfit"),
                                            *poi,
                                            *((RooRealVar*) tempvector_coeffs->At(0)),
                                            *((RooRealVar*) tempvector_coeffs->At(1))
                                            );

        //pdf->insert(
        //    pair<string,RooGaussian> (
        //        (type+to_string(order)),
        //        new RooGaussian(
        //            (TString)("sigfit"),
        //            (TString)("sigfit"),
        //            &poi,
        //            tempvector_coeffs[0],
        //            tempvector_coeffs[1]
        //            )
        //        )
        //    );


    }
    return;
}
void shape::fillDataset()
{
    //if (!ttree) return;
    eventweight = new RooRealVar("finalweight","finalweight",-10000.0,10000.0);
    //data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),RooArgSet(*poi),Import(*ttree));
    data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),
                        RooArgSet(*poi,*eventweight),
                        Import(*ttree),WeightVar("finalweight")
                    );
    cout<<"filled dataset dataset"<<endl;
    norm = new RooRealVar((shape_dist+"_norm").c_str(),(shape_dist+"_norm").c_str(),data->sumEntries(),0.0,10*data->sumEntries());
    
    return;
}
void shape::fitToData()
{
    for(auto const& x: pdf){
        cout<<"working on pdf "<<(x.first).c_str()<<endl;
        canvas[x.first] = new TCanvas(("canvas_"+shape_name).c_str(),("canvas_"+shape_name).c_str(),600,600);
        poi->setRange(rangeMin,rangeMax);
        plotframe[x.first] = poi->frame();
        plotframe[x.first]->SetAxisRange(rangeMin,rangeMax);
        data->plotOn(plotframe[x.first],Range(rangeMin,rangeMax),Binning(shape_binning));
        fitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        //fitresults[x.first]->Print();
    }
    

    return;
}
void shape::finalFitToData(string type)
{
    for(auto const& x: pdf){
        cout<<"working on pdf final fit"<<(x.first).c_str()<<endl;
        plotframe[x.first] = poi->frame();
        data->plotOn(plotframe[x.first],Binning(shape_binning));
        //TIter next(coeffs[x.first]); 
        //TObject * obj = 0;
        //while(obj = next()){
        if(type.compare("gaussian")!=0){
        for(const auto&& obj: *coeffs[x.first]){
                cout<<"object name "<<obj->GetName()<<endl; 
                cout<<"object class name "<<obj->ClassName()<<endl; 
                obj->Print();
                //lesson ... learned ...type  cast in paranthesis!!
                string varname = obj->GetName();
                cout<<"variable name "<<varname<<endl;
                    ((RooRealVar*) obj)->setRange( 
                        ((RooRealVar*)obj)->getValV()/1.50,
                        ((RooRealVar*)obj)->getValV()*1.50
                        ); 
                    }
            }

        if(type.compare("gaussian")==0){
        cout<<"running gaussian final fit"<<endl;
        for(const auto&& obj: *coeffs[x.first]){
            cout<<"object name "<<obj->GetName()<<endl; 
            cout<<"object class name "<<obj->ClassName()<<endl; 
            obj->Print();
            //lesson ... learned ...type  cast in paranthesis!!
            string varname = obj->GetName();
            cout<<"variable name "<<varname<<endl;
                if(varname.find("sigma") != string::npos)  {
                    cout<<"setting variable ranges for sigma"<<endl;
                    ((RooRealVar*) obj)->setRange( 
                        ((RooRealVar*)obj)->getValV()/2.0,
                        ((RooRealVar*)obj)->getValV()*2.0
                        ); 
                    }
                }
            }
        finalfitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        //finalfitresults[x.first]->Print();
    }
    

    return;
}
void shape::printParamsNLLs()
{
    ofstream nllresults;
    nllresults.open(outputdir+"/nllresults"+shape_dist+"_"+shape_name+"_"+systematic+".txt");
    nllresults << "fit results for "<<shape_name<<endl; 
    //for(auto& [key,value]: coeffs){ this is the future .... c++17
    for(auto const& x: coeffs){
        nllresults<<"finalfit values for the coeff of "<<x.first<<endl;
        //for(RooRealVar* var: *value){
        for(const auto&& obj: *coeffs[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            cout<<var->GetName()<<"  "<<to_string(var->getValV())<<endl; 
            nllresults<<var->GetName()<<"  "<<to_string(var->getValV())<<endl; 
        }
    } 
    vector<float> nllvals;
    for(auto const& x: finalfitresults){
        nllresults<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
        cout<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
        nllvals.push_back(x.second->minNll());
        
    }
    int nllsize = (int) nllvals.size();
    for (int i=0; i<(nllsize-1);i++){
        cout<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
        nllresults<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
    }

    return;
}
void shape::createPlots(string type,int order,string dist)
{
    gPad->SetLeftMargin(0.15);
    string key = type+to_string(order);//key for the map
    plotframe[key]->SetAxisRange(rangeMin,rangeMax);
    plotframe[key]->SetTitle(("#mu#mu for "+dist+" "+output).c_str());
    plotframe[key]->GetYaxis()->SetTitleOffset(1.6);
    plotframe[key]->GetXaxis()->SetTitle("M_{#mu#mu}");
    TGaxis::SetMaxDigits(2);
    data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
    
    //pdf[key]->plotOn(plotframe[key],
    //        VisualizeError(*finalfitresults[key],1.0,kFALSE),
    //        Range(rangeMin,rangeMax),
    //        //DrawOption("F"),
    //        FillColor(kOrange)
    //        );
    if(type.compare("gaussian")==0){
    pdf[key]->plotOn(plotframe[key],
            VisualizeError(*fitresults[key],0.4,kFALSE),
            Range(rangeMin,rangeMax),
            //DrawOption("F"),
            FillColor(kOrange)
            );
    }
    else{
    pdf[key]->plotOn(plotframe[key],
            VisualizeError(*finalfitresults[key],1.0,kFALSE),
            Range(rangeMin,rangeMax),
            //DrawOption("F"),
            FillColor(kOrange)
            );
    }
    pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
    data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
    pdf[key]->paramOn(plotframe[key],data);
    plotframe[key]->getAttText()->SetTextSize(0.02);
    plotframe[key]->Draw("same");
    canvas[key]->SaveAs((outputdir+"/"+to_string(order)+"_"+shape_name+"_"+systematic+".png").c_str()); 
    return;
}
Int_t shape::GetEntry(int entry)
{
    return ttree->GetEntry(entry);
}

Long64_t shape::GetEntries()
{
    return ttree->GetEntries();
}

TTree* shape::GetTree()
{
    return ttree;
}

void shape::write(TFile* tfile){
    tfile->cd();
   ttree->Write();

}
void shape::fill(inTuple* iTuple, int entry){
    iTuple->GetEntry(entry);
    value = iTuple->event;
    std::cout<<" event number "<<iTuple->event<<std::endl;
    newvar =  3.1415926534 * (value);

    //newvar = value;

    ttree->Fill();

}

class office{
    public:
    //Data members
    po::variables_map shape_vm;
    string           output="default";//string for naming
    string           outputdir="default_folder";
    string           channel="mmmt";
    int           year=2017;
    RooWorkspace * wsp; 
    vector<string> systematics = {
        "scale_e","scale_m_etalt1p2",
                       "scale_m_eta1p2to2p1","scale_m_etagt2p1",
                       "scale_t_1prong","scale_t_1prong1pizero",
                       "scale_t_3prong","scale_t_3prong1pizero"
    };
    //vector<shape*> signalshapes;


    //functions
    office(po::variables_map shape_vm);
    ~office();

    //void loadSignalShapes(vector<shape*> inputshapes);
    void interpolateParameters(map<string,shape*>& inputshapes);


};
office::office(po::variables_map shape_vm){
    output = shape_vm["output-file"].as<std::string>();
    outputdir = shape_vm["output-dir"].as<std::string>();
    channel = shape_vm["channel"].as<std::string>();
    year = stoi(output.substr(0,4));
    
}
office::~office(){}

void office::interpolateParameters(map<string,shape*>& inputshapes){
    cout<<"interpolating!"<<endl;
    TF1 * meanfit = new TF1("meanfit","pol1",16,66);
    TF1 * normfit = new TF1("normfit","pol3",16,66);
    TF1 * sigmafit = new TF1("sigmafit","pol3",16,66);
    TGraphErrors * meangraph = new TGraphErrors();
    TGraphErrors * normgraph = new TGraphErrors();
    TGraphErrors * sigmagraph = new TGraphErrors();
    int c = 0;
    //shape * sp = (shape*) inputshapes["a20"];
    cout<<"shape exist? name: "<<inputshapes["a20"]<<endl;
    cout<<"shape exist? name: "<<inputshapes["a20"]->shape_name<<endl;
    for(auto const& x: inputshapes){
        string mass = x.first.substr(1,2); 
        shape * sp = x.second;
        //cout<<"setting points "<<x.first<<endl; 
        //cout<<"point "<<to_string(c)<<endl;
        //cout<<"Mass "<<(mass)<<endl;
        //cout<<"shape name "<<sp->shape_name<<endl;
        //cout<<"var name "<<sp->coeffs["gaussian"+mass]->At(0)->GetName()<<endl;
        //cout<<"value "<<((RooRealVar*)sp->coeffs["gaussian"+mass]->At(0))->getVal()<<endl;

        meangraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["gaussian"+mass]->At(0))->getVal()
        );
        sigmagraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["gaussian"+mass]->At(1))->getVal()
        );
        normgraph->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()
        );

        c++;
    }
    TCanvas * can = new TCanvas("splinefits","splinefits",600,600);
    can->cd();
    gStyle->SetOptStat(1); 
    gStyle->SetOptFit(1); 
    gStyle->SetStatX(0.6);
    gStyle->SetStatY(0.9);

    bool doRatio = 0;
    TPaveText lumi=add_lumi(year,doRatio);
    TPaveText cms=add_CMS(doRatio);
    TPaveText pre=add_Preliminary(channel, doRatio);

    meangraph->SetName("mean");
    meangraph->SetTitle("");
    meangraph->GetXaxis()->SetTitle("Mass");
    meangraph->GetYaxis()->SetTitle("Mean Fit Parameter");
    meangraph->SetMarkerStyle(8);
    meangraph->Draw("AP");
    meangraph->Fit(meanfit);
    meanfit->Draw("same");
    lumi.Draw();
    cms.Draw();
    pre.Draw();


    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+output+".png"));
    can->Clear();

    normgraph->SetName("norm");
    normgraph->SetTitle("");
    normgraph->GetXaxis()->SetTitle("Mass");
    normgraph->GetYaxis()->SetTitle("Norm Fit Parameter");
    normgraph->SetMarkerStyle(8);
    normgraph->Draw("AP");
    normgraph->Fit(normfit);
    normfit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+output+".png"));
    can->Clear();

    sigmagraph->SetName("sigma");
    sigmagraph->SetTitle("");
    sigmagraph->GetXaxis()->SetTitle("Mass");
    sigmagraph->GetYaxis()->SetTitle("Sigma Fit Parameter");
    sigmagraph->SetMarkerStyle(8);
    sigmagraph->Draw("AP");
    sigmagraph->Fit(sigmafit);
    sigmafit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+output+".png"));
    can->Clear();
    cout<<"the year is "<<year<<endl;
    
    //doing the interpolation to a spline-like function to save to RooWorkspace 
    //signaltemplates = {};
    //signalnorms = {};
    //signaldatasets = {};
    //x = {};
    //m = {};
    //s = {};
    //MH = RooRealVar("MH","MH", 18, 63);
    //Mll = RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0);
    //MHerr = RooRealVar("MHerr","MHerr", 0, -4 , 4);
    //RooRealVar mean_c0 = meanfit->GetParameter(0);
    //RooRealVar mean_c1 = meanfit->GetParameter(1);
    //RooRealVar mean_c2 = meanfit->GetParameter(2);
    //RooRealVar mean_c3 = meanfit->GetParameter(3);
    //print(" mean formula ({0:f}+ {1:f}*@0)".format(mean_c0,mean_c1));
    //intMean = RooFormulaVar("intMean","intMean","({0:f}+ {1:f}*@0)".format(mean_c0,mean_c1),RooArgList(MH));
    //sigma_c0 = sigmafit.GetParameter(0);
    //sigma_c1 = sigmafit.GetParameter(1);
    //sigma_c2 = sigmafit.GetParameter(2);
    //sigma_c3 = sigmafit.GetParameter(3);
    //intSigma = RooFormulaVar("intSigma","intSigma","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(sigma_c0,sigma_c1,sigma_c2,sigma_c3),RooArgList(MH));
    //norm_c0 = normfit.GetParameter(0);
    //norm_c1 = normfit.GetParameter(1);
    //norm_c2 = normfit.GetParameter(2);
    //norm_c3 = normfit.GetParameter(3);
    //intNorm = RooFormulaVar("signal_norm","signal_norm","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(norm_c0,norm_c1,norm_c2,norm_c3),RooArgList(MH));
    //print("interpolated Mean formula ",intMean.Print());
    //intSignalTemplate = RooGaussian("signal"+"_"+str(args.channel),   "signal"+"_"+str(args.channel),Mll, intMean, intSigma );

    //for mass in range(16,66):;
    //    massEval = meanfit.Eval(mass);
    //    normEval = normfit.Eval(mass);
    //    sigmaEval = sigmafit.Eval(mass);
    //    signalnorms[str(mass)] = normEval;
    //    print("evaluation of mean at ",mass," is ",massEval, " generating signal template ");
    //    print("evaluation of norm at ",mass," is ",normEval);
    //    print("evaluation of sigma at ",mass," is ",sigmaEval);
    //    x[str(mass)] = RooRealVar("mll",    "mll",massEval-2.0,massEval+2.0);
    //    m[str(mass)] =RooRealVar("mean",    "mean",massEval,"GeV");
    //    s[str(mass)] = RooRealVar("sigma",    "sigma", sigmaEval,"GeV");


    return;
}






#endif
