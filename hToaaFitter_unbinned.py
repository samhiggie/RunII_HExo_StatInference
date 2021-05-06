#########################
#Author: Sam Higginbotham
'''

* File Name : hToaaFitter_unbinned.py

* Purpose : Fit the dimuon mass spectra of the pseudoscalar higgs candidate ... local version for now

* Creation Date : 23-10-2020

* Last Modified :

'''
#########################
import ROOT
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
#parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
parser.add_argument("-o",  "--output", default="",  help="postfix string")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
args = parser.parse_args()

#fIn = ROOT.TFile.Open("ggTo2mu2tau_40_2016.root","open")
#fIn = ROOT.TFile.Open("skimmed_mmmt.root,"open")
#filesloc = "histograms/"
#filesloc = "histograms_array_nw/"
filesloc = "histograms_2016_fits/"
#fIn2 = ROOT.TFile.Open("skimmed_2016_prompt_mmmt.root","READ")

fIn2 = ROOT.TFile.Open(args.input,"READ")

filesIn = {}


fIn2.cd()
hSignals = {}


filesIn["a40"] = [ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit0.root","open"),38.0,42.0,40.0,ROOT.kMagenta]

#getting the TTrees
datatree = fIn2.Get("data_obs")
sigtree = fIn2.Get("a40")
#hInSig = fIn2.Get("a40")   # signal distribution included above!
#tlist = ROOT.TList()
#tlist.Add(fIn2.Get("Bkg"))
##hInBkg.Add(fIn2.Get("vbf"))
#tlist.Add(fIn2.Get("irBkg"))
#tlist.Add(fIn2.Get("TrialphaBkg"))
#tlist.Add(fIn2.Get("rare"))

#merging the ttrees to make a single background TTree
#bkgtree = ROOT.TTree.MergeTrees(tlist)
bkgtree = fIn2.Get("Bkg")
bkgtree.SetName("bkg")

FFtree = fIn2.Get("Bkg")
ZZtree = fIn2.Get("irBkg")

fitParams = {}
dataHists = {}
fitModels = {}
pdfs = {}

for file in filesIn.keys():
    fitParams[file] = [
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", filesIn[file][1], filesIn[file][2]),#works for fine binning
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", 15.0, 60.0),
        ROOT.RooRealVar("mll",    "m_{#mu #mu}", 38.0, 42.0),
        ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    filesIn[file][3], filesIn[file][1], filesIn[file][2], "GeV"),
        ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#mu #mu}",1.0, 0.0,  10.0, "GeV"),
        ROOT.RooRealVar("lAlpha_"+str(file),   "#alpha of lorentz profile",     1.0, 0.0,      40.0)
        ]

for file in filesIn.keys():
    #dataHists[file] = [ROOT.RooDataHist("dh_"+str(file),"signal histo ",ROOT.RooArgList(fitParams[file][0]), hSignals[file])]
    fitModels[file] = [
        #ROOT.RooVoigtian(AMass,   "first voigtian PDF", MH, mean, alpha, sigma),
        ROOT.RooVoigtian("voigtian_"+str(file),   "first voigtian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][3], fitParams[file][2]),
        ROOT.RooRealVar("signalEvents_"+str(file), "",  sigtree.GetEntries(),   0.0, 1000.0),
        ]

#needed for binned fit ...
#for file in filesIn.keys():
#    pdfs[file] = ROOT.RooAddPdf("model_"+str(file), "", ROOT.RooArgList(fitModels[file][0]),ROOT.RooArgList(fitModels[file][1]))

#adding background and data RooDatasets
# should I include weights or not
finalweight = ROOT.RooRealVar("finalweight",   "finalweight",0.0,3.0) # is this the weight in the ttree or just a floating variable in the fit?

######################################################################################################
''' Obtaining Data
'''
######################################################################################################
#data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(datatree))
Mmm = ROOT.RooRealVar("mll","m_{#mu#mu}", 16, 66)
data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(datatree))
#bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(bkgtree))
#bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
#sig = ROOT.RooDataSet("sig","sig ",ROOT.RooArgSet(fitParams["a40"][0],finalweight), ROOT.RooFit.Import(sigtree),ROOT.RooFit.WeightVar("finalweight"))
varargset=ROOT.RooArgSet(fitParams["a40"][0],finalweight)
#tmpdst = ROOT.RooDataSet("tmpdataset","", varargset)
#tmpdst = ROOT.RooDataSet("a40","", varargset)
#tmpdst.read("skimmed_2016_prompt_mmmt.root",ROOT.RooArgSet(Mmm,finalweight))
#sig = ROOT.RooDataSet("sig","sig ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(sigtree))
#sig = ROOT.RooDataSet("sig","sig", varargset, ROOT.RooFit.Import(tmpdst),ROOT.RooFit.WeightVar(finalweight))
sig = ROOT.RooDataSet("sig","sig", varargset, ROOT.RooFit.Import(sigtree),ROOT.RooFit.WeightVar(finalweight))

FF = ROOT.RooDataSet("FF","FF ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(FFtree))
ZZ = ROOT.RooDataSet("ZZ","ZZ ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(ZZtree))

#setting up background model - normalization needed for binned models
#norm_bkg = ROOT.RooRealVar("bkg_norm",   "background Normalization",bkg.sumEntries(),0.0,bkg.sumEntries()*2)
#norm_sig = ROOT.RooRealVar("sig_norm",   "signal Normalization",sig.sumEntries(),0.0,sig.sumEntries()*2)


######################################################################################################
'''FF Fit
'''
######################################################################################################
#norm_bkg = ROOT.RooAbsReal("bkg_norm",   "background Normalization",0.0,10000.0)
#x_bkg = ROOT.RooRealVar("MH","m_{#mu#mu}", filesIn["a40"][1], filesIn["a40"][2]) # is this min and max?
c0_bkg = ROOT.RooRealVar("c0_bkg",   "coeff. of bernstein 0",0.1,-10.0,10.0)
c1_bkg = ROOT.RooRealVar("c1_bkg",  "coeff. of bernstein 1",5.0,-20.0,20.0)
c2_bkg = ROOT.RooRealVar("c2_bkg",   "coeff. of bernstein 2",5.0,-20.0,20.0)
c3_bkg = ROOT.RooRealVar("c3_bkg",   "coeff. of bernstein 3",5.0,-20.0,20.0)
c4_bkg = ROOT.RooRealVar("c4_bkg",   "coeff. of bernstein 4",5.0,-20.0,20.0)

c0_bkg_sq = ROOT.RooFormulaVar("c0_bkg_sq","@0*@1",ROOT.RooArgList(c0_bkg,c0_bkg))
c1_bkg_sq = ROOT.RooFormulaVar("c1_bkg_sq","@0*@1",ROOT.RooArgList(c1_bkg,c1_bkg))
c2_bkg_sq = ROOT.RooFormulaVar("c2_bkg_sq","@0*@1",ROOT.RooArgList(c2_bkg,c2_bkg))
c3_bkg_sq = ROOT.RooFormulaVar("c3_bkg_sq","@0*@1",ROOT.RooArgList(c3_bkg,c3_bkg))
c4_bkg_sq = ROOT.RooFormulaVar("c4_bkg_sq","@0*@1",ROOT.RooArgList(c4_bkg,c4_bkg))
#bkgHists = ROOT.RooDataHist("bkg","background histo ",ROOT.RooArgList(fitParams["a40"][0]), hInBkg)


#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",x_bkg,ROOT.RooArgList(c0_bkg,c1_bkg,c2_bkg)) #fully 2nd order polynominal
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",x_bkg,ROOT.RooArgList(c0_sq)) #constant only
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",fitParams["a40"][0],ROOT.RooArgList(c0_sq,c1_sq,c2_sq,c3_sq)) #constant only
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",fitParams["a40"][0],ROOT.RooArgList(c0_bkg_sq)) #constant only

FFfit = ROOT.RooBernstein("FFfit","FFfit",Mmm,ROOT.RooArgList(c0_bkg_sq,c1_bkg_sq,c2_bkg_sq,c3_bkg_sq)) #constant only

######################################################################################################
'''ZZ Fit
'''
######################################################################################################
#for ZZ
c0_ZZ = ROOT.RooRealVar("c0_ZZ",   "coeff. of bernstein 0",10.0,-1010.0,1010.0)
c1_ZZ = ROOT.RooRealVar("c1_ZZ",  "coeff. of bernstein 1",10.0,-1010.0,1010.0)
c2_ZZ = ROOT.RooRealVar("c2_ZZ",   "coeff. of bernstein 2",10.0,-1010.0,1010.0)
c3_ZZ = ROOT.RooRealVar("c3_ZZ",   "coeff. of bernstein 3",10.0,-1010.0,1010.0)
c4_ZZ = ROOT.RooRealVar("c4_ZZ",   "coeff. of bernstein 4",10.0,-1010.0,1010.0)

c0_ZZ_sq = ROOT.RooFormulaVar("c0_ZZ_sq","@0*@1",ROOT.RooArgList(c0_ZZ,c0_ZZ))
c1_ZZ_sq = ROOT.RooFormulaVar("c1_ZZ_sq","@0*@1",ROOT.RooArgList(c1_ZZ,c1_ZZ))
c2_ZZ_sq = ROOT.RooFormulaVar("c2_ZZ_sq","@0*@1",ROOT.RooArgList(c2_ZZ,c2_ZZ))
c3_ZZ_sq = ROOT.RooFormulaVar("c3_ZZ_sq","@0*@1",ROOT.RooArgList(c3_ZZ,c3_ZZ))
c4_ZZ_sq = ROOT.RooFormulaVar("c4_ZZ_sq","@0*@1",ROOT.RooArgList(c4_ZZ,c4_ZZ))

ZZfit = ROOT.RooBernstein("ZZfit","ZZfit",Mmm,ROOT.RooArgList(c0_ZZ_sq,c1_ZZ_sq,c2_ZZ_sq,c3_ZZ_sq,c4_ZZ_sq)) #constant only
######################################################################################################
'''Signal Fit
'''
######################################################################################################


sigfit = ROOT.RooVoigtian("sigfit",   "sigfit",
                    fitParams["a40"][0], #mll
                    fitParams["a40"][1], #mean
                    fitParams["a40"][3], #alpha
                    fitParams["a40"][2]) #sigma

# these may need to be RooABSPDFs ... not extended? for the unbinned fit...
#bkgfitModel = ROOT.RooExtendPdf("bkg","bkg",bkgfit, norm_bkg)
#sigfitModel = ROOT.RooExtendPdf("sig","sig",sigfit, norm_sig)
#making overall signal and bkg in single category
######################################################################################################
''' FF Plotting and Fitting Area
'''
######################################################################################################
overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0)
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()
#massFrame = fitParams["a40"][0].frame()
c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)
massFrame.Draw()
massFrame = overmass.frame()

#plotting data points
FF.plotOn(massFrame,ROOT.RooFit.Binning(16))

fitresult = FFfit.fitTo(FF,ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())

FFfit.paramOn(massFrame)
FFfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen),
             ROOT.RooFit.LineStyle(ROOT.kDashed),
             ROOT.RooFit.VisualizeError(fitresult,1,ROOT.kFALSE),
             ROOT.RooFit.FillColor(ROOT.kOrange))
FF.plotOn(massFrame,ROOT.RooFit.Binning(16))
massFrame.Draw()

c.SaveAs("DiMuonMass_FF_"+args.output+".pdf")
c.SaveAs("DiMuonMass_FF_"+args.output+".png")
c.Clear()

######################################################################################################
''' ZZ Plotting and Fitting Area
'''
######################################################################################################
overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0)
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()
#massFrame = fitParams["a40"][0].frame()
c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)
massFrame.Draw()
massFrame = overmass.frame()

#plotting data points
#bkg.plotOn(massFrame,ROOT.RooFit.Binning(16))
#FF.plotOn(massFrame,ROOT.RooFit.Binning(16))
ZZ.plotOn(massFrame,ROOT.RooFit.Binning(16))

#bkgfitModel.fitTo(bkg,ROOT.RooFit.Extended())
#bkgfitModel.paramOn(massFrame)
#bkgfit.fitTo(bkg,ROOT.RooFit.Extended())
#bkgfit.fitTo(bkg,ROOT.RooFit.Range(14,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
fitresult = ZZfit.fitTo(ZZ,ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())

ZZfit.paramOn(massFrame)
ZZfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen),
             ROOT.RooFit.LineStyle(ROOT.kDashed),
             ROOT.RooFit.VisualizeError(fitresult,1,ROOT.kFALSE),
             ROOT.RooFit.FillColor(ROOT.kOrange))
ZZ.plotOn(massFrame,ROOT.RooFit.Binning(16))
massFrame.Draw()

c.SaveAs("DiMuonMass_ZZ_"+args.output+".pdf")
c.SaveAs("DiMuonMass_ZZ_"+args.output+".png")
c.Clear()

######################################################################################################
''' Signal Plotting and Fitting Area
'''
######################################################################################################
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 11.0, 65.0)
overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()
massFrame.Draw()
massFrame = overmass.frame()

sig.plotOn(massFrame)
sigfit.fitTo(sig,ROOT.RooFit.Range(38,42), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
sigfit.paramOn(massFrame)
#cout<< rrv->getVal() <<"  +/-  "<<rrv->getError();

print  "signal fit value mean",fitParams["a40"][1].getVal()," error ",fitParams["a40"][1].getError()
print "signal fit value alpha",fitParams["a40"][3].getVal()," error ",fitParams["a40"][3].getError()
print "signal fit value sigma",fitParams["a40"][2].getVal()," error ",fitParams["a40"][2].getError()




sigfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))

massFrame.Draw()

c.SaveAs("DiMuonMass_sig_"+args.output+".pdf")
c.SaveAs("DiMuonMass_sig_"+args.output+".png")
c.Clear()

######################################################################################################
''' Saving the Workspace
'''
######################################################################################################
#import data and PDF into workspaco
workspace = ROOT.RooWorkspace("w")

getattr(workspace,'import')(data,ROOT.RooFit.RenameVariable("mll","MH"))
ZZfit.SetName("ZZ")
FFfit.SetName("FF")
sigfit.SetName("sig")
#try to rename mll to MH via
getattr(workspace,'import')(sigfit,ROOT.RooFit.RenameVariable("mll","MH"))
getattr(workspace,'import')(ZZfit,ROOT.RooFit.RenameVariable("mll","MH"))
getattr(workspace,'import')(FFfit,ROOT.RooFit.RenameVariable("mll","MH"))

workspace.Print()

workspace.writeToFile("HToAAWorkspace_"+args.output+".root")


del workspace
#sigCat[0].plotOn(massFrame[0])

##Make the datacard
outFile = open("datacard_"+args.output+".txt","w")#write mode

outFile.write("imax 1\n") #number of bins - only one category ... no control region
outFile.write("jmax 2\n") #number of processes minus 1
outFile.write("kmax *\n") #number of nuisance parameters
outFile.write("---------------\n")
#outFile.write("shapes * * HToAAWorkspace_combined.root w:$PROCESS\n")
outFile.write("shapes * bin1 HToAAWorkspace_"+args.output+".root w:$PROCESS\n")
outFile.write("---------------\n")

outFile.write("bin         bin1   \n")
#outFile.write("observation   "+str(bkgHists.sumEntries())+"\n")
outFile.write("observation   -1 \n") # for parametric fit this needs to be -1

outFile.write("------------------------------\n")
outFile.write("bin                     bin1     bin1    bin1\n")
outFile.write("process                 sig      ZZ      FF\n")
outFile.write("process                 0        1       2\n")
#outFile.write("rate                   "+str(sigtree.GetEntries())+"   "+str(bkgtree.GetEntries())+"   \n")
outFile.write("rate                   "+str(sig.sumEntries())+"   "+str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
outFile.write("------------------------------\n")
outFile.write("lumi     lnN              1.1    1.0    1.0\n")
outFile.write("g1Mean_a40  param          39.9432    0.0994924 \n") # form of shape paramters in fit include "name param mean std"
outFile.write("sigmaM_a40  param          0.2527    0.401558 \n") # form of shape paramters in fit include "name param mean std"
outFile.write("lAlpha_a40  param          0.528479    0.256726 \n") # form of shape paramters in fit include "name param mean std"

outFile.close()
