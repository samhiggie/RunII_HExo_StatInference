#########################
#Author: Sam Higginbotham
'''

* File Name : hToaaFitter.py

* Purpose : Fit the dimuon mass spectra of the pseudoscalar higgs candidate ... local version for now

* Creation Date : 21-06-2020

* Last Modified :

'''
#########################
import ROOT
from datetime import datetime

#fIn = ROOT.TFile.Open("ggTo2mu2tau_40_2016.root","open")
fIn = ROOT.TFile.Open("ggTo2mu2tau_40_2016.root","open")
tIn1 = fIn.Get("Events")
#filesloc = "histograms/"
#filesloc = "histograms_array_nw/"
filesloc = "histograms_2016_fits/"
#fIn2 = ROOT.TFile.Open(filesloc+"mll_m40.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_m40.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fine.root","open")
fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit0.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit1.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit2.root","open")
#fIn2 = ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit3.root","open")

filesIn = {}
#filesIn["a15"] = [ROOT.TFile.Open(filesloc+"mll_m15.root","open"),13.0,17.0,15.0,ROOT.kRed]
#filesIn["a20"] = [ROOT.TFile.Open(filesloc+"mll_m20.root","open"),18.0,22.0,20.0,ROOT.kOrange]
#filesIn["a25"] = [ROOT.TFile.Open(filesloc+"mll_m25.root","open"),23.0,27.0,25.0,ROOT.kYellow]
#filesIn["a30"] = [ROOT.TFile.Open(filesloc+"mll_m30.root","open"),28.0,32.0,30.0,ROOT.kGreen]
#filesIn["a35"] = [ROOT.TFile.Open(filesloc+"mll_m35.root","open"),33.0,37.0,35.0,ROOT.kBlue]
##filesIn["a40"] = [ROOT.TFile.Open(filesloc+"mll_m40.root","open"),38.0,42.0,40.0,ROOT.kMagenta]
#filesIn["a45"] = [ROOT.TFile.Open(filesloc+"mll_m45.root","open"),43.0,47.0,45.0,ROOT.kViolet]
#filesIn["a50"] = [ROOT.TFile.Open(filesloc+"mll_m50.root","open"),48.0,52.0,50.0,ROOT.kSpring]
#filesIn["a55"] = [ROOT.TFile.Open(filesloc+"mll_m55.root","open"),53.0,57.0,55.0,ROOT.kCyan]
#filesIn["a60"] = [ROOT.TFile.Open(filesloc+"mll_m60.root","open"),58.0,62.0,60.0,ROOT.kAzure]


#hInSig =  fIn2.Get("mmmt_inclusive/a40")
fIn2.cd()
#fIn2.Get("mmmt_inclusive").cd()
hSignals = {}

#for file in filesIn.keys():
#    filesIn[file][0].cd()
#    filesIn[file][0].Get("mmmt_inclusive").cd() # for the mll_m15.root like files
#    hSignals[file] = ROOT.gDirectory.Get(file).Clone()

filesIn["a40"] = [ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit0.root","open"),38.0,42.0,40.0,ROOT.kMagenta]
#filesIn["a40"] = [ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit1.root","open"),38.0,42.0,40.0,ROOT.kMagenta]
#filesIn["a40"] = [ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit2.root","open"),38.0,42.0,40.0,ROOT.kMagenta]
#filesIn["a40"] = [ROOT.TFile.Open(filesloc+"final_mmmt_inclusive_mll_fit3.root","open"),38.0,42.0,40.0,ROOT.kMagenta]
filesIn["a40"][0].cd()
filesIn["a40"][0].Get("mmmt_inclusive").cd() # for the mll_m15.root like files
hSignals["a40"] = ROOT.gDirectory.Get("a40").Clone()


#hInData = fIn2.Get("mmmt_inclusive/data_obs")
#hInBkg = fIn2.Get("mmmt_inclusive/DY")
#hInBkg.Add(fIn2.Get("mmmt_inclusive/W"))
#hInBkg.Add(fIn2.Get("mmmt_inclusive/TT"))
#hInBkg.Add(fIn2.Get("mmmt_inclusive/ST"))
#hInBkg.Add(fIn2.Get("mmmt_inclusive/WZ"))
#hInBkg.Add(fIn2.Get("mmmt_inclusive/ZZ"))

hInData = fIn2.Get("mmmt_inclusive/data_obs")
#hInSig = fIn2.Get("a40")   # signal distribution included above!
hInBkg = fIn2.Get("mmmt_inclusive/FF_1")
#hInBkg.Add(fIn2.Get("mmmt_inclusive/vbf"))
hInBkg.Add(fIn2.Get("mmmt_inclusive/WHTT"))
hInBkg.Add(fIn2.Get("mmmt_inclusive/HZJ"))
hInBkg.Add(fIn2.Get("mmmt_inclusive/Other"))

fitParams = {}
dataHists = {}
fitModels = {}
pdfs = {}

for file in filesIn.keys():
    fitParams[file] = [
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", filesIn[file][1], filesIn[file][2]),#works for fine binning
        ROOT.RooRealVar("MH",    "m_{#mu #mu}", 15.0, 60.0),
        ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    filesIn[file][3], filesIn[file][1], filesIn[file][2], "GeV"),
        ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#mu #mu}",  0.0,  2.0, "GeV"),
        ROOT.RooRealVar("lAlpha_"+str(file),   "#alpha of lorentz profile",     0.3, 0.0,      5.0)
        ]

for file in filesIn.keys():
    dataHists[file] = [ROOT.RooDataHist("dh_"+str(file),"signal histo ",ROOT.RooArgList(fitParams[file][0]), hSignals[file])]
    fitModels[file] = [
        ROOT.RooVoigtian("voigtian_"+str(file),   "first voigtian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][3], fitParams[file][2]),
        ROOT.RooRealVar("signalEvents_"+str(file), "",  hSignals[file].GetSumOfWeights(),   0.0, 1000.0),
        ]

for file in filesIn.keys():
    pdfs[file] = ROOT.RooAddPdf("model_"+str(file), "", ROOT.RooArgList(fitModels[file][0]),ROOT.RooArgList(fitModels[file][1]))

#adding background and data RooDatasets
data = ROOT.RooDataHist("data_obs","data histo ",ROOT.RooArgList(fitParams["a40"][0]), hInData)

#setting up background model
norm_bkg = ROOT.RooRealVar("bkg_norm",   "background Normalization",hInBkg.Integral(),0.0,hInBkg.Integral()*2)
norm_sig = ROOT.RooRealVar("sig_norm",   "signal Normalization",hSignals["a40"].Integral(),0.0,hSignals["a40"].Integral()*2)

#norm_bkg = ROOT.RooAbsReal("bkg_norm",   "background Normalization",0.0,10000.0)
x_bkg = ROOT.RooRealVar("MH","m_{#mu#mu}", filesIn["a40"][1], filesIn["a40"][2]) # is this min and max?
c0_bkg = ROOT.RooRealVar("c0_bkg",   "coeff. of bernstein 0",0.1,10.0)
c1_bkg = ROOT.RooRealVar("c1_bkg",  "coeff. of bernstein 1",0.1,10.0)
c2_bkg = ROOT.RooRealVar("c2_bkg",   "coeff. of bernstein 2",0.1,10.0)

bkgHists = ROOT.RooDataHist("bkg","background histo ",ROOT.RooArgList(fitParams["a40"][0]), hInBkg)


bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",x_bkg,ROOT.RooArgList(c0_bkg,c1_bkg,c2_bkg))
sigfit = ROOT.RooVoigtian("sigfit",   "sigfit", fitParams["a40"][0], fitParams["a40"][1], fitParams["a40"][3], fitParams["a40"][2])

bkgfitModel = ROOT.RooExtendPdf("bkg","bkg",bkgfit, norm_bkg)
sigfitModel = ROOT.RooExtendPdf("sig","sig",sigfit, norm_sig)
#making overall signal and bkg in single category

overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 10.0, 65.0)
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()


for file in filesIn.keys():
    pdfs[file].fitTo(dataHists[file][0],ROOT.RooFit.Extended())
    dataHists[file][0].plotOn(massFrame)
    pdfs[file].plotOn(massFrame, ROOT.RooFit.LineColor(filesIn[file][4]), ROOT.RooFit.LineStyle(ROOT.kDashed))

c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)
massFrame.Draw()

c.SaveAs("DiMuonMass_combined_sig.pdf")
c.SaveAs("DiMuonMass_combined_sig.png")
c.Clear()

print fitParams
print fitModels
print pdfs


c.cd()


massFrame = overmass.frame()

bkgfitModel.fitTo(bkgHists,ROOT.RooFit.Extended())
bkgfitModel.paramOn(massFrame)
sigfitModel.fitTo(dataHists["a40"][0],ROOT.RooFit.Extended())
sigfitModel.paramOn(massFrame)
dataHists["a40"][0].plotOn(massFrame)
sigfitModel.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))
bkgfitModel.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed))

massFrame.Draw()

c.SaveAs("DiMuonMass_sig_bkg.pdf")
c.SaveAs("DiMuonMass_sig_bkg.png")
c.Clear()

#import data and PDF into workspaco
workspace = ROOT.RooWorkspace("w")

#workspace.import(cat)
#for file in filesIn.keys():
#    getattr(workspace,'import')(dataHists[file][0])
#    getattr(workspace,'import')(pdfs[file])
#dataHists["a40"][0].SetName("data_obs")
getattr(workspace,'import')(data)
#bkgHists.SetName("data_obs")
dataHists["a40"][0].SetName("sig")
getattr(workspace,'import')(dataHists["a40"][0])
bkgHists.SetName("bkg")
getattr(workspace,'import')(bkgHists)

#getattr(workspace,'import')(pdfs["a40"])
getattr(workspace,'import')(sigfitModel)
getattr(workspace,'import')(bkgfitModel)
getattr(workspace,'import')(norm_bkg)

workspace.Print()

workspace.writeToFile("HToAAWorkspace_combined.root")


del workspace
#sigCat[0].plotOn(massFrame[0])

##Make the datacard
outFile = open("datacard.txt","w")#write mode

outFile.write("imax 1\n")
#outFile.write("jmax "+str(hInData.GetNbinsX()-1)+"\n")
outFile.write("jmax 1\n")
outFile.write("kmax *\n")
outFile.write("---------------\n")
outFile.write("shapes * * HToAAWorkspace_combined.root w:$PROCESS\n")
outFile.write("---------------\n")

outFile.write("bin         bin1   \n")
#outFile.write("observation   "+str(hSignals["a40"].Integral())+"\n")
#outFile.write("observation   "+str(bkgHists.sumEntries())+"\n")
outFile.write("observation   -1 \n")

outFile.write("------------------------------\n")
outFile.write("bin                      bin1       bin1 \n")
outFile.write("process                 sig     bkg\n")
outFile.write("process                 0     1 \n")
#outFile.write("rate                   "+str(hSignals["a40"].Integral())+"   "+str(hInBkg.Integral())+"   \n")
outFile.write("rate                   "+str(dataHists["a40"][0].sumEntries())+"   "+str(bkgHists.sumEntries())+"   \n")
outFile.write("------------------------------\n")
outFile.write("lumi     lnN              1.1    1.0   \n")
outFile.write("sigmaM_a40  param          1.0    0.1   \n")

outFile.close()
