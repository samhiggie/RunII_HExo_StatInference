#########################
#Author: Sam Higginbotham
'''

* File Name : fittingClasses.py

* Purpose : To interface with RooFit and use shapes and rooWorkspaces

* Creation Date : 16-12-2021

* Last Modified :

'''
#########################
import ROOT
import ROOT.RooFit
import numpy as np
from datetime import datetime
import argparse
parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
#parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
parser.add_argument("-id",  "--inputDir", default="",  help="postfix string from previous MakeDataCard step")
parser.add_argument("-o",  "--output", default="",  help="postfix string")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-cards",  "--cards", default=False,action='store_true',  help="create the datacards")
parser.add_argument("-workspace",  "--workspace", default=True,action='store_false',  help="create the workspace")
parser.add_argument("-opt",  "--optimize", default=False,action='store_true',  help="do double fit")
parser.add_argument("-seterror",  "--seterror", default=False,action='store_true',  help="save error on background to half of nominal fit value")
parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
args = parser.parse_args()

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
class office():
    def __init__(self):
        self.wsp = ROOT.RooWorkspace("w")
        self.info = {} #going to contain the year and channel
        #outFile = open("datacard_full_"+args.output+".txt","w")#write mode
        #self.txtfile = open("","w")

    def hireFunction(self,function):
        getattr(self.wsp,"import")(function)
        return

    def printCards(self):

        return

class shape():
    def __init__(self):
        self.name = ""
        self.variables = []
        self.vars={}
        self.newvariables = {}
        self.newvariablesbins = []
        self.cuts = {}
        self.newvarcuts = {}
        self.binning = []
        self.extraplots = {}
        self.systematic = "Nominal"
        self.tree = ROOT.TTree()
        self.tfile = ROOT.TFile()
        self.poi = ROOT.RooRealVar("mll","m_{#mu#mu}", 18, 62)
        self.pdf = {} # dictionary for various types of pdfs
        self.coeffs = {} # dictionary for RooRealVars
        self.data = ROOT.RooDataSet()
        self.fitresults = {} # dictionary to hold fit results from the PDFs
        self.plotframe = {}
        self.canvas = {}

    def connectFile(self,path):
        print("trying to open ",path)
        self.tfile = ROOT.TFile.Open(path,"read")
        self.tfile.ls()
        return

    def fillTree(self,pathtotree):
        self.tree = self.tfile.Get(pathtotree)
        self.tree.GetEntries()
        return

    #order is optional for signal
    def fillPDFs(self,type,order,dist):
        self.name = dist+"_"+self.systematic+"_"+args.output
        self.coeffs[type+"_"+str(order)]=[]
        if type=="bernstein":
            for coeff in range(0,order+1):
                ##self.coeffs[type][str(coeff)] = \
                self.coeffs[type+"_"+str(order)].append(\
                ROOT.RooRealVar("c"+str(coeff)+"_"+str(dist)+"_"+str(self.systematic),
                        "c"+str(coeff)+"_"+str(dist)+"_"+str(self.systematic),
                        0.0,
                        -1000000.0,
                        1000000.0))
            self.pdf[type+"_"+str(order)]= \
            ROOT.RooBernstein(str(dist)+"_"+str(type)+"_"+str(order)+"_fit_"+str(self.systematic),
                            str(dist)+"_"+str(type)+"_"+str(order)+"_fit_"+str(self.systematic),
                            self.poi,
                            ROOT.RooArgList(*self.coeffs[type+"_"+str(order)])
            )

        #if type=="voigtian":

        return

    def fillDataset(self,dist):
        self.data = ROOT.RooDataSet(dist,dist,ROOT.RooArgSet(self.poi),ROOT.RooFit.Import(self.tree))

        return

    def fitToData(self,dist):

        for pdfkey,pdf in self.pdf.items():
            self.fitresults[pdfkey]=pdf.fitTo(
                self.data,ROOT.RooFit.Range(18,62),
                ROOT.RooFit.Minimizer("Minuit2","migrad"),
                ROOT.RooFit.Save()
            )
            print("for type :",pdfkey)
            print(self.fitresults[pdfkey].Print('V'))
        return

    def finalFitToData(self,dist):
        for pdfkey,pdf in self.pdf.items():
            for coefnum in range(len(self.coeffs[pdfkey])): # this is a list
                self.coeffs[pdfkey][coefnum].setRange(
                            self.coeffs[pdfkey][coefnum].getValV()/1.50,
                            self.coeffs[pdfkey][coefnum].getValV()*1.50)

            self.fitresults[pdfkey]=pdf.fitTo(
                self.data,ROOT.RooFit.Range(18,62),
                #ROOT.RooFit.Minimizer("Minuit2","migrad"),
                ROOT.RooFit.Minimizer("Minuit2","minos"),
                ROOT.RooFit.Save()
            )
            print("for type :",pdfkey)
            print(self.fitresults[pdfkey].Print('V'))
        return

    def printParamsNLLs(self,dist):
        nllresults = open("nllresults_"+dist+"_"+self.name+".txt","w")
        nllresults.write("fit results for "+self.name+"\n")
        for coeffkey,coeff in self.coeffs.items():
            nllresults.write("final fit values for the coff of "+str(coeffkey)+"\n")
            for value in coeff:
                nllresults.write(value.GetName()+" "+str(value.getValV())+" with error "+str(value.getError())+"\n")


        for fitresultkey,fitresult in self.fitresults.items():
            print("nll for type ",fitresultkey,"   ",fitresult.minNll())
            nllresults.write("nll for type "+str(fitresultkey)+"   "+str(fitresult.minNll())+"\n")
        fitorder = list(self.fitresults.keys())
        for ordnum in range(len(fitorder)-1):
            print("\u0394nll  between  ",ordnum,"   and " ,ordnum+1,"   ",self.fitresults[fitorder[ordnum]].minNll()-self.fitresults[fitorder[ordnum+1]].minNll())
            nllresults.write("\u0394nll  between  "+str(ordnum)+"   and " +str(ordnum+1)+"   "+str(self.fitresults[fitorder[ordnum]].minNll()-self.fitresults[fitorder[ordnum+1]].minNll())+"\n")
        return

    def createPlots(self,type,order,dist):
        self.plotframe[type+"_"+str(order)] = self.poi.frame()
        self.canvas[type+"_"+str(order)] = ROOT.TCanvas("c","",600,600)
        self.canvas[type+"_"+str(order)].cd()
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gPad.SetLeftMargin(0.15)
        self.plotframe[type+"_"+str(order)].GetYaxis().SetTitleOffset(1.6)
        self.plotframe[type+"_"+str(order)].GetXaxis().SetTitle("M_{#mu#mu}")
        ROOT.TGaxis().SetMaxDigits(2)
        self.data.plotOn(self.plotframe[type+"_"+str(order)],ROOT.RooFit.Binning(25))
        self.pdf[type+"_"+str(order)].paramOn(self.plotframe[type+"_"+str(order)])
        self.pdf[type+"_"+str(order)].plotOn(self.plotframe[type+"_"+str(order)], ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Range(20,60))
        self.plotframe[type+"_"+str(order)].Draw()
        self.canvas[type+"_"+str(order)].Draw()
        self.canvas[type+"_"+str(order)].SaveAs(type+"_"+str(order)+"_"+str(self.name)+".png")



        return
