#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/TMVAMultiClassGui.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"

#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

TString outfileName("TMVA.root");

void Cuts(TString prefix = "")
{
    gErrorIgnoreLevel = kFatal;
    TString myMethodList = "";
    TMVA::Tools::Instance();
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    TFile* file = new TFile("/eos/user/b/bharikri/CERNBox_synced/Projects/HeavyIon/Results/2022_03_17_Photon_Final_Reproduction/1_Skimming_weighting/QCDPhoton.root");
    TFile* file1 = new TFile("/eos/user/b/bharikri/CERNBox_synced/Projects/HeavyIon/Results/2022_03_17_Photon_Final_Reproduction/1_Skimming_weighting/EMEnrichedDijet.root");

    TTree* SigTree = (TTree*)file->Get("pho_tree");
    TTree* BkgTree = (TTree*)file1->Get("pho_tree");

    float nEv = SigTree->GetEntries();

    TFile *outputFile = TFile::Open(prefix+outfileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("Cuts", outputFile, "V:!Silent:Color:Correlations:DrawProgressBar:AnalysisType=Classification:Transformations=D"); 
    // AnalysisType=Classification for binary 
    // Transformations transforms the input variables. I is identity. Other options are D;P;G

    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    dataloader->AddSignalTree(SigTree);
    dataloader->AddBackgroundTree(BkgTree);
    // dataloader->SetWeightExpression(prefix+"weight_pthat");

    dataloader->SetSignalWeightExpression("weight*weight_pthat*weight_cent");
    dataloader->SetBackgroundWeightExpression("weight*weight_pthat*weight_cent");

    dataloader->AddSpectator("phoEt");
    dataloader->AddSpectator("phoEta");
    dataloader->AddSpectator("phoPhi");
    dataloader->AddSpectator("phoSCEta");
    dataloader->AddSpectator("phoSCPhi");
    dataloader->AddSpectator("phoHoverE");
    dataloader->AddSpectator("pho_ecalClusterIsoR4");
    dataloader->AddSpectator("pho_hcalRechitIsoR4");
    dataloader->AddSpectator("pho_trackIsoR4PtCut20");
    dataloader->AddSpectator("pho_genMatchedIndex");
    dataloader->AddSpectator("mcPID");
    dataloader->AddSpectator("hiBin");
    dataloader->AddSpectator("rho");
    
    dataloader->AddVariable("phoHoverE",'F');
    dataloader->AddVariable("phoSigmaIEtaIEta_2012",'F');
    // dataloader->AddVariable("phoSumISO := pfcIso3subUE+pfnIso3subUE+pfpIso3subUE",'F');
    dataloader->AddVariable("phoSumISO := pho_ecalClusterIsoR3+pho_hcalRechitIsoR3+pho_trackIsoR3PtCut20",'F');
    
    std::cout<<"Add Variable Step"<<endl;

    // Cuts are applied before the training is performed. Only cuts applied are to remove NaN
    TCut Quality = "hiBin <= 180 && phoSigmaIEtaIEta_2012>0.002 && !(phoEta<-1.39 && phoPhi<-0.9 && phoPhi>-1.6)"; // Centrality for 2018 PbPb + noise/spike + HEM (already applied in skim sample)
    TCut GenPhoton = "pho_genMatchedIndex > -1 && mcPID == 22"; // True Photons
    TCut GenSig    = "mcCalIsoDR04 < 5 && (abs(mcMomPID) <= 22 || mcMomPID == -999)";
    TCut PhotonKin = "phoEtCorrected > 40 && pho_trackIsoR3PtCut20 > -500 && pho_trackIsoR3PtCut20 < 500";
    TCut SelEB = "fabs(phoSCEta) < 1.48";
    TCut SelEE = "fabs(phoSCEta) > 1.57 && fabs(phoSCEta) < 2.0";
    TCut mycut, mycutbkg; 
    mycut = Quality && GenPhoton && GenSig && SelEB && PhotonKin;
    mycutbkg = Quality && (!GenPhoton || !GenSig) && SelEB && PhotonKin;

    std::stringstream text; 
    Int_t n_sig= int(0.3333333*SigTree->GetEntries(mycut));
    Int_t n_bkg= int(0.3333333*BkgTree->GetEntries(mycutbkg));
    text << "SplitMode=Random:SplitSeed=12345:NormMode=EqualNumEvents:nTrain_Signal=" << n_sig<<":nTrain_Background="<<n_bkg;
        
    dataloader->PrepareTrainingAndTestTree(mycut,mycutbkg,(text.str()).c_str());
    
    factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts", "V:!H:FitMethod=GA:EffMethod=EffSel:VarProp=FMin:PopSize=600:Steps=60"); 

    factory->TrainAllMethods();
    // To move the xml file
    // TMVA::MethodBase* mBase =  dynamic_cast<TMVA::MethodBase*>(factory->GetMethod("dataset", "Cuts"));
    // gSystem->Exec(Form("cp -v %s %s", mBase->GetWeightFileName().Data(), "/home/llr/cms/bharikri/Projects/Template_fit/Calorimeter_Isolation/"));  

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

    outputFile->Close();

    delete factory;
    delete dataloader;

    // Uncomment the following lines for GUI
    // if (!gROOT->IsBatch())
    // TMVA::TMVAGui(outfileName);
}

int main(int argc, char* argv[]) {

    // if(argc!=2)  return -1;

    Cuts("EB_");
    return 0;
}
