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

#include <sstream>
#include <fstream>

using namespace std;

TString outfileName("TMVA.root");

void RegressionTraining(TString prefix = "")
{
    gErrorIgnoreLevel = kFatal;
    TString myMethodList = "";
    TMVA::Tools::Instance();
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    TFile* file = new TFile("/eos/user/b/bharikri/CERNBox_synced/Projects/HeavyIon/Results/2022_03_17_Photon_Final_Reproduction/1_Skimming_weighting/QCDPhoton_Bharad.root");

    TTree* SignalTree = (TTree*)file->Get("pho_tree");

    TFile *outputFile = TFile::Open(prefix+outfileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("tmvaTrainRegr_BDTG_pbpb18", outputFile, "V:!Silent:Color:DrawProgressBar:AnalysisType=Regression:Transformations=I"); 
    // AnalysisType=Classification for binary 
    // Transformations transforms the input variables. I is identity. Other options are D;P;G

    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    dataloader->AddRegressionTree(SignalTree);
    std::cout<<"Add Regression Tree Step"<<endl;
    // dataloader->SetWeightExpression("weight_pthat"); // Regression doesn't need weights

    dataloader->AddTarget("mcE");
    std::cout<<"Add Target Step"<<endl;

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
    dataloader->AddSpectator("hiBin");
    dataloader->AddSpectator("rho");
    dataloader->AddSpectator("phoE");

    dataloader->AddVariable("phoSCRawE",'F');
    dataloader->AddVariable("phoSCEta",'F');
    dataloader->AddVariable("phoSCPhi",'F');
    dataloader->AddVariable("phoSCEtaWidth",'F');
    dataloader->AddVariable("phoSCPhiWidth",'F');
    dataloader->AddVariable("phoE3x3_2012",'F');
    dataloader->AddVariable("phoMaxEnergyXtal_2012",'F');
    dataloader->AddVariable("phoE2nd_2012",'F');
    dataloader->AddVariable("phoE_LR := (phoELeft_2012-phoERight_2012)/(phoELeft_2012+phoERight_2012)",'F');
    dataloader->AddVariable("phoE_TB := (phoETop_2012-phoEBottom_2012)/(phoETop_2012+phoEBottom_2012)",'F');
    dataloader->AddVariable("phoSigmaIEtaIEta_2012",'F');
    dataloader->AddVariable("phoSigmaIEtaIPhi_2012",'F');
    dataloader->AddVariable("phoSigmaIPhiIPhi_2012",'F');
    dataloader->AddVariable("rho",'F');
    if(prefix.Contains("EE"))
        dataloader->AddVariable("phoESEn",'F');
    
    std::cout<<"Add Variable Step"<<endl;

    // Cuts are applied before the training is performed.
    TCut Quality = "hiBin <= 180 && phoSigmaIEtaIEta_2012>0.002 && !(phoEta<-1.39 && phoPhi<-0.9 && phoPhi>-1.6)"; // Centrality for 2018 PbPb + noise/spike + HEM (already applied in skim sample)
    TCut GenPhoton = "pho_genMatchedIndex > -1 && mcPID == 22"; // True Photons
    TCut GenSig    = "mcCalIsoDR04 < 5 && (abs(mcMomPID) <= 22 || mcMomPID == -999)";
    TCut SelEB = "fabs(phoSCEta) < 1.48";
    TCut SelEE = "fabs(phoSCEta) > 1.57 && fabs(phoSCEta) < 2.0";
    TCut mycut; 
    if(prefix.Contains("EE")){
        mycut = Quality && GenPhoton && GenSig && SelEE;
    }
    else{
        mycut = Quality && GenPhoton && GenSig && SelEB;
    }
    
    std::stringstream text; 
    Int_t n_reg= int(0.3333333*SignalTree->GetEntries(mycut));
    text << "SplitMode=Random:SplitSeed=12345:NormMode=EqualNumEvents:nTrain_Regression=" << n_reg;
    
    dataloader->PrepareTrainingAndTestTree(mycut, (text.str()).c_str());
    
    std::cout<<"PrepareTrainingAndTestTree Step"<<endl;

    factory->BookMethod(dataloader, TMVA::Types::kBDT, prefix.IsNull()?"EB":TString(prefix(0,2)), "V:!H:DoBoostMonitor=True:NTrees=2000:MinNodeSize=0.2%:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50");

    factory->TrainAllMethods();
    std::cout << "==> TMVA Training is done!" << std::endl;

    factory->TestAllMethods();
    std::cout << "==> TMVA Testing is done!" << std::endl;

    factory->EvaluateAllMethods();
    std::cout << "==> TMVA Evaluation is done!" << std::endl;

    outputFile->Close();

    std::cout << "==> Wrote root file: " <<prefix.Data() << outputFile->GetName()<< std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    std::cout <<std::endl<< "Weights.xml File is generated in the dataset folder. Copy it into input folder for implementation" << std::endl;

    delete factory;
    delete dataloader;

    // Uncomment the following lines for GUI
    // if (!gROOT->IsBatch())
    // TMVA::TMVARegGui(outfileName);
}

int main(){

    // gSystem->Exec("rm -rf dataset/ EE_TMVA.root EB_TMVA.root"); //Delete Previous training results

    RegressionTraining("EB_");
    RegressionTraining("EE_");

    return 0;
}