#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TMath.h>
#include <TSystem.h>

#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "phoERegression.h"
#endif

using namespace std;

void DrawLatex(Float_t , Float_t , Int_t , const char* , Float_t );
void SetPlotStyle();
void displayProgress(long, long);

Double_t RelativePhi(Double_t mphi,Double_t vphi) {

  if (vphi < -1 * TMath::Pi())
    vphi += (2 * TMath::Pi());
  else if (vphi > TMath::Pi())
    vphi -= (2 * TMath::Pi());
  if (mphi < -1 * TMath::Pi())
    mphi += (2 * TMath::Pi());
  else if (mphi > TMath::Pi())
    mphi -= (2 * TMath::Pi());
  double dphi = mphi - vphi;
  if (dphi < -1 * TMath::Pi())
    dphi += (2 * TMath::Pi());
  else if (dphi > TMath::Pi())
    dphi -= (2 * TMath::Pi());
  return dphi; // dphi in [-Pi, Pi]
}

void SkimFile_MC_allPho_Bharad(){

    TString type="QCDPhoton";
    TString type_dir = "/home/llr/cms/bharikri/Projects/Photon_Analysis/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/PhotonAnalysis/test/1_Skimming_weighting/";            

   gROOT->SetBatch();
   gErrorIgnoreLevel = kWarning;

   // Centrality weights 
   const float Ncoll[200] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
   const float SIGCUTOFF[6] = {0,29.808813095,49.885753632,79.932403564,119.952224731,169.969177246};
   const float SIGCUTOFF_MAX[6] = {29.808813095,49.885753632,79.932403564,119.952224731,169.969177246,0};
   Double_t SIGpthat_weight[6] = {1,0.206865924,0.047953729,0.009373522,0.002099855,0.000533636}; 
   const float BKGCUTOFF[6] = {0,49.908184052,79.949401855,119.969474792,169.972900391,219.981445312};
   const float BKGCUTOFF_MAX[6] = {49.908184052,79.949401855,119.969474792,169.972900391,219.981445312,0};
   Double_t BKGpthat_weight[6] = {1,0.470964235,0.205821138,0.081499942,0.023460744,0.007632521}; 
   
   SetPlotStyle();
   TChain data("data"), EventTree("EventTree"), HiTree("HiTree"),skimanalysis("skimanalysis"),hltanalysis("hltanalysis"),jets("jets");

   if(type.Contains("EMEnrichedDijet")){
   // EventTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_30/HiForestAOD_*.root/runAnalyzer/run");
   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_30/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_30/HiForestAOD_*.root/hltanalysis/HltTree");

   // EventTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_50/HiForestAOD_*.root/runAnalyzer/run");
   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_50/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_50/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_50/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_50/HiForestAOD_*.root/hltanalysis/HltTree");

   // EventTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_80/HiForestAOD_*.root/runAnalyzer/run");
   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_80/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_80/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_80/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_80/HiForestAOD_*.root/hltanalysis/HltTree");

   // EventTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_120/HiForestAOD_*.root/runAnalyzer/run");
   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_120/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_120/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_120/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_120/HiForestAOD_*.root/hltanalysis/HltTree");

   // EventTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_170/HiForestAOD_*.root/runAnalyzer/run");
   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_170/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_170/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_170/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_170/HiForestAOD_*.root/hltanalysis/HltTree");

   data.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_220/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_220/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_220/HiForestAOD_*.root/skimanalysis/HltTree");
   hltanalysis.Add("/data_CMS/cms/bharikri/EMEnrichedDijet-pthat_220/HiForestAOD_*.root/hltanalysis/HltTree");
   }

   if(type.Contains("QCDPhoton")){
      
   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_15.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15/HiForestAOD_*.root/hltanalysis/HltTree");

   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_15_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_15_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");

   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_30.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30/HiForestAOD_*.root/hltanalysis/HltTree");

   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_30_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_30_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");
   
   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_50.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50/HiForestAOD_*.root/hltanalysis/HltTree");

   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_50_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_50_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");

   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_80.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80/HiForestAOD_*.root/hltanalysis/HltTree");

   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_80_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_80_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");

   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_120.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120/HiForestAOD_*.root/hltanalysis/HltTree");

   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_120_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_120_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");

   data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170/HiForestAOD_*.root/skimanalysis/HltTree");
   jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_170.root/akCs2PFJetAnalyzer_substructure/t");
   hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170/HiForestAOD_*.root/hltanalysis/HltTree");
   
   // data.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170_Filter_30/HiForestAOD_*.root/ggHiNtuplizerGED/EventTree");
   // HiTree.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170_Filter_30/HiForestAOD_*.root/hiEvtAnalyzer/HiTree");
   // skimanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170_Filter_30/HiForestAOD_*.root/skimanalysis/HltTree");
   // jets.Add("/data_CMS/cms/bharikri/clusters/QCDPhoton_pthat_170_Filter_30.root/akCs2PFJetAnalyzer_substructure/t");
   // hltanalysis.Add("/data_CMS/cms/bharikri/QCDPhoton_pthat_170_Filter_30/HiForestAOD_*.root/hltanalysis/HltTree");
   }

   // data.AddFriend("EventTree");
   data.AddFriend("HiTree");
   data.AddFriend("skimanalysis");
   data.AddFriend("hltanalysis");
   // data.AddFriend("jets");

  // ----------------------------------------------------------------------------------------------------------------
  // Variables

   Int_t hiBin=0;
   float pthat=0;
   float weight = 0;
   
   data.SetBranchAddress("hiBin", &hiBin);
   data.SetBranchAddress("pthat", &pthat);
   data.SetBranchAddress("weight", &weight);

   // reco::GenParticle
   std::vector<int>*    mcPID;
   std::vector<int>*    mcMomPID;
   std::vector<float>*  mcPt;
   std::vector<float>*  mcEta;
   std::vector<float>*  mcPhi;
   std::vector<float>*  mcE;
   std::vector<float>*  mcEt;
   std::vector<float>*  mcMass;
   std::vector<float>*  mcCalIsoDR04;

   // reco::Photon
   std::vector<float>*  phoE;
   std::vector<float>*  phoEt;
   std::vector<float>*  phoEta;
   std::vector<float>*  phoPhi;
   std::vector<float>*  phoHoverE;
   std::vector<float>*  pho_ecalClusterIsoR3;
   std::vector<float>*  pho_hcalRechitIsoR3;
   std::vector<float>*  pho_trackIsoR3PtCut20;
   std::vector<float>*  pho_ecalClusterIsoR4;
   std::vector<float>*  pho_hcalRechitIsoR4;
   std::vector<float>*  pho_trackIsoR4PtCut20;
   std::vector<int>*    pho_genMatchedIndex;
   std::vector<float>*  phoSigmaEtaEta_2012;
   std::vector<float>*  pfcIso3subUE;
   std::vector<float>*  pfnIso3subUE;
   std::vector<float>*  pfpIso3subUE;
   
   std::vector<float>*  pho_swissCrx;
   std::vector<float>*  pho_seedTime;

   std::vector<float>* phoSCRawE;
   std::vector<float>* phoSCEta;
   std::vector<float>* phoSCPhi;
   std::vector<float>* phoSCEtaWidth;
   std::vector<float>* phoSCPhiWidth;
   std::vector<float>* phoE3x3_2012;
   std::vector<float>* phoMaxEnergyXtal_2012;
   std::vector<float>* phoE2nd_2012;
   std::vector<float>* phoELeft_2012;
   std::vector<float>* phoERight_2012;
   std::vector<float>* phoETop_2012;
   std::vector<float>* phoEBottom_2012;
   std::vector<float>* phoSigmaIEtaIEta_2012;
   std::vector<float>* phoSigmaIEtaIPhi_2012;
   std::vector<float>* phoSigmaIPhiIPhi_2012;
   Float_t rho;
   std::vector<float>* phoESEn;

   /// ----

   mcPID = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcMass = 0;
   mcCalIsoDR04 = 0;
   mcMomPID = 0;

   // reco::Photon
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoHoverE = 0;
   pho_ecalClusterIsoR3 = 0;
   pho_hcalRechitIsoR3 = 0;
   pho_trackIsoR3PtCut20 = 0;
   pho_ecalClusterIsoR4 = 0;
   pho_hcalRechitIsoR4 = 0;
   pho_trackIsoR4PtCut20 = 0;
   pho_genMatchedIndex = 0;
   phoSigmaEtaEta_2012 = 0;
   pfcIso3subUE = 0;
   pfnIso3subUE = 0;
   pfpIso3subUE = 0;
   
   pho_swissCrx = 0;
   pho_seedTime = 0;

   phoSCRawE = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoE3x3_2012 = 0;
   phoMaxEnergyXtal_2012 = 0;
   phoE2nd_2012 = 0;
   phoELeft_2012 = 0;
   phoERight_2012 = 0;
   phoETop_2012 = 0;
   phoEBottom_2012 = 0;
   phoSigmaIEtaIEta_2012 = 0;
   phoSigmaIEtaIPhi_2012 = 0;
   phoSigmaIPhiIPhi_2012 = 0;
   rho = 0;
   phoESEn = 0;

   data.SetBranchAddress("mcPID",        &mcPID);
   data.SetBranchAddress("mcPt",         &mcPt);
   data.SetBranchAddress("mcEta",        &mcEta);
   data.SetBranchAddress("mcPhi",        &mcPhi);
   data.SetBranchAddress("mcE",          &mcE);
   data.SetBranchAddress("mcEt",         &mcEt);
   data.SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
   data.SetBranchAddress("mcMass",       &mcMass);
   data.SetBranchAddress("mcMomPID",     &mcMomPID);

   data.SetBranchAddress("phoE",                  &phoE);
   data.SetBranchAddress("phoEt",                 &phoEt);
   data.SetBranchAddress("phoEta",                &phoEta);
   data.SetBranchAddress("phoPhi",                &phoPhi);
   data.SetBranchAddress("phoHoverE",             &phoHoverE);
   data.SetBranchAddress("pho_ecalClusterIsoR3",  &pho_ecalClusterIsoR3);
   data.SetBranchAddress("pho_hcalRechitIsoR3",   &pho_hcalRechitIsoR3);
   data.SetBranchAddress("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20);
   data.SetBranchAddress("pho_ecalClusterIsoR4",  &pho_ecalClusterIsoR4);
   data.SetBranchAddress("pho_hcalRechitIsoR4",   &pho_hcalRechitIsoR4);
   data.SetBranchAddress("pho_trackIsoR4PtCut20", &pho_trackIsoR4PtCut20);
   data.SetBranchAddress("pho_genMatchedIndex",   &pho_genMatchedIndex);
   data.SetBranchAddress("phoSigmaEtaEta_2012",   &phoSigmaEtaEta_2012);
   data.SetBranchAddress("pfcIso3subUE",          &pfcIso3subUE);
   data.SetBranchAddress("pfnIso3subUE",          &pfnIso3subUE);
   data.SetBranchAddress("pfpIso3subUE",          &pfpIso3subUE);
   
   data.SetBranchAddress("phoSCRawE",          &phoSCRawE);
   data.SetBranchAddress("phoSCEta",          &phoSCEta);
   data.SetBranchAddress("phoSCPhi",          &phoSCPhi);
   data.SetBranchAddress("phoSCEtaWidth",          &phoSCEtaWidth);
   data.SetBranchAddress("phoSCPhiWidth",          &phoSCPhiWidth);
   data.SetBranchAddress("phoE3x3_2012",          &phoE3x3_2012);
   data.SetBranchAddress("phoMaxEnergyXtal_2012",          &phoMaxEnergyXtal_2012);
   data.SetBranchAddress("phoE2nd_2012",          &phoE2nd_2012);
   data.SetBranchAddress("phoELeft_2012",          &phoELeft_2012);
   data.SetBranchAddress("phoERight_2012",          &phoERight_2012);
   data.SetBranchAddress("phoETop_2012",          &phoETop_2012);
   data.SetBranchAddress("phoEBottom_2012",          &phoEBottom_2012);
   data.SetBranchAddress("phoSigmaIEtaIEta_2012",          &phoSigmaIEtaIEta_2012);
   data.SetBranchAddress("phoSigmaIEtaIPhi_2012",          &phoSigmaIEtaIPhi_2012);
   data.SetBranchAddress("phoSigmaIPhiIPhi_2012",          &phoSigmaIPhiIPhi_2012);
   data.SetBranchAddress("rho",          &rho);
   data.SetBranchAddress("phoESEn",          &phoESEn);

   data.SetBranchAddress("pho_swissCrx",&pho_swissCrx);
   data.SetBranchAddress("pho_seedTime",&pho_seedTime);

   int pprimaryVertexFilter = 0;
   int pclusterCompatibilityFilter = 0;
   int phfCoincFilter2Th4 = 0;

   data.SetBranchAddress("pprimaryVertexFilter",  &pprimaryVertexFilter);
   data.SetBranchAddress("pclusterCompatibilityFilter",  &pclusterCompatibilityFilter);
   data.SetBranchAddress("phfCoincFilter2Th4",  &phfCoincFilter2Th4);

   Int_t L1_SingleEG21_BptxAND = 0;
   Int_t HLT_HIGEDPhoton40_v1 = 0;

   data.SetBranchAddress("L1_SingleEG21_BptxAND",  &L1_SingleEG21_BptxAND);
   data.SetBranchAddress("HLT_HIGEDPhoton40_v1",  &HLT_HIGEDPhoton40_v1);

// ----------------------------------------------------------------------------------------------------------------
// Histograms
    TTree *pho_tree;
   pho_tree = new TTree("pho_tree","pho_tree");

   Float_t train_weight = 0;
   Float_t train_weight_pthat = 0;
   Float_t train_weight_cent = 0;

   Float_t train_pthat = 0;
   Float_t train_phoSCRawE = 0;
   Float_t train_phoSCEta = 0;
   Float_t train_phoSCPhi = 0;
   Float_t train_phoSCEtaWidth = 0;
   Float_t train_phoSCPhiWidth = 0;
   Float_t train_phoE3x3_2012 = 0;
   Float_t train_phoMaxEnergyXtal_2012 = 0;
   Float_t train_phoE2nd_2012 = 0;
   Float_t train_phoELeft_2012 = 0;
   Float_t train_phoERight_2012 = 0;
   Float_t train_phoETop_2012 = 0;
   Float_t train_phoEBottom_2012 = 0;
   Float_t train_phoSigmaIEtaIEta_2012 = 0;
   Float_t train_phoSigmaIEtaIPhi_2012 = 0;
   Float_t train_phoSigmaIPhiIPhi_2012 = 0;
   Float_t train_rho = 0;
   Float_t train_phoESEn = 0;

   Float_t train_phoEt = 0;
   Float_t train_phoEtCorrected = 0;
   Float_t train_phoEta = 0;
   Float_t train_phoPhi = 0;
   Float_t train_phoHoverE = 0;
   Float_t train_pho_ecalClusterIsoR3 = 0;
   Float_t train_pho_hcalRechitIsoR3 = 0;
   Float_t train_pho_trackIsoR3PtCut20 = 0;
   Float_t train_pho_ecalClusterIsoR4 = 0;
   Float_t train_pho_hcalRechitIsoR4 = 0;
   Float_t train_pho_trackIsoR4PtCut20 = 0;
   Float_t train_phoE = 0;
   Float_t train_phoSigmaEtaEta_2012 = 0;
   
   Float_t train_pfcIso3subUE = 0;
   Float_t train_pfnIso3subUE = 0;
   Float_t train_pfpIso3subUE = 0;
   Float_t train_SumPFIso = 0;
   Float_t train_SumCalIso = 0;
   Int_t   train_hiBin = 0;

   Int_t   train_L1_SingleEG21_BptxAND = 0;
   Int_t   train_HLT_HIGEDPhoton40_v1 = 0;

   Float_t train_mcCalIsoDR04 = 0;
   Float_t train_mcE = 0;
   Float_t train_mcEt = 0;
   Float_t train_mcEta = 0;
   Float_t train_mcPhi = 0;

   Int_t   train_pho_genMatchedIndex = 0;
   Int_t   train_mcPID = 0;
   Int_t   train_mcMomPID = 0;

   pho_tree->Branch("weight",&train_weight);
   pho_tree->Branch("weight_pthat",&train_weight_pthat);
   pho_tree->Branch("weight_cent",&train_weight_cent);
   
   pho_tree->Branch("pthat",&train_pthat);
   pho_tree->Branch("phoSCRawE",&train_phoSCRawE);
   pho_tree->Branch("phoSCEta",&train_phoSCEta);
   pho_tree->Branch("phoSCPhi",&train_phoSCPhi);
   pho_tree->Branch("phoSCEtaWidth",&train_phoSCEtaWidth);
   pho_tree->Branch("phoSCPhiWidth",&train_phoSCPhiWidth);
   pho_tree->Branch("phoE3x3_2012",&train_phoE3x3_2012);
   pho_tree->Branch("phoMaxEnergyXtal_2012",&train_phoMaxEnergyXtal_2012);
   pho_tree->Branch("phoE2nd_2012",&train_phoE2nd_2012);
   pho_tree->Branch("phoELeft_2012",&train_phoELeft_2012);
   pho_tree->Branch("phoERight_2012",&train_phoERight_2012);
   pho_tree->Branch("phoETop_2012",&train_phoETop_2012);
   pho_tree->Branch("phoEBottom_2012",&train_phoEBottom_2012);
   pho_tree->Branch("phoSigmaIEtaIEta_2012",&train_phoSigmaIEtaIEta_2012);
   pho_tree->Branch("phoSigmaIEtaIPhi_2012",&train_phoSigmaIEtaIPhi_2012);
   pho_tree->Branch("phoSigmaIPhiIPhi_2012",&train_phoSigmaIPhiIPhi_2012);
   pho_tree->Branch("rho",&train_rho);
   pho_tree->Branch("phoESEn",&train_phoESEn);

   pho_tree->Branch("phoEt",&train_phoEt);
   pho_tree->Branch("phoEtCorrected",&train_phoEtCorrected);
   pho_tree->Branch("phoEta",&train_phoEta);
   pho_tree->Branch("phoPhi",&train_phoPhi);
   pho_tree->Branch("phoHoverE",&train_phoHoverE);
   pho_tree->Branch("pho_ecalClusterIsoR3",&train_pho_ecalClusterIsoR3);
   pho_tree->Branch("pho_hcalRechitIsoR3",&train_pho_hcalRechitIsoR3);
   pho_tree->Branch("pho_trackIsoR3PtCut20",&train_pho_trackIsoR3PtCut20);
   pho_tree->Branch("pho_ecalClusterIsoR4",&train_pho_ecalClusterIsoR4);
   pho_tree->Branch("pho_hcalRechitIsoR4",&train_pho_hcalRechitIsoR4);
   pho_tree->Branch("pho_trackIsoR4PtCut20",&train_pho_trackIsoR4PtCut20);
   pho_tree->Branch("phoE",&train_phoE);
   pho_tree->Branch("phoSigmaEtaEta_2012",&train_phoSigmaEtaEta_2012);
   
   pho_tree->Branch("pfcIso3subUE",&train_pfcIso3subUE);
   pho_tree->Branch("pfnIso3subUE",&train_pfnIso3subUE);
   pho_tree->Branch("pfpIso3subUE",&train_pfpIso3subUE);
   pho_tree->Branch("SumPFIso",&train_SumPFIso);
   pho_tree->Branch("SumCalIso",&train_SumCalIso);
   pho_tree->Branch("hiBin",&train_hiBin);
   pho_tree->Branch("L1_SingleEG21_BptxAND",&train_L1_SingleEG21_BptxAND);
   pho_tree->Branch("HLT_HIGEDPhoton40_v1",&train_HLT_HIGEDPhoton40_v1);

   pho_tree->Branch("mcCalIsoDR04",&train_mcCalIsoDR04);
   pho_tree->Branch("mcE",&train_mcE);
   pho_tree->Branch("mcEt",&train_mcEt);
   pho_tree->Branch("mcEta",&train_mcEta);
   pho_tree->Branch("mcPhi",&train_mcPhi);

   pho_tree->Branch("pho_genMatchedIndex",&train_pho_genMatchedIndex);
   pho_tree->Branch("mcPID",&train_mcPID);
   pho_tree->Branch("mcMomPID",&train_mcMomPID);

   phoERegression myPhoERegr;
   // myPhoERegr.initiliazeReaderEB("/home/llr/cms/bharikri/Projects/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/JetAnalysis/test/dataset/weights/tmvaTrainRegr_BDTG_pbpb18_EB.weights.xml");
   // myPhoERegr.initiliazeReaderEE("/home/llr/cms/bharikri/Projects/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/JetAnalysis/test/dataset/weights/tmvaTrainRegr_BDTG_pbpb18_EE.weights.xml");

   myPhoERegr.initiliazeReaderEB("/home/llr/cms/bharikri/Projects/Photon_Analysis/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/PhotonAnalysis/test/1_Skimming_weighting//tmvaTrainRegr_BDTG_pbpb18_EB.xml");
   myPhoERegr.initiliazeReaderEE("/home/llr/cms/bharikri/Projects/Photon_Analysis/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/PhotonAnalysis/test/1_Skimming_weighting//tmvaTrainRegr_BDTG_pbpb18_EE1.xml");

// ----------------------------------------------------------------------------------------------------------------
//        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
// ----------------------------------------------------------------------------------------------------------------

   Int_t nEv=data.GetEntries();
   std::cout<<nEv<<endl;
   float scale = 1;

   for(int iEntry=0; iEntry< nEv; iEntry++){
      displayProgress(iEntry,nEv);
      data.GetEntry(iEntry);
      scale =1;

      if(type.Contains("QCDPhoton")){
         if(pthat<SIGCUTOFF_MAX[0]){
            scale *=SIGpthat_weight[0];
        }
        else if(pthat<SIGCUTOFF_MAX[1]){
            scale *=SIGpthat_weight[1];
        }
        else if(pthat<SIGCUTOFF_MAX[2]){
            scale *=SIGpthat_weight[2];
        }
        else if(pthat<SIGCUTOFF_MAX[3]){
            scale *=SIGpthat_weight[3];
        }
        else if(pthat<SIGCUTOFF_MAX[4]){
            scale *=SIGpthat_weight[4];
        }
        else {
            scale *=SIGpthat_weight[5];
        }
      }
      else if(type.Contains("EMEnrichedDijet")){
         if(pthat<BKGCUTOFF_MAX[0]){
            scale *=BKGpthat_weight[0];
        }
        else if(pthat<BKGCUTOFF_MAX[1]){
            scale *=BKGpthat_weight[1];
        }
        else if(pthat<BKGCUTOFF_MAX[2]){
            scale *=BKGpthat_weight[2];
        }
        else if(pthat<BKGCUTOFF_MAX[3]){
            scale *=BKGpthat_weight[3];
        }
        else if(pthat<BKGCUTOFF_MAX[4]){
            scale *=BKGpthat_weight[4];
        }
        else {
            scale *=BKGpthat_weight[5];
        }
      }

      int pho_index=-1, genin=-1;;
      float Etmax=-1;
      float corrected_Et = -1;

      std::vector<float> phoVars;
      phoVars.reserve(17);

      if(pprimaryVertexFilter<=0) continue;
      if(pclusterCompatibilityFilter<=0) continue;
      if(phfCoincFilter2Th4<=0) continue;
      // if(hiBin>180) continue;

      for(int ipho=0; ipho<phoEt->size(); ipho++){
        // if(!(fabs(phoEta->at(ipho))<1.44)) continue; // || (fabs(phoEta->at(ipho))>1.6 && fabs(phoEta->at(ipho))<2.0)
        // if(pho_genMatchedIndex->at(ipho)==-1) continue;
        // if(mcPID->at(pho_genMatchedIndex->at(ipho)) != 22) continue;
        // if(mcCalIsoDR04->at(pho_genMatchedIndex->at(ipho))>=5) continue;
        if(phoEta->at(ipho)<-1.39 && phoPhi->at(ipho)<-0.9 && phoPhi->at(ipho)>-1.6) continue;  // HEM Failure
        // if(phoEt->at(ipho)<30) continue;                          // Minimum ET cut
        if(phoSigmaEtaEta_2012->at(ipho)<=0.002) continue;        // Spike Cuts   
        if(fabs(pho_seedTime->at(ipho))>=3 || pho_swissCrx->at(ipho)>=0.9) continue;


        phoVars[0] = phoSCRawE->at(ipho);
        phoVars[1] = phoSCEta->at(ipho);
        phoVars[2] = phoSCPhi->at(ipho);
        phoVars[3] = phoSCEtaWidth->at(ipho);
        phoVars[4] = phoSCPhiWidth->at(ipho);
        phoVars[5] = phoE3x3_2012->at(ipho);
        phoVars[6] = phoMaxEnergyXtal_2012->at(ipho);
        phoVars[7] = phoE2nd_2012->at(ipho);
        phoVars[8] = phoELeft_2012->at(ipho);
        phoVars[9] = phoERight_2012->at(ipho);
        phoVars[10] = phoETop_2012->at(ipho);
        phoVars[11] = phoEBottom_2012->at(ipho);
        phoVars[12] = phoSigmaIEtaIEta_2012->at(ipho);
        phoVars[13] = phoSigmaIEtaIPhi_2012->at(ipho);
        phoVars[14] = phoSigmaIPhiIPhi_2012->at(ipho);
        phoVars[15] = rho;
        phoVars[16] = phoESEn->at(ipho);

        corrected_Et = myPhoERegr.getCorrectedPt(phoVars, phoEt->at(ipho), phoEta->at(ipho), phoSCEta->at(ipho));

        Etmax = corrected_Et;
        pho_index = ipho;
        genin = pho_genMatchedIndex->at(ipho);

        train_weight = weight;
        train_weight_pthat = scale;
        train_weight_cent = Ncoll[hiBin];
        train_pthat = pthat;
        train_phoSCRawE = phoSCRawE->at(ipho);
        train_phoSCEta = phoSCEta->at(ipho);
        train_phoSCPhi = phoSCPhi->at(ipho);
        train_phoSCEtaWidth = phoSCEtaWidth->at(ipho);
        train_phoSCPhiWidth = phoSCPhiWidth->at(ipho);
        train_phoE3x3_2012 = phoE3x3_2012->at(ipho);
        train_phoMaxEnergyXtal_2012 = phoMaxEnergyXtal_2012->at(ipho);
        train_phoE2nd_2012 = phoE2nd_2012->at(ipho);
        train_phoELeft_2012 = phoELeft_2012->at(ipho);
        train_phoERight_2012 = phoERight_2012->at(ipho);
        train_phoETop_2012 = phoETop_2012->at(ipho);
        train_phoEBottom_2012 = phoEBottom_2012->at(ipho);
        train_phoSigmaIEtaIEta_2012 = phoSigmaIEtaIEta_2012->at(ipho);
        train_phoSigmaIEtaIPhi_2012 = phoSigmaIEtaIPhi_2012->at(ipho);
        train_phoSigmaIPhiIPhi_2012 = phoSigmaIPhiIPhi_2012->at(ipho);
        train_rho = rho;
        train_phoESEn = phoESEn->at(ipho);
        train_phoEt = phoEt->at(ipho);
        train_phoEtCorrected = Etmax;
        train_phoEta = phoEta->at(ipho);
        train_phoPhi = phoPhi->at(ipho);
        train_phoHoverE = phoHoverE->at(ipho);
        train_pho_ecalClusterIsoR3 = pho_ecalClusterIsoR3->at(ipho);
        train_pho_hcalRechitIsoR3 = pho_hcalRechitIsoR3->at(ipho);
        train_pho_trackIsoR3PtCut20 = pho_trackIsoR3PtCut20->at(ipho);
        train_pho_ecalClusterIsoR4 = pho_ecalClusterIsoR4->at(ipho);
        train_pho_hcalRechitIsoR4 = pho_hcalRechitIsoR4->at(ipho);
        train_pho_trackIsoR4PtCut20 = pho_trackIsoR4PtCut20->at(ipho);
        train_phoE = phoE->at(ipho);     
        train_phoSigmaEtaEta_2012 = phoSigmaEtaEta_2012->at(ipho);    
        train_pfcIso3subUE = pfcIso3subUE->at(ipho);
        train_pfnIso3subUE = pfnIso3subUE->at(ipho);
        train_pfpIso3subUE = pfpIso3subUE->at(ipho);
        train_SumPFIso = pfcIso3subUE->at(ipho)+pfnIso3subUE->at(ipho)+pfpIso3subUE->at(ipho);
        train_SumCalIso = pho_ecalClusterIsoR3->at(ipho)+pho_hcalRechitIsoR3->at(ipho)+pho_trackIsoR3PtCut20->at(ipho);
        train_hiBin = hiBin;
        train_L1_SingleEG21_BptxAND = L1_SingleEG21_BptxAND;
        train_HLT_HIGEDPhoton40_v1 = HLT_HIGEDPhoton40_v1;
        train_pho_genMatchedIndex = genin;
        if(genin!=-1){
        train_mcCalIsoDR04 = mcCalIsoDR04->at(genin);
        train_mcE = mcE->at(genin);
        train_mcEt = mcEt->at(genin);
        train_mcEta = mcEta->at(genin);
        train_mcPhi = mcPhi->at(genin);
        train_mcPID = mcPID->at(genin);
        train_mcMomPID = mcMomPID->at(genin);
        }
        else{
        train_mcCalIsoDR04 = -99999;
        train_mcE = -99999;
        train_mcEt = -99999;
        train_mcEta = -99999;
        train_mcPhi = -99999;
        train_mcPID = -99999;
        train_mcMomPID = -99999;
        }
        pho_tree->Fill();

      }


      // std::cout<<"i = "<< iEntry << " hiBin = "<<hiBin <<" Genin = "<<genin<<endl;

   } // End Event loop

   TFile *fout;
   fout = new TFile(type_dir + type+"_Kaya.root", "recreate");
   pho_tree->Write("",TObject::kOverwrite);

   fout->Close();

}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.042)                                             
{                                                                                                                                        
  TLatex* latex = new TLatex(x, y, text);                                                                                                
  latex->SetNDC();                                                                                                                       
  latex->SetTextSize(textSize);                                                                                                          
  latex->SetTextColor(color);                                                                                                            
  latex->SetTextFont(42);                                                                                                                
  latex->Draw();                                                                                                                         
}    


void SetPlotStyle()
{
   // from ATLAS plot style macro

   // use plain black on white colors
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetHistLineColor(1);

   gStyle->SetPalette(1);

   // set the paper & margin sizes
   gStyle->SetPaperSize(20, 26);
   gStyle->SetPadTopMargin(0.05);
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadBottomMargin(0.16);
   gStyle->SetPadLeftMargin(0.16);

   // set title offsets (for axis label)
   gStyle->SetTitleXOffset(1.4);
   gStyle->SetTitleYOffset(1.4);

   // use large fonts
   gStyle->SetTextFont(42);
   gStyle->SetTextSize(0.05);
   gStyle->SetLabelFont(42, "x");
   gStyle->SetTitleFont(42, "x");
   gStyle->SetLabelFont(42, "y");
   gStyle->SetTitleFont(42, "y");
   gStyle->SetLabelFont(42, "z");
   gStyle->SetTitleFont(42, "z");
   gStyle->SetLabelSize(0.05, "x");
   gStyle->SetTitleSize(0.05, "x");
   gStyle->SetLabelSize(0.05, "y");
   gStyle->SetTitleSize(0.05, "y");
   gStyle->SetLabelSize(0.05, "z");
   gStyle->SetTitleSize(0.05, "z");

   // use bold lines and markers
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(1.2);
   gStyle->SetHistLineWidth(2.);
   gStyle->SetLineStyleString(2, "[12 12]");

   // get rid of error bar caps
   gStyle->SetEndErrorSize(0.);

   // do not display any of the standard histogram decorations
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // put tick marks on top and RHS of plots
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
}

void displayProgress(long current, long max)
{
   using std::cerr;
   if (max < 2500)
      return;
   if (current % (max / 2500) != 0 && current < max - 1)
      return;

   int width = 52; // Hope the terminal is at least that wide.
   int barWidth = width - 2;
   cerr << "\x1B[2K";    // Clear line
   cerr << "\x1B[2000D"; // Cursor left
   cerr << '[';
   for (int i = 0; i < barWidth; ++i)
   {
      if (i < barWidth * current / max)
      {
         cerr << '=';
      }
      else
      {
         cerr << ' ';
      }
   }
   cerr << ']';
   cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
   cerr.flush();
}

int main(){
   SkimFile_MC_allPho_Bharad();
   return 0;
}
