#include <iostream>
#include <cstdlib>
#include <ctime>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "vector"
#include "map"

void sigShuffler() {
  std::vector<TFile*> inFiles;
  std::vector<TTree*> inTrees;

  std::vector<std::pair<unsigned short, unsigned long long> > eventList;
  
  unsigned short nInFiles = 40;
  std::string inFileNames[] =  {
    "culled_sig/culled_sig_QCD_HT500to700_0_0.root",
    "culled_sig/culled_sig_QCD_HT500to700_0_1.root",
    "culled_sig/culled_sig_QCD_HT500to700_0_2.root",
    "culled_sig/culled_sig_QCD_HT500to700_0_3.root",
    "culled_sig/culled_sig_QCD_HT500to700_1_0.root",
    "culled_sig/culled_sig_QCD_HT500to700_1_1.root",
    "culled_sig/culled_sig_QCD_HT500to700_1_2.root",
    "culled_sig/culled_sig_QCD_HT500to700_1_3.root",
    "culled_sig/culled_sig_QCD_HT500to700_2_0.root",
    "culled_sig/culled_sig_QCD_HT500to700_2_1.root",
    "culled_sig/culled_sig_QCD_HT500to700_2_2.root",
    "culled_sig/culled_sig_QCD_HT500to700_2_3.root",
    "culled_sig/culled_sig_QCD_HT500to700_3_0.root",
    "culled_sig/culled_sig_QCD_HT500to700_3_1.root",
    "culled_sig/culled_sig_QCD_HT500to700_3_2.root",
    "culled_sig/culled_sig_QCD_HT500to700_3_3.root",
    "culled_sig/culled_sig_QCD_HT700to1000_0_0.root",
    "culled_sig/culled_sig_QCD_HT700to1000_0_1.root",
    "culled_sig/culled_sig_QCD_HT700to1000_0_2.root",
    "culled_sig/culled_sig_QCD_HT700to1000_0_3.root",
    "culled_sig/culled_sig_QCD_HT700to1000_1_0.root",
    "culled_sig/culled_sig_QCD_HT700to1000_1_1.root",
    "culled_sig/culled_sig_QCD_HT700to1000_1_2.root",
    "culled_sig/culled_sig_QCD_HT700to1000_1_3.root",
    "culled_sig/culled_sig_QCD_HT700to1000_2_0.root",
    "culled_sig/culled_sig_QCD_HT700to1000_2_1.root",
    "culled_sig/culled_sig_QCD_HT700to1000_2_2.root",
    "culled_sig/culled_sig_QCD_HT700to1000_2_3.root",
    "culled_sig/culled_sig_QCD_HT1000to1500_0_0.root",
    "culled_sig/culled_sig_QCD_HT1000to1500_0_1.root",
    "culled_sig/culled_sig_QCD_HT1000to1500_0_2.root",
    "culled_sig/culled_sig_QCD_HT1000to1500_0_3.root",
    "culled_sig/culled_sig_QCD_HT1500to2000_0_0.root",
    "culled_sig/culled_sig_QCD_HT1500to2000_0_1.root",
    "culled_sig/culled_sig_QCD_HT1500to2000_0_2.root",
    "culled_sig/culled_sig_QCD_HT1500to2000_0_3.root",
    "culled_sig/culled_sig_QCD_HT2000toInf_0_0.root",
    "culled_sig/culled_sig_QCD_HT2000toInf_0_1.root",
    "culled_sig/culled_sig_QCD_HT2000toInf_0_2.root",
    "culled_sig/culled_sig_QCD_HT2000toInf_0_3.root",
  };
  // Set branch addresses
  const Int_t kMaxpassFilter_HBHE = 1;
  const Int_t kMaxpassFilter_HBHELoose = 1;
  const Int_t kMaxpassFilter_HBHETight = 1;
  const Int_t kMaxpassFilter_HBHEIso = 1;
  const Int_t kMaxpassFilter_CSCHalo = 1;
  const Int_t kMaxpassFilter_CSCTightHalo2015 = 1;
  const Int_t kMaxpassFilter_HCALlaser = 1;
  const Int_t kMaxpassFilter_ECALDeadCell = 1;
  const Int_t kMaxpassFilter_GoodVtx = 1;
  const Int_t kMaxpassFilter_TrkFailure = 1;
  const Int_t kMaxpassFilter_EEBadSc = 1;
  const Int_t kMaxpassFilter_ECALlaser = 1;
  const Int_t kMaxpassFilter_TrkPOG = 1;
  const Int_t kMaxpassFilter_TrkPOG_manystrip = 1;
  const Int_t kMaxpassFilter_TrkPOG_toomanystrip = 1;
  const Int_t kMaxpassFilter_TrkPOG_logError = 1;
  const Int_t kMaxpassFilter_METFilters = 1;
  const Int_t kMaxpassFilter_CSCTightHaloTrkMuUnvetoFilter = 1;
  const Int_t kMaxpassFilter_globalTightHalo2016 = 1;
  const Int_t kMaxpassFilter_HcalStripHalo = 1;
  const Int_t kMaxpassFilter_chargedHadronTrackResolution = 1;
  const Int_t kMaxpassFilter_muonBadTrack = 1;
  Int_t           ph_N;
  vector<int>     *ph_pdgId;
  vector<float>   *ph_charge;
  vector<float>   *ph_e;
  vector<float>   *ph_eta;
  vector<float>   *ph_phi;
  vector<float>   *ph_mass;
  vector<float>   *ph_pt;
  vector<float>   *ph_et;
  vector<float>   *ph_rho;
  vector<float>   *ph_superCluster_eta;
  vector<float>   *ph_superCluster_phi;
  vector<float>   *ph_sigmaIetaIeta;
  vector<float>   *ph_hOverE;
  vector<float>   *ph_isoGamma;
  vector<float>   *ph_isoCh;
  vector<bool>    *ph_passEleVeto;
  vector<int>     *ph_passLooseId;
  vector<int>     *ph_passMediumId;
  vector<int>     *ph_passTightId;
  vector<float>   *ph_mvaVal;
  vector<float>   *ph_mvaCat;
  Float_t         rho;
  Int_t           jetAK4_N;
  vector<float>   *jetAK4_pt;
  vector<float>   *jetAK4_eta;
  vector<float>   *jetAK4_mass;
  vector<float>   *jetAK4_phi;
  vector<float>   *jetAK4_e;
  vector<float>   *jetAK4_jec;
  vector<float>   *jetAK4_jecUp;
  vector<float>   *jetAK4_jecDown;
  vector<bool>    *jetAK4_IDLoose;
  vector<bool>    *jetAK4_IDTight;
  vector<bool>    *jetAK4_IDTightLepVeto;
  vector<int>     *jetAK4_charge;
  vector<float>   *jetAK4_csv;
  vector<float>   *jetAK4_vtxMass;
  vector<float>   *jetAK4_vtxNtracks;
  vector<float>   *jetAK4_vtx3DVal;
  vector<float>   *jetAK4_vtx3DSig;
  vector<int>     *jetAK4_partonFlavour;
  vector<int>     *jetAK4_hadronFlavour;
  vector<int>     *jetAK4_genParton_pdgID;
  vector<int>     *jetAK4_nbHadrons;
  vector<int>     *jetAK4_ncHadrons;
  vector<float>   *jetAK4_jer_sf;
  vector<float>   *jetAK4_jer_sf_up;
  vector<float>   *jetAK4_jer_sf_down;
  vector<float>   *jetAK4_jer_sigma_pt;
  Int_t           jetAK8_N;
  vector<float>   *jetAK8_pt;
  vector<float>   *jetAK8_eta;
  vector<float>   *jetAK8_mass;
  vector<float>   *jetAK8_phi;
  vector<float>   *jetAK8_e;
  vector<float>   *jetAK8_jec;
  vector<float>   *jetAK8_jecUp;
  vector<float>   *jetAK8_jecDown;
  vector<bool>    *jetAK8_IDLoose;
  vector<bool>    *jetAK8_IDTight;
  vector<bool>    *jetAK8_IDTightLepVeto;
  vector<int>     *jetAK8_charge;
  vector<int>     *jetAK8_partonFlavour;
  vector<int>     *jetAK8_hadronFlavour;
  vector<int>     *jetAK8_genParton_pdgID;
  vector<int>     *jetAK8_nbHadrons;
  vector<int>     *jetAK8_ncHadrons;
  vector<float>   *jetAK8_jer_sf;
  vector<float>   *jetAK8_jer_sf_up;
  vector<float>   *jetAK8_jer_sf_down;
  vector<float>   *jetAK8_jer_sigma_pt;
  vector<float>   *jetAK8Puppi_jer_sf;
  vector<float>   *jetAK8Puppi_jer_sf_up;
  vector<float>   *jetAK8Puppi_jer_sf_down;
  vector<float>   *jetAK8Puppi_jer_sigma_pt;
  vector<float>   *jetAK8_Hbbtag;
  vector<float>   *jetAK8_csv;
  vector<float>   *jetAK8_tau1;
  vector<float>   *jetAK8_tau2;
  vector<float>   *jetAK8_tau3;
  vector<float>   *jetAK8_pruned_mass;
  vector<float>   *jetAK8_pruned_massCorr;
  vector<float>   *jetAK8_pruned_jec;
  vector<float>   *jetAK8_pruned_jecUp;
  vector<float>   *jetAK8_pruned_jecDown;
  vector<float>   *jetAK8_softdrop_mass;
  vector<float>   *jetAK8_softdrop_massCorr;
  vector<float>   *jetAK8_softdrop_jec;
  vector<float>   *jetAK8_softdrop_jecUp;
  vector<float>   *jetAK8_softdrop_jecDown;
  vector<int>     *jetAK8_subjet_softdrop_N;
  vector<vector<float> > *jetAK8_subjet_softdrop_pt;
  vector<vector<float> > *jetAK8_subjet_softdrop_eta;
  vector<vector<float> > *jetAK8_subjet_softdrop_mass;
  vector<vector<float> > *jetAK8_subjet_softdrop_phi;
  vector<vector<float> > *jetAK8_subjet_softdrop_e;
  vector<vector<int> > *jetAK8_subjet_softdrop_charge;
  vector<vector<int> > *jetAK8_subjet_softdrop_genParton_pdgID;
  vector<vector<int> > *jetAK8_subjet_softdrop_nbHadrons;
  vector<vector<int> > *jetAK8_subjet_softdrop_ncHadrons;
  vector<vector<int> > *jetAK8_subjet_softdrop_partonFlavour;
  vector<vector<int> > *jetAK8_subjet_softdrop_hadronFlavour;
  vector<vector<float> > *jetAK8_subjet_softdrop_csv;
  vector<float>   *jetAK8_puppi_pt;
  vector<float>   *jetAK8_puppi_eta;
  vector<float>   *jetAK8_puppi_mass;
  vector<float>   *jetAK8_puppi_phi;
  vector<float>   *jetAK8_puppi_e;
  vector<float>   *jetAK8_puppi_pruned_mass;
  vector<float>   *jetAK8_puppi_pruned_massCorr;
  vector<float>   *jetAK8_puppi_pruned_jec;
  vector<float>   *jetAK8_puppi_softdrop_mass;
  vector<float>   *jetAK8_puppi_softdrop_massCorr;
  vector<float>   *jetAK8_puppi_softdrop_jec;
  vector<float>   *jetAK8_puppi_tau1;
  vector<float>   *jetAK8_puppi_tau2;
  vector<float>   *jetAK8_puppi_tau3;
  vector<int>     *jetAK8_subjet_puppi_softdrop_N;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_pt;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_eta;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_mass;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_phi;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_e;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_charge;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_genParton_pdgID;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_nbHadrons;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_ncHadrons;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_partonFlavour;
  vector<vector<int> > *jetAK8_subjet_puppi_softdrop_hadronFlavour;
  vector<vector<float> > *jetAK8_subjet_puppi_softdrop_csv;
  vector<int>     *jetAK8_subjet_pruned_N;
  vector<vector<float> > *jetAK8_subjet_pruned_pt;
  vector<vector<float> > *jetAK8_subjet_pruned_eta;
  vector<vector<float> > *jetAK8_subjet_pruned_mass;
  vector<vector<float> > *jetAK8_subjet_pruned_phi;
  vector<vector<float> > *jetAK8_subjet_pruned_e;
  vector<vector<int> > *jetAK8_subjet_pruned_charge;
  vector<vector<int> > *jetAK8_subjet_pruned_genParton_pdgID;
  vector<vector<int> > *jetAK8_subjet_pruned_nbHadrons;
  vector<vector<int> > *jetAK8_subjet_pruned_ncHadrons;
  vector<vector<int> > *jetAK8_subjet_pruned_partonFlavour;
  vector<vector<int> > *jetAK8_subjet_pruned_hadronFlavour;
  vector<vector<float> > *jetAK8_subjet_pruned_csv;
  map<string,bool> *HLT_isFired;
  Bool_t          passFilter_HBHE;
  Bool_t          passFilter_HBHELoose;
  Bool_t          passFilter_HBHETight;
  Bool_t          passFilter_HBHEIso;
  Bool_t          passFilter_CSCHalo;
  Bool_t          passFilter_CSCTightHalo2015;
  Bool_t          passFilter_HCALlaser;
  Bool_t          passFilter_ECALDeadCell;
  Bool_t          passFilter_GoodVtx;
  Bool_t          passFilter_TrkFailure;
  Bool_t          passFilter_EEBadSc;
  Bool_t          passFilter_ECALlaser;
  Bool_t          passFilter_TrkPOG;
  Bool_t          passFilter_TrkPOG_manystrip;
  Bool_t          passFilter_TrkPOG_toomanystrip;
  Bool_t          passFilter_TrkPOG_logError;
  Bool_t          passFilter_METFilters;
  Bool_t          passFilter_CSCTightHaloTrkMuUnvetoFilter;
  Bool_t          passFilter_globalTightHalo2016;
  Bool_t          passFilter_HcalStripHalo;
  Bool_t          passFilter_chargedHadronTrackResolution;
  Bool_t          passFilter_muonBadTrack;
  Int_t           EVENT_event;
  Int_t           EVENT_run;
  Int_t           EVENT_lumiBlock;
  Int_t           PV_N;
  Bool_t          PV_filter;
  vector<float>   *PV_chi2;
  vector<float>   *PV_ndof;
  vector<float>   *PV_rho;
  vector<float>   *PV_z;
  TBranch        *b_ph_N;   //!
  TBranch        *b_ph_pdgId;   //!
  TBranch        *b_ph_charge;   //!
  TBranch        *b_ph_e;   //!
  TBranch        *b_ph_eta;   //!
  TBranch        *b_ph_phi;   //!
  TBranch        *b_ph_mass;   //!
  TBranch        *b_ph_pt;   //!
  TBranch        *b_ph_et;   //!
  TBranch        *b_ph_rho;   //!
  TBranch        *b_ph_superCluster_eta;   //!
  TBranch        *b_ph_superCluster_phi;   //!
  TBranch        *b_ph_sigmaIetaIeta;   //!
  TBranch        *b_ph_hOverE;   //!
  TBranch        *b_ph_isoGamma;   //!
  TBranch        *b_ph_isoCh;   //!
  TBranch        *b_ph_passEleVeto;   //!
  TBranch        *b_ph_passLooseId;   //!
  TBranch        *b_ph_passMediumId;   //!
  TBranch        *b_ph_passTightId;   //!
  TBranch        *b_ph_mvaVal;   //!
  TBranch        *b_ph_mvaCat;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_jetAK4_N;   //!
  TBranch        *b_jetAK4_pt;   //!
  TBranch        *b_jetAK4_eta;   //!
  TBranch        *b_jetAK4_mass;   //!
  TBranch        *b_jetAK4_phi;   //!
  TBranch        *b_jetAK4_e;   //!
  TBranch        *b_jetAK4_jec;   //!
  TBranch        *b_jetAK4_jecUp;   //!
  TBranch        *b_jetAK4_jecDown;   //!
  TBranch        *b_jetAK4_IDLoose;   //!
  TBranch        *b_jetAK4_IDTight;   //!
  TBranch        *b_jetAK4_IDTightLepVeto;   //!
  TBranch        *b_jetAK4_charge;   //!
  TBranch        *b_jetAK4_csv;   //!
  TBranch        *b_jetAK4_vtxMass;   //!
  TBranch        *b_jetAK4_vtxNtracks;   //!
  TBranch        *b_jetAK4_vtx3DVal;   //!
  TBranch        *b_jetAK4_vtx3DSig;   //!
  TBranch        *b_jetAK4_partonFlavour;   //!
  TBranch        *b_jetAK4_hadronFlavour;   //!
  TBranch        *b_jetAK4_genParton_pdgID;   //!
  TBranch        *b_jetAK4_nbHadrons;   //!
  TBranch        *b_jetAK4_ncHadrons;   //!
  TBranch        *b_jetAK4_jer_sf;   //!
  TBranch        *b_jetAK4_jer_sf_up;   //!
  TBranch        *b_jetAK4_jer_sf_down;   //!
  TBranch        *b_jetAK4_jer_sigma_pt;   //!
  TBranch        *b_jetAK8_N;   //!
  TBranch        *b_jetAK8_pt;   //!
  TBranch        *b_jetAK8_eta;   //!
  TBranch        *b_jetAK8_mass;   //!
  TBranch        *b_jetAK8_phi;   //!
  TBranch        *b_jetAK8_e;   //!
  TBranch        *b_jetAK8_jec;   //!
  TBranch        *b_jetAK8_jecUp;   //!
  TBranch        *b_jetAK8_jecDown;   //!
  TBranch        *b_jetAK8_IDLoose;   //!
  TBranch        *b_jetAK8_IDTight;   //!
  TBranch        *b_jetAK8_IDTightLepVeto;   //!
  TBranch        *b_jetAK8_charge;   //!
  TBranch        *b_jetAK8_partonFlavour;   //!
  TBranch        *b_jetAK8_hadronFlavour;   //!
  TBranch        *b_jetAK8_genParton_pdgID;   //!
  TBranch        *b_jetAK8_nbHadrons;   //!
  TBranch        *b_jetAK8_ncHadrons;   //!
  TBranch        *b_jetAK8_jer_sf;   //!
  TBranch        *b_jetAK8_jer_sf_up;   //!
  TBranch        *b_jetAK8_jer_sf_down;   //!
  TBranch        *b_jetAK8_jer_sigma_pt;   //!
  TBranch        *b_jetAK8Puppi_jer_sf;   //!
  TBranch        *b_jetAK8Puppi_jer_sf_up;   //!
  TBranch        *b_jetAK8Puppi_jer_sf_down;   //!
  TBranch        *b_jetAK8Puppi_jer_sigma_pt;   //!
  TBranch        *b_jetAK8_Hbbtag;   //!
  TBranch        *b_jetAK8_csv;   //!
  TBranch        *b_jetAK8_tau1;   //!
  TBranch        *b_jetAK8_tau2;   //!
  TBranch        *b_jetAK8_tau3;   //!
  TBranch        *b_jetAK8_pruned_mass;   //!
  TBranch        *b_jetAK8_pruned_massCorr;   //!
  TBranch        *b_jetAK8_pruned_jec;   //!
  TBranch        *b_jetAK8_pruned_jecUp;   //!
  TBranch        *b_jetAK8_pruned_jecDown;   //!
  TBranch        *b_jetAK8_softdrop_mass;   //!
  TBranch        *b_jetAK8_softdrop_massCorr;   //!
  TBranch        *b_jetAK8_softdrop_jec;   //!
  TBranch        *b_jetAK8_softdrop_jecUp;   //!
  TBranch        *b_jetAK8_softdrop_jecDown;   //!
  TBranch        *b_jetAK8_subjet_softdrop_N;   //!
  TBranch        *b_jetAK8_subjet_softdrop_pt;   //!
  TBranch        *b_jetAK8_subjet_softdrop_eta;   //!
  TBranch        *b_jetAK8_subjet_softdrop_mass;   //!
  TBranch        *b_jetAK8_subjet_softdrop_phi;   //!
  TBranch        *b_jetAK8_subjet_softdrop_e;   //!
  TBranch        *b_jetAK8_subjet_softdrop_charge;   //!
  TBranch        *b_jetAK8_subjet_softdrop_genParton_pdgID;   //!
  TBranch        *b_jetAK8_subjet_softdrop_nbHadrons;   //!
  TBranch        *b_jetAK8_subjet_softdrop_ncHadrons;   //!
  TBranch        *b_jetAK8_subjet_softdrop_partonFlavour;   //!
  TBranch        *b_jetAK8_subjet_softdrop_hadronFlavour;   //!
  TBranch        *b_jetAK8_subjet_softdrop_csv;   //!
  TBranch        *b_jetAK8_puppi_pt;   //!
  TBranch        *b_jetAK8_puppi_eta;   //!
  TBranch        *b_jetAK8_puppi_mass;   //!
  TBranch        *b_jetAK8_puppi_phi;   //!
  TBranch        *b_jetAK8_puppi_e;   //!
  TBranch        *b_jetAK8_puppi_pruned_mass;   //!
  TBranch        *b_jetAK8_puppi_pruned_massCorr;   //!
  TBranch        *b_jetAK8_puppi_pruned_jec;   //!
  TBranch        *b_jetAK8_puppi_softdrop_mass;   //!
  TBranch        *b_jetAK8_puppi_softdrop_massCorr;   //!
  TBranch        *b_jetAK8_puppi_softdrop_jec;   //!
  TBranch        *b_jetAK8_puppi_tau1;   //!
  TBranch        *b_jetAK8_puppi_tau2;   //!
  TBranch        *b_jetAK8_puppi_tau3;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_N;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_pt;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_eta;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_mass;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_phi;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_e;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_charge;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_genParton_pdgID;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_nbHadrons;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_ncHadrons;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_partonFlavour;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_hadronFlavour;   //!
  TBranch        *b_jetAK8_subjet_puppi_softdrop_csv;   //!
  TBranch        *b_jetAK8_subjet_pruned_N;   //!
  TBranch        *b_jetAK8_subjet_pruned_pt;   //!
  TBranch        *b_jetAK8_subjet_pruned_eta;   //!
  TBranch        *b_jetAK8_subjet_pruned_mass;   //!
  TBranch        *b_jetAK8_subjet_pruned_phi;   //!
  TBranch        *b_jetAK8_subjet_pruned_e;   //!
  TBranch        *b_jetAK8_subjet_pruned_charge;   //!
  TBranch        *b_jetAK8_subjet_pruned_genParton_pdgID;   //!
  TBranch        *b_jetAK8_subjet_pruned_nbHadrons;   //!
  TBranch        *b_jetAK8_subjet_pruned_ncHadrons;   //!
  TBranch        *b_jetAK8_subjet_pruned_partonFlavour;   //!
  TBranch        *b_jetAK8_subjet_pruned_hadronFlavour;   //!
  TBranch        *b_jetAK8_subjet_pruned_csv;   //!
  TBranch        *b_HLT_isFired;   //!
  TBranch        *b_passFilter_HBHE_;   //!
  TBranch        *b_passFilter_HBHELoose_;   //!
  TBranch        *b_passFilter_HBHETight_;   //!
  TBranch        *b_passFilter_HBHEIso_;   //!
  TBranch        *b_passFilter_CSCHalo_;   //!
  TBranch        *b_passFilter_CSCTightHalo2015_;   //!
  TBranch        *b_passFilter_HCALlaser_;   //!
  TBranch        *b_passFilter_ECALDeadCell_;   //!
  TBranch        *b_passFilter_GoodVtx_;   //!
  TBranch        *b_passFilter_TrkFailure_;   //!
  TBranch        *b_passFilter_EEBadSc_;   //!
  TBranch        *b_passFilter_ECALlaser_;   //!
  TBranch        *b_passFilter_TrkPOG_;   //!
  TBranch        *b_passFilter_TrkPOG_manystrip_;   //!
  TBranch        *b_passFilter_TrkPOG_toomanystrip_;   //!
  TBranch        *b_passFilter_TrkPOG_logError_;   //!
  TBranch        *b_passFilter_METFilters_;   //!
  TBranch        *b_passFilter_CSCTightHaloTrkMuUnvetoFilter_;   //!
  TBranch        *b_passFilter_globalTightHalo2016_;   //!
  TBranch        *b_passFilter_HcalStripHalo_;   //!
  TBranch        *b_passFilter_chargedHadronTrackResolution_;   //!
  TBranch        *b_passFilter_muonBadTrack_;   //!
  TBranch        *b_EVENT_event;   //!
  TBranch        *b_EVENT_run;   //!
  TBranch        *b_EVENT_lumiBlock;   //!
  TBranch        *b_PV_N;   //!
  TBranch        *b_PV_filter;   //!
  TBranch        *b_PV_chi2;   //!
  TBranch        *b_PV_ndof;   //!
  TBranch        *b_PV_rho;   //!
  TBranch        *b_PV_z;   //!

  for (unsigned short i=0; i<nInFiles; ++i) { 
    inFiles.push_back(new TFile(inFileNames[i].c_str()));
    inFiles.back()->SetCompressionLevel(0);
    inTrees.push_back((TTree*) inFiles.back()->Get("tree"));
    inTrees.back()->LoadBaskets();
    cout << "inFile " << i << " is: " << inFiles.at(i)->GetName() << endl;
    cout << "inTree " << i << " is: " << inTrees.at(i) << endl;
    cout << "inTree " << i << " nEntries: " << inTrees.at(i)->GetEntries() << endl;
    
    for (Long64_t iEvent=0; iEvent < inTrees.at(i)->GetEntries(); ++iEvent) {
      eventList.push_back(std::pair<unsigned short, unsigned long long>(i, iEvent ));
    }
   // Set object pointer
   ph_pdgId = 0;
   ph_charge = 0;
   ph_e = 0;
   ph_eta = 0;
   ph_phi = 0;
   ph_mass = 0;
   ph_pt = 0;
   ph_et = 0;
   ph_rho = 0;
   ph_superCluster_eta = 0;
   ph_superCluster_phi = 0;
   ph_sigmaIetaIeta = 0;
   ph_hOverE = 0;
   ph_isoGamma = 0;
   ph_isoCh = 0;
   ph_passEleVeto = 0;
   ph_passLooseId = 0;
   ph_passMediumId = 0;
   ph_passTightId = 0;
   ph_mvaVal = 0;
   ph_mvaCat = 0;
   jetAK4_pt = 0;
   jetAK4_eta = 0;
   jetAK4_mass = 0;
   jetAK4_phi = 0;
   jetAK4_e = 0;
   jetAK4_jec = 0;
   jetAK4_jecUp = 0;
   jetAK4_jecDown = 0;
   jetAK4_IDLoose = 0;
   jetAK4_IDTight = 0;
   jetAK4_IDTightLepVeto = 0;
   jetAK4_charge = 0;
   jetAK4_csv = 0;
   jetAK4_vtxMass = 0;
   jetAK4_vtxNtracks = 0;
   jetAK4_vtx3DVal = 0;
   jetAK4_vtx3DSig = 0;
   jetAK4_partonFlavour = 0;
   jetAK4_hadronFlavour = 0;
   jetAK4_genParton_pdgID = 0;
   jetAK4_nbHadrons = 0;
   jetAK4_ncHadrons = 0;
   jetAK4_jer_sf = 0;
   jetAK4_jer_sf_up = 0;
   jetAK4_jer_sf_down = 0;
   jetAK4_jer_sigma_pt = 0;
   jetAK8_pt = 0;
   jetAK8_eta = 0;
   jetAK8_mass = 0;
   jetAK8_phi = 0;
   jetAK8_e = 0;
   jetAK8_jec = 0;
   jetAK8_jecUp = 0;
   jetAK8_jecDown = 0;
   jetAK8_IDLoose = 0;
   jetAK8_IDTight = 0;
   jetAK8_IDTightLepVeto = 0;
   jetAK8_charge = 0;
   jetAK8_partonFlavour = 0;
   jetAK8_hadronFlavour = 0;
   jetAK8_genParton_pdgID = 0;
   jetAK8_nbHadrons = 0;
   jetAK8_ncHadrons = 0;
   jetAK8_jer_sf = 0;
   jetAK8_jer_sf_up = 0;
   jetAK8_jer_sf_down = 0;
   jetAK8_jer_sigma_pt = 0;
   jetAK8Puppi_jer_sf = 0;
   jetAK8Puppi_jer_sf_up = 0;
   jetAK8Puppi_jer_sf_down = 0;
   jetAK8Puppi_jer_sigma_pt = 0;
   jetAK8_Hbbtag = 0;
   jetAK8_csv = 0;
   jetAK8_tau1 = 0;
   jetAK8_tau2 = 0;
   jetAK8_tau3 = 0;
   jetAK8_pruned_mass = 0;
   jetAK8_pruned_massCorr = 0;
   jetAK8_pruned_jec = 0;
   jetAK8_pruned_jecUp = 0;
   jetAK8_pruned_jecDown = 0;
   jetAK8_softdrop_mass = 0;
   jetAK8_softdrop_massCorr = 0;
   jetAK8_softdrop_jec = 0;
   jetAK8_softdrop_jecUp = 0;
   jetAK8_softdrop_jecDown = 0;
   jetAK8_subjet_softdrop_N = 0;
   jetAK8_subjet_softdrop_pt = 0;
   jetAK8_subjet_softdrop_eta = 0;
   jetAK8_subjet_softdrop_mass = 0;
   jetAK8_subjet_softdrop_phi = 0;
   jetAK8_subjet_softdrop_e = 0;
   jetAK8_subjet_softdrop_charge = 0;
   jetAK8_subjet_softdrop_genParton_pdgID = 0;
   jetAK8_subjet_softdrop_nbHadrons = 0;
   jetAK8_subjet_softdrop_ncHadrons = 0;
   jetAK8_subjet_softdrop_partonFlavour = 0;
   jetAK8_subjet_softdrop_hadronFlavour = 0;
   jetAK8_subjet_softdrop_csv = 0;
   jetAK8_puppi_pt = 0;
   jetAK8_puppi_eta = 0;
   jetAK8_puppi_mass = 0;
   jetAK8_puppi_phi = 0;
   jetAK8_puppi_e = 0;
   jetAK8_puppi_pruned_mass = 0;
   jetAK8_puppi_pruned_massCorr = 0;
   jetAK8_puppi_pruned_jec = 0;
   jetAK8_puppi_softdrop_mass = 0;
   jetAK8_puppi_softdrop_massCorr = 0;
   jetAK8_puppi_softdrop_jec = 0;
   jetAK8_puppi_tau1 = 0;
   jetAK8_puppi_tau2 = 0;
   jetAK8_puppi_tau3 = 0;
   jetAK8_subjet_puppi_softdrop_N = 0;
   jetAK8_subjet_puppi_softdrop_pt = 0;
   jetAK8_subjet_puppi_softdrop_eta = 0;
   jetAK8_subjet_puppi_softdrop_mass = 0;
   jetAK8_subjet_puppi_softdrop_phi = 0;
   jetAK8_subjet_puppi_softdrop_e = 0;
   jetAK8_subjet_puppi_softdrop_charge = 0;
   jetAK8_subjet_puppi_softdrop_genParton_pdgID = 0;
   jetAK8_subjet_puppi_softdrop_nbHadrons = 0;
   jetAK8_subjet_puppi_softdrop_ncHadrons = 0;
   jetAK8_subjet_puppi_softdrop_partonFlavour = 0;
   jetAK8_subjet_puppi_softdrop_hadronFlavour = 0;
   jetAK8_subjet_puppi_softdrop_csv = 0;
   jetAK8_subjet_pruned_N = 0;
   jetAK8_subjet_pruned_pt = 0;
   jetAK8_subjet_pruned_eta = 0;
   jetAK8_subjet_pruned_mass = 0;
   jetAK8_subjet_pruned_phi = 0;
   jetAK8_subjet_pruned_e = 0;
   jetAK8_subjet_pruned_charge = 0;
   jetAK8_subjet_pruned_genParton_pdgID = 0;
   jetAK8_subjet_pruned_nbHadrons = 0;
   jetAK8_subjet_pruned_ncHadrons = 0;
   jetAK8_subjet_pruned_partonFlavour = 0;
   jetAK8_subjet_pruned_hadronFlavour = 0;
   jetAK8_subjet_pruned_csv = 0;
   HLT_isFired = 0;
   PV_chi2 = 0;
   PV_ndof = 0;
   PV_rho = 0;
   PV_z = 0;
   inTrees.back()->SetBranchAddress("ph_N", &ph_N, &b_ph_N);
   inTrees.back()->SetBranchAddress("ph_pdgId", &ph_pdgId, &b_ph_pdgId);
   inTrees.back()->SetBranchAddress("ph_charge", &ph_charge, &b_ph_charge);
   inTrees.back()->SetBranchAddress("ph_e", &ph_e, &b_ph_e);
   inTrees.back()->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   inTrees.back()->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   inTrees.back()->SetBranchAddress("ph_mass", &ph_mass, &b_ph_mass);
   inTrees.back()->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   inTrees.back()->SetBranchAddress("ph_et", &ph_et, &b_ph_et);
   inTrees.back()->SetBranchAddress("ph_rho", &ph_rho, &b_ph_rho);
   inTrees.back()->SetBranchAddress("ph_superCluster_eta", &ph_superCluster_eta, &b_ph_superCluster_eta);
   inTrees.back()->SetBranchAddress("ph_superCluster_phi", &ph_superCluster_phi, &b_ph_superCluster_phi);
   inTrees.back()->SetBranchAddress("ph_sigmaIetaIeta", &ph_sigmaIetaIeta, &b_ph_sigmaIetaIeta);
   inTrees.back()->SetBranchAddress("ph_hOverE", &ph_hOverE, &b_ph_hOverE);
   inTrees.back()->SetBranchAddress("ph_isoGamma", &ph_isoGamma, &b_ph_isoGamma);
   inTrees.back()->SetBranchAddress("ph_isoCh", &ph_isoCh, &b_ph_isoCh);
   inTrees.back()->SetBranchAddress("ph_passEleVeto", &ph_passEleVeto, &b_ph_passEleVeto);
   inTrees.back()->SetBranchAddress("ph_passLooseId", &ph_passLooseId, &b_ph_passLooseId);
   inTrees.back()->SetBranchAddress("ph_passMediumId", &ph_passMediumId, &b_ph_passMediumId);
   inTrees.back()->SetBranchAddress("ph_passTightId", &ph_passTightId, &b_ph_passTightId);
   inTrees.back()->SetBranchAddress("ph_mvaVal", &ph_mvaVal, &b_ph_mvaVal);
   inTrees.back()->SetBranchAddress("ph_mvaCat", &ph_mvaCat, &b_ph_mvaCat);
   inTrees.back()->SetBranchAddress("rho", &rho, &b_rho);
   inTrees.back()->SetBranchAddress("jetAK4_N", &jetAK4_N, &b_jetAK4_N);
   inTrees.back()->SetBranchAddress("jetAK4_pt", &jetAK4_pt, &b_jetAK4_pt);
   inTrees.back()->SetBranchAddress("jetAK4_eta", &jetAK4_eta, &b_jetAK4_eta);
   inTrees.back()->SetBranchAddress("jetAK4_mass", &jetAK4_mass, &b_jetAK4_mass);
   inTrees.back()->SetBranchAddress("jetAK4_phi", &jetAK4_phi, &b_jetAK4_phi);
   inTrees.back()->SetBranchAddress("jetAK4_e", &jetAK4_e, &b_jetAK4_e);
   inTrees.back()->SetBranchAddress("jetAK4_jec", &jetAK4_jec, &b_jetAK4_jec);
   inTrees.back()->SetBranchAddress("jetAK4_jecUp", &jetAK4_jecUp, &b_jetAK4_jecUp);
   inTrees.back()->SetBranchAddress("jetAK4_jecDown", &jetAK4_jecDown, &b_jetAK4_jecDown);
   inTrees.back()->SetBranchAddress("jetAK4_IDLoose", &jetAK4_IDLoose, &b_jetAK4_IDLoose);
   inTrees.back()->SetBranchAddress("jetAK4_IDTight", &jetAK4_IDTight, &b_jetAK4_IDTight);
   inTrees.back()->SetBranchAddress("jetAK4_IDTightLepVeto", &jetAK4_IDTightLepVeto, &b_jetAK4_IDTightLepVeto);
   inTrees.back()->SetBranchAddress("jetAK4_charge", &jetAK4_charge, &b_jetAK4_charge);
   inTrees.back()->SetBranchAddress("jetAK4_csv", &jetAK4_csv, &b_jetAK4_csv);
   inTrees.back()->SetBranchAddress("jetAK4_vtxMass", &jetAK4_vtxMass, &b_jetAK4_vtxMass);
   inTrees.back()->SetBranchAddress("jetAK4_vtxNtracks", &jetAK4_vtxNtracks, &b_jetAK4_vtxNtracks);
   inTrees.back()->SetBranchAddress("jetAK4_vtx3DVal", &jetAK4_vtx3DVal, &b_jetAK4_vtx3DVal);
   inTrees.back()->SetBranchAddress("jetAK4_vtx3DSig", &jetAK4_vtx3DSig, &b_jetAK4_vtx3DSig);
   inTrees.back()->SetBranchAddress("jetAK4_partonFlavour", &jetAK4_partonFlavour, &b_jetAK4_partonFlavour);
   inTrees.back()->SetBranchAddress("jetAK4_hadronFlavour", &jetAK4_hadronFlavour, &b_jetAK4_hadronFlavour);
   inTrees.back()->SetBranchAddress("jetAK4_genParton_pdgID", &jetAK4_genParton_pdgID, &b_jetAK4_genParton_pdgID);
   inTrees.back()->SetBranchAddress("jetAK4_nbHadrons", &jetAK4_nbHadrons, &b_jetAK4_nbHadrons);
   inTrees.back()->SetBranchAddress("jetAK4_ncHadrons", &jetAK4_ncHadrons, &b_jetAK4_ncHadrons);
   inTrees.back()->SetBranchAddress("jetAK4_jer_sf", &jetAK4_jer_sf, &b_jetAK4_jer_sf);
   inTrees.back()->SetBranchAddress("jetAK4_jer_sf_up", &jetAK4_jer_sf_up, &b_jetAK4_jer_sf_up);
   inTrees.back()->SetBranchAddress("jetAK4_jer_sf_down", &jetAK4_jer_sf_down, &b_jetAK4_jer_sf_down);
   inTrees.back()->SetBranchAddress("jetAK4_jer_sigma_pt", &jetAK4_jer_sigma_pt, &b_jetAK4_jer_sigma_pt);
   inTrees.back()->SetBranchAddress("jetAK8_N", &jetAK8_N, &b_jetAK8_N);
   inTrees.back()->SetBranchAddress("jetAK8_pt", &jetAK8_pt, &b_jetAK8_pt);
   inTrees.back()->SetBranchAddress("jetAK8_eta", &jetAK8_eta, &b_jetAK8_eta);
   inTrees.back()->SetBranchAddress("jetAK8_mass", &jetAK8_mass, &b_jetAK8_mass);
   inTrees.back()->SetBranchAddress("jetAK8_phi", &jetAK8_phi, &b_jetAK8_phi);
   inTrees.back()->SetBranchAddress("jetAK8_e", &jetAK8_e, &b_jetAK8_e);
   inTrees.back()->SetBranchAddress("jetAK8_jec", &jetAK8_jec, &b_jetAK8_jec);
   inTrees.back()->SetBranchAddress("jetAK8_jecUp", &jetAK8_jecUp, &b_jetAK8_jecUp);
   inTrees.back()->SetBranchAddress("jetAK8_jecDown", &jetAK8_jecDown, &b_jetAK8_jecDown);
   inTrees.back()->SetBranchAddress("jetAK8_IDLoose", &jetAK8_IDLoose, &b_jetAK8_IDLoose);
   inTrees.back()->SetBranchAddress("jetAK8_IDTight", &jetAK8_IDTight, &b_jetAK8_IDTight);
   inTrees.back()->SetBranchAddress("jetAK8_IDTightLepVeto", &jetAK8_IDTightLepVeto, &b_jetAK8_IDTightLepVeto);
   inTrees.back()->SetBranchAddress("jetAK8_charge", &jetAK8_charge, &b_jetAK8_charge);
   inTrees.back()->SetBranchAddress("jetAK8_partonFlavour", &jetAK8_partonFlavour, &b_jetAK8_partonFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_hadronFlavour", &jetAK8_hadronFlavour, &b_jetAK8_hadronFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_genParton_pdgID", &jetAK8_genParton_pdgID, &b_jetAK8_genParton_pdgID);
   inTrees.back()->SetBranchAddress("jetAK8_nbHadrons", &jetAK8_nbHadrons, &b_jetAK8_nbHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_ncHadrons", &jetAK8_ncHadrons, &b_jetAK8_ncHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_jer_sf", &jetAK8_jer_sf, &b_jetAK8_jer_sf);
   inTrees.back()->SetBranchAddress("jetAK8_jer_sf_up", &jetAK8_jer_sf_up, &b_jetAK8_jer_sf_up);
   inTrees.back()->SetBranchAddress("jetAK8_jer_sf_down", &jetAK8_jer_sf_down, &b_jetAK8_jer_sf_down);
   inTrees.back()->SetBranchAddress("jetAK8_jer_sigma_pt", &jetAK8_jer_sigma_pt, &b_jetAK8_jer_sigma_pt);
   inTrees.back()->SetBranchAddress("jetAK8Puppi_jer_sf", &jetAK8Puppi_jer_sf, &b_jetAK8Puppi_jer_sf);
   inTrees.back()->SetBranchAddress("jetAK8Puppi_jer_sf_up", &jetAK8Puppi_jer_sf_up, &b_jetAK8Puppi_jer_sf_up);
   inTrees.back()->SetBranchAddress("jetAK8Puppi_jer_sf_down", &jetAK8Puppi_jer_sf_down, &b_jetAK8Puppi_jer_sf_down);
   inTrees.back()->SetBranchAddress("jetAK8Puppi_jer_sigma_pt", &jetAK8Puppi_jer_sigma_pt, &b_jetAK8Puppi_jer_sigma_pt);
   inTrees.back()->SetBranchAddress("jetAK8_Hbbtag", &jetAK8_Hbbtag, &b_jetAK8_Hbbtag);
   inTrees.back()->SetBranchAddress("jetAK8_csv", &jetAK8_csv, &b_jetAK8_csv);
   inTrees.back()->SetBranchAddress("jetAK8_tau1", &jetAK8_tau1, &b_jetAK8_tau1);
   inTrees.back()->SetBranchAddress("jetAK8_tau2", &jetAK8_tau2, &b_jetAK8_tau2);
   inTrees.back()->SetBranchAddress("jetAK8_tau3", &jetAK8_tau3, &b_jetAK8_tau3);
   inTrees.back()->SetBranchAddress("jetAK8_pruned_mass", &jetAK8_pruned_mass, &b_jetAK8_pruned_mass);
   inTrees.back()->SetBranchAddress("jetAK8_pruned_massCorr", &jetAK8_pruned_massCorr, &b_jetAK8_pruned_massCorr);
   inTrees.back()->SetBranchAddress("jetAK8_pruned_jec", &jetAK8_pruned_jec, &b_jetAK8_pruned_jec);
   inTrees.back()->SetBranchAddress("jetAK8_pruned_jecUp", &jetAK8_pruned_jecUp, &b_jetAK8_pruned_jecUp);
   inTrees.back()->SetBranchAddress("jetAK8_pruned_jecDown", &jetAK8_pruned_jecDown, &b_jetAK8_pruned_jecDown);
   inTrees.back()->SetBranchAddress("jetAK8_softdrop_mass", &jetAK8_softdrop_mass, &b_jetAK8_softdrop_mass);
   inTrees.back()->SetBranchAddress("jetAK8_softdrop_massCorr", &jetAK8_softdrop_massCorr, &b_jetAK8_softdrop_massCorr);
   inTrees.back()->SetBranchAddress("jetAK8_softdrop_jec", &jetAK8_softdrop_jec, &b_jetAK8_softdrop_jec);
   inTrees.back()->SetBranchAddress("jetAK8_softdrop_jecUp", &jetAK8_softdrop_jecUp, &b_jetAK8_softdrop_jecUp);
   inTrees.back()->SetBranchAddress("jetAK8_softdrop_jecDown", &jetAK8_softdrop_jecDown, &b_jetAK8_softdrop_jecDown);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_N", &jetAK8_subjet_softdrop_N, &b_jetAK8_subjet_softdrop_N);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_pt", &jetAK8_subjet_softdrop_pt, &b_jetAK8_subjet_softdrop_pt);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_eta", &jetAK8_subjet_softdrop_eta, &b_jetAK8_subjet_softdrop_eta);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_mass", &jetAK8_subjet_softdrop_mass, &b_jetAK8_subjet_softdrop_mass);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_phi", &jetAK8_subjet_softdrop_phi, &b_jetAK8_subjet_softdrop_phi);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_e", &jetAK8_subjet_softdrop_e, &b_jetAK8_subjet_softdrop_e);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_charge", &jetAK8_subjet_softdrop_charge, &b_jetAK8_subjet_softdrop_charge);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_genParton_pdgID", &jetAK8_subjet_softdrop_genParton_pdgID, &b_jetAK8_subjet_softdrop_genParton_pdgID);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_nbHadrons", &jetAK8_subjet_softdrop_nbHadrons, &b_jetAK8_subjet_softdrop_nbHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_ncHadrons", &jetAK8_subjet_softdrop_ncHadrons, &b_jetAK8_subjet_softdrop_ncHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_partonFlavour", &jetAK8_subjet_softdrop_partonFlavour, &b_jetAK8_subjet_softdrop_partonFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_hadronFlavour", &jetAK8_subjet_softdrop_hadronFlavour, &b_jetAK8_subjet_softdrop_hadronFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_softdrop_csv", &jetAK8_subjet_softdrop_csv, &b_jetAK8_subjet_softdrop_csv);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_pt", &jetAK8_puppi_pt, &b_jetAK8_puppi_pt);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_eta", &jetAK8_puppi_eta, &b_jetAK8_puppi_eta);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_mass", &jetAK8_puppi_mass, &b_jetAK8_puppi_mass);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_phi", &jetAK8_puppi_phi, &b_jetAK8_puppi_phi);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_e", &jetAK8_puppi_e, &b_jetAK8_puppi_e);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_pruned_mass", &jetAK8_puppi_pruned_mass, &b_jetAK8_puppi_pruned_mass);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_pruned_massCorr", &jetAK8_puppi_pruned_massCorr, &b_jetAK8_puppi_pruned_massCorr);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_pruned_jec", &jetAK8_puppi_pruned_jec, &b_jetAK8_puppi_pruned_jec);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_softdrop_mass", &jetAK8_puppi_softdrop_mass, &b_jetAK8_puppi_softdrop_mass);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_softdrop_massCorr", &jetAK8_puppi_softdrop_massCorr, &b_jetAK8_puppi_softdrop_massCorr);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_softdrop_jec", &jetAK8_puppi_softdrop_jec, &b_jetAK8_puppi_softdrop_jec);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_tau1", &jetAK8_puppi_tau1, &b_jetAK8_puppi_tau1);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_tau2", &jetAK8_puppi_tau2, &b_jetAK8_puppi_tau2);
   inTrees.back()->SetBranchAddress("jetAK8_puppi_tau3", &jetAK8_puppi_tau3, &b_jetAK8_puppi_tau3);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_N", &jetAK8_subjet_puppi_softdrop_N, &b_jetAK8_subjet_puppi_softdrop_N);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_pt", &jetAK8_subjet_puppi_softdrop_pt, &b_jetAK8_subjet_puppi_softdrop_pt);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_eta", &jetAK8_subjet_puppi_softdrop_eta, &b_jetAK8_subjet_puppi_softdrop_eta);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_mass", &jetAK8_subjet_puppi_softdrop_mass, &b_jetAK8_subjet_puppi_softdrop_mass);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_phi", &jetAK8_subjet_puppi_softdrop_phi, &b_jetAK8_subjet_puppi_softdrop_phi);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_e", &jetAK8_subjet_puppi_softdrop_e, &b_jetAK8_subjet_puppi_softdrop_e);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_charge", &jetAK8_subjet_puppi_softdrop_charge, &b_jetAK8_subjet_puppi_softdrop_charge);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_genParton_pdgID", &jetAK8_subjet_puppi_softdrop_genParton_pdgID, &b_jetAK8_subjet_puppi_softdrop_genParton_pdgID);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_nbHadrons", &jetAK8_subjet_puppi_softdrop_nbHadrons, &b_jetAK8_subjet_puppi_softdrop_nbHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_ncHadrons", &jetAK8_subjet_puppi_softdrop_ncHadrons, &b_jetAK8_subjet_puppi_softdrop_ncHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_partonFlavour", &jetAK8_subjet_puppi_softdrop_partonFlavour, &b_jetAK8_subjet_puppi_softdrop_partonFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_hadronFlavour", &jetAK8_subjet_puppi_softdrop_hadronFlavour, &b_jetAK8_subjet_puppi_softdrop_hadronFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_puppi_softdrop_csv", &jetAK8_subjet_puppi_softdrop_csv, &b_jetAK8_subjet_puppi_softdrop_csv);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_N", &jetAK8_subjet_pruned_N, &b_jetAK8_subjet_pruned_N);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_pt", &jetAK8_subjet_pruned_pt, &b_jetAK8_subjet_pruned_pt);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_eta", &jetAK8_subjet_pruned_eta, &b_jetAK8_subjet_pruned_eta);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_mass", &jetAK8_subjet_pruned_mass, &b_jetAK8_subjet_pruned_mass);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_phi", &jetAK8_subjet_pruned_phi, &b_jetAK8_subjet_pruned_phi);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_e", &jetAK8_subjet_pruned_e, &b_jetAK8_subjet_pruned_e);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_charge", &jetAK8_subjet_pruned_charge, &b_jetAK8_subjet_pruned_charge);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_genParton_pdgID", &jetAK8_subjet_pruned_genParton_pdgID, &b_jetAK8_subjet_pruned_genParton_pdgID);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_nbHadrons", &jetAK8_subjet_pruned_nbHadrons, &b_jetAK8_subjet_pruned_nbHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_ncHadrons", &jetAK8_subjet_pruned_ncHadrons, &b_jetAK8_subjet_pruned_ncHadrons);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_partonFlavour", &jetAK8_subjet_pruned_partonFlavour, &b_jetAK8_subjet_pruned_partonFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_hadronFlavour", &jetAK8_subjet_pruned_hadronFlavour, &b_jetAK8_subjet_pruned_hadronFlavour);
   inTrees.back()->SetBranchAddress("jetAK8_subjet_pruned_csv", &jetAK8_subjet_pruned_csv, &b_jetAK8_subjet_pruned_csv);
   inTrees.back()->SetBranchAddress("HLT_isFired", &HLT_isFired, &b_HLT_isFired);
   inTrees.back()->SetBranchAddress("passFilter_HBHE", &passFilter_HBHE, &b_passFilter_HBHE_);
   inTrees.back()->SetBranchAddress("passFilter_HBHELoose", &passFilter_HBHELoose, &b_passFilter_HBHELoose_);
   inTrees.back()->SetBranchAddress("passFilter_HBHETight", &passFilter_HBHETight, &b_passFilter_HBHETight_);
   inTrees.back()->SetBranchAddress("passFilter_HBHEIso", &passFilter_HBHEIso, &b_passFilter_HBHEIso_);
   inTrees.back()->SetBranchAddress("passFilter_CSCHalo", &passFilter_CSCHalo, &b_passFilter_CSCHalo_);
   inTrees.back()->SetBranchAddress("passFilter_CSCTightHalo2015", &passFilter_CSCTightHalo2015, &b_passFilter_CSCTightHalo2015_);
   inTrees.back()->SetBranchAddress("passFilter_HCALlaser", &passFilter_HCALlaser, &b_passFilter_HCALlaser_);
   inTrees.back()->SetBranchAddress("passFilter_ECALDeadCell", &passFilter_ECALDeadCell, &b_passFilter_ECALDeadCell_);
   inTrees.back()->SetBranchAddress("passFilter_GoodVtx", &passFilter_GoodVtx, &b_passFilter_GoodVtx_);
   inTrees.back()->SetBranchAddress("passFilter_TrkFailure", &passFilter_TrkFailure, &b_passFilter_TrkFailure_);
   inTrees.back()->SetBranchAddress("passFilter_EEBadSc", &passFilter_EEBadSc, &b_passFilter_EEBadSc_);
   inTrees.back()->SetBranchAddress("passFilter_ECALlaser", &passFilter_ECALlaser, &b_passFilter_ECALlaser_);
   inTrees.back()->SetBranchAddress("passFilter_TrkPOG", &passFilter_TrkPOG, &b_passFilter_TrkPOG_);
   inTrees.back()->SetBranchAddress("passFilter_TrkPOG_manystrip", &passFilter_TrkPOG_manystrip, &b_passFilter_TrkPOG_manystrip_);
   inTrees.back()->SetBranchAddress("passFilter_TrkPOG_toomanystrip", &passFilter_TrkPOG_toomanystrip, &b_passFilter_TrkPOG_toomanystrip_);
   inTrees.back()->SetBranchAddress("passFilter_TrkPOG_logError", &passFilter_TrkPOG_logError, &b_passFilter_TrkPOG_logError_);
   inTrees.back()->SetBranchAddress("passFilter_METFilters", &passFilter_METFilters, &b_passFilter_METFilters_);
   inTrees.back()->SetBranchAddress("passFilter_CSCTightHaloTrkMuUnvetoFilter", &passFilter_CSCTightHaloTrkMuUnvetoFilter, &b_passFilter_CSCTightHaloTrkMuUnvetoFilter_);
   inTrees.back()->SetBranchAddress("passFilter_globalTightHalo2016", &passFilter_globalTightHalo2016, &b_passFilter_globalTightHalo2016_);
   inTrees.back()->SetBranchAddress("passFilter_HcalStripHalo", &passFilter_HcalStripHalo, &b_passFilter_HcalStripHalo_);
   inTrees.back()->SetBranchAddress("passFilter_chargedHadronTrackResolution", &passFilter_chargedHadronTrackResolution, &b_passFilter_chargedHadronTrackResolution_);
   inTrees.back()->SetBranchAddress("passFilter_muonBadTrack", &passFilter_muonBadTrack, &b_passFilter_muonBadTrack_);
   inTrees.back()->SetBranchAddress("EVENT_event", &EVENT_event, &b_EVENT_event);
   inTrees.back()->SetBranchAddress("EVENT_run", &EVENT_run, &b_EVENT_run);
   inTrees.back()->SetBranchAddress("EVENT_lumiBlock", &EVENT_lumiBlock, &b_EVENT_lumiBlock);
   inTrees.back()->SetBranchAddress("PV_N", &PV_N, &b_PV_N);
   inTrees.back()->SetBranchAddress("PV_filter", &PV_filter, &b_PV_filter);
   inTrees.back()->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   inTrees.back()->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   inTrees.back()->SetBranchAddress("PV_rho", &PV_rho, &b_PV_rho);
   inTrees.back()->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
  }
  std::random_shuffle(eventList.begin(), eventList.end());

  std::vector<TFile*> outFiles;
  std::vector<TTree*> outTrees;
  for (int i=0; i<20; ++i) {
    std::string name = "sig_shuffled_" + std::to_string(i) + ".root";
    outFiles.push_back(new TFile(name.c_str(), "RECREATE"));
    outTrees.push_back(inTrees.back()->CloneTree(0));
    cout << "outFile " << i << " is: " << outFiles.at(i)->GetName() << endl;
    cout << "outTree " << i << " is: " << outTrees.at(i) << endl;
  }
  string::size_type totalEntries = eventList.size();
  cout << " Number of entries total is " << totalEntries << endl;
  
  for ( string::size_type iEvent = 0; iEvent < totalEntries; ++iEvent) {
    
    //cout << "Will get event " << eventList.at(iEvent).second << " from file " << eventList.at(iEvent).first << endl;
    inTrees.at(eventList.at(iEvent).first)->GetEntry(eventList.at(iEvent).second);
    outTrees[iEvent%20]->Fill();
    if (iEvent%5000==0) {
      cout << std::fixed << std::setw(3) << std::setprecision(1) << (float(iEvent)/float(totalEntries))*100 << "% done: Scanned " << iEvent << " events." << endl;
    }
  }
   for (int i=0; i<20; ++i) {
     outFiles[i]->cd();
     outTrees[i]->Write();
     outFiles[i]->Close();
   }
}
