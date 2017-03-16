#define sampleFactory_cxx
#include <iostream>
#include "sampleFactory.h"
#include <TROOT.h>
#include <TMath.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TH1F.h>



void sampleFactory::Loop(std::string mcName, std::string dataName, std::string sigName)
{
//   In a ROOT session, you can do:
//      root> .L sampleFactory.C
//      root> sampleFactory t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   /**
   * create output trees:
   * one is a straight copy of all branches specified in header
   * one is the fake data
   * one is the dummy signal
   */
   bool debug = true;
   TH1F* debugHist = new TH1F("debugHist", "dijet invariant mass", 200, 0, 1000);
   TFile* debugFile = new TFile("debug.root", "recreate");

   TFile* dataFile = new TFile(dataName.c_str(), "recreate");
   TTree* dataTree = fChain->CloneTree(0);
   TFile* sigFile = new TFile(sigName.c_str(), "recreate");
   TTree* sigTree = fChain->CloneTree(0);

   /**
   * also create some TLorentzVectors for doing dijet arithmetic
   */
   TLorentzVector jetA  = TLorentzVector();
   TLorentzVector jetB  = TLorentzVector();
   TLorentzVector dijet = TLorentzVector();

   /**
   * create a random number generator
   * we'll use it to try to create a "signal" looking distribution
   * where the dijet invariant mass is gaussian distributed around 750
   */
   TRandom* rand = new TRandom();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if (debug) {
        cout << "Generating a uniform random number for entry " << jentry << ": " << rand->Uniform(1) << endl;
      }
     
      if (jetAK4_N >= 2) {
        jetA.SetPtEtaPhiE(jetAK4_pt->at(0), jetAK4_eta->at(0), jetAK4_phi->at(0), jetAK4_e->at(0));
        jetB.SetPtEtaPhiE(jetAK4_pt->at(1), jetAK4_eta->at(1), jetAK4_phi->at(1), jetAK4_e->at(1));
        dijet = jetA + jetB;
        if (debug) {
          cout << "entry " << jentry << " has at invariant mass of two leading jets " << dijet.M() << endl;
          cout << "gaussian at " << dijet.M() << " is: " << TMath::Gaus(dijet.M(), 750, 750*0.05) << endl;
        }
        if (rand->Uniform(1) < TMath::Gaus(dijet.M(), 750, 750*0.05)) {
          sigTree->Fill();   
          if (debug) debugHist->Fill(dijet.M());
        }
        if (rand->Uniform(1) > 0.1) {
          dataTree->Fill();   
        }
      }
   }
   sigFile->cd();
   sigTree->Write();
   sigFile->Close();
   dataFile->cd();
   dataTree->Write();
   dataFile->Close();
   if (debug) {
     debugFile->cd();
     debugHist->Write();
     debugFile->Close();
   }
}
