/*
Macro for Making tree from Delphes ROOT File.
Open with
root -l tree.C'("input_file_name.root","output_file_name.root")'
*/


#include <iostream>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif


void ntuple(const char* input, const char* output)
{
  gSystem->Load("libDelphes");

  TFile *out = new TFile(output, "recreate");
  TTree* tree = new TTree("tree","tree");
  
  const unsigned short Electron_N = 100;
  unsigned short nElectron;
  float Electron_pt[Electron_N], Electron_eta[Electron_N], Electron_phi[Electron_N];

  const unsigned short Muon_N = 100;
  unsigned short nMuon;
  float Muon_pt[Muon_N], Muon_eta[Muon_N], Muon_phi[Muon_N];
  
  const unsigned short Jet_N = 100;
  unsigned short nJet;
  float Jet_pt[Jet_N], Jet_eta[Jet_N], Jet_phi[Jet_N];

  tree->Branch("nElectron",&nElectron,"nElectron/s");
  tree->Branch("Electron_pt",Electron_pt,"Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta",Electron_eta,"Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi",Electron_phi,"Electron_phi[nElectron]/F");
  
  tree->Branch("nMuon",&nMuon,"nMuon/s");
  tree->Branch("Muon_pt",Muon_pt,"Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta",Muon_eta,"Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi",Muon_phi,"Muon_phi[nMuon]/F");
  
  tree->Branch("nJet",&nJet,"nJet/s");
  tree->Branch("Jet_pt",Jet_pt,"Jet_pt[nJet]/F");
  tree->Branch("Jet_eta",Jet_eta,"Jet_eta[nJet]/F");
  tree->Branch("Jet_phi",Jet_phi,"Jet_phi[nJet]/F");
  
  TChain chain("Delphes");
  chain.Add(input);

  ExRootTreeReader * treeReader = new ExRootTreeReader(&chain);
  long nentries = chain.GetEntries();

  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  
  for(int entry = 0; entry < nentries; ++entry){
 
    treeReader->ReadEntry(entry);

    nElectron=0;
    for(int i=0; i<branchElectron->GetEntries(); ++i){
      
      if (electron->Jet
      const Electron* electron = (const Electron*) branchElectron->At(i);

      Electron_pt[nElectron] = electron->PT;
      Electron_eta[nElectron] = electron->Eta;
      Electron_phi[nElectron] = electron->Phi;

      ++nElectron;
      if (nElectron >= Electron_N ) break;
    }
    
    nMuon=0;
    for(int i=0; i<branchMuon->GetEntries(); ++i){

      const Muon* muon = (const Muon*) branchMuon->At(i);

      Muon_pt[nMuon] = muon->PT;
      Muon_eta[nMuon] = muon->Eta;
      Muon_phi[nMuon] = muon->Phi;

      ++nMuon;
      if (nMuon >= Muon_N ) break;
    }
    
    nJet=0;
    for(int i=0; i<branchJet->GetEntries(); ++i){

      const Jet* jet = (const Jet*) branchJet->At(i);

      if (jet->BTag == 1) continue;
      
      Jet_pt[nJet] = jet->PT;
      Jet_eta[nJet] = jet->Eta;
      Jet_phi[nJet] = jet->Phi;

      ++nJet;
      if (nJet >= Jet_N ) break;
    }
    tree->Fill();
  }
  tree->Write();
  out->Close();
}
