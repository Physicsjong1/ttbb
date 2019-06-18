/*
root -l -b -q ana.C'("step_1.root", "step_1_plots.root")'
*/

//------------------------------------------------------------------------------
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>

GenParticle* MotherParticle(const GenParticle* p, const TClonesArray* branchParticle, int i = 1){

  if( p == 0) return 0;

  int PID = p->PID;
  int index[2];
  index[i-1] = p->M1;
  index[i-1] = p->M2;

  GenParticle * m = 0;

  if( index[i] > 0){
    m = (GenParticle *) branchParticle->At( index[i] );
    double mPID = m->PID;
    cout << "particle " << PID << " (" << p->P4().Pt() << ") " << " is from " << " mother (M" << i << ") " << mPID << " ( " << m->P4().Pt() << " ) " << endl;
  }else{
    cout << "No mother particle " << "M" << i << " from  " << PID << "  !!!" << endl;
  }

  return m; 

}

GenParticle* DaughterParticle(const GenParticle* p, const TClonesArray* branchParticle, int i = 0){

  if( p == 0) return 0;

  int PID = p->PID;
  int index[2];
  index[i] = p->D1;
  index[i] = p->D2;

  GenParticle * m = 0;

  if( index[i] > 0){
    m = (GenParticle *) branchParticle->At( index[i] );
    double mPID = m->PID;
    cout << "particle " << PID << " (" << p->P4().Pt() << ") " << " decays to " << " particle (D" << i << ") " << mPID << " ( " << m->P4().Pt() << " ) " << endl;
  }else{
    cout << "No decaying particle " << "D" << i << " from  " << PID << "  !!!" << endl;
  }

  return m;

}

bool isFromTop(const GenParticle* p, const TClonesArray* branchParticle){

  bool output = false;
  bool debug = false;
  bool printout = false;
  double m1 = p->M1;
  double m2 = p->M2;
  int PID = p->PID;
  int mPID1 = -1;
  int mPID2 = -1;

  while( m1 >= 0 || m2 >= 0){
    if( m1 >= 0 && m2 < 0){
      GenParticle *mother = (GenParticle *) branchParticle->At(m1);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID1 = mother->PID;
      if(debug) cout << "case1 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) { 
        output = true;
        break;
      }
    }else if( m1 < 0 && m2 >= 0){
      GenParticle *mother = (GenParticle *) branchParticle->At(m2);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID2 = mother->PID;
      if(debug) cout << "case2 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) {
        output = true;
        break;
      }
    }else if( m1 == m2 ){
      GenParticle *mother = (GenParticle *) branchParticle->At(m1);
      m1 = mother->M1;
      m2 = mother->M2;
      mPID1 = mother->PID;
      mPID2 = mother->PID;
      if(debug) cout << "case 3 : mother = " << mother->PID << endl;
      if( abs(mother->PID) == 6 ) {
        output = true;
        break;
      }
    }else{
      GenParticle *mother1 = (GenParticle *) branchParticle->At(m1);
      GenParticle *mother2 = (GenParticle *) branchParticle->At(m2);

      mPID1 = mother1->PID;
      mPID2 = mother2->PID;
      if(debug) cout << "case4 : mother = " << mother1->PID << " " << mother2->PID << endl;

      double m11 = mother1->M1;
      double m12 = mother1->M2;
      double m21 = mother2->M1;
      double m22 = mother2->M2;
      double mo1_m1_PID = -1;
      double mo1_m2_PID = -1;
      double mo2_m1_PID = -1;
      double mo2_m2_PID = -1;
      if( m11 >= 0 ) {
        GenParticle * mo1_m1 = (GenParticle*) branchParticle->At(m11);
        mo1_m1_PID = mo1_m1->PID;
      }
      if( m12 >= 0 ) {
        GenParticle * mo1_m2 = (GenParticle*) branchParticle->At(m12);
        mo1_m2_PID =  mo1_m2->PID;
      }
      if( m21 >= 0 ) {
        GenParticle * mo2_m1 = (GenParticle*) branchParticle->At(m21);
        mo2_m1_PID = mo2_m1->PID;
      }
      if( m22 >= 0 ) {
        GenParticle * mo2_m2 = (GenParticle*) branchParticle->At(m22);
        mo2_m2_PID = mo2_m2->PID;
      }

      bool fromtop = abs(mother1->PID) == 6 || abs(mother2->PID) == 6; 
      bool fromtop2 = abs(mo1_m1_PID) == 6 || abs(mo1_m2_PID) == 6;
      bool fromtop3 = abs(mo2_m1_PID) == 6 || abs(mo2_m2_PID) == 6;
      if( fromtop ) {
        mPID1 = mother1->PID;
        mPID2 = mother2->PID;
        output = true;
        break;
      }
      if( fromtop2 ) {
        mPID1 = mo1_m1_PID;
        mPID2 = mo1_m2_PID;
        output = true;
        break;
      }
      if( fromtop3 ) {
        mPID1 = mo2_m1_PID;
        mPID2 = mo2_m2_PID;
        output = true;
        break;
      }

  if( (m11 >= 0 && m12 >= 0) && (m21 < 0 && m22 < 0) ) {
        m1 = m11;
        m2 = m12;
      }else if( (m11 < 0 && m12 < 0) && (m21 >= 0 && m22 >= 0) ){
        m1 = m21;
        m2 = m22;
      }else if( (m11 >= 0 && m12 < 0) && (m21 >= 0 && m22 < 0) ){
        m1 = m11;
        m2 = m21;
      }else if( (m11 < 0 && m12 >= 0) && (m21 < 0 && m22 >= 0) ){
        m1 = m12;
        m2 = m22;
      }else if( (m11 >= 0 && m12 < 0) && (m21 < 0 && m22 >= 0) ){
        m1 = m11;
        m2 = m22;
      }else if( (m11 < 0 && m12 >= 0) && (m21 >= 0 && m22 < 0) ){ 
        m1 = m12;
        m2 = m21; 
      }else if( (m11 < 0 && m12 < 0) && (m21 < 0 && m22 < 0)  ){
        break;
      }else{
        if( abs(mo1_m1_PID) == 5 || abs(mo1_m2_PID) == 5){
          m1 = m11;
          m2 = m12;
        }else if( abs(mo2_m1_PID) == 5 || abs(mo2_m2_PID) == 5)  {
          m1 = m21;
          m2 = m22;
        }else{
          if(debug){
            cout << "===== DEBUG =====" << endl;
            cout << m11 << " " << m12 << " " << m21 << " " << m22 << endl;
            cout << mo1_m1_PID << " " << mo1_m2_PID << " " << mo2_m1_PID << " " << mo2_m2_PID << endl;
            cout << "===== BREAK =====" << endl;
          }
          break;
        }
      }

    }
  }
  if(printout) cout << "Ancestor PID = " << mPID1 << " " << mPID2 << " test particle ID = " << PID << endl;

  return output;
}

void ana(const char *inputFile, const char *outputFile)
{
 gSystem->Load("libDelphes");

 // Create chain of root trees
 TChain chain("Delphes");
 chain.Add(inputFile);

 TString filename = inputFile;
 
 // Switch for single lepton or dilepton
 bool isdilepton = false;
 if( filename.Contains("di") == true ){
   isdilepton = true;
   cout<<"Dilepton Channel"<<endl;
 }
 else cout<<"Single Lepton Channel"<<endl;

 //Output
 TFile *fout = TFile::Open(outputFile,"RECREATE");
 fout->cd();

 double bjet1_pt;
 double bjet1_eta;
 double bjet1_phi;
 double bjet1_e;
 double bjet2_pt;
 double bjet2_eta;
 double bjet2_phi;
 double bjet2_e;
 double bjet3_pt;
 double bjet3_eta;
 double bjet3_phi;
 double bjet3_e;

 unsigned short nJet, nbJet, nMuon, nElectron;
 float Jet_pt, Jet_eta, Jet_phi;
 float Electron1_pt, Electron1_eta, Electron1_phi;
 float Electron2_pt, Electron2_eta, Electron2_phi;
 float Muon1_pt, Muon1_eta, Muon1_phi;
 float Muon2_pt, Muon2_eta, Muon2_phi;

 //Tree
 TTree * tree = new TTree( "tree", "tree for ttbb");
 tree->Branch("nbJet",&nbJet,"nbJet/s");
 tree->Branch("bjet1_pt",&bjet1_pt,"bjet1_pt/d");
 tree->Branch("bjet1_eta",&bjet1_eta,"bjet1_eta/d");
 tree->Branch("bjet1_phi",&bjet1_phi,"bjet1_phi/d");
 tree->Branch("bjet1_e",&bjet2_e,"bjet1_e/d");
 tree->Branch("bjet2_pt",&bjet2_pt,"bjet2_pt/d");
 tree->Branch("bjet2_eta",&bjet2_eta,"bjet2_eta/d");
 tree->Branch("bjet2_phi",&bjet2_phi,"bjet2_phi/d");
 tree->Branch("bjet2_e",&bjet2_e,"bjet2_e/d");
 tree->Branch("bjet3_pt",&bjet3_pt,"bjet3_pt/d");
 tree->Branch("bjet3_eta",&bjet3_eta,"bjet3_eta/d");
 tree->Branch("bjet3_phi",&bjet3_phi,"bjet3_phi/d");
 tree->Branch("bjet3_e",&bjet3_e,"bjet3_e/d");

 tree->Branch("nJet",&nJet,"nJet/s");
 tree->Branch("Jet_pt",&Jet_pt,"Jet_pt/f");
 tree->Branch("Jet_eta",&Jet_eta,"Jet_eta/f");
 tree->Branch("Jet_phi",&Jet_phi,"Jet_phi/f");

 tree->Branch("nElectron",&nElectron,"nElectron/s");
 tree->Branch("Electron1_pt",&Electron1_pt,"Electron1_pt/f");
 tree->Branch("Electron1_eta",&Electron1_eta,"Electron1_eta/f");
 tree->Branch("Electron1_phi",&Electron1_phi,"Electron1_phi/f");
 tree->Branch("Electron2_pt",&Electron2_pt,"Electron2_pt/f");
 tree->Branch("Electron2_eta",&Electron2_eta,"Electron2_eta/f");
 tree->Branch("Electron2_phi",&Electron2_phi,"Electron2_phi/f");
 
 tree->Branch("nMuon",&nMuon,"nMuon/s");
 tree->Branch("Muon1_pt",&Muon1_pt,"Muon1_pt/f");
 tree->Branch("Muon1_eta",&Muon1_eta,"Muon1_eta/f");
 tree->Branch("Muon1_phi",&Muon1_phi,"Muon1_phi/f"); 
 tree->Branch("Muon2_pt",&Muon2_pt,"Muon2_pt/f");
 tree->Branch("Muon2_eta",&Muon2_eta,"Muon2_eta/f");
 tree->Branch("Muon2_phi",&Muon2_phi,"Muon2_phi/f");

 // Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
 Long64_t numberOfEntries = treeReader->GetEntries();

 // Get pointers to branches used in this analysis
 TClonesArray *branchJet  = treeReader->UseBranch("Jet");
 TClonesArray *branchParticle  = treeReader->UseBranch("Particle");
 TClonesArray *branchElectron = treeReader->UseBranch("Electron");
 TClonesArray *branchMuon = treeReader->UseBranch("Muon");

 // Book histograms
 TH1 *histnjet = new TH1F("njet", "Number of jets", 14, 0.0, 14.0);
 TH1 *histnelectron = new TH1F("nelectron", "Number of electrons", 3, 0.0, 3.0);
 TH1 *histnmuon = new TH1F("nmuon", "Number of muons", 3, 0.0, 3.0);
 
 TH1 *histnbjet = new TH1F("nbjet", "Number of b-jets", 10, 0.0, 10.0);
 TH1 *histMbb = new TH1F("mbb", "M_{inv}(b, b)", 40, 20.0, 180.0);
 TH1 *histdRbb = new TH1F("dRbb", "dR(b, b)", 40, 0, 4.0);

 TH1 *hist_gennbjet = new TH1F("gennbjet", "Number of b-jets", 6, 0.0, 6.0);
 TH1 *hist_genMbb = new TH1F("genmbb", "M_{inv}(b, b)", 40, 20.0, 180.0);
 TH1 *hist_gendRbb = new TH1F("gendRbb", "dR(b, b)", 40, 0, 4.0);

 TH1 *hist_matchednbjet = new TH1F("matchednbjet", "Number of b-jets", 6, 0.0, 6.0);
 TH1 *hist_matchedMbb = new TH1F("matchedmbb", "M_{inv}(b, b)", 40, 20.0, 180.0);
 TH1 *hist_matcheddRbb = new TH1F("matcheddRbb", "dR(b, b)", 40, 0, 4.0);

 Int_t numberOfSelectedEvents = 0;
 Int_t numberOfMatchedEvents = 0;

 vector<Jet *> bJets;
 vector<Jet *> Jets;
 vector<Electron *> Electrons;
 vector<Muon *> Muons;

 TLorentzVector p4[2];
 Jet *jet;
 Electron *electron;
 Muon *muon;

 int entry, i, njet, nbjet, nelectron, nmuon;
 bool pass = false;

 // Loop over all events
 for(entry = 0; entry < numberOfEntries; ++entry)
 {
   if(entry%1000 == 0) cout << "event number: " << entry << endl;

   // Load selected branches with data from specified event
   treeReader->ReadEntry(entry);
   
   Jets.clear();
   Electrons.clear();
   Muons.clear(); 
   bJets.clear();
   
   Jet_pt = 999;
   Jet_eta = 999;
   Jet_phi = 999;

   Electron1_pt = 999;
   Electron1_eta = 999;
   Electron1_phi = 999;
   Electron2_pt = 999;
   Electron2_eta = 999;
   Electron2_phi = 999;

   Muon1_pt = 999;
   Muon1_eta = 999;
   Muon1_phi = 999;
   Muon2_pt = 999;
   Muon2_eta = 999;
   Muon2_phi = 999;

   bjet1_pt = 999;
   bjet1_eta = 999;
   bjet1_phi = 999;
   bjet1_e = 999;
   bjet2_pt = 999 ;
   bjet2_eta = 999;
   bjet2_phi = 999;
   bjet2_e = 999;
   bjet3_pt = 999 ;
   bjet3_eta = 999;
   bjet3_phi = 999;
   bjet3_e = 999;

   bool exitloop1 = false;
   bool exitloop2 = false;
   bool exitloop3 = false;

   // Jet and b-tag cuts
   njet = 0;
   nbjet = 0;
   for(i = 0; i < branchJet->GetEntriesFast(); ++i){
     jet = (Jet*) branchJet->At(i);
     if(jet->PT > 30 || fabs(jet->Eta) < 2.5){
       njet++;
       Jets.push_back(jet);
       if(jet->BTag) {
         nbjet++;
         bJets.push_back(jet);
       }
     }
     else exitloop1 = true;
   }
   //Electron cut
   nelectron = 0;
   for(i = 0; i < branchElectron->GetEntries(); ++i){
     electron = (Electron*) branchElectron->At(i);
     if( (electron->PT > 30) & (fabs(electron->Eta) <2.5)){
       nelectron++;
       Electrons.push_back(electron);
     }
     else exitloop2 = true;
   }
   //Muon cut
   nmuon = 0;
   for(i = 0; i < branchMuon->GetEntries(); ++i){
     muon = (Muon*) branchMuon->At(i);
     if( (muon->PT > 30) & (fabs(muon->Eta) < 2.5) ){
       nmuon++;
       Muons.push_back(muon);
     }
     else exitloop3 = true;
   }
   if (exitloop1 == true || exitloop2 == true || exitloop3 == true) continue;

   // Dilepton channel cuts
   if(isdilepton) {
     pass = ((nelectron >= 2) || (nmuon >= 2) || (nelectron >=1 & nmuon >= 1)) & (njet >= 2) & (nbjet >= 2);
     //cout<<"dicut"<<endl;
   }
   // Single lepton channel cuts
   else{
     pass = (nelectron == 1 || nmuon == 1) & (njet >=4) & (nbjet >= 2);
     //cout<<"singlecut"<<endl;
   }

   if(!pass) continue;

   // Fill the ntuples
   Jet_pt = Jets[0]->PT;
   Jet_eta = Jets[0]->Eta;
   Jet_phi = Jets[0]->Phi;
   if( nelectron == 1){
     Electron1_pt = Electrons[0]->PT;
     Electron1_eta = Electrons[0]->Eta;
     Electron1_phi = Electrons[0]->Phi;
   }
   if( nmuon ==1 ){
     Muon1_pt = Muons[0]->PT;
     Muon1_eta = Muons[0]->Eta;
     Muon1_phi = Muons[0]->Phi;
   }
   if(isdilepton){
     if(nelectron == 2){
       Electron2_pt = Electrons[1]->PT;
       Electron2_eta = Electrons[1]->Eta;
       Electron2_phi = Electrons[1]->Phi;
     }
     if(nmuon == 2){
       Muon2_pt = Muons[1]->PT;
       Muon2_eta = Muons[1]->Eta;
       Muon2_phi = Muons[1]->Phi;
     }
   }
   nJet = njet;
   nbJet = nbjet;
   nElectron = nelectron;
   nMuon = nmuon;

   histnjet->Fill( nJet );
   histnelectron->Fill( nelectron );
   histnmuon->Fill( nmuon );

   TLorentzVector addbjets[2]; 
   int nb = 0;
   int nbFromTop = 0;
   int nb_status3 = 0; 
   int ntop = 0;
   vector<GenParticle*> GenAddbJets;
   for(i = 0; i < branchParticle->GetEntriesFast(); ++i){
     GenParticle *genP = (GenParticle*) branchParticle->At(i);
     if( abs(genP->PID) == 6 ) {
       GenParticle *dauP1 = (GenParticle *) branchParticle->At(genP->D1);
       GenParticle *dauP2 = (GenParticle *) branchParticle->At(genP->D2);
       double dauPID1 = dauP1->PID;
       double dauPID2 = dauP2->PID;
       if( abs(dauPID1) != 6 && abs(dauPID2) != 6){
         ntop++;
       }
     }
     if( abs(genP->PID) == 5){
       if( genP->Status == 2) nb_status3++;
       //check if this is the last b quark 
       GenParticle *dauP1 = (GenParticle *) branchParticle->At(genP->D1);
       GenParticle *dauP2 = (GenParticle *) branchParticle->At(genP->D2);
       double dauPID1 = dauP1->PID;
       double dauPID2 = dauP2->PID;
       if( abs(dauPID1) == 5 || abs(dauPID2) == 5) continue; 
       //cout << "test b quark = " << genP->P4().Pt() << " decays to " << dauPID1 << " and " << dauPID2 << endl;
       nb++;
       bool fromTop = isFromTop(genP, branchParticle);
       if(fromTop) {
         nbFromTop++;
       }else{
         GenAddbJets.push_back(genP); 
       }
     }
   }

   //cout << "=========" << " Number of top = " << ntop << " number of b = " << nb << " (from top = " << nbFromTop << " )" << "=========" << endl;

   vector<Jet*> matchedbjets; 
   for( int i = 0; i < bJets.size(); i++){
     for(int j = 0 ; j < GenAddbJets.size(); j++){
       TLorentzVector recobjet = bJets[i]->P4();
       TLorentzVector genbjet = GenAddbJets[j]->P4();
       double dR = recobjet.DeltaR( genbjet );
       //cout << "test dR = " << dR << endl;
       if( dR < 0.4 ) matchedbjets.push_back( jet ) ;
     }
   }

   //cout << "matched = " << matchedbjets.size() << endl;
   histnbjet->Fill(bJets.size());
   hist_gennbjet->Fill(GenAddbJets.size());
   if( GenAddbJets.size() > 1){
     double genMbb = ( GenAddbJets[0]->P4() + GenAddbJets[1]->P4() ).M();
     double gendRbb = GenAddbJets[0]->P4().DeltaR( GenAddbJets[1]->P4() );
     hist_genMbb->Fill( genMbb );
     hist_gendRbb->Fill( gendRbb );
   }
   hist_matchednbjet->Fill( matchedbjets.size() );
   if( matchedbjets.size() > 1){
     double matched_mbb = (matchedbjets[0]->P4() + matchedbjets[1]->P4() ).M();
     double matched_dRbb = matchedbjets[0]->P4().DeltaR( matchedbjets[1]->P4() ); 
     hist_matchedMbb->Fill(matched_mbb);
     hist_matcheddRbb->Fill(matched_dRbb);
   }

   // select events with at least 2 b-jets and 2 opposite sign muons
   bjet1_pt = bJets[0]->P4().Pt();
   bjet1_eta = bJets[0]->P4().Eta();
   bjet1_phi = bJets[0]->P4().Phi();
   bjet1_e = bJets[0]->P4().E();
   bjet2_pt = bJets[1]->P4().Pt();
   bjet2_eta = bJets[1]->P4().Eta();
   bjet2_phi = bJets[1]->P4().Phi();
   bjet2_e = bJets[1]->P4().E();

   if(bJets.size() >=3){
     bjet3_pt = bJets[2]->P4().Pt();
     bjet3_eta = bJets[2]->P4().Eta();
     bjet3_phi = bJets[2]->P4().Phi();
     bjet3_e = bJets[2]->P4().E();
   }

   float mbb = 999;
   float dRbb = 999;

   // select two bjets with minimum dR
   TLorentzVector matchedRecoJets[2];
   for(int b1 = 0; b1 < bJets.size()-1; b1++){
     for(int b2 = b1+1; b2 < bJets.size(); b2++){
       p4[0] = bJets[b1]->P4();
       p4[1] = bJets[b2]->P4();

       float tmp_mbb = ((p4[0]) + (p4[1])).M();
       float tmp_dRbb = p4[0].DeltaR(p4[1]);
       if(tmp_dRbb < dRbb) {
          dRbb = tmp_dRbb;
          mbb = tmp_mbb; 
          matchedRecoJets[0] = p4[0];
          matchedRecoJets[1] = p4[1];
       }
     }
   }

   bool matched = false;
   bool matched1 = false;
   bool matched2 = false;
   for(int j = 0 ; j < GenAddbJets.size(); j++){
     if( matchedRecoJets[0].DeltaR( GenAddbJets[j]->P4() ) < 0.5 )  matched1 = true;
     if( matchedRecoJets[1].DeltaR( GenAddbJets[j]->P4() ) < 0.5 )  matched2 = true;
    }
   if( matched1 && matched2 ) matched = true;

   //cout<<"nelectron = "<< nelectron <<" nmuon = "<< nmuon <<" njet = "<< njet <<" nbjet = "<< nbjet << endl;

   ++numberOfSelectedEvents;
  if(matched) numberOfMatchedEvents++;

   histMbb->Fill(mbb);
   histdRbb->Fill(dRbb);
   tree->Fill();

 }

 cout << "Total number of events = " << numberOfSelectedEvents << endl;;
 cout << "Total number of matched events = " << numberOfMatchedEvents << endl;;
 double eff = (double) numberOfMatchedEvents/ (double) numberOfSelectedEvents;
 cout << "Matching eff. = " << eff << endl;
 fout->Write();
 fout->Close();

}



