// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TFitResultPtr.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TEfficiency.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include <TRandom3.h>
#include <TMarker3DBox.h>
#include "TSystem.h"
#include "TROOT.h"

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector> // Need this for generate dictionary for nested vectors

// Include customized functions and constants
#include "Helpers.h"

// Downstream fixed
void FDVtxPlot_z_down_fixed()
{
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 0. FD: read event from FD MC ntuple: before earth curvature rotation
  //
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // Ntuple path on FNAL dunegpvm machine
  // For FNAL machine:
  t->Add("/pnfs/dune/persistent/users/flynnguo/myFDntuples/myntuple_68092381_99all.root");

  // Define variables for FD event
  int FD_Run; // # of the run being processed
  int FD_SubRun; // # of the sub-run being processed
  int FD_Event; // # of the event being processed
  int FD_Sim_nNumu; // # of Sim muon neutrinos (numu and numubar)
  double FD_Gen_numu_E; // Energy of generator level neutrino [GeV]
  int FD_Sim_nMu; // # of Sim muons (mu+/mu-)
  int FD_CCNC_truth; // 0 =CC 1 =NC
  int FD_neuPDG; // Generator level neutrino PDG
  double FD_Sim_mu_start_vx; // Position of the muon trajectory at start point on the x-axis [cm]
  double FD_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double FD_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  double FD_Sim_mu_end_vx; // Position of the muon trajectory at end point on the x-axis [cm]
  double FD_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double FD_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double FD_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double FD_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double FD_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  double FD_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  double FD_Sim_mu_end_px; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double FD_Sim_mu_end_py; // Momentum of the muon trajectory at end point on the y-axis [GeV]
  double FD_Sim_mu_end_pz; // Momentum of the muon trajectory at end point on the z-axis [GeV]
  double FD_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  double FD_Sim_hadronic_Edep_b2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  int FD_Sim_n_hadronic_Edep_b; // # of hadronic energy deposits
  vector<float> *FD_Sim_hadronic_hit_Edep_b2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_b = 0; // Position of each energy deposit on the x-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_y_b = 0; // Position of each energy deposit on the y-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_z_b = 0; // Position of each energy deposit on the z-axis [cm]

  // Extract event info from ntuple
  t->SetBranchAddress("Run",                      &FD_Run);
  t->SetBranchAddress("SubRun",                   &FD_SubRun);
  t->SetBranchAddress("Event",                    &FD_Event);
  t->SetBranchAddress("Sim_nNumu",                &FD_Sim_nNumu);
  t->SetBranchAddress("Gen_numu_E",               &FD_Gen_numu_E);
  t->SetBranchAddress("Sim_nMu",                  &FD_Sim_nMu);
  t->SetBranchAddress("CCNC_truth",               &FD_CCNC_truth);
  t->SetBranchAddress("neuPDG",                   &FD_neuPDG);
  t->SetBranchAddress("Sim_mu_start_vx",          &FD_Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",          &FD_Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",          &FD_Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",            &FD_Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",            &FD_Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",            &FD_Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",          &FD_Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",          &FD_Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",          &FD_Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_start_E",           &FD_Sim_mu_start_E);
  t->SetBranchAddress("Sim_mu_end_px",            &FD_Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",            &FD_Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",            &FD_Sim_mu_end_pz);
  t->SetBranchAddress("Sim_mu_end_E",             &FD_Sim_mu_end_E);
  t->SetBranchAddress("Sim_hadronic_Edep_b2",     &FD_Sim_hadronic_Edep_b2);
  t->SetBranchAddress("Sim_n_hadronic_Edep_b",    &FD_Sim_n_hadronic_Edep_b);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_b2", &FD_Sim_hadronic_hit_Edep_b2);
  t->SetBranchAddress("Sim_hadronic_hit_x_b",     &FD_Sim_hadronic_hit_x_b);
  t->SetBranchAddress("Sim_hadronic_hit_y_b",     &FD_Sim_hadronic_hit_y_b);
  t->SetBranchAddress("Sim_hadronic_hit_z_b",     &FD_Sim_hadronic_hit_z_b);

  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Set Palette
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.85);		//Stat box x position (top right hand corner)
  gStyle->SetStatY(0.9); 		//Stat box y position
  gStyle->SetStatW(0.2);	 		//Stat box width as fraction of pad size
  gStyle->SetStatH(0.15);	 		//Size of each line in stat box

  // Set X,Y d_transverse
  vector<double> Dtransverse;
  Int_t Dtrans_steps = 0;
  Double_t Dtrans_stepsize = 10;

  Dtransverse.clear();

  Dtrans_steps = ( FDActiveVol_max[0] - 30 -10) / Dtrans_stepsize; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps + 1; i_Dtrans_step++ ){
    Dtransverse.emplace_back( i_Dtrans_step*Dtrans_stepsize );
    cout<< "Dtransverse_xy: " << i_Dtrans_step*Dtrans_stepsize << endl;
  }

  Int_t Dtransverse_size = Dtransverse.size();
  std::cout << "Dtransverse size: "<< Dtransverse_size << "\n\n";

  // Set Z d_transverse
  vector<double> Dtransverse_z;
  Int_t Dtrans_steps_z = 0;
  Double_t Dtrans_stepsize_z = 10;

  Dtransverse_z.clear();

  Dtrans_steps_z = 300 / Dtrans_stepsize_z; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps_z + 1; i_Dtrans_step++ ){
    Dtransverse_z.emplace_back( i_Dtrans_step*Dtrans_stepsize_z );
    cout<< "Dtransverse_z: " << i_Dtrans_step*Dtrans_stepsize_z << endl;
  }

  Int_t Dtransverse_z_size = Dtransverse_z.size();
  std::cout << "Dtransverse z size: "<< Dtransverse_z_size << "\n\n";

  Int_t Dtransverse_tot = Dtransverse_size*Dtransverse_z_size;
  std::cout << "Dtransverse tot size: "<< Dtransverse_tot << "\n\n";




  TFile * outFile = new TFile("Output_FDVtxPlots_Z_Downstream_Fixed.root", "RECREATE");



  //  If you don’t need ROOT’s object lifetime management, call TH1::AddDirectory(false) before creating the first histogram, and you’ll get “regular C++ behavior”.
  TH1::AddDirectory(kFALSE); // Avoid potential memory leak

  TGraph** g_dtransverse = new TGraph*[Dtransverse_z_size];


  Int_t z_counter = 0;
  for (Double_t i_Dtransverse_z : Dtransverse_z )
  {
    vector<double> x_list;
    vector<double> y_list;
    x_list.clear();
    y_list.clear();
    Int_t Dtransverse_counter = 0;



    for (Double_t i_Dtransverse : Dtransverse)
    {
      //
      // Declare variables used in this program
      //
      int nentries = 0; // Total input events
      float vetoEnergyFD; // Total hadron deposited energy in FD veto region
      int iwritten = 0; // Output event counter

      //
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //
      Double_t Eff_FD_HV = 0.; // Hardonic veto eff in FD
      Int_t Eff_FD_HV_nm = 0;
      Int_t Eff_FD_HV_dnm = 0;

      Double_t setFD_FV_min[] = { FDActiveVol_min[0]+30+i_Dtransverse, FDActiveVol_min[1]+30+i_Dtransverse, FDActiveVol_min[2]+30+i_Dtransverse_z};
      Double_t setFD_FV_max[] = { FDActiveVol_max[0]-30-i_Dtransverse, FDActiveVol_max[1]-30-i_Dtransverse, 960.};

      //------------------------------------------------------------------------------
      // Loop over FD events
      //
      nentries = t->GetEntries();
      std::cout << "Tot evts: " << nentries << std::endl;

      Int_t FV_iwritten = 0;
      for ( int ientry = 0; ientry < nentries; ientry++ )
      {
        t->GetEntry(ientry);
        if ( ientry%10000 == 0 )
        {
          std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;
        }

        //
        // Skip events without muon/hadronic deposits
        //
        if ( FD_Sim_nMu == 0 || FD_Sim_n_hadronic_Edep_b == 0 ) continue;
        if ( FD_CCNC_truth == 1) continue;   // only use CC events
        if ( abs(FD_neuPDG) != 14 ) continue;       // only use muon neu
        iwritten++;
        if(FD_Sim_mu_start_vx < setFD_FV_max[0] && FD_Sim_mu_start_vx > setFD_FV_min[0] && FD_Sim_mu_start_vy < setFD_FV_max[1] && FD_Sim_mu_start_vy > setFD_FV_min[1] && FD_Sim_mu_start_vz < setFD_FV_max[2] && FD_Sim_mu_start_vz > setFD_FV_min[2])
        {
          Eff_FD_HV_dnm++;

          vetoEnergyFD = 0.;
          // Loop over hadron E deposits
          for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_b; ihadronhit++ )
          {

            // Veto region size: 30 cm from the active volume
            if ( ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_min[0]      && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_min[1]      && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_min[2]      && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
                 ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_max[0] - 30 && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_max[0] ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_max[1] - 30 && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_max[1] ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_max[2] - 30 && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_max[2] )
               ){
                 vetoEnergyFD += FD_Sim_hadronic_hit_Edep_b2->at(ihadronhit);
            } // end if hadron deposit in FD veto region

          } // end loop over hadron E deposits

          // Skip FD event if the total hadron E in veto region exceeds vetoEnergy [MeV]
          //
          if (verbose) cout << "\n ientry: " << ientry << ", iwritten: " << iwritten << ", VetoEnergyFD[MeV]: " << vetoEnergyFD <<"\n\n";
          if ( vetoEnergyFD > 30 ) continue; // 30 MeV
          Eff_FD_HV_nm++;


        } // end if in FD FV

      } // end loop over events entries

      std::cout << "Tot entries: " << nentries << std::endl;
      std::cout << "Num of chosen events: " << iwritten << std::endl;
      std::cout << "Dtransverse: " << i_Dtransverse << ", Dtransverse_z: " << i_Dtransverse_z << endl;
      cout << "FD_FV_x: " << setFD_FV_min[0] << " - " << setFD_FV_max[0] << ", FD_FV_y: " << setFD_FV_min[1] << " - " << setFD_FV_max[1] << ", FD_FV_z: " << setFD_FV_min[2] << " - " << setFD_FV_max[2] << endl;

      // Calculate Eff_FD_HV
      if(Eff_FD_HV_dnm!=0)
      {
        Eff_FD_HV = Eff_FD_HV_nm*1.0/Eff_FD_HV_dnm;
      }
      else
      {
        Eff_FD_HV = 0;
      }


      std::cout << "# of evts in FV (Eveto<30MeV): " << Eff_FD_HV_nm << std::endl;
      std::cout << "# of evts in FV: " << Eff_FD_HV_dnm << std::endl;
      std::cout << "Eff_FD_HV: " << Eff_FD_HV << std::endl;
      if(Eff_FD_HV>1) cout<< "Eff_FD_HV OVER THAN 1, i_Dtransverse is: " << i_Dtransverse << endl;
      x_list.push_back(i_Dtransverse);
      y_list.push_back(Eff_FD_HV);

      Dtransverse_counter++;
      cout << "Dtransverse_counter: " << Dtransverse_counter << endl;

    } // end i_Dtransverse

    g_dtransverse[z_counter] = new TGraph(Dtransverse_size,&x_list[0],&y_list[0]);
    g_dtransverse[z_counter]->SetMinimum(0.);
    g_dtransverse[z_counter]->SetMaximum(1.0);
    g_dtransverse[z_counter]->SetMarkerStyle(20);
    g_dtransverse[z_counter]->SetMarkerColor(2);
    g_dtransverse[z_counter]->SetMarkerSize(2);
    g_dtransverse[z_counter]->Draw("AC*");
    TString d_transverse_z_name = Form("Eff vs Dtransverse, d_transverse_z=%.2f_cm, upstream_z=%.2f_cm", i_Dtransverse_z, i_Dtransverse_z+30);
    g_dtransverse[z_counter]->SetTitle(d_transverse_z_name);
    g_dtransverse[z_counter]->GetXaxis()->SetTitle("d_transverse [cm]");
    g_dtransverse[z_counter]->GetYaxis()->SetTitle("Eff_{FD}^{FV}");
    g_dtransverse[z_counter]->Write();
    z_counter++;
  } // end i_Dtransverse_z

  delete[] g_dtransverse;

  outFile->Close();
} // end main


// Upstream fixed
void FDVtxPlot_z_up_fixed()
{
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 0. FD: read event from FD MC ntuple: before earth curvature rotation
  //
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // Ntuple path on FNAL dunegpvm machine
  // For FNAL machine:
  t->Add("/pnfs/dune/persistent/users/flynnguo/myFDntuples/myntuple_68092381_99all.root");

  // Define variables for FD event
  int FD_Run; // # of the run being processed
  int FD_SubRun; // # of the sub-run being processed
  int FD_Event; // # of the event being processed
  int FD_Sim_nNumu; // # of Sim muon neutrinos (numu and numubar)
  double FD_Gen_numu_E; // Energy of generator level neutrino [GeV]
  int FD_Sim_nMu; // # of Sim muons (mu+/mu-)
  int FD_CCNC_truth; // 0 =CC 1 =NC
  int FD_neuPDG; // Generator level neutrino PDG
  double FD_Sim_mu_start_vx; // Position of the muon trajectory at start point on the x-axis [cm]
  double FD_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double FD_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  double FD_Sim_mu_end_vx; // Position of the muon trajectory at end point on the x-axis [cm]
  double FD_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double FD_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double FD_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double FD_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double FD_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  double FD_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  double FD_Sim_mu_end_px; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double FD_Sim_mu_end_py; // Momentum of the muon trajectory at end point on the y-axis [GeV]
  double FD_Sim_mu_end_pz; // Momentum of the muon trajectory at end point on the z-axis [GeV]
  double FD_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  double FD_Sim_hadronic_Edep_b2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  int FD_Sim_n_hadronic_Edep_b; // # of hadronic energy deposits
  vector<float> *FD_Sim_hadronic_hit_Edep_b2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_b = 0; // Position of each energy deposit on the x-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_y_b = 0; // Position of each energy deposit on the y-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_z_b = 0; // Position of each energy deposit on the z-axis [cm]

  // Extract event info from ntuple
  t->SetBranchAddress("Run",                      &FD_Run);
  t->SetBranchAddress("SubRun",                   &FD_SubRun);
  t->SetBranchAddress("Event",                    &FD_Event);
  t->SetBranchAddress("Sim_nNumu",                &FD_Sim_nNumu);
  t->SetBranchAddress("Gen_numu_E",               &FD_Gen_numu_E);
  t->SetBranchAddress("Sim_nMu",                  &FD_Sim_nMu);
  t->SetBranchAddress("CCNC_truth",               &FD_CCNC_truth);
  t->SetBranchAddress("neuPDG",                   &FD_neuPDG);
  t->SetBranchAddress("Sim_mu_start_vx",          &FD_Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",          &FD_Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",          &FD_Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",            &FD_Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",            &FD_Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",            &FD_Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",          &FD_Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",          &FD_Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",          &FD_Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_start_E",           &FD_Sim_mu_start_E);
  t->SetBranchAddress("Sim_mu_end_px",            &FD_Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",            &FD_Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",            &FD_Sim_mu_end_pz);
  t->SetBranchAddress("Sim_mu_end_E",             &FD_Sim_mu_end_E);
  t->SetBranchAddress("Sim_hadronic_Edep_b2",     &FD_Sim_hadronic_Edep_b2);
  t->SetBranchAddress("Sim_n_hadronic_Edep_b",    &FD_Sim_n_hadronic_Edep_b);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_b2", &FD_Sim_hadronic_hit_Edep_b2);
  t->SetBranchAddress("Sim_hadronic_hit_x_b",     &FD_Sim_hadronic_hit_x_b);
  t->SetBranchAddress("Sim_hadronic_hit_y_b",     &FD_Sim_hadronic_hit_y_b);
  t->SetBranchAddress("Sim_hadronic_hit_z_b",     &FD_Sim_hadronic_hit_z_b);

  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Set Palette
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.85);		//Stat box x position (top right hand corner)
  gStyle->SetStatY(0.9); 		//Stat box y position
  gStyle->SetStatW(0.2);	 		//Stat box width as fraction of pad size
  gStyle->SetStatH(0.15);	 		//Size of each line in stat box


  // Set X,Y d_transverse
  vector<double> Dtransverse;
  Int_t Dtrans_steps = 0;
  Double_t Dtrans_stepsize = 10;

  Dtransverse.clear();

  Dtrans_steps = ( FDActiveVol_max[0] - 30 -10) / Dtrans_stepsize; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps + 1; i_Dtrans_step++ ){
    Dtransverse.emplace_back( i_Dtrans_step*Dtrans_stepsize );
    cout<< "Dtransverse_xy: " << i_Dtrans_step*Dtrans_stepsize << endl;
  }

  Int_t Dtransverse_size = Dtransverse.size();
  std::cout << "Dtransverse size: "<< Dtransverse_size << "\n\n";

  // Set Z d_transverse
  vector<double> Dtransverse_z;
  Int_t Dtrans_steps_z = 0;
  Double_t Dtrans_stepsize_z = 10;

  Dtransverse_z.clear();

  Dtrans_steps_z = 600 / Dtrans_stepsize_z; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps_z + 1; i_Dtrans_step++ ){
    Dtransverse_z.emplace_back( i_Dtrans_step*Dtrans_stepsize_z );
    cout<< "Dtransverse_z: " << i_Dtrans_step*Dtrans_stepsize_z << endl;
  }

  Int_t Dtransverse_z_size = Dtransverse_z.size();
  std::cout << "Dtransverse z size: "<< Dtransverse_z_size << "\n\n";

  Int_t Dtransverse_tot = Dtransverse_size*Dtransverse_z_size;
  std::cout << "Dtransverse tot size: "<< Dtransverse_tot << "\n\n";



  TFile * outFile = new TFile("Output_FDVtxPlots_Z_Upstream_Fixed.root", "RECREATE");


  //  If you don’t need ROOT’s object lifetime management, call TH1::AddDirectory(false) before creating the first histogram, and you’ll get “regular C++ behavior”.
  TH1::AddDirectory(kFALSE); // Avoid potential memory leak

  TGraph** g_dtransverse = new TGraph*[Dtransverse_z_size];


  Int_t z_counter = 0;
  for (Double_t i_Dtransverse_z : Dtransverse_z )
  {
    vector<double> x_list;
    vector<double> y_list;
    x_list.clear();
    y_list.clear();
    Int_t Dtransverse_counter = 0;

    Double_t Eff_FD_HV = 0.; // Hardonic veto eff in FD
    Int_t Eff_FD_HV_nm = 0;
    Int_t Eff_FD_HV_dnm = 0;

    for (Double_t i_Dtransverse : Dtransverse)
    {
      //
      // Declare variables used in this program
      //
      int nentries = 0; // Total input events
      float vetoEnergyFD; // Total hadron deposited energy in FD veto region
      int iwritten = 0; // Output event counter
      //
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //
      Double_t Eff_FD_HV = 0.; // Hardonic veto eff in FD
      Int_t Eff_FD_HV_nm = 0;
      Int_t Eff_FD_HV_dnm = 0;

      Double_t setFD_FV_min[] = { FDActiveVol_min[0]+30+i_Dtransverse, FDActiveVol_min[1]+30+i_Dtransverse, 100.};
      Double_t setFD_FV_max[] = { FDActiveVol_max[0]-30-i_Dtransverse, FDActiveVol_max[1]-30-i_Dtransverse, FDActiveVol_max[2]-30-i_Dtransverse_z};

      //------------------------------------------------------------------------------
      // Loop over FD events
      //
      nentries = t->GetEntries();
      std::cout << "Tot evts: " << nentries << std::endl;

      Int_t FV_iwritten = 0;
      for ( int ientry = 0; ientry < nentries; ientry++ )
      {
        t->GetEntry(ientry);
        if ( ientry%10000 == 0 )
        {
          std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;
        }

        //
        // Skip events without muon/hadronic deposits
        //
        if ( FD_Sim_nMu == 0 || FD_Sim_n_hadronic_Edep_b == 0 ) continue;
        if ( FD_CCNC_truth == 1) continue;   // only use CC events
        if ( abs(FD_neuPDG) != 14 ) continue;       // only use muon neu
        iwritten++;
        if(FD_Sim_mu_start_vx < setFD_FV_max[0] && FD_Sim_mu_start_vx > setFD_FV_min[0] && FD_Sim_mu_start_vy < setFD_FV_max[1] && FD_Sim_mu_start_vy > setFD_FV_min[1] && FD_Sim_mu_start_vz < setFD_FV_max[2] && FD_Sim_mu_start_vz > setFD_FV_min[2])
        {
          Eff_FD_HV_dnm++;

          vetoEnergyFD = 0.;
          // Loop over hadron E deposits
          for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_b; ihadronhit++ )
          {

            // Veto region size: 30 cm from the active volume
            if ( ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_min[0]      && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_min[1]      && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_min[2]      && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
                 ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_max[0] - 30 && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_max[0] ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_max[1] - 30 && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_max[1] ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_max[2] - 30 && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_max[2] )
               ){
                 vetoEnergyFD += FD_Sim_hadronic_hit_Edep_b2->at(ihadronhit);
            } // end if hadron deposit in FD veto region

          } // end loop over hadron E deposits

          // Skip FD event if the total hadron E in veto region exceeds vetoEnergy [MeV]
          //
          if (verbose) cout << "\n ientry: " << ientry << ", iwritten: " << iwritten << ", VetoEnergyFD[MeV]: " << vetoEnergyFD <<"\n\n";
          if ( vetoEnergyFD > 30 ) continue; // 30 MeV
          Eff_FD_HV_nm++;


        } // end if in FD FV

      } // end loop over events entries

      std::cout << "Tot entries: " << nentries << std::endl;
      std::cout << "Num of chosen events: " << iwritten << std::endl;
      std::cout << "Dtransverse: " << i_Dtransverse << ", Dtransverse_z: " << i_Dtransverse_z << endl;
      cout << "FD_FV_x: " << setFD_FV_min[0] << " - " << setFD_FV_max[0] << ", FD_FV_y: " << setFD_FV_min[1] << " - " << setFD_FV_max[1] << ", FD_FV_z: " << setFD_FV_min[2] << " - " << setFD_FV_max[2] << endl;

      // Calculate Eff_FD_HV
      Eff_FD_HV = Eff_FD_HV_nm*1.0/Eff_FD_HV_dnm;
      std::cout << "# of evts in FV (Eveto<30MeV): " << Eff_FD_HV_nm << std::endl;
      std::cout << "# of evts in FV: " << Eff_FD_HV_dnm << std::endl;
      std::cout << "Eff_FD_HV: " << Eff_FD_HV << std::endl;
      if(Eff_FD_HV>1) cout<< "Eff_FD_HV OVER THAN 1, i_Dtransverse is: " << i_Dtransverse << endl;
      x_list.push_back(i_Dtransverse);
      y_list.push_back(Eff_FD_HV);

      Dtransverse_counter++;
      cout << "Dtransverse_counter: " << Dtransverse_counter << endl;

    } // end i_Dtransverse

    g_dtransverse[z_counter] = new TGraph(Dtransverse_size,&x_list[0],&y_list[0]);
    g_dtransverse[z_counter]->SetMinimum(0.);
    g_dtransverse[z_counter]->SetMaximum(1.0);
    g_dtransverse[z_counter]->SetMarkerStyle(20);
    g_dtransverse[z_counter]->SetMarkerColor(2);
    g_dtransverse[z_counter]->SetMarkerSize(2);
    g_dtransverse[z_counter]->Draw("AC*");
    TString d_transverse_z_name = Form("Eff vs Dtransverse, d_transverse_z=%.2f_cm, downstream=%.2f_cm", i_Dtransverse_z, FDActiveVol_max[2]-30-i_Dtransverse_z);
    g_dtransverse[z_counter]->SetTitle(d_transverse_z_name);
    g_dtransverse[z_counter]->GetXaxis()->SetTitle("d_transverse [cm]");
    g_dtransverse[z_counter]->GetYaxis()->SetTitle("Eff_{FD}^{FV}");
    g_dtransverse[z_counter]->Write();
    z_counter++;
  } // end i_Dtransverse_z

  delete[] g_dtransverse;

  outFile->Close();
} // end main

// Both sides (downstream+upstream) fixed
void FDVtxPlot_z_fixed_both()
{
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 0. FD: read event from FD MC ntuple: before earth curvature rotation
  //
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // Ntuple path on FNAL dunegpvm machine
  // For FNAL machine:
  t->Add("/pnfs/dune/persistent/users/flynnguo/myFDntuples/myntuple_68092381_99all.root");

  // Define variables for FD event
  int FD_Run; // # of the run being processed
  int FD_SubRun; // # of the sub-run being processed
  int FD_Event; // # of the event being processed
  int FD_Sim_nNumu; // # of Sim muon neutrinos (numu and numubar)
  double FD_Gen_numu_E; // Energy of generator level neutrino [GeV]
  int FD_Sim_nMu; // # of Sim muons (mu+/mu-)
  int FD_CCNC_truth; // 0 =CC 1 =NC
  int FD_neuPDG; // Generator level neutrino PDG
  double FD_Sim_mu_start_vx; // Position of the muon trajectory at start point on the x-axis [cm]
  double FD_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double FD_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  double FD_Sim_mu_end_vx; // Position of the muon trajectory at end point on the x-axis [cm]
  double FD_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double FD_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double FD_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double FD_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double FD_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  double FD_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  double FD_Sim_mu_end_px; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double FD_Sim_mu_end_py; // Momentum of the muon trajectory at end point on the y-axis [GeV]
  double FD_Sim_mu_end_pz; // Momentum of the muon trajectory at end point on the z-axis [GeV]
  double FD_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  double FD_Sim_hadronic_Edep_b2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  int FD_Sim_n_hadronic_Edep_b; // # of hadronic energy deposits
  vector<float> *FD_Sim_hadronic_hit_Edep_b2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_b = 0; // Position of each energy deposit on the x-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_y_b = 0; // Position of each energy deposit on the y-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_z_b = 0; // Position of each energy deposit on the z-axis [cm]

  // Extract event info from ntuple
  t->SetBranchAddress("Run",                      &FD_Run);
  t->SetBranchAddress("SubRun",                   &FD_SubRun);
  t->SetBranchAddress("Event",                    &FD_Event);
  t->SetBranchAddress("Sim_nNumu",                &FD_Sim_nNumu);
  t->SetBranchAddress("Gen_numu_E",               &FD_Gen_numu_E);
  t->SetBranchAddress("Sim_nMu",                  &FD_Sim_nMu);
  t->SetBranchAddress("CCNC_truth",               &FD_CCNC_truth);
  t->SetBranchAddress("neuPDG",                   &FD_neuPDG);
  t->SetBranchAddress("Sim_mu_start_vx",          &FD_Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",          &FD_Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",          &FD_Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",            &FD_Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",            &FD_Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",            &FD_Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",          &FD_Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",          &FD_Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",          &FD_Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_start_E",           &FD_Sim_mu_start_E);
  t->SetBranchAddress("Sim_mu_end_px",            &FD_Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",            &FD_Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",            &FD_Sim_mu_end_pz);
  t->SetBranchAddress("Sim_mu_end_E",             &FD_Sim_mu_end_E);
  t->SetBranchAddress("Sim_hadronic_Edep_b2",     &FD_Sim_hadronic_Edep_b2);
  t->SetBranchAddress("Sim_n_hadronic_Edep_b",    &FD_Sim_n_hadronic_Edep_b);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_b2", &FD_Sim_hadronic_hit_Edep_b2);
  t->SetBranchAddress("Sim_hadronic_hit_x_b",     &FD_Sim_hadronic_hit_x_b);
  t->SetBranchAddress("Sim_hadronic_hit_y_b",     &FD_Sim_hadronic_hit_y_b);
  t->SetBranchAddress("Sim_hadronic_hit_z_b",     &FD_Sim_hadronic_hit_z_b);

  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Set Palette
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.85);		//Stat box x position (top right hand corner)
  gStyle->SetStatY(0.9); 		//Stat box y position
  gStyle->SetStatW(0.2);	 		//Stat box width as fraction of pad size
  gStyle->SetStatH(0.15);	 		//Size of each line in stat box





  // Set X,Y d_transverse
  vector<double> Dtransverse;
  Int_t Dtrans_steps = 0;
  Double_t Dtrans_stepsize = 10;

  Dtransverse.clear();

  Dtrans_steps = ( FDActiveVol_max[0] - 30 -10) / Dtrans_stepsize; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps + 1; i_Dtrans_step++ ){
    Dtransverse.emplace_back( i_Dtrans_step*Dtrans_stepsize );
    cout<< "Dtransverse_xy: " << i_Dtrans_step*Dtrans_stepsize << endl;
  }

  Int_t Dtransverse_size = Dtransverse.size();
  std::cout << "Dtransverse size: "<< Dtransverse_size << "\n\n";

  // Set Z d_transverse
  vector<double> Dtransverse_z;
  Int_t Dtrans_steps_z = 0;
  Double_t Dtrans_stepsize_z = 10;

  Dtransverse_z.clear();

  Dtrans_steps_z = 300 / Dtrans_stepsize_z; // 30 is the veto size
  for ( int i_Dtrans_step = 0; i_Dtrans_step < Dtrans_steps_z + 1; i_Dtrans_step++ ){
    Dtransverse_z.emplace_back( i_Dtrans_step*Dtrans_stepsize_z );
    cout<< "Dtransverse_z: " << i_Dtrans_step*Dtrans_stepsize_z << endl;
  }

  Int_t Dtransverse_z_size = Dtransverse_z.size();
  std::cout << "Dtransverse z size: "<< Dtransverse_z_size << "\n\n";

  Int_t Dtransverse_tot = Dtransverse_size*Dtransverse_z_size;
  std::cout << "Dtransverse tot size: "<< Dtransverse_tot << "\n\n";




  TFile * outFile = new TFile("Output_FDVtxPlots_Z_ivy.root", "RECREATE");



  //  If you don’t need ROOT’s object lifetime management, call TH1::AddDirectory(false) before creating the first histogram, and you’ll get “regular C++ behavior”.
  TH1::AddDirectory(kFALSE); // Avoid potential memory leak

  TGraph** g_dtransverse = new TGraph*[Dtransverse_z_size];


    Int_t z_counter = 0;
    vector<double> x_list;
    vector<double> y_list;
    x_list.clear();
    y_list.clear();
    Int_t Dtransverse_counter = 0;



    for (Double_t i_Dtransverse : Dtransverse)
    {
      //
      // Declare variables used in this program
      //
      int nentries = 0; // Total input events
      float vetoEnergyFD; // Total hadron deposited energy in FD veto region
      int iwritten = 0; // Output event counter

      //
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //
      Double_t Eff_FD_HV = 0.; // Hardonic veto eff in FD
      Int_t Eff_FD_HV_nm = 0;
      Int_t Eff_FD_HV_dnm = 0;

      Double_t setFD_FV_min[] = { FDActiveVol_min[0]+30+i_Dtransverse, FDActiveVol_min[1]+30+i_Dtransverse, 100.};
      Double_t setFD_FV_max[] = { FDActiveVol_max[0]-30-i_Dtransverse, FDActiveVol_max[1]-30-i_Dtransverse, 1244.};

      //------------------------------------------------------------------------------
      // Loop over FD events
      //
      nentries = t->GetEntries();
      std::cout << "Tot evts: " << nentries << std::endl;

      Int_t FV_iwritten = 0;
      for ( int ientry = 0; ientry < nentries; ientry++ )
      {
        t->GetEntry(ientry);
        if ( ientry%10000 == 0 )
        {
          std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;
        }

        //
        // Skip events without muon/hadronic deposits
        //
        if ( FD_Sim_nMu == 0 || FD_Sim_n_hadronic_Edep_b == 0 ) continue;
        if ( FD_CCNC_truth == 1) continue;   // only use CC events
        if ( abs(FD_neuPDG) != 14 ) continue;       // only use muon neu
        iwritten++;

        if(FD_Sim_mu_start_vx < setFD_FV_max[0] && FD_Sim_mu_start_vx > setFD_FV_min[0] && FD_Sim_mu_start_vy < setFD_FV_max[1] && FD_Sim_mu_start_vy > setFD_FV_min[1] && FD_Sim_mu_start_vz < setFD_FV_max[2] && FD_Sim_mu_start_vz > setFD_FV_min[2])
        {
          Eff_FD_HV_dnm++;

          vetoEnergyFD = 0.;
          // Loop over hadron E deposits
          for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_b; ihadronhit++ )
          {

            // Veto region size: 30 cm from the active volume
            if ( ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_min[0]      && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_min[1]      && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_min[2]      && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
                 ( FD_Sim_hadronic_hit_x_b->at(ihadronhit) > FDActiveVol_max[0] - 30 && FD_Sim_hadronic_hit_x_b->at(ihadronhit) < FDActiveVol_max[0] ) ||
                 ( FD_Sim_hadronic_hit_y_b->at(ihadronhit) > FDActiveVol_max[1] - 30 && FD_Sim_hadronic_hit_y_b->at(ihadronhit) < FDActiveVol_max[1] ) ||
                 ( FD_Sim_hadronic_hit_z_b->at(ihadronhit) > FDActiveVol_max[2] - 30 && FD_Sim_hadronic_hit_z_b->at(ihadronhit) < FDActiveVol_max[2] )
               ){
                 vetoEnergyFD += FD_Sim_hadronic_hit_Edep_b2->at(ihadronhit);
            } // end if hadron deposit in FD veto region

          } // end loop over hadron E deposits

          // Skip FD event if the total hadron E in veto region exceeds vetoEnergy [MeV]
          //
          if (verbose) cout << "\n ientry: " << ientry << ", iwritten: " << iwritten << ", VetoEnergyFD[MeV]: " << vetoEnergyFD <<"\n\n";
          if ( vetoEnergyFD > 30 ) continue; // 30 MeV
          Eff_FD_HV_nm++;


        } // end if in FD FV

      } // end loop over events entries

      std::cout << "Tot entries: " << nentries << std::endl;
      std::cout << "Num of chosen events: " << iwritten << std::endl;
      cout << "FD_FV_x: " << setFD_FV_min[0] << " - " << setFD_FV_max[0] << ", FD_FV_y: " << setFD_FV_min[1] << " - " << setFD_FV_max[1] << ", FD_FV_z: " << setFD_FV_min[2] << " - " << setFD_FV_max[2] << endl;

      // Calculate Eff_FD_HV
      Eff_FD_HV = Eff_FD_HV_nm*1.0/Eff_FD_HV_dnm;
      std::cout << "# of evts in FV (Eveto<30MeV): " << Eff_FD_HV_nm << std::endl;
      std::cout << "# of evts in FV: " << Eff_FD_HV_dnm << std::endl;
      std::cout << "Eff_FD_HV: " << Eff_FD_HV << std::endl;
      if(Eff_FD_HV>1) cout<< "Eff_FD_HV OVER THAN 1, i_Dtransverse is: " << i_Dtransverse << endl;
      x_list.push_back(i_Dtransverse);
      y_list.push_back(Eff_FD_HV);

      Dtransverse_counter++;
      cout << "Dtransverse_counter: " << Dtransverse_counter << endl;

    } // end i_Dtransverse

    g_dtransverse[z_counter] = new TGraph(Dtransverse_size,&x_list[0],&y_list[0]);
    g_dtransverse[z_counter]->SetMinimum(0.);
    g_dtransverse[z_counter]->SetMaximum(1.0);
    g_dtransverse[z_counter]->SetMarkerStyle(20);
    g_dtransverse[z_counter]->SetMarkerColor(2);
    g_dtransverse[z_counter]->SetMarkerSize(2);
    g_dtransverse[z_counter]->Draw("AC*");
    g_dtransverse[z_counter]->GetXaxis()->SetTitle("d_transverse [cm]");
    g_dtransverse[z_counter]->GetYaxis()->SetTitle("Eff_{FD}^{FV}");
    g_dtransverse[z_counter]->Write();
    z_counter++;

  delete[] g_dtransverse;

  outFile->Close();
} // end main
