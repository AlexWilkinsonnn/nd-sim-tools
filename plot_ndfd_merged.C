#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"

void plot_ndfd_merged()
{
  gStyle->SetOptStat(0);

  const std::string fname = "/dune/app/users/awilkins/nd-sim-tools/merge_nd-fd_test_more.root";

  TFile* f = TFile::Open(fname.c_str());
  TTree* t = (TTree*)f->Get("nd_fd_reco");

  std::cout << t->GetEntries() << " events in file. Tree looks like:\n";
  t->Show(0);

  TH1F* hFDNumuNuE = new TH1F("hFDNumuNuE", "FD ND Ev Reco", 50, 0, 10);
  TH1F* hEvReco = new TH1F("hEvReco", "FD ND Ev Reco", 50, 0, 10);
  t->Draw("Ev_reco>>+hEvReco", "", "goff");
  t->Draw("fd_numu_nu_E>>+hFDNumuNuE", "", "goff");
  TCanvas* c = new TCanvas();
  auto leg = new TLegend(0.50, 0.7, 0.85, 0.85);
  leg->SetTextSize(0.03);
  hEvReco->SetLineColor(kRed);
  leg->AddEntry(hEvReco, "ND EvReco");
  hEvReco->Draw("hist");
  hFDNumuNuE->SetLineColor(kBlue);
  leg->AddEntry(hFDNumuNuE, "FD EvReco");
  hFDNumuNuE->Draw("hist same");
  leg->Draw();

  // hFDNumuNuE = new TH1F("hFDNumuNuE", "FD ND Ev Reco muon_contained", 50, 0, 10);
  // hEvReco = new TH1F("hEvReco", "FD ND Ev Reco muon_contained", 50, 0, 10);
  // t->Draw("Ev_reco>>+hEvReco", "muon_contained==1", "goff");
  // t->Draw("fd_numu_nu_E>>+hFDNumuNuE", "muon_contained==1", "goff");
  // c = new TCanvas();
  // leg = new TLegend(0.50, 0.7, 0.85, 0.85);
  // leg->SetTextSize(0.03);
  // hEvReco->SetLineColor(kRed);
  // leg->AddEntry(hEvReco, "ND EvReco");
  // hEvReco->Draw("hist");
  // hFDNumuNuE->SetLineColor(kBlue);
  // leg->AddEntry(hFDNumuNuE, "FD EvReco");
  // hFDNumuNuE->Draw("hist same");
  // leg->Draw();

  TH1F* hEvDiff = new TH1F("hEvDiff", "ND - FD Ev Reco (muon contained)", 50, -3, 3);
  t->Draw("Ev_reco - fd_numu_nu_E>>+hEvDiff", "", "goff");
  c = new TCanvas();
  hEvDiff->Draw("hist");

  TH1F* hNumuCVNNDReco = new TH1F("hNumuCVNNDReco", "CVN Scores", 50, 0, 1);
  TH1F* hNumuCVNNDNoReco = new TH1F("hNumuCVNNDNoReco", "CVN Scores", 50, 0, 1);
  t->Draw("fd_numu_score>>+hNumuCVNNDReco", "reco_numu==1", "goff");
  t->Draw("fd_numu_score>>+hNumuCVNNDNoReco", "reco_numu==0", "goff");
  c = new TCanvas();
  leg = new TLegend(0.50, 0.7, 0.85, 0.85);
  leg->SetTextSize(0.03);
  hNumuCVNNDNoReco->SetLineColor(kBlue);
  leg->AddEntry(hNumuCVNNDNoReco, "reco_numu==0 at ND");
  hNumuCVNNDNoReco->Draw("hist");
  hNumuCVNNDReco->SetLineColor(kRed);
  leg->AddEntry(hNumuCVNNDReco, "reco_numu==1 at ND");
  hNumuCVNNDReco->Draw("hist same");
  leg->Draw();
}


