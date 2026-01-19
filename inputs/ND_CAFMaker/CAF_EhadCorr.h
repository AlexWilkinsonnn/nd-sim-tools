#define CAF_EhadCorr_cxx
#ifdef CAF_EhadCorr_cxx

#include "TFile.h"
#include "TTree.h"
#include "Ntuple/NtpMCEventRecord.h"

class CAF_EhadCorr {

public:
  CAF_EhadCorr( std::string filename, bool isGas );
  ~CAF_EhadCorr();
  void fill();
  void fillPOT();
  void write();
  void addRWbranch( int parId, std::string name, std::string wgt_var, std::vector<double> &vars );
  void Print();
  void setToBS();

  // Make ntuple variables public so they can be set from other file
  
    // configuration variables
  int isFD, isFHC;
  // event accounting
  int run, subrun, event;
  // Truth information
  int isCC, neutrinoPDG, neutrinoPDGunosc, mode, LepPDG; 
  double Ev, Q2, W, X, Y, NuMomX, NuMomY, NuMomZ, LepMomX, LepMomY, LepMomZ, LepE, LepNuAngle;
  // True particle counts
  int nP, nN, nipip, nipim, nipi0, nikp, nikm, nik0, niem, niother, nNucleus, nUNKNOWN;
  double eP, eN, ePip, ePim, ePi0, eOther;
  double eRecoP, eRecoN, eRecoPip, eRecoPim, eRecoPi0, eRecoOther;
  double CorreRecoP, CorreRecoN, CorreRecoPip, CorreRecoPim, CorreRecoPi0, CorreRecoOther;
  double CorrP, CorrN, CorrPip, CorrPim, CorrPi0, CorrOther;
  // Particle energy info
  int nFSP;
  int pdg[100];
  float PrimEtrue[100];

  // vertex -- smear it?
  double vtx_x, vtx_y, vtx_z;
  double det_x;

  // Reco information CV
  double Ev_reco, CorrEv_reco, Elep_reco, Ehad_reco, CorrEhad_reco, \
	theta_reco;
  int reco_numu, reco_nue, reco_nc, reco_q;
  int muon_contained, muon_tracker, muon_ecal, muon_exit, reco_lepton_pdg;
  float muon_endpoint[3];
  std::string * muon_endVolName;
  double Ehad_veto, CorrEhad_veto;
  double TotCorr, TotCorrCollar;
  double pileup_energy;
  float PrimEreco[100];
  float PrimErecoCorr[100];
  float CorrPrimEreco[100];

  // Gas TPC variables
  int gastpc_pi_min_mult, gastpc_pi_pl_mult;
  double trkLen[100], trkLenPerp[100], ptrue[100], partEvReco[100];

  // reweights -- make sure big enough to hold all the variations for each knob, and all the knobs
  // the names, and what they actually mean, are determined automatically from the fhicl input file
  int nwgt[100];
  double cvwgt[100];
  double wgt[100][100];
  bool iswgt[100];

  // store the GENIE record as a branch
  genie::NtpMCEventRecord * mcrec;

  // Event-by-event geometric efficiency throw results
  std::vector< std::vector < std::vector < uint64_t > > > * geoEffThrowResults;

  // Infill info
  int vtxInGap;
  double hadEFracInGap;
  double lepEFracInGap;

  // meta
  double pot;
  int meta_run, meta_subrun;
  int version;

  TFile * cafFile;
  TTree * cafMVA;
  TTree * cafPOT;
  TTree * genie;
};

#endif  
