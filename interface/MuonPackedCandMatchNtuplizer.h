/**
 * \file MuonPackedCandMatchNtuplizer.h
 * \brief Declaration of the standalone MiniAOD muon-to-packed-candidate study analyzer.
 */

#ifndef HeavyFlavorAnalysis_TPSMuonPackedMatch_MuonPackedCandMatchNtuplizer_h
#define HeavyFlavorAnalysis_TPSMuonPackedMatch_MuonPackedCandMatchNtuplizer_h

#include <string>
#include <vector>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TTree;

/**
 * \class MuonPackedCandMatchNtuplizer
 * \brief Studies how MiniAOD muons can be associated with packed PF candidates.
 *
 * The analyzer writes one event-level tree for all studied muons and one retained
 * `(muon, packed candidate)` tree for disagreement cases so metric-based and
 * pointer-based matching views can be compared without modifying TPS-Onia2MuMu.
 */
class MuonPackedCandMatchNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  /// Track observables and propagated momentum uncertainties used by the match metrics.
  struct TrackKinematicsWithErrors {
    float pt = -9999.f;
    float eta = -9999.f;
    float phi = -9999.f;
    float px = -9999.f;
    float py = -9999.f;
    float pz = -9999.f;
    float ptErr = -9999.f;
    float etaErr = -9999.f;
    float phiErr = -9999.f;
    float dxyErr = -9999.f;
    float dzErr = -9999.f;
    float pxErr = -9999.f;
    float pyErr = -9999.f;
    float pzErr = -9999.f;
    bool hasTrackErrors = false;
  };

  /// Per-candidate comparison metrics stored in the retained-row tree.
  struct MatchMetrics {
    float deltaPRelVec = -9999.f;
    float chi2MomentumDiag = -9999.f;
    float deltaDzSelectedPV = -9999.f;
    float deltaDzAssocPV = -9999.f;
    float chi2MomentumDzSelectedPV = -9999.f;
    float chi2MomentumDzAssocPV = -9999.f;
    int vectorPass = 0;
    int chi2Pass = 0;
    int dzPvPass = 0;
    int dzAssocPass = 0;
  };

  explicit MuonPackedCandMatchNtuplizer(const edm::ParameterSet &iConfig);
  ~MuonPackedCandMatchNtuplizer() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) override;
  void endJob() override;

  /// Book the event-level and retained-row TTrees.
  void bookTrees();
  /// Reset per-event branch buffers before processing a new event.
  void clearEventBranches();
  /// Select the PV used for selected-PV observables and dz/dxy comparisons.
  int selectPrimaryVertexIndex(const reco::VertexCollection &vertices) const;
  /// Compute the scalar sum of squared track transverse momenta for a vertex.
  double vertexSumPt2(const reco::Vertex &vertex) const;
  /// Prefer muon.track() and fall back to muon.muonBestTrack() when available.
  TrackKinematicsWithErrors extractMuonTrackKinematics(const pat::Muon &muon, int &trackSource) const;
  /// Extract pseudoTrack-based observables for packed candidates with track details.
  TrackKinematicsWithErrors extractPackedTrackKinematics(const pat::PackedCandidate &cand) const;
  /// Copy reco::Track observables and uncertainties into the ntuple-friendly struct.
  TrackKinematicsWithErrors extractTrackKinematics(const reco::Track &track) const;
  /// Propagate the track covariance from (q/p, lambda, phi) to Cartesian momentum errors.
  void propagateCartesianMomentumErrors(const reco::Track &track, float &pxErr, float &pyErr, float &pzErr) const;

  // Consumes tokens and input tags configured from the Python cfg.
  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertexToken_;

  edm::InputTag muonTag_;
  edm::InputTag packedTag_;
  edm::InputTag primaryVertexTag_;

  // Runtime switches and thresholds controlling the matching study.
  bool studyAllMuons_;
  std::string muonSelectionStr_;
  StringCutObjectSelector<pat::Muon> muonSelector_;
  std::string pvSelectionMode_;
  double vectorRelPThreshold_;
  double momentumChi2Threshold_;
  double momentumDzPvChi2Threshold_;
  double momentumDzAssocChi2Threshold_;
  bool requireChargeMatch_;
  bool storeAllPrimaryVertices_;
  bool storePointerDiagnostics_;
  bool storeDetailedRowsOnlyOnDisagreement_;
  int loserRowsPerMethod_;
  std::string disagreementRowMode_;
  bool debugUnmatchedSoftMuons_;

  // TFileService products written by the analyzer.
  TTree *eventTree_;
  TTree *candidateTree_;

  // Event identifiers and summary counters.
  UInt_t run_;
  UInt_t lumi_;
  ULong64_t event_;
  Int_t selectedPVIndex_;
  Int_t nPV_;
  Int_t nPVStored_;
  Int_t nSlimmedMuons_;
  Int_t nMuStored_;

  // Primary-vertex payload stored in the event tree.
  std::vector<float> pv_x_;
  std::vector<float> pv_y_;
  std::vector<float> pv_z_;
  std::vector<float> pv_xErr_;
  std::vector<float> pv_yErr_;
  std::vector<float> pv_zErr_;
  std::vector<float> pv_chi2_;
  std::vector<float> pv_ndof_;
  std::vector<float> pv_sumPt2_;
  std::vector<int> pv_nTracks_;
  std::vector<int> pv_isValid_;

  // Muon-level event-tree payload: kinematics, IDs, pointer diagnostics, and match summaries.
  std::vector<int> mu_index_;
  std::vector<int> mu_passSelection_;
  std::vector<int> mu_charge_;
  std::vector<int> mu_hasTrack_;
  std::vector<int> mu_trackSource_;
  std::vector<float> mu_pt_;
  std::vector<float> mu_eta_;
  std::vector<float> mu_phi_;
  std::vector<float> mu_px_;
  std::vector<float> mu_py_;
  std::vector<float> mu_pz_;
  std::vector<float> mu_ptErr_;
  std::vector<float> mu_etaErr_;
  std::vector<float> mu_phiErr_;
  std::vector<float> mu_dxyErr_;
  std::vector<float> mu_dzErr_;
  std::vector<float> mu_pxErr_;
  std::vector<float> mu_pyErr_;
  std::vector<float> mu_pzErr_;
  std::vector<float> mu_dzSelectedPV_;
  std::vector<float> mu_dxySelectedPV_;
  std::vector<int> mu_passedCutBasedIdLoose_;
  std::vector<int> mu_passedCutBasedIdMedium_;
  std::vector<int> mu_passedCutBasedIdMediumPrompt_;
  std::vector<int> mu_passedCutBasedIdTight_;
  std::vector<int> mu_passedCutBasedIdGlobalHighPt_;
  std::vector<int> mu_passedCutBasedIdTrkHighPt_;
  std::vector<int> mu_passedCutBasedIdSoft_;
  std::vector<float> mu_mvaIDValue_;
  std::vector<int> mu_passedMvaIDwpSoft_;
  std::vector<int> mu_passedMvaIDwpMedium_;
  std::vector<int> mu_passedMvaIDwpTight_;
  std::vector<int> mu_nSourceCandidatePtrs_;
  std::vector<int> mu_hasPfCandidateRef_;
  std::vector<int> mu_pfCandidateRefResolved_;
  std::vector<int> mu_sourceCandidateResolvedCount_;
  std::vector<int> mu_pointerResolvedCount_;
  std::vector<int> mu_pointerMultiplicityAnomaly_;
  std::vector<int> mu_nMethodsConsidered_;
  std::vector<int> mu_nMethodsMatched_;
  std::vector<int> mu_nDistinctMatchedPackedIdx_;
  std::vector<int> mu_hasAnyMethodMatch_;
  std::vector<int> mu_isAgreementOnMatch_;
  std::vector<int> mu_isAgreementOnMismatch_;
  std::vector<int> mu_hasMethodDisagreement_;
  std::vector<int> mu_nPackedPassingAnyKinematic_;
  std::vector<int> mu_nPackedPointerResolved_;
  std::vector<int> mu_nPackedPassingAnyCriterion_;
  std::vector<float> mu_bestVectorRelP_;
  std::vector<float> mu_bestMomentumChi2_;
  std::vector<float> mu_bestMomentumDzPvChi2_;
  std::vector<float> mu_bestMomentumDzAssocChi2_;
  std::vector<int> mu_matchVectorPackedIdx_;
  std::vector<int> mu_matchChi2PackedIdx_;
  std::vector<int> mu_matchDzPvPackedIdx_;
  std::vector<int> mu_matchDzAssocPackedIdx_;
  std::vector<int> mu_matchPointerPackedIdx_;
  std::vector<int> mu_nPassVector_;
  std::vector<int> mu_nPassChi2_;
  std::vector<int> mu_nPassDzPv_;
  std::vector<int> mu_nPassDzAssoc_;

  // Candidate-tree payload for the current retained `(event, muon, packed candidate)` row.
  UInt_t cand_run_;
  UInt_t cand_lumi_;
  ULong64_t cand_event_;
  Int_t cand_selectedPVIndex_;
  Int_t cand_muonIndex_;
  Int_t cand_packedIndex_;
  Int_t cand_muPassSelection_;
  Int_t cand_muCharge_;
  Int_t cand_muPointerResolvedCount_;
  Float_t cand_muPt_;
  Float_t cand_muEta_;
  Float_t cand_muPhi_;
  Float_t cand_muPx_;
  Float_t cand_muPy_;
  Float_t cand_muPz_;
  Float_t cand_muPtErr_;
  Float_t cand_muPxErr_;
  Float_t cand_muPyErr_;
  Float_t cand_muPzErr_;
  Float_t cand_muDzSelectedPV_;
  Float_t cand_muDxySelectedPV_;
  Float_t cand_muDzErr_;

  Int_t cand_charge_;
  Int_t cand_pdgId_;
  Int_t cand_hasTrackDetails_;
  Int_t cand_trackHighPurity_;
  Int_t cand_vertexRefKey_;
  Int_t cand_pvAssocQuality_;
  Int_t cand_fromPV_;
  Float_t cand_pt_;
  Float_t cand_eta_;
  Float_t cand_phi_;
  Float_t cand_px_;
  Float_t cand_py_;
  Float_t cand_pz_;
  Float_t cand_mass_;
  Float_t cand_ptErr_;
  Float_t cand_etaErr_;
  Float_t cand_phiErr_;
  Float_t cand_dxyErr_;
  Float_t cand_dzErr_;
  Float_t cand_pxErr_;
  Float_t cand_pyErr_;
  Float_t cand_pzErr_;
  Float_t cand_dzAssociatedPV_;
  Float_t cand_dzSelectedPV_;
  Float_t cand_dxyAssociatedPV_;
  Float_t cand_dxySelectedPV_;
  Float_t cand_deltaPRelVec_;
  Float_t cand_chi2MomentumDiag_;
  Float_t cand_deltaDzSelectedPV_;
  Float_t cand_deltaDzAssocPV_;
  Float_t cand_chi2MomentumDzSelectedPV_;
  Float_t cand_chi2MomentumDzAssocPV_;
  Int_t cand_vectorPass_;
  Int_t cand_chi2Pass_;
  Int_t cand_dzPvPass_;
  Int_t cand_dzAssocPass_;
  Int_t cand_finalVector_;
  Int_t cand_finalChi2_;
  Int_t cand_finalDzPv_;
  Int_t cand_finalDzAssoc_;
  Int_t cand_pointerMatchAny_;
  Int_t cand_pointerMatchSource_;
  Int_t cand_pointerMatchPfRef_;
  Int_t cand_finalPointer_;
};

#endif
