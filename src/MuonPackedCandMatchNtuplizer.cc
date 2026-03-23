#include "../interface/MuonPackedCandMatchNtuplizer.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"

namespace {

constexpr float kInvalidFloat = -9999.f;
constexpr int kInvalidInt = -1;

inline double sqr(double value) { return value * value; }

struct PointerSummary {
  int nSourceCandidatePtrs = 0;
  int hasPfCandidateRef = 0;
  int pfCandidateRefResolved = 0;
  int sourceCandidateResolvedCount = 0;
  int pointerResolvedCount = 0;
  int pointerMultiplicityAnomaly = 0;
  int firstPointerPackedIdx = kInvalidInt;
  std::set<int> allMatches;
  std::set<int> sourceMatches;
  std::set<int> pfMatches;
};

struct LocalCandidateStudy {
  int packedIndex = kInvalidInt;
  int pointerMatchAny = 0;
  int pointerMatchSource = 0;
  int pointerMatchPfRef = 0;
  MuonPackedCandMatchNtuplizer::TrackKinematicsWithErrors trackKin;
  MuonPackedCandMatchNtuplizer::MatchMetrics metrics;
};

struct DebugCandidateRank {
  int packedIndex = kInvalidInt;
  double rawMomentumChi2 = std::numeric_limits<double>::infinity();
  float deltaPRelVec = kInvalidFloat;
  int chargeMatch = 0;
  int pointerMatchSource = 0;
  int pointerMatchPfRef = 0;
  MuonPackedCandMatchNtuplizer::TrackKinematicsWithErrors trackKin;
};

std::string formatFloatForDebug(double value) {
  if (!std::isfinite(value) || value == static_cast<double>(kInvalidFloat)) {
    return "invalid";
  }

  std::ostringstream stream;
  stream << std::fixed << std::setprecision(6) << value;
  return stream.str();
}

std::string formatProductIdForDebug(const edm::ProductID &productId, bool valid) {
  if (!valid) {
    return "invalid";
  }

  std::ostringstream stream;
  stream << productId;
  return stream.str();
}

}  // namespace

MuonPackedCandMatchNtuplizer::MuonPackedCandMatchNtuplizer(const edm::ParameterSet &iConfig)
    : muonToken_(consumes<edm::View<pat::Muon>>(
          iConfig.getUntrackedParameter<edm::InputTag>("muons", edm::InputTag("slimmedMuons")))),
      packedToken_(consumes<edm::View<pat::PackedCandidate>>(
          iConfig.getUntrackedParameter<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates")))),
      primaryVertexToken_(consumes<reco::VertexCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices",
                                                       edm::InputTag("offlineSlimmedPrimaryVertices")))),
      muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muons", edm::InputTag("slimmedMuons"))),
      packedTag_(iConfig.getUntrackedParameter<edm::InputTag>("packedCandidates",
                                                              edm::InputTag("packedPFCandidates"))),
      primaryVertexTag_(iConfig.getUntrackedParameter<edm::InputTag>(
          "primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"))),
      studyAllMuons_(iConfig.getUntrackedParameter<bool>("studyAllMuons", true)),
      muonSelectionStr_(iConfig.getUntrackedParameter<std::string>("MuonSelection",
                                                                   "pt > 2.5 && abs(eta) < 2.4")),
      muonSelector_(muonSelectionStr_),
      pvSelectionMode_(iConfig.getUntrackedParameter<std::string>("PVSelectionMode", "firstVertex")),
      legacyBoxThreshold_(iConfig.getUntrackedParameter<double>("LegacyBoxThreshold", 0.5)),
      vectorRelPThreshold_(iConfig.getUntrackedParameter<double>("VectorRelPThreshold", 0.01)),
      momentumChi2Threshold_(iConfig.getUntrackedParameter<double>("MomentumChi2Threshold", 25.0)),
      momentumDzPvChi2Threshold_(iConfig.getUntrackedParameter<double>("MomentumDzPvChi2Threshold", 25.0)),
      momentumDzAssocChi2Threshold_(
          iConfig.getUntrackedParameter<double>("MomentumDzAssocChi2Threshold", 25.0)),
      requireChargeMatch_(iConfig.getUntrackedParameter<bool>("RequireChargeMatch", true)),
      storeAllPrimaryVertices_(iConfig.getUntrackedParameter<bool>("StoreAllPrimaryVertices", true)),
      storePointerDiagnostics_(iConfig.getUntrackedParameter<bool>("StorePointerDiagnostics", true)),
      debugUnmatchedSoftMuons_(iConfig.getUntrackedParameter<bool>("DebugUnmatchedSoftMuons", false)),
      eventTree_(nullptr),
      candidateTree_(nullptr),
      run_(0),
      lumi_(0),
      event_(0),
      selectedPVIndex_(kInvalidInt),
      nPV_(0),
      nSlimmedMuons_(0),
      nMuStored_(0),
      cand_run_(0),
      cand_lumi_(0),
      cand_event_(0),
      cand_selectedPVIndex_(kInvalidInt),
      cand_muonIndex_(kInvalidInt),
      cand_packedIndex_(kInvalidInt),
      cand_muPassSelection_(0),
      cand_muCharge_(0),
      cand_muPointerResolvedCount_(0),
      cand_muPt_(kInvalidFloat),
      cand_muEta_(kInvalidFloat),
      cand_muPhi_(kInvalidFloat),
      cand_muPx_(kInvalidFloat),
      cand_muPy_(kInvalidFloat),
      cand_muPz_(kInvalidFloat),
      cand_muPtErr_(kInvalidFloat),
      cand_muPxErr_(kInvalidFloat),
      cand_muPyErr_(kInvalidFloat),
      cand_muPzErr_(kInvalidFloat),
      cand_muDzSelectedPV_(kInvalidFloat),
      cand_muDxySelectedPV_(kInvalidFloat),
      cand_muDzErr_(kInvalidFloat),
      cand_charge_(0),
      cand_pdgId_(0),
      cand_hasTrackDetails_(0),
      cand_trackHighPurity_(0),
      cand_vertexRefKey_(kInvalidInt),
      cand_pvAssocQuality_(kInvalidInt),
      cand_fromPV_(kInvalidInt),
      cand_pt_(kInvalidFloat),
      cand_eta_(kInvalidFloat),
      cand_phi_(kInvalidFloat),
      cand_px_(kInvalidFloat),
      cand_py_(kInvalidFloat),
      cand_pz_(kInvalidFloat),
      cand_mass_(kInvalidFloat),
      cand_ptErr_(kInvalidFloat),
      cand_etaErr_(kInvalidFloat),
      cand_phiErr_(kInvalidFloat),
      cand_dxyErr_(kInvalidFloat),
      cand_dzErr_(kInvalidFloat),
      cand_pxErr_(kInvalidFloat),
      cand_pyErr_(kInvalidFloat),
      cand_pzErr_(kInvalidFloat),
      cand_dzAssociatedPV_(kInvalidFloat),
      cand_dzSelectedPV_(kInvalidFloat),
      cand_dxyAssociatedPV_(kInvalidFloat),
      cand_dxySelectedPV_(kInvalidFloat),
      cand_deltaPRelVec_(kInvalidFloat),
      cand_chi2MomentumDiag_(kInvalidFloat),
      cand_deltaDzSelectedPV_(kInvalidFloat),
      cand_deltaDzAssocPV_(kInvalidFloat),
      cand_chi2MomentumDzSelectedPV_(kInvalidFloat),
      cand_chi2MomentumDzAssocPV_(kInvalidFloat),
      cand_legacyBoxPass_(0),
      cand_vectorPass_(0),
      cand_chi2Pass_(0),
      cand_dzPvPass_(0),
      cand_dzAssocPass_(0),
      cand_finalLegacy_(0),
      cand_finalVector_(0),
      cand_finalChi2_(0),
      cand_finalDzPv_(0),
      cand_finalDzAssoc_(0),
      cand_pointerMatchAny_(0),
      cand_pointerMatchSource_(0),
      cand_pointerMatchPfRef_(0),
      cand_finalPointer_(0) {
  usesResource("TFileService");
}

void MuonPackedCandMatchNtuplizer::beginJob() { bookTrees(); }

void MuonPackedCandMatchNtuplizer::endJob() {}

void MuonPackedCandMatchNtuplizer::bookTrees() {
  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Events", "Muon to packed-candidate matching study");
  candidateTree_ = fs->make<TTree>("MuonCandMatch", "Retained muon-candidate comparison rows");

  eventTree_->Branch("run", &run_, "run/i");
  eventTree_->Branch("lumi", &lumi_, "lumi/i");
  eventTree_->Branch("event", &event_, "event/l");
  eventTree_->Branch("selectedPVIndex", &selectedPVIndex_, "selectedPVIndex/I");
  eventTree_->Branch("nPV", &nPV_, "nPV/I");
  eventTree_->Branch("nSlimmedMuons", &nSlimmedMuons_, "nSlimmedMuons/I");
  eventTree_->Branch("nMuStored", &nMuStored_, "nMuStored/I");
  eventTree_->Branch("pv_x", &pv_x_);
  eventTree_->Branch("pv_y", &pv_y_);
  eventTree_->Branch("pv_z", &pv_z_);
  eventTree_->Branch("pv_xErr", &pv_xErr_);
  eventTree_->Branch("pv_yErr", &pv_yErr_);
  eventTree_->Branch("pv_zErr", &pv_zErr_);
  eventTree_->Branch("pv_chi2", &pv_chi2_);
  eventTree_->Branch("pv_ndof", &pv_ndof_);
  eventTree_->Branch("pv_sumPt2", &pv_sumPt2_);
  eventTree_->Branch("pv_nTracks", &pv_nTracks_);
  eventTree_->Branch("pv_isValid", &pv_isValid_);
  eventTree_->Branch("mu_index", &mu_index_);
  eventTree_->Branch("mu_passSelection", &mu_passSelection_);
  eventTree_->Branch("mu_charge", &mu_charge_);
  eventTree_->Branch("mu_hasTrack", &mu_hasTrack_);
  eventTree_->Branch("mu_trackSource", &mu_trackSource_);
  eventTree_->Branch("mu_pt", &mu_pt_);
  eventTree_->Branch("mu_eta", &mu_eta_);
  eventTree_->Branch("mu_phi", &mu_phi_);
  eventTree_->Branch("mu_px", &mu_px_);
  eventTree_->Branch("mu_py", &mu_py_);
  eventTree_->Branch("mu_pz", &mu_pz_);
  eventTree_->Branch("mu_ptErr", &mu_ptErr_);
  eventTree_->Branch("mu_etaErr", &mu_etaErr_);
  eventTree_->Branch("mu_phiErr", &mu_phiErr_);
  eventTree_->Branch("mu_dxyErr", &mu_dxyErr_);
  eventTree_->Branch("mu_dzErr", &mu_dzErr_);
  eventTree_->Branch("mu_pxErr", &mu_pxErr_);
  eventTree_->Branch("mu_pyErr", &mu_pyErr_);
  eventTree_->Branch("mu_pzErr", &mu_pzErr_);
  eventTree_->Branch("mu_dzSelectedPV", &mu_dzSelectedPV_);
  eventTree_->Branch("mu_dxySelectedPV", &mu_dxySelectedPV_);
  eventTree_->Branch("mu_passedCutBasedIdLoose", &mu_passedCutBasedIdLoose_);
  eventTree_->Branch("mu_passedCutBasedIdMedium", &mu_passedCutBasedIdMedium_);
  eventTree_->Branch("mu_passedCutBasedIdMediumPrompt", &mu_passedCutBasedIdMediumPrompt_);
  eventTree_->Branch("mu_passedCutBasedIdTight", &mu_passedCutBasedIdTight_);
  eventTree_->Branch("mu_passedCutBasedIdGlobalHighPt", &mu_passedCutBasedIdGlobalHighPt_);
  eventTree_->Branch("mu_passedCutBasedIdTrkHighPt", &mu_passedCutBasedIdTrkHighPt_);
  eventTree_->Branch("mu_passedCutBasedIdSoft", &mu_passedCutBasedIdSoft_);
  eventTree_->Branch("mu_passedMvaIDwpSoft", &mu_passedMvaIDwpSoft_);
  eventTree_->Branch("mu_mvaIDValue", &mu_mvaIDValue_);
  eventTree_->Branch("mu_passedMvaIDwpMedium", &mu_passedMvaIDwpMedium_);
  eventTree_->Branch("mu_passedMvaIDwpTight", &mu_passedMvaIDwpTight_);
  eventTree_->Branch("mu_nSourceCandidatePtrs", &mu_nSourceCandidatePtrs_);
  eventTree_->Branch("mu_hasPfCandidateRef", &mu_hasPfCandidateRef_);
  eventTree_->Branch("mu_pfCandidateRefResolved", &mu_pfCandidateRefResolved_);
  eventTree_->Branch("mu_sourceCandidateResolvedCount", &mu_sourceCandidateResolvedCount_);
  eventTree_->Branch("mu_pointerResolvedCount", &mu_pointerResolvedCount_);
  eventTree_->Branch("mu_pointerMultiplicityAnomaly", &mu_pointerMultiplicityAnomaly_);
  eventTree_->Branch("mu_matchLegacyPackedIdx", &mu_matchLegacyPackedIdx_);
  eventTree_->Branch("mu_matchVectorPackedIdx", &mu_matchVectorPackedIdx_);
  eventTree_->Branch("mu_matchChi2PackedIdx", &mu_matchChi2PackedIdx_);
  eventTree_->Branch("mu_matchDzPvPackedIdx", &mu_matchDzPvPackedIdx_);
  eventTree_->Branch("mu_matchDzAssocPackedIdx", &mu_matchDzAssocPackedIdx_);
  eventTree_->Branch("mu_matchPointerPackedIdx", &mu_matchPointerPackedIdx_);
  eventTree_->Branch("mu_nPassLegacyBox", &mu_nPassLegacyBox_);
  eventTree_->Branch("mu_nPassVector", &mu_nPassVector_);
  eventTree_->Branch("mu_nPassChi2", &mu_nPassChi2_);
  eventTree_->Branch("mu_nPassDzPv", &mu_nPassDzPv_);
  eventTree_->Branch("mu_nPassDzAssoc", &mu_nPassDzAssoc_);

  candidateTree_->Branch("run", &cand_run_, "run/i");
  candidateTree_->Branch("lumi", &cand_lumi_, "lumi/i");
  candidateTree_->Branch("event", &cand_event_, "event/l");
  candidateTree_->Branch("selectedPVIndex", &cand_selectedPVIndex_, "selectedPVIndex/I");
  candidateTree_->Branch("muonIndex", &cand_muonIndex_, "muonIndex/I");
  candidateTree_->Branch("packedIndex", &cand_packedIndex_, "packedIndex/I");
  candidateTree_->Branch("mu_passSelection", &cand_muPassSelection_, "mu_passSelection/I");
  candidateTree_->Branch("mu_charge", &cand_muCharge_, "mu_charge/I");
  candidateTree_->Branch("mu_pointerResolvedCount", &cand_muPointerResolvedCount_, "mu_pointerResolvedCount/I");
  candidateTree_->Branch("mu_pt", &cand_muPt_, "mu_pt/F");
  candidateTree_->Branch("mu_eta", &cand_muEta_, "mu_eta/F");
  candidateTree_->Branch("mu_phi", &cand_muPhi_, "mu_phi/F");
  candidateTree_->Branch("mu_px", &cand_muPx_, "mu_px/F");
  candidateTree_->Branch("mu_py", &cand_muPy_, "mu_py/F");
  candidateTree_->Branch("mu_pz", &cand_muPz_, "mu_pz/F");
  candidateTree_->Branch("mu_ptErr", &cand_muPtErr_, "mu_ptErr/F");
  candidateTree_->Branch("mu_pxErr", &cand_muPxErr_, "mu_pxErr/F");
  candidateTree_->Branch("mu_pyErr", &cand_muPyErr_, "mu_pyErr/F");
  candidateTree_->Branch("mu_pzErr", &cand_muPzErr_, "mu_pzErr/F");
  candidateTree_->Branch("mu_dzSelectedPV", &cand_muDzSelectedPV_, "mu_dzSelectedPV/F");
  candidateTree_->Branch("mu_dxySelectedPV", &cand_muDxySelectedPV_, "mu_dxySelectedPV/F");
  candidateTree_->Branch("mu_dzErr", &cand_muDzErr_, "mu_dzErr/F");
  candidateTree_->Branch("cand_charge", &cand_charge_, "cand_charge/I");
  candidateTree_->Branch("cand_pdgId", &cand_pdgId_, "cand_pdgId/I");
  candidateTree_->Branch("cand_hasTrackDetails", &cand_hasTrackDetails_, "cand_hasTrackDetails/I");
  candidateTree_->Branch("cand_trackHighPurity", &cand_trackHighPurity_, "cand_trackHighPurity/I");
  candidateTree_->Branch("cand_vertexRefKey", &cand_vertexRefKey_, "cand_vertexRefKey/I");
  candidateTree_->Branch("cand_pvAssocQuality", &cand_pvAssocQuality_, "cand_pvAssocQuality/I");
  candidateTree_->Branch("cand_fromPV", &cand_fromPV_, "cand_fromPV/I");
  candidateTree_->Branch("cand_pt", &cand_pt_, "cand_pt/F");
  candidateTree_->Branch("cand_eta", &cand_eta_, "cand_eta/F");
  candidateTree_->Branch("cand_phi", &cand_phi_, "cand_phi/F");
  candidateTree_->Branch("cand_px", &cand_px_, "cand_px/F");
  candidateTree_->Branch("cand_py", &cand_py_, "cand_py/F");
  candidateTree_->Branch("cand_pz", &cand_pz_, "cand_pz/F");
  candidateTree_->Branch("cand_mass", &cand_mass_, "cand_mass/F");
  candidateTree_->Branch("cand_ptErr", &cand_ptErr_, "cand_ptErr/F");
  candidateTree_->Branch("cand_etaErr", &cand_etaErr_, "cand_etaErr/F");
  candidateTree_->Branch("cand_phiErr", &cand_phiErr_, "cand_phiErr/F");
  candidateTree_->Branch("cand_dxyErr", &cand_dxyErr_, "cand_dxyErr/F");
  candidateTree_->Branch("cand_dzErr", &cand_dzErr_, "cand_dzErr/F");
  candidateTree_->Branch("cand_pxErr", &cand_pxErr_, "cand_pxErr/F");
  candidateTree_->Branch("cand_pyErr", &cand_pyErr_, "cand_pyErr/F");
  candidateTree_->Branch("cand_pzErr", &cand_pzErr_, "cand_pzErr/F");
  candidateTree_->Branch("cand_dzAssociatedPV", &cand_dzAssociatedPV_, "cand_dzAssociatedPV/F");
  candidateTree_->Branch("cand_dzSelectedPV", &cand_dzSelectedPV_, "cand_dzSelectedPV/F");
  candidateTree_->Branch("cand_dxyAssociatedPV", &cand_dxyAssociatedPV_, "cand_dxyAssociatedPV/F");
  candidateTree_->Branch("cand_dxySelectedPV", &cand_dxySelectedPV_, "cand_dxySelectedPV/F");
  candidateTree_->Branch("deltaPRelVec", &cand_deltaPRelVec_, "deltaPRelVec/F");
  candidateTree_->Branch("chi2MomentumDiag", &cand_chi2MomentumDiag_, "chi2MomentumDiag/F");
  candidateTree_->Branch("deltaDzSelectedPV", &cand_deltaDzSelectedPV_, "deltaDzSelectedPV/F");
  candidateTree_->Branch("deltaDzAssocPV", &cand_deltaDzAssocPV_, "deltaDzAssocPV/F");
  candidateTree_->Branch("chi2MomentumDzSelectedPV", &cand_chi2MomentumDzSelectedPV_,
                         "chi2MomentumDzSelectedPV/F");
  candidateTree_->Branch("chi2MomentumDzAssocPV", &cand_chi2MomentumDzAssocPV_,
                         "chi2MomentumDzAssocPV/F");
  candidateTree_->Branch("legacyBoxPass", &cand_legacyBoxPass_, "legacyBoxPass/I");
  candidateTree_->Branch("vectorPass", &cand_vectorPass_, "vectorPass/I");
  candidateTree_->Branch("chi2Pass", &cand_chi2Pass_, "chi2Pass/I");
  candidateTree_->Branch("dzPvPass", &cand_dzPvPass_, "dzPvPass/I");
  candidateTree_->Branch("dzAssocPass", &cand_dzAssocPass_, "dzAssocPass/I");
  candidateTree_->Branch("finalLegacy", &cand_finalLegacy_, "finalLegacy/I");
  candidateTree_->Branch("finalVector", &cand_finalVector_, "finalVector/I");
  candidateTree_->Branch("finalChi2", &cand_finalChi2_, "finalChi2/I");
  candidateTree_->Branch("finalDzPv", &cand_finalDzPv_, "finalDzPv/I");
  candidateTree_->Branch("finalDzAssoc", &cand_finalDzAssoc_, "finalDzAssoc/I");
  candidateTree_->Branch("pointerMatchAny", &cand_pointerMatchAny_, "pointerMatchAny/I");
  candidateTree_->Branch("pointerMatchSource", &cand_pointerMatchSource_, "pointerMatchSource/I");
  candidateTree_->Branch("pointerMatchPfRef", &cand_pointerMatchPfRef_, "pointerMatchPfRef/I");
  candidateTree_->Branch("finalPointer", &cand_finalPointer_, "finalPointer/I");
}

void MuonPackedCandMatchNtuplizer::clearEventBranches() {
  selectedPVIndex_ = kInvalidInt;
  nPV_ = 0;
  nSlimmedMuons_ = 0;
  nMuStored_ = 0;

  pv_x_.clear();
  pv_y_.clear();
  pv_z_.clear();
  pv_xErr_.clear();
  pv_yErr_.clear();
  pv_zErr_.clear();
  pv_chi2_.clear();
  pv_ndof_.clear();
  pv_sumPt2_.clear();
  pv_nTracks_.clear();
  pv_isValid_.clear();

  mu_index_.clear();
  mu_passSelection_.clear();
  mu_charge_.clear();
  mu_hasTrack_.clear();
  mu_trackSource_.clear();
  mu_pt_.clear();
  mu_eta_.clear();
  mu_phi_.clear();
  mu_px_.clear();
  mu_py_.clear();
  mu_pz_.clear();
  mu_ptErr_.clear();
  mu_etaErr_.clear();
  mu_phiErr_.clear();
  mu_dxyErr_.clear();
  mu_dzErr_.clear();
  mu_pxErr_.clear();
  mu_pyErr_.clear();
  mu_pzErr_.clear();
  mu_dzSelectedPV_.clear();
  mu_dxySelectedPV_.clear();
  mu_passedCutBasedIdLoose_.clear();
  mu_passedCutBasedIdMedium_.clear();
  mu_passedCutBasedIdMediumPrompt_.clear();
  mu_passedCutBasedIdTight_.clear();
  mu_passedCutBasedIdGlobalHighPt_.clear();
  mu_passedCutBasedIdTrkHighPt_.clear();
  mu_passedCutBasedIdSoft_.clear();
  mu_passedMvaIDwpSoft_.clear();
  mu_mvaIDValue_.clear();
  mu_passedMvaIDwpMedium_.clear();
  mu_passedMvaIDwpTight_.clear();
  mu_nSourceCandidatePtrs_.clear();
  mu_hasPfCandidateRef_.clear();
  mu_pfCandidateRefResolved_.clear();
  mu_sourceCandidateResolvedCount_.clear();
  mu_pointerResolvedCount_.clear();
  mu_pointerMultiplicityAnomaly_.clear();
  mu_matchLegacyPackedIdx_.clear();
  mu_matchVectorPackedIdx_.clear();
  mu_matchChi2PackedIdx_.clear();
  mu_matchDzPvPackedIdx_.clear();
  mu_matchDzAssocPackedIdx_.clear();
  mu_matchPointerPackedIdx_.clear();
  mu_nPassLegacyBox_.clear();
  mu_nPassVector_.clear();
  mu_nPassChi2_.clear();
  mu_nPassDzPv_.clear();
  mu_nPassDzAssoc_.clear();
}

int MuonPackedCandMatchNtuplizer::selectPrimaryVertexIndex(const reco::VertexCollection &vertices) const {
  if (vertices.empty()) {
    return kInvalidInt;
  }

  if (pvSelectionMode_ == "firstVertex") {
    return 0;
  }

  int bestIndex = 0;
  double bestValue = -std::numeric_limits<double>::infinity();
  for (size_t idx = 0; idx < vertices.size(); ++idx) {
    const auto &vertex = vertices[idx];
    const double value = (pvSelectionMode_ == "mostTracks") ? static_cast<double>(vertex.tracksSize())
                                                            : vertexSumPt2(vertex);
    if (value > bestValue) {
      bestValue = value;
      bestIndex = static_cast<int>(idx);
    }
  }
  return bestIndex;
}

double MuonPackedCandMatchNtuplizer::vertexSumPt2(const reco::Vertex &vertex) const {
  double sumPt2 = 0.0;
  for (auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
    if (it->isNonnull()) {
      sumPt2 += (*it)->pt() * (*it)->pt();
    }
  }
  return sumPt2;
}

MuonPackedCandMatchNtuplizer::TrackKinematicsWithErrors
MuonPackedCandMatchNtuplizer::extractTrackKinematics(const reco::Track &track) const {
  TrackKinematicsWithErrors result;
  result.pt = track.pt();
  result.eta = track.eta();
  result.phi = track.phi();
  result.px = track.px();
  result.py = track.py();
  result.pz = track.pz();
  result.ptErr = track.ptError();
  result.etaErr = track.etaError();
  result.phiErr = track.phiError();
  result.dxyErr = track.dxyError();
  result.dzErr = track.dzError();

  propagateCartesianMomentumErrors(track, result.pxErr, result.pyErr, result.pzErr);
  result.hasTrackErrors = (result.pxErr > 0.f && result.pyErr > 0.f && result.pzErr > 0.f &&
                           std::isfinite(result.pxErr) && std::isfinite(result.pyErr) &&
                           std::isfinite(result.pzErr));
  return result;
}

void MuonPackedCandMatchNtuplizer::propagateCartesianMomentumErrors(const reco::Track &track,
                                                                    float &pxErr,
                                                                    float &pyErr,
                                                                    float &pzErr) const {
  pxErr = kInvalidFloat;
  pyErr = kInvalidFloat;
  pzErr = kInvalidFloat;

  const double qop = track.qoverp();
  if (!std::isfinite(qop) || std::abs(qop) < 1e-12) {
    return;
  }

  const double lambda = track.lambda();
  const double phi = track.phi();
  const double p = std::abs(1.0 / qop);
  const double signQop = (qop > 0.0) ? 1.0 : -1.0;
  const double dpdqop = -signQop / (qop * qop);

  const double cosLam = std::cos(lambda);
  const double sinLam = std::sin(lambda);
  const double cosPhi = std::cos(phi);
  const double sinPhi = std::sin(phi);

  const double covQQ = track.covariance(reco::TrackBase::i_qoverp, reco::TrackBase::i_qoverp);
  const double covQL = track.covariance(reco::TrackBase::i_qoverp, reco::TrackBase::i_lambda);
  const double covQF = track.covariance(reco::TrackBase::i_qoverp, reco::TrackBase::i_phi);
  const double covLL = track.covariance(reco::TrackBase::i_lambda, reco::TrackBase::i_lambda);
  const double covLF = track.covariance(reco::TrackBase::i_lambda, reco::TrackBase::i_phi);
  const double covFF = track.covariance(reco::TrackBase::i_phi, reco::TrackBase::i_phi);

  const double dpxdq = dpdqop * cosLam * cosPhi;
  const double dpxdl = -p * sinLam * cosPhi;
  const double dpxdf = -p * cosLam * sinPhi;

  const double dpydq = dpdqop * cosLam * sinPhi;
  const double dpydl = -p * sinLam * sinPhi;
  const double dpydf = p * cosLam * cosPhi;

  const double dpzdq = dpdqop * sinLam;
  const double dpzdl = p * cosLam;
  const double dpzdf = 0.0;

  const auto variance = [&](double dq, double dl, double df) {
    return dq * dq * covQQ + dl * dl * covLL + df * df * covFF +
           2.0 * (dq * dl * covQL + dq * df * covQF + dl * df * covLF);
  };

  const double pxVar = variance(dpxdq, dpxdl, dpxdf);
  const double pyVar = variance(dpydq, dpydl, dpydf);
  const double pzVar = variance(dpzdq, dpzdl, dpzdf);

  if (pxVar >= 0.0) {
    pxErr = std::sqrt(pxVar);
  }
  if (pyVar >= 0.0) {
    pyErr = std::sqrt(pyVar);
  }
  if (pzVar >= 0.0) {
    pzErr = std::sqrt(pzVar);
  }
}

MuonPackedCandMatchNtuplizer::TrackKinematicsWithErrors
MuonPackedCandMatchNtuplizer::extractMuonTrackKinematics(const pat::Muon &muon, int &trackSource) const {
  trackSource = 0;
  if (!muon.track().isNull() && muon.track().isAvailable()) {
    trackSource = 1;
    return extractTrackKinematics(*muon.track());
  }

  try {
    const reco::TrackRef bestTrack = muon.muonBestTrack();
    if (bestTrack.isNonnull() && bestTrack.isAvailable()) {
      trackSource = 2;
      return extractTrackKinematics(*bestTrack);
    }
  } catch (...) {
  }

  return TrackKinematicsWithErrors();
}

MuonPackedCandMatchNtuplizer::TrackKinematicsWithErrors
MuonPackedCandMatchNtuplizer::extractPackedTrackKinematics(const pat::PackedCandidate &cand) const {
  if (!cand.hasTrackDetails()) {
    return TrackKinematicsWithErrors();
  }

  try {
    return extractTrackKinematics(cand.pseudoTrack());
  } catch (...) {
    return TrackKinematicsWithErrors();
  }
}

void MuonPackedCandMatchNtuplizer::analyze(const edm::Event &iEvent, const edm::EventSetup &) {
  clearEventBranches();

  edm::Handle<edm::View<pat::Muon>> muons;
  edm::Handle<edm::View<pat::PackedCandidate>> packedCands;
  edm::Handle<reco::VertexCollection> vertices;

  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(packedToken_, packedCands);
  iEvent.getByToken(primaryVertexToken_, vertices);

  run_ = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  nSlimmedMuons_ = static_cast<int>(muons->size());

  const bool hasVertices = vertices.isValid() && !vertices->empty();
  selectedPVIndex_ = hasVertices ? selectPrimaryVertexIndex(*vertices) : kInvalidInt;
  const reco::Vertex *selectedPV =
      (hasVertices && selectedPVIndex_ >= 0 && selectedPVIndex_ < static_cast<int>(vertices->size()))
          ? &vertices->at(selectedPVIndex_)
          : nullptr;

  if (vertices.isValid()) {
    if (storeAllPrimaryVertices_) {
      for (size_t i = 0; i < vertices->size(); ++i) {
        const auto &pv = vertices->at(i);
        pv_x_.push_back(pv.x());
        pv_y_.push_back(pv.y());
        pv_z_.push_back(pv.z());
        pv_xErr_.push_back(pv.xError());
        pv_yErr_.push_back(pv.yError());
        pv_zErr_.push_back(pv.zError());
        pv_chi2_.push_back(pv.chi2());
        pv_ndof_.push_back(pv.ndof());
        pv_sumPt2_.push_back(vertexSumPt2(pv));
        pv_nTracks_.push_back(static_cast<int>(pv.tracksSize()));
        pv_isValid_.push_back((!pv.isFake() && pv.ndof() > 0.) ? 1 : 0);
      }
    } else if (selectedPV != nullptr) {
      pv_x_.push_back(selectedPV->x());
      pv_y_.push_back(selectedPV->y());
      pv_z_.push_back(selectedPV->z());
      pv_xErr_.push_back(selectedPV->xError());
      pv_yErr_.push_back(selectedPV->yError());
      pv_zErr_.push_back(selectedPV->zError());
      pv_chi2_.push_back(selectedPV->chi2());
      pv_ndof_.push_back(selectedPV->ndof());
      pv_sumPt2_.push_back(vertexSumPt2(*selectedPV));
      pv_nTracks_.push_back(static_cast<int>(selectedPV->tracksSize()));
      pv_isValid_.push_back((!selectedPV->isFake() && selectedPV->ndof() > 0.) ? 1 : 0);
    }
  }
  nPV_ = static_cast<int>(pv_x_.size());

  std::vector<int> legacyAvailable;
  legacyAvailable.reserve(packedCands->size());
  for (size_t i = 0; i < packedCands->size(); ++i) {
    legacyAvailable.push_back(static_cast<int>(i));
  }

  for (size_t muIdx = 0; muIdx < muons->size(); ++muIdx) {
    const auto &muon = muons->at(muIdx);
    const bool passSelection = muonSelector_(muon);
    if (!studyAllMuons_ && !passSelection) {
      continue;
    }

    PointerSummary pointerSummary;
    pointerSummary.nSourceCandidatePtrs = muon.numberOfSourceCandidatePtrs();
    pointerSummary.hasPfCandidateRef = (storePointerDiagnostics_ && muon.pfCandidateRef().isNonnull()) ? 1 : 0;

    int trackSource = 0;
    const TrackKinematicsWithErrors muTrackKin = extractMuonTrackKinematics(muon, trackSource);
    const bool muHasTrack = (trackSource != 0);
    const double muMomentum = muon.p();
    const int passedCutBasedIdLoose = muon.passed(reco::Muon::CutBasedIdLoose) ? 1 : 0;
    const int passedCutBasedIdMedium = muon.passed(reco::Muon::CutBasedIdMedium) ? 1 : 0;
    const int passedCutBasedIdMediumPrompt = muon.passed(reco::Muon::CutBasedIdMediumPrompt) ? 1 : 0;
    const int passedCutBasedIdTight = muon.passed(reco::Muon::CutBasedIdTight) ? 1 : 0;
    const int passedCutBasedIdGlobalHighPt = muon.passed(reco::Muon::CutBasedIdGlobalHighPt) ? 1 : 0;
    const int passedCutBasedIdTrkHighPt = muon.passed(reco::Muon::CutBasedIdTrkHighPt) ? 1 : 0;
    const int passedCutBasedIdSoft = muon.passed(reco::Muon::SoftCutBasedId) ? 1 : 0;
    const int passedMvaIDwpSoft = muon.passed(reco::Muon::SoftMvaId) ? 1 : 0;
    const int passedMvaIDwpMedium = muon.passed(reco::Muon::MvaIDwpMedium) ? 1 : 0;
    const int passedMvaIDwpTight = muon.passed(reco::Muon::MvaIDwpTight) ? 1 : 0;
    const bool debugSoftMuon = debugUnmatchedSoftMuons_ && storePointerDiagnostics_ && (passedCutBasedIdSoft == 1);

    const float muDzSelectedPV =
        (selectedPV != nullptr && muHasTrack && std::isfinite(muTrackKin.dzErr))
            ? ((trackSource == 1) ? muon.track()->dz(selectedPV->position()) : muon.muonBestTrack()->dz(selectedPV->position()))
            : kInvalidFloat;
    const float muDxySelectedPV =
        (selectedPV != nullptr && muHasTrack)
            ? ((trackSource == 1) ? muon.track()->dxy(selectedPV->position())
                                  : muon.muonBestTrack()->dxy(selectedPV->position()))
            : kInvalidFloat;

    const auto pfPtr = (storePointerDiagnostics_ && muon.pfCandidateRef().isNonnull())
                           ? edm::refToPtr(muon.pfCandidateRef())
                           : reco::CandidatePtr();

    std::vector<LocalCandidateStudy> retainedRows;
    retainedRows.reserve(16);
    std::vector<DebugCandidateRank> debugRanks;
    if (debugSoftMuon) {
      debugRanks.reserve(packedCands->size());
    }

    int nPassLegacyBox = 0;
    int nPassVector = 0;
    int nPassChi2 = 0;
    int nPassDzPv = 0;
    int nPassDzAssoc = 0;

    int bestVectorIdx = kInvalidInt;
    double bestVectorValue = std::numeric_limits<double>::infinity();
    int bestChi2Idx = kInvalidInt;
    double bestChi2Value = std::numeric_limits<double>::infinity();
    int bestDzPvIdx = kInvalidInt;
    double bestDzPvValue = std::numeric_limits<double>::infinity();
    int bestDzAssocIdx = kInvalidInt;
    double bestDzAssocValue = std::numeric_limits<double>::infinity();

    for (size_t packedIdx = 0; packedIdx < packedCands->size(); ++packedIdx) {
      const auto &cand = packedCands->at(packedIdx);

      int pointerMatchSource = 0;
      int pointerMatchPfRef = 0;
      if (storePointerDiagnostics_) {
        const auto packedPtr = packedCands->ptrAt(packedIdx);
        for (size_t srcIdx = 0; srcIdx < muon.numberOfSourceCandidatePtrs(); ++srcIdx) {
          const auto sourcePtr = muon.sourceCandidatePtr(srcIdx);
          if (sourcePtr.isNonnull() && sourcePtr.id() == packedPtr.id() && sourcePtr.key() == packedPtr.key()) {
            pointerMatchSource = 1;
            pointerSummary.sourceMatches.insert(static_cast<int>(packedIdx));
            pointerSummary.allMatches.insert(static_cast<int>(packedIdx));
          }
        }
        if (pfPtr.isNonnull() && pfPtr.id() == packedPtr.id() && pfPtr.key() == packedPtr.key()) {
          pointerMatchPfRef = 1;
          pointerSummary.pfMatches.insert(static_cast<int>(packedIdx));
          pointerSummary.allMatches.insert(static_cast<int>(packedIdx));
        }
      }

      MatchMetrics metrics;
      const bool chargePass = (!requireChargeMatch_ || cand.charge() == muon.charge());
      const bool legacyChargePass = (cand.charge() == muon.charge());
      if (muMomentum > 0.) {
        const double dpx = cand.px() - muon.px();
        const double dpy = cand.py() - muon.py();
        const double dpz = cand.pz() - muon.pz();
        metrics.deltaPRelVec = std::sqrt(dpx * dpx + dpy * dpy + dpz * dpz) / muMomentum;
      }

      if (debugSoftMuon) {
        DebugCandidateRank rank;
        rank.packedIndex = static_cast<int>(packedIdx);
        rank.deltaPRelVec = metrics.deltaPRelVec;
        rank.chargeMatch = (cand.charge() == muon.charge()) ? 1 : 0;
        rank.pointerMatchSource = pointerMatchSource;
        rank.pointerMatchPfRef = pointerMatchPfRef;
        rank.trackKin = extractPackedTrackKinematics(cand);
        if (muTrackKin.hasTrackErrors) {
          const double dpx = cand.px() - muon.px();
          const double dpy = cand.py() - muon.py();
          const double dpz = cand.pz() - muon.pz();
          rank.rawMomentumChi2 =
              dpx * dpx / sqr(muTrackKin.pxErr) + dpy * dpy / sqr(muTrackKin.pyErr) + dpz * dpz / sqr(muTrackKin.pzErr);
        }
        debugRanks.push_back(rank);
      }

      if (muMomentum > 0. && muon.track().isNonnull() && legacyChargePass &&
          std::abs(cand.px() - muon.px()) < legacyBoxThreshold_ * muMomentum &&
          std::abs(cand.py() - muon.py()) < legacyBoxThreshold_ * muMomentum &&
          std::abs(cand.pz() - muon.pz()) < legacyBoxThreshold_ * muMomentum) {
        metrics.legacyBoxPass = 1;
        ++nPassLegacyBox;
      }

      if (chargePass && metrics.deltaPRelVec >= 0.f && metrics.deltaPRelVec < vectorRelPThreshold_) {
        metrics.vectorPass = 1;
        ++nPassVector;
        if (metrics.deltaPRelVec < bestVectorValue) {
          bestVectorValue = metrics.deltaPRelVec;
          bestVectorIdx = static_cast<int>(packedIdx);
        }
      }

      if (chargePass && muTrackKin.hasTrackErrors) {
        const double dpx = cand.px() - muon.px();
        const double dpy = cand.py() - muon.py();
        const double dpz = cand.pz() - muon.pz();
        metrics.chi2MomentumDiag =
            dpx * dpx / sqr(muTrackKin.pxErr) + dpy * dpy / sqr(muTrackKin.pyErr) + dpz * dpz / sqr(muTrackKin.pzErr);
        if (metrics.chi2MomentumDiag < momentumChi2Threshold_) {
          metrics.chi2Pass = 1;
          ++nPassChi2;
          if (metrics.chi2MomentumDiag < bestChi2Value) {
            bestChi2Value = metrics.chi2MomentumDiag;
            bestChi2Idx = static_cast<int>(packedIdx);
          }
        }
      }

      const bool selectedPvAvailable = (selectedPV != nullptr && muHasTrack && muTrackKin.dzErr > 0.f);
      if (selectedPvAvailable) {
        try {
          const double candDzSelectedPV = cand.dz(selectedPV->position());
          metrics.deltaDzSelectedPV = std::abs(candDzSelectedPV - muDzSelectedPV);
          if (metrics.chi2MomentumDiag >= 0.f) {
            metrics.chi2MomentumDzSelectedPV =
                metrics.chi2MomentumDiag + sqr(metrics.deltaDzSelectedPV / muTrackKin.dzErr);
            if (chargePass && metrics.chi2MomentumDzSelectedPV < momentumDzPvChi2Threshold_) {
              metrics.dzPvPass = 1;
              ++nPassDzPv;
              if (metrics.chi2MomentumDzSelectedPV < bestDzPvValue) {
                bestDzPvValue = metrics.chi2MomentumDzSelectedPV;
                bestDzPvIdx = static_cast<int>(packedIdx);
              }
            }
          }
        } catch (...) {
        }
      }

      const auto vertexRef = cand.vertexRef();
      const bool assocPvAvailable = (vertexRef.isNonnull() && vertices.isValid() &&
                                     vertexRef.key() < vertices->size() && muHasTrack && muTrackKin.dzErr > 0.f);
      if (assocPvAvailable) {
        try {
          const auto &assocPV = vertices->at(vertexRef.key());
          const double muDzAssocPV = (trackSource == 1) ? muon.track()->dz(assocPV.position())
                                                        : muon.muonBestTrack()->dz(assocPV.position());
          metrics.deltaDzAssocPV = std::abs(cand.dzAssociatedPV() - muDzAssocPV);
          if (metrics.chi2MomentumDiag >= 0.f) {
            metrics.chi2MomentumDzAssocPV =
                metrics.chi2MomentumDiag + sqr(metrics.deltaDzAssocPV / muTrackKin.dzErr);
            if (chargePass && metrics.chi2MomentumDzAssocPV < momentumDzAssocChi2Threshold_) {
              metrics.dzAssocPass = 1;
              ++nPassDzAssoc;
              if (metrics.chi2MomentumDzAssocPV < bestDzAssocValue) {
                bestDzAssocValue = metrics.chi2MomentumDzAssocPV;
                bestDzAssocIdx = static_cast<int>(packedIdx);
              }
            }
          }
        } catch (...) {
        }
      }

      if (metrics.legacyBoxPass || metrics.vectorPass || metrics.chi2Pass || metrics.dzPvPass ||
          metrics.dzAssocPass || pointerMatchSource || pointerMatchPfRef) {
        LocalCandidateStudy study;
        study.packedIndex = static_cast<int>(packedIdx);
        study.pointerMatchSource = pointerMatchSource;
        study.pointerMatchPfRef = pointerMatchPfRef;
        study.pointerMatchAny = (pointerMatchSource || pointerMatchPfRef) ? 1 : 0;
        study.trackKin = extractPackedTrackKinematics(cand);
        study.metrics = metrics;
        retainedRows.push_back(study);
      }
    }

    pointerSummary.sourceCandidateResolvedCount = static_cast<int>(pointerSummary.sourceMatches.size());
    pointerSummary.pfCandidateRefResolved = pointerSummary.pfMatches.empty() ? 0 : 1;
    pointerSummary.pointerResolvedCount = static_cast<int>(pointerSummary.allMatches.size());
    pointerSummary.pointerMultiplicityAnomaly = (pointerSummary.pointerResolvedCount > 1) ? 1 : 0;
    if (!pointerSummary.allMatches.empty()) {
      pointerSummary.firstPointerPackedIdx = *pointerSummary.allMatches.begin();
    }

    if (debugSoftMuon && pointerSummary.pointerResolvedCount == 0) {
      const auto debugLess = [](const DebugCandidateRank &lhs, const DebugCandidateRank &rhs) {
        const bool lhsValidChi2 = std::isfinite(lhs.rawMomentumChi2);
        const bool rhsValidChi2 = std::isfinite(rhs.rawMomentumChi2);
        if (lhsValidChi2 != rhsValidChi2) {
          return lhsValidChi2;
        }
        if (lhsValidChi2 && lhs.rawMomentumChi2 != rhs.rawMomentumChi2) {
          return lhs.rawMomentumChi2 < rhs.rawMomentumChi2;
        }

        const bool lhsValidVector = std::isfinite(lhs.deltaPRelVec) && lhs.deltaPRelVec >= 0.f;
        const bool rhsValidVector = std::isfinite(rhs.deltaPRelVec) && rhs.deltaPRelVec >= 0.f;
        if (lhsValidVector != rhsValidVector) {
          return lhsValidVector;
        }
        if (lhsValidVector && lhs.deltaPRelVec != rhs.deltaPRelVec) {
          return lhs.deltaPRelVec < rhs.deltaPRelVec;
        }
        return lhs.packedIndex < rhs.packedIndex;
      };

      const size_t topN = std::min<size_t>(5, debugRanks.size());
      if (topN > 0) {
        std::partial_sort(debugRanks.begin(), debugRanks.begin() + topN, debugRanks.end(), debugLess);
      }

      std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon] run:lumi:event=" << run_ << ":" << lumi_ << ":" << event_
                << " muIdx=" << muIdx
                << " q=" << muon.charge()
                << " pt=" << formatFloatForDebug(muon.pt())
                << " eta=" << formatFloatForDebug(muon.eta())
                << " phi=" << formatFloatForDebug(muon.phi())
                << " ptErr=" << formatFloatForDebug(muTrackKin.ptErr)
                << " etaErr=" << formatFloatForDebug(muTrackKin.etaErr)
                << " phiErr=" << formatFloatForDebug(muTrackKin.phiErr)
                << " trackSource=" << trackSource
                << " cutBasedIds={Loose:" << passedCutBasedIdLoose
                << ",Medium:" << passedCutBasedIdMedium
                << ",MediumPrompt:" << passedCutBasedIdMediumPrompt
                << ",Tight:" << passedCutBasedIdTight
                << ",GlobalHighPt:" << passedCutBasedIdGlobalHighPt
                << ",TrkHighPt:" << passedCutBasedIdTrkHighPt
                << ",Soft:" << passedCutBasedIdSoft
                << "}"
                << " mva={Soft:" << passedMvaIDwpSoft
                << ",Medium:" << passedMvaIDwpMedium
                << ",Tight:" << passedMvaIDwpTight
                << ",Value:" << formatFloatForDebug(muon.mvaIDValue())
                << "}"
                << " ptrDiag={nSourcePtrs:" << pointerSummary.nSourceCandidatePtrs
                << ",hasPfRef:" << pointerSummary.hasPfCandidateRef
                << ",srcResolved:" << pointerSummary.sourceCandidateResolvedCount
                << ",pfResolved:" << pointerSummary.pfCandidateRefResolved
                << "}"
                << '\n';

      if (!muTrackKin.hasTrackErrors) {
        std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][TopChi2] raw momentum chi2 unavailable; using deltaPRelVec tie-break ordering"
                  << '\n';
      }

      for (size_t rankIdx = 0; rankIdx < topN; ++rankIdx) {
        const auto &rank = debugRanks[rankIdx];
        const auto &cand = packedCands->at(rank.packedIndex);
        std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][TopChi2] rank=" << (rankIdx + 1)
                  << " packedIdx=" << rank.packedIndex
                  << " q=" << cand.charge()
                  << " pdgId=" << cand.pdgId()
                  << " pt=" << formatFloatForDebug(cand.pt())
                  << " eta=" << formatFloatForDebug(cand.eta())
                  << " phi=" << formatFloatForDebug(cand.phi())
                  << " ptErr=" << formatFloatForDebug(rank.trackKin.ptErr)
                  << " etaErr=" << formatFloatForDebug(rank.trackKin.etaErr)
                  << " phiErr=" << formatFloatForDebug(rank.trackKin.phiErr)
                  << " rawMomentumChi2=" << formatFloatForDebug(rank.rawMomentumChi2)
                  << " deltaPRelVec=" << formatFloatForDebug(rank.deltaPRelVec)
                  << " hasTrackDetails=" << (cand.hasTrackDetails() ? 1 : 0)
                  << " trackHighPurity=" << (cand.trackHighPurity() ? 1 : 0)
                  << " chargeMatch=" << rank.chargeMatch
                  << " pointerMatchSource=" << rank.pointerMatchSource
                  << " pointerMatchPfRef=" << rank.pointerMatchPfRef
                  << '\n';
      }

      const edm::ProductID packedCollectionId = packedCands->empty() ? edm::ProductID() : packedCands->ptrAt(0).id();
      if (muon.numberOfSourceCandidatePtrs() == 0) {
        std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] no sourceCandidatePtr entries on muon" << '\n';
      } else {
        for (size_t srcIdx = 0; srcIdx < muon.numberOfSourceCandidatePtrs(); ++srcIdx) {
          const auto sourcePtr = muon.sourceCandidatePtr(srcIdx);
          const bool sourcePtrNonnull = sourcePtr.isNonnull();
          const bool sourcePtrAvailable = sourcePtrNonnull && sourcePtr.isAvailable();
          const bool idMatchesPackedCollection = sourcePtrNonnull && !packedCands->empty() && sourcePtr.id() == packedCollectionId;
          const bool keyInRange = sourcePtrNonnull && sourcePtr.key() < packedCands->size();

          std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] srcIdx=" << srcIdx
                    << " nonnull=" << (sourcePtrNonnull ? 1 : 0)
                    << " available=" << (sourcePtrAvailable ? 1 : 0)
                    << " id=" << formatProductIdForDebug(sourcePtrNonnull ? sourcePtr.id() : edm::ProductID(), sourcePtrNonnull)
                    << " key=" << (sourcePtrNonnull ? std::to_string(sourcePtr.key()) : std::string("invalid"))
                    << " packedCollectionId=" << formatProductIdForDebug(packedCollectionId, !packedCands->empty())
                    << '\n';

          if (!sourcePtrNonnull) {
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] classification=null sourceCandidatePtr" << '\n';
            continue;
          }

          if (!idMatchesPackedCollection && !keyInRange) {
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] classification=id mismatch with packedPFCandidates collection and key out of range"
                      << '\n';
            continue;
          }

          if (!idMatchesPackedCollection) {
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] classification=id mismatch with packedPFCandidates collection"
                      << '\n';
          } else if (!keyInRange) {
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] classification=key out of range for packedPFCandidates"
                      << '\n';
          } else {
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] classification=id matches and key points to packed candidate"
                      << '\n';
          }

          if (keyInRange) {
            const auto packedPtrAtKey = packedCands->ptrAt(sourcePtr.key());
            const auto &candAtKey = packedCands->at(sourcePtr.key());
            std::cout << "[MuonPackedCandMatch][UnmatchedSoftMuon][SourcePtr] packedAtKey"
                      << " packedIdx=" << sourcePtr.key()
                      << " id=" << packedPtrAtKey.id()
                      << " key=" << packedPtrAtKey.key()
                      << " q=" << candAtKey.charge()
                      << " pdgId=" << candAtKey.pdgId()
                      << " pt=" << formatFloatForDebug(candAtKey.pt())
                      << " eta=" << formatFloatForDebug(candAtKey.eta())
                      << " phi=" << formatFloatForDebug(candAtKey.phi())
                      << '\n';
          }
        }
      }
    }

    int legacyWinnerIdx = kInvalidInt;
    if (muon.track().isNonnull()) {
      for (auto availIt = legacyAvailable.begin(); availIt != legacyAvailable.end(); ++availIt) {
        const auto &cand = packedCands->at(*availIt);
        if (cand.charge() != muon.charge()) {
          continue;
        }
        if (muMomentum > 0. && std::abs(cand.px() - muon.px()) < legacyBoxThreshold_ * muMomentum &&
            std::abs(cand.py() - muon.py()) < legacyBoxThreshold_ * muMomentum &&
            std::abs(cand.pz() - muon.pz()) < legacyBoxThreshold_ * muMomentum) {
          legacyWinnerIdx = *availIt;
          legacyAvailable.erase(availIt);
          break;
        }
      }
    }

    mu_index_.push_back(static_cast<int>(muIdx));
    mu_passSelection_.push_back(passSelection ? 1 : 0);
    mu_charge_.push_back(muon.charge());
    mu_hasTrack_.push_back(muHasTrack ? 1 : 0);
    mu_trackSource_.push_back(trackSource);
    mu_pt_.push_back(muon.pt());
    mu_eta_.push_back(muon.eta());
    mu_phi_.push_back(muon.phi());
    mu_px_.push_back(muon.px());
    mu_py_.push_back(muon.py());
    mu_pz_.push_back(muon.pz());
    mu_ptErr_.push_back(muTrackKin.ptErr);
    mu_etaErr_.push_back(muTrackKin.etaErr);
    mu_phiErr_.push_back(muTrackKin.phiErr);
    mu_dxyErr_.push_back(muTrackKin.dxyErr);
    mu_dzErr_.push_back(muTrackKin.dzErr);
    mu_pxErr_.push_back(muTrackKin.pxErr);
    mu_pyErr_.push_back(muTrackKin.pyErr);
    mu_pzErr_.push_back(muTrackKin.pzErr);
    mu_dzSelectedPV_.push_back(muDzSelectedPV);
    mu_dxySelectedPV_.push_back(muDxySelectedPV);
    mu_passedCutBasedIdLoose_.push_back(passedCutBasedIdLoose);
    mu_passedCutBasedIdMedium_.push_back(passedCutBasedIdMedium);
    mu_passedCutBasedIdMediumPrompt_.push_back(passedCutBasedIdMediumPrompt);
    mu_passedCutBasedIdTight_.push_back(passedCutBasedIdTight);
    mu_passedCutBasedIdGlobalHighPt_.push_back(passedCutBasedIdGlobalHighPt);
    mu_passedCutBasedIdTrkHighPt_.push_back(passedCutBasedIdTrkHighPt);
    mu_passedCutBasedIdSoft_.push_back(passedCutBasedIdSoft);
    mu_mvaIDValue_.push_back(muon.mvaIDValue());
    mu_passedMvaIDwpSoft_.push_back(passedMvaIDwpSoft);
    mu_passedMvaIDwpMedium_.push_back(passedMvaIDwpMedium);
    mu_passedMvaIDwpTight_.push_back(passedMvaIDwpTight);
    mu_nSourceCandidatePtrs_.push_back(pointerSummary.nSourceCandidatePtrs);
    mu_hasPfCandidateRef_.push_back(pointerSummary.hasPfCandidateRef);
    mu_pfCandidateRefResolved_.push_back(pointerSummary.pfCandidateRefResolved);
    mu_sourceCandidateResolvedCount_.push_back(pointerSummary.sourceCandidateResolvedCount);
    mu_pointerResolvedCount_.push_back(pointerSummary.pointerResolvedCount);
    mu_pointerMultiplicityAnomaly_.push_back(pointerSummary.pointerMultiplicityAnomaly);
    mu_matchLegacyPackedIdx_.push_back(legacyWinnerIdx);
    mu_matchVectorPackedIdx_.push_back(bestVectorIdx);
    mu_matchChi2PackedIdx_.push_back(bestChi2Idx);
    mu_matchDzPvPackedIdx_.push_back(bestDzPvIdx);
    mu_matchDzAssocPackedIdx_.push_back(bestDzAssocIdx);
    mu_matchPointerPackedIdx_.push_back(pointerSummary.firstPointerPackedIdx);
    mu_nPassLegacyBox_.push_back(nPassLegacyBox);
    mu_nPassVector_.push_back(nPassVector);
    mu_nPassChi2_.push_back(nPassChi2);
    mu_nPassDzPv_.push_back(nPassDzPv);
    mu_nPassDzAssoc_.push_back(nPassDzAssoc);
    ++nMuStored_;

    for (const auto &retained : retainedRows) {
      const auto &cand = packedCands->at(retained.packedIndex);
      cand_run_ = run_;
      cand_lumi_ = lumi_;
      cand_event_ = event_;
      cand_selectedPVIndex_ = selectedPVIndex_;
      cand_muonIndex_ = static_cast<int>(muIdx);
      cand_packedIndex_ = retained.packedIndex;
      cand_muPassSelection_ = passSelection ? 1 : 0;
      cand_muCharge_ = muon.charge();
      cand_muPointerResolvedCount_ = pointerSummary.pointerResolvedCount;
      cand_muPt_ = muon.pt();
      cand_muEta_ = muon.eta();
      cand_muPhi_ = muon.phi();
      cand_muPx_ = muon.px();
      cand_muPy_ = muon.py();
      cand_muPz_ = muon.pz();
      cand_muPtErr_ = muTrackKin.ptErr;
      cand_muPxErr_ = muTrackKin.pxErr;
      cand_muPyErr_ = muTrackKin.pyErr;
      cand_muPzErr_ = muTrackKin.pzErr;
      cand_muDzSelectedPV_ = muDzSelectedPV;
      cand_muDxySelectedPV_ = muDxySelectedPV;
      cand_muDzErr_ = muTrackKin.dzErr;

      cand_charge_ = cand.charge();
      cand_pdgId_ = cand.pdgId();
      cand_hasTrackDetails_ = cand.hasTrackDetails() ? 1 : 0;
      cand_trackHighPurity_ = cand.trackHighPurity() ? 1 : 0;
      cand_vertexRefKey_ = cand.vertexRef().isNonnull() ? static_cast<int>(cand.vertexRef().key()) : kInvalidInt;
      cand_pvAssocQuality_ = cand.vertexRef().isNonnull() ? static_cast<int>(cand.pvAssociationQuality()) : kInvalidInt;
      cand_fromPV_ = cand.vertexRef().isNonnull() ? static_cast<int>(cand.fromPV()) : kInvalidInt;
      cand_pt_ = cand.pt();
      cand_eta_ = cand.eta();
      cand_phi_ = cand.phi();
      cand_px_ = cand.px();
      cand_py_ = cand.py();
      cand_pz_ = cand.pz();
      cand_mass_ = cand.mass();
      cand_ptErr_ = retained.trackKin.ptErr;
      cand_etaErr_ = retained.trackKin.etaErr;
      cand_phiErr_ = retained.trackKin.phiErr;
      cand_dxyErr_ = retained.trackKin.dxyErr;
      cand_dzErr_ = retained.trackKin.dzErr;
      cand_pxErr_ = retained.trackKin.pxErr;
      cand_pyErr_ = retained.trackKin.pyErr;
      cand_pzErr_ = retained.trackKin.pzErr;
      cand_dzAssociatedPV_ = cand.vertexRef().isNonnull() ? cand.dzAssociatedPV() : kInvalidFloat;
      cand_dzSelectedPV_ = (selectedPV != nullptr) ? cand.dz(selectedPV->position()) : kInvalidFloat;
      cand_dxyAssociatedPV_ = cand.vertexRef().isNonnull() ? cand.dxy() : kInvalidFloat;
      cand_dxySelectedPV_ = (selectedPV != nullptr) ? cand.dxy(selectedPV->position()) : kInvalidFloat;
      cand_deltaPRelVec_ = retained.metrics.deltaPRelVec;
      cand_chi2MomentumDiag_ = retained.metrics.chi2MomentumDiag;
      cand_deltaDzSelectedPV_ = retained.metrics.deltaDzSelectedPV;
      cand_deltaDzAssocPV_ = retained.metrics.deltaDzAssocPV;
      cand_chi2MomentumDzSelectedPV_ = retained.metrics.chi2MomentumDzSelectedPV;
      cand_chi2MomentumDzAssocPV_ = retained.metrics.chi2MomentumDzAssocPV;
      cand_legacyBoxPass_ = retained.metrics.legacyBoxPass;
      cand_vectorPass_ = retained.metrics.vectorPass;
      cand_chi2Pass_ = retained.metrics.chi2Pass;
      cand_dzPvPass_ = retained.metrics.dzPvPass;
      cand_dzAssocPass_ = retained.metrics.dzAssocPass;
      cand_finalLegacy_ = (retained.packedIndex == legacyWinnerIdx) ? 1 : 0;
      cand_finalVector_ = (retained.packedIndex == bestVectorIdx) ? 1 : 0;
      cand_finalChi2_ = (retained.packedIndex == bestChi2Idx) ? 1 : 0;
      cand_finalDzPv_ = (retained.packedIndex == bestDzPvIdx) ? 1 : 0;
      cand_finalDzAssoc_ = (retained.packedIndex == bestDzAssocIdx) ? 1 : 0;
      cand_pointerMatchAny_ = retained.pointerMatchAny;
      cand_pointerMatchSource_ = retained.pointerMatchSource;
      cand_pointerMatchPfRef_ = retained.pointerMatchPfRef;
      cand_finalPointer_ = retained.pointerMatchAny;
      candidateTree_->Fill();
    }
  }

  eventTree_->Fill();
}

DEFINE_FWK_MODULE(MuonPackedCandMatchNtuplizer);
