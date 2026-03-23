# MuonPackedCandMatchNtuplizer

`MuonPackedCandMatchNtuplizer` is a standalone validation `EDAnalyzer` for studying how `pat::Muon` objects in `slimmedMuons` can be associated to `pat::PackedCandidate` objects in `packedPFCandidates` without modifying `TPS-Onia2MuMu`.

## Package layout

- Analyzer header: `interface/MuonPackedCandMatchNtuplizer.h`
- Analyzer implementation: `src/MuonPackedCandMatchNtuplizer.cc`
- CMSSW config: `test/runMuonPackedCandMatch_cfg.py`
- Analysis workbook: `notebooks/MuonPackedCandMatchWorkbook.ipynb`
- This note: `doc/MuonPackedCandMatchStudy.md`

## Implemented matching views

The analyzer stores diagnostics for six matching views.

1. Legacy sequential comparator: reproduces the currently active inline matching logic in `MultiLepPAT`.
2. Vector threshold comparator: minimum `|p_cand - p_mu| / |p_mu|`.
3. Momentum chi2 comparator: minimum Cartesian momentum chi2.
4. Momentum+dz comparator with respect to the selected primary vertex.
5. Momentum+dz comparator with respect to the candidate-associated primary vertex.
6. Pointer diagnostics using `sourceCandidatePtr()` and `pfCandidateRef()`.

The candidate tree keeps a packed candidate if it passes at least one preselection flag or if it is directly resolved by a PAT pointer diagnostic.

## Validated debug status

The current implementation has been smoke-tested on one MC event and one data event.

- MC input: `/eos/user/c/chiw/JpsiJpsiUps/MC_samples/miniAOD/DPS-JpsiJpsi-Y/filter_JpsiPtMin4p0_YPtMin6p0/HO_DPS_JpsiJpsi_Y_Run3Summer22_miniAOD_292.root`
- MC run mode: `runOnMC=True era=Run2022 maxEvents=1`
- Output: `test_muonPackedMatch_mc_numEvent1.root`
- Observed content: `muonPackedCandMatch/Events` with 1 entry and `muonPackedCandMatch/MuonCandMatch` with 93 entries

- Data input: `/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root`
- Data run mode: `runOnMC=False era=Run2023D maxEvents=1`
- Output: `test_muonPackedMatch_data_numEvent1.root`
- Observed content: `muonPackedCandMatch/Events` with 1 entry and `muonPackedCandMatch/MuonCandMatch` with 172 entries

Two workflow details are worth calling out explicitly.

- Data access requires a valid X509 user proxy.
- The `TFileService` output lives inside the ROOT directory `muonPackedCandMatch`, so the tree paths are `muonPackedCandMatch/Events` and `muonPackedCandMatch/MuonCandMatch` rather than top-level `Events` and `MuonCandMatch`.

## Build and run

From the CMSSW release area:

```bash
cd /eos/user/c/chiw/JpsiJpsiPhi/CMSSW_15_0_15_JpsiJpsiPhi_refactor
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(scramv1 runtime -sh)"
scram b -j 4
```

### MC example

```bash
cmsRun src/HeavyFlavorAnalysis/TPSMuonPackedMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=file:/eos/user/c/chiw/JpsiJpsiUps/MC_samples/miniAOD/DPS-JpsiJpsi-Y/filter_JpsiPtMin4p0_YPtMin6p0/HO_DPS_JpsiJpsi_Y_Run3Summer22_miniAOD_292.root     outputFile=muonPackedMatch_mc.root     runOnMC=True     era=Run2022     maxEvents=1000
```

### Data example

```bash
cmsRun src/HeavyFlavorAnalysis/TPSMuonPackedMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root     outputFile=muonPackedMatch_data.root     runOnMC=False     era=Run2023D     maxEvents=1000
```

### Single-event debugging

```bash
cmsRun src/HeavyFlavorAnalysis/TPSMuonPackedMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root     outputFile=muonPackedMatch_debug.root     runOnMC=False     era=Run2023D     maxEvents=1     eventRange=369873:446440
```

### Output filename note

Because this cfg uses `VarParsing('analysis')`, the final ROOT filename may receive suffixes when command-line overrides are present. For example,

- `outputFile=test_muonPackedMatch_mc.root maxEvents=1` produced `test_muonPackedMatch_mc_numEvent1.root`
- `outputFile=test_muonPackedMatch_data.root maxEvents=1` produced `test_muonPackedMatch_data_numEvent1.root`

This is expected behavior in the current cfg style inherited from `ConfFile_cfg.py`.

## ROOT output layout

The ntuple is written under the ROOT directory `muonPackedCandMatch`.

- `muonPackedCandMatch/Events`
- `muonPackedCandMatch/MuonCandMatch`

### `Events` tree

One row per event. Muon content is stored as jagged arrays indexed by the per-event muon ordering.

```text
run
lumi
event
selectedPVIndex
nPV
nSlimmedMuons
nMuStored
pv_x
pv_y
pv_z
pv_xErr
pv_yErr
pv_zErr
pv_chi2
pv_ndof
pv_sumPt2
pv_nTracks
pv_isValid
mu_index
mu_passSelection
mu_charge
mu_hasTrack
mu_trackSource
mu_pt
mu_eta
mu_phi
mu_px
mu_py
mu_pz
mu_ptErr
mu_etaErr
mu_phiErr
mu_dxyErr
mu_dzErr
mu_pxErr
mu_pyErr
mu_pzErr
mu_dzSelectedPV
mu_dxySelectedPV
mu_nSourceCandidatePtrs
mu_hasPfCandidateRef
mu_pfCandidateRefResolved
mu_sourceCandidateResolvedCount
mu_pointerResolvedCount
mu_pointerMultiplicityAnomaly
mu_matchLegacyPackedIdx
mu_matchVectorPackedIdx
mu_matchChi2PackedIdx
mu_matchDzPvPackedIdx
mu_matchDzAssocPackedIdx
mu_matchPointerPackedIdx
mu_nPassLegacyBox
mu_nPassVector
mu_nPassChi2
mu_nPassDzPv
mu_nPassDzAssoc
```

### `MuonCandMatch` tree

One row per retained `(event, muonIndex, packedIndex)` pair.

```text
run
lumi
event
selectedPVIndex
muonIndex
packedIndex
mu_passSelection
mu_charge
mu_pointerResolvedCount
mu_pt
mu_eta
mu_phi
mu_px
mu_py
mu_pz
mu_ptErr
mu_pxErr
mu_pyErr
mu_pzErr
mu_dzSelectedPV
mu_dxySelectedPV
mu_dzErr
cand_charge
cand_pdgId
cand_hasTrackDetails
cand_trackHighPurity
cand_vertexRefKey
cand_pvAssocQuality
cand_fromPV
cand_pt
cand_eta
cand_phi
cand_px
cand_py
cand_pz
cand_mass
cand_ptErr
cand_etaErr
cand_phiErr
cand_dxyErr
cand_dzErr
cand_pxErr
cand_pyErr
cand_pzErr
cand_dzAssociatedPV
cand_dzSelectedPV
cand_dxyAssociatedPV
cand_dxySelectedPV
deltaPRelVec
chi2MomentumDiag
deltaDzSelectedPV
deltaDzAssocPV
chi2MomentumDzSelectedPV
chi2MomentumDzAssocPV
legacyBoxPass
vectorPass
chi2Pass
dzPvPass
dzAssocPass
finalLegacy
finalVector
finalChi2
finalDzPv
finalDzAssoc
pointerMatchAny
pointerMatchSource
pointerMatchPfRef
finalPointer
```

## Error extraction and propagation

Muon and packed-candidate native track errors are stored directly through the standard `reco::Track` accessors.

```cpp
ptErr  = track.ptError();
etaErr = track.etaError();
phiErr = track.phiError();
dxyErr = track.dxyError();
dzErr  = track.dzError();
```

Cartesian momentum uncertainties are not taken directly from `track.covariance()(0,0..2)`. The `reco::TrackBase` covariance basis is `(qoverp, lambda, phi, dxy, dsz)`, so the analyzer propagates the `(qoverp, lambda, phi)` block through the Jacobian:

```cpp
const auto C = track.covariance();
const double qop = track.qoverp();
const double lam = track.lambda();
const double phi = track.phi();
const double p = std::abs(1.0 / qop);
const double dpdqop = -std::copysign(1.0, qop) / (qop * qop);

// px = p * cos(lam) * cos(phi)
// py = p * cos(lam) * sin(phi)
// pz = p * sin(lam)
```

```cpp
J = d(px, py, pz) / d(qoverp, lambda, phi)
Cov_xyz = J * Cov_(qoverp, lambda, phi) * J^T
pxErr = sqrt(Cov_xyz(0,0))
pyErr = sqrt(Cov_xyz(1,1))
pzErr = sqrt(Cov_xyz(2,2))
```

For `pat::PackedCandidate`, track-based errors are only filled when `hasTrackDetails()` is true. Otherwise, the analyzer keeps the kinematics and PV-association information but writes sentinel values for track-error quantities.

## Matching metrics stored per candidate row

Each retained candidate row stores the following quantities needed for downstream notebook studies.

- Kinematics and native errors for both the muon and the candidate
- `cand_pdgId`, `cand_fromPV`, `cand_pvAssocQuality`, `cand_vertexRefKey`
- `cand_dzAssociatedPV`, `cand_dzSelectedPV`, `cand_dxyAssociatedPV`, `cand_dxySelectedPV`
- `deltaPRelVec`
- `chi2MomentumDiag`
- `deltaDzSelectedPV`
- `deltaDzAssocPV`
- `chi2MomentumDzSelectedPV`
- `chi2MomentumDzAssocPV`
- Preselection flags: `legacyBoxPass`, `vectorPass`, `chi2Pass`, `dzPvPass`, `dzAssocPass`
- Final-choice flags: `finalLegacy`, `finalVector`, `finalChi2`, `finalDzPv`, `finalDzAssoc`, `finalPointer`
- Pointer diagnostics: `pointerMatchAny`, `pointerMatchSource`, `pointerMatchPfRef`

## Workbook usage

The Jupyter workbook lives at `notebooks/MuonPackedCandMatchWorkbook.ipynb`.

It is designed for a Python environment with:

- `uproot`
- `awkward`
- `numpy`
- `pandas`
- `matplotlib`

> Example kernel availble at `/afs/cern.ch/user/c/chiw/.local/share/jupyter/kernels/eos-jetflav-kernel/kernel.json`.

The notebook is organized to:

1. open the ntuple from `muonPackedCandMatch/Events` and `muonPackedCandMatch/MuonCandMatch`
2. flatten the event-level jagged muon content into a muon table
3. summarize candidate multiplicities per muon for each matching method
4. build agreement matrices among legacy, vector, chi2, dz-selected-PV, dz-associated-PV, and pointer winners
5. inspect pointer-resolved vs pointer-unresolved categories
6. attach winner-level candidate metadata such as `fromPV`, `pvAssociationQuality`, `vertexRefKey`, and `pdgId`
7. produce side-by-side data and MC plots with the same helper functions

The workbook intentionally uses the real tree paths under `muonPackedCandMatch/...`, which was confirmed during the smoke tests.

## Configuration details

The cfg in `test/runMuonPackedCandMatch_cfg.py` follows the same style as `TPS-Onia2MuMu/test/ConfFile_cfg.py` for:

- `VarParsing`
- GlobalTag selection from data/MC era dictionaries
- `TFileService`
- optional `eventRange`

The default data era for the provided file is `Run2023D`, which maps to `130X_dataRun3_PromptAnalysis_v1`. The validated MC era for the provided Summer22 sample is `Run2022`, which maps to `130X_mcRun3_2022_realistic_v5`.
