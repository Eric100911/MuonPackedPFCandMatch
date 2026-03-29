# MuonPackedCandMatchNtuplizer

`MuonPackedCandMatchNtuplizer` is a standalone validation `EDAnalyzer` for studying how `pat::Muon` objects in `slimmedMuons` can be associated to `pat::PackedCandidate` objects in `packedPFCandidates` without modifying `TPS-Onia2MuMu`.

## Package layout

- Analyzer header: `interface/MuonPackedCandMatchNtuplizer.h`
- Analyzer implementation: `src/MuonPackedCandMatchNtuplizer.cc`
- CMSSW config: `test/runMuonPackedCandMatch_cfg.py`
- Analysis workbook: `notebooks/MuonPackedCandMatchWorkbook.ipynb`
- This note: `doc/MuonPackedCandMatchStudy.md`

## Implemented matching views

The analyzer stores diagnostics for five matching views.

1. Vector threshold comparator: minimum $|p_{\mathrm{cand}} - p_{\mu}| / |p_{\mu}|$.
2. Momentum $\chi^2$ comparator: minimum Cartesian momentum $\chi^2$.
3. Momentum+$\Delta z$ comparator with respect to the selected primary vertex.
4. Momentum+$\Delta z$ comparator with respect to the candidate-associated primary vertex.
5. Pointer diagnostics using `sourceCandidatePtr()` and `pfCandidateRef()`.

The event tree stores all studied muons with full kinematics, IDs, pointer diagnostics, and method-summary branches. In the compact default mode, the event tree keeps the total PV count together with only the selected PV payload. The candidate tree keeps retained rows only for disagreement cases, where the enabled methods either pick different packed candidates or mix match and mismatch outcomes. In that mode, the retained set is the union of distinct winners, pointer-resolved candidates, and a configurable per-method loser shortlist. Set `storeDetailedRowsOnlyOnDisagreement=False` in the cfg to recover the broader candidate-row dump.

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
cmsRun src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=file:/eos/user/c/chiw/JpsiJpsiUps/MC_samples/miniAOD/DPS-JpsiJpsi-Y/filter_JpsiPtMin4p0_YPtMin6p0/HO_DPS_JpsiJpsi_Y_Run3Summer22_miniAOD_292.root     outputFile=muonPackedMatch_mc.root     runOnMC=True     era=Run2022     maxEvents=1000
```

### Data example

```bash
cmsRun src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root     outputFile=muonPackedMatch_data.root     runOnMC=False     era=Run2023D     maxEvents=1000
```

### Single-event debugging

```bash
cmsRun src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py     inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root     outputFile=muonPackedMatch_debug.root     runOnMC=False     era=Run2023D     maxEvents=1     eventRange=369873:446440
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

Key branch groups are:

- Event and PV summary: `run`, `lumi`, `event`, `selectedPVIndex`, `nPV`, `nPVStored`, `nSlimmedMuons`, `nMuStored`, plus the `pv_*` arrays. In compact mode, `nPV` is the total reconstructed PV count while `nPVStored` is the number of PV rows actually written to the tree.
- Muon kinematics and propagated uncertainties: `mu_pt`, `mu_eta`, `mu_phi`, `mu_px`, `mu_py`, `mu_pz`, `mu_*Err`, `mu_dzSelectedPV`, `mu_dxySelectedPV`.
- Muon IDs: `mu_passedCutBasedId*`, `mu_mvaIDValue`, `mu_passedMvaIDwp*`.
- Pointer diagnostics: `mu_nSourceCandidatePtrs`, `mu_hasPfCandidateRef`, `mu_pfCandidateRefResolved`, `mu_sourceCandidateResolvedCount`, `mu_pointerResolvedCount`, `mu_pointerMultiplicityAnomaly`.
- Final winner indices for the surviving methods: `mu_matchVectorPackedIdx`, `mu_matchChi2PackedIdx`, `mu_matchDzPvPackedIdx`, `mu_matchDzAssocPackedIdx`, `mu_matchPointerPackedIdx`.
- Per-method pass counts: `mu_nPassVector`, `mu_nPassChi2`, `mu_nPassDzPv`, `mu_nPassDzAssoc`.
- Agreement summary branches: `mu_nMethodsConsidered`, `mu_nMethodsMatched`, `mu_nDistinctMatchedPackedIdx`, `mu_hasAnyMethodMatch`, `mu_isAgreementOnMatch`, `mu_isAgreementOnMismatch`, `mu_hasMethodDisagreement`, `mu_nPackedPassingAnyKinematic`, `mu_nPackedPointerResolved`, `mu_nPackedPassingAnyCriterion`, `mu_bestVectorRelP`, `mu_bestMomentumChi2`, `mu_bestMomentumDzPvChi2`, `mu_bestMomentumDzAssocChi2`.

### `MuonCandMatch` tree

One row per retained `(event, muonIndex, packedIndex)` pair. With the default cfg, these rows are written only for disagreement muons.

For a disagreement muon, the compact retained set is:

- the distinct winning candidates from the vector, momentum-$\chi^2$, selected-PV, and associated-PV methods
- the unique pointer winner, when one exists
- all pointer-resolved candidates when pointer resolution is ambiguous
- up to `LoserRowsPerMethod` non-winning candidates for each considered kinematic method, ranked by the corresponding raw score

Each row stores:

- Muon summary fields duplicated from the event tree: event ids, `muonIndex`, `mu_passSelection`, charge, kinematics, propagated uncertainties, selected-PV impact parameters, and `mu_pointerResolvedCount`.
- Candidate kinematics and metadata: charge, `pdgId`, `hasTrackDetails`, `trackHighPurity`, `vertexRefKey`, `pvAssocQuality`, `fromPV`, momentum, mass, and track-error fields.
- Matching metrics: `deltaPRelVec`, `chi2MomentumDiag`, `deltaDzSelectedPV`, `deltaDzAssocPV`, `chi2MomentumDzSelectedPV`, `chi2MomentumDzAssocPV`.
- Surviving kinematic pass flags: `vectorPass`, `chi2Pass`, `dzPvPass`, `dzAssocPass`.
- Final-choice flags: `finalVector`, `finalChi2`, `finalDzPv`, `finalDzAssoc`, `finalPointer`.
- Pointer diagnostics: `pointerMatchAny`, `pointerMatchSource`, `pointerMatchPfRef`.

## Error extraction and propagation

Muon and packed-candidate native track errors are stored directly through the standard `reco::Track` accessors.

$$
\begin{aligned}
\sigma_{p_T} &= \texttt{track.ptError()} \\
\sigma_{\eta} &= \texttt{track.etaError()} \\
\sigma_{\phi} &= \texttt{track.phiError()} \\
\sigma_{d_{xy}} &= \texttt{track.dxyError()} \\
\sigma_{d_z} &= \texttt{track.dzError()}
\end{aligned}
$$

Cartesian momentum uncertainties are not taken directly from `track.covariance()(0,0..2)`. The `reco::TrackBase` covariance basis is `(qoverp, lambda, phi, dxy, dsz)`, so the analyzer propagates the `(qoverp, lambda, phi)` block through the Jacobian:

$$
\begin{aligned}
p_x &= p \cos \lambda \cos \phi \\
p_y &= p \cos \lambda \sin \phi \\
p_z &= p \sin \lambda
\end{aligned}
$$

$$
\begin{aligned}
J &= \frac{\partial (p_x, p_y, p_z)}{\partial (q/p, \lambda, \phi)} \\
\mathrm{Cov}_{xyz} &= J \, \mathrm{Cov}_{(q/p,\lambda,\phi)} \, J^T \\
\sigma_{p_x} &= \sqrt{\mathrm{Cov}_{xyz}(0,0)} \\
\sigma_{p_y} &= \sqrt{\mathrm{Cov}_{xyz}(1,1)} \\
\sigma_{p_z} &= \sqrt{\mathrm{Cov}_{xyz}(2,2)}
\end{aligned}
$$

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
- Preselection flags: `vectorPass`, `chi2Pass`, `dzPvPass`, `dzAssocPass`
- Final-choice flags: `finalVector`, `finalChi2`, `finalDzPv`, `finalDzAssoc`, `finalPointer`
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
4. build agreement matrices among vector, chi2, dz-selected-PV, dz-associated-PV, and pointer winners
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

The default slimming mode is controlled by `storeDetailedRowsOnlyOnDisagreement`. Leave it at `True` for compact ntuples, or set it to `False` when you need the broader candidate-row dump for debugging.

The disagreement-row selection is controlled by:

- `DisagreementRowMode`, which defaults to `winnersPlusPerMethodLosers`
- `LoserRowsPerMethod`, which defaults to `3`

The default compact cfg also sets `StoreAllPrimaryVertices=False`, so only the selected PV payload is written even though `nPV` still reports the full reconstructed PV multiplicity.

The default data era for the provided file is `Run2023D`, which maps to `130X_dataRun3_PromptAnalysis_v1`. The validated MC era for the provided Summer22 sample is `Run2022`, which maps to `130X_mcRun3_2022_realistic_v5`.
