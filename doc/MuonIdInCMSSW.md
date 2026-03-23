# Muon ID Internals In CMSSW

This note records the CMSSW-side findings that matter for reading muon quality information from `slimmedMuons` in `MuonPackedCandMatchNtuplizer`.

## Where muon IDs live in CMSSW

### `reco::Muon::Selector` bits

The modern cut-based, isolation, timing, trigger, and some MVA working points are represented as selector bits on `reco::Muon` in:

- `DataFormats/MuonReco/interface/Muon.h`
- `DataFormats/MuonReco/interface/MuonSelectors.h`

Representative enum entries include:

```cpp
CutBasedIdLoose = 0,
CutBasedIdMedium = 1,
CutBasedIdMediumPrompt = 2,
CutBasedIdTight = 3,
CutBasedIdGlobalHighPt = 4,
CutBasedIdTrkHighPt = 5,
...
SoftCutBasedId = 13,
SoftMvaId = 14,
...
MvaIDwpMedium = 34,
MvaIDwpTight = 35,
```

The generic access pattern is:

```cpp
mu.passed(reco::Muon::CutBasedIdLoose);
mu.passed(reco::Muon::CutBasedIdMedium);
mu.passed(reco::Muon::CutBasedIdMediumPrompt);
mu.passed(reco::Muon::CutBasedIdTight);
mu.passed(reco::Muon::SoftCutBasedId);
mu.passed(reco::Muon::SoftMvaId);
```

### `pat::Muon` accessors

`pat::Muon` inherits the selector-bit interface and adds PAT-level score accessors in:

- `DataFormats/PatCandidates/interface/Muon.h`
- `DataFormats/PatCandidates/src/Muon.cc`

Relevant score accessors are:

```cpp
float mvaIDValue() const;
float softMvaValue() const;
float softMvaRun3Value() const;
```

The current prompt-MVA quantities used in this package are:

```cpp
mu.mvaIDValue();
mu.passed(reco::Muon::MvaIDwpMedium);
mu.passed(reco::Muon::MvaIDwpTight);
```

### Legacy string-based selector family

`pat::Muon::muonID(const std::string&)` is a different API. It maps strings such as `TMLastStationLoose` and `GlobalMuonPromptTight` to the older `muon::SelectionType` enum.

```cpp
mu.muonID("TMLastStationLoose");
```

This is not the same mechanism as the modern `CutBasedId*` selector bits.

## Difference between `CutBasedId`, `MvaId`, and `LowPtMva`

### `CutBasedId*`

`CutBasedId*` working points are boolean selector bits set from explicit logic in:

- `DataFormats/MuonReco/src/MuonSelectors.cc`

The key entry point is `muon::makeSelectorBitset(...)`. In CMSSW 15_0_15 it explicitly sets the cut-based muon ID bits, for example:

```cpp
if (muon::isLooseMuon(muon))
  selectors |= (1UL << reco::Muon::CutBasedIdLoose);
if (vertex) {
  if (muon::isTightMuon(muon, *vertex))
    selectors |= (1UL << reco::Muon::CutBasedIdTight);
  if (muon::isHighPtMuon(muon, *vertex))
    selectors |= (1UL << reco::Muon::CutBasedIdGlobalHighPt);
  if (muon::isTrackerHighPtMuon(muon, *vertex))
    selectors |= (1UL << reco::Muon::CutBasedIdTrkHighPt);
}
if (muon::isMediumMuon(muon, run2016_hip_mitigation)) {
  selectors |= (1UL << reco::Muon::CutBasedIdMedium);
  if (vertex && std::abs(muon.muonBestTrack()->dz(vertex->position())) < 0.1 &&
      std::abs(muon.muonBestTrack()->dxy(vertex->position())) < 0.02)
    selectors |= (1UL << reco::Muon::CutBasedIdMediumPrompt);
}
```

So `CutBasedId*` means stored booleans derived from explicit cuts, not ML scores.

### Soft-muon working points

The soft-muon working points split into one cut-based selector and one MVA selector:

- `SoftCutBasedId` is a selector bit filled by `muon::makeSelectorBitset(...)` when `muon::isSoftMuon(muon, *vertex, run2016_hip_mitigation)` passes.
- `SoftMvaId` is a selector bit filled in `PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc` when the PAT producer computes the soft-muon MVA and applies its working-point cut. In CMSSW 15_0_15 that cut is `muon.softMvaValue() > 0.58`.

Representative PAT-side logic:

```cpp
float mva = globalCache()->softMuonMvaEstimator().computeMva(muon);
muon.setSoftMvaValue(mva);
muon.setSelector(reco::Muon::SoftMvaId, muon.softMvaValue() > 0.58);
```

For this package, the working-point bits we store are:

```cpp
mu.passed(reco::Muon::SoftCutBasedId);
mu.passed(reco::Muon::SoftMvaId);
```

### `MvaId*`

In current MiniAOD production, the prompt muon MVA is filled in:

- `PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc`
- `PhysicsTools/PatAlgos/src/MuonMvaIDEstimator.cc`

The producer stores the floating-point score and then assigns working-point selector bits:

```cpp
float mvaID = 0.0;
if (computeMuonIDMVA_) {
  if (muon.isLooseMuon()) {
    mvaID = globalCache()->muonMvaIDEstimator().computeMVAID(muon)[1];
  } else {
    mvaID = -99;
  }
  muon.setMvaIDValue(mvaID);
  muon.setSelector(reco::Muon::MvaIDwpMedium, muon.mvaIDValue() > mvaIDmediumCut_);
  muon.setSelector(reco::Muon::MvaIDwpTight, muon.mvaIDValue() > mvaIDtightCut_);
}
```

The default thresholds configured in `PhysicsTools/PatAlgos/python/producersLayer1/muonProducer_cfi.py` are:

```cpp
mvaIDwpMedium = cms.double(0.08)
mvaIDwpTight = cms.double(0.20)
```

The underlying estimator uses an ONNX model with inputs such as:

- global-muon flag
- valid fraction
- normalized chi2
- local chi2
- kink
- segment compatibility
- valid muon hits
- matched stations
- valid pixel hits
- tracker layers with measurement
- `pt`
- `eta`

So `MvaId*` here means one PAT-stored prompt-muon MVA score plus two prompt-MVA working-point selector bits.

### `LowPtMva`

`LowPtMva` is not a built-in `pat::Muon` data member in standard MiniAOD. In current CMSSW it appears as a separate NanoAOD-era score produced by a dedicated module in:

- `PhysicsTools/NanoAOD/python/muons_cff.py`

The relevant configuration is:

```cpp
muonMVALowPt = muonPROMPTMVA.clone(
    weightFile = cms.FileInPath("PhysicsTools/NanoAOD/data/mu_BDTG_lowpt.weights.xml"),
    name = cms.string("muonMVALowPt"),
    variables = _legacy_muon_BDT_variable,
)
```

NanoAOD then exposes it as an external variable:

```cpp
mvaLowPt = ExtVar(cms.InputTag("muonMVALowPt"), float, doc="Low pt muon ID score")
```

This makes `LowPtMva` a separate score with a different model and feature set, not the same thing as `pat::Muon::mvaIDValue()`.

## What is actually safe to read from `slimmedMuons`

For this package, the robust MiniAOD-native quantities are:

```cpp
mu.passed(reco::Muon::CutBasedIdLoose);
mu.passed(reco::Muon::CutBasedIdMedium);
mu.passed(reco::Muon::CutBasedIdMediumPrompt);
mu.passed(reco::Muon::CutBasedIdTight);
mu.passed(reco::Muon::CutBasedIdGlobalHighPt);
mu.passed(reco::Muon::CutBasedIdTrkHighPt);
mu.passed(reco::Muon::SoftCutBasedId);
mu.passed(reco::Muon::SoftMvaId);

mu.mvaIDValue();
mu.passed(reco::Muon::MvaIDwpMedium);
mu.passed(reco::Muon::MvaIDwpTight);
```

These are the quantities added to the ntuple in this update.

## What we intentionally do not store

This pass intentionally excludes:

- `LowPtMva*`, because it is not a standard built-in `pat::Muon` quantity on `slimmedMuons`
- `softMvaValue()` and `softMvaRun3Value()` scores, because this update adds only the soft working-point bits and not the raw soft-MVA score payload
- legacy string-based IDs such as `mu.muonID("TMLastStationLoose")`
- `MvaLoose`, `MvaMedium`, `MvaTight`, `MvaVTight`, `MvaVVTight`, because the enums exist but the standard MiniAOD producer path inspected here does not populate them as the current prompt-MVA working-point interface

## Code references

CMSSW release inspected: `CMSSW_15_0_15`

Main files examined:

- `DataFormats/MuonReco/interface/Muon.h`
- `DataFormats/MuonReco/interface/MuonSelectors.h`
- `DataFormats/MuonReco/src/MuonSelectors.cc`
- `DataFormats/PatCandidates/interface/Muon.h`
- `DataFormats/PatCandidates/src/Muon.cc`
- `PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc`
- `PhysicsTools/PatAlgos/src/MuonMvaIDEstimator.cc`
- `PhysicsTools/NanoAOD/python/muons_cff.py`
- `PhysicsTools/PatAlgos/python/producersLayer1/muonProducer_cfi.py`
- `PhysicsTools/PatAlgos/python/slimming/miniAOD_tools.py`
