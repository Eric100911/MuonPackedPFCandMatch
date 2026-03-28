# MuonPackedPFCandMatch

`MuonPackedCandMatchNtuplizer` is a standalone CMSSW validation `EDAnalyzer` for studying how `pat::Muon` objects in `slimmedMuons` can be associated with `pat::PackedCandidate` objects in `packedPFCandidates` without changing `TPS-Onia2MuMu`.

## Package Layout

- `interface/MuonPackedCandMatchNtuplizer.h`: analyzer interface and ntuple payload definitions
- `src/MuonPackedCandMatchNtuplizer.cc`: analyzer implementation
- `test/runMuonPackedCandMatch_cfg.py`: `cmsRun` entry point with `runOnMC`, `era`, `eventRange`, and debug switches
- `doc/MuonPackedCandMatchStudy.md`: extended technical note and validation details
- `notebooks/MuonPackedCandMatchWorkbook.ipynb`: offline analysis notebook

## Matching Views

The analyzer stores diagnostics for six matching strategies:

1. Legacy sequential comparator
2. Vector relative-momentum comparator
3. Cartesian momentum chi2 comparator
4. Momentum-plus-`dz` comparator against the selected primary vertex
5. Momentum-plus-`dz` comparator against the candidate-associated primary vertex
6. Direct PAT pointer diagnostics using `sourceCandidatePtr()` and `pfCandidateRef()`

## Build

From the CMSSW work area:

```bash
cd /eos/home-c/chiw/JpsiJpsiPhi/CMSSW_15_0_15_JpsiJpsiPhi_refactor
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(scramv1 runtime -sh)"
scram b -j 4 HeavyFlavorAnalysis/MuonPackedPFCandMatch
```

## Run

MC example:

```bash
cmsRun src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py \
    inputFiles=file:/eos/user/c/chiw/JpsiJpsiUps/MC_samples/miniAOD/DPS-JpsiJpsi-Y/filter_JpsiPtMin4p0_YPtMin6p0/HO_DPS_JpsiJpsi_Y_Run3Summer22_miniAOD_292.root \
    outputFile=muonPackedMatch_mc.root \
    runOnMC=True era=Run2022 maxEvents=1000
```

Data example:

```bash
cmsRun src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py \
    inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root \
    outputFile=muonPackedMatch_data.root \
    runOnMC=False era=Run2023D maxEvents=1000
```

`VarParsing('analysis')` may append tags such as `_numEvent1` to the final ROOT filename. For remote data, ensure a valid X509 proxy is available before running.

## Output

The ntuple is written under the ROOT directory `muonPackedCandMatch`:

- `muonPackedCandMatch/Events`
- `muonPackedCandMatch/MuonCandMatch`

See `doc/MuonPackedCandMatchStudy.md` for validated sample inputs, observed event counts, and the full branch inventory.
