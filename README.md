# MuonPackedPFCandMatch

`MuonPackedCandMatchNtuplizer` is a standalone CMSSW validation `EDAnalyzer` for studying how `pat::Muon` objects in `slimmedMuons` can be associated with `pat::PackedCandidate` objects in `packedPFCandidates` without changing `TPS-Onia2MuMu`.

## Fresh CMSSW Setup

Create a new CMSSW release area, clone this package directly under `src/HeavyFlavorAnalysis/`, and build it there:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_15_0_15
cd CMSSW_15_0_15/src
eval "$(scramv1 runtime -sh)"

mkdir -p HeavyFlavorAnalysis
git clone git@github.com:Eric100911/MuonPackedPFCandMatch.git HeavyFlavorAnalysis/MuonPackedPFCandMatch

scram b -j 4 HeavyFlavorAnalysis/MuonPackedPFCandMatch
```

## Package Layout

- `interface/MuonPackedCandMatchNtuplizer.h`: analyzer interface and ntuple payload definitions
- `src/MuonPackedCandMatchNtuplizer.cc`: analyzer implementation
- `test/runMuonPackedCandMatch_cfg.py`: `cmsRun` entry point with `runOnMC`, `era`, `eventRange`, and debug switches
- `doc/MuonPackedCandMatchStudy.md`: extended technical note and validation details
- `notebooks/MuonPackedCandMatchWorkbook.ipynb`: offline analysis notebook

## Repository Structure

The package is expected to live at `src/HeavyFlavorAnalysis/MuonPackedPFCandMatch` inside the CMSSW release area:

```text
CMSSW_15_0_15/
└── src/
    └── HeavyFlavorAnalysis/
        └── MuonPackedPFCandMatch/
            ├── interface/
            ├── src/
            ├── test/
            ├── doc/
            ├── notebooks/
            └── BuildFile.xml
```

## Matching Views

The analyzer stores diagnostics for six matching strategies:

1. Legacy sequential comparator
2. Vector relative-momentum comparator
3. Cartesian momentum chi2 comparator
4. Momentum-plus-`dz` comparator against the selected primary vertex
5. Momentum-plus-`dz` comparator against the candidate-associated primary vertex
6. Direct PAT pointer diagnostics using `sourceCandidatePtr()` and `pfCandidateRef()`

## Build

From the CMSSW release area:

```bash
cd $CMSSW_BASE/src
scram b -j 4 HeavyFlavorAnalysis/MuonPackedPFCandMatch
```

## Run

From `$CMSSW_BASE/src`.

MC example:

```bash
cmsRun HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py \
    inputFiles=file:/eos/user/c/chiw/JpsiJpsiUps/MC_samples/miniAOD/DPS-JpsiJpsi-Y/filter_JpsiPtMin4p0_YPtMin6p0/HO_DPS_JpsiJpsi_Y_Run3Summer22_miniAOD_292.root \
    outputFile=muonPackedMatch_mc.root \
    runOnMC=True era=Run2022 maxEvents=1000
```

Data example:

```bash
cmsRun HeavyFlavorAnalysis/MuonPackedPFCandMatch/test/runMuonPackedCandMatch_cfg.py \
    inputFiles=/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root \
    outputFile=muonPackedMatch_data.root \
    runOnMC=False era=Run2023D maxEvents=1000
```

`VarParsing('analysis')` may append tags such as `_numEvent1` to the final ROOT filename. For remote data, ensure a valid X509 proxy is available before running, for example with `voms-proxy-info --timeleft`.

## Output

The ntuple is written under the ROOT directory `muonPackedCandMatch`:

- `muonPackedCandMatch/Events`
- `muonPackedCandMatch/MuonCandMatch`

See `doc/MuonPackedCandMatchStudy.md` for validated sample inputs, observed event counts, and the full branch inventory.
