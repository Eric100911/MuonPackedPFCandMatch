"""Run the MuonPackedCandMatchNtuplizer on MiniAOD inputs.

This cfg selects the GlobalTag from `runOnMC` and `era`, exposes the main
matching-study switches through VarParsing, and writes the ntuple with
TFileService.
"""

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# --- Command-line arguments exposed to cmsRun ---
ivars = VarParsing.VarParsing('analysis')

# These options control dataset interpretation, event filtering, and optional
# debug printouts without changing the analyzer source.
ivars.register(
    'runOnMC',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='Whether to run on MC (default: False)'
)
ivars.register(
    'era',
    default='Run2023D',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Data-taking era for GlobalTag selection'
)
ivars.register(
    'eventRange',
    default='',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Optional event range to process'
)
ivars.register(
    'studyAllMuons',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='Store all slimmedMuons instead of only analysis-selected muons'
)
ivars.register(
    'debugUnmatchedSoftMuons',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='Print stdout debug for SoftCutBasedId muons with zero pointer matches'
)
ivars.register(
    'storeDetailedRowsOnlyOnDisagreement',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info='Store MuonCandMatch rows only for muons where the enabled methods disagree'
)
ivars.register(
    'loserRowsPerMethod',
    default=3,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info='Maximum number of non-winning loser rows to keep per considered kinematic method'
)
ivars.register(
    'disagreementRowMode',
    default='winnersPlusPerMethodLosers',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info='Detailed-row retention mode for disagreement muons'
)

# --- Default input/output values for local testing ---
ivars.inputFiles = ('/store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/369/873/00000/33e0e861-ddbc-4afe-a76b-31be5057dff1.root',)
ivars.outputFile = 'muonPackedCandMatch.root'
# Note: VarParsing('analysis') may append tags such as _numEvent1 to the
# final ROOT file name when command-line overrides are used.
ivars.era = 'Run2023D'

# --- Era-to-GlobalTag lookup for data and MC ---
globalTagDict = {
    'data': {
        'Run2022C': '124X_dataRun3_PromptAnalysis_v1',
        'Run2022D': '124X_dataRun3_PromptAnalysis_v1',
        'Run2022E': '124X_dataRun3_Prompt_v10',
        'Run2022F': '124X_dataRun3_PromptAnalysis_v2',
        'Run2022G': '124X_dataRun3_PromptAnalysis_v2',
        'Run2023C': '130X_dataRun3_PromptAnalysis_v1',
        'Run2023D': '130X_dataRun3_PromptAnalysis_v1',
        'Run2024C': '150X_dataRun3_v2',
        'Run2024D': '150X_dataRun3_v2',
        'Run2024E': '150X_dataRun3_v2',
        'Run2024F': '150X_dataRun3_v2',
        'Run2024G': '150X_dataRun3_v2',
        'Run2024H': '150X_dataRun3_v2',
        'Run2024I': '150X_dataRun3_v2',
        'Run2025C': '150X_dataRun3_v2',
        'Run2025D': '150X_dataRun3_v2',
        'Run2025E': '150X_dataRun3_v2',
        'Run2025F': '150X_dataRun3_v2',
        'Run2025G': '150X_dataRun3_v2',
    },
    'MC': {
        'Run2022': '130X_mcRun3_2022_realistic_v5',
        'Run2022EE': '130X_mcRun3_2022_realistic_postEE_v6',
        'Run2023': '130X_mcRun3_2023_realistic_v14',
        'Run2023BPix': '130X_mcRun3_2023_realistic_postBPix_v2',
        'Run2024': '150X_mcRun3_2024_realistic_v2',
        'Run2025': '150X_mcRun3_2024_realistic_v2',
    }
}

# --- Parse and validate runtime arguments ---
ivars.parseArguments()

if ivars.runOnMC and ivars.era not in globalTagDict['MC']:
    raise ValueError(f"Invalid MC era '{ivars.era}'. Available options: {list(globalTagDict['MC'].keys())}")
if not ivars.runOnMC and ivars.era not in globalTagDict['data']:
    raise ValueError(f"Invalid data era '{ivars.era}'. Available options: {list(globalTagDict['data'].keys())}")

# --- Core CMSSW process setup ---
process = cms.Process('MUONMATCH')

# --- Standard services and detector conditions ---
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

# Select the GlobalTag from the validated era table instead of hardcoding it in the analyzer cfg.
if ivars.runOnMC:
    myGlobalTag = globalTagDict['MC'][ivars.era]
else:
    myGlobalTag = globalTagDict['data'][ivars.era]
process.GlobalTag = GlobalTag(process.GlobalTag, myGlobalTag, '')

# --- Input source and optional event-range filtering ---
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(ivars.maxEvents))
process.source = cms.Source(
    'PoolSource',
    skipEvents=cms.untracked.uint32(0),
    fileNames=cms.untracked.vstring(ivars.inputFiles),
)
if ivars.eventRange != '':
    process.source.eventRange = cms.untracked.VEventRange(ivars.eventRange)

# --- Matching-study analyzer configuration ---
process.muonPackedCandMatch = cms.EDAnalyzer(
    'MuonPackedCandMatchNtuplizer',

    # MiniAOD collections studied by the standalone validation analyzer.
    muons=cms.untracked.InputTag('slimmedMuons'),
    packedCandidates=cms.untracked.InputTag('packedPFCandidates'),
    primaryVertices=cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),

    # Runtime selection, compact-mode row retention, and matching thresholds mirrored in the technical docs.
    studyAllMuons=cms.untracked.bool(ivars.studyAllMuons),
    StoreDetailedRowsOnlyOnDisagreement=cms.untracked.bool(ivars.storeDetailedRowsOnlyOnDisagreement),
    LoserRowsPerMethod=cms.untracked.int32(ivars.loserRowsPerMethod),
    DisagreementRowMode=cms.untracked.string(ivars.disagreementRowMode),
    MuonSelection=cms.untracked.string('pt > 2.5 && abs(eta) < 2.4'),
    PVSelectionMode=cms.untracked.string('firstVertex'),
    VectorRelPThreshold=cms.untracked.double(0.01),
    MomentumChi2Threshold=cms.untracked.double(25.0),
    MomentumDzPvChi2Threshold=cms.untracked.double(25.0),
    MomentumDzAssocChi2Threshold=cms.untracked.double(25.0),
    RequireChargeMatch=cms.untracked.bool(True),
    StoreAllPrimaryVertices=cms.untracked.bool(False),
    StorePointerDiagnostics=cms.untracked.bool(True),
    DebugUnmatchedSoftMuons=cms.untracked.bool(ivars.debugUnmatchedSoftMuons),
)

# --- ROOT output ---
process.TFileService = cms.Service(
    'TFileService',
    fileName=cms.string(ivars.outputFile),
)

# --- Schedule ---
process.p = cms.Path(process.muonPackedCandMatch)
process.schedule = cms.Schedule(process.p)
