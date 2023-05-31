import ROOT
import sys
from helper_functions import *
import math

################################################################################################################
# Simple script demonstrating how to apply acceptance selection and parameterized efficiencies
# found in: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-13/hepdata_info.pdf
# Authors: Emily Anne Thompson (LBNL), Catherine Xu (LBNL)
# Contact: emily.anne.thompson@cern.ch
# Date: 05/30/2023
#
# Input to this script are ROOT ntuples with the following branches:
#
# mcEventWeight: event weight from monte-carlo generator
# 
# (see document for more information on formation of truth jets and matching to LLP decays)
# truthJet_Pt: Truth jet pT in GeV
# truthJet_Eta: Truth jet eta
# truthJet_hasMatchedLLP: Boolean stored as true if truth jet has been matched to LLP decay
# truthJet_matchedLLP_x/y/z: Position of LLP decay which truth jet has been matched to
# truthJet_matchedLLP_pdgId: Pdg ID () of the LLP which the truth jet has been matched to
#
# (truthSparticle refers to a real SUSY particle in the truth record)
# truthSparticle_PVCut: Boolean stored as true if the location of the LLP decay is > 4mm in x-y direction from all PVs
# truthSparticle_VtxX/Y/Z: Position of LLP decay in x-y-z directions
# truthSparticle_Barcode: Barcode of LLP in MC truth record (used to identify its decay products)
#
# (truth_LLPChild refers to particle in truth record which originates from LLP decay)
# truth_LLPChild_Pt: Pt of LLP child in GeV
# truth_LLPChild_Eta: Eta of LLP child
# truth_LLPChild_Phi: Phi of LLP child
# truth_LLPChild_Parent_barcode: Barcode of the parent particle which this child originates from
# truth_LLPChild_Parent_pdgID: Pdg ID of the parent particle which this child originates from
# truth_LLPChild_Charge: Charge of LLP child
# truth_LLPChild_Status: Status of LLP child in MC truth record (status=1 means particle is stable in the generator we test with)
# truth_LLPChild_d0: Approximate d0 of LLP child, see documentation for how this is estimated
#
################################################################################################################


# Empty histograms to implement jet selection
h_highpt_filters  = ROOT.TH1F("HighptFilters", "HighptFilters; Jet Multiplicity; Entries;", 10, 0, 10);
h_lowpt_filters   = ROOT.TH1F("LowptFilters", "LowptFilters; Jet Multiplicity; Entries;", 10, 0, 10);
h_trackless_filters   = ROOT.TH1F("TracklessFilters", "TracklessFilters; Jet Multiplicity; Entries;", 4, 0, 4);

# Efficiencies from HEPdata: https://www.hepdata.net/record/137762
f_eventeff_highpt_R1150 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_mm.root", "READ")
h_eventeff_highpt_R1150 = f_eventeff_highpt_R1150.Get("event_efficiency_HighPt_R_1150_mm/Hist1D_y1")
f_eventeff_trackless_R1150 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_Trackless_R_1150_mm.root", "READ")
h_eventeff_trackless_R1150 = f_eventeff_trackless_R1150.Get("event_efficiency_Trackless_R_1150_mm/Hist1D_y1")
f_eventeff_highpt_R1150_R3870 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_HighPt_R_1150_3870_mm.root", "READ")
h_eventeff_highpt_R1150_R3870 = f_eventeff_highpt_R1150_R3870.Get("event_efficiency_HighPt_R_1150_3870_mm/Hist1D_y1")
f_eventeff_trackless_R1150_R3870 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_Trackless_R_1150_3870_mm.root", "READ")
h_eventeff_trackless_R1150_R3870 = f_eventeff_trackless_R1150_R3870.Get("event_efficiency_Trackless_R_1150_3870_mm/Hist1D_y1")
f_eventeff_highpt_R3870 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_HighPt_R_3870_mm.root", "READ")
h_eventeff_highpt_R3870 = f_eventeff_highpt_R3870.Get("event_efficiency_HighPt_R_3870_mm/Hist1D_y1")
f_eventeff_trackless_R3870 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-event_efficiency_Trackless_R_3870_mm.root", "READ")
h_eventeff_trackless_R3870 = f_eventeff_trackless_R3870.Get("event_efficiency_Trackless_R_3870_mm/Hist1D_y1")

f_vertexeff_R22 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_22_mm.root", "READ")
h_vertexeff_R22 = f_vertexeff_R22.Get("vertex_efficiency_R_22_mm/Hist2D_y1")
f_vertexeff_R22_R25 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_22_25_mm.root", "READ")
h_vertexeff_R22_R25 = f_vertexeff_R22_R25.Get("vertex_efficiency_R_22_25_mm/Hist2D_y1")
f_vertexeff_R25_R29 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_25_29_mm.root", "READ")
h_vertexeff_R25_R29 = f_vertexeff_R25_R29.Get("vertex_efficiency_R_25_29_mm/Hist2D_y1")
f_vertexeff_R29_R38 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_29_38_mm.root", "READ")
h_vertexeff_R29_R38 = f_vertexeff_R29_R38.Get("vertex_efficiency_R_29_38_mm/Hist2D_y1")
f_vertexeff_R38_R46 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_38_46_mm.root", "READ")
h_vertexeff_R38_R46 = f_vertexeff_R38_R46.Get("vertex_efficiency_R_38_46_mm/Hist2D_y1")
f_vertexeff_R46_R73 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_46_73_mm.root", "READ")
h_vertexeff_R46_R73 = f_vertexeff_R46_R73.Get("vertex_efficiency_R_46_73_mm/Hist2D_y1")
f_vertexeff_R73_R84 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_73_84_mm.root", "READ")
h_vertexeff_R73_R84 = f_vertexeff_R73_R84.Get("vertex_efficiency_R_73_84_mm/Hist2D_y1")
f_vertexeff_R84_R111 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_84_111_mm.root", "READ")
h_vertexeff_R84_R111 = f_vertexeff_R84_R111.Get("vertex_efficiency_R_84_111_mm/Hist2D_y1")
f_vertexeff_R111_R120 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_111_120_mm.root", "READ")
h_vertexeff_R111_R120 = f_vertexeff_R111_R120.Get("vertex_efficiency_R_111_120_mm/Hist2D_y1")
f_vertexeff_R120_R145 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_120_145_mm.root", "READ")
h_vertexeff_R120_R145 = f_vertexeff_R120_R145.Get("vertex_efficiency_R_120_145_mm/Hist2D_y1")
f_vertexeff_R145_R180 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_145_180_mm.root", "READ")
h_vertexeff_R145_R180 = f_vertexeff_R145_R180.Get("vertex_efficiency_R_145_180_mm/Hist2D_y1")
f_vertexeff_R180_R300 = ROOT.TFile("efficiencies/HEPData-ins2628398-v1-vertex_efficiency_R_180_300_mm.root", "READ")
h_vertexeff_R180_R300 = f_vertexeff_R180_R300.Get("vertex_efficiency_R_180_300_mm/Hist2D_y1")

# Function to apply truth jet acceptance selection
def applyJetSelection(truthJet_Pt, truthJet_Eta, truthJet_matchedLLP_x, truthJet_matchedLLP_y, truthJet_matchedLLP_z, truthJet_matchedLLP_pdgId, truthJet_hasMatchedLLP):

  h_highpt_filters.Reset()
  h_lowpt_filters.Reset()
  h_trackless_filters.Reset()

  # Pt jet thresholds for each jet multiplicity
  high_pt_threshold = {4: 250,
                       5: 195,
                       6: 116,
                       7: 90}
  low_pt_threshold = {4: 137,
                      5: 101,
                      6: 83,
                      7: 55}
  trackless_pt_threshold = {1: 70,
                            2: 50}

  truthJets_sumpt = 0
  for i in range(0, len(truthJet_Pt)):
  
    # Calculate sum of pT
    truthJets_sumpt += truthJet_Pt.at(i)

    # Don't use jets that have decayed outside of the calorimeter
    if (truthJet_hasMatchedLLP.at(i) == 1):

      # Make sure we are looking at correct particle
      if (abs(truthJet_matchedLLP_pdgId.at(i))<1000022 and abs(truthJet_matchedLLP_pdgId.at(i))>1000024): continue

      # Apply selection on R 
      truthJet_matchedLLP_R = math.sqrt(pow(truthJet_matchedLLP_x.at(i), 2.) + pow(truthJet_matchedLLP_y.at(i), 2.) )
      if( truthJet_matchedLLP_R > 3870. ): continue;  

    # Count jets above thresholds
    for j in range(4, 8):
      if (truthJet_Pt.at(i) > high_pt_threshold.get(j)): 
        h_highpt_filters.AddBinContent(j)
      if (truthJet_Pt.at(i) > low_pt_threshold.get(j)): h_lowpt_filters.AddBinContent(j);
       
    # Count displaced jets
    if (abs(truthJet_matchedLLP_pdgId.at(i)) >= 1000022 and abs(truthJet_matchedLLP_pdgId.at(i)) <= 1000024 ):
      if (truthJet_Pt.at(i) >  trackless_pt_threshold.get(1) and abs(truthJet_Eta.at(i))<2.5 ): h_trackless_filters.AddBinContent(1);
      if (truthJet_Pt.at(i) >  trackless_pt_threshold.get(2) and abs(truthJet_Eta.at(i))<2.5 ): h_trackless_filters.AddBinContent(2);


  passHighptSR = False
  passTracklessSR = False
  passLowptJet = False
  passTracklessJet = False
  for j in range(4, 8):
    if (h_highpt_filters.GetBinContent(j) >= j): passHighptSR = True
    if (h_lowpt_filters.GetBinContent(j) >= j): passLowptJet = True

  if (h_trackless_filters.GetBinContent(1) >= 1): passTracklessJet = True
  if (h_trackless_filters.GetBinContent(2) >= 2): passTracklessJet = True
  passTracklessSR= passLowptJet and passTracklessJet

  return passHighptSR, passTracklessSR, truthJets_sumpt

# Function to find maximum R_decay of all LLP decays in event
def getRmax(truthSparticle_VtxX, truthSparticle_VtxY, truthSparticle_VtxZ):
  Rmax = 0.
  for j in range(0, truthSparticle_VtxX.size()):
    R = math.hypot(truthSparticle_VtxX.at(j), truthSparticle_VtxY.at(j))
    if (R>Rmax): Rmax = R
  return Rmax

if len(sys.argv)<2:
  print("Please specify dsid!")
  exit(0)

dsid = sys.argv[1]
if (len(sys.argv)>2):
  tag = sys.argv[2]
else: tag = "v1"

# Load atlas style so plots look nice
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasStyle.C")
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasUtils.C")
ROOT.gROOT.LoadMacro("atlasstyle-00-04-02/AtlasLabels.C")
ROOT.SetAtlasStyle()

# Gather inputs and cross-section
inputs = get_files(dsid)
info = get_sample_info(dsid)
print("Analyzing DSID ", dsid, " with sample info = ", info)
xsec = get_cross_section(info)

# Set counters
allevents = 0
passHighPtJetAcc = 0
passHighPtAcc = 0
passTracklessJetAcc = 0
passTracklessAcc = 0
passHighPtEff = 0
passTracklessEff = 0

# Load ROOT TTree and loop through events in ntuple
for f in inputs:
   inFile = ROOT.TFile(f, 'READ')

   # Load TTree
   t = inFile.Get("trees_SRDV_")
   nentries = t.GetEntries()

   print("Total number of entries in file: ", nentries)
   # Get weight information
   metadata = inFile.Get("MetaData_EventCount")
   total_sum_weights = metadata.GetBinContent(3)

   # Luminosity (dummy value - doesn't matter for efficiency)
   lumi = 1.0

   # Calculate event weight
   weight_nomc = lumi * xsec / total_sum_weights

   # Loop through entries in TTree
   for i in range(0, nentries):
      t.GetEntry(i)

      # Set weight
      weight = weight_nomc * t.mcEventWeight 
      allevents += weight

      # Apply event-level acceptance
      passHighPtJet, passTracklessJet, truthJets_sumpt = applyJetSelection(t.truthJet_Pt, t.truthJet_Eta, t.truthJet_matchedLLP_x, t.truthJet_matchedLLP_y, t.truthJet_matchedLLP_z, t.truthJet_matchedLLP_pdgId, t.truthJet_hasMatchedLLP)

      if (passHighPtJet): passHighPtJetAcc += weight
      if (passTracklessJet): passTracklessJetAcc += weight

      # Skip events failing acceptance to save time
      if (not passHighPtJet and not passTracklessJet): continue

      # Apply vertex-level acceptance
      passes_vertex_acceptance = False
      vertex_efficiency = 1.0

      # Loop through all LLPs in event
      for j in range(0, t.truthSparticle_PVCut.size()):

        # Determine if decay happened in fiducial volume
        truthSparticle_fiducial = False
        R = math.hypot(t.truthSparticle_VtxX.at(j), t.truthSparticle_VtxY.at(j))
        if (R<300. and abs(t.truthSparticle_VtxZ.at(j))<300.): truthSparticle_fiducial = True

   
        p_tracks = ROOT.TLorentzVector(0., 0., 0., 0.);
        n_tracks = 0.
        flag_pass_truth_LLPChild_d0 = False 
        # Loop through LLP decay products
        for k in range(0, t.truth_LLPChild_Pt.size()):

          # Skip decay products from different LLP
          if (t.truthSparticle_Barcode.at(j) != t.truth_LLPChild_Parent_barcode.at(k)): continue
     
          # Store four-vector 
          tmp_track = ROOT.TLorentzVector();
          tmp_track.SetPtEtaPhiM(t.truth_LLPChild_Pt.at(k), t.truth_LLPChild_Eta.at(k), t.truth_LLPChild_Phi.at(k), 0.135) # Assumes charged pion mass
          
          # Count track if is passes selection
          if (abs(t.truth_LLPChild_Charge.at(k))>0. and t.truth_LLPChild_Pt.at(k)>1. and t.truth_LLPChild_Status.at(k)>0. and abs(t.truth_LLPChild_Parent_pdgID.at(k))>=1000022 and abs(t.truth_LLPChild_Parent_pdgID.at(k))<=1000024):
            p_tracks+=tmp_track
            n_tracks+=1

          # Determine if LLP decay has at least 1 charged "track" with d0>2 mm
          if (abs(t.truth_LLPChild_Charge.at(k))>0. and abs(t.truth_LLPChild_d0.at(k))>2.): flag_pass_truth_LLPChild_d0 = True

        # Apply vertex acceptance and find efficiency
        if (truthSparticle_fiducial and t.truthSparticle_PVCut.at(j)==1 and n_tracks>=5 and p_tracks.M()>10. and flag_pass_truth_LLPChild_d0):
          passes_vertex_acceptance = True
          if (R<22.): tmp_weight = h_vertexeff_R22.GetBinContent(h_vertexeff_R22.FindBin(p_tracks.M(), n_tracks))
          elif (R<25.): tmp_weight = h_vertexeff_R22_R25.GetBinContent(h_vertexeff_R22_R25.FindBin(p_tracks.M(), n_tracks))
          elif (R<29.): tmp_weight = h_vertexeff_R25_R29.GetBinContent(h_vertexeff_R25_R29.FindBin(p_tracks.M(), n_tracks))
          elif (R<38.): tmp_weight = h_vertexeff_R29_R38.GetBinContent(h_vertexeff_R29_R38.FindBin(p_tracks.M(), n_tracks))
          elif (R<46.): tmp_weight = h_vertexeff_R38_R46.GetBinContent(h_vertexeff_R38_R46.FindBin(p_tracks.M(), n_tracks))
          elif (R<73.): tmp_weight = h_vertexeff_R46_R73.GetBinContent(h_vertexeff_R46_R73.FindBin(p_tracks.M(), n_tracks))
          elif (R<84.): tmp_weight = h_vertexeff_R73_R84.GetBinContent(h_vertexeff_R73_R84.FindBin(p_tracks.M(), n_tracks))
          elif (R<111.): tmp_weight = h_vertexeff_R84_R111.GetBinContent(h_vertexeff_R84_R111.FindBin(p_tracks.M(), n_tracks))
          elif (R<120.): tmp_weight = h_vertexeff_R111_R120.GetBinContent(h_vertexeff_R111_R120.FindBin(p_tracks.M(), n_tracks))
          elif (R<145.): tmp_weight = h_vertexeff_R120_R145.GetBinContent(h_vertexeff_R120_R145.FindBin(p_tracks.M(), n_tracks))
          elif (R<180.): tmp_weight = h_vertexeff_R145_R180.GetBinContent(h_vertexeff_R145_R180.FindBin(p_tracks.M(), n_tracks))
          elif (R<300.): tmp_weight = h_vertexeff_R180_R300.GetBinContent(h_vertexeff_R180_R300.FindBin(p_tracks.M(), n_tracks))
          else: print("Warning! Vertex that passes acceptance with R = ", R, " ?? ")

          vertex_efficiency = vertex_efficiency * (1. - tmp_weight )

      vertex_efficiency = 1. - vertex_efficiency     

      if (passHighPtJet and passes_vertex_acceptance): passHighPtAcc += weight
      if (passTracklessJet and passes_vertex_acceptance): passTracklessAcc += weight

      # Prepare to calculate event-level efficiency
      event_efficiency_highpt = 0.
      event_efficiency_trackless = 0.
      Rmax = getRmax(t.truthSparticle_VtxX, t.truthSparticle_VtxY, t.truthSparticle_VtxZ)

      # Converte truthJets_sumpt from GeV to TeV
      truthJets_sumpt = truthJets_sumpt/1000.      

      # Find event-level efficiency
      if (passHighPtJet and Rmax<1150):
        event_efficiency_highpt = h_eventeff_highpt_R1150.GetBinContent(h_eventeff_highpt_R1150.FindBin(truthJets_sumpt))
      elif (passHighPtJet and Rmax < 3870):
        event_efficiency_highpt = h_eventeff_highpt_R1150_R3870.GetBinContent(h_eventeff_highpt_R1150_R3870.FindBin(truthJets_sumpt))
      elif (passHighPtJet):
        event_efficiency_highpt = h_eventeff_highpt_R3870.GetBinContent(h_eventeff_highpt_R3870.FindBin(truthJets_sumpt))

      if (passTracklessJet and Rmax<1150):
        event_efficiency_trackless = h_eventeff_trackless_R1150.GetBinContent(h_eventeff_trackless_R1150.FindBin(truthJets_sumpt))
      elif (passTracklessJet and Rmax < 3870):
        event_efficiency_trackless = h_eventeff_trackless_R1150_R3870.GetBinContent(h_eventeff_trackless_R1150_R3870.FindBin(truthJets_sumpt))
      elif (passTracklessJet):
        event_efficiency_trackless = h_eventeff_trackless_R3870.GetBinContent(h_eventeff_trackless_R3870.FindBin(truthJets_sumpt))

      passHighPtEff += weight * (event_efficiency_highpt * vertex_efficiency)
      passTracklessEff += weight * (event_efficiency_trackless * vertex_efficiency)

print("High-pt jet acceptance: ", passHighPtJetAcc/allevents*100., " %")
print("High-pt full acceptance: ", passHighPtAcc/allevents*100., " %")
print("High-pt acc x eff: ", passHighPtEff/allevents*100., " %")
print("Trackless jet acceptance: ", passTracklessJetAcc/allevents*100., " %")
print("Trackless full acceptance: ", passTracklessAcc/allevents*100., " %")
print("Trackless acc x eff: ", passTracklessEff/allevents*100., " %") 

