__doc__ = """List of samples to be used in the analysis"""

##
## Baseline samples
##
lZjets   = range(167797,167844+1)+ range(167749,167760+1) # SherpaMassiveCB
lYjets   = [] # SherpaMassiveCB
lWjets   = range(167740,167748+1)+range(167761,167796+1)+range(180534,180542+1) # SherpaMassiveCB

lttbar   = [110401] # PowhegPythia_P2012_ttbar_nonallhad
lsinglet = [110302,110305,110070,110071, # PowhegPythia
            ]       
lttbarX  = [119583,119353,119355,174830,174831,174832,174833] # MadGraph5+Pythia
lTop     = lttbar + lsinglet + lttbarX

# diboson samples   NB: WZ !!! dilepton filter !!!
lDiBoson = range(187150,187158+1) + range(187160,187178+1) + range(187180,187188+1)
lZgamma  = [] 
lWgamma  = []
lVgamma  = lZgamma + lWgamma
lTriBoson = []
lDiBoson = lDiBoson + lVgamma + lTriBoson

lQCDMC   = range(147910,147917+1) # Pythia8 dijet

lbaseline = lZjets + lWjets + lYjets + lTop + lQCDMC + lDiBoson

##
## Alternative MC estimate samples
##
lZlljetsAlpgen = []
lZvvjetsAlpgen = []
lZjetsAlpgen   = lZlljetsAlpgen + lZvvjetsAlpgen
lYjetsAlpgen   = []

# Sherpa with massless C and B quarks
lZlljetsSherpa = []
lZvvjetsSherpa = []
lZDYSherpa     = []
lZjetsSherpa   = lZlljetsSherpa + lZvvjetsSherpa + lZDYSherpa
lYjetsSherpa   = []
# Massive B quarks samples
lWjetsMassiveB  = []  # Sherpa MassiveB samples

# As long as no b-tagging is used, the following samples are fine
lWjetsMasslessB = []

# Be sure that top_hfor_d3pd is recomputed with HforToolD3PD::BBONLY option
lWLjetsAlpgen  = []
lWbbjetsAlpgen = []
lWjetsAlpgen   = lWLjetsAlpgen + lWbbjetsAlpgen
# Be sure that top_hfor_d3pd is computed properly. For Pythia overlap is removed wrt both b- and c- quarks.
lWLjetsAlpgenPythia = []
lWHjetsAlpgenPythia = []
lWjetsAlpgenPythia  = lWLjetsAlpgenPythia + lWHjetsAlpgenPythia
                 
lttbarMcAtNlo  = [] # McAtNlo+Herwig+Jimmy
lTopMcAtNlo    = lttbarMcAtNlo + lsinglet + lttbarX
lttbarSherpa   = []
lTopSherpa     = lttbarSherpa + lsinglet + lttbarX
lttbarAlpgen   = []
lTopAlpgen     = lttbarAlpgen + lsinglet + lttbarX 
lttbarPowhegP6 = [] # PowhegPythia_AUET2BCT10_ttbar
lTopPowhegP6   = lttbarPowhegP6 + lsinglet + lttbarX 
lsingletMorePS = []
lsingletLessPS = []

# Sherpa diboson with massless C,B
lDiBosonSherpa = []
lDiBosonSherpa = lDiBosonSherpa + lVgamma + lTriBoson
lWWAlpgen = []
lWgammaAlgpen = []

lalt = lZjetsAlpgen + lYjetsAlpgen + lWjetsAlpgen + lTopMcAtNlo
