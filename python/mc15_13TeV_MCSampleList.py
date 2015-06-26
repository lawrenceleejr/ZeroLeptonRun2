__doc__ = """List of samples to be used in the analysis"""

##
## Baseline samples
##
lZjets   = range(361372,361467+1) # Sherpa
lYjets   = range(361039,361056+1) # Sherpa
lWjets   = range(361300,361371+1) # Sherpa

lttbar   = [410000] # PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
lsinglet = [410011, 410012, # .PowhegPythia t-channel
            410013,410014,410015,410016, # PowhegPythia Wt
            ]       
lttbarX  = [] # MadGraph5+Pythia
lTop     = lttbar + lsinglet + lttbarX

# diboson samples
lDiBoson = [361084,361085]
lZgamma  = [] 
lWgamma  = []
lVgamma  = lZgamma + lWgamma
lTriBoson = []
lDiBoson = lDiBoson + lVgamma + lTriBoson

lQCDMC   = range(361020,361031+1) # Pythia8 dijet JZxW

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
