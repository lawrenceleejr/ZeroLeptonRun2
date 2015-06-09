__doc__ = """2015 week1 samples for first tests"""

##
## Baseline samples
##
lZjets   = range(361444,361449+1) # Sherpa Znunu
lZjets  +=  [361106,361107,361108] # Pythia8 Zll

lYjets   = range(361039,361047+1) # Sherpa
lWjets   = range(361100,361108+1) # PowHeg+Pythia8

lttbar   = [410000] # PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
lsinglet = []       
lttbarX  = []
lTop     = lttbar + lsinglet + lttbarX

# diboson samples
lDiBoson = [361084,361085]
lZgamma  = [] 
lWgamma  = []
lVgamma  = lZgamma + lWgamma
lTriBoson = []
lDiBoson = lDiBoson + lVgamma + lTriBoson

lQCDMC   = range(361020,361032+1) # Pythia8 dijet JZxW

lbaseline = lZjets + lWjets + lYjets + lTop + lQCDMC + lDiBoson

