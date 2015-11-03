__doc__ = """List of samples to be used in the analysis"""

##
## Baseline samples   built from 50ns tags [r6630, r6647, r6655, r6767, r6793, r6802, r6828]
##
lZjets   = range(361444,361467+1) # Sherpa Znunu 
lZjets  += range(361372,361395+1) # Sherpa Zee
lZjets  += range(361396,361419+1) # Sherpa Zmumu
lZjets  += range(361420,361443+1) # Sherpa Ztautau

lZjetsLO = range(407100,407183+1) # L0 Z+jets

lYjets   = range(361039,361061+1) # Sherpa

lWjets   = range(361300,361323+1) # Sherpa Wenu missing pt<70
lWjets  += range(361324,361347+1) # Sherpa Wmunu  missing pt<70 
lWjets  += range(361348,361371+1) # Sherpa Wtaunu


lttbar   = [410000] # PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
lsinglet = [410011, 410012, # PowhegPythia t-channel
            410013, 410014, 410015, 410016, # PowhegPythia Wt
            410025, 410026, # PowhegPythia s-channel
            ]       
lttbarX  = [] # MadGraph5+Pythia
lTop     = lttbar + lsinglet + lttbarX

# diboson samples
lDiBoson =  range(361063,361089+1)
lTriBoson = []

lZgamma  = [301908,301909] # Sherpa nunugamma 35<pt<140
lWgamma  = [301893,301894,301902] # Sherpa  several missing samples
lVgamma  = lZgamma + lWgamma
lDiBoson = lDiBoson  + lTriBoson + lVgamma

lQCDMC   = range(361020,361032+1) # Pythia8 dijet JZxW

lbaseline = lZjets + lWjets + lYjets + lTop + lQCDMC + lDiBoson

##
## Alternative MC samples
##

##[187150,187158+1]+[187160,187168+1]+[187170,187178+1]+[187180,187188+1]+[206618,206620+1] # diboson PowhegPythia

lZjetsAlt = [361106,361107,361108] # PowhegPythia8 samples
lZjetsAlt += range(361500,361519+1)  # Madgraph+Pythia8

lWjetsAlt = [361100,361101,361102,361103,361104,361105]  # PowhegPythia8 samples
lWjetsAlt += range(361520,361534+1)  # Madgraph+Pythia8

lttbarAlt = [410001,410002,410003,410004,410006]
lsingletAlt = [ 410015, 410016, # PowhegPythia Wt dilepton (overlap with 13+14
               ]       

lalt = lZjetsLO + lZjetsAlt + lWjetsAlt + lttbarAlt + lsingletAlt

#print lbaseline
