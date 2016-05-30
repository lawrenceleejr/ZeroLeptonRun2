__doc__ = """List of samples to be used in the analysis"""

##
## Baseline samples   built from 50ns tags [r6630, r6647, r6655, r6767, r6793, r6802, r6828]
##
lZjets   = range(361444,361467+1) # Sherpa Znunu 
lZjets  += range(361372,361395+1) # Sherpa Zee
lZjets  += range(361396,361419+1) # Sherpa Zmumu
lZjets  += range(361420,361443+1) # Sherpa Ztautau

#lZjetsLO = range(407100,407183+1) # L0 Z+jets

lYjets   = range(361039,361061+1) # Sherpa

lWjets   = range(361300,361323+1) # Sherpa Wenu missing pt<70
lWjets  += range(361324,361347+1) # Sherpa Wmunu  missing pt<70 
lWjets  += range(361348,361371+1) # Sherpa Wtaunu


lttbar   = [410000] # PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
lsinglet = [410011, 410012, # PowhegPythia t-channel
            410013, 410014, # PowhegPythia Wt
            410025, 410026, # PowhegPythia s-channel
            ]       
lttbarX   = range(410066,410068+1)  # MadGraph5+Pythia ttW
lttbarX  += range(410069,410075+1)  # MadGraph5+Pythia ttZ
lttbarX  += [410081] # MadGraph5+Pythia ttWW

lTop     = lttbar + lsinglet + lttbarX

# diboson samples
lDiBoson =  range(361063,361073+1)+[361077,]+range(361081,361089+1) 

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

lZjetsPowheg = [361106,361107,361108] # PowhegPythia8 samples
lZjetsMadgraph = range(361500,361519+1)  # Madgraph+Pythia8
lZtautaujetsSherpa22 = range(363102,363122+1) + range(363361,363363+1)
lZmumujetsSherpa22 = range(363364,363387+1)
lZeejetsSherpa22 = range(363388,363411+1)
lZnunujetsSherpa22 = range(363412,363435+1)
lZjetsSherpa22 = lZtautaujetsSherpa22+lZmumujetsSherpa22+lZeejetsSherpa22+lZnunujetsSherpa22
lZjetsAlt = lZjetsPowheg + lZjetsMadgraph + lZjetsSherpa22

lWjetsPowheg = [361100,361101,361102,361103,361104,361105]  # PowhegPythia8 samples
lWjetsMadgraph = range(361520,361534+1)  # Madgraph+Pythia
lWjtaunuetsSherpa22 = range(363331,363354+1) 
lWmunujetsSherpa22 = range(363436,363459+1)
lWenujetsSherpa22 = range(363460,363483+1)
lWjetsSherpa22 = lWjtaunuetsSherpa22+lWmunujetsSherpa22+lWenujetsSherpa22
lWjetsAlt = lWjetsPowheg + lWjetsMadgraph + lWjetsSherpa22

lttbarPowhegSyst = [410001,410002]
lttbarMCatNLO = [410003]
lttbarPowhegHerwing = [410004]
lttbarPowhegP8 = [410006]
lttbarSherpa = [410021,410022,410023]
lttbarAlt = lttbarPowhegSyst + lttbarMCatNLO + lttbarPowhegHerwing + lttbarPowhegP8 + lttbarSherpa
lsingletAlt = [ ]

lDiBosonPowheg =  range(361600,361611+1) 
lDiBosonAlt = lDiBosonPowheg


lalt = lZjetsAlt + lWjetsAlt + lttbarAlt + lsingletAlt + lDiBosonAlt


##
## truth level samples
##

lZjetsTruth  = range(362000,362190+1) # Sherpa Z->nunu syst variations
lZjetsTruth += range(361444,361467+1) # Sherpa Z->nunu
lZjetsTruth += range(361515,361519+1) # Madgraph+Pythia8 Z->nunu
lYjetsTruth = range(361039,361061+1) # Sherpa gamma+jets

ltruth = lZjetsTruth + lYjetsTruth

