__doc__ = """List of samples to be used in the analysis"""
import makeSignalPointPickle

##
## Baseline samples   built from 50ns tags [r6630, r6647, r6655, r6767, r6793, r6802, r6828]
##
#lZMassiveCB   = range(361444,361467+1) # Sherpa Znunu
#lZMassiveCB  += range(361372,361395+1) # Sherpa Zee
#lZMassiveCB  += range(361396,361419+1) # Sherpa Zmumu
#lZMassiveCB  += range(361420,361443+1) # Sherpa Ztautau

ZMassiveCB   = range(363412,363435+1) # Sherpa 2.2 Znunu
ZMassiveCB  += range(363391,363411+1) # Sherpa 2.2 Zee
ZMassiveCB  += range(363364,363387+1) # Sherpa 2.2 Zmumu
ZMassiveCB  += range(363102,363122+1) # + range(363361,363363+1) # Sherpa 2.2 Ztautau


#ZMassiveCBLO = range(407100,407183+1) # L0 Z+jets

GammaMassiveCB   = range(361039,361061+1) # Sherpa

#Wjets   = range(361300,361323+1) # Sherpa Wenu missing pt<70
#Wjets  += range(361324,361347+1) # Sherpa Wmunu  missing pt<70
#Wjets  += range(361348,361371+1) # Sherpa Wtaunu

Wjets   = range(363463,363483+1) # Sherpa 2.2 Wenu missing pt<70
Wjets  += range(363439,363459+1) # Sherpa 2.2 Wmunu  missing pt<70
Wjets  += range(363334,363354+1) # Sherpa 2.2 Wtaunu


ttbar   = [410000] # PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonahad
singlet = [410011, 410012, # PowhegPythia t-channel
            410013, 410014, # PowhegPythia Wt
            410025, 410026, # PowhegPythia s-channel
            ]
ttbarX   = range(410066,410068+1)  # MadGraph5+Pythia ttW
ttbarX  += range(410069,410075+1)  # MadGraph5+Pythia ttZ
ttbarX  += [410081] # MadGraph5+Pythia ttWW

Top     = ttbar + singlet + ttbarX

# diboson sampes
#DibosonMassiveCB =  range(361063,361073+1)+[361077,]+range(361081,361089+1)
DibosonMassiveCB =  range(361063,361073+1)+[361077,]+range(361087,361088+1)+range(361091,361096+1) # use of 361087 is temporary to be repaced with 361097

TriBoson = []

Zgamma  = [301908,301909] # Sherpa nunugamma 35<pt<140
Wgamma  = [301893,301894,301902] # Sherpa  several missing sampes
Vgamma  = Zgamma + Wgamma
DibosonMassiveCB = DibosonMassiveCB  + TriBoson + Vgamma

QCD   = range(361020,361032+1) # Pythia8 dijet JZxW

baseline = ZMassiveCB + Wjets + GammaMassiveCB + Top + QCD + DibosonMassiveCB

##
## Alternative MC samples
##

##[187150,187158+1]+[187160,187168+1]+[187170,187178+1]+[187180,187188+1]+[206618,206620+1] # diboson PowhegPythia

ZMassiveCBPowheg = [361106,361107,361108] # PowhegPythia8 sampes
ZMassiveCBMadgraph = range(361500,361519+1)  # Madgraph+Pythia8
ZtautaujetsSherpa22 = range(363102,363122+1) + range(363361,363363+1)
ZmumujetsSherpa22 = range(363364,363387+1)
ZeejetsSherpa22 = range(363388,363411+1)
ZnunujetsSherpa22 = range(363412,363435+1)
ZMassiveCBSherpa22 = ZtautaujetsSherpa22+ZmumujetsSherpa22+ZeejetsSherpa22+ZnunujetsSherpa22
ZMassiveCBAlt = ZMassiveCBPowheg + ZMassiveCBMadgraph + ZMassiveCBSherpa22

WjetsPowheg = [361100,361101,361102,361103,361104,361105]  # PowhegPythia8 sampes
WjetsMadgraph = range(361520,361534+1)  # Madgraph+Pythia
WjtaunuetsSherpa22 = range(363331,363354+1)
WmunujetsSherpa22 = range(363436,363459+1)
WenujetsSherpa22 = range(363460,363483+1)
WjetsSherpa22 = WjtaunuetsSherpa22+WmunujetsSherpa22+WenujetsSherpa22
WjetsAlt = WjetsPowheg + WjetsMadgraph + WjetsSherpa22

ttbarPowhegSyst = [410001,410002]
ttbarMCatNLO = [410003]
ttbarPowhegHerwing = [410004]
ttbarPowhegP8 = [410006]
ttbarSherpa = [410021,410022,410023]
ttbarAlt = ttbarPowhegSyst + ttbarMCatNLO + ttbarPowhegHerwing + ttbarPowhegP8 + ttbarSherpa
singletAlt = [ ]

DibosonMassiveCBPowheg =  range(361600,361611+1)
DibosonMassiveCBAlt = DibosonMassiveCBPowheg


alt = ZMassiveCBAlt + WjetsAlt + ttbarAlt + singletAlt + DibosonMassiveCBAlt

SS_direct = list(makeSignalPointPickle.pointdict["SS_direct"].keys())
GG_direct = list(makeSignalPointPickle.pointdict["GG_direct"].keys())

##
## truth level samples
##

# ZMassiveCBTruth  = range(362000,362190+1) # Sherpa Z->nunu syst variations
# ZMassiveCBTruth += range(361444,361467+1) # Sherpa Z->nunu
# ZMassiveCBTruth += range(361515,361519+1) # Madgraph+Pythia8 Z->nunu
# GammaMassiveCBTruth = range(361039,361061+1) # Sherpa gamma+jets

# truth = ZMassiveCBTruth + GammaMassiveCBTruth

