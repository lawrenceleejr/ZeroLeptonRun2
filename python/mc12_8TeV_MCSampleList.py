__doc__ = """List of samples to be used in the analysis"""

##
## Baseline samples
##
lZjets   = range(167749,167760+1)+range(167797,167844+1) # SherpaMassiveCB
lYjets   = range(177574,177585+1)                        # SherpaMassiveCB
lWjets   = range(167740,167748+1)+range(167761,167796+1) # SherpaMassiveCB

lttbar   = [117050] # PowhegPythia_P2011C_ttbar
lsinglet = [108343,108344,108345,108346, # McAtNLO+FHerwig for s-chan and Wt-chan
            117360,117361,117362]        # AcerMC+P6 for t-chan
lttbarX  = [119353,119354,119355,119356,119583] # MadGraph5+P6
lTop     = lttbar + lsinglet + lttbarX

# SherpaMassiveCB diboson samples
lDiBoson = [177997,183734,183736,183738,               # WW
            179974,179975,183735,183737,183739,        # WZ
            183585,183587,183589,183591,               # ZW
            126894,177999,183586,183588,183590,183592] # ZZ
lZgamma  = [126854,145161,145162,146828,164438] # Massless CB
lWgamma  = [126739,126742,126856]               # Massless CB
lVgamma  = lZgamma + lWgamma
lTriBoson = [167006,167007,167008]
lDiBoson = lDiBoson + lVgamma + lTriBoson

lQCDMC   = range(147900,147907+1) # Pythia8 dijet, no weight

lbaseline = lZjets + lWjets + lYjets + lTop + lQCDMC + lDiBoson

##
## Alternative MC estimate samples
##
lZlljetsAlpgen = range(107650,107655+1)+range(107660,107665+1)+range(107670,107675+1) # Alpgen (no HFOR)
lZvvjetsAlpgen = range(156803,156828+1) # Alpgen with ptZ(truth) slices
lZjetsAlpgen   = lZlljetsAlpgen + lZvvjetsAlpgen
lYjetsAlpgen   = range(156839,156863+1) # Alpgen with ptgamma(truth) slices

# Sherpa with massless C and B quarks
lZlljetsSherpa = [147770,147771,147772] # Sherpa NLO
lZvvjetsSherpa = [157537,157538,157539,157540] # Sherpa with ptZ(truth) slices
lZDYSherpa     = range(173041,173046+1) # Sherpa DY samples
lZjetsSherpa   = lZlljetsSherpa + lZvvjetsSherpa + lZDYSherpa
lYjetsSherpa   = [113715,113716,113717,126371,126955,126956] # Sherpa with ptgamma(truth) slices
# Massive B quarks samples
lWjetsMassiveB  = range(167150,167161+1)  # Sherpa MassiveB samples

# As long as no b-tagging is used, the following samples are fine
lWjetsMasslessB = [147774,147775,147776, # Sherpa Inclusive samples
                   144992,144993,144994, # >= 3 parton-level jets
                   157534,157535,157536] # ptW(truth) > 200 GeV

# Be sure that top_hfor_d3pd is recomputed with HforToolD3PD::BBONLY option
lWLjetsAlpgen  = range(107680,107685+1)+range(107690,107695+1)+range(107700,107705+1)+range(172001,172006+1)+range(172011,172016+1)+range(172021,172026+1) # W+light jets
lWbbjetsAlpgen = range(107280,107283+1) # Wbb
lWjetsAlpgen   = lWLjetsAlpgen + lWbbjetsAlpgen
# Be sure that top_hfor_d3pd is computed properly. For Pythia overlap is removed wrt both b- and c- quarks.
lWLjetsAlpgenPythia = range(117680,117685+1)+range(117690,117695+1)+range(117700,117705+1) # W+light jets
lWHjetsAlpgenPythia = range(110801,110804+1)+range(126601,126605+1)+range(126606,126609+1) # Wbb+Wc+Wcc
lWjetsAlpgenPythia  = lWLjetsAlpgenPythia + lWHjetsAlpgenPythia
                 
lttbarMcAtNlo  = [105200] # McAtNlo+Herwig+Jimmy
lTopMcAtNlo    = lttbarMcAtNlo + lsinglet + lttbarX
lttbarSherpa   = range(117800,117809+1)
lTopSherpa     = lttbarSherpa + lsinglet + lttbarX
lttbarAlpgen   = range(164440,164443+1)+range(164450,164453+1)
lTopAlpgen     = lttbarAlpgen + lsinglet + lttbarX 
lttbarPowhegP6 = [105861] # PowhegPythia_AUET2BCT10_ttbar
lTopPowhegP6   = lttbarPowhegP6 + lsinglet + lttbarX 
lsingletMorePS = [117213,117215,117217,117219,117221,117223,117245]
lsingletLessPS = [117214,117216,117218,117220,117222,117224,117246]

# Sherpa diboson with massless C,B
lDiBosonSherpa = [143065,157814,157815,157816,157817,157818,157819,
                  126892,126893,126894,126895]
lDiBosonSherpa = lDiBosonSherpa + lVgamma + lTriBoson
lWWAlpgen = [110829,110830,110831,110832]
lWgammaAlgpen = range(146430,146435+1)

lalt = lZjetsAlpgen + lYjetsAlpgen + lWjetsAlpgen + lTopMcAtNlo
