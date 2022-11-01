import matplotlib.pyplot as plt

import boolean2
import numpy as np
import pylab
import seaborn as sns
import pandas as pd
import seaborn as sb
from boolean2 import util
from util import create_heatmap





Initial_conditions='''
FANCM = False
FAcore = False
FANCD2I = False
MUS81 = False
FANCJBRCA1 = False
XPF = False
FAN1 = False
PCNATLS = False
MRN = False
BRCA1 = False
ssDNARPA = False
FANCD1N = False
RAD51 = False
HRR = False
USP1 = False
KU = False
DNAPK = False
NHEJ = False
ATR = False
ATM = False
p53 = False
CHK1 = False
CHK2 = False
H2AX = False
CHKREC = False
'''



Rules='''
ICL* = ICL and not DSB
FANCM* = ICL and not CHKREC
FAcore* = FANCM and (ATR or ATM) and not CHKREC
FANCD2I* = FAcore and ((ATM or ATR) or (H2AX and DSB)) and not USP1
MUS81* = ICL
FANCJBRCA1* = (ICL or ssDNARPA) and (ATM or ATR)
XPF* = (MUS81 and not FANCM) or (MUS81 and p53 and not (FAcore and FANCD2I and FAN1))
FAN1* = MUS81 and FANCD2I
ADD* = (ADD or (MUS81 and (FAN1 or XPF))) and not PCNATLS
DSB* = (DSB or FAN1 or XPF) and not (NHEJ or HRR)
PCNATLS* = (ADD or (ADD and FAcore)) and not (USP1 or FAN1)
MRN* = DSB and ATM and not ((KU and FANCD2I) or RAD51 or CHKREC)
BRCA1* = DSB and (ATM or CHK2 or ATR) and not CHKREC
ssDNARPA* = DSB and ((FANCD2I and FANCJBRCA1) or MRN) and not (RAD51 or KU)
FANCD1N* = (ssDNARPA and BRCA1) or (FANCD2I and ssDNARPA) and not CHKREC
RAD51* = ssDNARPA and FANCD1N and not CHKREC
HRR* = DSB and RAD51 and FANCD1N and BRCA1 and not CHKREC
USP1* = ((FANCD1N and FANCD2I) or PCNATLS) and not FANCM
KU* = DSB and not (MRN or FANCD2I or CHKREC)
DNAPK* = (DSB and KU) and not CHKREC
NHEJ* = (DSB and DNAPK and XPF and not ((FANCJBRCA1 and ssDNARPA) or CHKREC)) or ((DSB and DNAPK and KU) and not (ATM and ATR))
ATR* = (ssDNARPA or FANCM or ATM) and not CHKREC
ATM* = (ATR or DSB) and not CHKREC
p53* = (((ATM and CHK2) or (ATR and CHK1)) or DNAPK) and not CHKREC
CHK1* = (ATM or ATR or DNAPK) and not CHKREC
CHK2* = (ATM or ATR or DNAPK) and not CHKREC
H2AX* = DSB and (ATM or ATR or DNAPK) and not CHKREC
CHKREC* = ((PCNATLS or NHEJ or HRR) and not DSB) or ((not ADD) and (not ICL) and (not DSB) and not CHKREC)
'''




for i in range(3):
    if i==0:
        initial_model = Initial_conditions + '''DSB=True\nADD=False\nICL=False''' + Rules
    elif i==1:
        initial_model = Initial_conditions + '''DSB=False\nADD=True\nICL=False''' + Rules
    else:
        initial_model = Initial_conditions + '''DSB=False\nADD=False\nICL=True''' + Rules





    time_iteration=20
    #coll = util.Collector()
    #step one


    model= boolean2.Model(initial_model,mode='sync')

    model.initialize()


    model.iterate(steps=time_iteration)
    #coll.collect(states=model.states, nodes=model.nodes)
    states=[]
    for state in model.states:
       states.append([int(x) for x in state.values()])
    nodes=state.keys()
    print (type(nodes))

    create_heatmap(states,nodes)





