from util import *

# construct the model
Initial_conditions = '''
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

Rules = '''
1: ICL* = ICL and not DSB
1: FANCM* = ICL and not CHKREC
1: FAcore* = FANCM and (ATR or ATM) and not CHKREC
1: FANCD2I* = FAcore and ((ATM or ATR) or (H2AX and DSB)) and not USP1
1: MUS81* = ICL
1: FANCJBRCA1* = (ICL or ssDNARPA) and (ATM or ATR)
1: XPF* = (MUS81 and not FANCM) or (MUS81 and p53 and not (FAcore and FANCD2I and FAN1))
1: FAN1* = MUS81 and FANCD2I
1: ADD* = (ADD or (MUS81 and (FAN1 or XPF))) and not PCNATLS
1: DSB* = (DSB or FAN1 or XPF) and not (NHEJ or HRR)
1: PCNATLS* = (ADD or (ADD and FAcore)) and not (USP1 or FAN1)
1: MRN* = DSB and ATM and not ((KU and FANCD2I) or RAD51 or CHKREC)
1: BRCA1* = DSB and (ATM or CHK2 or ATR) and not CHKREC
1: ssDNARPA* = DSB and ((FANCD2I and FANCJBRCA1) or MRN) and not (RAD51 or KU)
1: FANCD1N* = (ssDNARPA and BRCA1) or (FANCD2I and ssDNARPA) and not CHKREC
1: RAD51* = ssDNARPA and FANCD1N and not CHKREC
1: HRR* = DSB and RAD51 and FANCD1N and BRCA1 and not CHKREC
1: USP1* = ((FANCD1N and FANCD2I) or PCNATLS) and not FANCM
1: KU* = DSB and not (MRN or FANCD2I or CHKREC)
1: DNAPK* = (DSB and KU) and not CHKREC
1: NHEJ* = (DSB and DNAPK and XPF and not ((FANCJBRCA1 and ssDNARPA) or CHKREC)) or ((DSB and DNAPK and KU) and not (ATM and ATR))
1: ATR* = (ssDNARPA or FANCM or ATM) and not CHKREC
1: ATM* = (ATR or DSB) and not CHKREC
1: p53* = (((ATM and CHK2) or (ATR and CHK1)) or DNAPK) and not CHKREC
1: CHK1* = (ATM or ATR or DNAPK) and not CHKREC
1: CHK2* = (ATM or ATR or DNAPK) and not CHKREC
1: H2AX* = DSB and (ATM or ATR or DNAPK) and not CHKREC
1: CHKREC* = ((PCNATLS or NHEJ or HRR) and not DSB) or (not ADD and not ICL and not DSB and not CHKREC)
'''
# CHKREC* = ((PCNATLS or NHEJ or HRR) and not DSB) or ((not ADD) and (not ICL) and (not DSB) and not CHKREC)

# DSB = double-strand break
# ADD = DNA adduct
# ICL = interstrand crosslink
for i in range(3):  # TODO: we should generalize this
    if i == 0:
        model_text = Initial_conditions + 'DSB=True\nADD=False\nICL=False\n' + Rules
    elif i == 1:
        model_text = Initial_conditions + 'DSB=False\nADD=True\nICL=False\n' + Rules
    else:
        model_text = Initial_conditions + 'DSB=False\nADD=False\nICL=True\n' + Rules

    # run a synchronous updating simulation of the BooleanNet model
    n_iterations = 20
    run_FA_Boolean_synch(model_text, n_iterations)

    # run an asynchronous updating simulation of the BooleanNet model
    n_iterations = 20
    n_runs = 100
    run_FA_Boolean_asynch(model_text, n_iterations, n_runs)

    plt.show()
