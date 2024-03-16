from util import *
from fanconi_pysb import model as model_pysb
import re

# construct the BooleanNet model
Initial_conditions = '''
ICL = False
FANCM = False
FAcore = False
FANCD2I = False
MUS81 = False
FANCJBRCA1 = False
XPF = False
FAN1 = False
ADD = False
DSB = False
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

CONDITIONS = [
    {'mutations': None, 'initials': [("DSB", True)]},  # FIG 2A
    {'mutations': None, 'initials': [("ADD", True)]},  # FIG 2B
    {'mutations': None, 'initials': [("ICL", True)]},  # FIG 2C
    {'mutations': [("ICL", True)], 'initials': None},  # FIG 3
    {'mutations': [("FAcore", False)], 'initials': [("ICL", True)]},  # FIG 4A
    {'mutations': [("FAcore", False), ("ICL", True)], 'initials': None},  # FIG 4B
    {'mutations': [("FANCD1N", False)], 'initials': [("ICL", True)]},  # FIG 5A
    {'mutations': [("FANCD1N", False), ("ICL", True)], 'initials': None}  # FIG 5B
]

for cond in CONDITIONS:

    print(cond)
    prefix = ""

    model_text = Initial_conditions + Rules  # for BooleanNet model
    param_values = {}  # for PySB model

    # loop over mutations
    if cond['mutations'] is not None:
        prefix = "MUT"
        for mut in cond['mutations']:
            prefix += "_%s_%s" % (mut[0], mut[1])

            # === BooleanNet model ===

            # set the initial state of the node
            m = re.search(r'%s\s*=\s*\w+' % mut[0], model_text)
            if m is not None:
                model_text = model_text.replace(m.group(0), "%s = %s" % (mut[0], mut[1]))
            else:
                print("Error: can't find node '%s' in the initial conditions." % mut[0])
                quit()

            # remove the update rule
            m = re.search(r'\s*\d+:\s*%s\*.*' % mut[0], model_text)
            if m is not None:
                model_text = model_text.replace(m.group(0), "")
            else:
                print("Error: can't find node '%s' in the rules." % mut[0])
                quit()

            # === PySB model ===

            # set the initial state of the node
            param_values['%s_%s_init' % (mut[0], mut[1])] = 1.0
            param_values['%s_%s_init' % (mut[0], not mut[1])] = 0.0

            # remove (deactivate) the update rule
            param_values['k_rate_%s' % mut[0]] = 0.0

    # loop over initial conditions
    if cond['initials'] is not None:
        prefix += "INIT" if prefix == "" else "_INIT"
        for init in cond['initials']:
            prefix += "_%s_%s" % (init[0], init[1])

            # === BooleanNet model ===

            # set the initial state of the node
            m = re.search(r'%s\s*=\s*\w+' % init[0], model_text)
            if m is not None:
                model_text = model_text.replace(m.group(0), "%s = %s" % (init[0], init[1]))
            else:
                print("Error: can't find node '%s' in initial conditions." % init[0])
                quit()

            # === PySB model ===

            # set the initial state of the node
            param_values['%s_%s_init' % (init[0], init[1])] = 1.0
            param_values['%s_%s_init' % (init[0], not init[1])] = 0.0

    # if no mutations or initial conditions, this is the base model
    if prefix == "":
        prefix = "base"

    # synchronous updating
    n_iterations = 30
    run_FA_Boolean_synch(model_text, n_iterations, outfile='%s_synch.pdf' % prefix)

    # asynchronous updating
    n_iterations = 30
    n_runs = 100
    run_FA_Boolean_asynch(model_text, n_iterations, n_runs, verbose=False, outfile='%s_asynch.pdf' % prefix)

    # PySB general asynchronous (Gillespie) simulations
    t_end = 30
    n_runs = 100
    run_FA_Boolean_pysb(model_pysb, t_end, n_runs, param_values=param_values, verbose=False,
                        outfile='%s_pysb.pdf' % prefix)

# Export model generated by Boolean2Rules so don't have to regenerate it each time (very time-consuming)
# from pysb.export import export
# pysb_model = export(model, 'pysb_flat')
# with open('fanconi_pysb.py', 'w') as f:
#     f.write(pysb_model)

plt.show()
