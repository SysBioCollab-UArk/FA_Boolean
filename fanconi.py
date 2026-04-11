from util import *
from fanconi_pysb import model as model_pysb
import re
import csv
import os

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

# CONDITIONS = [
#     {'mutations': None, 'initials': [("DSB", True)]},  # FIG 2A
#     {'mutations': None, 'initials': [("ADD", True)]},  # FIG 2B
#     {'mutations': None, 'initials': [("ICL", True)]},  # FIG 2C
#     {'mutations': [("ICL", True)], 'initials': None},  # FIG 3
#     {'mutations': [("FAcore", False)], 'initials': [("ICL", True)]},  # FIG 4A
#     {'mutations': [("FAcore", False), ("ICL", True)], 'initials': None},  # FIG 4B
#     {'mutations': [("FANCD1N", False)], 'initials': [("ICL", True)]},  # FIG 5A
#     {'mutations': [("FANCD1N", False), ("ICL", True)], 'initials': None}  # FIG 5B
# ]

initials = [("ICL", True)]
exclude = {'ICL', 'DSB', 'ADD'}  # set
nodes = sorted(list(set(re.findall(r'(\w+)\s*=', Initial_conditions)) - exclude))

CONDITIONS = [{'mutations': None, 'initials': initials}]  # wild-type
for node in nodes:
    CONDITIONS.append({'mutations': [(node, True)], 'initials': initials})
    CONDITIONS.append({'mutations': [(node, False)], 'initials': initials})

SYNCH = False
ASYNCH = False
PYSB = True

if ASYNCH or PYSB:
    OVERWRITE = True

    species_to_plot = [
        ['ICL', 'ADD', 'DSB'],
        ['FANCM', 'FAcore', 'FANCD2I', 'MUS81', 'FANCJBRCA1', 'XPF', 'FAN1'],
        ['PCNATLS', 'MRN', 'BRCA1', 'ssDNARPA', 'FANCD1N', 'RAD51', 'HRR', 'USP1'],
        ['KU', 'DNAPK', 'NHEJ'],
        ['ATR', 'ATM', 'p53', 'CHK1', 'CHK2', 'H2AX', 'CHKREC']
    ]

    stoch_kwargs = {
        'species': ['ICL', 'ADD', 'DSB'],  # must be a list
        'threshold': 0.05
    }
    if len(stoch_kwargs.get('species', [])) > 0:
        csvwriters = []
        if ASYNCH:
            filename = 'REPAIR_TIMES_asynch.csv'
            mode = 'a' if os.path.exists(filename) and not OVERWRITE else 'w'
            csvfile = open(filename, mode)
            csvwriter_asynch = csv.writer(csvfile, delimiter=',')
            csvwriters.append(csvwriter_asynch)
        if PYSB:
            filename = 'REPAIR_TIMES_pysb.csv'
            mode = 'a' if os.path.exists(filename) and not OVERWRITE else 'w'
            csvfile = open(filename, mode)
            csvwriter_psyb = csv.writer(csvfile, delimiter=',')
            csvwriters.append(csvwriter_psyb)
        if mode == 'w':
            for csvwriter in csvwriters:
                csvwriter.writerow(['mutations', 'initials'] + stoch_kwargs['species'])

# TODO: for synch updating, find attractors for all 2^28 possible initial conditions

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
            if SYNCH or ASYNCH:
                # set the initial state of the node
                # NOTE: to ensure that we match the whole node name and not just a part of it (e.g., BRCA1 vs
                # FANCJBRCA1), use the ^ character together with the MULTILINE flag
                pattern = r'^(\s*%s\s*=\s*).*' % mut[0]
                if re.search(pattern, model_text, flags=re.MULTILINE) is not None:
                    model_text = re.sub(pattern, r'\1%s' % mut[1], model_text, flags=re.MULTILINE)
                else:
                    print("Error: can't find node '%s' in the initial conditions." % mut[0])
                    quit()
                # remove the update rule
                pattern = r'\s*\d+:\s*%s\*.*' % mut[0]
                if re.search(pattern, model_text) is not None:
                    model_text = re.sub(pattern, '', model_text)
                else:
                    print("Error: can't find node '%s' in the rules." % mut[0])
                    quit()

            # === PySB model ===
            if PYSB:
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
            if SYNCH or ASYNCH:
                # set the initial state of the node
                # NOTE: to ensure that we match the whole node name and not just a part of it (e.g., BRCA1 vs
                # FANCJBRCA1), use the ^ character together with the MULTILINE flag
                pattern = r'^(\s*%s\s*=\s*).*' % init[0]
                # m = re.search(r'%s\s*=\s*\w+' % init[0], model_text)
                if re.search(pattern, model_text, flags=re.MULTILINE) is not None:
                    model_text = re.sub(pattern, r'\1%s' % init[1], model_text, flags=re.MULTILINE)
                else:
                    print("Error: can't find node '%s' in initial conditions." % init[0])
                    quit()

            # === PySB model ===
            if PYSB:
                # set the initial state of the node
                param_values['%s_%s_init' % (init[0], init[1])] = 1.0
                param_values['%s_%s_init' % (init[0], not init[1])] = 0.0

    # if no mutations or initial conditions, this is the base model
    if prefix == "":
        prefix = "base"

    # synchronous updating
    if SYNCH:
        n_iter = 100
        run_FA_Boolean_synch(model_text, n_iter, detect_cycles=True, outfile='%s_synch.pdf' % prefix)

    # asynchronous updating
    if ASYNCH:
        n_iter = 50
        n_stoch_sims = 100
        run_FA_Boolean_asynch(model_text, n_iter, n_stoch_sims, species_to_plot, verbose=False,
                              outfile='%s_asynch.pdf' % prefix, **stoch_kwargs)

    # PySB general asynchronous (Gillespie) simulations
    if PYSB:
        t_end = 50
        n_stoch_sims = 1000
        n_reps = 100
        for n in range(n_reps):
            print(n, end='\n' if (n + 1) % 20 == 0 or n == n_reps - 1 else ' ')
            outfile = None if n == 0 else '%s_pysb.pdf' % prefix
            pvals = {k: v for k, v in param_values.items()}
            t_thresh = run_FA_Boolean_pysb(model_pysb, t_end, n_stoch_sims, species_to_plot, param_values=pvals,
                                           verbose=False, outfile=outfile, **stoch_kwargs)
            mutations = [''] if cond['mutations'] is None else [str(m) for mut in cond['mutations'] for m in mut]
            initials = [''] if cond['initials'] is None else [str(m) for mut in cond['initials'] for m in mut]
            if t_thresh[0] == np.inf:
                t_thresh[1:] = [np.nan] * len(t_thresh[1:])
                # print(t_thresh)
            csvwriter_psyb.writerow(['_'.join(mutations), '_'.join(initials)] + t_thresh)

# Export model generated by Boolean2Rules so don't have to regenerate it each time (very time-consuming)
# from pysb.export import export
# pysb_model = export(model, 'pysb_flat')
# with open('fanconi_pysb.py', 'w') as f:
#     f.write(pysb_model)

plt.show()
