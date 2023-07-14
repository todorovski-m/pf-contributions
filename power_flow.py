import time
import os
import pandas as pd
import numpy as np
import scipy.io as sio


MPC_COLUMN_NAMES = {
    "BUS": ["BUS_I", "BUS_TYPE", "PD", "QD", "GS", "BS", "BUS_AREA", "VM", "VA", "BASE_KV", "ZONE", "VMAX", "VMIN"],
    "BRANCH": ["F_BUS", "T_BUS", "BR_R", "BR_X", "BR_B", "RATE_A", "RATE_B", "RATE_C", "TAP", "SHIFT", "BR_STATUS",
               "ANGMIN", "ANGMAX", "PF", "QF", "PT", "QT"],
    "GEN": ["GEN_BUS", "PG", "QG", "QMAX", "QMIN", "VG", "MBASE", "GEN_STATUS", "PMAX", "PMIN", "PC1", "PC2",
            "QC1MIN", "QC1MAX", "QC2MIN", "QC2MAX", "RAMP_AGC", "RAMP_10", "RAMP_30", "RAMP_Q", "APF"],
    "BUS_SELECTED": ["BUS_I", "BUS_TYPE", "PD", "QD", "VM", "VA"],
    "BRANCH_SELECTED": ["F_BUS", "T_BUS", "BR_STATUS", "TAP", "PF", "QF", "PT", "QT"],
    "GEN_SELECTED": ["GEN_BUS", "GEN_STATUS", "PG", "QG"],
}


def runpf(case_file):
    start = time.time()
    case_file = "'" + case_file + "'"
    command = f'octave --eval "my_runpf({case_file});"'
    os.system(command)
    results = sio.loadmat("results/results.mat")
    if results["success"] == 0:
        raise Exception("Newton's method failed.")
    print("matpower time = {:.3f}\n".format(time.time()-start))

    bus = results["bus"]
    bus = pd.DataFrame(np.array(bus), columns=MPC_COLUMN_NAMES["BUS"])
    bus = bus[MPC_COLUMN_NAMES["BUS_SELECTED"]]
    for c in ["BUS_I", "BUS_TYPE"]:
        bus[c] = bus[c].astype(int)
    bus["TYPE"] = "L"
    bus["ID"] = ["L"+str(i) for i in bus["BUS_I"]]
    bus["PS"] = np.array(results["PS"])
    bus["QS"] = np.array(results["QS"])

    branch = results["branch"]
    branch = pd.DataFrame(np.array(branch), columns=MPC_COLUMN_NAMES["BRANCH"])
    idx = branch["TAP"] == 0
    branch.loc[idx, "TAP"] = 1
    branch = branch[MPC_COLUMN_NAMES["BRANCH_SELECTED"]]
    for c in ["F_BUS", "T_BUS", "BR_STATUS"]:
        branch[c] = branch[c].astype(int)
    branch["QSF"] = np.array(results["QSF"])
    branch["QST"] = np.array(results["QST"])
    branch["TYPE"] = "E"
    n = len(branch)
    branch["ID"] = ["E"+str(i) for i in range(1, n+1)]

    gen = results["gen"]
    gen = pd.DataFrame(np.array(gen), columns=MPC_COLUMN_NAMES["GEN"])
    gen = gen[MPC_COLUMN_NAMES["GEN_SELECTED"]]
    for c in ["GEN_BUS", "GEN_STATUS"]:
        gen[c] = gen[c].astype(int)
    gen["TYPE"] = "G"
    gen["ID"] = ["G"+str(i) for i in gen["GEN_BUS"]]

    if "buslocation" in results:
        xy = results["buslocation"]
    else:
        xy = None

    return bus, gen, branch, xy
