import time
import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import logging
from power_flow import runpf


def replaceDuplicates(names):
    """
    Function to append duplicate strings by alphanumeric strings to make all
    strings in the array unique
    """
    count = {}
    for i in range(len(names)):
        if names[i] not in count:
            count[names[i]] = 1
        else:
            count[names[i]] += 1
            n = count[names[i]]
            # add suffix .2, .3 etc for duplicate names
            names[i] += "." + str(n)
    for i in range(len(names)):
        if names[i] in count:
            if count[names[i]] > 1:
                # add suffix .1 to the first instance that was present multiple times
                names[i] += ".1"
    return names


def makeActiveFlows(bus, gen, branch):
    # since we are going to swap columns we make a copy of branch
    # otherwise changes will made to the original variable as well
    branch = branch.copy()
    # Swap branch orientation when both PF and PT are negative
    idx = (branch["PF"] < 0) & (branch["PT"] < 0)
    if any(idx):
        c = ["F_BUS", "T_BUS", "PF", "PT"]
        logging.info(f"Swap branch orientation when both PF and PT are negative ({branch[idx].shape[0]})")
        logging.debug(f"Before swap\n{branch.loc[idx, c]}")
        # Swap F_BUS and T_BUS
        temp = branch.loc[idx, "F_BUS"]
        branch.loc[idx, "F_BUS"] = branch.loc[idx, "T_BUS"]
        branch.loc[idx, "T_BUS"] = temp
        # Swap PF and PT and change sign
        temp = branch.loc[idx, "PF"]
        branch.loc[idx, "PF"] = -branch.loc[idx, "PT"]
        branch.loc[idx, "PT"] = -temp
        logging.debug(f"After swap\n{branch.loc[idx, c]}")

    # Active power flows - network elements
    cols = ["F_BUS", "T_BUS", "PF", "PT", "TYPE", "TAP", "ID"]
    flows = branch[cols].copy()

    # Active power flows - loads with positive power
    idx = bus["PD"] > 0
    loads = bus.loc[idx, ["BUS_I", "PD", "TYPE", "ID"]]
    loads.insert(1, "T_BUS", 0)
    loads = loads.rename(columns={"BUS_I": "F_BUS", "PD": "PF"})
    loads["PT"] = loads["PF"]
    loads["TAP"] = 1

    # Active power flows - loads with negative power, converted to generators
    idx = bus["PD"] < 0
    loads_to_gens = bus.loc[idx, ["BUS_I", "PD", "TYPE", "ID"]]
    loads_to_gens.insert(0, "F_BUS", 0)
    loads_to_gens = loads_to_gens.rename(
        columns={"BUS_I": "T_BUS", "PD": "PF"})
    loads_to_gens["PF"] = -loads_to_gens["PF"]
    loads_to_gens["PT"] = loads_to_gens["PF"]
    loads_to_gens["TYPE"] = "G"
    loads_to_gens["TAP"] = 1
    # loads_to_gens["ID"] = [i.replace("L", "G") for i in loads_to_gens["ID"]]
    if any(idx):
        logging.info("Loads with negative power converted to generators ({})".format(
            loads_to_gens.shape[0]))

    # Active power flows - generators with positive power
    idx = gen["PG"] > 0
    gens = gen.loc[idx, ["GEN_BUS", "PG", "TYPE", "ID"]]
    gens.insert(0, "F_BUS", 0)
    gens = gens.rename(columns={"GEN_BUS": "T_BUS", "PG": "PT"})
    gens["PF"] = gens["PT"]
    gens["TAP"] = 1

    # Active power flows - generators with negative power, converted to loads
    idx = gen["PG"] < 0
    gens_to_loads = gen.loc[idx, ["GEN_BUS", "PG", "TYPE", "ID"]]
    gens_to_loads.insert(1, "T_BUS", 0)
    gens_to_loads = gens_to_loads.rename(
        columns={"GEN_BUS": "F_BUS", "PG": "PT"})
    gens_to_loads["PT"] = -gens_to_loads["PT"]
    gens_to_loads["PF"] = gens_to_loads["PT"]
    gens_to_loads["TYPE"] = "L"
    gens_to_loads["TAP"] = 1
    # gens_to_loads["ID"] = [i.replace("G", "L") for i in gens_to_loads["ID"]]
    if any(idx):
        logging.info("Generators with negative power converted to loads ({})".format(
            gens_to_loads.shape[0]))
    flows = pd.concat([flows, loads, loads_to_gens, gens, gens_to_loads])

    # Replace possible duplicates in ID
    ID = list(flows["ID"])
    flows["ID"] = replaceDuplicates(ID)
    flows = flows.set_index("ID")
    return flows[cols[:-1]]


def makeReactiveFlows(bus, gen, branch):
    # since we are going to swap columns we make a copy of branch
    # otherwise changes will made to the original variable as well
    branch = branch.copy()
    # Swap branch orientation when both QF and QT are negative
    idx = (branch["QF"] < 0) & (branch["QT"] < 0)
    if any(idx):
        c = ["F_BUS", "T_BUS", "QF", "QT", "QSF", "QST", "TAP"]
        logging.info(f"Swap branch orientation when both QF and QT are negative ({branch[idx].shape[0]})")
        logging.debug(f"Before swap\n{branch.loc[idx, c]}")
        # Swap F_BUS and T_BUS
        temp = branch.loc[idx, "F_BUS"]
        branch.loc[idx, "F_BUS"] = branch.loc[idx, "T_BUS"]
        branch.loc[idx, "T_BUS"] = temp
        # Swap QF and QT and change sign
        temp = branch.loc[idx, "QF"]
        branch.loc[idx, "QF"] = -branch.loc[idx, "QT"]
        branch.loc[idx, "QT"] = -temp
        # Swap QSF and QST
        temp = branch.loc[idx, "QSF"]
        branch.loc[idx, "QSF"] = branch.loc[idx, "QST"]
        branch.loc[idx, "QST"] = temp
        logging.debug(f"After swap\n{branch.loc[idx, c]}")

    # Reactive power flows - network elements
    cols = ["F_BUS", "T_BUS", "QF", "QF1", "QT1",
            "QT", "QSF", "QST", "TAP", "TYPE", "ID"]
    branch["QF1"] = branch["QF"] + branch["QSF"]
    branch["QT1"] = branch["QT"] - branch["QST"]
    idx = branch["TAP"] == 0
    branch.loc[idx, "TAP"] = 1
    flows = branch[cols].copy()

    # Reactive power flows - loads with positive power
    idx = bus["QD"] > 0
    loads = bus.loc[idx, ["BUS_I", "QD", "TYPE", "ID"]]
    loads.insert(1, "T_BUS", 0)
    loads = loads.rename(columns={"BUS_I": "F_BUS", "QD": "QF"})
    loads["QT"] = loads["QF"]
    loads["QF1"] = loads["QF"]
    loads["QT1"] = loads["QF"]
    loads["QSF"] = 0
    loads["QST"] = 0
    loads["TAP"] = 1

    # Reactive power flows - loads with negative power converted to generators
    idx = bus["QD"] < 0
    loads_to_gens = bus.loc[idx, ["BUS_I", "QD", "TYPE", "ID"]]
    loads_to_gens.insert(0, "F_BUS", 0)
    loads_to_gens = loads_to_gens.rename(
        columns={"BUS_I": "T_BUS", "QD": "QF"})
    loads_to_gens["QF"] = -loads_to_gens["QF"]
    loads_to_gens["QT"] = loads_to_gens["QF"]
    loads_to_gens["QF1"] = loads_to_gens["QF"]
    loads_to_gens["QT1"] = loads_to_gens["QF"]
    loads_to_gens["QSF"] = 0
    loads_to_gens["QST"] = 0
    loads_to_gens["TYPE"] = "G"
    # loads_to_gens["ID"] = [i.replace("L", "G") for i in loads_to_gens["ID"]]
    loads_to_gens["TAP"] = 1
    if any(idx):
        logging.info(f"Loads with negative power converted to generators ({loads_to_gens.shape[0]})")
        logging.debug(f"Loads with negative power converted to generators\n{loads_to_gens}")

    # Reactive power flows - shunts with positive power converted to generators
    idx = bus["QS"] > 0
    shunts_to_gens = bus.loc[idx, ["BUS_I", "QS", "TYPE", "ID"]]
    shunts_to_gens.insert(0, "F_BUS", 0)
    shunts_to_gens = shunts_to_gens.rename(
        columns={"BUS_I": "T_BUS", "QS": "QF"})
    shunts_to_gens["QT"] = shunts_to_gens["QF"]
    shunts_to_gens["QF1"] = shunts_to_gens["QF"]
    shunts_to_gens["QT1"] = shunts_to_gens["QF"]
    shunts_to_gens["QSF"] = 0
    shunts_to_gens["QST"] = 0
    shunts_to_gens["TYPE"] = "G"
    shunts_to_gens["ID"] = [i.replace("L", "S") for i in shunts_to_gens["ID"]]
    shunts_to_gens["TAP"] = 1
    if any(idx):
        logging.info(f"Shunts with positive power converted to generators ({shunts_to_gens.shape[0]})")
        logging.debug(f"Shunts with positive power converted to generators\n{shunts_to_gens}")

    # Reactive power flows - shunts with negative power converted to loads
    idx = bus["QS"] < 0
    shunts_to_loads = bus.loc[idx, ["BUS_I", "QS", "TYPE", "ID"]]
    shunts_to_loads.insert(1, "T_BUS", 0)
    shunts_to_loads = shunts_to_loads.rename(
        columns={"BUS_I": "F_BUS", "QS": "QT"})
    shunts_to_loads["QT"] = -shunts_to_loads["QT"]
    shunts_to_loads["QF"] = shunts_to_loads["QT"]
    shunts_to_loads["QF1"] = shunts_to_loads["QT"]
    shunts_to_loads["QT1"] = shunts_to_loads["QT"]
    shunts_to_loads["QSF"] = 0
    shunts_to_loads["QST"] = 0
    shunts_to_loads["TYPE"] = "L"
    shunts_to_loads["ID"] = [i.replace("L", "S") for i in shunts_to_loads["ID"]]
    shunts_to_loads["TAP"] = 1
    if any(idx):
        logging.info(f"Shunts with negative power converted to loads ({shunts_to_loads.shape[0]})")
        logging.debug(f"Shunts with negative power converted to loads\n{shunts_to_loads}")

    # Reactive power flows - generators with positive power
    idx = gen["QG"] > 0
    gens = gen.loc[idx, ["GEN_BUS", "QG", "TYPE", "ID"]]
    gens.insert(0, "F_BUS", 0)
    gens = gens.rename(columns={"GEN_BUS": "T_BUS", "QG": "QT"})
    gens["QF"] = gens["QT"]
    gens["QF1"] = gens["QF"]
    gens["QT1"] = gens["QF"]
    gens["QSF"] = 0
    gens["QST"] = 0
    gens["TAP"] = 1

    # Reactive power flows - generators with negative power converted to loads
    idx = gen["QG"] < 0
    gens_to_loads = gen.loc[idx, ["GEN_BUS", "QG", "TYPE", "ID"]]
    gens_to_loads.insert(1, "T_BUS", 0)
    gens_to_loads = gens_to_loads.rename(
        columns={"GEN_BUS": "F_BUS", "QG": "QT"})
    gens_to_loads["QT"] = -gens_to_loads["QT"]
    gens_to_loads["QF"] = gens_to_loads["QT"]
    gens_to_loads["QF1"] = gens_to_loads["QT"]
    gens_to_loads["QT1"] = gens_to_loads["QT"]
    gens_to_loads["QSF"] = 0
    gens_to_loads["QST"] = 0
    gens_to_loads["TYPE"] = "L"
    # gens_to_loads["ID"] = [i.replace("G", "L") for i in gens_to_loads["ID"]]
    gens_to_loads["TAP"] = 1
    if any(idx):
        logging.info(f"Generators with negative power converted to loads ({gens_to_loads.shape[0]})")
        logging.debug(f"Generators with negative power converted to loads\n{gens_to_loads}")

    # Branches that are receiving reactive powers from both ends, converted to loads
    idx = (flows["TYPE"] == "E") & (flows["QF"] > 0) & (flows["QT"] < 0)
    if any(idx):
        F_BUS = list(flows.loc[idx, "F_BUS"])
        T_BUS = list(flows.loc[idx, "T_BUS"])
        QF = list(flows.loc[idx, "QF"])
        QT = list(-flows.loc[idx, "QT"])
        ID = list(flows.loc[idx, "ID"])
        logging.info(
            f"Branches that are receiving reactive powers from both ends converted to loads ({flows[idx].shape[0]})")
        logging.debug(f"Branches that are receiving reactive powers from both ends converted to loads\n{flows[idx]}")
        flows = flows.drop(flows[idx].index)
        branches_to_loads = []
        for i in range(len(F_BUS)):
            d = {
                "F_BUS": F_BUS[i],
                "T_BUS": 0,
                "QF": QF[i],
                "ID": ID[i] + "_QF",
            }
            branches_to_loads.append(d)
            d = {
                "F_BUS": T_BUS[i],
                "T_BUS": 0,
                "QF": QT[i],
                "ID": ID[i] + "_QT",
            }
            branches_to_loads.append(d)
        branches_to_loads = pd.DataFrame(branches_to_loads)
        Q = branches_to_loads["QF"]
        branches_to_loads["QT"] = Q
        branches_to_loads["QF1"] = Q
        branches_to_loads["QT1"] = Q
        branches_to_loads["QSF"] = 0
        branches_to_loads["QST"] = 0
        branches_to_loads["TYPE"] = "L"
        branches_to_loads["TAP"] = 1
        logging.debug(
            f"Loads from branches that are receiving reactive powers from both ends\n{branches_to_loads}")
    else:
        branches_to_loads = pd.DataFrame()

    # Transformers with negative QSF converted to loads
    idx = flows["QSF"] < 0
    qsf_to_loads = flows.loc[idx, ["F_BUS", "QSF"]]
    flows.loc[idx, "QSF"] = 0
    flows.loc[idx, "QF1"] = flows.loc[idx, "QF"]
    qsf_to_loads["T_BUS"] = 0
    Q = -qsf_to_loads["QSF"]
    qsf_to_loads["QF"] = Q
    qsf_to_loads["QT"] = Q
    qsf_to_loads["QF1"] = Q
    qsf_to_loads["QT1"] = Q
    qsf_to_loads["QSF"] = 0
    qsf_to_loads["QST"] = 0
    qsf_to_loads["TYPE"] = "L"
    qsf_to_loads["ID"] = flows.loc[idx, "ID"] + \
        ".S" + flows.loc[idx, "F_BUS"].astype(str)
    qsf_to_loads["TAP"] = 1
    if any(idx):
        logging.info(f"Transformers with negative QSF converted to loads ({qsf_to_loads.shape[0]})")
        logging.debug(f"Transformers with negative QSF converted to loads\n{qsf_to_loads}")

    # Transformers with negative QST converted to loads
    idx = flows["QST"] < 0
    qst_to_loads = flows.loc[idx, ["T_BUS", "QST"]]
    qst_to_loads = qst_to_loads.rename(columns={"T_BUS": "F_BUS"})
    flows.loc[idx, "QST"] = 0
    flows.loc[idx, "QT1"] = flows.loc[idx, "QT"]
    qst_to_loads["T_BUS"] = 0
    Q = -qst_to_loads["QST"]
    qst_to_loads["QF"] = Q
    qst_to_loads["QT"] = Q
    qst_to_loads["QF1"] = Q
    qst_to_loads["QT1"] = Q
    qst_to_loads["QSF"] = 0
    qst_to_loads["QST"] = 0
    qst_to_loads["TYPE"] = "L"
    qst_to_loads["ID"] = flows.loc[idx, "ID"] + \
        ".S" + flows.loc[idx, "T_BUS"].astype(str)
    qst_to_loads["TAP"] = 1
    if any(idx):
        logging.info(f"Transformers with negative QST converted to loads ({qst_to_loads.shape[0]})")
        logging.debug(f"Transformers with negative QST converted to loads\n{qst_to_loads}")
    flows = pd.concat([flows, loads, loads_to_gens, gens, gens_to_loads,
                       qsf_to_loads, qst_to_loads, branches_to_loads,
                       shunts_to_gens, shunts_to_loads])

    # Replace possible duplicates in ID
    ID = list(flows["ID"])
    flows["ID"] = replaceDuplicates(ID)
    flows = flows.set_index("ID")
    return flows[cols[:-1]]


def unpack(dict, power_key):
    branches = list(dict.keys())
    val = list(dict.values())
    P = [item[power_key] for item in val]
    side = [item["side"] for item in val]
    sign = [item["sign"] for item in val]
    return branches, P, side, sign


def expand(x, side):
    return x + "_" + side


def other(side):
    opposite = {"PF": "PT", "PT": "PF", "QF": "QT", "QT": "QF"}
    return opposite[side]


def makeTable(A, branch_names, gen_names, F_BUS, T_BUS, TYPE, DONE, TAP):
    A = pd.DataFrame(A, index=branch_names, columns=gen_names)
    A["F_BUS"] = F_BUS
    A["T_BUS"] = T_BUS
    A["TYPE"] = TYPE
    A["DONE"] = DONE
    if TAP:
        A["TAP"] = TAP
    else:
        A["TAP"] = 1
    return A


def drawGraph(DATA, bus_pos, title=""):
    idx = (DATA["F_BUS"] > 0) & (DATA["T_BUS"] > 0) & (DATA["TAP"] == 1)
    LINES = nx.from_pandas_edgelist(
        DATA[idx], source="F_BUS", target="T_BUS", create_using=nx.Graph)

    idx = (DATA["F_BUS"] > 0) & (DATA["T_BUS"] > 0) & (DATA["TAP"] != 1)
    TRANSF = nx.from_pandas_edgelist(
        DATA[idx], source="F_BUS", target="T_BUS", edge_attr=True, create_using=nx.Graph)
    TRANSF_LABELS = nx.get_edge_attributes(TRANSF, "TAP")

    idx = (DATA["F_BUS"] > 0) & (DATA["T_BUS"] > 0) & DATA["DONE"]
    NOT_DONE = nx.from_pandas_edgelist(
        DATA[idx], source="F_BUS", target="T_BUS", create_using=nx.Graph)
    nl = len(LINES.edges()) + len(TRANSF.edges())
    nd = len(NOT_DONE.edges())
    if not bus_pos:
        bus_pos = nx.spring_layout(LINES)
    plt.figure()
    plt.title(title + " ({}/{})".format(nd, nl))
    nx.draw(LINES, bus_pos, edge_color="#C0C0C0", node_size=35)
    nx.draw(TRANSF, bus_pos, edge_color="#C0C0C0", node_size=35)
    nx.draw_networkx_edge_labels(TRANSF, bus_pos, edge_labels=TRANSF_LABELS)
    nx.draw(NOT_DONE, bus_pos, edge_color="#FF0000", node_size=35, width=2)
    node_ind = nx.draw_networkx_nodes(LINES, bus_pos, node_color="#C0C0C0")
    node_ind.set_edgecolor("k")
    nx.draw_networkx_labels(LINES, bus_pos)
    plt.show()


def makeAPFC(bus, gen, branch, bus_pos, itermax, graph):
    logging.info("Starting ...")
    # Active power flows
    flows = makeActiveFlows(bus, gen, branch)
    m, n = flows.shape
    flows_mat = flows.values
    flows_row = {v: k for v, k in zip(flows.index, range(m))}
    flows_col = {v: k for v, k in zip(flows.columns, range(n))}

    # Bus names
    bus_names = pd.concat([flows["F_BUS"], flows["T_BUS"]]).unique()
    bus_names = list(bus_names)
    bus_names.sort()
    bus_names.pop(0)

    # Sources of active power
    idx = flows["TYPE"] == "G"
    gen_names = list(flows[idx].index)

    # A[k]: set of branches incident to bus k where each branch total active power flow direction is away from bus k
    # B[k]: set of branches incident to bus k where each branch total active power flow direction is towards bus k
    A, B = {}, {}
    for k in bus_names:
        A[k] = {}
        B[k] = {}

    # make directed graph as a convenient way to explore bus-branch connections
    # reset index so that ID becomes column that will be used in edge_attr
    G = nx.from_pandas_edgelist(
        flows.reset_index(), source="F_BUS", target="T_BUS", edge_attr=True, create_using=nx.MultiDiGraph)

    # find branches incident to each bus and categorize them into sets A and B
    # save the active power P that leaves or enters bus k, branch side PF
    # "from" or PT "to" and the multiple sign (1 or -1) that is used to obtain
    # the correct power flow according to definition of sets A and B
    for k in bus_names:
        # branches that are oriented "to" buses
        for (_, _, data) in G.in_edges(k, data=True):
            ID, PT = data["ID"], data["PT"]
            # branch power is PT taken as is if positive or multiplied by -1 if negative (branch reversal)
            if PT > 0:
                B[k][ID] = {"P": PT, "side": "PT", "sign": 1}
            else:
                A[k][ID] = {"P": -PT, "side": "PT", "sign": -1}
        # branches that are oriented "from" buses
        for (_, _, data) in G.out_edges(k, data=True):
            ID, PF = data["ID"], data["PF"]
            # branch power is PF taken as is if positive or multiplied by -1 if negative (branch reversal)
            if PF > 0:
                A[k][ID] = {"P": PF, "side": "PF", "sign": 1}
            else:
                B[k][ID] = {"P": -PF, "side": "PF", "sign": -1}

    logging.debug(f"{15*'='} Sets A and B {15*'='}")
    for k in bus_names:
        logging.debug(f"A[{k}] = {A[k]}")
        logging.debug(f"B[{k}] = {B[k]}\n")
    logging.debug(55*"=")

    # Components of active power flows
    # Each branch is listed twice with suffixes "_PF" and "_PT" for both ends
    branch_names = [
        i + suffix for i in flows.index for suffix in ["_PF", "_PT"]]
    # List "from" and "to" buses twice
    F_BUS = [i for i in flows["F_BUS"] for _ in range(2)]
    T_BUS = [i for i in flows["T_BUS"] for _ in range(2)]
    TYPE = [i for i in flows["TYPE"] for _ in range(2)]
    TAP = [i for i in flows["TAP"] for _ in range(2)]
    # Interleave PF and PT
    PFT = [val for pair in zip(flows["PF"], flows["PT"]) for val in pair]
    # numpy array that will hold the results
    nr, nc = len(branch_names), len(gen_names)
    APFC = np.zeros((nr, nc))
    # Dictionaries that hold row and column indices for each element in branch_names and gen_names
    row_i = {v: k for v, k in zip(branch_names, range(nr))}
    col_i = {v: k for v, k in zip(gen_names, range(nc))}
    # List that indicates whether APFC have been calculated
    DONE = nr * [False]

    # Generator branches are supplied by respective generators only,
    # therefore we know their active power flow components (APFC).
    # They are starting points in the procedure.
    idx = flows["TYPE"] == "G"
    PG = flows.loc[idx, "PF"]
    ig = flows[idx].index
    for g, P in zip(ig, PG):
        i, j = [row_i[g + suffix] for suffix in ["_PF", "_PT"]]
        gg = col_i[g]
        APFC[i, gg] = P
        APFC[j, gg] = P
        DONE[i] = True
        DONE[j] = True
    known_components = list(gen_names)
    logging.debug(f"Initial known components: {gen_names}")

    # Calculate APFC of other branches
    processed_buses = []
    repeat = True
    iter = 0
    start = time.time()
    while repeat and iter < itermax:
        iter += 1
        logging.debug(f"{10*'='} Start of iteration {iter} {10*'='}")
        iter_start = time.time()
        # Find bus k for which APFC are known for all branches in B[k]
        for k in bus_names:
            if k not in processed_buses:
                logging.debug(f"Processing bus k = {k}")
                check = all(b in known_components for b in B[k].keys())
                if check:
                    # If all APFC for branches in B[k] are known we may calculate
                    # the APFC for the branches in A[k].
                    branches_a, P_a, side_a, sign_a = unpack(A[k], "P")
                    branches_b, P_b, side_b, sign_b = unpack(B[k], "P")
                    logging.debug(f"B[{k}]: {branches_b}, {P_b}, {side_b}, {sign_b}, sum(P) = {sum(P_b)}")
                    logging.debug(f"A[{k}]: {branches_a}, {P_a}, {side_a}, {sign_a}, sum(P) = {sum(P_a)}")
                    # Loop through all branches in A[k] and calculate their APFC
                    for i, P, side_i, sign_i in zip(branches_a, P_a, side_a, sign_a):
                        ii = row_i[expand(i, side_i)]
                        logging.debug(f"outflow: {i}_{side_i} ({F_BUS[ii]}-{T_BUS[ii]})")
                        # autopep8: off
                        # if abs(P) <= 1e-4:
                            # DONE[row_i[expand(i, "PF")]] = True
                            # DONE[row_i[expand(i, "PT")]] = True
                            # known_components.append(i)
                            # continue
                        # f = P/sum(P_b)
                        if abs(P) <= 1e-4 and sum(P_b) == 0:
                            f = 0
                            DONE[row_i[expand(i, "PF")]] = True
                            DONE[row_i[expand(i, "PT")]] = True
                        else:
                            f = P/sum(P_b)
                        # autopep8: on
                        logging.debug(f"f = {P}/{sum(P_b)} = {f}")
                        # logging.info("f = {}/{} = {}".format(P, sum(P_b), f))
                        # Loop through all components of all branches in B[k]
                        for g in gen_names:
                            for j, side_j, sign_j in zip(branches_b, side_b, sign_b):
                                # APFC for branch j (from B[k]) and generator g for the side (side_j) that brings power into bus k
                                jj = row_i[expand(j, side_j)]
                                gg = col_i[g]
                                Pjg = APFC[jj, gg]
                                # APFC for branch i (from A[k]) and generator g for the side (side_i) that takes out power from bus k
                                # for both branches we must multiple with the power flow orientation sign
                                CP = sign_i * sign_j * f * Pjg
                                # logging.info("in: {}.{} = {}, CP = {}".format(j, g, Pjg, CP))
                                if CP != 0:
                                    logging.debug(f"inflow: {j}.{g} = {Pjg}, CP = {CP} " +
                                                  f"--> component for: {expand(i, side_i)}.{g} ({F_BUS[ii]}-{T_BUS[ii]})")
                                ii = row_i[expand(i, side_i)]
                                APFC[ii, gg] += CP
                                DONE[ii] = True
                                # APFC for the other side of branch i
                                # P_other_side = flows.loc[i, other(side_i)] # <-- slow
                                P_other_side = flows_mat[flows_row[i], flows_col[other(side_i)]]
                                ii = row_i[expand(i, other(side_i))]
                                CP_other_side = sign_i * P_other_side/P * CP
                                APFC[ii, gg] += CP_other_side
                                if CP_other_side != 0:
                                    logging.debug(
                                        f"P_other_side = {CP_other_side} --> component for: {expand(i, other(side_i))}.{g}")
                                DONE[ii] = True
                        known_components.append(i)
                    processed_buses.append(k)
                else:
                    logging.debug("    not all APFC known, skipping ...")
        repeat = not all(DONE)
        logging.info("iter = {}, time = {:.3f}".format(iter, time.time()-iter_start))
        if graph:
            T = makeTable(APFC, branch_names, gen_names,
                          F_BUS, T_BUS, TYPE, DONE, TAP)
            drawGraph(T, bus_pos, "iter = {}".format(iter))
    APFC = makeTable(APFC, branch_names, gen_names,
                     F_BUS, T_BUS, TYPE, DONE, TAP)
    APFC["Total"] = APFC[gen_names].sum(axis=1)
    APFC["Error"] = APFC["Total"] - PFT
    logging.info("total time = {:.3f}".format(time.time()-start))
    logging.info("APFC: nr = {}, iters = {}, DONE = {}".format(nr, iter, not repeat))

    if repeat:
        for k in bus_names:
            if k not in processed_buses:
                logging.info("\nbus: {}".format(k))
                logging.info("IN")
                for b in B[k]:
                    f = flows.loc[b, "F_BUS"]
                    t = flows.loc[b, "T_BUS"]
                    ii = row_i[expand(b, B[k][b]["side"])]
                    logging.info("  branch: {} ({}-{}), {}, {}".format(b,
                                                                       f, t, B[k][b], DONE[ii]))
                logging.info("OUT")
                for b in A[k]:
                    f = flows.loc[b, "F_BUS"]
                    t = flows.loc[b, "T_BUS"]
                    ii = row_i[expand(b, A[k][b]["side"])]
                    logging.info("  branch: {} ({}-{}), {}, {}".format(b,
                                                                       f, t, A[k][b], DONE[ii]))

    logging.info("End of calculations\n")
    return flows, APFC


def makeRPFC(bus, gen, branch, bus_pos, itermax, graph):
    logging.info("Starting ...")
    # Reactive power flows
    flows = makeReactiveFlows(bus, gen, branch)
    m, n = flows.shape
    flows_mat = flows.values
    flows_row = {v: k for v, k in zip(flows.index, range(m))}
    flows_col = {v: k for v, k in zip(flows.columns, range(n))}

    # Bus names
    bus_names = pd.concat([flows["F_BUS"], flows["T_BUS"]]).unique()
    bus_names = list(bus_names)
    bus_names.sort()
    bus_names.pop(0)

    # Sources of reactive power
    idx = flows["TYPE"] == "G"
    gen_names = list(flows[idx].index)
    idx = flows["TYPE"] == "E"
    reactive_sources = gen_names + list(flows[idx].index)

    # C[k]: set of branches incident to bus k where each branch total reactive power flow direction is away from bus k
    # D[k]: set of branches incident to bus k where each branch total reactive power flow direction is towards bus k
    C, D = {}, {}
    for k in bus_names:
        C[k] = {}
        D[k] = {}

    # make directed graph as a convenient way to explore bus-branch connections
    # reset index so that ID becomes column that will be used in edge_attr
    G = nx.from_pandas_edgelist(
        flows.reset_index(), source="F_BUS", target="T_BUS", edge_attr=True, create_using=nx.MultiDiGraph)

    # find branches incident to each bus and categorize them into sets C and D
    # save the reactive power Q that leaves or enters bus k, branch side QF
    # "from" or QT "to" and the multiple sign (1 or -1) that is used to obtain
    # the correct power flow according to definition of sets C and D
    for k in bus_names:
        # branches that are oriented "to" buses
        for (_, _, data) in G.in_edges(k, data=True):
            ID, QT = data["ID"], data["QT"]
            # branch power is QT taken as is if positive or multiplied by -1 if negative (branch reversal)
            if QT > 0:
                D[k][ID] = {"Q": QT, "side": "QT", "sign": 1}
            else:
                C[k][ID] = {"Q": -QT, "side": "QT", "sign": -1}
        # branches that are oriented "from" buses
        for (_, _, data) in G.out_edges(k, data=True):
            ID, QF = data["ID"], data["QF"]
            # branch power is QF taken as is if positive or multiplied by -1 if negative (branch reversal)
            if QF > 0:
                C[k][ID] = {"Q": QF, "side": "QF", "sign": 1}
            else:
                D[k][ID] = {"Q": -QF, "side": "QF", "sign": -1}

    logging.debug(f"{15*'='} Sets C and D {15*'='}")
    for k in bus_names:
        logging.debug(f"C[{k}] = {C[k]}")
        logging.debug(f"D[{k}] = {D[k]}\n")
    logging.debug(55*"=")

    # Components of reactive power flows
    # Each branch is listed twice with suffixes "_QF" and "_QT" for both ends
    branch_names = [
        i + suffix for i in flows.index for suffix in ["_QF", "_QT"]]
    # List "from" and "to" buses twice
    F_BUS = [i for i in flows["F_BUS"] for _ in range(2)]
    T_BUS = [i for i in flows["T_BUS"] for _ in range(2)]
    TYPE = [i for i in flows["TYPE"] for _ in range(2)]
    TAP = [i for i in flows["TAP"] for _ in range(2)]
    # Interleave QF and QT
    QFT = [val for pair in zip(flows["QF"], flows["QT"]) for val in pair]
    # numpy array that will hold the results
    nr, nc = len(branch_names), len(reactive_sources)
    RPFC = np.zeros((nr, nc))
    # Dictionaries that hold row and column indices for each element in branch_names and gen_names
    row_i = {v: k for v, k in zip(branch_names, range(nr))}
    col_i = {v: k for v, k in zip(reactive_sources, range(nc))}
    # List that indicates whether APFC have been calculated
    DONE = nr * [False]

    # Generator branches are supplied by respective generators only,
    # therefore we know their reactive power flow components (RPFC).
    # They are starting points in the procedure.
    idx = flows["TYPE"] == "G"
    QG = flows.loc[idx, "QF"]
    ig = flows[idx].index
    for g, Q in zip(ig, QG):
        i, j = [row_i[g + suffix] for suffix in ["_QF", "_QT"]]
        gg = col_i[g]
        RPFC[i, gg] = Q
        RPFC[j, gg] = Q
        DONE[i] = True
        DONE[j] = True
    known_components = list(gen_names)
    logging.debug(f"Initial known components: {gen_names}")

    # Branches that are inserting reactive powers at both ends
    # their RPFC are from themselfs only
    idx = (flows["TYPE"] == "E") & (flows["QF"] < 0) & (flows["QT"] > 0)
    ID = flows.loc[idx].index
    for i in flows.loc[idx].index:
        known_components.append(i)
        f = flows.loc[i, "F_BUS"]
        t = flows.loc[i, "T_BUS"]
        logging.debug(f"Initial known components: {i} ({f}-{t}) (inserts reactive power at both ends)")
        for c in ["QF", "QT"]:
            Q = flows.loc[i, c]
            ii = row_i[expand(i, c)]
            RPFC[ii, col_i[i]] = Q
            DONE[ii] = True

    # Calculate RPFC of other branches
    processed_buses = []
    repeat = True
    iter = 0
    start = time.time()
    while repeat and iter < itermax:
        iter += 1
        logging.debug(f"{10*'='} Start of iteration {iter} {10*'='}")
        iter_start = time.time()
        for k in bus_names:
            # Find bus k for which RPFC are known for all branches in D[k]
            if k not in processed_buses:
                logging.debug(f"Processing bus k = {k}")
                check = all(b in known_components for b in D[k].keys())
                if check:
                    # If all RPFC for branches in D[k] are known we may calculate
                    # the RPFC for the branches in C[k].
                    branches_c, Q_c, side_c, sign_c = unpack(C[k], "Q")
                    branches_d, Q_d, side_d, sign_d = unpack(D[k], "Q")
                    logging.debug(f"D[{k}]: {branches_d}, {Q_d}, {side_d}, {sign_d}, sum(Q) = {sum(Q_d)}")
                    logging.debug(f"C[{k}]: {branches_c}, {Q_c}, {side_c}, {sign_c}, sum(Q) = {sum(Q_c)}")
                    # Loop through all branches in C[k] and calculate their RPFC
                    for i, Q, side_i, sign_i in zip(branches_c, Q_c, side_c, sign_c):
                        ii = row_i[expand(i, side_i)]
                        logging.debug(f"outflow: {i}_{side_i} ({F_BUS[ii]}-{T_BUS[ii]})")
                        # autopep8: off
                        # if abs(Q) <= 1e-4:
                            # DONE[row_i[expand(i, "QF")]] = True
                            # DONE[row_i[expand(i, "QT")]] = True
                            # known_components.append(i)
                            # continue
                        # f = Q/sum(Q_d)
                        if abs(Q) <= 1e-4 and sum(Q_d) == 0:
                            f = 0
                            DONE[row_i[expand(i, "QF")]] = True
                            DONE[row_i[expand(i, "QT")]] = True
                        else:
                            f = Q/sum(Q_d)
                        # autopep8: on
                        logging.debug(f"f = {Q}/{sum(Q_d)} = {f}")
                        # Loop through all components of all branches in D[k]
                        for g in reactive_sources:
                            for j, side_j, sign_j in zip(branches_d, side_d, sign_d):
                                # RPFC for branch j (from D[k]) and generator g for the side (side_j) that brings power into bus k
                                jj = row_i[expand(j, side_j)]
                                gg = col_i[g]
                                Qjg = RPFC[jj, gg]
                                # RPFC for branch i (from C[k]) and generator g for the side (side_i) that takes out power from bus k
                                # for both branches we must multiple with the power flow orientation sign
                                CQ = sign_i * sign_j * f * Qjg
                                ii = row_i[expand(i, side_i)]
                                RPFC[ii, gg] += CQ
                                if CQ != 0:
                                    logging.debug(f"inflow: {j}.{g} = {Qjg}, CQ = {CQ} " +
                                                  f"--> component for: {expand(i, side_i)}.{g} ({F_BUS[ii]}-{T_BUS[ii]})")
                                DONE[ii] = True
                                # RPFC for the other side of branch i
                                ii = row_i[expand(i, other(side_i))]
                                # Q1 = flows.loc[i, side_i+"1"] # <-- slow
                                # Q1_other_side = flows.loc[i, other(side_i)+"1"] # <-- slow
                                Q1 = flows_mat[flows_row[i], flows_col[side_i+"1"]]
                                Q1_other_side = flows_mat[flows_row[i], flows_col[other(side_i)+"1"]]
                                CQ_other_side = sign_i * Q1_other_side/Q1 * CQ
                                RPFC[ii, gg] += CQ_other_side
                                if CQ_other_side != 0:
                                    logging.debug(
                                        f"CQ_other_side = {CQ_other_side} --> component for: {expand(i, other(side_i))}.{g}")
                                DONE[ii] = True
                        # Components due to shunt elements of branch i
                        if flows.loc[i, "TYPE"] == "E":
                            ii = row_i[expand(i, side_i)]
                            logging.debug(f"Components due to shunt elements of branch {i} ({F_BUS[ii]}-{T_BUS[ii]})")
                            shunt = side_i[0] + "S" + side_i[1]
                            shunt_other_side = side_i[0] + \
                                "S" + other(side_i)[1]
                            # QS = flows.loc[i, shunt] # <-- slow
                            # QS_other_side = flows.loc[i, shunt_other_side] # <-- slow
                            QS = flows_mat[flows_row[i], flows_col[shunt]]
                            QS_other_side = flows_mat[flows_row[i], flows_col[shunt_other_side]]
                            logging.debug(
                                f"QS = {QS} (side: {side_i}), QS_other_side = {QS_other_side} (side: {other(side_i)})")
                            ii = row_i[expand(i, other(side_i))]
                            gg = col_i[i]
                            # Q1 = flows.loc[i, side_i+"1"] # <-- slow
                            # Q1_other_side = flows.loc[i, other(side_i)+"1"] # <-- slow
                            Q1 = flows_mat[flows_row[i], flows_col[side_i+"1"]]
                            Q1_other_side = flows_mat[flows_row[i], flows_col[other(side_i)+"1"]]
                            logging.debug(
                                f"Q1 = {Q1} (side: {side_i}), Q1_other_side = {Q1_other_side} (side: {other(side_i)})")
                            if abs(Q1_other_side) < 1e-4 and abs(Q1) < 1e-4:
                                h = 0
                            else:
                                h = Q1_other_side/Q1
                            RPFC[ii, gg] = h * QS + QS_other_side
                            logging.debug(f"h = {Q1_other_side}/{Q1} = {Q1_other_side/Q1}")
                            logging.debug(
                                f"Qh = {h} * {QS} + {QS_other_side} = {h * QS + QS_other_side} --> component for: {expand(i, other(side_i))}.{i}")
                            DONE[ii] = True
                        else:
                            ii = row_i[expand(i, other(side_i))]
                            DONE[ii] = True
                        known_components.append(i)
                        for ik in known_components:
                            jk = row_i[ik + "_QF"]
                            logging.debug(f"    known components {ik} ({F_BUS[jk]}-{T_BUS[jk]})")
                    processed_buses.append(k)
                else:
                    logging.debug("    not all RPFC known, skipping ...")
        repeat = not all(DONE)
        logging.info("iter = {}, time = {:.3f}".format(iter, time.time()-iter_start))
        if graph:
            T = makeTable(RPFC, branch_names, reactive_sources,
                          F_BUS, T_BUS, TYPE, DONE, TAP)
            drawGraph(T, bus_pos, "iter = {}".format(iter))

    RPFC = makeTable(RPFC, branch_names, reactive_sources,
                     F_BUS, T_BUS, TYPE, DONE, TAP)
    RPFC["Total"] = RPFC[reactive_sources].sum(axis=1)
    RPFC["Error"] = RPFC["Total"] - QFT
    logging.info("total time = {:.3f}".format(time.time()-start))
    logging.info("RPFC: nr = {}, iters = {}, DONE = {}".format(nr, iter, not repeat))

    if repeat:
        for k in bus_names:
            if k not in processed_buses:
                logging.info("\nbus: {}".format(k))
                logging.info("IN")
                for b in D[k]:
                    f = flows.loc[b, "F_BUS"]
                    t = flows.loc[b, "T_BUS"]
                    ii = row_i[expand(b, D[k][b]["side"])]
                    logging.info("  branch: {} ({}-{}), {}, {}".format(b,
                                                                       f, t, D[k][b], DONE[ii]))
                logging.info("OUT")
                for b in C[k]:
                    f = flows.loc[b, "F_BUS"]
                    t = flows.loc[b, "T_BUS"]
                    ii = row_i[expand(b, C[k][b]["side"])]
                    logging.info("  branch: {} ({}-{}), {}, {}".format(b,
                                                                       f, t, C[k][b], DONE[ii]))

    logging.info("End of calculations\n")
    return flows, RPFC


def cleanRows(TABLE, power):
    idx = TABLE["TYPE"] == "E"
    E = TABLE[idx]
    LG = TABLE[~idx]
    F = "_" + power + "F"
    T = "_" + power + "T"
    LG.index = [i.replace(F, "").replace(T, "") for i in LG.index]
    LG = LG.drop_duplicates()
    TABLE = pd.concat([E, LG.drop_duplicates()])
    return TABLE


def components(case_file, itermax, graph):
    pd.set_option("display.precision", 3)
    graph = graph == 1

    logging.info(f"case_file = {case_file}")

    logging.info("Solve the power flow\n")
    bus, gen, branch, xy = runpf(case_file)
    logging.info(f"Total P gens   {sum(gen['PG']):9.3f} MW")
    logging.info(f"Total P load   {sum(bus['PD']):9.3f} MW")
    logging.info(f"Total P losses {sum(branch['PF']-branch['PT']):9.3f} MW\n")
    logging.info(f"Total Q gens   {sum(gen['QG']):9.3f} Mvar")
    logging.info(f"Total Q shunts {sum(bus['QS']):9.3f} Mvar")
    logging.info(f"Total Q load   {sum(bus['QD']):9.3f} Mvar")
    logging.info(f"Total Q losses {sum(branch['QF']-branch['QT']):9.3f} Mvar\n")

    logging.info("Read bus coordinates")
    if xy is not None:
        bus_pos = {}
        for i, x, y in xy:
            bus_pos[int(i)] = [x, y]
    else:
        bus_pos = None

    logging.info("Calculate active power flow components")
    flows_p, APFC = makeAPFC(bus, gen, branch, bus_pos, itermax, graph)
    APFC = cleanRows(APFC, "P")

    logging.info("Calculate reactive power flow components")
    flows_q, RPFC = makeRPFC(bus, gen, branch, bus_pos, itermax, graph)
    RPFC = cleanRows(RPFC, "Q")

    logging.info("Save the results")
    flows_p.to_csv(os.path.join("results", case_file + "-flows_p.csv"))
    APFC.to_csv(os.path.join("results", case_file + "-apfc.csv"))
    flows_q.to_csv(os.path.join("results", case_file + "-flows_q.csv"))
    RPFC.to_csv(os.path.join("results", case_file + "-rpfc.csv"))


if __name__ == "__main__":
    if not os.path.exists("results"):
        os.makedirs("results")

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="casefile", type=str,
                        default="case9", help="Matpower case file name")
    parser.add_argument("-i", dest="itermax", type=int,
                        default=20, help="Maximum number of iteration")
    parser.add_argument("-g", dest="graph", type=int,
                        default=0, help="Draw graph: 0 = no, 1 = yes")
    args = parser.parse_args()

    # Set loggers
    logging.getLogger().setLevel(logging.NOTSET)
    formater = logging.Formatter("%(funcName)s %(levelname)s: %(message)s")
    # Add stdout handler, with level INFO
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(formater)
    logging.getLogger().addHandler(console)
    # Add file handler, with level DEBUG
    file = logging.FileHandler(os.path.join("results", args.casefile + ".log"), mode="w")
    file.setLevel(logging.DEBUG)
    file.setFormatter(formater)
    logging.getLogger().addHandler(file)

    components(args.casefile, args.itermax, args.graph)
