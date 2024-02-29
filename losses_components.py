import argparse
import os
import pandas as pd


def gen_load(case_file):
    file_name = os.path.join("results", case_file + "-apfc.csv")
    APFC = pd.read_csv(file_name)
    APFC = APFC.rename(columns={"Unnamed: 0": "name"})

    gen_name = [c for c in APFC.columns if c.startswith("G")]

    idx = (APFC["TYPE"] == "E") & (APFC["name"].str.endswith("_PF"))
    APFC_PF = APFC.loc[idx, ["name", "F_BUS", "T_BUS"] + gen_name]
    APFC_PF["name"] = APFC_PF["name"].apply(lambda x: x.replace("_PF", ""))

    idx = (APFC["TYPE"] == "E") & (APFC["name"].str.endswith("_PT"))
    APFC_PT = APFC.loc[idx, ["name", "F_BUS", "T_BUS"] + gen_name]
    APFC_PT["name"] = APFC_PT["name"].apply(lambda x: x.replace("_PT", ""))

    LOSSES = pd.merge(APFC_PF, APFC_PT, how="left", on=["F_BUS", "T_BUS"])
    for g in gen_name:
        LOSSES[g] = (LOSSES[f"{g}_x"] - LOSSES[f"{g}_y"]) * 1000
    LOSSES = LOSSES.rename(columns={"name_x": "name"})
    LOSSES = LOSSES[["name", "F_BUS", "T_BUS"] + gen_name]
    LOSSES["total"] = LOSSES[gen_name].sum(axis=1)
    file_name = os.path.join("results", case_file + "-line-losses.csv")
    LOSSES.to_csv(file_name, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        dest="casefile",
        type=str,
        default="case33",
        help="Matpower case file name",
    )
    args = parser.parse_args()

    gen_load(args.casefile)
