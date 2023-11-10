#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import errno
from itertools import product
import os

import uproot
import utils

from utils import TEMPLATES
from config import Configurator



config = Configurator()
VARS = ["MH","MET","2DEllipses"]
CHANNELS = ["inv"]
REGIONS = list(utils.REGION_ID_MAP.keys())
ALL_PROCESSES = ["AtoZH"] + config.backgrounds

RATE_NP = {
    # experimental 
    "lumi_13TeV_2017": { p: 1.02 for p in ALL_PROCESSES },
    "Norm_bjet": {"DYJets_bjet": 1.4, "WJets_bjet": 1.4}, #not on light jets -> freely floating
    
    # theory AN2019_094_v12
    "QCDscale_ttbar": {"TTW": "0.836/1.255", "TTZ": "0.907/1.081"},  #not on TT -> freely floating
    "QCDscale_singletop": {"SingleTop": "0.979/1.031"},
    "QCDscale_VV": {"VV": 1.03},
    "pdf_gg": {"TTZ": 1.035},  #not on TT -> freely floating
    "pdf_qqbar": {"VV": 1.05, "TTW": 1.036}, 
    "pdf_qg": {"SingleTop": 1.028},
}

SHAPES_NP = {
    # b tagging efficiency 
    "CMS_btag_bc_YEAR": ALL_PROCESSES,
    "CMS_btag_light_YEAR": ALL_PROCESSES,


    "CMS_pileup_YEAR": ALL_PROCESSES,
    "CMS_pileupJetID_YEAR": ALL_PROCESSES,
    "CMS_toppt_a": ["TT"],
    "CMS_toppt_b": ["TT"],

    # factorization scale - shape effect only 
    "CMS_muf_AtoZH": ["AtoZH"],
    "CMS_muf_DYJets_bjet": ["DYJets_bjet"],
    "CMS_mur_DYJets_ljet": ["DYJets_ljet"],
    "CMS_muf_TT": ["TT"],
    "CMS_muf_TTZ": ["TTZ"],
    "CMS_muf_SingleTop": ["SingleTop"],
    "CMS_muf_TTW": ["TTW"],
    "CMS_muf_VV": ["VV"],
    "CMS_muf_WJets_bjet": ["WJets_bjet"],
    "CMS_muf_WJets_ljet": ["WJets_ljet"],
    "CMS_muf_QCD": ["QCD"],

    # renormilization - shape effect only
    "CMS_mur_AtoZH": ["AtoZH"],
    "CMS_mur_DYJets_bjet": ["DYJets_bjet"],
    "CMS_muf_DYJets_ljet": ["DYJets_ljet"],
    "CMS_mur_TT": ["TT"],
    "CMS_mur_TTZ": ["TTZ"],
    "CMS_mur_SingleTop": ["SingleTop"],
    "CMS_mur_TTW": ["TTW"],
    "CMS_mur_VV": ["VV"],
    "CMS_mur_WJets_bjet": ["WJets_bjet"],
    "CMS_mur_WJets_ljet": ["WJets_ljet"],
    "CMS_mur_QCD": ["QCD"],

    # radiation
    "fsr_2": ALL_PROCESSES,
    "isr_2": ALL_PROCESSES,

    # pdf uncertainty 
    "pdf_AtoZH": ["AtoZH"],
    "pdf_DYJets_bjet": ["DYJets_bjet"],
    "pdf_DYJets_ljet": ["DYJets_ljet"],
    "pdf_SingleTop": ["SingleTop"],
    "pdf_TT": ["TT"],
    "pdf_TTW": ["TTW"],
    "pdf_TTZ": ["TTZ"],
    "pdf_VV": ["VV"],
    "pdf_WJets_bjet": ["WJets_bjet"],
    "pdf_WJets_ljet": ["WJets_ljet"],

    # ?
    "CMS_vjets_EWK_d1K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],
    "CMS_vjets_EWK_d2K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],
    "CMS_vjets_EWK_d3K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],
    "CMS_vjets_QCD_NLO_d1K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],
    "CMS_vjets_QCD_NLO_d2K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],
    "CMS_vjets_QCD_NLO_d3K": ["DYJets_bjet","DYJets_ljet","WJets_bjet","WJets_ljet"],

    # jet energy correction and resolution
    "JER": ALL_PROCESSES,
    "JES": ALL_PROCESSES,
}

# Chars per column
N1 = 79
N2 = 8
N3 = max([len(x) for x in REGIONS]) + 10


class Datacard():

    def __init__(self, year, mass_point, svar, channel, region, args):
        self.year = year
        self.mass_point = mass_point
        self.svar = svar
        self.channel = channel
        self.region = region
        self.fname = f"{self.mass_point}_{self.svar}_{self.channel}_{self.region}"
        self.processes = None
        self.args = args
        try:
            self.filter_processes()
            self.generate_datacard()
        except FileNotFoundError as e:
            print(f"File {e.filename} not found! Will skip.")

    def filter_processes(self):
        root_fpath = f"{TEMPLATES}/{self.year}/{self.fname}.root"
        if not os.path.isfile(root_fpath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), root_fpath)
        with uproot.open(root_fpath) as f:
            keys = f.keys()
        self.processes = [x for x in config.backgrounds if any([(f"{x}" in k) for k in keys])]
        self.processes = ["AtoZH"] + self.processes

    def write_block_header(self, f, block_name: str):
        f.write(f"# {block_name.capitalize()}\n")
        f.write(N1 * "-" + "\n")

    def write_parameters(self, f):
        self.write_block_header(f, "parameters")
        f.write(f"imax 1\njmax {len(self.processes)-1}\nkmax *\n")
        f.write(f"shapes * {self.region} {self.fname}.root "
                f"{self.region}/$PROCESS {self.region}/$PROCESS_$SYSTEMATIC\n\n")

    def write_channels(self, f):
        self.write_block_header(f, "channels")
        f.write(f"bin          {self.region}\n")
        f.write("observation  -1\n\n")

    def pad(self, s, n_pad):
        s = str(s)
        n_pad = n_pad - len(s)
        return s + n_pad * ' '

    def write_processes(self, f):
        padded_processes = [self.pad(x, N3) for x in self.processes]
        id_offset = 0  # int(self.region != "SignalRegion")
        padded_ids = [self.pad(i + id_offset, N3) for i in range(len(self.processes))]
        self.write_block_header(f, "processes")
        f.write(self.pad("bin", N1 + N2) + len(self.processes) * self.pad(self.region, N3) + "\n")
        f.write(self.pad("process", N1 + N2) + "".join(padded_processes) + "\n")
        f.write(self.pad("process", N1 + N2) + "".join(padded_ids) + "\n")
        f.write(self.pad("rate", N1 + N2) + len(self.processes) * self.pad("-1", N3))
        f.write("\n\n")

    def shape_np_applicable(self, shape_np):
        vetos = {
            shape_np == "CMS_sfelec" and self.channel == "diMuon",
            shape_np == "CMS_sfmu" and self.channel == "diElectron",
            shape_np == "CMS_trigger_sf_ee" and self.channel == "diMuon",
            shape_np == "CMS_trigger_sf_mumu" and self.channel == "diElectron"
        }
        return not any(vetos)

    def write_lnN_systematics(self, f):
        self.write_block_header(f, "systematics")
        for nuisance in RATE_NP:
            f.write(self.pad(nuisance, N1) + self.pad("lnN", 8))
            for process in self.processes:
                if process in RATE_NP[nuisance]:
                    np_val = RATE_NP[nuisance][process]
                    if isinstance(np_val, dict):
                        np_val = np_val[self.year]
                    f.write(self.pad(np_val, N3))
                else:
                    f.write(self.pad("-", N3))
            f.write("\n")

    def write_shape_systematics(self, f):
        for shape_np, applicable_processes in SHAPES_NP.items():
            shape_np = shape_np.replace("YEAR", self.year)
            f.write(self.pad(shape_np, N1) + self.pad("shape", 8))
            for process in self.processes:
                is_applicable = (
                    self.shape_np_applicable(shape_np)
                    and process in applicable_processes
                )
                if is_applicable:
                    f.write(self.pad(1, N3))
                else:
                    f.write(self.pad("-", N3))
            f.write("\n")

    def write_auto_mc_stats(self, f):
        f.write(f"\n* autoMCStats 0")

    def write_scaling(self, f):
        f.write("\nrate_BR_vv rateParam * AtoZH 1") #0.2 for inclusive xsec in AZH -> ttZ
        f.write("\nnuisance edit freeze rate_BR_vv ifexists")
        f.write("\nlumiscale rateParam * * 1") #3.33 for full run 2
        f.write("\nnuisance edit freeze lumiscale ifexists")

    def write_bkg_rate_params(self, f):
        if "TT" in self.processes:
            f.write("rate_TT rateParam * TT 1 [-2,2]")
        if "DYJets_bjet" in self.processes:
            f.write("\nrate_Vbjet rateParam * DYJets_bjet 1 [-2,2]")
        if "DYJets_ljet" in self.processes:
            f.write("\nrate_Vljet rateParam * DYJets_ljet 1 [-2,2]")
        if "WJets_bjet" in self.processes:
            f.write("\nrate_Vbjet rateParam * WJets_bjet 1 [-2,2]")
        if "WJets_ljet" in self.processes:
            f.write("\nrate_Vljet rateParam * WJets_ljet 1 [-2,2]")
        

    def generate_datacard(self):
        print(f"{self.year}/{self.fname}.dat")
        with open(f"{TEMPLATES}/{self.year}/{self.fname}.dat", 'w') as f:
            self.write_parameters(f)
            self.write_channels(f)
            self.write_processes(f)
            self.write_lnN_systematics(f)
            self.write_shape_systematics(f)
            self.write_bkg_rate_params(f)
            self.write_scaling(f)
            self.write_auto_mc_stats(f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=str, help="all/UL16/UL17/UL18/combined")
    parser.add_argument("--nojes", action="store_true")
    parser.add_argument("--inclusivejes", action="store_true")
    parser.add_argument("--nosamplevars", action="store_true")
    parser.add_argument("--autoMCStats", type=int, default=-1,
                        help="an integer for autoMCStats")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()

    # Years
    if args.year in ["UL16", "UL17", "UL18"]:
        YEARS = [args.year]
    if args.year == "all":
        YEARS = ["UL16", "UL17", "UL18"]

    ll = [VARS, YEARS, config.signals, CHANNELS, REGIONS]
    i = 0
    for svar, year, mass_point, channel, region in product(*ll):
        if not utils.is_valid_set(channel, region, svar):
            continue
        if (region == "CRDiffLeptonFlavours") != (channel == "ElectronMuon"):
            continue
        datacard = Datacard(year, mass_point, svar, channel, region, args)