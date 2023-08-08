#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import errno
from itertools import product
import os

import uproot

from config import Configurator
import utils


config = Configurator()
VARS = ["MTA","MTH","MH","MET","2DEllipses"]
CHANNELS = ["inv"]
REGIONS = list(utils.REGION_ID_MAP.keys())
ALL_PROCESSES = ["AtoZH"] + config.backgrounds

RATE_NP = {
    "lumi_13TeV_2016_uncorrelated": {
        p: {"UL16": 1.01, "UL17": '-', "UL18": '-'} for p in ALL_PROCESSES
    },
    "lumi_13TeV_2017_uncorrelated": {
        p: {"UL16": '-', "UL17": 1.02, "UL18": '-'} for p in ALL_PROCESSES
    },
    "lumi_13TeV_2018_uncorrelated": {
        p: {"UL16": '-', "UL17": '-', "UL18": 1.015} for p in ALL_PROCESSES
    },
    "lumi_13TeV_correlated_16_17_18": {
        p: {"UL16": 1.006, "UL17": 1.009, "UL18": 1.02} for p in ALL_PROCESSES
    },
    "lumi_13TeV_correlated_17_18": {
        p: {"UL16": '-', "UL17": 1.006, "UL18": 1.002} for p in ALL_PROCESSES
    },
    # qcd scales: AN2019_094_v12
    "QCDscale_ttbar": {"TTW": "0.836/1.255", "TTZ": "0.907/1.081"},  # "TT": "0.965/1.024",
    "QCDscale_singletop": {"SingleTop": "0.979/1.031"},
    "QCDscale_VV": {"VV": 1.03},
    "QCDscale_V": {"WJets": 1.038},  # "DYJets": 1.02},
    # pdf scales: AN2019_094_v12
    "pdf_gg": {"TTZ": 1.035},  # "TT": 1.042,
    "pdf_qqbar": {"WJets": "1.008/0.996", "VV": 1.05, "TTW": 1.036},  # "DYJets": 1.02,
    "pdf_qg": {"SingleTop": 1.028},
    # hdamp
    #"hdamp_5j": {"TT": {"UL16": "1.061/0.929", "UL17": "1.063/0.928", "UL18": "1.064/0.924"}},
    #"hdamp_6j": {"TT": {"UL16": "1.085/0.926", "UL17": "1.068/0.930", "UL18": "1.087/0.914"}},
    # TuneCP5
    #"TuneCP5_0b": {"TT": {"UL16": "1.045/0.974", "UL17": "1.049/0.974", "UL18": "0.989/0.969"}},
    #"TuneCP5_1b": {"TT": {"UL16": "0.988/1.005", "UL17": "0.999/0.999", "UL18": "1.003/0.990"}},
    #"TuneCP5_2b": {"TT": {"UL16": "0.999/0.998", "UL17": "1.003/0.990", "UL18": "1.009/1.005"}},
}

SHAPES_NP = {
    "CMS_btag_bc_YEAR": ALL_PROCESSES,
    "CMS_btag_light_YEAR": ALL_PROCESSES,
    "CMS_pileup_YEAR": ALL_PROCESSES,
    "CMS_pileupJetID_YEAR": ALL_PROCESSES,
    #"CMS_trigger_sf_ee": ALL_PROCESSES,
    #"CMS_trigger_sf_mumu": ALL_PROCESSES,
    #"CMS_sfelec_reco": ALL_PROCESSES,
    #"CMS_sfelec_tight_id": ALL_PROCESSES,
    #"CMS_sfmu_isolation": ALL_PROCESSES,
    #"CMS_sfmu_tight_id": ALL_PROCESSES,
    "CMS_toppt_a": ["TT"], #or ["AtoZH", "TT", "TTW", "TTZ"],?
    "CMS_toppt_b": ["TT"], #or ["AtoZH", "TT", "TTW", "TTZ"],?
    #"CMS_scale_j_Absolute": ALL_PROCESSES,
    #"CMS_scale_j_BBEC1": ALL_PROCESSES,
    #"CMS_scale_j_EC2": ALL_PROCESSES,
    #"CMS_scale_j_FlavorQCD": ALL_PROCESSES,
    #"CMS_scale_j_HF": ALL_PROCESSES,
    #"CMS_scale_j_RelativeBal": ALL_PROCESSES,
    #"CMS_scale_j_Absolute_YEAR": ALL_PROCESSES,
    #"CMS_scale_j_BBEC1_YEAR": ALL_PROCESSES,
    #"CMS_scale_j_EC2_YEAR": ALL_PROCESSES,
    #"CMS_scale_j_HF_YEAR": ALL_PROCESSES,
    #"CMS_scale_j_RelativeSample_YEAR": ALL_PROCESSES,
    #"CMS_scale_j_Total": ALL_PROCESSES,
    #"CMS_res_j_YEAR": ALL_PROCESSES,
    "CMS_mur_AtoZH": ["AtoZH"],
    "CMS_muf_AtoZH": ["AtoZH"],
    "CMS_mur_ZJets": ["ZJets"],
    "CMS_muf_ZJets": ["ZJets"],
    "CMS_mur_TT": ["TT"],
    "CMS_muf_TT": ["TT"],
    "CMS_mur_TTZ": ["TTZ"],
    "CMS_muf_TTZ": ["TTZ"],
    "CMS_mur_SingleTop": ["SingleTop"],
    "CMS_muf_SingleTop": ["SingleTop"],
    "CMS_mur_TTW": ["TTW"],
    "CMS_muf_TTW": ["TTW"],
    "CMS_mur_VV": ["VV"],
    "CMS_muf_VV": ["VV"],
    "CMS_mur_WJets": ["WJets"],
    "CMS_muf_WJets": ["WJets"],
    "CMS_mur_QCD": ["QCD"],
    "CMS_muf_QCD": ["QCD"],
    "fsr_2": ALL_PROCESSES,
    "isr_2": ALL_PROCESSES,
    "pdf_AtoZH": ["AtoZH"],
    "pdf_ZJets": ["ZJets"],
    "pdf_SingleTop": ["SingleTop"],
    "pdf_TT": ["TT"],
    "pdf_TTW": ["TTW"],
    "pdf_TTZ": ["TTZ"],
    "pdf_VV": ["VV"],
    "pdf_WJets": ["WJets"],
    "CMS_vjets_EWK_d1K": ["ZJets","WJets"],
    "CMS_vjets_EWK_d2K": ["ZJets","WJets"],
    "CMS_vjets_EWK_d3K": ["ZJets","WJets"],
    "CMS_vjets_QCD_NLO_d1K": ["ZJets","WJets"],
    "CMS_vjets_QCD_NLO_d2K": ["ZJets","WJets"],
    "CMS_vjets_QCD_NLO_d3K": ["ZJets","WJets"],
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
        print(svar)
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
        combine_fpath = f"{self.year}/{self.fname}.root"
        if not os.path.isfile(combine_fpath):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), combine_fpath)
        with uproot.open(combine_fpath) as f:
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
        if self.args.autoMCStats > -1:
            N = self.args.autoMCStats
            f.write(f"* autoMCStats {N}\n")

    def add_bkg_rate_params(self, f, tt: bool, vj: bool):
        if tt and "TT" in self.processes:
            f.write(f"rate_ttbar rateParam * TT 1 [-2,2]")
        if vj and "ZJets" in self.processes:
            f.write(f"\nrate_vjets rateParam * ZJets 1 [-2,2]")
        if vj and "WJets" in self.processes:
            f.write(f"\nrate_vjets rateParam * WJets 1 [-2,2]")

    def generate_datacard(self):
        print(f"{self.year}/{self.fname}.dat")
        with open(f"{self.year}/{self.fname}.dat", 'w') as f:
            self.write_parameters(f)
            self.write_channels(f)
            self.write_processes(f)
            self.write_lnN_systematics(f)
            self.write_shape_systematics(f)
            self.write_auto_mc_stats(f)
            self.add_bkg_rate_params(f, True, True)


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

