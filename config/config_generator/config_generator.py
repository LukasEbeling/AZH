#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import collections.abc
import glob
import itertools
import os
import sys
import yaml
import re

from lxml import etree
import uproot


CMSSW_BASE = os.environ.get("CMSSW_BASE")
UHH_DATASETS_PATH = os.path.join(CMSSW_BASE, "src/UHH2/common/UHH2-datasets")
sys.path.append(UHH_DATASETS_PATH)
from CrossSectionHelper import MCSampleValuesHelper


def parse_cmd_arguments():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "analysis_step",
        help="Preselection/Reconstruction"
    )
    parser.add_argument(
        "year",
        help="UL16preVFP/UL16postVFP/UL17/UL18"
    )
    group.add_argument(
        "--data",
        help="Set this flag to generate a config for data",
        action="store_true"
    )
    return parser.parse_args()


class Config():

    def __init__(self):
        self.element_tree = None

        # XML Components
        self.xml_doctype = None
        self.xml_config = None
        self.xml_job_config = None
        self.xml_cycle = None

    def set_element_tree(self):
        self.element_tree = etree.ElementTree(self.xml_job_config)


class ConfigFactory():

    def __init__(self, args, var):
        self.config = Config()
        self.args = args
        self.year = args.year
        self.step = args.analysis_step
        self.is_data = args.data
        self.is_uhh_dataset = (self.step == "Preselection")
        self.data_type = "DATA" if self.is_data else "MC"
        self.mc_helper = MCSampleValuesHelper()
        self.jes_variation = var["jes"]
        self.jer_variation = var["jer"]
        self.params = None
        self.samples = None
        self.output_fname = None
        self.target_lumis = None
        self.output_dirs = None
        self.n_events_break = None
        self.file_split = None

        # Fill Fields
        self.load_params()
        self.load_samples()
        self.set_output_fname()

        print(
            self.data_type,
            self.step,
            "jes" + self.jes_variation.title(),
            "jer" + self.jer_variation.title(),
            self.year,
        )

    @property
    def param_file_dir(self): return "param_files"

    def generate_config_elements(self):
        self.config.xml_doctype = self.generate_doctype()
        self.build_xml_structure()
        self.generate_input_data_blocks()
        self.generate_user_config()
        self.config.set_element_tree()
        self.save_xml_file()
        self.postprocess_entities()

    def set_output_fname(self):
        jet_var_suffix = (f"_JES{self.jes_variation}"
                          f"_JER{self.jer_variation}")
        self.output_fname = (
            f"{self.params['analysis_name']}{self.step}"
            f"_{self.data_type}"
            f"_{self.year}"
            f"{jet_var_suffix}.xml"
        )

    def build_xml_structure(self):
        job_config_attrs = {
            "JobName": self.step + self.data_type + self.year,
            "OutputLevel": "Info"
        }

        step_dir = self.output_dirs[self.step]
        output_dir = os.path.join(
            CMSSW_BASE,
            "src/UHH2/2HDM/data",
            step_dir,
            self.data_type,
            self.year
        ) + '/'
        cycle_attrs = {
            'Name': self.params['cycle']['name'],
            'OutputDirectory': output_dir,
            'PostFix': "",
            'TargetLumi': self.target_lumis[self.year]
        }
        lib_element = {'Name': 'libSUHH2AtoZH'}
        pkg_element = {'Name': 'SUHH2AtoZH.par'}
        comment = self.generate_config_comment()

        self.config.xml_job_config = etree.Element("JobConfiguration", job_config_attrs)
        self.config.xml_job_config.addprevious(comment)
        etree.SubElement(self.config.xml_job_config, "Library", lib_element)
        etree.SubElement(self.config.xml_job_config, "Package", pkg_element)
        self.config.xml_cycle = etree.SubElement(self.config.xml_job_config, "Cycle", cycle_attrs)

    @property
    def input_path(self):
        if self.step == "Preselection":
            path_common = "/nfs/dust/cms/user/hundhad/AtoZH_Common/signal_samples/"
            return os.path.join(path_common, self.year)

        # else the step is "Reconstruction"
        input_path = os.path.join(
            CMSSW_BASE,
            "src/UHH2/2HDM/data",
            "output_01_preselection",
            self.data_type,
            self.year
        )

        return input_path

    def load_params(self):
        with open(f"{self.param_file_dir}/params_common.yaml", 'r') as f:
            self.params = yaml.safe_load(f)
        with open(f"{self.param_file_dir}/job_splits.yaml", 'r') as f:
            p = yaml.safe_load(f)

        self.n_events_break = p[self.step][self.data_type]["n_events_break"]
        self.file_split = p[self.step][self.data_type]["file_split"]
        self.target_lumis = self.params["target_lumis"]
        self.output_dirs = self.params["cycle"]["output_subdirs"]

    def resolve_sample_variations(self, x):
        if isinstance(x, dict):
                return [v for v in x.values()]
        else:
            return x

    def flatten_irregular_list(self, x):
        if isinstance(x, list):
            return [a for i in x for a in self.flatten_irregular_list(i)]
        else:
            return [x]

    def filter_samples_for_combined_variations(self, samples):
        """
        Filters out combinations like hdampUP + JESdown.
        """
        if (self.jes_variation == self.jer_variation):
            return samples
        samples = [s for s in samples if not re.search(r"(hdamp|TuneCP5)(UP|DOWN|up|down)", s)]
        return samples

    def load_uhh2_dataset_samples(self):
        with open("../signals.txt", 'r') as f:
            signal_samples = ["AToZHToLLTTbar_" + x.removesuffix("\n") for x in f]

        with open(f"{self.param_file_dir}/uhh_samples_{self.data_type.lower()}.yaml", 'r') as f:
            samples_by_process = yaml.safe_load(f)

        samples = [list(map(self.resolve_sample_variations, y)) for y in samples_by_process.values()]
        samples = self.flatten_irregular_list(samples)
        #samples += signal_samples

        if self.is_data:
            samples = []
            for x in samples_by_process.values():
                for sample, years in x.items():
                    if self.year in years:
                        samples.append(sample)

        samples = self.filter_samples_for_combined_variations(samples)

        return samples

    def load_root_samples(self):
        def root_file_not_empty(fpath):
            if os.path.getsize(fpath) > 50000:
                return True
            with uproot.open(fpath) as f:
                return any(["AnalysisTree" in x for x in f.keys()])

        if self.step != "Reconstruction" or self.is_data:
            suffix = ""
        else:
            suffix = (f"_JES{self.jes_variation}"
                      f"_JER{self.jer_variation}")
        path = os.path.join(
            self.input_path,
            f"*{suffix}_{self.year}*.root"
        )
        print(path)
        samples = glob.glob(path)
        samples = list(filter(root_file_not_empty, samples))
        samples = [s.split('/')[-1] for s in samples]
        samples = self.filter_samples_for_combined_variations(samples)
        return samples

    def load_samples(self):
        if self.is_uhh_dataset:
            self.samples = self.load_uhh2_dataset_samples()
        else:
            self.samples = self.load_root_samples()

    def get_bare_sample_name(self, sample_fname):
        bare_sample_name = sample_fname.removeprefix("uhh2.AnalysisModuleRunner.MC.")
        bare_sample_name = bare_sample_name.removeprefix("uhh2.AnalysisModuleRunner.DATA.")
        bare_sample_name = bare_sample_name.removesuffix(".root")
        bare_sample_name = bare_sample_name.removesuffix(".xml")
        bare_sample_name = bare_sample_name.removesuffix(f"_{self.year}")
        bare_sample_name = bare_sample_name.removesuffix(f"_JES{self.jes_variation}_JER{self.jer_variation}")
        return bare_sample_name

    def get_lumi_from_sample_name(self, sample_name):
        if self.is_data:
            return "1"

        with open(f"{self.param_file_dir}/private_sample_lumis.yaml") as f:
            p = yaml.safe_load(f)[self.year]
        try:
            return p[f"{sample_name}_{self.year}"]
        except KeyError:
            return str(self.mc_helper.get_lumi(sample_name, "13TeV", self.year))

    def get_version_from_sample_name(self, sample_name):
        sample_version = sample_name
        sample_version = sample_name
        sample_version += f"_JES{self.jes_variation}_JER{self.jer_variation}"
        sample_version += f"_{self.year}"
        return sample_version

    def generate_root_input_blocks(self):
        for sample in self.samples:
            bare_sample_name = self.get_bare_sample_name(sample)
            sample_version = self.get_version_from_sample_name(bare_sample_name)
            lumi = self.get_lumi_from_sample_name(re.sub(r"_\d+$", '', bare_sample_name))
            attrs = {
                'Lumi': lumi,
                'NEventsMax': "-1",
                'Type': self.data_type,
                'Version': sample_version,
                'Cacheable': "False"
            }
            in_attrs = {
                'FileName': os.path.join(self.input_path, sample),
                'Lumi': '0.0'
            }

            input_data = etree.SubElement(self.config.xml_cycle, "InputData", attrs)
            etree.SubElement(input_data, "In", in_attrs)
            etree.SubElement(input_data, "InputTree", {"Name": "AnalysisTree"})
            etree.SubElement(input_data, "OutputTree", {"Name": "AnalysisTree"})

    def generate_uhh_dataset_input_blocks(self):
        for sample in self.samples:
            bare_sample_name = self.get_bare_sample_name(sample)
            sample_version = self.get_version_from_sample_name(bare_sample_name)
            lumi = "1"
            if not self.is_data:
                lumi = str(self.mc_helper.get_lumi(sample, "13TeV", self.year))
            attrs = {
                "Lumi": lumi,
                "NEventsMax": "-1",
                "Type": self.data_type,
                "Version": sample_version,
                "Cacheable": "False"
            }
            entity_ref = etree.fromstring(f"<temp>{sample}</temp>")

            input_data = etree.SubElement(self.config.xml_cycle, "InputData", attrs)
            input_data.append(entity_ref)
            etree.SubElement(input_data, "InputTree", {"Name": "AnalysisTree"})
            etree.SubElement(input_data, "OutputTree", {"Name": "AnalysisTree"})

    def generate_input_data_blocks(self):
        if self.is_uhh_dataset:
            return self.generate_uhh_dataset_input_blocks()
        else:
            return self.generate_root_input_blocks()

    def set_custom_user_config_values(self, name, value):
        btag_sf_filenames = {
            "UL16preVFP": "wp_deepJet_106XUL16preVFP_v2.csv",
            "UL16postVFP": "wp_deepJet_106XUL16postVFP_v3.csv",
            "UL17": "wp_deepJet_106XUL17_v3.csv",
            "UL18": "wp_deepJet_106XUL18_v2.csv"
        }
        lumi_filenames = {
            "UL16preVFP": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root",
            "UL16postVFP": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root",
            "UL17": "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root",
            "UL18": "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root"
        }

        uhh_data_path = "src/UHH2/common/UHH2-data/"
        if name == "additionalBranches":
            value = ' '.join(value)
        if name == "jecsmear_direction":
            value = self.jes_variation
        if name == "jersmear_direction":
            value = self.jer_variation
        if name == "BTagMCEffFile":
            mc_path = "src/UHH2/2HDM/data/output_01_preselection/MC/"
            fname = "BTagMCEfficiencyHists.root"
            value = os.path.join(CMSSW_BASE, mc_path, self.year, fname)
        if name == "BTagScaleFactorCSV":
            uhh_btag_path = uhh_data_path + "btagging_SFs_UL/"
            fname = btag_sf_filenames[self.year]
            value = os.path.join(CMSSW_BASE, uhh_btag_path, fname)
        if name == "lumi_file":
            fname = lumi_filenames[self.year]
            value = os.path.join(CMSSW_BASE, uhh_data_path, self.year, fname)

        return value

    @property
    def analysis_module(self): return f"InvAtoZH{self.step}"

    def generate_user_config(self):
        user_config = etree.SubElement(self.config.xml_cycle, 'UserConfig')
        for name, value in self.params['user_config'].items():
            value = self.set_custom_user_config_values(name, value)
            etree.SubElement(user_config, "Item", {'Name': name, 'Value': value})

        pu_base_dir = os.path.join(CMSSW_BASE, "src/UHH2/common/UHH2-data/")
        mc_path = pu_base_dir + f"{self.year}/MyMCPileupHistogram_{self.year}.root"
        data_path = pu_base_dir + f"{self.year}/MyDataPileupHistogram_{self.year}.root"
        data_path_down = pu_base_dir + f"{self.year}/MyDataPileupHistogram_{self.year}_66017.root"
        data_path_up = pu_base_dir + f"{self.year}/MyDataPileupHistogram_{self.year}_72383.root"

        etree.SubElement(user_config, "Item", {'Name': 'AnalysisModule', 'Value': self.analysis_module})
        etree.SubElement(user_config, "Item", {'Name': 'pileup_directory', 'Value': mc_path})
        etree.SubElement(user_config, "Item", {'Name': 'pileup_directory_data', 'Value': data_path})
        etree.SubElement(user_config, "Item", {'Name': 'pileup_directory_data_down', 'Value': data_path_down})
        etree.SubElement(user_config, "Item", {'Name': 'pileup_directory_data_up', 'Value': data_path_up})

    def generate_config_comment(self) -> str:
        config_parse = etree.Element(
            "ConfigParse",
            {
                "NEventsBreak": self.n_events_break,
                "LastBreak": "0",
                "FileSplit": self.file_split,
                "AutoResubmit": "0"
            }
        )

        mail = self.params["config_parse"]["mail"]
        notification = self.params["config_parse"]["notification"]
        workdir = (f"workdir_{self.step}"
                   f"_{self.data_type}"
                   f"_JES{self.jes_variation}"
                   f"_JER{self.jer_variation}"
                   f"_{self.year}")
        config_sge = etree.Element(
            "ConfigSGE",
            {
                "RAM": "2",
                "DISK": "2",
                "Mail": mail,
                "Notification": notification,
                "Workdir": workdir
            }
        )

        comment = (f'\n  {etree.tostring(config_parse, encoding="unicode")}'
                   f'\n  {etree.tostring(config_sge, encoding="unicode")}\n')

        return etree.Comment(comment)

    def generate_uhh_dataset_entities(self):
        entities = ''

        for uhh_sample in self.samples:
            sample_path = self.mc_helper.get_xml(uhh_sample, "13TeV", str(self.year))
            absolut_xml_path = os.path.join(UHH_DATASETS_PATH, sample_path)
            entities += f'  <!ENTITY {uhh_sample} SYSTEM "{absolut_xml_path}">\n'

        return entities

    def generate_entities(self):
        if self.is_uhh_dataset:
            entities = self.generate_uhh_dataset_entities()
        else:
            entities = ""

        return entities

    def generate_doctype(self):
        entities = self.generate_entities()
        doctype = (f'<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd"[\n'
                   f'{entities}'
                   f']>')
        return doctype

    def save_xml_file(self):
        self.config.element_tree.write(
            self.output_fname,
            xml_declaration=True,
            encoding='UTF-8',
            pretty_print=True,
            doctype=self.config.xml_doctype
        )

    def postprocess_entities(self):
        with open(self.output_fname, 'r') as f:
            _f = f.read()
        with open(self.output_fname, 'w') as f:
            _f = _f.replace('<temp>', '&')
            _f = _f.replace('</temp>', ';')
            f.write(_f)


if __name__ == "__main__":
    args = parse_cmd_arguments()
    # Build list of the relevant JES and JER variations
    variations = [
        {"jes": "", "jer": ""},
    ]

    if args.data:
        variations = [
            {"jes": "nominal", "jer": "nominal"},
        ]
    else:
        variations = ["nominal", "up", "down"]
        variations = [
            {"jes": v1, "jer": v2}
            for v1, v2 in itertools.product(variations, repeat=2)
            if not (v1 != "nominal" and v2 != "nominal")
        ]

    # Build config files
    for var in variations:
        factory = ConfigFactory(args, var)
        factory.generate_config_elements()
