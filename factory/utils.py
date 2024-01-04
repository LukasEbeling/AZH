import copy
import os

import matplotlib.patches as mpatches
import numpy as np
import scipy.integrate as integrate


CMSSW_BASE = os.environ.get('CMSSW_BASE')
CACHE = os.path.join(CMSSW_BASE,"src/UHH2/AZH/data/output_03_templates/cache")
TEMPLATES = os.path.join(CMSSW_BASE,"src/UHH2/AZH/data/output_03_templates")

ELLIPSE_X = "MET"
ELLIPSE_Y = "m_H"

REGION_ID_MAP = {
    "SR_6J": 1,
    "SR_5J": 2,
    "SR_1B_5J": 3,
    "SR_1B_6J": 4,
    "IR_0B_5J": 10,
    "IR_0B_6J": 11,
    "LR_2B_5J": 12,
    "LR_2B_6J": 13,
    "LR_1B_5J": 14,
    "LR_1B_6J": 15,
    "LR_0B_5J": 16,
    "LR_0B_6J": 17,
    #"SR_DNN": 100,
    #"TT_DNN": 101,
    #"WJ_DNN": 102,
    #"DY_DNN": 103,
    #"QCD_DNN": 104,
}

OUTPUT_MAP = {
    "z_pt_reco": "ZPT",
    "ellipses": "2DEllipses",
    "jetsAk4CHS.m_phi_1": "Jet1Phi",
    "jetsAk4CHS.m_eta_1": "Jet1Eta",
    "jetsAk4CHS.m_pt_1": "Jet1Pt",
    "jetsAk4CHS.m_pt_2": "Jet2Pt",
    "jetsAk4CHS.m_pt_3": "Jet3Pt",
    "jetsAk4CHS.m_pt_4": "Jet4Pt",
    "slimmedElectronsUSER.m_pt_1": "Elec1Pt",
    "slimmedElectronsUSER.m_pt_2": "Elec2Pt",
    "slimmedMuonsUSER.m_pt_1": "Muon1Pt",
    "slimmedMuonsUSER.m_pt_2": "Muon2Pt",
    "z_mass_reco": "ZMass",
    "mt_A": "MTA",
    "mt_H": "MTH",
    "m_H": "MH",
    "MET": "MET",
    "HT":"HT",
}

class Ellipse():

    def __init__(self, mean_x, mean_y, width, height, angle, n_std):
        self.mean_x = mean_x
        self.mean_y = mean_y
        self.width = width
        self.height = height
        self.angle = angle
        self.rad_angle = np.radians(angle)
        self.semi_w = width / 2
        self.semi_h = height / 2
        self.n_std = n_std

    def set_width(self, value):
        self.width = value
        self.semi_w = value / 2

    def set_height(self, value):
        self.height = value
        self.semi_h = value / 2

    def rescale_axes(self, factor):
        width = self.width * factor
        height = self.height * factor
        self.set_width(width)
        self.set_height(height)

    def get_plt_patch(self, color):
        mean = (self.mean_x, self.mean_y)
        ell = mpatches.Ellipse(xy=mean, width=self.width, height=self.height,
                               angle=self.angle, edgecolor=color, fc='None',
                               lw=2, zorder=4)
        return ell

    def pct_points_in(self, x_vals, y_vals, weights):
        p_in = 0
        norm = np.sum(weights)
        for x, y, w in zip(x_vals, y_vals, weights):
            if self.is_point_included(x, y):
                p_in += w / norm
        return p_in

    def is_point_included(self, x, y):
        cos = np.cos(self.rad_angle)
        sin = np.sin(self.rad_angle)
        n1 = (cos * (x - self.mean_x) + sin * (y - self.mean_y)) ** 2
        n2 = (sin * (x - self.mean_x) - cos * (y - self.mean_y)) ** 2
        v = n1 / self.semi_w ** 2 + n2 / self.semi_h ** 2
        return v < 1

    def rescale_to_nstd(self, x_vals, y_vals, weights):
        target = integrate.quad(normal_distribution, -self.n_std, self.n_std)[0]
        pct_in = self.pct_points_in(x_vals, y_vals, weights)
        required_accuracy = 1.5  # in percent
        f_scale = 1
        i = 1
        while abs(target - pct_in) * 100 > required_accuracy:
            self.rescale_axes(f_scale)
            pct_in = self.pct_points_in(x_vals, y_vals, weights)
            if i > 1900:
                print("Over 1900: ", target - pct_in)
            if pct_in < 1e-4:
                f_scale = 2
            else:
                f_scale = (target / pct_in) ** 0.5

            assert (i < 2000), "Rescaling got stuck"
            i += 1


def normal_distribution(x: float):
    return 1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x ** 2)


def weighted_median(data, weights):
    """
    Args:
      data (list or numpy.array): data
      weights (list or numpy.array): weights
    """
    data = np.array(data).squeeze()
    weights = np.array(weights).squeeze()

    s_data, s_weights = map(np.array, zip(*sorted(zip(data, weights))))
    midpoint = 0.5 * sum(s_weights)
    if any(weights > midpoint):
        w_median = (data[weights == np.max(weights)])[0]
    else:
        cs_weights = np.cumsum(s_weights)
        idx = np.where(cs_weights <= midpoint)[0][-1]
        if cs_weights[idx] == midpoint:
            w_median = np.mean(s_data[idx:idx + 2])
        else:
            w_median = s_data[idx + 1]
    return w_median


def split_ul16(_years: list):
    years = copy.deepcopy(_years)
    if "UL16" in years:
        years.remove("UL16")
        years.append("UL16preVFP")
        years.append("UL16postVFP")
    return years


def pt_of_z_edge(m_a, m_h, m_z=91.187):
    """
    Returns the value of the characteristic p_t(Z) edge in the AZH
    decay for a given combination of m_A and m_H.
    """
    x = (m_a**2 - m_h**2 - m_z**2) ** 2 - 4 * m_h**2 * m_z**2
    return 1 / (2 * m_a) * np.sqrt(x)


def collection_key(x):
    """
    Utility to make a key out of a collection plus
    position tuple provided in config.yaml
    """
    if isinstance(x, list):
        return f"{x[0]}_{x[1]}"
    else:
        return x


def is_valid_set(channel, region, svar):
    return True