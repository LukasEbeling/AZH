#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import mplhep as hep
import matplotlib.pyplot as plt
import re
import os


plt.style.use(hep.style.CMS)

#Config
#SUSHI = "/nfs/dust/cms/user/ebelingl/SusHi-1.7.0/sushi_outputs/"
SUSHI = "/nfs/dust/cms/user/ebelingl/sushi_outputs/"
THDMC = "/nfs/dust/cms/user/ebelingl/2HDMC-1.7.0/2HDMC_out.txt"
ANALYSIS = "/nfs/dust/cms/user/ebelingl/uhh2_106X_v2/CMSSW_10_6_28/src/UHH2/AZH/"

BR_T = 0.665
BR_Z = 0.2

BRAZILIAN_GREEN = "#009C3B"
BRAZILIAN_GOLD = "#FFDF00"

class xsecLoader():

    def __init__(self):
        self.load_2HDMC()    
    
    def load_2HDMC(self):
        with open(THDMC,"r") as f: file = list(f)
        pattern = r"MASSPOINT|A\s+->\s+Z\s+H|H\s+->\s+t\s+t"
        lines = [l for l in file if re.search(pattern,l)]
        if os.path.exists('2HDMC.txt'): 
            print("existing 2HDMC.txt used")
            return False
        with open('2HDMC.txt', 'w') as f:
            for i in range(len(lines)-3):
                if re.search(r"MASSPOINT",lines[i]) and re.search(r"H\s+->\s+t\s+t",lines[i+1]) and re.search(r"A\s+->\s+Z\s+H",lines[i+2]):
                    f.write(lines[i])
                    f.write(lines[i+1])
                    f.write(lines[i+2])

    def get_xsec(self,mA,mH,tan=1):
        path = SUSHI+ f"mA{mA}_mH{mH}_tanb{tan}.out"
        if not os.path.exists(path): return 0.0
        with open(path, "r") as f: lines = list(f)
        match = re.findall(r"\d\.\d+E.\d+", lines[18])[0]
        print(mA,mH)
        return float(match)
    
    def get_br_H(self,mA,mH,tan=1):
        with open("2HDMC.txt", "r") as f: lines = list(f)
        for i in range(len(lines)):
            if re.search(f"MASSPOINT,{tan},{mA},{mH}",lines[i]):
                match = re.findall(r"\s+\d+.\d+e.\d+", lines[i+1])[1]
                return float(match)
        return 0.0
            
    def get_br_A(self,mA,mH,tan=1):
        with open("2HDMC.txt", "r") as f: lines = list(f)
        for i in range(len(lines)):
            if re.search(f"MASSPOINT,{tan},{mA},{mH}",lines[i]):
                match = re.findall(r"\s+\d+.\d+e.\d+", lines[i+2])[1]
                return float(match)
        return 0.0
    
    def get_total(self,mA,mH,tan=1):
        return self.get_xsec(mA,mH,tan=1)*self.get_br_A(mA,mH,tan=1)*self.get_br_H(mA,mH,tan=1)*BR_T**2*BR_Z
    

class limitLoader():

    def get_observed(self,var,mA,mH): pass
    
    def get_expected(self,var,mA,mH):
        path = ANALYSIS+f"combine/UL17/expected_{mA}_{mH}_{var}.log"
        with open(path,"r") as f: lines = list(f)
        return [float(re.findall(r"\d+.\d+",line)[1]) for line in lines]
    
    def get_2s_down(self,var,mA,mH):
        return self.get_expected(var,mA,mH)[0]
    
    def get_1s_down(self,var,mA,mH):
        return self.get_expected(var,mA,mH)[1]
    
    def get_median(self,var,mA,mH):
        return self.get_expected(var,mA,mH)[2]
    
    def get_1s_up(self,var,mA,mH):
        return self.get_expected(var,mA,mH)[3]
    
    def get_2s_up(self,var,mA,mH):
        return self.get_expected(var,mA,mH)[4]

    


def plot_limits():
    masses = [(600,400),(750,400),(800,400),(1000,400)]
    
    theory = xsecLoader()
    limits = limitLoader()
    
    tan1 = [theory.get_total(mA,mH,1) for mA,mH in masses]
    down2s = [limits.get_2s_down("met",mA,mH) for mA,mH in masses]
    down1s = [limits.get_1s_down("met",mA,mH) for mA,mH in masses]
    median = [limits.get_median("met",mA,mH) for mA,mH in masses]
    up1s = [limits.get_1s_up("met",mA,mH) for mA,mH in masses]
    up2s = [limits.get_2s_up("met",mA,mH) for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Work in progress',data=True, lumi=41.48, year="2017")
            
    points = [m[0] for m in masses]
    axes.plot(points,median,color="black",marker = "o",label="expected")
    axes.plot(points,tan1,color="red",label=r"tan$\beta$=1")
    axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label=r"$2\sigma$")
    axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label=r"$1\sigma$")

    axes.legend()
    plt.xlabel('$m_T$ of A [GeV] for H@400GeV')
    plt.ylabel('$\sigma \cdot BR$ [pb]')
    axes.set_yscale('log')
    fig.savefig(ANALYSIS+"plots/limits/limit.png")
    plt.close()

def plot_plane():
    samples = [(600,400),(750,400),(800,400),(1000,400),(1000,850),(700,450),(750,650)]
    masses = [(mA,mH) for mA in range(130,1200,10) for mH in range(130,mA,10)]

    theory = xsecLoader()
    tan1 = [theory.get_total(mA,mH,1) for mA,mH in masses]
    A = [m[0] for m in masses]
    H = [m[1] for m in masses]

    fig, axes = plt.subplots(figsize=(14, 8))
    hep.cms.label(ax=axes,llabel='Work in progress',data=True, lumi=41.48, year="2017")
    plt.scatter(A, H, c=tan1, cmap='viridis')
    plt.scatter([s[0] for s in samples],[s[1] for s in samples], marker = "x", c="black")
    colorbar = plt.colorbar()
    colorbar.set_label('theory $\sigma \cdot BR$ [pb]')
    plt.xlabel('$m_T$ of A [GeV]')
    plt.ylabel('$m_T$ of H [GeV]')
    fig.savefig(ANALYSIS+"plots/limits/plane.png")
    plt.close()


    
if __name__ == "__main__":
    plot_limits()
    plot_plane()
