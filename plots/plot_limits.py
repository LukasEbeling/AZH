#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import mplhep as hep
import matplotlib.pyplot as plt

from utilities import Theory, Limit, load_masses

plt.style.use(hep.style.CMS)

BRAZILIAN_GREEN = "#009C3B"
BRAZILIAN_GOLD = "#FFDF00"  


def plot_limits(channel="inv", var="MET", region="all", year="UL17"):
    masses = [(mA,mH) for mA,mH in load_masses() if mH==400]
    masses = sorted(masses, key=lambda x: x[0])
    
    theory = Theory(tanb=1)
    limit = Limit(channel,var,region,year)
    
    plot_range = range(500,2000,10)
    tan = [theory.get_inclusive(mA,400) for mA in plot_range]

    down2s = [limit.load(mA,mH,pct="2.5%") for mA,mH in masses]
    down1s = [limit.load(mA,mH,pct="16.0%") for mA,mH in masses]
    median = [limit.load(mA,mH,pct="50.0%") for mA,mH in masses]
    up1s = [limit.load(mA,mH,pct="84.0%") for mA,mH in masses]
    up2s = [limit.load(mA,mH,pct="97.5%") for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Work in progress',data=True, lumi=41.48, year="2017")
            
    points = [m[0] for m in masses]
    axes.plot(list(plot_range),tan,color="red",label=r"theory tan$\beta$=1")
    axes.plot(points,median,color="black",marker = "o",label="expected 95% CL")
    axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label=r"$2\sigma$")
    axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label=r"$1\sigma$")

    axes.legend()
    plt.xlabel('$m_A$ [GeV] for H@400GeV')
    plt.ylabel('$\sigma \cdot BR$ [pb]')
    axes.set_yscale('log')
    fig.savefig(f"limit_{var}.png")
    fig.savefig(f"limit_{var}.pdf")
    plt.close()  

    
if __name__ == "__main__":
    plot_limits()

