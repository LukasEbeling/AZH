#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import mplhep as hep
import matplotlib.pyplot as plt
import os

from plot_utils import Theory, Limit, BestLimit, load_masses

plt.style.use(hep.style.CMS)

BRAZILIAN_GREEN = "#009C3B"
BRAZILIAN_GOLD = "#FFDF00"  

LUMI = 41.48
#LUMI = 138


def plot_limits(channel="inv", region="all", year="UL17", MH=400, tanb=1):
    masses = [(mA,mH) for mA,mH in load_masses() if mH==MH]
    masses = sorted(masses, key=lambda x: x[0])
    
    theory = Theory(tanb)
    limit = BestLimit(channel,region,year)
    
    points = [m[0] for m in masses]
    plot_range = range(min(points),max(points),10)

    tan = [theory.get_inclusive(mA,MH) for mA in plot_range]
    #tan = [theory.get_invisible(mA,400) for mA in plot_range]

    down2s = [limit.load(mA,mH,pct="2.5%") for mA,mH in masses]
    down1s = [limit.load(mA,mH,pct="16.0%") for mA,mH in masses]
    median = [limit.load(mA,mH,pct="50.0%") for mA,mH in masses]
    up1s = [limit.load(mA,mH,pct="84.0%") for mA,mH in masses]
    up2s = [limit.load(mA,mH,pct="97.5%") for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Private work',data=True, lumi=LUMI, year="2017")
            
    axes.plot(list(plot_range),tan,color="red",label=r"theory $\sigma \times$BR"+"\n"+r"(tan$\beta$="+str(tanb)+")")
    axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label="2 std. deviation")
    axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label="1 std. deviation")
    axes.plot(points,median,color="black",marker = "o",label="median expected")

    handles, labels = axes.get_legend_handles_labels()
    axes.legend(
        reversed(handles),
        reversed(labels),
        title=f"$m_H$ = {MH}GeV",
        fontsize=18,
        title_fontsize=18
    )

    plt.xlabel(f'$m_A$ [GeV]',fontsize=25)
    plt.ylabel(r'$\sigma(pp\rightarrow A) \times BR(A\rightarrow ZH\rightarrow Zt\bar t)$ [pb]',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axes.set_yscale('log')
    
    os.makedirs("limits",exist_ok=True)
    fig.savefig(f"limits/limit_{MH}.png",bbox_inches='tight')
    fig.savefig(f"limits/limit_{MH}.pdf",bbox_inches='tight')
    plt.close()  

    
def compare_regions(channel="inv", year="UL17", MH=400):
    masses = [(mA,mH) for mA,mH in load_masses() if mH==MH]
    masses = sorted(masses, key=lambda x: x[0])

    limit = BestLimit(channel,"all",year)
    down2s = [limit.load(mA,mH,pct="2.5%") for mA,mH in masses]
    down1s = [limit.load(mA,mH,pct="16.0%") for mA,mH in masses]
    median = [limit.load(mA,mH,pct="50.0%") for mA,mH in masses]
    up1s = [limit.load(mA,mH,pct="84.0%") for mA,mH in masses]
    up2s = [limit.load(mA,mH,pct="97.5%") for mA,mH in masses]

    SR1 = [BestLimit(channel,"SR_2B_6J",year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    SR2 = [BestLimit(channel,"SR_2B_5J",year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    SR3 = [BestLimit(channel,"SR_1B_6J",year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    SR4 = [BestLimit(channel,"SR_1B_5J",year).load(mA,mH,pct="50.0%") for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Private work',data=True, lumi=LUMI, year="2017")
    points = [m[0] for m in masses]

    axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label="2 std. deviation")
    axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label="1 std. deviation")
    axes.plot(points,median,color="black",marker = "o",label="median expected")
    axes.plot(points,SR1,color="red",marker = "o",label="median expected (b=2 j=6)")
    axes.plot(points,SR2,color="blue",marker = "o",label="median expected (b=2 j=5)")
    axes.plot(points,SR3,color="green",marker = "o",label="median expected (b=1 j=6)")
    axes.plot(points,SR4,color="purple",marker = "o",label="median expected (b=1 j=5)")

    handles, labels = axes.get_legend_handles_labels()
    axes.legend(
        reversed(handles),
        reversed(labels),
        title=f"$m_H$ = {MH}GeV",
        fontsize=18,
        title_fontsize=18
    )

    plt.xlabel(f'$m_A$ [GeV]',fontsize=25)
    plt.ylabel(r'$\sigma(pp\rightarrow A) \times BR(A\rightarrow ZH\rightarrow Zt\bar t)$ [pb]',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axes.set_yscale('log')
    
    os.makedirs("limits",exist_ok=True)
    fig.savefig(f"limits/regions_{MH}.png",bbox_inches='tight')
    fig.savefig(f"limits/regions_{MH}.pdf",bbox_inches='tight')
    plt.close()  

def compare_vars(channel="inv", region="all", year="UL17", MH=400):
    masses = [(mA,mH) for mA,mH in load_masses() if mH==MH]
    masses = sorted(masses, key=lambda x: x[0])

    limit = BestLimit(channel,"all",year)
    down2s = [limit.load(mA,mH,pct="2.5%") for mA,mH in masses]
    down1s = [limit.load(mA,mH,pct="16.0%") for mA,mH in masses]
    median = [limit.load(mA,mH,pct="50.0%") for mA,mH in masses]
    up1s = [limit.load(mA,mH,pct="84.0%") for mA,mH in masses]
    up2s = [limit.load(mA,mH,pct="97.5%") for mA,mH in masses]

    limit_mta = [Limit(channel,"MTA",region,year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    limit_met = [Limit(channel,"MET",region,year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    limit_mh = [Limit(channel,"MH",region,year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    limit_twd = [Limit(channel,"2DEllipses",region,year).load(mA,mH,pct="50.0%") for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Private work',data=True, lumi=LUMI, year="2017")
    points = [m[0] for m in masses]

    #axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label="2 std. deviation")
    #axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label="1 std. deviation")
    #axes.plot(points,median,color="black",marker = "o",label="median expected")
    axes.plot(points,limit_mta,color="red",marker = "o",label=r"median expected ($m_{T,A}$)")
    axes.plot(points,limit_met,color="purple",marker = "o",label=r"median expected ($p_T^{miss}$)")
    axes.plot(points,limit_mh,color="blue",marker = "o",label=r"median expected ($m_H$)")
    axes.plot(points,limit_twd,color="green",marker = "o",label="median expected (2D)")

    handles, labels = axes.get_legend_handles_labels()
    axes.legend(
        reversed(handles),
        reversed(labels),
        title=f"$m_H$ = {MH}GeV",
        fontsize=18,
        title_fontsize=18
    )

    plt.xlabel(f'$m_A$ [GeV]',fontsize=25)
    plt.ylabel(r'$\sigma(pp\rightarrow A) \times BR(A\rightarrow ZH\rightarrow Zt\bar t)$ [pb]',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axes.set_yscale('log')
    
    os.makedirs("limits",exist_ok=True)
    fig.savefig(f"limits/variables_{MH}.png",bbox_inches='tight')
    fig.savefig(f"limits/variables_{MH}.pdf",bbox_inches='tight')
    plt.close()  


if __name__ == "__main__":
    compare_vars(MH=400)
    compare_vars(MH=600)
    #compare_regions()
    plot_limits(tanb=0.5,MH=400)
    plot_limits(tanb=0.5,MH=600)
    plot_limits(tanb=0.5,MH=800)
    plot_limits(tanb=0.5,MH=1000)



