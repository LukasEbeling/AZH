#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import mplhep as hep
import matplotlib.pyplot as plt
import os

from plot_utils import Theory, Limit, BestLimit, load_masses

plt.style.use(hep.style.CMS)

BRAZILIAN_GREEN = "#009C3B"
BRAZILIAN_GOLD = "#FFDF00"  


def plot_combination(year="UL17"):
    masses = [(600,400),(750,400),(800,400),(1000,400)]
    limit = Limit("combined","MET","all",year) 
    down2s = [limit.load(mA,mH,pct="2.5%") for mA,mH in masses]
    down1s = [limit.load(mA,mH,pct="16.0%") for mA,mH in masses]
    median = [limit.load(mA,mH,pct="50.0%") for mA,mH in masses]
    up1s = [limit.load(mA,mH,pct="84.0%") for mA,mH in masses]
    up2s = [limit.load(mA,mH,pct="97.5%") for mA,mH in masses]

    inv = [Limit("inv","MET","all",year).load(mA,mH,pct="50.0%") for mA,mH in masses]
    lep = [Limit("lep","2DEllipses","all",year).load(mA,mH,pct="50.0%") for mA,mH in masses]

    fig, axes = plt.subplots(figsize=(12, 8))
    hep.cms.label(ax=axes,llabel='Private work',data=True, lumi=41.48, year="2017")
    points = [m[0] for m in masses]

    axes.fill_between(points, down2s, up2s, color=BRAZILIAN_GOLD, label="2 std. deviation")
    axes.fill_between(points, down1s, up1s, color=BRAZILIAN_GREEN, label="1 std. deviation")
    axes.plot(points,median,color="black",marker = "o",label="median expected")

    axes.plot(points,lep,color="blue",marker = "o",label=r"median expected ($Z\rightarrow ll$)")
    axes.plot(points,inv,color="red",marker = "o",label=r"median expected ($Z\rightarrow\nu\nu$)")

    handles, labels = axes.get_legend_handles_labels()
    axes.legend(
        reversed(handles),
        reversed(labels),
        title=f"$m_H$ = 400 GeV",
        fontsize=18,
        title_fontsize=18
    )


    plt.xlabel(f'$m_A$ [GeV]',fontsize=25)
    plt.ylabel(r'$\sigma(pp\rightarrow A) \times BR(A\rightarrow ZH\rightarrow Zt\bar t)$ [pb]',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axes.set_yscale('log')
    
    os.makedirs("limits",exist_ok=True)
    fig.savefig(f"limits/combination.png")
    fig.savefig(f"limits/combination.pdf")
    plt.close() 



if __name__ == "__main__":
    plot_combination()




