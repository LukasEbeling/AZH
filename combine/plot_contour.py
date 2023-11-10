#!/nfs/dust/cms/user/ebelingl/anaconda3/envs/py311/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import os

from matplotlib.tri import Triangulation, LinearTriInterpolator
from plot_utils import Theory, Limit, load_masses

plt.style.use(hep.style.CMS)


def plot_plane(tanb):

    year = "UL17"
    mass_points = load_masses()

    grid = np.array(
        [
            [MA, MH]
            for MA in range(500, 2000, 2)
            for MH in range(400, int(MA*4/5)+50, 2)
        ]
    )

    fig, ax = plt.subplots(figsize=(11, 8))

    limit = Limit()
    theory = Theory(tanb)

    expected_values = np.array(
        [limit.load(MA, MH) for MA, MH in mass_points]
    )

    triObj = Triangulation(mass_points[:, 0], mass_points[:, 1])
    f_interpolate = LinearTriInterpolator(triObj, expected_values)
    expected_values_interpolated = f_interpolate(grid[:, 0], grid[:, 1])
    mask = expected_values_interpolated.mask
    expected_values_interpolated = expected_values_interpolated.data[~mask]
    grid=grid[~mask]


    theory_values = np.array(
        #[theory.get_inclusive(MA, MH) for MA, MH in grid]
        [theory.get_invisible(MA, MH) for MA, MH in grid]
    )

    theory_limit_ratio = theory_values / expected_values_interpolated
    
    contour = ax.tricontourf(
        grid[:, 0], 
        grid[:, 1], 
        theory_limit_ratio,
        levels=10,
            cmap="viridis"
    )

    ax.tricontour(
        grid[:, 0],
        grid[:, 1],
        theory_limit_ratio,
        levels=[0, 1],
        colors="red",
    )

    ax.scatter(mass_points[:, 0], mass_points[:, 1], color="orange")

    fig.colorbar(contour, ax=ax, label=r"$\sigma_{Theory}/\sigma_{Expected}$")
    ax.set_xlabel(r"m$_{A}$ [GeV]")
    ax.set_ylabel(r"m$_{H}$ [GeV]")
    ax.set_xlim(500, 2000)
    ax.set_ylim(400, 1600)

    os.makedirs("limits",exist_ok=True)
    plt.savefig(f"limits/mass_plane_{tanb}.png")
    plt.savefig(f"limits/mass_plane_{tanb}.pdf")
    plt.close()


if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("--year", type=str, help="UL16/UL17/UL18/combined")
    #args = parser.parse_args()

    for tanb in [0.5,1,3]:
        plot_plane(tanb)