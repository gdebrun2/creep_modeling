import os
from config import MICROSTRUCTURE_DIR, RESULTS_DIR
from plotting.common import new_fig_ax, material_params_suptitle, style_axes, save_fig
import numpy as np
from data_utils import CreepData, SimResults
from plotting.mechanical import EXP_COLOR, plot_mean_with_std, SIM_COLOR
from pathlib import Path

micro_dir = "/Users/gtdebru/creep_modeling/microstructures"
run_path = "/Users/gtdebru/creep_modeling/utils/run.py"
res_path = RESULTS_DIR / "creep_calibration"
plot_dir = Path("/Users/gtdebru/creep_modeling/figures")

load = 525
igas = 0

micros = os.listdir(micro_dir)
micros = sorted(
    [file for file in micros if file.endswith(".dat") and len(file.split("_")) == 1],
    key=lambda x: int(x.strip("micro").strip(".dat")),
)


for micro in micros:
    cmd = f"python {run_path}"
    cmd += f" --sim_case creep_calibration --load {load} --igas {igas}"
    cmd += f" --micro {micro_dir + "/" + micro}"
    cmd += f" --fftw_save {micro_dir + "/" + micro.split('.')[0] + "_fftw.save"}"
    os.system(cmd)

nruns = len(micros)
results_dirs = os.listdir(res_path)
results_dirs = [dir for dir in results_dirs if not dir.startswith(".")]
results_dirs = sorted(results_dirs, key=lambda x: int(x))[-nruns:]
results = [str(res_path) + "/" + dir + "/evpfft_outputs" for dir in results_dirs]
results = [
    SimResults.load(res, micro_path=micros[i], skip_vtk=True)
    for i, res in enumerate(results)
]


def plot_strain_vs_time(
    sim: SimResults,
    ax,
):

    sim_strain = sim.epav33
    sim_time = sim.sim_time.copy() / 3600

    ax.plot(
        sim_time,
        sim_strain * 100,
        zorder=10,
        linestyle="--",
    )

    return ax


fig, ax = new_fig_ax()
ax.plot([1], [1], alpha=1, color="black", linestyle="--", label="Sim.")
title = "Strain vs Time"
material_params_suptitle(fig, start_dir=results_dirs[0], gas=igas == 2)
title += f" {load} MPa"
ylabel = "Plastic Strain (%)"
xlabel = "Time (hrs)"

exp = CreepData.load(skip_roughness=True).get_load_vals(load)
exp_strain = exp["mean_plastic_strain"]
std_strain = exp["std_plastic_strain"]
exp_time = exp["time"][1:].copy() / 3600

plot_mean_with_std(
    ax, exp_time, exp_strain * 100, std_strain * 100, label="Exp.", color=EXP_COLOR
)

for res in results:
    plot_strain_vs_time(res, ax)

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)
ax.legend()
style_axes(ax)

save_fig(fig, Path(plot_dir) / f"micro_sweep_{load}.png")
