import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator
from SNSPlottingWidget import SNSPlottingWidget
pw = SNSPlottingWidget()
S_I_Glua2 = []
I_S_Glua2 = []
glua2_maps = {}
num_cells = {
    "MAX1": 36,
    "MAX2": 23,
    "MAX3": 35
}
# scale = 
op_folder = "./Figures/HaloGluA2_analysis/"
def measure_cell(file_path,signal_name):
    df = pd.read_csv(file_path)

    df.columns = df.columns.str.strip()
    # Extract ROI from Label
    df["roi"] = df["Label"].astype(str).str.split(":").str[-1].str.strip()
    # Find bg mean
    bg_rows = df[df["roi"].str.lower() == "bg"]
    
    # breakpoint()
    if bg_rows.empty:
        print(f"No bg row found in {cell_name} {signal_name}")
        bg_mean = 0
    else:
        bg_mean = bg_rows["Mean"].mean()
    
    # Keep only soma and numeric rows
    is_soma = df["roi"].str.lower() == "soma"
    is_numeric = df["roi"].str.match(r"^\d+$", na=False)

    df_keep = df.loc[is_soma | is_numeric].copy()

    df_keep.loc[df_keep["roi"].str.lower() == "soma", "roi"] = 0
    df_keep["roi"] = pd.to_numeric(df_keep["roi"])

    df_keep["compartment"] = "Dendrite"
    df_keep.loc[df_keep["roi"] == 0, "compartment"] = "Soma"

    # Background subtraction
    df_keep["Mean_bg_subtracted"] = df_keep["Mean"] - bg_mean
    df_keep["IntDen_bg_subtracted"] = df_keep["IntDen"] - df_keep["Area"]*bg_mean

    # Metadata
    df_keep["cell_id"] = "{}_{}".format(rep,j)
    df_keep["cell_name"] = cell_name
    df_keep["signal"] = signal_name
    # df_keep["source_folder"] = folder_name

    # Keep only useful columns first, then everything else if wanted
    cols_to_front = ["cell_id", "cell_name", "signal", "roi", "compartment", "Mean", "IntDen","Mean_bg_subtracted","IntDen_bg_subtracted"]
    remaining_cols = [c for c in df_keep.columns if c not in cols_to_front]
    df_keep = df_keep[cols_to_front + remaining_cols]
    return df_keep
def bin_roi_data(data,bin_width):
    x = data[:, 0]  # column 1
    y = data[:, 1]  # column 2
    bins = np.arange(x.min(), x.max() + bin_width, bin_width)
    bin_indices = np.digitize(x, bins)
    bin_means = []
    bin_centers = []

    for i in range(1, len(bins)):
        mask = bin_indices == i
        if np.any(mask):
            bin_means.append(y[mask].mean())
            bin_centers.append((bins[i-1] + bins[i]) / 2)
    return np.array([np.array(bin_centers),np.array(bin_means)])

def read_roi_csv_file(file_path1, bg_value1,file_path2, bg_value2,bin_width=10):
    data1 = np.loadtxt(file_path1, delimiter=',', skiprows=1, usecols=(1, 2))
    data1[:, 1] = data1[:, 1] - bg_value1
    binned_data1 = bin_roi_data(data1, bin_width=bin_width)
    data2 = np.loadtxt(file_path2, delimiter=',', skiprows=1, usecols=(1, 2))
    data2[:, 1] = data2[:, 1] - bg_value2
    binned_data2 = bin_roi_data(data2, bin_width=bin_width)
    # binned_data1[1] = binned_data1[1] 
    # binned_data2[1] = binned_data2[1] - bg_value2
    return binned_data1,binned_data2    

for rep in ["MAX1", "MAX2", "MAX3"]:
    base_path = f"./Protein/ROIS/{rep}/"
    glua2_maps[rep] = {}
    for j in range(1, num_cells[rep]+1):
        cell_name = f"cell_{j}"
        cell_path = os.path.join(base_path, cell_name)

        # Reference file used only to determine the GluA2 mapping
        ref_file = os.path.join(cell_path, "SUM", "PSD95", "measures", "all_data.csv")
        if not os.path.exists(ref_file):
            print(f"Missing reference file for {cell_name}: {ref_file}")
            continue

        ref_df = pd.read_csv(ref_file)
        ref_df.columns = ref_df.columns.str.strip()

        has_JF635i = ref_df["Label"].astype(str).str.contains("JF635i", case=False, na=False).any()

        if has_JF635i:
            glua2_map = {
                "surf_GluA2": "GluA2_1",
                "int_GluA2": "GluA2_2",
                "MAP2": "MAP2"
            }
        else:
            glua2_map = {
                "surf_GluA2": "GluA2_2",
                "int_GluA2": "GluA2_1",
                "MAP2": "MAP2"
            }

        glua2_maps[rep][cell_name] = glua2_map

        for signal_name, folder_name in glua2_map.items():
            file_path = os.path.join(cell_path, "SUM", folder_name, "measures", "all_data.csv")
            # MAP2_path = os.path.join(cell_path, "MAX", "MAP2", "measures", "all_data.csv")
            if not os.path.exists(file_path):
                print(f"Missing file for {cell_name}, {signal_name}: {file_path}")
                continue

            df_keep = measure_cell(file_path,signal_name)
            # breakpoint()
            if has_JF635i:
                S_I_Glua2.append(df_keep)
            else:
                I_S_Glua2.append(df_keep)


def geometric_mean_bootstrap(x, y, ci=95, n_boot=1000,seed=2026):
    x = np.asarray(x)
    y = np.asarray(y)
    
    G_samples = []
    np.random.seed(seed)
    for _ in range(n_boot):
        xb = np.random.choice(x, size=len(x), replace=True)
        yb = np.random.choice(y, size=len(y), replace=True)
        
        G_samples.append(np.sqrt(np.mean(xb) * np.mean(yb)))
    
    G_samples = np.array(G_samples)
    
    G = np.mean(G_samples)
    sd = np.std(G_samples, ddof=1)
    alpha = (100 - ci) / 2
    ci_low = np.percentile(G_samples, alpha, axis=0)
    ci_high = np.percentile(G_samples, 100 - alpha, axis=0)
    return G_samples, G, sd, ci_low, ci_high

def geometric_mean_bootstrap2D(x, y, n_boot=10000, ci=95, seed=2026):
    """
    Bootstrap estimate of sqrt(mean(x) * mean(y)) when x and y are matrices.

    Parameters
    ----------
    x : array-like, shape (n_samples_x, n_features)
        Each row is a sample.
    y : array-like, shape (n_samples_y, n_features)
        Each row is a sample.
    n_boot : int
        Number of bootstrap samples.
    ci : float
        Confidence interval percentage, e.g. 95 for a 95% CI.
    seed : int
        Random seed.

    Returns
    -------
    G : ndarray, shape (n_features,)
        Bootstrap estimate of sqrt(mean(x) * mean(y)).
    ci_low : ndarray, shape (n_features,)
        Lower bound of the confidence interval.
    ci_high : ndarray, shape (n_features,)
        Upper bound of the confidence interval.
    G_samples : ndarray, shape (n_boot, n_features)
        Bootstrap distribution.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if x.ndim != 2 or y.ndim != 2:
        raise ValueError("x and y must both be 2D arrays with rows as samples.")
    if x.shape[1] != y.shape[1]:
        raise ValueError("x and y must have the same number of columns/features.")

    rng = np.random.default_rng(seed)

    n_x, n_feat = x.shape
    n_y, _ = y.shape

    G_samples = np.empty((n_boot, n_feat), dtype=float)

    for i in range(n_boot):
        xb = x[rng.choice(n_x, size=n_x, replace=True), :]
        yb = y[rng.choice(n_y, size=n_y, replace=True), :]

        mx = np.mean(xb, axis=0)
        my = np.mean(yb, axis=0)

        G_samples[i] = np.sqrt(mx * my)

    G = np.mean(G_samples, axis=0)
    G_std = np.std(G_samples, axis=0, ddof=1)
    G_sem = G_std / np.sqrt(n_boot)
    alpha = (100 - ci) / 2
    ci_low = np.percentile(G_samples, alpha, axis=0)
    ci_high = np.percentile(G_samples, 100 - alpha, axis=0)

    return G, G_std,G_sem, ci_low, ci_high, G_samples
# Combine all cells/signals
S_I_Glua2 = pd.concat(S_I_Glua2, ignore_index=True)
I_S_Glua2 = pd.concat(I_S_Glua2, ignore_index=True)

print(S_I_Glua2.head())
measure = "IntDen"
# breakpoint()
# Save combined result
# S_I_Glua2.to_csv("all_cells_S_I_glua2_processed.csv", index=False)
# I_S_Glua2.to_csv("all_cells_I_S_glua2_processed.csv", index=False)


# Pivot so surf/int become columns
ratio1_df = S_I_Glua2.pivot_table(
    index=["cell_id", "roi", "compartment"],
    columns="signal",
    values="{}_bg_subtracted".format(measure)
).reset_index()

# Remove negatives
ratio1_df = ratio1_df[
    (ratio1_df["surf_GluA2"] >= 0) &
    (ratio1_df["int_GluA2"] > 0)   # avoids divide-by-zero too
]
# Compute ratio
ratio1_df = (
    ratio1_df.groupby(['cell_id', 'compartment'], as_index=False)
      .agg({
          'MAP2': 'sum',
          'int_GluA2': 'sum',
          'surf_GluA2': 'sum'
      })
)
ratio1_df["surf_to_int_ratio"] = (
    (ratio1_df["surf_GluA2"] / ratio1_df["int_GluA2"])
)
# breakpoint()
plt.hist(ratio1_df["surf_to_int_ratio"])
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.savefig(os.path.join(op_folder,"Glua2_surface_to_internal_histogram_S_I_{}".format(measure)))
plt.show()
cell_avg1 = (
    ratio1_df
    .groupby(["cell_id", "compartment"])["surf_to_int_ratio"]
    .mean()
    .reset_index()
)

paired = (
    cell_avg1
    .pivot(index="cell_id", columns="compartment", values="surf_to_int_ratio")
    .dropna(subset=["Soma", "Dendrite"])
    .reset_index()
)

# Convert back to long format for seaborn/statannotations
cell_avg_paired1 = paired.melt(
    id_vars="cell_id",
    value_vars=["Soma", "Dendrite"],
    var_name="compartment",
    value_name="surf_to_int_ratio"
)

print(cell_avg_paired1.groupby("compartment").agg(['min', 'max','mean', 'median','std']))
print(f"Number of paired cells: {paired.shape[0]}")


ratio2_df = I_S_Glua2.pivot_table(
    index=["cell_id", "roi", "compartment"],
    columns="signal",
    values="{}_bg_subtracted".format(measure)
).reset_index()

ratio2_df = ratio2_df[
    (ratio2_df["surf_GluA2"] >= 0) &
    (ratio2_df["int_GluA2"] > 0)   # avoids divide-by-zero too
]
ratio2_df = (
    ratio2_df.groupby(['cell_id', 'compartment'], as_index=False)
      .agg({
          'MAP2': 'sum',
          'int_GluA2': 'sum',
          'surf_GluA2': 'sum'
      })
)

# Compute ratio
ratio2_df["surf_to_int_ratio"] = (
    (ratio2_df["surf_GluA2"]/ ratio2_df["int_GluA2"])
)
# import matplotlib.pyplot as plt

# values = [1, 2, 2, 3, 3, 3, 4, 5]

plt.hist(ratio2_df["surf_to_int_ratio"])
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.savefig(os.path.join(op_folder,"Glua2_surface_to_internal_histogram_I_S_{}".format(measure)))

plt.show()
cell_avg2 = (
    ratio2_df
    .groupby(["cell_id", "compartment"])["surf_to_int_ratio"]
    .mean()
    .reset_index()
)
paired = (
    cell_avg2
    .pivot(index="cell_id", columns="compartment", values="surf_to_int_ratio")
    .dropna(subset=["Soma", "Dendrite"])
    .reset_index()
)

# Convert back to long format for seaborn/statannotations
cell_avg_paired2 = paired.melt(
    id_vars="cell_id",
    value_vars=["Soma", "Dendrite"],
    var_name="compartment",
    value_name="surf_to_int_ratio"
)
print(cell_avg_paired2.groupby("compartment").agg(['min', 'max','mean', 'median','std']))
print(f"Number of paired cells: {paired.shape[0]}")

x = "compartment"
y = "surf_to_int_ratio"
hue = None#["Soma","Dendrite"]
xlab = "Compartment"
ylab1 = r"$\frac{sGluA2(JF635i)}{iGluA2(JF549)}$"
title = ""
stat_test = "Wilcoxon"
fsize = 22
lsize= 22

fig_file = os.path.join(op_folder,"Glua2_surface_to_internal_soma_vs_dendrite_S_I_{}".format(measure))
pairs = [("Soma","Dendrite")]
my_pal1 = {"Soma":"#277DA1","Dendrite": '#43AA8B'}
ax_label=1
x_pos = [0,1]
y_pos = np.stack((cell_avg_paired1.groupby("compartment")["surf_to_int_ratio"].apply(list)[::-1])).T
pw.Swarmboxplotcombo(data=cell_avg_paired1,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab1,
                             pairs=pairs,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 1
                             )

ylab2 = r"$\frac{sGluA2(JF549i)}{iGluA2(JF635)}$"
fig_file = os.path.join(op_folder,"Glua2_surface_to_internal_soma_vs_dendrite_I_S_{}".format(measure))
# breakpoint()
y_pos = np.stack((cell_avg_paired2.groupby("compartment")["surf_to_int_ratio"].apply(list)[::-1])).T
pw.Swarmboxplotcombo(data=cell_avg_paired2,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab2,
                             pairs=pairs,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 1
                             )
# breakpoint()
G_soma, G_soma_mean, G_soma_std, G_soma_clow,G_soma_chigh = geometric_mean_bootstrap(cell_avg_paired1[cell_avg_paired1["compartment"] == "Soma"]["surf_to_int_ratio"],
                           cell_avg_paired2[cell_avg_paired2["compartment"] == "Soma"]["surf_to_int_ratio"])
G_dend, G_dend_mean, G_dend_std, G_dend_clow,G_dend_chigh = geometric_mean_bootstrap(cell_avg_paired1[cell_avg_paired1["compartment"] == "Dendrite"]["surf_to_int_ratio"],
                           cell_avg_paired2[cell_avg_paired2["compartment"] == "Dendrite"]["surf_to_int_ratio"])

print("Soma geometric mean ratio: {:.3f} ± {:.3f} (SEM), ci = {:.3f}, {:.3f}".format(G_soma_mean, G_soma_std/np.sqrt(G_soma.shape[0]),G_soma_clow, G_soma_chigh))
print("Dendrite geometric mean ratio: {:.3f} ± {:.3f} (SEM), ci = {:.3f},{:.3f}".format(G_dend_mean, G_dend_std/np.sqrt(G_dend.shape[0]), G_dend_clow, G_dend_chigh))
boots_trapped_data = pd.concat([
    pd.DataFrame({'surf_to_int_ratio': G_soma, 'compartment': 'Soma'}),
    pd.DataFrame({'surf_to_int_ratio': G_dend, 'compartment': 'Dendrite'})
], ignore_index=True)
fig_file = os.path.join(op_folder,"Glua2_surface_to_internal_soma_vs_dendrite_geometric_{}".format(measure))
ylab3 = r"Estimated $\frac{sGluA2}{iGluA2}$"
pw.Swarmboxplotcombo(data=boots_trapped_data,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab3,
                             pairs=pairs,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=None,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 0,
                             circle_size= 6
                             )

# breakpoint()
LOIS = [75, 100, 125, 150]
final_ratio = {}
d1_I_S = {}
d2_I_S = {}
d1_S_I = {}
d2_S_I = {}
op_data_folder = "./Scripts/Protein_data/HaloGluA2_analysis/"
os.makedirs(op_data_folder,exist_ok=True)
w,h= 8,4
for LOI in LOIS:
    filtered_df1 = S_I_Glua2[
        (S_I_Glua2["compartment"] == "Dendrite") & (S_I_Glua2["Length"] >= LOI)
    ]
    xlab1 = r"Dendritic distance ($\mu m$)"
    bin_size = 5 #in micron
    scale = 17.2308 #pixels per micron
    d1_S_I[LOI] = []
    d2_S_I[LOI] = []
    ratios_100_S_I = []
    for idx, row in filtered_df1.iterrows():
        rep_folder = row["cell_id"].split("_")[0]
        cell_f  = row["cell_name"]
        s_folder = glua2_maps[rep_folder][cell_f]["surf_GluA2"]
        i_folder = glua2_maps[rep_folder][cell_f]["int_GluA2"]
        f1_path = os.path.join("./Protein/ROIS",rep_folder ,cell_f, "SUM", s_folder, "data", "{}.csv".format(row["roi"]))
        f2_path = os.path.join("./Protein/ROIS",rep_folder ,cell_f, "SUM", i_folder, "data", "{}.csv".format(row["roi"]))
        bg1 = row["Mean"] - row["Mean_bg_subtracted"]
        bg2 = row["Mean"] - row["Mean_bg_subtracted"]
        d1,d2 = read_roi_csv_file(f1_path,bg1,f2_path,bg2,bin_width=bin_size)
        if (np.all(d1[1] >= 0)) and np.all(d2[1] > 0):
            # breakpoint()
            last_i = np.where(d1[0] <= LOI+1)[0][-1]+1
            d1_S_I[LOI].append(d1[1, :last_i])
            d2_S_I[LOI].append(d2[1,:last_i])
            ratios_100_S_I.append((d1[1]/d2[1])[d1[0] <= LOI+1])
        # print(idx, row["cell_id"], row["Length"])
    # for idx, row in filtered_df1.iterrows():
    ratios_100_S_I = np.array(ratios_100_S_I)
    mean_ratios_100_S_I = ratios_100_S_I.mean(axis=0)
    std_ratios_100_S_I = ratios_100_S_I.std(axis=0)
    sem_ratios_100_S_I = std_ratios_100_S_I / np.sqrt(ratios_100_S_I.shape[0])
    # std = data.std(axis=0)

    plt.figure(figsize=(w,h))

    # plot each row with low alpha
    x = np.arange(0,LOI,bin_size)  # assuming all have the same x values
    # for row in ratios_100_S_I:
    #     plt.plot(x, row, alpha=0.1)

    # plot mean with std as error bars
    plt.errorbar(x, mean_ratios_100_S_I,yerr=sem_ratios_100_S_I, 
                        color="#277DA1",label='data',marker='o',
                        markersize=10,linestyle='None' ,zorder=6)

    plt.xlabel(xlab1)
    plt.ylabel(ylab1)
    plt.xlim(-2, LOI+1)
    plt.xticks(np.arange(0, LOI+1, 25))
    plt.ylim(0., 1.2)
    plt.tight_layout()
    ax = plt.gca()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.title("Rows with mean ± std")
    pw.SaveFigures(os.path.join(op_folder,"Glua2_surface_to_internal_dendrite_S_I_profile_{}".format(LOI)))

    # plt.show()

    filtered_df2 = I_S_Glua2[
        (I_S_Glua2["compartment"] == "Dendrite") & (I_S_Glua2["Length"] >= LOI)
    ]
    d1_I_S[LOI] = []
    d2_I_S[LOI] = []
    ratios_100_I_S = []
    for idx, row in filtered_df2.iterrows():
        rep_folder = row["cell_id"].split("_")[0]
        cell_f  = row["cell_name"]
        s_folder = glua2_maps[rep_folder][cell_f]["surf_GluA2"]
        i_folder = glua2_maps[rep_folder][cell_f]["int_GluA2"]
        f1_path = os.path.join("./Protein/ROIS",rep_folder ,cell_f, "SUM", s_folder, "data", "{}.csv".format(row["roi"]))
        f2_path = os.path.join("./Protein/ROIS",rep_folder ,cell_f, "SUM", i_folder, "data", "{}.csv".format(row["roi"]))
        bg1 = row["Mean"] - row["Mean_bg_subtracted"]
        bg2 = row["Mean"] - row["Mean_bg_subtracted"]
        d1,d2 = read_roi_csv_file(f1_path,bg1,f2_path,bg2,bin_width=bin_size)
        if (np.all(d1[1] >= 0)) and np.all(d2[1] > 0):
            # breakpoint()
            last_i = np.where(d1[0] <= LOI+1)[0][-1]+1
            d1_I_S[LOI].append(d1[1,:last_i])
            d2_I_S[LOI].append(d2[1,:last_i])
            ratios_100_I_S.append((d1[1]/d2[1])[d1[0] <= LOI+1])
    ratios_100_I_S = np.array(ratios_100_I_S)
    mean_ratios_100_I_S = ratios_100_I_S.mean(axis=0)
    std_ratios_100_I_S = ratios_100_I_S.std(axis=0)
    sem_ratios_100_S_I = std_ratios_100_I_S / np.sqrt(ratios_100_I_S.shape[0])
    d1_I_S[LOI] = np.array(d1_I_S[LOI])
    d2_I_S[LOI] = np.array(d2_I_S[LOI])
    d1_S_I[LOI] = np.array(d1_S_I[LOI])
    d2_S_I[LOI] = np.array(d2_S_I[LOI])
    plt.figure(figsize=(w,h))
    plt.errorbar(x, mean_ratios_100_I_S,yerr=sem_ratios_100_S_I, 
                color="#43AA8B",label='data',marker='o',
                        markersize=10,linestyle='None' ,zorder=6)

    plt.xlabel(xlab1)
    plt.ylabel(ylab2)
    plt.xlim(-2, LOI+1)
    plt.xticks(np.arange(0, LOI+1, 25))
    # plt.ylim(1, 9)
    # plt.title("Rows with mean ± std")
    plt.tight_layout()
    ax = plt.gca()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pw.SaveFigures(os.path.join(op_folder,"Glua2_surface_to_internal_dendrite_I_S_profile_{}".format(LOI)))
    plt.show()
    # breakpoint()
    G, G_sem,G_std, ci_low,ci_high, G_samples = geometric_mean_bootstrap2D(ratios_100_S_I, ratios_100_I_S, n_boot=10000, ci=95, seed=2026)
    final_ratio[LOI] = [LOI,x,G, ci_low,ci_high]
    # print(np.sqrt(mean_ratios_100_I_S*mean_ratios_100_S_I))
    
    concat_data = np.vstack((d1_S_I[LOI],d1_I_S[LOI]))
    concat_data_norm = concat_data / concat_data[:,[0]]
    mean2 = np.mean(concat_data_norm,axis=0)
    std2 = np.std(concat_data_norm, axis=0, ddof=1)
    sem2 = std2 / np.sqrt(concat_data_norm.shape[0])
    x = np.arange(0, LOI, bin_size)
    # breakpoint()
    plt.figure(figsize=(w,h))
    plt.errorbar(x, mean2, yerr=sem2, color="#277DA1", label='data', marker='o', markersize=10, linestyle='None', zorder=6)
    plt.xlabel(xlab1)
    plt.ylabel("Normalized IntDen")
    plt.xlim(-2, LOI+1)
    plt.xticks(np.arange(0, LOI, 25))
    # plt.ylim(0., 1.2)
    plt.tight_layout()
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pw.SaveFigures(os.path.join(op_folder,"Glua2_intden_surf_dendrite_profile_{}".format(LOI)))
    plt.show()
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_surf_mean_profile_{}".format(LOI)),mean2)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_surf_std_profile_{}".format(LOI)),std2)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_surf_sem_profile_{}".format(LOI)),sem2)
    
    concat_data = np.vstack((d2_S_I[LOI],d2_I_S[LOI]))
    concat_data_norm = concat_data / concat_data[:,[0]]
    print("Num samples for ratio at LOI {}: {}".format(LOI,concat_data_norm.shape[0]))
    print("Num samples for ratio at LOI {}: {}".format(LOI,d1_I_S[LOI].shape[0]))
    print("Num samples for ratio at LOI {}: {}".format(LOI,d1_S_I[LOI].shape[0]))
    mean1 = np.mean(concat_data_norm,axis=0)
    std1 = np.std(concat_data_norm, axis=0, ddof=1) 
    sem1 = std1 / np.sqrt(concat_data_norm.shape[0])
    # breakpoint()
    x = np.arange(0, LOI, bin_size)
    plt.figure(figsize=(w,h))
    plt.errorbar(x, mean1, yerr=sem1, color="#43AA8B", label='data', marker='o', markersize=10, linestyle='None', zorder=6)
    plt.xlabel(xlab1)
    plt.ylabel("Normalized IntDen")
    plt.xlim(-2, LOI+1)
    plt.xticks(np.arange(0, LOI+1, 25))
    # plt.ylim(0., )
    plt.tight_layout()
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pw.SaveFigures(os.path.join(op_folder,"Glua2_intden_int_dendrite_profile_{}".format(LOI)))
    plt.show()
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_int_mean_profile_{}".format(LOI)),mean1)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_int_std_profile_{}".format(LOI)),std1)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_int_sem_profile_{}".format(LOI)),sem1)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_ratio_mean_profile_{}".format(LOI)),G)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_ratio_std_profile_{}".format(LOI)),G_std)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_ratio_sem_profile_{}".format(LOI)),G_sem)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_ratio_cilow_profile_{}".format(LOI)),ci_low)
    np.save(os.path.join(op_data_folder,"Glua2_dendrite_ratio_cihigh_profile_{}".format(LOI)),ci_high)
plt.figure(figsize=(8,6))
colors = ["#43AA8B", "#277DA1", "#F94144", "#F3722C"]
c1 = 0
for loi,x1, G, ci_low,ci_high, in final_ratio.values():
    plt.plot(x1, G, color=colors[c1],label=loi)
    plt.fill_between(x1, ci_low, ci_high, alpha=0.3,color=colors[c1])
    c1 += 1
plt.legend(title="Dendritic Length (µm)")
plt.xlabel(xlab1)
plt.ylabel(r"$\frac{sGluA2}{iGluA2}$")
plt.ylim(0, 4)
# plt.title("Rows with mean ± std")
plt.tight_layout()
ax = plt.gca()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
pw.SaveFigures(os.path.join(op_folder,"Glua2_surface_to_internal_dendrite_inteden_total_profile"))
plt.show()

# breakpoint()
# for LOI in LOIS:
    

#     breakpoint()
