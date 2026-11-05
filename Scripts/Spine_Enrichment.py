import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import re
import seaborn as sns
import statannotations
from statannotations.Annotator import Annotator
sns.set(style="whitegrid")



# def load_spine_dendrite_dataframe(excel_path,sheet_name = 0):
#     """
#     Reads the ROI-summary Excel sheet and returns a curated wide DataFrame
#     with spine vs dendrite columns per (cell_id, roi_index).

#     Output columns:
#       cell_id, roi_index,
#       Mean_AF647_spine, Mean_AF647_dendrite,
#       BG_corrected_spine, BG_corrected_dendrite
#     """

#     raw = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)

#     # Identify the "start" columns for each cell block from row 0
#     cell_start_cols = []
#     for c in range(raw.shape[1]):
#         v = raw.iloc[0, c]
#         if isinstance(v, str) and "cell" in v.lower():
#             cell_start_cols.append(c)

#     if not cell_start_cols:
#         raise ValueError("Could not find any cell columns in row 0 (no header containing 'cell').")

#     records = []

#     for c in cell_start_cols:
#         cell_id = raw.iloc[0, c]

#         # In YOUR sheet, ROI name is in the SAME column as the cell header (row 1),
#         # and the numeric columns are to the right.
#         roi_col = c
#         mean_col = c + 1
#         bg_col = c + 2

#         # Data starts from row 2 (row 1 contains headers like 'ROI name')
#         for r in range(2, raw.shape[0]):
#             roi = raw.iloc[r, roi_col]
#             if pd.isna(roi):
#                 continue

#             m = re.match(r"^\s*(spine|dendrite)\s*[_-]\s*(\d+)\s*$", str(roi), flags=re.IGNORECASE)
#             if not m:
#                 continue

#             roi_type = m.group(1).lower()
#             roi_index = int(m.group(2))

#             records.append({
#                 "cell_id": str(cell_id).strip(),
#                 "roi_type": roi_type,
#                 "roi_index": roi_index,
#                 "Mean_AF647": raw.iloc[r, mean_col],
#                 "BG_corrected": raw.iloc[r, bg_col],
#             })

#     long_df = pd.DataFrame(records)
#     if long_df.empty:
#         raise ValueError(
#             "No spine/dendrite ROIs found. "
#             "Tip: check ROI labels (expected spine_1/dendrite_1 etc.) and that data starts at row 2."
#         )

#     # Pivot to wide (spine vs dendrite columns)
#     wide_df = (
#         long_df.pivot_table(
#             index=["cell_id", "roi_index"],
#             columns="roi_type",
#             values=["Mean_AF647", "BG_corrected"],
#             aggfunc="first"
#         )
#         .reset_index()
#     )

#     # Flatten MultiIndex columns
#     wide_df.columns = [
#         f"{a}_{b}" if b else a
#         for a, b in wide_df.columns.to_flat_index()
#     ]

#     # Nice consistent order (keep only what exists)
#     preferred = [
#         "cell_id", "roi_index",
#         "Mean_AF647_spine", "Mean_AF647_dendrite",
#         "BG_corrected_spine", "BG_corrected_dendrite",
#     ]
#     wide_df = wide_df[[c for c in preferred if c in wide_df.columns]]

#     return wide_df

# def DumpDict(fname,df):
#     df.to_pickle(fname)

def ReadDataDict(fname):
    return pd.read_pickle(fname)

def ReadData():
    # dend_cell_sum = {}
    curr_wd = os.getcwd()
    print(curr_wd)
    return [
        ReadDataDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df1_sub.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/Attila_Spine/df1_div.pickle")),
        ReadDataDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df2_sub.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/Attila_Spine/df2_div.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/MK_Spine/surf_df.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/MK_Spine/int_df.pickle"))
    ]


# df1_sub = load_spine_dendrite_dataframe(
#     "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="AF647 GluA2 mean fluor BG-subtr"
# )
# df1_div = load_spine_dendrite_dataframe(
#     "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="AF647 GluA2 mean fluor BG-divid"
# )
# df2_sub = load_spine_dendrite_dataframe(
#     "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="Xph20 PSD95 mean fluor"
# )

# df2_div = load_spine_dendrite_dataframe(
#     "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="Xph20 PSD95 mean fluor BG-divid"
# )

# curr_wd = os.getcwd()
# DumpDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df1_sub.pickle"),df1_sub)
# DumpDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df1_div.pickle"),df1_div)
# DumpDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df2_sub.pickle"),df2_sub)
# DumpDict(os.path.join(curr_wd,"Protein_data/Attila_Spine/df2_div.pickle"),df2_div)
df1_sub,df1_div,df2_sub,df2_div,surf_df,int_df=ReadData()
def add_bg_ratio(df):
    df = df.copy()
    df["SE_ratio"] = df["BG_corrected_spine"] / df["BG_corrected_dendrite"]
    return df

df1_sub = add_bg_ratio(df1_sub)
df2_sub = add_bg_ratio(df2_sub)
df1_div = add_bg_ratio(df1_div)
df2_div = add_bg_ratio(df2_div)
# breakpoint()

def GetSEvals(stat,s_df,i_df,norm_by_area):
    ratio_SE = pd.DataFrame()
    if norm_by_area == True:
        ratio_SE["surf"] = ((s_df['sp_{}'.format(stat)] - s_df['sp_bg_{}'.format(stat)]) / s_df["sp_area"]) / (
                (s_df['dend_{}'.format(stat)] - s_df['dend_bg_{}'.format(stat)]) / s_df["dend_area"])
        ratio_SE["int"] = ((i_df['sp_{}'.format(stat)] - i_df['sp_bg_{}'.format(stat)]) / s_df["sp_area"]) / (
                (i_df['dend_{}'.format(stat)] - i_df['dend_bg_{}'.format(stat)]) / s_df["dend_area"])
    else:
        ratio_SE["surf"] = ((s_df['sp_{}'.format(stat)] - s_df['sp_bg_{}'.format(stat)])) / (
            (s_df['dend_{}'.format(stat)] - s_df['dend_bg_{}'.format(stat)]))
        ratio_SE["int"] = ((i_df['sp_{}'.format(stat)] - i_df['sp_bg_{}'.format(stat)])) / (
            (i_df['dend_{}'.format(stat)] - i_df['dend_bg_{}'.format(stat)]))
    ratio_SE["id"] = s_df.index
    return ratio_SE

def calculateSE(stat,s_df,i_df,x,y,pairs,labs,order,hue_order,fname,x_pos = [0,1],norm_by_area = True):
    ratio_SE = GetSEvals(stat,s_df,i_df,norm_by_area)
    y_pos = ratio_SE[["surf", "int"]].to_numpy()
    ratio_SE_melted = ratio_SE.melt(id_vars=["id"], var_name="localization", value_name="SE")
    return ratio_SE_melted

x="localization"
y = "SE"
pairs = [("surf","int")]
labs = ["Spine","Shaft"]
order = ["surf","int"]
hue = order
hue_order = order
# breakpoint()
# breakpoint()
x_pos = [0,1]


stat = "mean"
fig_file = "./Spine_Enrichment_comparisor.svg"
MK_data = calculateSE(stat,surf_df,int_df,x,y,pairs,labs,order,hue_order,fig_file,norm_by_area=False)

plot_df = pd.concat([
    df1_sub.assign(condition="BG-subtracted", dataset="GluA2 data"),
    # df2_sub.assign(condition="BG-subtracted", dataset="Our data"),
    df1_div.assign(condition="BG-divided", dataset="GluA2 data"),
    # df2_div.assign(condition="BG-divided", dataset="Our data"),
], ignore_index=True)
plot_df_cell = (plot_df.groupby(["cell_id", "condition", "dataset"], as_index=False).agg(SE_ratio=("SE_ratio", "mean"),n_rois=("SE_ratio", "size")))

fig, ax = plt.subplots(figsize=(5, 5))


pairs = [
    (("BG-subtracted", "Getz et al., 2022"), ("BG-subtracted", "Our data")),
    (("BG-divided", "Getz et al., 2022"), ("BG-divided", "Our data")),
]

# Boxplot (background)
sns.boxplot(
    data=plot_df_cell,
    x="condition",
    y="SE_ratio",
    # hue="dataset",
    showfliers=False,   # hide outliers (points shown in strip)
    width=0.5,
    color="lightgray",
    ax = ax
)

# Stripplot (foreground)
sns.stripplot(
    data=plot_df_cell,
    x="condition",
    y="SE_ratio",
    # hue="dataset",
    dodge=True,
    jitter=0.25,
    alpha=0.7,
    size=5,
    ax=ax
)

plt.ylabel("Synaptic enrichment")
plt.xlabel("")
plt.title("Spine vs dendrite enrichment")
plt.ylim(0, 9)
# plt.show()
plt.legend(title="Dataset")
plt.tight_layout()
plt.savefig("./Spine_Enrichment_attila.svg", dpi=300)
plt.show()


plt.figure(figsize=(5, 5))
sns.boxplot(
    data=MK_data[(MK_data["SE"]>=1.0)],
    x="localization",
    y="SE",
    # hue="dataset",
    showfliers=False,   # hide outliers (points shown in strip)
    width=0.5,
    color="k"
)

# Stripplot (foreground)
sns.stripplot(
    data=MK_data[(MK_data["SE"]>=1.0)],
    x="localization",
    y="SE",
    # hue="dataset",
    dodge=True,
    jitter=0.25,
    alpha=0.7,
    size=5
)

plt.ylabel("Synaptic enrichment")
plt.xlabel("")
plt.title("Spine vs dendrite enrichment")
plt.ylim(0, 9)
# plt.show()
plt.legend(title="Dataset")
plt.tight_layout()
plt.savefig("./Spine_Enrichment_MK.svg", dpi=300)
plt.show()

MK_data = MK_data.copy() 
MK_data["cell_id"] = ( MK_data["id"] .str.extract(r"(cell_\d+)") )

new_MKdata = MK_data[MK_data["SE"] >= 1.0].copy()
order = ["Attila data", "MK data"]
MK_data_cell = ( 
                new_MKdata .groupby(["cell_id", "localization"], as_index=False) .agg( SE=("SE", "mean"), # or "median" 
                                                                                   n_spines=("SE", "size") ) )
pairs = [("Attila data", "MK data")]
# --- Attila: keep only BG-subtracted + GluA2 data (already cell-averaged) ---
att_sel = (
    plot_df_cell
    .query("condition == 'BG-subtracted' and dataset == 'GluA2 data'")
    .assign(group="Attila data", value=lambda d: d["SE_ratio"])
    [["cell_id", "group", "value"]]
)

# --- MK: keep only surf (cell-averaged) ---
mk_sel = (
    MK_data_cell
    .query("localization == 'surf'")
    .assign(group="MK data", value=lambda d: d["SE"])
    [["cell_id", "group", "value"]]
)

combined = pd.concat([att_sel, mk_sel], ignore_index=True)

fig, ax = plt.subplots(figsize=(5, 5))

order = ["Attila data", "MK data"]

sns.boxplot(
    data=combined,
    x="group",
    y="value",
    order=order,
    showfliers=False,
    width=0.5,
    color="lightgray",
    ax=ax
)

sns.stripplot(
    data=combined,
    x="group",
    y="value",
    order=order,
    jitter=0.25,
    alpha=0.7,
    size=5,
    ax=ax
)
# --- Statistical annotation ---
annot = Annotator(ax, pairs, data=combined, x="group", y="value", order=order)
annot.configure(
    test="Mann-Whitney",   # good default for ratios / non-normal
    text_format="star",    # or "simple" to print p-values
    loc="outside",
    verbose=2
)
annot.apply_and_annotate()
ax.set_ylabel("Synaptic enrichment")
ax.set_xlabel("")
ax.set_ylim(0, 9)
ax.set_title("Attila data vs MK data")
plt.tight_layout()
plt.savefig("./Enrichment_BGsub_vs_MKsurf.svg", dpi=300)
plt.show()

summary_stats = (
    combined
    .groupby("group")["value"]
    .agg(
        n="count",
        mean="mean",
        std="std",
        median="median"
    )
)

print(summary_stats)

breakpoint()
