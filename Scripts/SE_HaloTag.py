import os
from pathlib import Path
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator
from SNSPlottingWidget import SNSPlottingWidget
pw = SNSPlottingWidget()
import numpy as np
from scipy.stats import wilcoxon

def load_spine_csv(file_path, signal):
    """
    Read the first 3 columns of the csv, assign ROI from Branch_i rows,
    remove branch rows, and add the signal column.
    """
    valid_signals = {"PSD95", "sGluA2", "iGluA2"}
    if signal not in valid_signals:
        raise ValueError(f"signal must be one of {valid_signals}")

    df = pd.read_csv(file_path, header=None, usecols=[0, 1, 2])
    df = df.dropna(how="all")
    # first row is header
    df.columns = df.iloc[0]
    df = df.iloc[1:].copy()

    first_col_name = df.columns[0]
    first_col = df[first_col_name].astype(str)

    # rows like Branch_1, Branch_2, ...
    is_branch = first_col.str.match(r"(?i)^branch_\d+$", na=False)

    # extract ROI from branch rows and forward-fill
    df["ROI"] = first_col.str.extract(r"(?i)^branch_(\d+)$")[0]
    df["ROI"] = df["ROI"].ffill()

    # remove branch rows themselves
    df = df[~is_branch].copy()

    # remove anything before first branch
    df = df[df["ROI"].notna()].copy()

    # Give each spine an index within ROI
    df["spine_id"] = df.groupby("ROI").cumcount() + 1
    
    df["ROI"] = df["ROI"].astype(int)
    df["signal"] = signal
    df = df.reset_index(drop=True)
    
    for col in ["PSD_intensity_BG_corr", "Dendritic_intensity_BG_corr","Spine_Dendrite_Ratio"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df

def load_area_csv(file_path):
    """
    Read the first 2 columns of the csv, assign ROI from Branch_i rows,
    remove branch rows, and add the signal column.
    """

    df = pd.read_csv(file_path, header=None, usecols=[0, 1])

    # first row is header
    df.columns = df.iloc[0]
    df = df.iloc[1:].copy()

    first_col_name = df.columns[0]
    first_col = df[first_col_name].astype(str)

    # rows like Branch_1, Branch_2, ...
    is_branch = first_col.str.match(r"(?i)^branch_\d+$", na=False)

    # extract ROI from branch rows and forward-fill
    df["ROI"] = first_col.str.extract(r"(?i)^branch_(\d+)$")[0]
    df["ROI"] = df["ROI"].ffill()

    # remove branch rows themselves
    df = df[~is_branch].copy()

    # remove anything before first branch
    df = df[df["ROI"].notna()].copy()

    # Give each spine an index within ROI
    df["spine_id"] = df.groupby("ROI").cumcount() + 1
    
    df["ROI"] = df["ROI"].astype(int)
    # df["signal"] = signal
    df = df.reset_index(drop=True)
    
    for col in ["Spine_Area", "Dendrite_Area"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df

def extract_roi_metadata(file_path, signal):
    """
    Extract ROI-level metadata and attach signal label.

    Returns:
        pd.DataFrame with:
        ROI, Branch_length, Background_value, PSD_number, Area, signal
    """

    valid_signals = {"PSD95", "sGluA2", "iGluA2"}
    if signal not in valid_signals:
        raise ValueError(f"signal must be one of {valid_signals}")

    df = pd.read_csv(file_path, header=None)

    # first row is header
    df.columns = df.iloc[0]
    df = df.iloc[1:].copy()

    # select rows where metadata exists
    meta = df[df["Branch_ID"].notna()].copy()

    # extract ROI
    meta["ROI"] = meta["Branch_ID"].astype(str).str.extract(r"(?i)^branch_(\d+)$")[0]
    meta = meta[meta["ROI"].notna()].copy()
    meta["ROI"] = meta["ROI"].astype(int)

    # keep relevant columns
    meta = meta[["ROI", "Branch_length", "Background_value"]].copy()
    
    # convert to numeric
    for col in ["Branch_length", "Background_value"]:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")

    # add signal label 
    meta["signal"] = signal

    meta = meta.drop_duplicates(subset="ROI").reset_index(drop=True)

    return meta



def build_master_dataframe(max1_folder,max_f = "3"):
    """
    Iterate through MAX1/cell_*/spines/*.csv,
    assign files to PSD95 / sGluA2 / iGluA2,
    and combine everything into one big dataframe.
    """
    max1_folder = Path(max1_folder)
    all_dfs = []
    all_meta = []
    all_areas = []
    for cell_dir in sorted(max1_folder.glob("cell_*")):
        spines_dir = cell_dir / "spines"
        if not spines_dir.exists():
            continue

        csv_files = list(spines_dir.glob("*.csv"))

        c1 = next((f for f in csv_files if f.name.endswith("_c1.csv")), None)
        c2 = next((f for f in csv_files if f.name.endswith("_c2.csv")), None)
        c3 = next((f for f in csv_files if f.name.endswith("_c3.csv")), None)
        area_file = next((f for f in csv_files if f.name.endswith("_areas.csv")), None)
        # breakpoint()
        if not all([c1, c2, c3]):
            print(f"Skipping {cell_dir.name}: missing one or more csv files")
            continue

        # assignment rule
        jf549i_present = any("JF549i" in f.name for f in [c1, c2, c3])

        file_map = {
            "PSD95": c3,
            "sGluA2": c2 if jf549i_present else c1,
            "iGluA2": c1 if jf549i_present else c2,
        }


        for signal, file_path in file_map.items():
            # spine-level data
            df = load_spine_csv(file_path, signal)
            df["cell"] = cell_dir.name + "_{}".format(max_f)
            df["file"] = file_path.name
            all_dfs.append(df)
            # print(f"Loaded {len(df)} rows for {signal} from {cell_dir.name,file_path.name}")
            # ROI metadata WITH signal
            meta_df = extract_roi_metadata(file_path, signal)
            meta_df["cell"] = cell_dir.name + "_{}".format(max_f)
            all_meta.append(meta_df)
        
        areas = load_area_csv(area_file)
        areas["cell"] = cell_dir.name + "_{}".format(max_f)
        all_areas.append(areas)
        # print(f"Loaded {len(areas)} rows for areas from {cell_dir.name,area_file.name}")
    master_df = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
    meta_df = pd.concat(all_meta, ignore_index=True) if all_meta else pd.DataFrame()
    all_area_df = pd.concat(all_areas, ignore_index=True) if all_areas else pd.DataFrame()

    return master_df, meta_df,all_area_df

master_df3, roi_meta_df3,all_area_df3 = build_master_dataframe("./Protein/MAX3","3")
master_df2, roi_meta_df2,all_area_df2 = build_master_dataframe("./Protein/MAX2","2")
master_df1, roi_meta_df1,all_area_df1 = build_master_dataframe("./Protein/MAX1","1")

master_df = pd.concat([master_df1,master_df2, master_df3], ignore_index=True)
roi_meta_df = pd.concat([roi_meta_df1,roi_meta_df2, roi_meta_df3], ignore_index=True)
all_area_df = pd.concat([all_area_df1,all_area_df2,all_area_df3], ignore_index=True)
print("master_df rows before merge:", len(master_df))
print("roi_meta_df rows:", len(roi_meta_df))
print("all_area_df rows:", len(all_area_df))
# breakpoint()
print(
    roi_meta_df.groupby(["cell", "ROI", "signal"]).size().max()
)

merged_df = master_df.merge(
    roi_meta_df,
    on=["cell", "ROI", "signal"],
    how="left",
    validate="many_to_one"
)
merged_df = merged_df.merge(
    all_area_df,
    on=["cell", "ROI", "spine_id"],
    how="left",
    validate="many_to_one"
)
breakpoint()

# breakpoint()
print("merged_df rows after merge:", len(merged_df))
sig_per =  0
sglua2_pass = merged_df[
    (merged_df["signal"] == "sGluA2") &
    (merged_df["PSD_intensity_BG_corr"].notna()) &
    (merged_df["Background_value"].notna()) &
    (merged_df["PSD_intensity_BG_corr"] >= merged_df["Background_value"] * (sig_per/100.))
].copy()

valid_spines = sglua2_pass[["cell", "ROI", "spine_id"]].drop_duplicates()

filtered_df = merged_df.merge(
    valid_spines,
    on=["cell", "ROI", "spine_id"],
    how="inner",
    validate="many_to_many"
)
filtered_df = filtered_df[
    filtered_df["signal"].isin(["sGluA2", "iGluA2", "PSD95"])
].copy()
# breakpoint()
filtered_df = filtered_df[
    (filtered_df["PSD_intensity_BG_corr"].notna()) &
    (filtered_df["Dendritic_intensity_BG_corr"].notna()) &
    (filtered_df["PSD_intensity_BG_corr"] >= 0) &
    (filtered_df["Dendritic_intensity_BG_corr"] >= 0)
].copy()
# breakpoint()
filtered_df.groupby(["cell", "ROI", "spine_id"])["signal"].nunique().value_counts()
print("Original rows:", len(merged_df))
print("Filtered rows:", len(filtered_df))
print("Unique spines:", len(valid_spines))
# pivot_df = filtered_df.pivot_table(
#     index=["cell", "ROI", "spine_id"],
#     columns="signal",
#     values="Spine_Dendrite_Ratio"
# ).reset_index()
# pivot_df = df.melt(
#     id_vars=["cell_id"],                 # keep identifier columns
#     value_vars=["PSD95", "sGluA2", "iGluA2"],  # columns to unpivot
#     var_name="signal",                  # new column for names
#     value_name="SE_ratio"               # new column for values
# )
# breakpoint()
fsize = 26
lsize= 26
x = "signal"
y = "Spine_Dendrite_Ratio"
hue = None#["Soma","Dendrite"]
order = ["sGluA2", "iGluA2", "PSD95"]
xlab = "Protein"
ylab = "Synaptic Enrichment"
title = ""
stat_test = "Wilcoxon"
op_folder = "./Figures/HaloGluA2_analysis/"
fig_file = os.path.join(op_folder,"SE_iGLuA2_sGluA2_PSD95_{}".format(sig_per))
pairs = [("sGluA2","iGluA2"), ("sGluA2","PSD95"), ("iGluA2","PSD95")]
my_pal1 = {"sGluA2":"#FF00FF","iGluA2": '#00FFFF',"PSD95": '#008000'}
ax_label=1
# x_pos = [0,1]
# y_pos = np.stack((cell_avg_paired1.groupby("compartment")["surf_to_int_ratio"].apply(list)[::-1])).T
cell_avg = (
    filtered_df
    .groupby(["cell", "signal"])[["Spine_Dendrite_Ratio","Spine_Area","Dendrite_Area"]]
    .mean()
    .reset_index()
)
# breakpoint()
cols = ['cell', 'Spine_Dendrite_Ratio']
df_subset = cell_avg.loc[cell_avg['signal'] == 'sGluA2', cols].copy()
df_subset = df_subset[cols]
df_subset.to_pickle('./Scripts/Protein_data/Halo_data_spines/sGluA2_surface_only_ratios.pkl')

pw.Swarmboxplotcombo(data=cell_avg,
                     fig_size = (8,6),
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab,
                             pairs=pairs,
                             order = order,
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
fig_file = os.path.join(op_folder,"PSD_roi_area")
pw.Swarmboxplotcombo(data=cell_avg,
                     fig_size = (8,6),
                             x=x,
                             y="Spine_Area",
                             hue=hue,
                             xlab=xlab,
                             ylab=r"Spine Area ($\mu m^2$)",
                             pairs=pairs,
                             order = order,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 0
                             )
fig_file = os.path.join(op_folder,"Dendritic_roi_area")
pw.Swarmboxplotcombo(data=cell_avg,
                     fig_size = (8,6),
                             x=x,
                             y="Dendrite_Area",
                             hue=hue,
                             xlab=xlab,
                             ylab=r"Dendritic Area ($\mu m^2$)",
                             pairs=pairs,
                             order = order,
                             title=title,
                             color_pal=my_pal1,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 0
                             )
print(cell_avg.groupby("signal").agg(['min', 'max','mean', 'median','std']))
# breakpoint()
# import pandas as pd

# df = your dataframe

# --- 1. Spine count per ROI per cell ---
spine_counts_roi = (
    filtered_df.groupby(['cell', 'ROI'])['spine_id']
    .nunique()
    .reset_index(name='n_spines')
)

# --- 2. Get branch length per ROI ---
# (assuming Branch_length is same for all rows within an ROI)
branch_lengths_roi = (
    filtered_df.groupby(['cell', 'ROI'])['Branch_length']
    .first()
    .reset_index()
)

# --- 3. Merge counts + lengths ---
roi_density1 = spine_counts_roi.merge(
    branch_lengths_roi,
    on=['cell', 'ROI']
)

# --- 4. Compute density per ROI ---
roi_density1['spine_density'] = (
    roi_density1['n_spines'] / roi_density1['Branch_length']
)

# ================================
# 📊 PER CELL DENSITY
# ================================

# --- 5. Total spines per cell ---
# Sum spines across ROIs within each cell
spine_counts_cell = (
    spine_counts_roi.groupby('cell')['n_spines']
    .sum()
    .reset_index(name='total_spines')
)

# --- 6. Total dendritic length per cell ---

branch_lengths_cell = (
    branch_lengths_roi.groupby('cell')['Branch_length']
    .sum()
    .reset_index(name='total_length')
)

# --- 7. Merge + compute density ---
cell_density1 = spine_counts_cell.merge(
    branch_lengths_cell,
    on='cell'
)

cell_density1['spine_density'] = (
    cell_density1['total_spines'] / cell_density1['total_length']
)
# ================================
# Outputs
# ================================
print("ROI-level density:")
print(roi_density1.head())

print("\nCell-level density:")
print(cell_density1.head())

print("cell_density = ",cell_density1["spine_density"].mean(),cell_density1["spine_density"].std(),len(cell_density1))
print("ROI density = ",roi_density1["spine_density"].mean(),roi_density1["spine_density"].std(),len(roi_density1))
# breakpoint()


# breakpoint()

# --- 1. Spine count per ROI per cell ---
spine_counts_roi = (
    merged_df.groupby(['cell', 'ROI'])['spine_id']
    .nunique()
    .reset_index(name='n_spines')
)

# --- 2. Get branch length per ROI ---
# (assuming Branch_length is same for all rows within an ROI)
branch_lengths_roi = (
    merged_df.groupby(['cell', 'ROI'])['Branch_length']
    .first()
    .reset_index()
)

# --- 3. Merge counts + lengths ---
roi_density2 = spine_counts_roi.merge(
    branch_lengths_roi,
    on=['cell', 'ROI']
)

# --- 4. Compute density per ROI ---
roi_density2['spine_density'] = (
    roi_density2['n_spines'] / roi_density2['Branch_length']
)

# ================================
# 📊 PER CELL DENSITY
# ================================

# --- 5. Total spines per cell ---
# Sum spines across ROIs within each cell
spine_counts_cell = (
    spine_counts_roi.groupby('cell')['n_spines']
    .sum()
    .reset_index(name='total_spines')
)
# --- 6. Total dendritic length per cell ---
# (sum across ROIs)
branch_lengths_cell = (
    branch_lengths_roi.groupby('cell')['Branch_length']
    .sum()
    .reset_index(name='total_length')
)

# --- 7. Merge + compute density ---
cell_density2 = spine_counts_cell.merge(
    branch_lengths_cell,
    on='cell'
)

cell_density2['spine_density'] = (
    cell_density2['total_spines'] / cell_density2['total_length']
)

print("cell_density = ",cell_density2["spine_density"].mean(),cell_density2["spine_density"].std(),len(cell_density2))
print("ROI density = ",roi_density2["spine_density"].mean(),roi_density2["spine_density"].std(),len(roi_density2))

roi_density1["dataset"] = "After"
roi_density2["dataset"] = "Before"
cell_density1["dataset"] = "After"
cell_density2["dataset"] = "Before"
cell_density_paired = cell_density2[['cell', 'spine_density']].merge(
    cell_density1[['cell', 'spine_density']],
    on='cell', how='inner', suffixes=('_before', '_after'), validate='one_to_one'
).sort_values('cell')
if not np.isfinite(cell_density_paired[['spine_density_before', 'spine_density_after']].to_numpy()).all():
    raise ValueError('Non-finite paired cell densities')
cell_density_differences = (
    cell_density_paired['spine_density_before'].to_numpy()
    - cell_density_paired['spine_density_after'].to_numpy()
)
cell_density_test = None
if np.any(cell_density_differences != 0):
    cell_density_test = wilcoxon(
        cell_density_differences, alternative='two-sided',
        zero_method='wilcox', correction=False
    )
    print(
        'S15B: two-sided paired Wilcoxon signed-rank test; '
        f'n={len(cell_density_paired)} cell pairs; '
        f'nonzero differences={np.count_nonzero(cell_density_differences)}; '
        f'W={cell_density_test.statistic:g}; P={cell_density_test.pvalue:.10g}'
    )
else:
    print('S15B: all paired cell densities are identical; Wilcoxon test not calculated')
# Preserve identical cell order in the two groups passed to statannotations.
cell_density_combined = cell_density_paired.melt(
    id_vars='cell', value_vars=['spine_density_before', 'spine_density_after'],
    var_name='dataset', value_name='spine_density'
)
cell_density_combined['dataset'] = cell_density_combined['dataset'].map({
    'spine_density_before': 'Before', 'spine_density_after': 'After'
})
roi_density_paired = roi_density2[['cell', 'ROI', 'spine_density']].merge(
    roi_density1[['cell', 'ROI', 'spine_density']],
    on=['cell', 'ROI'], how='inner', suffixes=('_before', '_after'),
    validate='one_to_one'
).sort_values(['cell', 'ROI'])
if not np.isfinite(roi_density_paired[['spine_density_before', 'spine_density_after']].to_numpy()).all():
    raise ValueError('Non-finite paired dendrite densities')
roi_density_differences = (
    roi_density_paired['spine_density_before'].to_numpy()
    - roi_density_paired['spine_density_after'].to_numpy()
)
roi_density_test = None
if np.any(roi_density_differences != 0):
    roi_density_test = wilcoxon(
        roi_density_differences, alternative='two-sided',
        zero_method='wilcox', correction=False
    )
    print(
        'S15C: exploratory two-sided paired Wilcoxon signed-rank test; '
        f'n={len(roi_density_paired)} dendrite pairs; '
        f'nonzero differences={np.count_nonzero(roi_density_differences)}; '
        f'W={roi_density_test.statistic:g}; P={roi_density_test.pvalue:.10g}. '
        'This test does not account for dendrites nested within cells.'
    )
else:
    print('S15C: all paired dendrite densities are identical; Wilcoxon test not calculated')
roi_density_combined = roi_density_paired.melt(
    id_vars=['cell', 'ROI'],
    value_vars=['spine_density_before', 'spine_density_after'],
    var_name='dataset', value_name='spine_density'
)
roi_density_combined['dataset'] = roi_density_combined['dataset'].map({
    'spine_density_before': 'Before', 'spine_density_after': 'After'
})
stat_test = 'Wilcoxon' if cell_density_test is not None else None
x="dataset"
y = "spine_density"
xlab = "filtering"
color_pal = {"Before": "#FF00FF", "After": "#008000"}
ylab = r"Spine Density (spines/$\mu m$)"
pairs = [("Before", "After")]
order = ["Before", "After"]
fig_file = os.path.join(op_folder,"cell_density_comparison")
pw.Swarmboxplotcombo(data=cell_density_combined,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab,
                             pairs=pairs,
                             order = order,
                             title=title,
                             color_pal=color_pal,
                             stat_test=stat_test,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 0
                             )
fig_file = os.path.join(op_folder,"roi_density_comparison")
# Exploratory S15C comparison; the test does not account for cell clustering.
pw.Swarmboxplotcombo(data=roi_density_combined,
                             x=x,
                             y=y,
                             hue=hue,
                             xlab=xlab,
                             ylab=ylab,
                             pairs=pairs,
                             order = order,
                             title=title,
                             color_pal=color_pal,
                             stat_test='Wilcoxon' if roi_density_test is not None else None,
                             xfsize=fsize,
                             yfsize=1.2 * fsize,
                             fname=fig_file,
                             save_it=1,
                             ax_lab=ax_label,
                             log_y = 0,
                             circle_size= 8
                             )
breakpoint()
