import pandas as pd
import re
import seaborn as sns
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

sns.set(style="whitegrid")



def load_spine_dendrite_dataframe(excel_path,sheet_name = 0):
    """
    Reads the ROI-summary Excel sheet and returns a curated wide DataFrame
    with spine vs dendrite columns per (cell_id, roi_index).

    Output columns:
      cell_id, roi_index,
      Mean_AF647_spine, Mean_AF647_dendrite,
      BG_corrected_spine, BG_corrected_dendrite
    """

    raw = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)

    # Identify the "start" columns for each cell block from row 0
    cell_start_cols = []
    for c in range(raw.shape[1]):
        v = raw.iloc[0, c]
        if isinstance(v, str) and "cell" in v.lower():
            cell_start_cols.append(c)

    if not cell_start_cols:
        raise ValueError("Could not find any cell columns in row 0 (no header containing 'cell').")

    records = []

    for c in cell_start_cols:
        cell_id = raw.iloc[0, c]

        # In YOUR sheet, ROI name is in the SAME column as the cell header (row 1),
        # and the numeric columns are to the right.
        roi_col = c
        mean_col = c + 1
        bg_col = c + 2

        # Data starts from row 2 (row 1 contains headers like 'ROI name')
        for r in range(2, raw.shape[0]):
            roi = raw.iloc[r, roi_col]
            if pd.isna(roi):
                continue

            m = re.match(r"^\s*(spine|dendrite)\s*[_-]\s*(\d+)\s*$", str(roi), flags=re.IGNORECASE)
            if not m:
                continue

            roi_type = m.group(1).lower()
            roi_index = int(m.group(2))

            records.append({
                "cell_id": str(cell_id).strip(),
                "roi_type": roi_type,
                "roi_index": roi_index,
                "Mean_AF647": raw.iloc[r, mean_col],
                "BG_corrected": raw.iloc[r, bg_col],
            })

    long_df = pd.DataFrame(records)
    if long_df.empty:
        raise ValueError(
            "No spine/dendrite ROIs found. "
            "Tip: check ROI labels (expected spine_1/dendrite_1 etc.) and that data starts at row 2."
        )

    # Pivot to wide (spine vs dendrite columns)
    wide_df = (
        long_df.pivot_table(
            index=["cell_id", "roi_index"],
            columns="roi_type",
            values=["Mean_AF647", "BG_corrected"],
            aggfunc="first"
        )
        .reset_index()
    )

    # Flatten MultiIndex columns
    wide_df.columns = [
        f"{a}_{b}" if b else a
        for a, b in wide_df.columns.to_flat_index()
    ]

    # Nice consistent order (keep only what exists)
    preferred = [
        "cell_id", "roi_index",
        "Mean_AF647_spine", "Mean_AF647_dendrite",
        "BG_corrected_spine", "BG_corrected_dendrite",
    ]
    wide_df = wide_df[[c for c in preferred if c in wide_df.columns]]

    return wide_df

df1_sub = load_spine_dendrite_dataframe(
    "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="AF647 GluA2 mean fluor BG-subtr"
)
df1_div = load_spine_dendrite_dataframe(
    "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="AF647 GluA2 mean fluor BG-divid"
)
df2_sub = load_spine_dendrite_dataframe(
    "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="Xph20 PSD95 mean fluor"
)

df2_div = load_spine_dendrite_dataframe(
    "../Attila/AI GluA2 surface intensity by Fabs/AI GluA2 manual ROI placement/Manual ROI analysis summary.xlsx",sheet_name ="Xph20 PSD95 mean fluor BG-divid"
)

def add_bg_ratio(df):
    df = df.copy()
    df["SE_ratio"] = df["BG_corrected_spine"] / df["BG_corrected_dendrite"]
    return df

df1_sub = add_bg_ratio(df1_sub)
df2_sub = add_bg_ratio(df2_sub)
df1_div = add_bg_ratio(df1_div)
df2_div = add_bg_ratio(df2_div)


plot_df = pd.concat([
    df1_sub.assign(condition="BG-subtracted", dataset="df1"),
    df2_sub.assign(condition="BG-subtracted", dataset="df2"),
    df1_div.assign(condition="BG-divided", dataset="df1"),
    df2_div.assign(condition="BG-divided", dataset="df2"),
], ignore_index=True)

plt.figure(figsize=(5, 5))
# Boxplot (background)
sns.boxplot(
    data=plot_df,
    x="condition",
    y="SE_ratio",
    hue="dataset",
    showfliers=False,   # hide outliers (points shown in strip)
    width=0.5,
    color="lightgray"
)

# Stripplot (foreground)
sns.stripplot(
    data=plot_df,
    x="condition",
    y="SE_ratio",
    hue="dataset",
    dodge=True,
    jitter=0.25,
    alpha=0.7,
    size=5
)

plt.ylabel("BG-corrected spine / dendrite ratio")
plt.xlabel("")
plt.title("Spine vs dendrite enrichment")

plt.legend(title="Dataset")
plt.tight_layout()
plt.show()


breakpoint()