import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Input file: all-vs-all fastANI results from Galaxy
input_file = "/mnt/c/Users/prisc/OneDrive/Desktop/GalaxyAnalyses/ani_all_vs_all.tsv"

# Read fastANI output
df = pd.read_csv(input_file, sep="\t", header=None)
df.columns = ["query", "reference", "ani", "fragments_matched", "total_fragments"]

# Clean names for plotting
def clean_name(path):
    name = str(path).split("/")[-1]

    # remove extensions
    for ext in [".fasta", ".fna", ".fa", ".zip"]:
        name = name.replace(ext, "")

    # remove fastANI suffixes
    name = name.replace("_query", "")
    name = name.replace("_ref", "")

    # botulinum -> boto
    if "botulinum" in name:
        name = name.replace("botulinum_", "")
        name = name.replace("botulinum", "")
        name = "boto " + name

    # sporogenes -> sporo
    elif "sporo" in name or "sporogenes" in name:
        name = name.replace("sporo_", "")
        name = name.replace("sporogenes_", "")
        name = name.replace("sporogenes", "")
        name = "sporo " + name

    # scaffold label
    elif "filtered_scaffold" in name:
        name = "SRR29831631"

    # general cleanup
    name = name.replace("_", " ").replace("  ", " ")
    name = " ".join(name.split())

    return name.strip()

df["query"] = df["query"].apply(clean_name)
df["reference"] = df["reference"].apply(clean_name)

# Keep highest ANI if duplicate pairs exist
df = df.groupby(["query", "reference"], as_index=False)["ani"].max()

# Build matrix
matrix = df.pivot(index="query", columns="reference", values="ani")

# Fill from reverse direction where available
matrix = matrix.combine_first(matrix.T)

# Make sure rows/columns contain same genomes
all_names = sorted(set(matrix.index).union(set(matrix.columns)))
matrix = matrix.reindex(index=all_names, columns=all_names)

# Fill diagonal
for name in matrix.index:
    matrix.loc[name, name] = 100.0

# Fill any remaining missing values with a low ANI floor
# fastANI often omits very distant pairs; 75 keeps them finite for clustering
matrix = matrix.fillna(75.0)

# Clustered heatmap
g = sns.clustermap(
    matrix,
    cmap="viridis",
    vmin=75,
    vmax=100,
    linewidths=0.5,
    figsize=(12, 12),
    annot=True,
    fmt=".1f"
)

# Make axis labels bold
for label in g.ax_heatmap.get_xticklabels():
    label.set_fontweight("bold")

for label in g.ax_heatmap.get_yticklabels():
    label.set_fontweight("bold")


g.fig.suptitle("Clustered ANI heatmap of reference genomes and SRR29831631", y=1.02)

output_file = "/mnt/c/Users/prisc/OneDrive/Desktop/GalaxyAnalyses/ani_heatmap_galaxy/ani_clustered_heatmap.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")

print(f"Saved: {output_file}")
