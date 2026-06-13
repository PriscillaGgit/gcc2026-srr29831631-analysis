import pandas as pd
import matplotlib.pyplot as plt

# Input file
input_file = "/mnt/c/Users/prisc/OneDrive/Desktop/GalaxyAnalyses/ani_results.tsv"

# Read fastANI output
df = pd.read_csv(input_file, sep="\t", header=None)

# fastANI has 5 columns
df.columns = ["query", "reference", "ani", "fragments_matched", "total_fragments"]

# Clean reference names
def clean_name(path):
    name = path.split("/")[-1]

    # remove extensions
    for ext in [".fasta", ".fna", ".fa", ".zip"]:
        name = name.replace(ext, "")

    # botulinum → boto
    if "botulinum" in name:
        name = name.replace("botulinum_", "")
        name = name.replace("botulinum", "")
        name = "boto " + name

    # sporogenes → sporo
    elif "sporo" in name or "sporogenes" in name:
        name = name.replace("sporo_", "")
        name = name.replace("sporogenes_", "")
        name = name.replace("sporogenes", "")
        name = "sporo " + name

    # remove generic junk
    name = name.replace("_ref", "")
    name = name.replace("ref", "")

    # clean formatting
    name = name.replace("_", " ")
    name = " ".join(name.split())

    return name.strip()

df["label"] = df["reference"].apply(clean_name)

# Sort by ANI (correct biological ordering)
df = df.sort_values("ani", ascending=False)

# Build 1-row heatmap matrix
heatmap_df = pd.DataFrame(
    [df["ani"].values],
    columns=df["label"].values,
    index=["SRR29831631"]
)

# Figure size scales with number of references
width = max(12, len(heatmap_df.columns) * 0.6)
fig, ax = plt.subplots(figsize=(width, 3))

# Plot heatmap
im = ax.imshow(heatmap_df.values, aspect="auto", vmin=90, vmax=100)

# Axis labels
ax.set_xticks(range(len(heatmap_df.columns)))
ax.set_xticklabels(heatmap_df.columns, rotation=60, ha="right", fontsize=9)
ax.set_yticks([0])
ax.set_yticklabels(["SRR29831631"], fontsize=10, fontweight="bold")

# Add ANI values inside boxes
for j, value in enumerate(heatmap_df.iloc[0]):
    color = "white" if value < 93 else "black"
    ax.text(j, 0, f"{value:.1f}", ha="center", va="center", fontsize=8, color=color)

# Title
ax.set_title("ANI heatmap of SRR29831631 against reference genomes")

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("ANI (%)")

# Layout
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save output
output_file = "/mnt/c/Users/prisc/OneDrive/Desktop/GalaxyAnalyses/ani_heatmap_galaxy/ani_heatmap.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")

print(f"Saved: {output_file}")
