import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

ani_file = "ani_all_vs_all.tsv"

df = pd.read_csv(ani_file, sep="\t", header=None)
df.columns = ["query", "reference", "ani", "fragments", "total_fragments"]

def clean_name(path):
    name = os.path.basename(str(path))
    name = name.replace(".fasta", "")
    name = name.replace(".fna", "")
    name = name.replace("_query", "")
    name = name.replace("_ref", "")
    name = name.replace("filtered_scaffold", "SRR29831631")
    return name

df["query_clean"] = df["query"].apply(clean_name)
df["reference_clean"] = df["reference"].apply(clean_name)

genomes = sorted(set(df["query_clean"]).union(set(df["reference_clean"])))

ani_matrix = pd.DataFrame(np.nan, index=genomes, columns=genomes)

for genome in genomes:
    ani_matrix.loc[genome, genome] = 100.0

for _, row in df.iterrows():
    ani_matrix.loc[row["query_clean"], row["reference_clean"]] = row["ani"]

for i in genomes:
    for j in genomes:
        a = ani_matrix.loc[i, j]
        b = ani_matrix.loc[j, i]

        if pd.isna(a) and not pd.isna(b):
            ani_matrix.loc[i, j] = b
        elif pd.isna(b) and not pd.isna(a):
            ani_matrix.loc[j, i] = a
        elif not pd.isna(a) and not pd.isna(b):
            mean_value = (a + b) / 2
            ani_matrix.loc[i, j] = mean_value
            ani_matrix.loc[j, i] = mean_value

ani_matrix = ani_matrix.fillna(75)

distance_matrix = 100 - ani_matrix

pca = PCA(n_components=3)
coords = pca.fit_transform(distance_matrix)

pca_df = pd.DataFrame(coords, columns=["PC1", "PC2", "PC3"])
pca_df["genome"] = genomes

def assign_group(name):
    lower = name.lower()
    if "srr29831631" in lower:
        return "SRR29831631"
    if "sporo" in lower or "sporogenes" in lower:
        return "C. sporogenes"
    if "tagluense" in lower:
        return "Outgroup"
    if "botulinum" in lower or "boto" in lower:
        return "C. botulinum"
    return "Other"

pca_df["group"] = pca_df["genome"].apply(assign_group)

plt.figure(figsize=(8, 6))

for group in pca_df["group"].unique():
    sub = pca_df[pca_df["group"] == group]
    plt.scatter(sub["PC1"], sub["PC2"], label=group, s=70)

for _, row in pca_df.iterrows():
    if row["group"] == "SRR29831631":
        plt.text(row["PC1"], row["PC2"], row["genome"], fontsize=10, fontweight="bold")
    else:
        plt.text(row["PC1"], row["PC2"], row["genome"], fontsize=7)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA of ANI distances")
plt.legend()
plt.tight_layout()
plt.savefig("ani_pca_2d.png", dpi=300)

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")

for group in pca_df["group"].unique():
    sub = pca_df[pca_df["group"] == group]
    ax.scatter(sub["PC1"], sub["PC2"], sub["PC3"], label=group, s=70)

for _, row in pca_df.iterrows():
    if row["group"] == "SRR29831631":
        ax.text(row["PC1"], row["PC2"], row["PC3"], row["genome"], fontsize=10, fontweight="bold")
    else:
        ax.text(row["PC1"], row["PC2"], row["PC3"], row["genome"], fontsize=7)

ax.set_title("3D PCA of ANI distances")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_zlabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)")
ax.legend()

plt.tight_layout()
plt.savefig("ani_pca_3d.png", dpi=300)

ani_matrix.to_csv("ani_matrix_for_pca.csv")
pca_df.to_csv("pca_coordinates.csv", index=False)

print("Saved ani_pca_2d.png")
print("Saved ani_pca_3d.png")
print("Saved ani_matrix_for_pca.csv")
print("Saved pca_coordinates.csv")
print("\nGenomes included:")
print(pca_df[["genome", "group"]])
