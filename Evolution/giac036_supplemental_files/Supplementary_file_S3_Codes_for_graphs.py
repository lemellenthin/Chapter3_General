#!/usr/bin/env python
# coding: utf-8

# # The state of Medusozoa genomics: current evidence and future challanges
# ## Supplementary file 
# This notebook was used to construct plots from Santander et al. 2021.  
# Graphs were further modified to improve visualization using Inkskape 0.92.
# 
# This code was written in Python 3.6.13
# This code is open under GNU General Public License GPLv3 
# (licence details available at https://www.gnu.org/licenses/gpl-3.0.en.html):

# In[ ]:


#Dependencies
import os
import re
from collections import Counter

import matplotlib.pyplot as plt # version 3.3.4
import numpy as np #version 1.19.5 
import pandas as pd # version 1.1.5
import seaborn as sns #version 0.11.
from ete3 import (
    AttrFace,
    CircleFace,
    ClusterTree,
    FaceContainer,
    ImgFace,
    NodeStyle,
    PhyloTree,
    RectFace,
    StackedBarFace,
    TextFace,
    Tree,
    TreeStyle,
    faces,
) #version 3.1.2


# ## FIGURE 1

# In[ ]:


# load tree
t = PhyloTree("Data/Medusozoa.tre", quoted_node_names=True, format=1)
# load information for graph
stats = pd.read_csv("Data/Figure1_scriptData.csv", sep="\t")
stats = stats.fillna(0)


# In[ ]:


# define a factor for normalization for bar size
# for Genome size
factor = stats["Genome size max."].max()
# for species number
factor2 = stats.SpeciesN.max()


# normalize
stats["SeqSpecies"] = stats["Sequenced species"].div(factor2)
stats["kariotyped species"] = stats["kariotyped species"].div(factor2)
stats["GsN"] = stats["Genome size number"].div(factor2)
stats["SpeciesN"] = stats["SpeciesN"].div(factor2)
stats["Genome size max."] = stats["Genome size max."].div(factor)
stats["Genome size min."] = stats["Genome size min."].div(factor)
stats["Allinfo"] = stats["Allinfo"].div(factor2)
stats = stats.set_index("Family")
stats["Allinfo"] = stats["Allinfo"].div(factor2)


# In[ ]:


# definition of functions for tree plotting


def addTreeFace(
    node_name, number, type_graph, column, color, function=False, factor=None
):
    node = t.get_leaves_by_name(node_name)
    # if the variable is a float or an int
    # this code graphs the variable taking into account
    # the factor used for normalization

    if type(number) in [float, int]:
        if number > 0:
            orinumber = number
            if function:
                number = number * 6
                orinumber = orinumber * factor
            if type_graph == "rect":
                graph = RectFace(
                    width=number * 30,
                    height=20,
                    fgcolor=color,
                    bgcolor=color,
                    label={
                        "text": str(int(orinumber)),
                        "font": "inter",
                        "color": "black",
                        "fontsize": 14,
                    },
                )
            elif type_graph == "circle":
                graph = CircleFace()
            if len(node) >= 1:
                node[0].add_face(graph, column=column, position="aligned")
            else:
                print(node_name, " Face not found: ", number)
        else:
            print(node_name + " incorrect value or missing value: ", number)
    elif type_graph == "str":
        graph = TextFace(number, fsize=14, ftype="inter")
        if len(node) >= 1:
            node[0].add_face(graph, column=column, position="aligned")
        else:
            print(node_name, "Incorrect Str or terminal not found: ", number)


# Start plotting
# plot Total species in first column
stats.apply(
    lambda x: addTreeFace(
        x.Family, x["SpeciesN"], "rect", 1, "#ffb000", function=True, factor=factor2
    ),
    axis=1,
)
# plot total number of species with at least  1 information
stats.apply(
    lambda x: addTreeFace(
        x.Family, x["Allinfo"], "rect", 2, "#000000", function=True, factor=factor2
    ),
    axis=1,
)

# plot text of Chromosome number range in third column
stats.apply(
    lambda x: addTreeFace(x.Family, x["Chromosome number range"], "str", 3, "#fe6100"),
    axis=1,
)

# Plot genome size min for each family in fourth column
stats.apply(
    lambda x: addTreeFace(
        x.Family,
        x["Genome size min."],
        "rect",
        4,
        "#dc267f",
        function=True,
        factor=factor,
    ),
    axis=1,
)
# Plot genome size max for each family in column 5
stats.apply(
    lambda x: addTreeFace(
        x.Family,
        x["Genome size max."],
        "rect",
        5,
        "#648fff",
        function=True,
        factor=factor,
    ),
    axis=1,
)

# Plot genome size range (text) for each family in column 6
stats.apply(
    lambda x: addTreeFace(x.Family, str(x["Genome size number"]), "str", 6, "#648fff"),
    axis=1,
)
# Plot number of unique genomes in column 7
stats.apply(
    lambda x: addTreeFace(x.Family, x.Unique_genomes, "rect", 7, "#fe6100"), axis=1
)
# plot total number of sequenced species in column 8
stats.apply(
    lambda x: addTreeFace(x.Family, x["Sequenced species"], "rect", 8, "#785ef0"),
    axis=1,
)

# add family name
for n in t.traverse():
    if not n.is_leaf():
        n.add_face(
            TextFace(n.name, fsize=14, ftype="inter"), column=0, position="float"
        )
    if n.is_leaf():
        N = AttrFace("name", fsize=14)
        n.add_face(N, 0, position="aligned")


# change tree appearence
t.convert_to_ultrametric()
t.ladderize(direction=2)

# prevent using labels and symbols in internal nodes when drawing tree
for n in t.traverse():
    n.img_style["size"] = 0
    n.img_style["vt_line_width"] = 2
    n.img_style["hz_line_width"] = 2

# define treestyle
ts = TreeStyle()
ts.min_leaf_separation = 30
ts.show_scale = False
ts.force_topology = False

t.show(tree_style=ts)
t.render(
    "figuras/tree_genomes.svg",
    tree_style=ts,
    units="px",
    h=500,
    dpi=300,
)


# ## Figure 2
# ### Figure 2A

# In[ ]:


# Read csv with BBmap results and csv with biogenome score manually calculated

df = pd.read_csv("Data/BBMAP_statistics_for_graph.csv", sep="\t", decimal=",")
bdf = pd.read_csv("Data/Biogenome_score.csv", sep="\t")
df = df.set_index("Species")

# modify statistics values for better visualization
# (i.e. bp to Mbp and counts in thousands)
df["Scaffold_Mbp"] = df.Scaffold_bp / 100000
df["Contig_n"] = df.Contig_n / 1000
df["Scaffold_n"] = df.Scaffold_n / 1000
df["Contig_N50"] = df.Contig_N50 / 1000
df["Scaffold_N50"] = df.Scaffold_N50 / 1000

# sort statistics by class and Scaffold N50
order = {"Hydrozoa": 0, "Staurozoa": 1, "Cubozoa": 2, "Scyphozoa": 3}
df["order"] = df.Class.apply(lambda x: order[x])
df = df.sort_values(["order", "Scaffold_N50"])

# Plot fig 2A.
# Create grid
fig = plt.figure(figsize=(18, 5))
ax1 = plt.subplot(1, 5, 1)
ax2 = plt.subplot(1, 5, 2)
ax3 = plt.subplot(4, 5, 3)
ax4 = plt.subplot2grid((4, 5), (1, 2), rowspan=3)
ax5 = plt.subplot(4, 5, 4)
ax6 = plt.subplot2grid((4, 5), (1, 3), rowspan=3)
ax7 = plt.subplot(1, 5, 5)
axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

# determine palette
pal = ["#648fff", "#fe6100", "#785ef0", "#dc267f"]

# Plot mean assembly length
totalplot = sns.barplot(
    y=df.Scaffold_Mbp,
    x=df.Class,
    hue=df.Class,
    dodge=False,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[0],
    alpha=1,
)
# Plot CG composition
gcplot = sns.violinplot(
    y=df["GC_avg"] * 100,
    inner="point",
    x=df.Class,
    hue=df.Class,
    dodge=False,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    cut=2,
    palette=pal,
    ax=axs[1],
    alpha=1,
)

# Plot contig and scaffold count in two different axis for better visualization
axs[2].vlines(x=df.index, ymin=df.Scaffold_n, ymax=df.Contig_n, linewidths=2)
sns.scatterplot(
    x=df.index,
    y=df["Scaffold_n"],
    s=100,
    hue=df.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[2],
    alpha=1,
)
countplot = sns.scatterplot(
    x=df.index,
    y=df["Contig_n"],
    s=50,
    hue=df.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[2],
    alpha=1,
)
axs[3].vlines(x=df.index, ymin=df.Scaffold_n, ymax=df.Contig_n, linewidths=2)
sns.scatterplot(
    x=df.index,
    y=df["Scaffold_n"],
    s=100,
    hue=df.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[3],
    alpha=1,
)
countplot2 = sns.scatterplot(
    x=df.index,
    y=df["Contig_n"],
    s=50,
    hue=df.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[3],
    alpha=1,
)

axs[2].set_ylim(bottom=650, top=850)  # outliers only
axs[3].set_ylim(bottom=-5, top=500)
axs[2].set_xticklabels(range(700, 800, 100))

# remove axis and ticks from the redundant axis
axs[2].spines["bottom"].set_visible(False)
axs[3].spines["top"].set_visible(False)
axs[2].xaxis.tick_top()
axs[2].tick_params(labeltop="off")
axs[3].xaxis.tick_bottom()
axs[2].tick_params(labeltop="off", bottom=False, top=False)

d = 0.015  # width of diagonal lines in axes coordinates

# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=axs[2].transAxes, color="k", clip_on=False)
axs[2].plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
kwargs.update(transform=axs[3].transAxes)  # switch to the bottom axes
axs[3].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal


# plot contig and scaffold N50 (same procedure than before)

axs[4].vlines(
    x=df.index, ymin=df.Contig_N50, ymax=df.Scaffold_N50, linewidths=2, alpha=1
)
sns.scatterplot(
    x=df.index,
    y=df["Contig_N50"],
    hue=df.Class,
    palette=pal,
    s=50,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[4],
    alpha=1,
)
N50plot = sns.scatterplot(
    x=df.index,
    y=df["Scaffold_N50"],
    hue=df.Class,
    s=100,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[4],
    alpha=1,
)
axs[5].vlines(x=df.index, ymin=df.Contig_N50, ymax=df.Scaffold_N50, linewidths=2)
sns.scatterplot(
    x=df.index,
    y=df["Contig_N50"],
    hue=df.Class,
    palette=pal,
    s=50,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[5],
    alpha=1,
)

N50plot2 = sns.scatterplot(
    x=df.index,
    y=df["Scaffold_N50"],
    hue=df.Class,
    s=100,
    palette=pal,
    hue_order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    ax=axs[5],
    alpha=1,
)

axs[4].set_ylim(bottom=11800, top=13300)  # outliers only
axs[4].set_xticklabels(range(1200, 1300, 1000))
axs[5].set_ylim(bottom=-5, top=5000)

axs[4].spines["bottom"].set_visible(False)
axs[5].spines["top"].set_visible(False)
axs[4].tick_params(labeltop="off", bottom=False, top=False)

# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=axs[4].transAxes, color="k", clip_on=False)
axs[4].plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal

kwargs.update(transform=axs[5].transAxes)  # switch to the bottom axes
axs[5].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

# Plot XYZ score
bar2 = sns.barplot(
    x=bdf[bdf.Class == "Staurozoa"]["X.Y.Z score"],
    y=bdf[bdf.Class == "Staurozoa"]["Count"],
    color=pal[1],
    ax=axs[6],
    alpha=1,
)
biogenomeplot = sns.barplot(
    x=bdf[bdf.Class == "Scyphozoa"]["X.Y.Z score"],
    y=bdf[bdf.Class == "Scyphozoa"]["Count"],
    estimator=sum,
    ci=None,
    color=pal[3],
    ax=axs[6],
    alpha=1,
)
biogenomeplot = sns.barplot(
    x=bdf[bdf.Class == "Hydrozoa"]["X.Y.Z score"],
    y=bdf[bdf.Class == "Hydrozoa"]["Count"],
    estimator=sum,
    ci=None,
    color=pal[0],
    ax=axs[6],
    alpha=1,
)
biogenomeplot = sns.barplot(
    x=bdf[bdf.Class == "Cubozoa"]["X.Y.Z score"],
    y=bdf[bdf.Class == "Cubozoa"]["Count"],
    estimator=sum,
    ci=None,
    color=pal[2],
    ax=axs[6],
    alpha=1,
)
biogenomeplot = sns.barplot(
    x=bdf[bdf.Class == "Staurozoa"]["X.Y.Z score"],
    y=bdf[bdf.Class == "Staurozoa"]["Count"],
    estimator=sum,
    ci=None,
    color=pal[1],
    ax=axs[6],
    alpha=1,
)

# set plot titles
plots_titles = {
    totalplot: ["Mean assembly length", "(Mbp)"],
    gcplot: ["GC content", "Percentage (%)"],
    countplot: ["Contig & scaffold number", "Count in thousands"],
    N50plot: ["Scaffold & contig N50", "Length (Kbp)"],
    biogenomeplot: ["Assembly quality", "Count"],
    N50plot2: ["", ""],
    countplot2: ["", ""],
}

# set plot titles and axis labels
for plot, name in plots_titles.items():
    if plot != biogenomeplot:
        plot.get_legend().remove()
    yaxis = plot.get_yaxis()
    plot.set_title(name[0], size=15)
    if plot in [countplot, N50plot]:
        plot.set_xticklabels([""])
    yaxis.set_label_text(name[1], size=12)


for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticklabels("")
    ax.xaxis.label.set_fontsize(12)


plt.rcParams["svg.fonttype"] = "none"
fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
fig.autofmt_xdate(rotation=90)
# save
# fig.savefig("Assemblystats.svg")


# ## Figure 2B 

# In[ ]:


# load annotation statistics
annotationdf = pd.read_csv("Data/Annotation_stats.csv", sep="\t")

# sort by class
annotationdf["Order"] = annotationdf.Class.apply(
    lambda x: {"Hydrozoa": 1, "Staurozoa": 2, "Cubozoa": 3, "Scyphozoa": 4}[x]
)
annotationdf = annotationdf.sort_values("Order")

# define palette
pal = ["#648fff", "#785ef0", "#dc267f"]

# define subplots
fig, axs = plt.subplots(ncols=5, nrows=1, figsize=[18, 5])


annotationdf.replace("not_informed", np.nan, inplace=True)

# define a second df to handle cases when Staurozoa values are missing
genedf = annotationdf.copy()
annotationdf = annotationdf[annotationdf.Class != "Staurozoa"]

# plot repeat content
repeat_f = sns.barplot(
    y=annotationdf["Total repetitive length"].div(1000000),
    x=annotationdf.Class,
    order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    palette=pal,
    hue_order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    dodge=False,
    ci=100,
    ax=axs[0],
    alpha=1,
)
# plot number of genes
totaln_f = sns.barplot(
    y=genedf["Number of genes"],
    x=genedf.Class,
    order=["Hydrozoa", "Staurozoa", "Cubozoa", "Scyphozoa"],
    dodge=False,
    ci=100,
    palette=["#648fff", "#fe6100", "#785ef0", "#dc267f"],
    ax=axs[1],
    alpha=1,
)

# plot exon count
exon = sns.barplot(
    x=annotationdf.Class,
    y=annotationdf["Mean exons per cds"],
    order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    palette=pal,
    dodge=False,
    ci=100,
    ax=axs[2],
    alpha=1,
)

# drop cases when no info about exon intron and gene length are missing
annotationdf2 = annotationdf.dropna(
    subset=["Mean exon length", "Mean gene length", "Mean intron length"], how="all"
)
x1 = annotationdf2.Species
axs[3].vlines(
    x=x1,
    ymin=annotationdf2["Mean exon length"].div(1000),
    ymax=annotationdf2["Mean gene length"].div(1000),
    linewidths=2,
    alpha=1,
)
lengtp = sns.scatterplot(
    x=x1,
    y=annotationdf2["Mean exon length"].div(1000),
    hue=annotationdf2.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    s=50,
    ax=axs[3],
    alpha=1,
)
lengtp = sns.scatterplot(
    x=x1,
    y=annotationdf2["Mean gene length"].div(1000),
    hue=annotationdf2.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    s=120,
    ax=axs[3],
    alpha=1,
)
lengtp = sns.scatterplot(
    x=x1,
    y=annotationdf2["Mean intron length"].div(1000),
    hue=annotationdf2.Class,
    palette=pal,
    hue_order=["Hydrozoa", "Cubozoa", "Scyphozoa"],
    s=90,
    ax=axs[3],
    alpha=1,
)

# remove some axis to improve visualization
for ax in axs:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticklabels("")
    ax.xaxis.label.set_fontsize(12)

# set plot titles
plots = [totaln_f, repeat_f, exon, lengtp]
titles = {
    totaln_f: [
        "Total number of genes",
        "Count n=" + str(genedf["Number of genes"].count()),
    ],
    exon: [
        "Mean exon number",
        "Count per gene n=" + str(annotationdf["Mean exons per cds"].count()),
    ],
    repeat_f: [
        "Total repeat length",
        "(Mbp) n=" + str(annotationdf["Total repetitive length"].count()),
    ],
    lengtp: ["Gene feature length", "(Kbp)", ""],
}

for plott in plots:
    plott.set_title(titles[plott][0], size=15)
    plott.axes.text(
        0.5,
        0.98,
        titles[plott][1],
        horizontalalignment="center",
        transform=plott.axes.transAxes,
        size=12,
    )
    plott.legend().set_visible(False)


fig.savefig("Genestats.svg")


# ## Figure 3
# This code was used to plot BUSCO Metazoa distribution.
# The same code was run with results from BUSCO eukaryota from supplementary files but is not indicated here. To run with eukaryota files must change the path variable to path_eukarya

# In[ ]:


#define path where BUSCO outputs are found
buscodf = pd.DataFrame()
path_metazoa = "Data/BUSCO_metazoa/"
path_eukarya = "Data/BUSCO_eukarya/"
path = path_metazoa

# import tree
t = Tree("Data/BUSCO_TREE.tre", quoted_node_names=True)
# import taxonomy info.
taxonomydf = pd.read_csv("Data/Taxonomy.csv", sep="\t")


# read BUSCO outputs and concatenate in a single df
for element in os.listdir(path):
    if element.split(".")[-1] == "BUSCO":
        df = pd.read_csv(path + element, skiprows=2, sep="\t")
        df["Assembly"] = element.split(".")[0]
        buscodf = buscodf.append(df)

ref = pd.DataFrame()
ref["# Busco id"] = buscodf["# Busco id"]

buscodf["# Busco id"] = "BUSCO" + buscodf["# Busco id"].apply(
    lambda x: x.split("at")[0]
)


# In[ ]:


def codif(x):
    if x == "Complete":
        return 3
    elif x == "Duplicated":
        return 1
    elif x == "Fragmented":
        return 2
    elif x == "Missing":
        return 0


# codify presence and absence with numbers
df = buscodf.groupby(["# Busco id", "Assembly"]).agg({"Status": list})
df.Status = df.Status.apply(lambda x: x[0])
df.Status = df.Status.apply(lambda x: codif(x))
df = df.unstack()

# sort values by "completeness" .. those genes that are frequently present
# complete or fragmented will be plotted in the left side of the plot
df["order"] = df.apply(lambda x: sum(x.Status), axis=1)
df.sort_values(by="order", ascending=False, inplace=True)
df = df.drop(columns="order")
df.columns = df.columns.droplevel()


# In[ ]:


def append_att(x, t=t):
    node = t.get_leaves_by_name(x.name)[0]
    node.add_feature("assembly", x.name)
    for gene in x.index:
        node.add_feature(gene, x[gene])
    globals()[x.name] = node
    return globals()[x.name]


def get_sp(x, tree):
    for node in t.traverse():
        if node.is_leaf():
            if re.search(x.Assembly_accession.split(".")[0], node.assembly) != None:
                node.name = x.Species


def set_color(att):
    if att == 1:  # duplicated pink
        color = "#dc267f"
    elif att == 2:  # fragmented violet
        color = "#785ef0"
    elif att == 3:
        color = "#648fff"  # complete blue
    elif att == 0:  # missing white
        color = "white"
    return color

##add gene information to each tree node (leaf)
df.apply(lambda x: append_att(x))


# In[ ]:


##add squares corresponding to genes for each terminal (heron referred as nodes)
nodes = list(df.columns)
realnodes = []
for node in nodes:
    realnode = globals()[node]
    genes = list(df.index)
    for indx in range(len(genes)):
        att = getattr(realnode, genes[indx])
        color = set_color(att)
        rect = RectFace(10, 10, color, bgcolor=color, label=None)
        realnode.add_face(rect, column=indx, position="aligned")
    realnodes.append(realnode)

att_list = list(df.index)
taxonomydf.apply(lambda x: get_sp(x, t), axis=1)

##build headers
ts = TreeStyle()
headers_p = (1, 10, 0)
for indx in range(len(att_list)):
    header = TextFace(att_list[indx], ftype="Inter Medium", fsize=8)
    header.rotable = True
    header.rotation = 270
    ts.aligned_foot.add_face(header, column=indx)


###define style and add assemblies names
ts.show_leaf_name = False
ts.min_leaf_separation = 2
ts.show_scale = False
ts.allow_face_overlap = True
ts.draw_guiding_lines = True
ts.complete_branch_lines_when_necessary = True
ts.extra_branch_line_type = 0
ts.extra_branch_line_color = "black"
ts.force_topology = True
ts.guiding_lines_type = 0
ts.guiding_lines_color = "black"


#add assembly accessions
for n in t.traverse():
    n.img_style["size"] = 0
    n.img_style["vt_line_width"] = 3
    n.img_style["hz_line_width"] = 3
    if n.is_leaf():
        name = TextFace(n.name, ftype="inter", fsize=12, fgcolor="black")
        name.margin_left = 0
        name.margin_right = 0
        name.fstyle = "italic"
        name.hz_align = 0
        n.add_face(name, column=-1, position="aligned")
        assembly = TextFace(n.assembly, ftype="inter", fsize=8, fgcolor="black")
        assembly.margin_left, assembly.margin_right, assembly.hz_align = 0, 0, 0
        n.add_face(assembly, column=-1, position="float")


# In[ ]:


#Create references
C = RectFace(50, 50, "#648fff", bgcolor="#648fff")
Cl = TextFace("Complete", ftype="inter", fsize=15, fgcolor="black")
M = RectFace(50, 50, "black", bgcolor="white")
Ml = TextFace("Missing", ftype="inter", fsize=15, fgcolor="black")
F = RectFace(50, 50, "#785ef0", bgcolor="#785ef0")
Fl = TextFace("Fragmented", ftype="inter", fsize=15, fgcolor="black")
D = RectFace(50, 50, "#dc267f", bgcolor="#dc267f")
Dl = TextFace("Duplicated", ftype="inter", fsize=15, fgcolor="black")

# add references
ts.legend.add_face(C, column=1)
ts.legend.add_face(Cl, column=1)
ts.legend.add_face(F, column=2)
ts.legend.add_face(Fl, column=2)
ts.legend.add_face(D, column=3)
ts.legend.add_face(Dl, column=3)
ts.legend.add_face(M, column=4)
ts.legend.add_face(Ml, column=4)
ts.legend_position = 4
D.margin_right = 60
D.margin_top = 10
C.margin_right = 60
C.margin_top = 10
F.margin_right = 60
F.margin_top = 10
M.margin_right = 60
M.margin_top = 10

Dl.margin_right = 60
Cl.margin_right = 60
Fl.margin_right = 60
Ml.margin_right = 60

t.show(tree_style=ts)
t.render(
    "BUSCO_tree_eukaria.svg", layout=None, tree_style=ts, units="px", h=500, dpi=200
)

