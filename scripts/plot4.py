# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 21:35:32 2024

@author: aaron
"""

from scripts.common_imports import *

def plot1(data, Dataset, allele):
    
    data = data.plot2data()
    
    # Filter
    
    data = (
        data[data["Dataset"] == Dataset]
    )
    
    data = data[["Marker", "L1", "UMI", allele]]
    
    data["Marker"] = data["Marker"].astype(str)
    
    data["L1"] = data["L1"].astype(str)
    
    data["ID"] = data["Marker"] + "." + data["L1"]
    
    data["med"] = list(map(np.median, data[allele]))
    
    data["UMI"] = pd.Categorical(data["UMI"], categories=["D2", "D1", "D0"], ordered=True)
    
    # Statistics
    
    df = data.pivot(index=["ID", "Marker", "L1"], columns="UMI", values="med").reset_index()
    
    # MMRp CI
    
    MMRp = df[df["L1"] == "MMRp"]
    
    MMRp["D0-D2"] = MMRp["D0"] - MMRp["D2"]
    
    # MMRd
    
    MMRd = df[df["L1"] == "MMRd"]
    
    MMRd["D0-D2"] = MMRd["D0"] - MMRd["D2"]
    
    statistic, pval = stats.wilcoxon(MMRp["D0-D2"], MMRd["D0-D2"], alternative="greater")
    
    print(pval)
    
    plot = (
        p9.ggplot(data, p9.aes(x="UMI", y="med", group="ID", color="L1"))
        + p9.geom_line(size=0.8, alpha=0.5)  # Connect lines for each sample
        + p9.scale_color_manual(values=["red", "blue"])  # Red for Group 1, Blue for Group 2
        + p9.labs(
            title=Dataset + " (" + allele + ")",
            x="UMI",
            y="Relative Frequency",
            color="Groups"
        )
        + p9.theme_gray()
    )
    
    plot.show()