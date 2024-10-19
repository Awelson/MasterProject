# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:34:36 2024

@author: aaron
"""

import polars as pl
import pandas as pd
import numpy as np
from scipy import stats
import plotnine as p9
import statsmodels.api as sm

def tablegen1(data, type1, markerlist, samplelist, var1, var2):
    # Filter by type
    if type1 == "all":
        filter1 = data
    elif type1 == "MMRp":
        filter1 = data.filter(pl.col("MMR") == "MMRp")
    elif type1 == "MMRd":
        filter1 = data.filter(pl.col("MMR") == "MMRd")
    elif type1 == "lynch":
        filter1 = data.filter(pl.col("Type") == "lynch")
    elif type1 == "sporadic":
        filter1 = data.filter(pl.col("Type") == "sporadic")
    else:
        filter1 = data.filter((pl.col("Gene") == type1))
        
    # Filter by markerlist
    filter2 = filter1.filter(pl.col("Marker").is_in(markerlist))
    
    # Filter by samplelist
    filter3 = filter2.filter(pl.col("SampleName").is_in(samplelist))
    
    # Filter by allele
    filter4 = filter3.filter(pl.col("Allele").is_in([var1, var2]))
    
    pivot = filter4.pivot(
        values="RelativeFrequency", 
        index=["SampleName", "Marker", "MMR", "Type", "Gene"], 
        columns="Allele"
    ).drop_nulls()
    
    pivot = pivot.filter(
    (pl.col(pivot.columns[-1]) != 0) & (pl.col(pivot.columns[-2]) != 0))
    
    return pivot

def spearman(table):
    result = (
        table.group_by("Marker")
        .agg(
            pl.map_groups(
                exprs = table.columns[-2:],
                function = lambda x: stats.spearmanr(x[0], x[1]).correlation
            ).alias("Spearman_Correlation")
        )
    )
    return result


def plot1(data, type1, type2, markerlist, samplelist, var1, var2):
    
    table1 = spearman(tablegen1(data, type1, markerlist, samplelist, var1, var2))
    table2 = spearman(tablegen1(data, type2, markerlist, samplelist, var1, var2))
    
    data1 = table1["Spearman_Correlation"].to_numpy()
    data2 = table2["Spearman_Correlation"].to_numpy()
    
    df = pd.DataFrame({
        'Type': [type1] * len(data1) + [type2] * len(data2),
        'Spearman_Correlation': list(data1) + list(data2)
    })

    plot = (
        p9.ggplot(df, p9.aes(x='Type', y='Spearman_Correlation', fill='Type')) +
        p9.geom_boxplot(width=0.4) +
        # p9.geom_point(color='black', fill='black', alpha=0.5) +
        p9.geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, stackratio=1.3, alpha=0.6, binwidth=0.05) +
        p9.scale_fill_brewer(type="qualitative", palette=1) +
        p9.geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        p9.labs(y='Spearman Correlation') +
        p9.ggtitle(var1 + ' vs ' + var2 + ' (Allele Frequency)') +
        p9.theme_gray()
    )

    plot.show()
    
def plot3(data, type1, type2, markerlist, samplelist, var1, var2):
    
    table1 = spearman(tablegen1(data, type1, markerlist, samplelist, var1, var2))
    table2 = spearman(tablegen1(data, type2, markerlist, samplelist, var1, var2))
    
    join = table1.join(table2, on="Marker", how="inner")
    
    join = join.rename({join.columns[1]: type1, join.columns[2]: type2})
    
    # remove outlier
    # join = join.filter(pl.col(type2) > -0.3)
    
    spcor = stats.spearmanr(join[type1].to_numpy(), join[type2].to_numpy()).correlation
    
    plot = (
        p9.ggplot(join, p9.aes(x=type1, y=type2)) +
        p9.geom_point() +
        p9.scale_x_continuous(breaks=np.arange(-0.4, 0.6, 0.2), limits = (-0.4, 0.6)) +
        # scale for 0 vs -1
        # p9.scale_x_continuous(breaks=np.arange(-0.9, 0.4, 0.2), limits = (-0.9, 0.4)) +
        p9.geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        p9.geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
        p9.labs(x = type1, y = type2) +
        p9.stat_smooth(method="rlm", se=False, linetype = "dotted", color="blue") +
        p9.annotate("text", x=0.4, y=0.1, label=f"Spearman Correlation: {spcor:.2f}", 
              size=10, color="red") +
        p9.ggtitle('Spearman Correlations ' + '(' + var1 + ' vs ' + var2 + ')') +
        p9.theme_gray()
    )
    
    plot.show()