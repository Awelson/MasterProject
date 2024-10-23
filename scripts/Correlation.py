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
        
# Marker Lengths

Markerdata = pl.read_csv("RawData/Positions.20190515.csv")
Markerdata = (Markerdata.select(["Gene", "Ms.Length"])).rename({"Gene": "Marker"})

lengths = {row['Marker']: row['Ms.Length'] for row in Markerdata.to_dicts()}


# Define the Marker class

class Marker:
    def __init__(self, MarkerName, allele_frequencies):
        """
        marker_id: ID or name of the marker (e.g., 1, 2, 3)
        allele_frequencies: A dictionary with allele as key and frequency as value
                            e.g., {'-1': 0.3, '0': 0.5, '1': 0.2}
        """
        self.MarkerName = MarkerName
        self.Allele = allele_frequencies

    def __repr__(self):
        return f"Marker(MarkerName={self.MarkerName}, allele_frequencies={self.Allele})"

# Define the Sample class

class Sample:
    def __init__(self, SampleName, MMR, Type, Gene):
        """
        sample_id: ID or name of the sample
        markers: Dictionary of marker objects {marker_id: Marker object}
        """
        self.SampleName = SampleName
        self.MMR = MMR
        self.Type = Type
        self.Gene = Gene
        self.Marker = {}

    def __repr__(self):
        return f"Sample(SampleName={self.SampleName}, MMR={self.MMR}, Type={self.Type}, Gene={self.Gene}, markers={self.Marker})"

# Function to turn sample data (from a dataframe) into objects (stored in a dictionary)

def to_obj(data):
    samples = {}

    for row in data.iter_rows(named=True):
        SampleName = row['SampleName']
        MarkerName = row['Marker']
        Allele = row['Allele']
        RelativeFrequency = row['RelativeFrequency']
        MMR = row['MMR']
        Type = row['Type']
        Gene = row['Gene']

        # Create or get the existing sample
        if SampleName not in samples:
            samples[SampleName] = Sample(SampleName, MMR, Type, Gene)

        sample = samples[SampleName]

        # Create or get the existing marker for this sample
        if MarkerName not in sample.Marker:
            sample.Marker[MarkerName] = Marker(MarkerName, {})

        # Add allele frequency to the marker
        sample.Marker[MarkerName].Allele[Allele] = RelativeFrequency

    return samples


# Loading CRCdata

CRCdata = pl.read_csv("RawData/CRCdata2.csv", separator = ",", infer_schema_length = 10000, schema_overrides={"Allele": pl.Categorical})

CRCfullmarkerlist = CRCdata["Marker"].unique().to_list()

CRCfullsamplelist = CRCdata["SampleName"].unique().to_list()

CRCsample = to_obj(CRCdata)

# Extractor function. Given a triple (Sample, Marker, Allele), return its relative frequency (zero if not found)

def extractor(sample, marker, allele):
    try: 
        return CRCsample[sample].Marker[marker].Allele[allele]
    except Exception:
        return 0
    
extractor("1003", "GM09", "-1")
    

# Spearman Correlation

def spearman(type1, var1, var2):
    
    # Loading data
    
    data = CRCsample
    
    markerlist = CRCfullmarkerlist
    
    # Creation of filtered samplelist
    
    samplelist = CRCfullsamplelist
    samplelist = [sample for sample in samplelist if (type1 == data[sample].MMR or type1 == data[sample].Type or type1 == data[sample].Gene)]
    
    # spearman correlation (var1 vs var2) is calculated on a per marker basis
    
    temp = {}
    
    for marker in markerlist:
        allele1 = [extractor(sample, marker, var1) for sample in samplelist]
        allele2 = [extractor(sample, marker, var2) for sample in samplelist]
        spcor = stats.spearmanr(allele1, allele2).correlation
        temp[marker] = spcor
    
    return temp

list(spearman("MMRd", "0", "1").values())

def plot1(type1, type2, var1, var2):
    
    x = list(spearman(type1, var1, var2).values())
    y = list(spearman(type2, var1, var2).values())
    
    df = pd.DataFrame({
        'Type': [type1] * len(x) + [type2] * len(y),
        'Spearman_Correlation': x + y
    })
    
    plot = (
        p9.ggplot(df, p9.aes(x='Type', y='Spearman_Correlation', fill='Type')) +
        p9.geom_boxplot(width=0.4) +
        p9.geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, stackratio=1.3, alpha=0.6, binwidth=0.05) +
        p9.scale_fill_brewer(type="qualitative", palette=1) +
        p9.geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        p9.labs(y='Spearman Correlation') +
        p9.ggtitle(var1 + ' vs ' + var2 + ' (Allele Frequency)') +
        p9.theme_gray()
    )

    plot.show()
    
def plot2(type1, var1):
    
    # Loading data
    data = CRCsample
    markerlist = CRCfullmarkerlist
    
    # Creation of filtered samplelist
    samplelist = CRCfullsamplelist
    samplelist = [sample for sample in samplelist if (type1 == data[sample].MMR or type1 == data[sample].Type or type1 == data[sample].Gene)]
    
    df = []
    
    # marker loop to get df with columns (marker, avg relative frequency, length)
    for marker in markerlist:
        d = {
            'Marker': marker,
            'Frequency': stats.tmean([extractor(sample, marker, var1) for sample in samplelist]),
            'Length': lengths[marker]
        }
        df.append(d)
    
    df = pd.DataFrame(df)
    
    # spearman coefficient
    
    spcor = stats.spearmanr(df['Frequency'], df['Length']).correlation
    
    print(spcor)
    
    plot = (
        p9.ggplot(df, p9.aes(x='Length', y='Frequency')) +
        p9.geom_point(size=2) +
        p9.stat_smooth(method="rlm", se=False, linetype = "dotted", color="blue") +
        p9.labs(y= 'Relative Allele Frequency: ' + var1, x= "Marker Length") +
        p9.ggtitle(type1) +
        p9.theme_gray()
    )
    
    plot.show()


# def tablegen1(data, type1, markerlist, samplelist, var):
#     # Filter by type
#     if type1 == "all":
#         filter1 = data
#     elif type1 == "MMRp":
#         filter1 = data.filter(pl.col("MMR") == "MMRp")
#     elif type1 == "MMRd":
#         filter1 = data.filter(pl.col("MMR") == "MMRd")
#     elif type1 == "lynch":
#         filter1 = data.filter(pl.col("Type") == "lynch")
#     elif type1 == "sporadic":
#         filter1 = data.filter(pl.col("Type") == "sporadic")
#     else:
#         filter1 = data.filter((pl.col("Gene") == type1))
        
#     # Filter by markerlist
#     filter2 = filter1.filter(pl.col("Marker").is_in(markerlist))
    
#     # Filter by samplelist
#     filter3 = filter2.filter(pl.col("SampleName").is_in(samplelist))
    
#     # Filter by allele
#     filter4 = filter3.filter(pl.col("Allele").cast(str).is_in(var))
    
#     pivot = filter4.pivot(
#         values="RelativeFrequency", 
#         index=["SampleName", "Marker", "MMR", "Type", "Gene"], 
#         columns="Allele"
#     ).drop_nulls()
    
#     n = len(var)
    
#     last_n_columns = pivot.columns[-n:]
    
#     filter_condition = pl.lit(True)
#     for col in last_n_columns:
#         filter_condition &= (pl.col(col) != 0)

#     return pivot.filter(filter_condition)

# def spearman(table):
#     result = (
#         table.group_by("Marker")
#         .agg(
#             pl.map_groups(
#                 exprs = table.columns[-2:],
#                 function = lambda x: stats.spearmanr(x[0], x[1]).correlation
#             ).alias("Spearman_Correlation")
#         )
#     )
#     return result


# def plot1(data, type1, type2, markerlist, samplelist, var1, var2):
    
#     table1 = spearman(tablegen1(data, type1, markerlist, samplelist, var1, var2))
#     table2 = spearman(tablegen1(data, type2, markerlist, samplelist, var1, var2))
    
#     data1 = table1["Spearman_Correlation"].to_numpy()
#     data2 = table2["Spearman_Correlation"].to_numpy()
    
#     df = pd.DataFrame({
#         'Type': [type1] * len(data1) + [type2] * len(data2),
#         'Spearman_Correlation': list(data1) + list(data2)
#     })

#     plot = (
#         p9.ggplot(df, p9.aes(x='Type', y='Spearman_Correlation', fill='Type')) +
#         p9.geom_boxplot(width=0.4) +
#         # p9.geom_point(color='black', fill='black', alpha=0.5) +
#         p9.geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, stackratio=1.3, alpha=0.6, binwidth=0.05) +
#         p9.scale_fill_brewer(type="qualitative", palette=1) +
#         p9.geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
#         p9.labs(y='Spearman Correlation') +
#         p9.ggtitle(var1 + ' vs ' + var2 + ' (Allele Frequency)') +
#         p9.theme_gray()
#     )

#     plot.show()
    
# def plot2(data, type1, markerlist, samplelist, markerlength, var):
#     table = tablegen1(data, type1, markerlist, samplelist, var)
#     table = table.join(markerlength, on="Marker", how="left")
    
#     return table
    
# def plot3(data, type1, type2, markerlist, samplelist, var1, var2):
    
#     table1 = spearman(tablegen1(data, type1, markerlist, samplelist, var1, var2))
#     table2 = spearman(tablegen1(data, type2, markerlist, samplelist, var1, var2))
    
#     join = table1.join(table2, on="Marker", how="inner")
    
#     join = join.rename({join.columns[1]: type1, join.columns[2]: type2})
    
#     # remove outlier
#     # join = join.filter(pl.col(type2) > -0.3)
    
#     spcor = stats.spearmanr(join[type1].to_numpy(), join[type2].to_numpy()).correlation
    
#     plot = (
#         p9.ggplot(join, p9.aes(x=type1, y=type2)) +
#         p9.geom_point() +
#         p9.scale_x_continuous(breaks=np.arange(-0.4, 0.6, 0.2), limits = (-0.4, 0.6)) +
#         # scale for 0 vs -1
#         # p9.scale_x_continuous(breaks=np.arange(-0.9, 0.4, 0.2), limits = (-0.9, 0.4)) +
#         p9.geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
#         p9.geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
#         p9.labs(x = type1, y = type2) +
#         p9.stat_smooth(method="rlm", se=False, linetype = "dotted", color="blue") +
#         p9.annotate("text", x=0.4, y=0.1, label=f"Spearman Correlation: {spcor:.2f}", 
#               size=10, color="red") +
#         p9.ggtitle('Spearman Correlations ' + '(' + var1 + ' vs ' + var2 + ')') +
#         p9.theme_gray()
#     )
    
#     plot.show()