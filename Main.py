# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:02:55 2024

@author: aaron
"""

from scripts.common_imports import *
pl.enable_string_cache()
import scripts.Class as C
import scripts.plot1 as plot1
import scripts.plot2 as plot2

# Loading data

Alldata = pl.read_csv("Data/All.csv", separator = ",", infer_schema_length = 10000, schema_overrides={"SampleName": pl.String, "Allele (length)": pl.String, "Marker": pl.String, "UMI": pl.String, "L1": pl.String, "Dataset": pl.String, "Type": pl.String, "Allele (SNP)": pl.String})

All = C.Class(Alldata)

# QC step

toremove = All.BadSamples(1000, 100, 10).pandas

toremove = list(toremove["SampleName"].unique())

Allf = All.Filter(toremove)

# plot2

Alleles = ["Deletion", "Insertion"]

UMIs = ["D0", "D1", "D2"]

test = Allf.byMut().noZero().Freq().Pivot().byType().plot2data(Alleles, UMIs).pandas

test = Allf.byLen().Freq().Pivot().byType().plot1data().pandas

plot2.plot(test)

plot1.plot(test)

