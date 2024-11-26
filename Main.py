# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:02:55 2024

@author: aaron
"""

from scripts.common_imports import *
pl.enable_string_cache()
import scripts.Class as C
import scripts.plot4 as plot4
import scripts.plot3 as plot3
import scripts.plot2 as plot2
import scripts.plot1 as plot1


# Loading withSNP data

Alldata = pl.read_csv("Data/All.csv", separator = ",", infer_schema_length = 10000, schema_overrides={"SampleName": pl.String, "Allele (length)": pl.String, "Marker": pl.String, "UMI": pl.String, "L1": pl.String, "L2": pl.String, "L3": pl.String, "Dataset": pl.String, "Type": pl.String, "Allele (SNP)": pl.String})

All = C.Sample(Alldata)

# Creating a filtered version of "All"

Allfdata = All.fil(1000, 100, 100)

# Loading Allfdata as an object

Allfdata = pl.from_pandas(Allfdata)

Allf = C.Sample(Allfdata)

plot2.plot2(Allf, "Type", "AF")

