# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:02:55 2024

@author: aaron
"""

import polars as pl


CRCdata = pl.read_csv("RawData/CRCdata1.csv", separator = ",", infer_schema_length = 10000, schema_overrides={"Allele": pl.Categorical})

ECdata = pl.read_csv("RawData/ECdata1.csv", separator = ",", infer_schema_length = 10000, schema_overrides={"Allele": pl.Categorical})

CRCfullmarkerlist = CRCdata["Marker"].unique().to_list()

CRCfullsamplelist = CRCdata["SampleName"].unique().to_list()

import scripts.Correlation as C

C.plot1(CRCdata, "MMRp", "MMRd", CRCfullmarkerlist, CRCfullsamplelist, "0", "1")

C.plot3(CRCdata, "MMRp", "MMRd", CRCfullmarkerlist, CRCfullsamplelist, "0", "1")
