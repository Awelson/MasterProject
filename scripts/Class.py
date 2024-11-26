# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:09:06 2024

@author: aaron
"""

from scripts.common_imports import *

class Sample:
    def __init__(self, data):
        
        self.pandas = data.to_pandas()
        self.polars= data
        
    def __repr__(self):
        return f"data:{self.data}"
    
    def withoutSNP1(self):
        
        df = self.polars
        
        group_cols = ['SampleName', 'Allele (length)', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type']
        
        result = (df.group_by(group_cols)
              .agg([
                  pl.col('Counts').sum()
              ]))
    
        return(result.to_pandas())
    
    def withoutSNP2(self):
        
        df = self.withoutSNP1()
        
        df = pl.from_pandas(df)
        
        group_cols = ['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type']
        
        result = (df.group_by(group_cols)
              .agg([
                  pl.col('Counts'),
                  pl.col('Allele (length)'),
                  pl.col('Counts').sum().alias("Total")
              ]))
        
        result = result.explode(['Counts', 'Allele (length)'])        
        
        return(result.to_pandas())
    
    def withoutSNP3(self):
        
        df = self.withoutSNP1()
        
        df = df.loc[df['Allele (length)'] != "0"]
        
        df = pl.from_pandas(df)
        
        group_cols = ['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type']
        
        result = (df.group_by(group_cols)
              .agg([
                  pl.col('Counts'),
                  pl.col('Allele (length)'),
                  pl.col('Counts').sum().alias("Total")
              ]))
        
        result = result.explode(['Counts', 'Allele (length)'])        
        
        return(result.to_pandas())
    
    def totalmarkercounts1(self):
        
        df = self.withoutSNP2()
        
        df = pl.from_pandas(df)
        
        df = df.select(["SampleName", "Marker", "UMI", "Total"]).unique()
        
        return(df.to_pandas())
    
    def totalmarkercounts2(self, D0, D1, D2):
        
        df = self.totalmarkercounts1()
        
        df1 = df.loc[df["UMI"] == "D0"].loc[df["Total"] <= D0]
        df2 = df.loc[df["UMI"] == "D1"].loc[df["Total"] <= D1]
        df3 = df.loc[df["UMI"] == "D2"].loc[df["Total"] <= D2]
        
        df = pd.concat([df1,df2,df3])
        
        return(df)
    
    def fil(self, D0, D1, D2):
        
        toremove = self.totalmarkercounts2(D0, D1, D2)
        
        toremove = pl.from_pandas(toremove)
        
        df = self.polars
        
        result = df.join(toremove, on=["SampleName", "Marker", "UMI"], how="anti")
        
        return(result.to_pandas())
    
    def withoutSNP_meltcounts(self):
        
        df = self.withoutSNP2()
        
        df = pl.from_pandas(df)
        
        melted = df.pivot(
            index=['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type'],  # Columns to index
            columns='Allele (length)',  # Column to pivot
            values='Counts'
        )
        
        melted = melted.fill_null(0)
        
        # melted = df.pivot_table(index=['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset'], columns='Allele (length)', values='Counts', fill_value=0)
        
        return(melted.to_pandas())
    
    def withoutSNP_meltfrequencies(self, type1):
        
        if type1 == "AF":
            
            df = self.withoutSNP2()
        
            df["RelativeFrequency"] = df["Counts"] / df["Total"]
            
            df = pl.from_pandas(df)
            
            melted = df.pivot(
                index=['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type'],  # Columns to index
                columns='Allele (length)',  # Column to pivot
                values='RelativeFrequency'
            )
            
            melted = melted.fill_null(0)
            
            return(melted.to_pandas())
        
        elif type1 == "VF":
            
            df = self.withoutSNP3()
            
            df["VariantFrequency"] = df["Counts"] / df["Total"]
    
            df = pl.from_pandas(df)
            
            melted = df.pivot(
                index=['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type'],  # Columns to index
                columns='Allele (length)',  # Column to pivot
                values='VariantFrequency'
            )
            
            melted = melted.fill_null(0)
            
            return(melted.to_pandas())
    
    def plot2data(self, by, type1):
        
        df = self.withoutSNP_meltfrequencies(type1)
        
        df = pl.from_pandas(df)
        
        if by == "Dataset":
            group_cols = ['Marker', 'UMI', 'L1', 'Dataset']
        elif by == "Type":
            group_cols = ['Marker', 'UMI', 'L1', 'Type']
            
        if type1 == "AF":
        
            result = (df.group_by(group_cols)
                  .agg([
                      pl.col('-2'),
                      pl.col('-1'),
                      pl.col('0'),
                      pl.col('1')
                  ]))
            
        if type1 == "VF":
        
            result = (df.group_by(group_cols)
                  .agg([
                      pl.col('-2'),
                      pl.col('-1'),
                      pl.col('1')
                  ]))
            
        data = result.to_pandas()
        
        if type1 == "AF":
            varlist = ['-2', '-1', '0', '1']
        elif type1 == "VF":
            varlist = ['-2', '-1', '1']
        
        # Melt
        
        data = pd.melt(data, id_vars=['Marker', 'L1', 'UMI', by], 
                       value_vars= varlist,
                       var_name='Allele (length)', 
                       value_name='List')
        
        # Pivot
        
        data = data.pivot(index=['Marker', 'Allele (length)', 'UMI', by], 
                                  columns='L1', 
                                  values='List').reset_index()

        # pvalue with mannwhitney

        pvalue = [stats.mannwhitneyu(data['MMRd'][i], data['MMRp'][i], alternative='two-sided')[1] for i in range(len(data))]
        
        data["pvalue"] = pvalue
        
        data["-log_10p"] = -np.log10(data["pvalue"])
        
        # AUC
        
        auc = []
        
        for i in range(len(data)):
            
            row = data.loc[i]
            
            MMRd = row["MMRd"]
            MMRp = row["MMRp"]
        
            x = np.array(["1", "0"])
            y_true = np.repeat(x, [len(MMRd), len(MMRp)], axis = 0)
            y_scores = list(MMRd) + list(MMRp)
            
            aucscore = roc_auc_score(y_true, y_scores)
            
            auc.append(aucscore)
            
        data["AUC"] = auc
        
        return(data)
    
    # JUNK
    
    def stutterratio(self):
        
        df = self.data

        pivot_df = df.pivot_table(index=['SampleName', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset'], columns='Allele (length)', values='Counts', fill_value=0)

        # Calculate ratios
        pivot_df['-2'] = pivot_df["-2"] / pivot_df["0"]
        pivot_df['-1'] = pivot_df["-1"] / pivot_df["0"]
        pivot_df['1'] = pivot_df["1"] / pivot_df["0"]
                
        # Reset index to make it a flat DataFrame
        ratios_df = pivot_df[['-2', '-1', '1']].reset_index()
        
        return(ratios_df)
    
    def errorratio(self):
        
        df = self.stutterratio()
        
        pivot_df = df.pivot_table(index=['SampleName', 'Marker', 'L1', 'L2', 'L3'], columns='UMI', values=['-2','-1','1'], fill_value=0)
        
        pivot_df.columns = ['_'.join(map(str, col)).strip() for col in pivot_df.columns.values]
        
        pivot_df = pivot_df.reset_index()
        
        def temp(D0, D2):
            x = (D0-D2)/D0
            return(x)
        
        pivot_df['error_-2'] = temp(pivot_df['-2_D0'], pivot_df['-2_D2']) 
        
        pivot_df['error_-1'] = temp(pivot_df['-1_D0'], pivot_df['-1_D2']) 
        
        pivot_df['error_1'] = temp(pivot_df['1_D0'], pivot_df['1_D2']) 
        
        pivot_df = pivot_df.drop(columns = ['-2_D0', '-1_D0', '1_D0', '-2_D1', '-1_D1', '1_D1', '-2_D2', '-1_D2', '1_D2'])
        
        return(pivot_df)
    
class Sample2:
    
    def __init__(self, data):
        
        self.data = data.to_pandas()
        
    def __repr__(self):
        return f"data:{self.data}"
    
    def meltSNPcounts(self):
        return (self.data.pivot_table(
            index=['SampleName', 'Marker', 'Allele (length)', 'UMI', 'L1', 'L2', 'L3', 'Dataset'],
            columns='Allele (SNP)',
            values='Counts',
            fill_value=0
        ).reset_index())

    def commonSNPcounts(self):
        df = self.meltSNPcounts()
        metadata = df[['SampleName', 'L1', 'L2', 'L3', 'Dataset']].drop_duplicates()
        
        alleles = ['C', 'T', '-', 'A', 'G']
        
        # Total counts of SNP allele by summing along Allele (length)
        grouped = df.groupby(["SampleName", "Marker", "UMI"])[alleles].sum()
        
        values = grouped.values
        
        # Get indices of top 2 counts for each row
        top_2_indices = np.argsort(values, axis=1)[:, -2:][:, ::-1]
        
        # Get the actual values and corresponding alleles
        counts = np.take_along_axis(values, top_2_indices, axis=1)
        allele_indices = top_2_indices
        
        # Create result DataFrame
        result = pd.DataFrame({
            'SampleName': grouped.index.get_level_values(0),
            'Marker': grouped.index.get_level_values(1),
            'UMI': grouped.index.get_level_values(2),
            'Counts 1': counts[:, 0],
            'Allele 1': [alleles[i] for i in allele_indices[:, 0]],
            'Counts 2': counts[:, 1],
            'Allele 2': [alleles[i] for i in allele_indices[:, 1]]
        })
        
        # Merge with metadata
        return pd.merge(result, metadata, on='SampleName', how='left')
    
    