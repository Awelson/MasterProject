# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:38:30 2024

@author: aaron
"""

from scripts.common_imports import *

# data is a an object from the class 'Sample2'

def plot1(data, Dx, alpha = 0):
    
    data = data.commonSNPcounts()
    
    # Filter
    
    data = (
        data.loc[data["UMI"] == Dx]
            .drop(columns=['UMI', 'L2', 'L3', 'Allele 1', 'Allele 2'])
    )
    
    # Ratio Counts 2 / Counts 1
    
    data["MAF"] = data["Counts 2"] / (data["Counts 1"] + data["Counts 2"])
    
    
    # alpha filter
    
    data = data.loc[data["MAF"] > alpha]
    
    # Plot
    
    plot = (
        p9.ggplot(data, p9.aes(x='MAF', fill='L1'))
        + p9.geom_histogram(p9.aes(y=p9.after_stat('density'), binwidth=0.05))
        + p9.facet_wrap(['Dataset', 'L1'], ncol=2)
        + p9.theme_gray()
        + p9.labs(title='-', x='MAF', y='Density')
    )
    
    plot.show()
    
def plot2(data, dataset, alpha = 0):
    
    data = data.commonSNPcounts()
    
    # Filter
    
    data = (
        data.loc[data["Dataset"] == dataset]
            .drop(columns=['Dataset', 'L2', 'L3', 'Allele 1', 'Allele 2'])
    )
    
    # Ratio Counts 2 / Counts 1
    
    data["MAF"] = data["Counts 2"] / (data["Counts 1"] + data["Counts 2"])
    
    # alpha filter
    
    data = data.loc[data["MAF"] > alpha]
    
    plot = (
        p9.ggplot(data, p9.aes(x='MAF', fill='L1'))
        + p9.geom_histogram(p9.aes(y=p9.after_stat('density')), alpha=0.5, position="identity", binwidth=0.05)
        + p9.facet_wrap('UMI', nrow=3)
        + p9.theme_gray()
        + p9.labs(title='-', x='MAF', y='Density')
    )
        
    plot.show()
    
def thresholds(data, Dx, dataset, alpha = 0.05, ran = (1000, 50)):
    
    data = data.commonSNPcounts()
    
    # Filter
    
    data = (
        data.loc[data["UMI"] == Dx]
            .loc[data["Dataset"] == dataset]
            .drop(columns=['UMI', 'Dataset', 'L2', 'L3', 'Allele 1', 'Allele 2'])
    )
    
    # Ratio Counts 2 / Counts 1
    
    data["MAF"] = data["Counts 2"] / (data["Counts 1"] + data["Counts 2"])
    
    filtered1 = data.loc[data["MAF"] > alpha]
    
    tlist = list(range(0, ran[0] - 1, ran[1]))
    
    list1 = []
    
    for beta in tlist:
    
        filtered2 = filtered1.loc[filtered1["Counts 1"] + filtered1["Counts 2"] > beta]
        
        markerlist = list(filtered2["Marker"].drop_duplicates())
        
        d = []
        
        for marker in markerlist:
            
            filtered3 = filtered2.loc[filtered2["Marker"] == marker]
            
            MMRp = filtered3.loc[filtered3["L1"] == "MMRp"]
            
            MMRd = filtered3.loc[filtered3["L1"] == "MMRd"]
            MAF1 = list(MMRp["MAF"])
            MAF2 = list(MMRd["MAF"])
            
            # Wilcoxon p-val
    
            pval = stats.mannwhitneyu(MAF1, MAF2, alternative='less')[1]
            
            # AUC
            
            x = np.array(["0", "1"])
            y_true = np.repeat(x, [len(MAF1), len(MAF2)], axis = 0)
            
            y_scores = MAF1 + MAF2
            
            auc = roc_auc_score(y_true, y_scores)
            
            # dictionary
            
            d1 = {
                
                "Marker": marker,
                "pval" : pval,
                "AUC" : auc
                
                }
            
            d.append(d1)
            
        df = pd.DataFrame(d)
        
        fdf = df.loc[df["pval"] < 0.05]
        
        x = len(fdf)
        
        print(fdf)
        
        list1.append((beta, x, fdf["AUC"].sum()/x))
        
    return(list1)
        
# def plot2(data, Dx, dataset, alpha = 0.1):
    
#     data = plot2data(data, Dx, dataset, alpha)
    
#     data["logpval"] = -np.log10(data["pval"])
    
#     plot = (
#         p9.ggplot(data, p9.aes(x='AUC', y='logpval'))
#         + p9.geom_point()
#         + p9.geom_vline(xintercept=0.5, color='red', linetype='dotted')
#         + p9.geom_hline(yintercept=1.3, color='red', linetype='dotted')
#         + p9.theme_gray()
#         )
    
#     plot.show()
    
def IM16plot(data, Dx, dataset, alpha = 0.05):
    
    data = data.commonSNPcounts()
    
    # Filter
    
    data = (
        data.loc[data["UMI"] == Dx]
            .loc[data["Dataset"] == dataset]
            .loc[data["Marker"] == "IM16_SNP1"]
            .drop(columns=['UMI', 'Dataset', 'L2', 'L3', 'Allele 1', 'Allele 2'])
    )
    
    # Ratio Counts 2 / Counts 1
    
    data["MAF"] = data["Counts 2"] / (data["Counts 1"] + data["Counts 2"])
    
    # alpha filter
    
    data = data.loc[data["MAF"] > alpha]
    
    # Plot
    
    plot = (
        p9.ggplot(data, p9.aes(x='MAF', fill='L1'))
        + p9.geom_histogram(p9.aes(y='..count..'), alpha=0.5, position="identity", binwidth=0.05)
        + p9.theme_gray()
        + p9.labs(title='-', x='MAF', y='Counts')
    )
    
    plot.show()
    
    