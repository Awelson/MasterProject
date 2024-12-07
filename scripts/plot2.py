# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:28:28 2024

@author: aaron

# Score plots
"""

from scripts.common_imports import *

def plot(data):
    
    data["-log_10p"] = -np.log10(data["pvalue"])
    
    Alleles = list(data["Allele"].unique())
    
    UMIs = list(data["UMI"].unique())
    
    by = data.columns[2]
    
    if len(UMIs) > 1:
        
        plot = (
            p9.ggplot(data, p9.aes(y="Allele", x="AUC", color="-log_10p"))
            + p9.geom_point()
            + p9.scale_y_discrete(limits = Alleles)
            + p9.scale_color_gradientn(
                colors=["blue", "orange", "red"],
                values=[0,0.1,1]
            )
            + p9.geom_vline(xintercept=0.5, color="red", linetype="dashed", size=0.4)
            + p9.facet_grid(f'{by} ~ UMI')
            + p9.labs(title="Groupings")
            + p9.theme_gray()
            + p9.theme(panel_spacing=0.05) 
        )
        
        plot.show()
        
    else:
        
        data = data.drop(columns = ["UMI"])
        
        plot = (
            p9.ggplot(data, p9.aes(y="Allele", x="AUC", color="-log_10p"))
            + p9.geom_point()
            + p9.scale_y_discrete(limits = Alleles)
            + p9.scale_color_gradientn(
                colors=["blue", "orange", "red"],
                values=[0,0.1,1]
            )
            + p9.geom_vline(xintercept=0.5, color="red", linetype="dashed", size=0.4)
            + p9.facet_wrap(by, ncol=3)
            + p9.theme_gray()
            + p9.theme(panel_spacing=0.05) 
        )
        
        plot.show()
    
    # Filter
    
    # data = (
    #     data.loc[data["UMI"] == Dx]
    #         .drop(columns=['UMI'])
    # )

    # plot = (
    #     p9.ggplot(data, p9.aes(y="Allele (length)", x="AUC", color="-log_10p"))
    #     + p9.geom_point()
    #     + p9.scale_y_discrete(limits = varlist)
    #     + p9.scale_color_gradientn(
    #         colors=["blue", "orange", "red"],
    #         values=[0,0.1,1]
    #     )
    #     + p9.geom_vline(xintercept=0.5, color="red", linetype="dashed", size=0.4)
    #     + p9.facet_wrap(by, ncol=3)
    #     + p9.theme_gray()
    #     + p9.theme(panel_spacing=0.05) 
    # )
    
    # plot.show()
    
def plot2(data, by, type1):
    
    data = data.plot2data(by, type1)
    
    if type1 == "AF":
        varlist = ['-2', '-1', '0', '1']
    elif type1 == "VF":
        varlist = ['-2', '-1', '1']

    plot = (
        p9.ggplot(data, p9.aes(y="Allele (length)", x="AUC", color="-log_10p"))
        + p9.geom_point()
        + p9.scale_y_discrete(limits = varlist)
        + p9.scale_color_gradientn(
            colors=["blue", "orange", "red"],
            values=[0,0.1,1]
        )
        + p9.geom_vline(xintercept=0.5, color="red", linetype="dashed", size=0.4)
        + p9.facet_grid(f'{by} ~ UMI')
        + p9.labs(title="by UMI")
        + p9.theme_gray()
        + p9.theme(panel_spacing=0.05) 
    )
    
    plot.show()