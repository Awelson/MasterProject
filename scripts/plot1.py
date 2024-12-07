# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 00:32:36 2024

@author: aaron
"""

from scripts.common_imports import *

def plot(data):
    
    by = data.columns[2]
    
    # def rgb_to_hex(rgb_tuple):
    #     r, g, b = [int(c * 255) for c in rgb_tuple]
    #     return "#{:02x}{:02x}{:02x}".format(r, g, b)
    
    colors = sns.color_palette("Accent", n_colors=3)
    
    colors = list(map(rgb_to_hex, colors))
    
    plot = (
        p9.ggplot(data, p9.aes(y="Ratio", x=by, fill=by))
        + p9.geom_boxplot()
        + p9.facet_wrap("UMI", ncol=3)
        + p9.scale_fill_manual(values=colors)
        + p9.theme_gray()
        )
    
    plot.show()