# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 14:17:45 2024

@author: aaron

# PCR error plots
"""

# PCR error rate of -1 vs +1 allele in CMMRD and Control (per marker)

from scripts.common_imports import *

def plot1(data, name):
    
    data = data.errorratio()
    
    markerlist = list(data['Marker'].unique())
    
    data2 = data.loc[data['L1'] == "MMRd"]
    
    data3 = data.loc[data['L1'] == "MMRp"]
    
    def temp(data):
    
        d = []
        
        for marker in markerlist:
            
            fdata = data.loc[data['Marker'] == marker]
            
            fdata = fdata.dropna(subset=['error_-1', 'error_1'])
            
            list1 = list(fdata['error_-1'])
            list2 = list(fdata['error_1'])
            
            t_statistic, p_value = stats.ttest_rel(list1, list2)
            
            x = np.average(list1)
            y = np.average(list2)
            
            d.append({
                "Marker": marker,
                "p_value": p_value,
                "average_error_-1": x,
                "average_error_1": y
            })
            
        df = pd.DataFrame(d)
        
        return(df)
    
    df2, df3 = list(map(temp, [data2, data3]))
    
    df2['Subset'], df3['Subset'] = ["MMRd", "MMRp"]
    
    df = pd.concat([df2, df3]).reset_index().drop(columns=['index'])
    
    # This filters out GM09 which is acting weirdly
    df = df[(0 <= df["average_error_-1"]) & (df["average_error_-1"] <= 1) & 
        (0 <= df["average_error_1"]) & (df["average_error_1"] <= 1)]
    
    df['Legend'] = df['p_value'].apply(lambda p: 'green' if p < 0.05 else 'gray')
    
    df = pd.melt(df, 
                  id_vars=['Marker', 'p_value', 'Subset', 'Legend'],
                  value_vars=['average_error_-1', 'average_error_1'],
                  var_name='error_type',
                  value_name='average_error')
    
    error_colors = {"average_error_-1": "red", "average_error_1": "blue"}
    
    plot = (
        p9.ggplot(df, p9.aes(x='average_error', y='Marker')) +
        p9.geom_segment(p9.aes(x=0, xend=1, y='Marker', yend='Marker', group='Subset', color='Legend')) +  # Use line_color for coloring
        p9.geom_point(p9.aes(color='error_type'), size=3) +
        p9.facet_wrap('~Subset') +
        p9.scale_color_manual(values={**error_colors, 'gray': 'gray', 'green': 'green'}, labels=["NS", "p<0.05", "error_-1", "error_1"]) +  # Ensure both colors are defined
        p9.scale_x_continuous(breaks=[0.0, 0.5, 1.0]) +
        p9.labs(x="Error Value", y="Marker", title="PCR Errors by Subset " + "(" + name +")" ) +
        p9.theme_gray() +
        p9.theme(panel_spacing=0.05)
    )
    
    plot.show()
    

# PCR error rate CMMRD vs Control in -1 and +1 allele (per marker)
def plot2(data, name):
    
    data = data.errorratio()
    
    markerlist = list(data['Marker'].unique())
    
    def temp(val):
        
        d = []
        
        for marker in markerlist:
            
            fdata = data.loc[data['Marker'] == marker]
            
            CMMRD = fdata.loc[fdata["L1"] == "MMRd"]
            Control = fdata.loc[fdata["L1"] == "MMRp"]
        
            list1 = CMMRD["error_"+val].dropna().tolist()
            list2 = Control["error_"+val].dropna().tolist()
            
            t_statistic, p_value = stats.ttest_ind(list1, list2)
            
            x = np.average(list1)
            y = np.average(list2)
            
            d.append({
                
                "Marker": marker,
                "p_value": p_value,
                "average_MMRd" : x,
                "average_MMRp": y
                
                })
            
        df = pd.DataFrame(d)
        
        return(df)
        
    df1, df2 = list(map(temp, ["-1", "1"]))
    
    df1['Subset'], df2['Subset'] = ["-1", "1"]
    
    df = pd.concat([df1, df2]).reset_index().drop(columns=['index'])
    
    # This filters out GM09 which acts weirdly 
    df = df[(0 <= df["average_MMRd"]) & (df["average_MMRd"] <= 1) & 
        (0 <= df["average_MMRp"]) & (df["average_MMRp"] <= 1)]
    
    df['Legend'] = df['p_value'].apply(lambda p: 'green' if p < 0.05 else 'gray')
        
    df = pd.melt(df, 
                  id_vars=['Marker', 'p_value', 'Subset', 'Legend'],
                  value_vars=['average_MMRd', 'average_MMRp'],
                  var_name='type',
                  value_name='average_error')
    
    error_colors = {"average_MMRd": "red", "average_MMRp": "blue"}
    
    plot = (
        p9.ggplot(df, p9.aes(x='average_error', y='Marker')) +
        p9.geom_segment(p9.aes(x=0, xend=1, y='Marker', yend='Marker', group='Subset', color='Legend')) +  # Use line_color for coloring
        p9.geom_point(p9.aes(color='type'), size=3) +
        p9.facet_wrap('~Subset') +
        p9.scale_color_manual(values={**error_colors, 'gray': 'gray', 'green': 'green'}, labels=["NS", "p<0.05", "MMRd", "MMRp"]) +  # Ensure both colors are defined
        p9.scale_x_continuous(breaks=[0.0, 0.5, 1.0]) +
        p9.labs(x="Error Value", y="Marker", title="PCR Errors by Allele " + "(" + name + ")") +
        p9.theme_gray() +
        p9.theme(panel_spacing=0.05)
    )
    
    plot.show()