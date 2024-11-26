# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:34:36 2024

@author: aaron
"""

from scripts.common_imports import *

# Marker Lengths

Markerdata = pl.read_csv("RawData/MarkerLengths.csv")
Markerdata = (Markerdata.select(["Gene", "Ms.Length"])).rename({"Gene": "Marker"})

lengths = {row['Marker']: row['Ms.Length'] for row in Markerdata.to_dicts()}

# Define the Marker class

class Marker:
    def __init__(self):

        self.Allelef = {}
        self.Allelec = {}

    def __repr__(self):
        return f"allelef={self.Allelef}, allelec={self.Allelec})"

# Define the Sample class

class Sample:
    def __init__(self, SampleName, L1, L2, L3):
        
        self.SampleName = SampleName
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.Marker = {}

    def __repr__(self):
        return f"Sample(SampleName={self.SampleName}, L1={self.L1}, L2={self.L2}, L3={self.L3}"

# Function to turn sample data (from a dataframe) into objects (stored in a dictionary)

def to_obj(data):
    samples = {}

    for row in data.iter_rows(named=True):
        SampleName = row['SampleName']
        MarkerName = row['Marker']
        Allele = row['Allele']
        RelativeFrequency = row['RelativeFrequency']
        L1 = row['L1']
        L2 = row['L2']
        L3 = row['L3']
        Count = row['Count']

        # Create or get the existing sample
        if SampleName not in samples:
            samples[SampleName] = Sample(SampleName, L1, L2, L3)

        sample = samples[SampleName]

        # Create or get the existing marker for this sample
        if MarkerName not in sample.Marker:
            sample.Marker[MarkerName] = Marker()

        # Add allele frequency to the marker
        sample.Marker[MarkerName].Allelef[Allele] = RelativeFrequency
        
        # Add allele count to the marker
        sample.Marker[MarkerName].Allelec[Allele] = Count

    return samples

# set_data

def set_data(input_name):
    
    # Check if the input_name exists as a global variable
    if input_name in globals():
        temp = globals()[input_name]
        data = to_obj(temp)
        fullsamplelist = temp["SampleName"].unique().to_list()
        fullmarkerlist = temp["Marker"].unique().to_list()
    
        return(data, fullsamplelist, fullmarkerlist)
    else:
        print(f"Warning: '{input_name}' is not defined.")
        return None
    
def plot2(dataset, type1, var1):
    
    # Loading data
    
    data = set_data(dataset)
    
    data, fullsamplelist, fullmarkerlist = data
    
    def fextractor(sample, marker, var):
        try:
            return data[sample].Marker[marker].Allelef[var]
        except Exception:
            return None
    
    # Creation of filtered samplelist
    if type1 == "all":
        filteredsamplelist = fullsamplelist
    else:
        filteredsamplelist = [sample for sample in fullsamplelist if (type1 == data[sample].L1 or type1 == data[sample].L2 or type1 == data[sample].L3)]
    
    df = []
    
    # marker loop to get df with columns (marker, avg relative frequency, length)
    for marker in fullmarkerlist:
        
        frequencies = [fextractor(sample, marker, var1) for sample in filteredsamplelist]
        frequencies = [x for x in frequencies if x is not None]
        
        d = {
            'Marker': marker,
            'Frequency': np.median(frequencies),
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

def dataframe1(D0, D2, type1, var1, var2):
    
    # Loading dataset1
    
    D0, fullsamplelist, fullmarkerlist = set_data(D0)
    
    # Loading dataset2
    
    D2 = set_data(D2)[0]
    
    # Filtering the sample list
    
    if type1 == "all":
        filteredsamplelist = fullsamplelist
    else:
        filteredsamplelist = [sample for sample in fullsamplelist if (type1 == D0[sample].L1 or type1 == D0[sample].L2 or type1 == D0[sample].L3)]
    
    # Function to calculate PCR error ratio of var1 at a (sample, marker)
    
    def errorrat(sample, marker, var1):
        
        # cextractor
        
        def cextractor(data, sample, marker, var):
            try:
                return data[sample].Marker[marker].Allelec[var]
            except Exception:
                return None
            
        try:
           val1 = cextractor(D2, sample, marker, var1) * cextractor(D0, sample, marker, "0")
           val2 = cextractor(D0, sample, marker, var1) * cextractor(D2, sample, marker, "0")
           return val1 / val2
       
        except Exception:
           return None
   
    df = []
   
    for marker in fullmarkerlist:
        
        list1 = [(errorrat(x, marker, var1), errorrat(x, marker, var2)) for x in filteredsamplelist]
        list1 = [(x, y) for x, y in list1 if x is not None and y is not None]
        
        list2 = [x[0] for x in list1]
        list3 = [x[1] for x in list1]
        
        statistic, pvalue = stats.ttest_rel(list2, list3)
        
        avar1 = np.average(list2)
        avar2 = np.average(list3)
        
        d = {
            'Marker': marker,
            var1 : avar1,
            var2 : avar2,
            "p-value": pvalue
        }
        
        df.append(d)
        
    df = pd.DataFrame(df)
    
    return df

test = dataframe1("D0data", "D2data", "all", "-1", "1")

def plot3(D0, D2, type1, var1, var2):
    
    data = dataframe1(D0, D2, type1, var1, var2)
    
    # Reshape data for plotnine
    data_long = pd.melt(data, id_vars=['Marker', 'p-value'], value_vars=[var1, var2], 
                        var_name='Point', value_name='Value')
    
    # Convert Point to categorical to ensure separate colors
    data_long['Point'] = data_long['Point'].astype('category')
    
    print(data_long)
    
    # Define a threshold for significance (e.g., 0.05)
    significance_threshold = 0.05
    data['significance'] = data['p-value'].apply(lambda x: '*' if x < significance_threshold else '')
    
    # Plot
    plot = (
        p9.ggplot(data_long, p9.aes(x='Value', y='Marker', group='Marker')) +
        p9.geom_line() +  # Connects points for each marker
        p9.geom_point(p9.aes(color='Point'), size=3) +  # Dots for point 1 and point 2
        p9.geom_text(
        data=data[data['p-value'] < significance_threshold],
        mapping=p9.aes(x=1.2, y='Marker', label='significance'),  # Display significance symbol
        color='red',
        size=10) +
        p9.labs(x='Error Ratio', y='Marker', color='Point') +
        p9.theme_gray()
    )
    
    plot.show()
    
def dataframe2(D0, D2, type1, type2, var1):
    
    # Loading dataset1
    
    D0, fullsamplelist, fullmarkerlist = set_data(D0)
    
    # Loading dataset2
    
    D2 = set_data(D2)[0]
    
    # Filtering the sample list 1
    
    if type1 == "all":
        filteredsamplelist1 = fullsamplelist
    else:
        filteredsamplelist1 = [sample for sample in fullsamplelist if (type1 == D0[sample].L1 or type1 == D0[sample].L2 or type1 == D0[sample].L3)]
        
    # Filtering the sample list 2
    
    if type2 == "all":
        filteredsamplelist2 = fullsamplelist
    else:
        filteredsamplelist2 = [sample for sample in fullsamplelist if (type2 == D0[sample].L1 or type2 == D0[sample].L2 or type2 == D0[sample].L3)]
    
    # Function to calculate PCR error ratio of var1 at a (sample, marker)
    
    def errorrat(sample, marker, var1):
        
        # cextractor
        
        def cextractor(data, sample, marker, var):
            try:
                return data[sample].Marker[marker].Allelec[var]
            except Exception:
                return None
            
        try:
           val1 = cextractor(D2, sample, marker, var1) * cextractor(D0, sample, marker, "0")
           val2 = cextractor(D0, sample, marker, var1) * cextractor(D2, sample, marker, "0")
           return val1 / val2
       
        except Exception:
           return None
   
    df = []
   
    for marker in fullmarkerlist:
        
        list1 = [errorrat(x, marker, var1) for x in filteredsamplelist1]
        list1 = [x for x in list1 if x is not None]
        
        list2 = [errorrat(x, marker, var1) for x in filteredsamplelist2]
        list2 = [x for x in list2 if x is not None]
        
        avar1 = np.average(list1)
        avar2 = np.average(list2)
        
        statistic, pvalue = stats.ttest_ind(list1, list2)
        
        d = {
            'Marker': marker,
            type1 : avar1,
            type2 : avar2,
            "p-value": pvalue
        }
        
        df.append(d)
        
    df = pd.DataFrame(df)
    
    return df

def plot4(D0, D2, type1, type2, var1):
    
    data = dataframe2(D0, D2, type1, type2, var1)
    
    # Reshape data for plotnine
    data_long = pd.melt(data, id_vars=['Marker', 'p-value'], value_vars=[type1, type2], 
                        var_name='Point', value_name='Value')
    
    # Convert Point to categorical to ensure separate colors
    data_long['Point'] = data_long['Point'].astype('category')
    
    print(data_long)
    
    # Define a threshold for significance (e.g., 0.05)
    significance_threshold = 0.05
    data['significance'] = data['p-value'].apply(lambda x: '*' if x < significance_threshold else '')
    
    # Plot
    plot = (
        p9.ggplot(data_long, p9.aes(x='Value', y='Marker', group='Marker')) +
        p9.geom_line() +  # Connects points for each marker
        p9.geom_point(p9.aes(color='Point'), size=3) +  # Dots for point 1 and point 2
        p9.geom_text(
        data=data[data['p-value'] < significance_threshold],
        mapping=p9.aes(x=1.2, y='Marker', label='significance'),  # Display significance symbol
        color='red',
        size=10) +
        p9.labs(x='Error Ratio', y='Marker', color='Point') +
        p9.ggtitle("Average error ratio of" + var1) +
        p9.theme_gray()
    )
    
    plot.show()
    
dataframe2("D0data", "D2data", "CMMRD", "Control", "1")

plot4("D0data", "D2data", "CMMRD", "Control", "-1")
    
plot3("D0data", "D2data", "CMMRD", "-1", "1")

    


