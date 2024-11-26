# MasterProject

### What is in RawData?

- 3 CRC
    - CRC1
        - genotype : NEWCRC/20170501
        - phenotype : NEWCRC/Tumours.20190923.csv
    - CRC2
        - genotype : NEWCRC/20170809
        - phenotype : NEWCRC/Tumours.20190923.csv
    - Jenny
        - genotype : OGCRC/UMImsiMNRsJenny210421
        - phenotype : OGCRC/Data.Jenny.Coding.csv

- 1 CMMRD
    - CMMRD
    	- genotype : CMMRD/byUMIGenotypesRearranged
	    - phenotype : CMMRD/SampleData.20220204.v2.rds

- 2 EC
    - Man
        - genotype : EC/ManUMI/20220622
        - phenotype : EC/ManUMI/Manchester.EC.Samples.20240111.csv
        - Notes : Samples were renamed accordingly to match "EC/Man/20220425.24" (scripts/Renaming.R)
    - Ohio
        - genotype : EC/OhioUMI/20220318, EC/UMIOhio/20220401 
        - phenotype : EC/OhioUMI/Ohio.SampleData.20240111.csv

Marker Lengths : MarkerLengths.csv

### Pre-Processing :

All datasets are grouped into one dataframe

- withSNP
    - Columns : SampleName, Marker, Allele (length), Allele (SNP), Counts, UMI, L1, L2, L3, Dataset, Type
    - Scripts
        - scripts/Main.R + ReadFolder2.R
    - csv output : Data/withSNP/All.csv

This csv is then loaded as a python object, I have made various python "methods" to transform the dataframe in various ways depending on the purpose

## Methods : 

### withoutSNP1

Procedure :
    - Grouping : 'SampleName', 'Allele (length)', 'Marker', 'UMI', 'L1', 'L2', 'L3', 'Dataset', 'Type'
    - "Counts" (new column) : the sum of all the counts across each "Allele (SNP)", this becomes the new counts for each group

<!--
Groups all observations (rows) which have the same values in the "Grouping" column together
-->

### withoutSNP2

From withoutSNP1, create a new column called "Total", the total counts across "Allele (length)" for each group (marker-sample combination)

### totalmarkercounts1

This gives the total number of reads/counts observed in a particular sample-marker-UMI

### totalmarkercounts2

Same as totalmarkercounts, but filters are present

filter parameters : D0, D1, D2

Filter :

Total < D0 in the D0 UMI group
Total < D1 in the D1 UMI group
Total < D2 in the D2 UMI group

### fil

Parameters : D0, D1, D2

The original dataset, but with observations from Sample-Marker-UMI groups present in "totalmarkercounts2" removed