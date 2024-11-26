# MasterProject

## 21/10/2024

### New Tasks

- [X] Make improvements to Powerpoint :
    - Motivation : Detection of +1 allele in samples, more frequent in MMRp than MMRd
    - Question : How did the +1 allele arise?
    - Hypothesis : 
        - +1 Allele is a PCR artefact
        - +1 Allele is caused by MMRd
            - unlikely due to previous results
        - +1 Allele is caused by MMR
            - aligns with previous result
            - MMR is known to be involved in trinucleotide expansion diseases
    - Significance :
        - Assuming Hyp 3 is true, +1 Allele frequency can be used as a marker to evaluate MMR function, esp in constitutional sample
        - Assuming Hyp 1 is true, +1 Allele frequency can be used to predict PCR error rate, and this can be used as an extra normalization factor
    - Ideas :
        - Different markers have different sensitivies to MMR effects
        - We also assume that there is a relationship between the +1 allele and 0 allele freqeuncy
            - Since we assume that the +1 allele is generated from insertion mutations occuring in the 0 allele 


## 25/10/2024

### To-Do

Mauro said that the new powerpoint was an improvement but some minor changes are necessary

- Powerpoint stuff

    - Remove "Caused by MMRd" hypothesis (it is already obvious that is not the case)
    - Explain the other hypotheses (what reason do I have to suggest these hypotheses?)
    - change "can be used as a marker" to "can be used as a feature"
    - "what about insertion?" should be on a higher (bullet point) level
    - Individual "to investigate" section for each set of graphs

- Data analysis

    - Blood / CMMRD samples
        - Compare $\mathtt{D_0(n)/D_0(0)}$ and $\mathtt{D_2(n)/D_2(0)}$ for $\mathtt{n\in \{-2,-1,1\}}$
        - Idea? : +1 counts ~ B_0 + B_1 (0 counts) + B_2 (with/without UMI)

## 4/11/2024

### To-Do

- [X] Create R script to make dataframe of allele counts of a marker at a sample, 
this time separated by SNP allele 
    - Finished 08/11/2024

- [X] Figure out code for plot 1 (VER 1)
    - Compare allele frequencies (which corresponds to mutation rate) between MMRp and MMRd samples on a per marker basis
    - Use Wilcoxon for statistical testing
    - x-axis is Wilcoxon effect size, y-axis is the different length alleles (-2,-1,0, and 1) color to indicate p-value
    - Finished 09/11/2024

- New Idea :
    - Microsatellites were chosen such that for each sample there exists a heterozygous SNP downstream of the microsatellite
    - PCR reads also includes this data, i.e. what is the downstream SNP of each read
    - Theoretically PCR reads at a sample-marker should be balanced between these two SNPs
    - An imbalance of counts between the two SNPs suggests homogeneiety in the initial sample
        - i.e. most of the cells in the initial sample have one allele preferentially deleted (via clonal expansion?)
    - [X] For each sample-marker, determine total counts (across all length alleles) for the most common and second most common SNP allele.
        - Add genotype of SNP1 (most common SNP) and SNP2 (second most common SNP) to the data
        - remove NAs, they are not alleles
        - calculate minor allele frequency (MAF) SNP2 / (SNP1 + SNP2)
        - Figure out why there is a missing marker in Sample2.commonSNPcounts
        - No missing marker on my end, the problem is on Mauro's end?
        - Sent the files again to Mauro just in case
        - Turns out the Ratios in LR17 were so low that they got filtered out
        - Finished 13/11/2024

## 11/11/2024

### To-Do

- [X] Change code for plot 1 (VER 2)
    - use AUC instead of Wilcoxon effect size
    - Finished 23/11/2024

- [X] Import Old CRC dataset and the 2 EC datasets
    - Finished 23/11/2024

- [X] Figure out statistics for Sample2.commonSNPcounts (plot1.plot2)
    - use MAF (minor/rare allele frequency) instead of ratios; but this shouldn't really make much of a difference
    - we checked that filtering out the total counts (including only observations with above 100) did not change anything 
    - Finished 15/11/2024

- [] Based on the results from plot2, figure out if, when going from D2 to D0, if the +1 allele frequency : 1) increases in MMRp, 2) decreases in MMRd, or 3) both (this will be plot 3)

- [X] Begin writing report
    - Finished 23/11/2024

## 25/11/2024

### Comments about report from Mauro, stuff to fix (main.20241123_msk.docx)

- [] Include citations in report
- [] It is more accurate to say that the MMR system identifies subtle changes to the DNA helical structure, this includes base mismatches, indel loops, etc...
- [] Give an overview (steps) of the MMR system before going into more details
- [] Why is there more focus on the MLH and MSH homologs? What about POLD/LIG1/POLE/EXO1 deficiencies?
- [] Endogenous alkylating agents -> DNA modification -> MMR system futile repair loop -> double strand break -> MMR is pro-mutagenic
    - quantify double strand break to measure MMR function (rocket and comet assays)
    - research the above topic a bit more
- [] No need to add "as determined by a supFG1 reporter assay"
- [] Revamp section on MMR in cancer
- [] MSI testing and IHC are only highly concordant in CRC, but less so in other types of cancer. MSI testing is bad in cancers with low levels of MSI
- [] MSI testing also uses gel electrophoresis, include that
- [] Challenges the "perception" is better than Challenges the "assumption"
- [] Add methodology section
- [] Explain AUC
- [] Why pattern is not observed in CMMRD? Hypothesize, then test :
    - Initial pool size of zero allele ~ +1 allele gained by PCR error
    - Thus, larger difference in zero allele frequency between MMRp and MMRd groups -> larger difference in +1 frequency in MMRp and MMRd group
    - Maybe this difference is high in CRC, lower in EC, and very low in CMMRD?
    - y variable : gain of +1 allele due to PCR error, given by the formula D0(+1) - D2(+1)
    - x variable : initial pool of 0 allele, given by D2(0) 

### Improvements to Code

- [] General cleanup
    - Dataset labels, change : "CMMRD" -> "CMMRD (blood)", "Jenny" -> "CRC1", "Man" -> "EC1", "Ohio" -> "EC2"



### Python Methods

- Methods :

    Transformations :
    - bySNP : collapses the length, Counts is now total of counts across all lengths at that SNP allele
        - noNA : removes "NA" SNP allele
        - MajorMinor : keeps only the 2 most major alleles
        - Freq1 : calculates the RelativeFrequency of each allele, count / total counts
                - MeltFreqs : melts frequencies

    - byLen : collapses the SNP, Counts is now total of counts across all SNPs at for that particular length allele
        - byMut : collapses into three alleles, -, 0, and +, the count of - is the sum of all the negative alleles, etc...
        - Freq1 : calculates the RelativeFrequency of each allele, count / total counts
            - MeltFreqs : melts frequencies
        - Freq2 : removes 0 alleles, then calculates the RelativeFrequency of each length allele
            - MeltFreqs : melts frequencies
        - MeltCounts : melts counts

    Plots

    - Plot1data_bydataset (ver2) : input is any byLen.MeltFreqs
        - for each Marker-UMI-Allele-Dataset -> a list of relativefrequencies in the MMRp and the MMRd group
    - Plot1data_bytype : input is any byLen.MeltFreqs
        - for each Marker-UMI-Allele-Type -> a list of relativefrequencies in the MMRp and the MMRd group
    - Plot1.plot1 : Filter by UMI, then plot
    - Plot1.plot2 : No UMI filter (UMI is in columns, dataset/type is in rows)

    Cleanup?
    - totalmarkercounts : In a sample-marker-UMI group, obtain the total number of reads mapping to that marker 
    - missingmarker : In a sample-UMI group, outputs list of markers with completely no data, i.e. missing markers
    - Filter1 : removes Sample-Marker-UMI group from the dataset if the totalmarkercount is below a certain threshold