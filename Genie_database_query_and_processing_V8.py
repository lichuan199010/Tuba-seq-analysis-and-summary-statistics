"""GENIE_data_processing.py """
'''
1. Input the panel data as dictionary
3. For each patient, find the corresponding sample and gene tested
4. Keep each patient if we found that is has the KRAS or EGFR mutations
5. For each gene, have a counter and count the number of mutations for each gene
Compile a list of tumor suppressor genes
Filter out the result for the tumor suppressor gene list
'''

## Import libraries ##
import pandas as pd
import numpy as np
import sys
import os
import collections

## Main working directory ##
p_d = 'C:/Users/lichuan/Downloads/CANCER/DATA/GENIE/GENIE_V8'
os.chdir(p_d)

# read in sample info
patient_info = "data_clinical_patient.txt"
sample_info = "data_clinical_sample.txt"
mutation_info = "data_mutations_extended.txt"

# readin a standard .txt dataframe
def df_readin(filename, use_cols=None, skip_rows=0):  
    df = pd.read_csv(
        filename,
        usecols=use_cols,
        skiprows=skip_rows,
        sep='\t'
        )
    return df

# somatic mutations from GENIE
genie_df = df_readin(mutation_info, use_cols=['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode', 'HGVSp_Short'])  # each row a tumor-mutation combo, multiple rows per tumour
genie_df = genie_df.rename(columns={'Hugo_Symbol': 'gene', 'Variant_Classification': 'category', 'Tumor_Sample_Barcode': 'sample', 'HGVSp_Short': 'mutant'})

# sample info - links SAMPLE_ID (sample) to PATIENT_ID (patient), gives cancer type
clinical_df_sample = df_readin(sample_info, skip_rows=4, use_cols=['PATIENT_ID', 'SAMPLE_ID', 'AGE_AT_SEQ_REPORT', 'ONCOTREE_CODE', 'SAMPLE_TYPE_DETAILED', 'SAMPLE_TYPE', 'SEQ_ASSAY_ID'])
clinical_df_sample = clinical_df_sample.rename(columns={'SAMPLE_ID': 'sample'})
genie_df_full = genie_df.merge(clinical_df_sample, on = ["sample"], how = "left")

panels = set(clinical_df_sample["SEQ_ASSAY_ID"])

assay_info = df_readin("assay_information.txt")
panels2 = set(assay_info["SEQ_ASSAY_ID"])

### let's focus on panels2 only!

# get a list for each panel - gene mapping
unique_panel_gene = genie_df_full.drop_duplicates(subset = ["gene", "SEQ_ASSAY_ID"])
panel_mut_dict = collections.defaultdict(set)

for idx, row in unique_panel_gene.iterrows():
    panel_mut_dict[row["SEQ_ASSAY_ID"]].add(row["gene"])

# get a list for each panel - sample mapping
unique_panel_sample = genie_df_full.drop_duplicates(subset = ["sample", "SEQ_ASSAY_ID"])
panel_sample_dict = collections.defaultdict(int)

for idx, row in unique_panel_gene.iterrows():
    panel_sample_dict[row["SEQ_ASSAY_ID"]] += 1

# whether or not we do luad only
clinical_df_sample = clinical_df_sample[clinical_df_sample['ONCOTREE_CODE'] == 'LUAD'] # filter for LUAD samples only 11107 samples

# run through each PATIENT_ID and save only those that have a single sample
# UNLESS there is a difference in age, sample type or panel
# keep: samples from younger patient, samples from primary tumour (in that order of priority)
patients_to_save = []
patients_to_save_d = {}
patients_to_chuck = []
samples_to_save = []
patient_grouped = clinical_df_sample.groupby('PATIENT_ID')

for each_patient, samples in patient_grouped:
    # check that they're all from LUAD
    cancers = samples['ONCOTREE_CODE'].tolist()
    if not all(x == 'LUAD' for x in cancers):
        raise Exception('cancer type filtering failed')
    # check that there is only one
    number_samples = len(samples)
    if number_samples == 1:
        patients_to_save.append(each_patient)
        sub_id = samples['sample'].tolist()[0]
        patients_to_save_d[sub_id] = each_patient
        continue
    if number_samples > 1:
        # check to see if there are any differences in the two entries
        dup_drop = samples.drop(columns=['sample']).drop_duplicates()
        if len(dup_drop) == 1:
            patients_to_chuck.append(each_patient)
            continue
        # if there are differences:
        # always take earlier sample if the ages differ of samples differ
        ages = samples['AGE_AT_SEQ_REPORT'].tolist()
        if not all(x==ages[0] for x in ages):
            samples = samples.sort_values(by=['AGE_AT_SEQ_REPORT']) # first entry from youngest patient
            sub_id = samples['sample'].tolist()[0]  # take first entry
            samples_to_save.append(sub_id)
            patients_to_save_d[sub_id] = each_patient
            continue
        # take primary tumour instead of metastasis if the samples are the same age
        sample_types = samples['SAMPLE_TYPE'].tolist()
        num_primary = sample_types.count('Primary')
        if num_primary == 1:  # if there is one primary sample
            to_take = samples[samples['SAMPLE_TYPE'] == 'Primary']
            sub_id = to_take['sample'].tolist()[0]
            samples_to_save.append(sub_id)
            patients_to_save_d[sub_id] = each_patient
            continue
        else:
            patients_to_chuck.append(each_patient)

# filter the sample file for only those patients that have a single LUAD sample
clinical_df_sample = clinical_df_sample.loc[clinical_df_sample['PATIENT_ID'].isin(patients_to_save)]
samples_to_use = list(set(clinical_df_sample['sample'].tolist()))
samples_to_use = samples_to_use + samples_to_save  #add in the ids from patients with multiple samples

genie_df_luad = genie_df.loc[genie_df['sample'].isin(samples_to_use)]
genie_df_luad = genie_df_luad.merge(clinical_df_sample, on = ["sample"], how = "left")
genie_df_luad.head()

# for each KRAS and EGFR samples, count the number of mutations belong to each sample
import collections


# find all samples that has mutation belonging to my category of interests
KRAS = ["p.G12", "p.G13", "p.Q61"]
EGFR = df_readin("EGFRList.txt")  # each row a tumor-mutation combo, multiple rows per tumour

subsetEGFR = genie_df_luad.loc[genie_df_luad["gene"].isin(["EGFR"]) & genie_df_luad["mutant"].isin(EGFR["Vars"])]
sampleEGFR = set(subsetEGFR["sample"])

def filterStr(s):
    s = str(s)
    return s.startswith("p.G12") or s.startswith("p.G13") or s.startswith("p.G61")
    
subsetKRAS = genie_df_luad.loc[genie_df_luad["gene"].isin(["KRAS"]) & genie_df_luad["mutant"].apply(filterStr)]
sampleKRAS = set(subsetKRAS["sample"])




### map subset to panel
sample_panel_dict = {}
for idx, row in subsetEGFR.iterrows():
    sample_panel_dict[row["sample"]] = row["SEQ_ASSAY_ID"]

for idx, row in subsetKRAS.iterrows():
    sample_panel_dict[row["sample"]] = row["SEQ_ASSAY_ID"]

## for each sample
# find its corresponding panel, add everything in the panel by 1 to total
# add positive findings by 1 to True
AllSamp = sampleKRAS.copy()
AllSamp = AllSamp.union(sampleEGFR)

# drop all cases with silent mutation and amplifications
KRAScntTot = collections.Counter()
KRAScntMut = collections.Counter()
EGFRcntTot = collections.Counter()
EGFRcntMut = collections.Counter()

for samp in AllSamp:
    
    subset = genie_df_luad[genie_df_luad["sample"] == samp]
    if samp in sampleEGFR and samp in sampleKRAS: continue
    if samp in sampleEGFR:
        tot = EGFRcntTot
        pos = EGFRcntMut
    elif samp in sampleKRAS:
        tot = KRAScntTot
        pos = KRAScntMut
    
    # get all samples belong to this samp
    
    ## In addition, we require that P53 genes were mutated
    # comment out this line if I want to have all mutants regardless of TP53 status
    if "TP53" not in subset["gene"].values: continue
    
    if samp not in sample_panel_dict: continue
    panel = sample_panel_dict[samp]
    if panel not in panel_mut_dict: continue
    
    for gene in panel_mut_dict[panel]:
        tot[gene] += 1
        
    for gene in set(subset["gene"]):
        pos[gene] += 1
        

freq = []
for key in pos:
    freq.append([pos[key]/tot[key], pos[key], tot[key], key])
    
freq.sort(reverse = True)
        
### read in a list of candidate genes
driver_genes_1 = 'Sanchez_Vega_Cell_2018_Supp_Tables_clean.xlsx'  # https://www.ncbi.nlm.nih.gov/pubmed/29625050
driver_genes_2 = 'mmc1.xlsx'  # https://www.ncbi.nlm.nih.gov/pubmed/29625053
# create lists of passengers and drivers from two sources
drivers1_df = pd.read_excel(driver_genes_1, index_col='Gene')
drivers2_df = pd.read_excel(
    driver_genes_2,
    skiprows = [0,1,2],
    sheet_name = 'Table S1',
    usecols=['Gene', 'Cancer', 'Tumor suppressor or oncogene prediction (by 20/20+)', 'KEY'],
    index_col='Gene')
drivers2_df = drivers2_df.rename({'Tumor suppressor or oncogene prediction (by 20/20+)': 'TS_or_OG'}, axis=1)

# tsg and og from paper 1
drivers1 = list(set(drivers1_df.index.tolist()))
drivers1_tsg = list(set(drivers1_df[drivers1_df['TS_or_OG'] == 'TSG'].index.tolist()))
drivers1_og = list(set(drivers1_df[drivers1_df['TS_or_OG'] == 'OG'].index.tolist()))

# tsg and og from paper 2
drivers2 = list(set(drivers2_df.index.tolist()))
drivers2_tsg = list(set(drivers2_df[(drivers2_df['TS_or_OG'] == 'tsg') | (drivers2_df['TS_or_OG'] == 'possible tsg')].index.tolist())) 
drivers2_og = list(set(drivers2_df[(drivers2_df['TS_or_OG'] == 'oncogene') | (drivers2_df['TS_or_OG'] == 'possible oncogene')].index.tolist())) 

#TSG are across any cancer type, combine from both papers
combined_drivers = list(set(drivers1 + drivers2))
combined_tsg = set(drivers1_tsg + drivers2_tsg)
combined_og = list(set(drivers1_og + drivers2_og))

len(combined_drivers)
len(combined_og)
len(combined_tsg)

### plot correlation for tumor suppressor genes

tsg_to_remove = ['MET', 'WT1', 'JAK1', 'MAP3K1', 'NRAS', 'PIK3R3']  # manual curation for non-TS

combined_tsg = [x for x in combined_tsg if x not in tsg_to_remove]  # remove these


krastot = []
egfrtot = []
krascnt = []
egfrcnt = []
    
for gene in combined_tsg:
    krastot.append(KRAScntTot[gene] if gene in KRAScntTot else 0)
    egfrtot.append(EGFRcntTot[gene] if gene in EGFRcntTot else 0)
    krascnt.append(KRAScntMut[gene] if gene in KRAScntMut else 0)
    egfrcnt.append(EGFRcntMut[gene] if gene in EGFRcntMut else 0)
    
output = pd.DataFrame(list(zip(combined_tsg, krastot, krascnt, egfrtot, egfrcnt)), columns = ["gene", "KrasTot", "KrasMut", "EGFRTot", "EGFRMut"])
output = output.loc[output["KrasTot"] > 0] 

output.to_csv("CntsForTSGsP53Mut.txt", index = False, sep = "\t")
