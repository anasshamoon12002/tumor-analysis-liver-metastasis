
# %%
import pandas as pd 
import numpy as np


# %%

pd.set_option('display.max_columns', None)

df_1 = pd.read_excel("../dataset/d1.xlsx")

print(list(df_1))

df_1.head()


# %%
df_2 = pd.read_excel("../dataset/d2.xlsx")

print(list(df_2))
print(df_2.shape)

df_2.head()

# %% [markdown]
# ## Get patients with surgery on liver and find unique biomarkers

# %%
liver_patients = df_1.loc[df_1["sites of surgery"].isin([1, "1;2", "1;4", "1;5"])]["AccessionNumber"]
all_patients = df_1['AccessionNumber']
print(len(liver_patients))
all_biomarkers = df_2.loc[df_2["AccessionNumber"].isin(all_patients)]["Biomarker"]
biomarkers = df_2.loc[df_2["AccessionNumber"].isin(liver_patients)]["Biomarker"]
unique_biomarkers = list(set(biomarkers))
all_biomarkers = list(set(all_biomarkers))
len(unique_biomarkers)
len(all_biomarkers)
# unique_biomarkers

# %%
unique_biomarkers.insert(0, 'AccessionNumber')
all_biomarkers.insert(0, 'AccessionNumber')
new_df = pd.DataFrame(columns=unique_biomarkers)
AllPatients = pd.DataFrame(columns=all_biomarkers)
all_biomarkers.remove('AccessionNumber')
unique_biomarkers.remove('AccessionNumber')

# %%
new_df["AccessionNumber"] = liver_patients
AllPatients["AccessionNumber"] = all_patients
AllPatients
# %%
# for index, row in df_2.iterrows():
    # print (row)

# %%
countGenes = df_2['Biomarker'].value_counts()
commonGenes = countGenes[countGenes > 592]
print(len(commonGenes))
#%%

df_2

#%%
#print(df_2[df_2['Technology'] == 'CNA']['TestResult'].unique())
#print(df_2[df_2['Technology'] == 'CNA']['Conclusion'].unique())

#TestResult
#'Amplification Not Detected' 'Indeterminate' 'Amplified' 'Intermediate']
#        0                            0            1              0

df_2.loc[df_2['TestResult'] == 'Amplification Not Detected', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Indeterminate', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Intermediate', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Amplified', 'TestResult'] = 1
df_2
#%%
# for access_num in liver_patients:

print(df_2[df_2['Technology'] == 'NGS Q3']['TestResult'].unique())

print(df_2[df_2['Technology'] == 'NGS Q3']['Conclusion'].unique())

#MUTATION AMPL RESULT
# 0        0      0
# 1        0      1
# 0        1      2
# 1        1      3

#not present -> -1



#TECHNOLOGY ==> CNA
#'Amplification Not Detected' 'Indeterminate' 'Amplified' 'Intermediate']
#        0                            0            1              0



#TECHNOLOGY ==> NGS Q3
# 'variantnotdetected' #0
# 'indeterminate' #0 
# 'Wild Type' #0
# 'Indeterminate' #0
# 'Intermediate' #0
# 'Stable' #0
# 'Deficient' #0

# 'Mutated, Variant of Unknown Significance'  #1
# 'Mutated, Presumed Pathogenic' #1
# 'variantdetected' #1
# 'Mutated, Pathogenic' #1 
# 'Mutated, Presumed Benign' #1
# 'Mutated, Unclassified' #1
# 'Low' #1
# 'High' #1



df_2.loc[df_2['TestResult'] == 'variantnotdetected', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'indeterminate', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Wild Type', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Indeterminate', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Stable', 'TestResult'] = 0
df_2.loc[df_2['TestResult'] == 'Deficient', 'TestResult'] = 0

df_2.loc[df_2['TestResult'] == 'Mutated, Pathogenic', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'Mutated, Variant of Unknown Significance', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'Mutated, Presumed Pathogenic', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'Mutated, Presumed Benign', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'Low', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'High', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'Mutated, Unclassified', 'TestResult'] = 1
df_2.loc[df_2['TestResult'] == 'variantdetected', 'TestResult'] = 1


'''Is MSI Stable high or low?
CRCs can be classified into 3 groups according to the MSI status: 
MSI-high (MSI-H), which exhibit ≥ 30 to 40% microsatellite marker instability,
MSI-Low (MSI-L), which exhibit instability at < 30 to 40% of loci,
and microsatellite stable (MSS), which exhibit no unstable markers.'''
# %%
df_2
# %%
print(df_2[df_2['Technology'] == 'NGS Q3']['TestResult'].unique())
print(df_2[df_2['Technology'] == 'CNA']['TestResult'].unique())
# %%
for col in new_df.columns:
    new_df[col].values[:] = -1

new_df

for col in AllPatients.columns:
    AllPatients[col].values[:] = -1
AllPatients["AccessionNumber"] = all_patients
AllPatients

#%%

df_NGS = df_2[df_2['Technology'] == 'NGS Q3']
#df_NGS = df_NGS[~df_NGS['NGS_PercentMutated'].isna()]
df_NGS

#%%
## PROBLEM TO SOLVE: there are duplicated biomarkers for patient because they used two different technologies ##
## the following code just overwrites the first occurence with the second ##

#pseudocode



for index, row in df_NGS.iterrows():
  #read the biomarker of that row
  Biomarker = row['Biomarker']
  #read the result of the biomarker in that row
  mutationResult = row['NGS_PercentMutated']
  #Technology
  if( np.isnan(mutationResult)):
    mutationResult = 0;
    
  #read patient
  AccessN = row['AccessionNumber']
  #if the patient has changed we need to update a new row in the other dataset
  #update the right row, the right column -> new_df[ROW_NUMBER][COLUMN_NAME]
  AllPatients.loc[ (AllPatients['AccessionNumber'] == AccessN), Biomarker] = mutationResult
# %%
AllPatients
#%%

#MUTATION AMPL RESULT
# 0        0      0
# 1        0      1
# 0        1      2
# 1        1      3

#not present -> -1

#%%
AllPatients[]
# %%
AllPatientsNaN = AllPatients.copy()
AllPatientsNaN.replace(-1, np.NaN, inplace=True)
AllPatientsNaN
#%%
if -1 in AllPatientsNaN.values:
   print("Yes")
# %%
CommonGenesPatients
# %%
AllPatientsNaN=AllPatientsNaN.dropna(axis='columns')
# %%
print(AllPatientsNaN.max())
AllPatientsNaN.loc[:, (AllPatientsNaN != -1).any(axis=0)]

# %%
AllPatientsNaN
# %%
