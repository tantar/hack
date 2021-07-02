# Loads libraries required for the app
from google.colab import drive
import streamlit as st
import pandas as pd
import numpy as np
import datetime as dt
import os
import altair as alt
import plotly.graph_objects as go
import plotly.express as px
import tableone as tab1
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# # The current version relies on data stored in Tarek's google drive
# drive.mount('/content/gdrive')
# datafolder = f'/content/gdrive/MyDrive/Hackathon/data'
# workfolder = f'/content/gdrive/MyDrive/Hackathon/'

# os.chdir(workfolder)\

#Set up the general format for the app 
st.title("GP2 Data Visualization Tool")
st.sidebar.title("Options")

CST = st.selectbox("Cohort Selection Type", ["Preset Cohorts", "Choose-Your-Own Cohorts"])

if CST == "Preset Cohorts":
  select_data = st.multiselect("Select Datasets",
        ['PPMI_PD', 'PPMI_HC', 'PPMI_PRODROMA', 'PPMI_SWEDD', 'PPMI_GENUN','PPMI_GENPD', 'PPMI_REGUN','PPMI_REGPD', 
       'BioFIND_PD', 'BioFIND_HC', 
       'PDBP_202_Control', 'PDBP_202_Other', 'PDBP_202_PD', 'PDBP_204_Control',
       'PDBP_204_Not Reported', 'PDBP_204_PD', 'PDBP_205_Control', 'PDBP_205_PD',
       'PDBP_206_Control', 'PDBP_206_PD', 'PDBP_207_Control', 'PDBP_207_Other', 'PDBP_207_PD',
       'PDBP_208_Control', 'PDBP_208_PD', 'PDBP_214_Control', 'PDBP_214_Other', 'PDBP_214_PD',
       'PDBP_220_PD', 'PDBP_224_Control', 'PDBP_224_PD', 'PDBP_225_Control', 'PDBP_225_PD',
       'PDBP_229_Control', 'PDBP_229_Other', 'PDBP_229_PD', 'PDBP_231_Other', 'PDBP_231_PD',
       'PDBP_232_Other', 'PDBP_232_PD', 'PDBP_233_Other', 'PDBP_233_PD', 'PDBP_235_Control',
       'PDBP_235_Not Reported', 'PDBP_235_Other', 'PDBP_235_PD', 'PDBP_236_Other', 'PDBP_236_PD',
       'PDBP_245_Other', 'PDBP_245_PD', 'PDBP_247_PD', 'PDBP_254_Control', 'PDBP_254_Other',
       'PDBP_254_PD',],
        ['PPMI_PD', 'BioFIND_PD', 'PDBP_202_PD','PDBP_204_PD', 'PDBP_205_PD',  
          'PDBP_206_PD', 'PDBP_207_PD', 'PDBP_208_PD', 'PDBP_245_PD','PDBP_247_PD'])
else:
  select_data = st.multiselect("What study arms are in your cohort of interest (COI)?",   
        ['PPMI_PD', 'PPMI_HC', 'PPMI_PRODROMA', 'PPMI_SWEDD', 'PPMI_GENUN','PPMI_GENPD', 'PPMI_REGUN','PPMI_REGPD', 
       'BioFIND_PD', 'BioFIND_HC', 
       'PDBP_202_Control', 'PDBP_202_Other', 'PDBP_202_PD', 'PDBP_204_Control',
       'PDBP_204_Not Reported', 'PDBP_204_PD', 'PDBP_205_Control', 'PDBP_205_PD',
       'PDBP_206_Control', 'PDBP_206_PD', 'PDBP_207_Control', 'PDBP_207_Other', 'PDBP_207_PD',
       'PDBP_208_Control', 'PDBP_208_PD', 'PDBP_214_Control', 'PDBP_214_Other', 'PDBP_214_PD',
       'PDBP_220_PD', 'PDBP_224_Control', 'PDBP_224_PD', 'PDBP_225_Control', 'PDBP_225_PD',
       'PDBP_229_Control', 'PDBP_229_Other', 'PDBP_229_PD', 'PDBP_231_Other', 'PDBP_231_PD',
       'PDBP_232_Other', 'PDBP_232_PD', 'PDBP_233_Other', 'PDBP_233_PD', 'PDBP_235_Control',
       'PDBP_235_Not Reported', 'PDBP_235_Other', 'PDBP_235_PD', 'PDBP_236_Other', 'PDBP_236_PD',
       'PDBP_245_Other', 'PDBP_245_PD', 'PDBP_247_PD', 'PDBP_254_Control', 'PDBP_254_Other',
       'PDBP_254_PD',],
        ['PPMI_PD', 'BioFIND_PD', 'PDBP_202_PD','PDBP_204_PD', 'PDBP_205_PD',  
          'PDBP_206_PD', 'PDBP_207_PD', 'PDBP_208_PD', 'PDBP_245_PD','PDBP_247_PD'])
  COI = st.text_input("Cohort Name")
  OD = st.selectbox("What will this cohort be compared to?", ["No Comparison","All Others", "All Other PD", "All Other HC"])
                               

###----------------------------------------------------------------------------------------------###
###This section of code performs data cleaning for output in streamlit
###----------------------------------------------------------------------------------------------'''
cols = ["participant_id","visit_name","visit_month",'date_visit', 'date_baseline',
       'age_at_baseline',"age_at_diagnosis", 'date_birth', 'sex', 'ethnicity', 'race',
       'education_level','education_years', 'family_hx_1st', 'family_hx_2nd',
       'family_hx_generations', 'study_arm', 'cohort', 'Phenotype',
       'primary_diagnosis','handedness', 'dominant_side','BMI', 
        
       'hx_hypertension','hx_hyperlipidemia', 'hx_diabetes', 'hx_heart_disease',
       'hx_dementia_mci', 'hx_stroke', 'hx_autoimmune_disease',
       'hx_depression', 'hx_bipolar', 'hx_restless_legs_syndrome',
       'hx_head_trauma', 'hx_constipation', 'hx_sleep', 'smoking_status',
       
       'levodopa_mg_daily', 'ledd_daily', 'levodopa_use','dopamine_agonist_use',
       'anticholinergics_use','cholin_esterase_inhibitor_use', 'anti_depressant_use',
       'sleeping_pills_use', 'stomach_medicines_use', 'laxatives_use',
       'NSAIDs_use', 'morphines_use', 'anti_psychotics_use', 'insulin_use',
       'hormone_replacement_therapy', 'hoehn_and_yahr_stage',
        
        "mds_updrs_part_i_summary_score","mds_updrs_part_ii_summary_score",
        "mds_updrs_part_iii_summary_score","mds_updrs_part_iv_summary_score",
        'depress_test_score','depress_test_name','gds15_total_score','moca_total_score',
        'rbd_summary_score','ess_total_score','schwab_england_pct_adl_score', 
        'smell_test_results', 'smell_test_name','smell_test_score','pdq39_mobility_score',
        'pdq39_adl_score', 
        'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
        'pdq39_cognition_score', 'pdq39_communication_score',
        'pdq39_discomfort_score', 'mod_schwab_england_pct_adl_score',

        'dx_essential_bradykinesia','dx_essential_rest_tremor',
        'dx_essential_rigidity', 'clinical_dx_depression',
        'clinical_dx_insomnia','hx_dementia','insuline_use'

        ]

# PPMI
study_arms = [x.replace('PPMI_', '') for x in select_data if 'PPMI' in x]
d = pd.read_csv("data/PPMI_processed.csv")
d1 = d[d.study_arm.isin(study_arms)]

if CST == "Preset Cohorts":
  d1["cohort"] = "PPMI"
else:
  d1["cohort"] = COI

  if OD == "All Others":
    d2 = d[~d.study_arm.isin(study_arms)]
    d2["cohort"] = OD
    d1 = d1.append(d2)
  elif OD == "All Other PD":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("PD"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)    
  elif OD == "All Other HC":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("HC"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)  

d1 = d1.loc[:,d1.columns.isin(cols)]

t1 = d1.copy()

# PDBP
study_arms = [x.replace('PDBP_', '') for x in select_data if 'PDBP' in x]
d = pd.read_csv("data/PDBP_processed.csv")
d1 = d[d.study_arm.isin(study_arms)]

if CST == "Preset Cohorts":
  d1["cohort"] = "PDBP"
else:
  d1["cohort"] = COI

  if OD == "All Others":
    d2 = d[~d.study_arm.isin(study_arms)]
    d2["cohort"] = OD
    d1 = d1.append(d2)

  elif OD == "All Other PD":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("PD"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)   

  elif OD == "All Other HC":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("HC"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)  

d1 = d1.loc[:,d1.columns.isin(cols)]
t2 = d1.copy()

# BioFIND
study_arms = [x.replace('BioFIND_', '') for x in select_data if 'BioFIND' in x]
study_arms = ['Control' if arm=='HC' else arm for arm in study_arms]
d = pd.read_csv("data/BIOFIND_processed.csv")
d1 = d[d.study_arm.isin(study_arms)]

if CST == "Preset Cohorts":
  d1["cohort"] = "BIOFIND"
else:
  d1["cohort"] = COI

  if OD == "All Others":
    d2 = d[~d.study_arm.isin(study_arms)]
    d2["cohort"] = OD
    d1 = d1.append(d2)
  elif OD == "All Other PD":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("PD"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)    
  elif OD == "All Other HC":
    d2 = d[(~d.study_arm.isin(study_arms))&(d.study_arm.str.contains("HC"))]
    d2["cohort"] = OD
    d1 = d1.append(d2)  



d = d[d.study_arm.isin(study_arms)]

d1 = d1.loc[:, d1.columns.isin(cols)].copy()
t3 = d1.copy()

# combine all
d = pd.concat([t1,t2,t3], ignore_index=True, axis=0)

if d.empty:
  st.warning("Please Select a Cohort")
  st.stop()

if st.sidebar.checkbox("Show Substudy Info"):
  st.sidebar.markdown("""
  We can use PPMI, BioFIND and PDBP in the current version. Please pick the study_arms of interest.    
  **PPMI substudies**
  * PPMI_PD/HC/SWEDD: original cohorts - PD patients, SWEDDS and contros
  * PPMI_PRODROMAL - people with prodromal symptoms (hyposmia or DATSCAN deficit)
  * PPMI_GENPD/GENUN: genetic cohort - affected/unaffected
  * PPMI_REGPD/REGUN: genetic registry - affected/unafected

  **PDBP substudies**   
  * 202: N = 200,  1-year, PD and other mimicking diseases plus controls
  * 204: N = 600,  cross-sectional, PD and controls
  * 205: N = 120, 5-year, 
  * 206: N = 240, 4-year
  * 207: N = 270, 3-year. PD and other mimicking diseases
  * 208: N = 75, 3-year
  """)
  


cats = ['sex','Phenotype','race','education_level','primary_diagnosis',
       'handedness', 'dominant_side','hx_depression','family_hx_1st', 'family_hx_2nd',
       'family_hx_generations','ethnicity','cohort','study_arm',

       'hx_hypertension','hx_hyperlipidemia', 'hx_diabetes', 'hx_heart_disease',
       'hx_dementia_mci', 'hx_stroke', 'hx_autoimmune_disease',
       'hx_depression', 'hx_bipolar', 'hx_restless_legs_syndrome',
       'hx_head_trauma', 'hx_constipation', 'hx_sleep', 'smoking_status',
       
       'levodopa_use',
       'dopamine_agonist_use', 'anticholinergics_use',
       'cholin_esterase_inhibitor_use', 'anti_depressant_use',
       'sleeping_pills_use', 'stomach_medicines_use', 'laxatives_use',
       'NSAIDs_use', 'morphines_use', 'anti_psychotics_use', 'insulin_use',
       'hormone_replacement_therapy', 
        
        'hoehn_and_yahr_stage','HY','HY3','depress_test_name','smell_test_name',

        'dx_essential_bradykinesia','dx_essential_rest_tremor',
        'dx_essential_rigidity','bradykinesia','resting_tremor','rigidity',
        'clinical_dx_depression','clinical_dx_insomnia','hx_dementia','insuline_use']

nnml = ["BL_durdx","daysB","daysDx","yearsDx","yearsB",
  
        'age_at_baseline','age_at_diagnosis','education_years','BL_durdx',
        'levodopa_mg_daily', 'ledd_daily','BMI','smell_test_score',

        "mds_updrs_part_i_summary_score","UPDRS1","mds_updrs_part_ii_summary_score","UPDRS2",
        "mds_updrs_part_iii_summary_score","UPDRS3","mds_updrs_part_iv_summary_score","UPDRS4",
        
        'depress_test_score','gds15_total_score','moca_total_score','MoCA',
        'rbd_summary_score','ess_total_score','pdq39_mobility_score', 'pdq39_adl_score',
       'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
       'pdq39_cognition_score', 'pdq39_communication_score','schwab_england_pct_adl_score', 
       'pdq39_discomfort_score','mod_schwab_england_pct_adl_score']



#Creating checkboxes and slide bars for dataset, diagram choice, and variable of choice
viz = st.sidebar.selectbox("Choose a Visualization", ["Tables","Line Plot","Bar Graph",
                                                         "Baseline Histogram","Scatter Plot",
                                                         "Sankey Diagram","Kaplan-Meier Survival"])

d.loc[:,"HY"] = d.hoehn_and_yahr_stage
d = d.rename({"mds_updrs_part_i_summary_score" : "UPDRS1",
              "mds_updrs_part_ii_summary_score": "UPDRS2",
              "mds_updrs_part_iii_summary_score":"UPDRS3",
              "mds_updrs_part_iv_summary_score":"UPDRS4", 
              "moca_total_score":"MoCA",
              'dx_essential_bradykinesia':'bradykinesia', 
              'dx_essential_rest_tremor':'resting_tremor',
              'dx_essential_rigidity':'rigidity'}, axis = 1)

#Essential Time Variables are Created
d.loc[:,"BL_durdx"] = d.age_at_baseline - d.age_at_diagnosis
d.loc[:,"daysB"] = (d.date_visit - d.date_baseline).fillna(d.visit_month * 30)
d.loc[:,"daysDx"] = (d.daysB + d.BL_durdx * 365.25).round(2)
d.loc[:,"yearsDx"] = (d.daysDx / 365.25).round(1)
d.loc[:,"yearsB"] = (d.daysB / 365.25).round(1)

yearsDx = "Years Since Diagnosis"


if "HY" in d.columns:
  d.loc[:,"HY3"] = ["HY<3" if i < 3 
                          else "HY>=3" for i in d.HY]

d.loc[:,"age_bin"] = ["<50" if i < 50
                            else ">50" for i in d.age_at_diagnosis]


for var in ["UPDRS1","UPDRS2","UPDRS3","UPDRS4"]:
  max = d.loc[:,var].max()
  d.loc[:,(var + "_bin")] = [("<" + str(round((max * .2),0)))  if i < round((max * .2),0)
                          else ("<" + str(round((max * .4),0))) if i < round((max * .4),0)
                          else ("<" + str(round((max * .6),0)))  if i < round((max * .6),0)
                          else ("<" + str(round((max * .8), 0)))  if i < round((max * .8), 0)
                          else ("<" + str(round(max, 0)))  if i <= round(max, 0)
                          else np.nan for i in d.loc[:,var]]

                        
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the scatterplots
###----------------------------------------------------------------------------------------------'''
if viz == "Scatter Plot":
  xax = st.sidebar.selectbox("Select an x-axis variable", d.columns[d.columns.isin(nnml)])
  yax = st.sidebar.selectbox("Select a y-axis variable", d.columns[d.columns.isin(nnml)])

  if st.sidebar.checkbox("Include only Baseline?"):
    d = d[d.visit_month == 0]

  d = d.dropna(subset = [xax,yax])

  tl = None
  if st.sidebar.checkbox("Add Trendline?"):
    tl = st.sidebar.selectbox("Select Trendline Type", ["ols", "lowess"])
    

  if st.sidebar.checkbox("Stratify Based on Third Variable?"):
    cs = st.sidebar.selectbox("Select a stratifying variable", d.columns[d.columns.isin(cats)])
    d = d.dropna(subset = [cs])    
    fig = px.scatter(d, x = xax, y = yax, color = cs, opacity = 0.5, trendline = tl)

  else:
    fig = px.scatter(d, x = xax, y = yax, opacity = 0.5, trendline = tl)


  st.write(fig)

###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the demographic summarization table
###----------------------------------------------------------------------------------------------'''
if viz == "Tables":

  foc = st.sidebar.selectbox("Select a Focus", ["Overall","Demographic","Motor","Non-Motor","Choose-Your-Own"])
  gby = st.sidebar.selectbox("Stratifying Variable", ["None","study_arm","cohort",'Phenotype'])

  if foc == "Overall":
    # columns to summarize
    columns = ['sex','Phenotype','race','education_level','primary_diagnosis','age_at_baseline',
              'age_at_diagnosis','education_years','BL_durdx','UPDRS1', 'UPDRS2', 'UPDRS3', 'UPDRS4']

    # columns containing categorical variables
    categorical = [i for i in cats if i in columns] 

    # non-normal variables
    nonnormal = [i for i in nnml if i in columns]
  
  if foc == "Demographic":
    # columns to summarize
    columns = ['sex','Phenotype','race','education_level','primary_diagnosis','age_at_baseline',
              'age_at_diagnosis','education_years','BL_durdx']

    # columns containing categorical variables
    categorical = [i for i in cats if i in columns] 

    # non-normal variables
    nonnormal = [i for i in nnml if i in columns]

  if foc == "Motor":
    # columns to summarize
    columns = ['BL_durdx', 'UPDRS3', 'UPDRS4','bradykinesia', 'resting_tremor',
       'rigidity']
    
    # columns containing categorical variables
    categorical = [i for i in cats if i in columns] 

    # non-normal variables
    nonnormal = [i for i in nnml if i in columns]
  
  if foc == "Non-Motor":
    columns = ["gds15_total_score",'depress_test_name','depress_test_score',"MoCA",
               'clinical_dx_depression','clinical_dx_insomnia',
               'pdq39_mobility_score', 'pdq39_adl_score',
               'smell_test_name','smell_test_score',
               'hx_dementia',
               'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
               'pdq39_cognition_score', 'pdq39_communication_score',
               'pdq39_discomfort_score',]
    
    # columns containing categorical variables
    categorical = [i for i in cats if i in columns] 

    # non-normal variables
    nonnormal = [i for i in nnml if i in columns]

  if foc == "Choose-Your-Own":
      columns = st.multiselect("Select Variables", list(d.columns[2:]),['sex','Phenotype','age_at_baseline'])

      # columns containing categorical variables
      categorical = [i for i in cats if i in columns] 

      # non-normal variables
      nonnormal = [i for i in nnml if i in columns]
  
  # limit the binary variable "death" to a single row
  # limit = {"death": 1}

  # set the order of the categorical variables
  order = dict()
  for col in columns:
    if "Yes" in d.loc[:, col].unique():
      order[col] = ["Yes","No","Unknown"]
      
  # alternative labels
  #labels={'death': 'Mortality'}

  # set decimal places for age to 0
  #decimals = {"Age": 0}

  # optionally, a categorical variable for stratification
  #groupby = ['visit_month']

  # rename the death column
  #labels={'death': 'Mortality'}

  # display minimum and maximum for listed variables
  #min_max = ['Height']

  if not columns:
      st.warning("No Variable Has Been Selected")
  elif gby == "None":
    dt = tab1.TableOne(d.drop_duplicates(subset = "participant_id",keep="first"), columns=columns, categorical=categorical,
                      nonnormal=nonnormal, order = order)
    st.write(dt)
  else:
    dt = tab1.TableOne(d.drop_duplicates(subset = "participant_id",keep="first"), columns=columns, categorical=categorical,
                      nonnormal=nonnormal, groupby = gby, order = order)
    st.write(dt)

  


###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing histograms for baseline variables
###----------------------------------------------------------------------------------------------'''
if viz == "Baseline Histogram":
  var = st.sidebar.selectbox("Pick a variable to visualize", d.columns[3:])
  counts = st.sidebar.checkbox("Convert to Percent Histogram")
  d = d[d.visit_month == 0]

  if var not in d.columns:
    st.write("This variable is not supported in the selected dataset")

  strat = st.sidebar.selectbox("Pick a Stratifier", d.columns[7:])

  if var in nnml:
    nb = st.sidebar.slider("Number of Bins", 1, 25, 50, step = 1)
  else:
    nb = 1
  
  if strat == "None":
    fig = px.histogram(d, x = var, y = "participant_id", nbins = nb)

    if counts == True:
          fig = px.histogram(d, x = var, y = "participant_id", barnorm = "percent",
                          nbins = nb)
          
  if counts == True:
    d = d.dropna(subset = [var, strat])
    fig = px.histogram(d, x = var, y = "participant_id", barnorm = "percent",color = strat, nbins = nb)   
  else:
    d = d.dropna(subset = [strat])
    fig = px.histogram(d, x = var, y = "participant_id", color = strat, nbins = nb)
  
  st.write(fig)
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the new stacked bar graphs
###----------------------------------------------------------------------------------------------'''
if viz == "Bar Graph":
  var = st.sidebar.selectbox("Pick a y-axis", ["HY3","UPDRS1","UPDRS2","UPDRS3","UPDRS4",
                                               'bradykinesia', 'resting_tremor','rigidity'])

  if var in ["UPDRS1","UPDRS2","UPDRS3","UPDRS4"]: 
    var = var + "_bin"

  bl = st.sidebar.selectbox("Pick an x-axis",["yearsDx","yearsB","visit_month",
                                              "HY3","UPDRS1","UPDRS2","UPDRS3",
                                              "UPDRS4"])


  strat = st.sidebar.selectbox("Pick a Stratifier", ["None","sex","age_bin"])
  counts = st.sidebar.selectbox("Pick a Histogram Type", ["Count","Percent"])
  nb = st.sidebar.slider("Number of Bins", 1, 25, 50, step = 1)
  d = d.dropna(subset = [var])

  if bl in ["yearsDx","yearsB"]:
    d = d.drop_duplicates(subset = ["participant_id",var])


  if var not in d.columns:
    st.write("This variable is not supported in the selected dataset")

  if counts == "Percent":
    if strat == "None":
      fig = px.histogram(d, x = bl, y = "participant_id", color = var, barnorm = "percent",
                          nbins = nb)
    else:
      d = d.dropna(subset = [strat])
      fig = px.histogram(d.dropna(subset = [strat]), facet_col = strat, x = bl, y = "participant_id", color = var, barnorm = "percent",
                        nbins = nb)
  if counts == "Count":
    if strat == "None":
      fig = px.histogram(d, x = bl, y = "participant_id", color = var, 
                          nbins = nb)
    else:
      d = d.dropna(subset = [strat])
      fig = px.histogram(d.dropna(subset = [strat]), facet_col = strat, x = bl, y = "participant_id", color = var,
                        nbins = nb)

  st.write(fig)
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the sankey diagrams
###----------------------------------------------------------------------------------------------'''
if viz == "Sankey Diagram":
  ny = st.slider("Select Number of Years Since Diagnosis", 1, 10, 5, step = 1)
  voi = st.selectbox("Select a Variable of Interest", ["HY"])
  d = d[(d.yearsDx >= 0) & (d.yearsDx <= ny)]

  ydx = d.loc[:,["participant_id","yearsDx","HY"]].dropna().drop_duplicates(keep="first", subset = ["participant_id","yearsDx"])
  ydx.loc[:,"HY"] = ydx.HY.round(0)
  ydx.loc[:,"yearsDx"] = ydx.yearsDx.round(0)
  ydx.loc[:,"HY3"] = ["0" if i <= 3 
                        else "1" for i in ydx.HY]

  SNK = ydx.pivot_table(values = "HY3", index = "participant_id", columns = "yearsDx", aggfunc = 'max')
  SNK.columns = range(0,ny+1)

  miss = st.selectbox("How are missing values handled?", ["Include Censoring", "Interpolate Values"])

  if miss == "Include Censoring":
    SNK = SNK.fillna(2)
    nl = 3
  if miss == "Interpolate Values":
    SNK = SNK.fillna(method = 'pad', axis = 1).fillna(method = 'bfill', axis = 1)
    nl = 2

  count = SNK.reset_index().pivot_table(values = "participant_id", columns = range(0,ny+1), aggfunc = "count")

  count = count.melt()

  MC2 = pd.DataFrame()
  for i in range(0,ny):
    MC = count.iloc[:,[i,(i+1),ny+1]].groupby([i,(i+1)]).agg(sum).reset_index()
    MC.columns = ["SRC","TRGT","VAL"]
    MC = MC.astype(int)
    MC.iloc[:,1] = [int(x)+(nl*i)+nl for x in list(MC.iloc[:,1])]
    if i > 0:
      MC.iloc[:,0] = [int(x)+(nl*i) for x in list(MC.iloc[:,0])]
    MC2 = MC2.append(MC)

  MC2 = MC2.sort_values(["SRC","TRGT"])

  if nl == 2:
    labl = ["HY<3", "HY>=3"]* int(len(MC2) / nl)
    node_y = [0,1]* int(len(MC2) / nl)
  elif nl == 3:
    labl = ["HY<3", "HY>=3", "CENS"]* int(len(MC2) / nl)
    node_y = [0,1,2]* int(len(MC2) / nl)

  fig = go.Figure(data=[go.Sankey(
      node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = labl,
        y = node_y,
        color = "blue"
      ),
      link = dict(
        source = np.array(MC2.SRC), # indices correspond to labels, eg A1, A2, A1, B1, ...
        target = np.array(MC2.TRGT),
        value = np.array(MC2.VAL)
    ))])
  
  fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
  st.write(fig)
  
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the Kaplan Meier Survival plots
###----------------------------------------------------------------------------------------------'''
if viz == "Kaplan-Meier Survival":

  var = st.sidebar.selectbox("Pick a variable", ["HY","UPDRS1","UPDRS2","UPDRS3","UPDRS4",
                                                "depress_test_score",'gds15_total_score',
                                                'moca_total_score','MoCA',
                                                'rbd_summary_score','ess_total_score','pdq39_mobility_score', 'pdq39_adl_score',
                                                'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
                                                'pdq39_cognition_score', 'pdq39_communication_score',
                                                'pdq39_discomfort_score'
                                                'bradykinesia', 'resting_tremor','rigidity'])
  

  

  if var in ['bradykinesia', 'resting_tremor','rigidity']:
    d.loc[:,"KM_var"] = d.loc[:,var].map({"Yes":1,"No":0,"Unknown":np.nan}).astype(int)
  else: 
    cutoff = st.sidebar.slider("Select a Cut-off Value", d.loc[:,var].min(), d.loc[:,var].median(), d.loc[:,var].max(), step = 1.0)
    if var == "MoCA":
      d.loc[:,"KM_var"] = [1 if i <= cutoff 
                      else 0 for i in d.loc[:,var]]
    else:
      d.loc[:,"KM_var"] = [1 if i >= cutoff 
                      else 0 for i in d.loc[:,var]]

  
  asdm = d.loc[:,"KM_var"].unique()
  d = d.dropna(subset = ["daysB",var])

  kmv = d.loc[:,["participant_id","KM_var"]].groupby("participant_id").agg(func="max").reset_index()

  kmv.columns = ["participant_id","KM_var_ev"]

  d = d.merge(kmv,how = "outer")

  d_0 = d[d.KM_var_ev == 0].drop_duplicates(subset = "participant_id",keep='last')
  d_1 = d[(d.KM_var_ev == 1) & (d.KM_var == 1)].drop_duplicates(subset = "participant_id",keep="first")

  dn = d_0.append(d_1).loc[:,["participant_id","daysB","KM_var"]]

  ## create a kmf object
  kmf = KaplanMeierFitter() 

  if st.sidebar.checkbox("Stratify Curve?"):
    ## Fit the data into the model
    kmf.fit(dn.loc[:,"daysB"], dn.loc[:,"KM_var"], label='Kaplan Meier Estimate')
    
    ## Create an estimate
    kmf.plot(ci_show=True, at_risk_counts = True)
  else:
    stvar = st.sidebar.selectbox("Pick a variable",['Phenotype',"study_arm",'sex','cohort','age_bin','ethnicity'])

    dn = dn.merge(d.loc[:,["participant_id",stvar]].drop_duplicates(subset="participant_id"), how = 'left') 
    items = dn.loc[:,stvar].unique()
    leni = len(items)

    ## fit the model for 1st cohort
    kmf.fit(dn.loc[dn.loc[:,stvar] == items[0],"daysB"], dn.loc[dn.loc[:,stvar] == items[0],"KM_var"], label=items[0])
    a1 = kmf.plot()

    for i in range(1,leni):
      ## fit the model for 1st cohort
      kmf.fit(dn.loc[dn.loc[:,stvar] == items[i],"daysB"], dn.loc[dn.loc[:,stvar] == items[i],"KM_var"], label=items[i])
      a1 = kmf.plot(ax = a1, at_risk_counts = True)


  if st.sidebar.checkbox("Stretch y-axis?",value=True):
    plt.ylim(0,1.05)

  st.set_option('deprecation.showPyplotGlobalUse', False)
  st.pyplot()
  
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the line plot
###----------------------------------------------------------------------------------------------'''
alt.data_transformers.enable('default', max_rows=None)

if viz == "Line Plot":

  var = st.sidebar.selectbox("Pick a Variable", ["HY","MoCA","UPDRS1","UPDRS2","UPDRS3","All UPDRS"])
  bl = st.sidebar.selectbox("Select a Baseline", ["Study Baseline","Diagnosis"])

  bl_name = "Baseline" if bl == "Study Baseline" else "Diagnosis"
  bl_val =  "yearsB" if bl == "Study Baseline" else "yearsDx"

  # if var == "HY":
  #     # Plot for HY
  
  if var == "All UPDRS":
          alt.data_transformers.enable('default', max_rows=None)

          source = pd.DataFrame({
              'Years Since '+ bl_name: d[bl_val],
              'UPDRS1': d["UPDRS1"],
              'UPDRS2': d["UPDRS2"],
              'UPDRS3': d["UPDRS3"],
          })


          base = alt.Chart(source).mark_circle(opacity=0, color='white').transform_fold(
              fold=['UPDRS1', 'UPDRS2', 'UPDRS3'],
              as_=['category', 'UPDRS Score']
          ).encode(
            alt.X(('Years Since '+ bl_name + ':Q'),
                    axis=alt.Axis(grid=False)),
              alt.Y('UPDRS Score:Q',
                    axis=alt.Axis(grid=False)),
              alt.Color('category:N')
          ).properties(width = 750, height = 500)
          # Creating 95 CI
          band = alt.Chart(source).mark_errorband(extent='ci').transform_fold(
            fold=['UPDRS1', 'UPDRS2', 'UPDRS3'],
              as_=['category', 'UPDRS Score']
          ).encode(
              alt.X(('Years Since '+ bl_name + ':Q'),
                    axis=alt.Axis(grid=False)),
              alt.Y('UPDRS Score:Q',
                    axis=alt.Axis(grid=False)),
              alt.Color('category:N')
          )

          base = base + band + base.transform_loess(('Years Since '+ bl_name), 'UPDRS Score', 
                                    groupby=['category']).mark_line(size=4)

          st.write(base)
  else:
    alt.data_transformers.enable('default', max_rows=None)

    source = pd.DataFrame({
                ('Years Since ' + bl_name): d[bl_val],
                var: d.loc[:,var]
            })

    base = alt.Chart(source).mark_circle(opacity=0, color='white'

            ).encode(
                alt.X(('Years Since '+ bl_name +':Q'),
                      axis=alt.Axis(grid=False)),
                alt.Y((var + ':Q'),
                      axis=alt.Axis(grid=False))
            ).properties(width = 750, height = 500)

            # Creating 95 CI
    band = alt.Chart(source).mark_errorband(extent='ci', color ='green').encode(
                alt.X(('Years Since '+ bl_name +':Q'),
                      axis=alt.Axis(grid=False)),
                alt.Y((var+':Q'),
                      axis=alt.Axis(grid=False)),   
            )

    base = base + band+ base.transform_loess(('Years Since'+ bl_name), var,
                                        groupby=['category']).mark_line(size=4, color="red")

    st.altair_chart(base)
