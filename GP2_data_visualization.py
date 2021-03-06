# Loads libraries required for the app
from warnings import simplefilter
simplefilter(action='ignore', category=DeprecationWarning)

import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
import plotly.graph_objects as go
import plotly.express as px
import re
import glob
import tableone as tab1
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
# pd.set_option('display.max_rows', 1000)


@st.cache
def MeanSD(series) :
    mean = series.mean().round(1)
    SD = series.std().round(1)
    out = str(mean) + "(" + str(SD) + ")"
    return out
@st.cache
def MedQ13(series) :
    median = series.median().round(1)
    SD = series.quantile([.25,.75])
    out = str(median) + "[" + str(series.quantile(.25).round(1)) + "," + str(series.quantile(.75).round(1))  + "]"
    return out
@st.cache
def LoadReqs():
    # Creates list of file names, datasets, and dataset names
    files = glob.glob(f'*.csv')
    d_sets = glob.glob(f'data/*.csv')
    d_names = [re.sub("[procesdvat.\\\\/]", "", i) for i in d_sets]

    # The Reference sheet contains info on all possible columns in GP2, facilitating easier categorization for any cohort
    ref = pd.read_csv("Reference.csv")

    # If study_arms and info_table need to be produced from scratch, the data files are initialized here
    study_arms = list()
    info = pd.DataFrame()

    if any(i not in files for i in ["study_arms.csv", "info_table.csv"]):
        st.warning("Missing Files Detected, Rewriting Study Arms and Info Table")
        for dn in d_sets:
            cohort_name = re.sub("[_procesdvat.\\/]", "", dn)
            cohort = pd.read_csv(dn)
            cohort = cohort.astype(str)

            # LS1 is unique in that it's study_arm refers to an experimental arm rather than phenotypic one
            if "LS1" in cohort_name:
                cohort.loc[:, "study_arm"] = [cohort_name + "_" + i for i in cohort.Phenotype]
            else:
                cohort.loc[:, "study_arm"] = [cohort_name + "_" + i for i in cohort.study_arm]

            study_arms = study_arms + list(cohort.study_arm.unique())

            study_arms.sort()
            tinfo2 = pd.DataFrame()

            for j in cohort.study_arm.unique():
                tinfo = pd.DataFrame(data={'Study Arm': j,
                                           'N': len(cohort.loc[(cohort.study_arm == j), "participant_id"].unique()),
                                           'N_OBS': cohort.loc[(cohort.study_arm == j), "participant_id"].count(),
                                           'Mean BL Age': MeanSD(
                                               cohort.loc[(cohort.study_arm == j), "age_at_baseline"]),
                                           'Median BL DurDx': MedQ13(pd.Series([i - j for i, j in
                                                                                zip(cohort.age_at_baseline,
                                                                                    cohort.age_at_diagnosis)])),
                                           'FU Time': (cohort.drop_duplicates(subset="participant_id", keep="last").loc[
                                                           (cohort.study_arm == j), "visit_month"].mean() / 12),
                                           'Visit Date': ["X" if "date_visit" in cohort.columns else " "],
                                           'Family History': ["X" if "family_hx_1st" in cohort.columns else " "],
                                           'Vital': ["X" if all(elem in cohort.columns for elem in list(
                                               ref.loc[ref.Modality == "Vital", "Item"])) else " "],
                                           'UPDRS': [
                                               "X" if "mds_updrs_part_i_summary_score" in cohort.columns else " "],
                                           'MoCA': ["X" if "moca_total_score" in cohort.columns else " "],
                                           'SCOPA-COG': ["X" if "scopa_cog_total_score" in cohort.columns else " "],
                                           'RBDSQ': ["X" if "rbd_summary_score" in cohort.columns else " "],
                                           'ESS': ["X" if "ess_total_score" in cohort.columns else " "],
                                           'PDQ39': ["X" if "pdq39_discomfort_score" in cohort.columns else " "],
                                           'Depression': ["X" if any(i in cohort.columns for i in ["gds15_total_score",
                                                                                                   "depress_test_score"]) else " "],
                                           'Olfactory': ["X" if "smell_test_results" in cohort.columns else " "],
                                           },

                                     index=[0])
                tinfo2 = tinfo2.append(tinfo)

            info = info.append(tinfo2)

        pd.Series(study_arms).to_csv("study_arms.csv", index=False)
        info.to_csv("info_table.csv", index=False)
    else:
        study_arms = list(pd.read_csv("study_arms.csv").iloc[:, 1])
        info = pd.read_csv("info_table.csv")

    return study_arms, info, ref, d_sets
@st.cache
def LoadData():
    d = pd.DataFrame()
    for dn in d_sets:
        cohort_name = re.sub("[_procesdvat.\\\\/]", "", dn)
        cohort = pd.read_csv(dn,low_memory=False) #usecols = ['participant_id',"visit_month","study_arm","Phenotype",'primary_diagnosis',"age_at_baseline","age_at_diagnosis"])

        if "LS1" in cohort_name:
            cohort.loc[:,"study_arm"] = [cohort_name + "_" + i for i in cohort.Phenotype]
        else:
            cohort.loc[:,"study_arm"] = [cohort_name + "_" + i for i in cohort.study_arm]


        if (mode == "Advanced") & (CST == "Preset Cohorts"):
            cohort.loc[:,"cohort"] = cohort_name
        else:
            if COI == "":
                cohort.loc[:,"cohort"] = "COI"
            else:
                cohort.loc[:,"cohort"] = cohort_name

            if OD == "All Others":
                d2 = dat[~dat.study_arm.isin(select_arms)]
                d2["cohort"] = OD
                d1 = d1.append(d2)
            elif OD == "All Other PD":
                d2 = dat[(~dat.study_arm.isin(select_arms)) & (dat.study_arm.str.contains("_PD"))]
                d2["cohort"] = OD
                d1 = d1.append(d2)
            elif OD == "All Other HC":
                d2 = dat[(~dat.study_arm.isin(select_arms)) & (dat.study_arm.str.contains("_HC"))]
                d2["cohort"] = OD
                d1 = d1.append(d2)

        d = pd.concat([d,cohort], ignore_index=True, axis=0)
    return d


# Set up the headers for the app
st.title("GP2 Data Visualization Tool")
#mode = st.sidebar.radio("Mode",["Basic","Advanced"])
mode = "Basic"
st.sidebar.title("Options")

study_arms, info, ref, d_sets = LoadReqs()

#In Advanced mode, users are given the ability to create and name their own cohort, then compare it to other cohorts
if mode == "Advanced":
    CST = st.selectbox("Cohort Selection Type", ["Preset Cohorts", "Choose-Your-Own Cohorts"])
#In basic mode, users do not select their own cohorts and cannot compare.
else:
    COI = "cohort"
    CST = "Preset Cohorts"
    OD = "No Comparison"

#Loads all possible columns from reference sheet
cols = list(ref.Item)

d = LoadData()


#In Basic Mode, the first condition will always be met.
# In Advanced Mode, the user is given the option to use preset cohorts or create their own
if CST == "Preset Cohorts":
        if st.checkbox("All PD"):
            select_data = st.multiselect("Select Datasets", study_arms,
                                         [x for x in study_arms if "_PD" in x])
        else:
            select_data = st.multiselect("Select Datasets", study_arms,
                                         ["PPMI_PD","BIOFIND_PD","LS1_PD"])
else:
        select_data = st.multiselect("What study arms are in your cohort of interest (COI)?", study_arms,
                                     ["PPMI_PD","BIOFIND_PD","LS1_PD"])

        COI = st.text_input("Cohort Name")
        OD = st.selectbox("What will this cohort be compared to?",
                          ["No Comparison", "All Others", "All Other PD", "All Other HC"])


# select_arms = [x.replace(re.sub("[procesdvat.\\\\]", "", dn), '') for x in select_data if
#                re.sub("[procesdvat.\\\\]", "", dn) in x]

# if cohort_name not in ["LS1", "PROBAND"]:
#     d1 = dat[dat.study_arm.isin(select_arms)]
# else:
#     d1 = dat.copy()

d = d[d.study_arm.isin(select_data)]

#If no cohorts have been selected, the apps halt
if not select_data:
    st.warning("Please Select a Cohort")
    st.stop()

###----------------------------------------------------------------------------------------------
###This section of code performs data cleaning for output in streamlit
###----------------------------------------------------------------------------------------------



# #Only loads d_sets already processed
# d_set_target = pd.Series(i[0] for i in pd.Series(select_data).str.split("_")).unique()
#
# d_sets2 = list()
# for j in d_set_target:
#     dst = [d for d in d_sets if j in d]
#     d_sets2 = d_sets2 + dst
#
#
# d1 = d1.loc[:, d1.columns.isin(cols)]
#
# t1 = d1.copy()

cats = list(ref.loc[ref.ItemType == "string", "Item"])

if mode == "Basic":
    # Creating checkboxes and slide bars for dataset, diagram choice, and variable of choice
    viz = st.sidebar.selectbox("Choose a Visualization", ["Tables", "Scatter Plot","Bar Graph","Baseline Histogram", "Kaplan-Meier Survival"])#,"Cohort Compare"
if mode == "Advanced":
    # Creating checkboxes and slide bars for dataset, diagram choice, and variable of choice
    viz = st.sidebar.selectbox("Choose a Visualization", ["Tables", "Scatter Plot","Line Plot",
                                                          "Bar Graph","Baseline Histogram",
                                                          "Sankey Diagram", "Kaplan-Meier Survival"])


d = d.rename({"mds_updrs_part_i_summary_score": "UPDRS1",
              "mds_updrs_part_ii_summary_score": "UPDRS2",
              "mds_updrs_part_iii_summary_score": "UPDRS3",
              "mds_updrs_part_iv_summary_score": "UPDRS4",
              "moca_total_score": "MoCA",
              "hoehn_and_yahr_stage":"HY"}, axis=1)

# # This code was used to subset PDBP patients based on those with genetic data
# PDBPs = pd.read_csv("PDBP_sample.csv")
# PDBPs = PDBPs[PDBPs.Original_genotypeID.notna()]
# d = d[d.participant_id.isin(PDBPs.Original_clinicalID)]

# Essential Time Variables are Created
if ("date_visit" in d.columns):
    d.loc[:, "daysB"] = (d.date_visit - d.date_baseline).fillna(d.visit_month * 30)
    d.loc[:, "yearsB"] = (d.daysB / 365.25).round(1)
else:
    d.loc[:,"daysB"] = d.visit_month * 30
    d.loc[:, "yearsB"] = (d.daysB / 365.25).round(1)

if ("age_at_baseline" in d.columns) & ("age_at_diagnosis" in d.columns):
    d.loc[:, "BL_durdx"] = d.age_at_baseline - d.age_at_diagnosis
    d.loc[:, "daysDx"] = (d.daysB + d.BL_durdx * 365.25).round(2)
    d.loc[:, "yearsDx"] = (d.daysDx / 365.25).round(1)

yearsDx = "Years Since Diagnosis"
#
# # Creating new variables
# if "HY" in d.columns:
#     d.loc[:, "HY3"] = ["HY<3" if i < 3
#                        else "HY>=3" for i in d.HY]
#
# if "age_at_diagnosis" in d.columns:
#     d.loc[:, "age_bin"] = ["<50" if i < 50
#                                  else ">50" for i in d.age_at_diagnosis]
#
# flx_cols = ["code_upd2403_time_spent_in_the_off_state","code_upd2404_functional_impact_of_fluctuations",
#               "code_upd2405_complexity_of_motor_fluctuations"]
#
# if all([c in d.columns for c in flx_cols]):
#     d.loc[:, 'motor_fluctuations'] = ["Yes" if i > 0 else "No" for i in d.loc[:, flx_cols].sum(axis=1, skipna=True)]
#
# if "code_upd2314_body_bradykinesia" in d.columns:
#     d.loc[:, 'bradykinesia'] = ["Yes" if i > 0 else "No" for i in d.code_upd2314_body_bradykinesia]
#
# rgd_cols = ["code_upd2303a_rigidity_neck", "code_upd2303b_rigidity_rt_upper_extremity",
#             "code_upd2303c_rigidity_left_upper_extremity", "code_upd2303d_rigidity_rt_lower_extremity",
#             "code_upd2303e_rigidity_left_lower_extremity"]
#
# if all([c in d.columns for c in rgd_cols]):
#     d.loc[:, 'rigidity'] = ["Yes" if i > 0 else "No" for i in d.loc[:, rgd_cols].sum(axis=1, skipna=True)]
#
# trem_cols = ["code_upd2317a_rest_tremor_amplitude_right_upper_extremity","code_upd2317b_rest_tremor_amplitude_left_upper_extremity",
#                                        "code_upd2317c_rest_tremor_amplitude_right_lower_extremity","code_upd2317d_rest_tremor_amplitude_left_lower_extremity",
#                                        "code_upd2317e_rest_tremor_amplitude_lip_or_jaw", "code_upd2318_consistency_of_rest_tremor"]
#
# if all([c in d.columns for c in trem_cols]):
#     d.loc[:, 'resting_tremor'] = ["Yes" if i > 0 else "No" for i in d.loc[:, trem_cols].sum(axis=1, skipna=True)]
#
cats = cats + ["bradykinesia","rigidity","resting_tremor","hx_dementia"]

nnml = d.columns[~d.columns.isin(cats)]

if st.checkbox("Show Study Arm Info Table"):
    st.write(info)
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the scatterplots
###----------------------------------------------------------------------------------------------'''
if viz == "Scatter Plot":
    xax = st.sidebar.selectbox("Select an x-axis variable", d.columns[d.columns.isin(nnml)])
    yax = st.sidebar.selectbox("Select a y-axis variable", d.columns[d.columns.isin(nnml)])

    if st.sidebar.checkbox("Include only Baseline?"):
        d = d[d.visit_month == 0]

    d = d.dropna(subset=[xax, yax])

    tl = None
    if st.sidebar.checkbox("Add Trendline?"):
        tl = st.sidebar.selectbox("Select Trendline Type", ["ols", "lowess"])

    if st.sidebar.checkbox("Stratify Based on Third Variable?"):
        cs = st.sidebar.selectbox("Select a stratifying variable", d.columns[d.columns.isin(cats)])
        d = d.dropna(subset=[cs])
        fig = px.scatter(d, x=xax, y=yax, color=cs, opacity=0.5, trendline=tl)

    else:
        fig = px.scatter(d, x=xax, y=yax, opacity=0.5, trendline=tl)

    st.write(fig)

###----------------------------------------------------------------------------------------------###
###The following code is responsible for producing the demographic summarization table
###----------------------------------------------------------------------------------------------'''
if viz == "Tables":

    foc = st.sidebar.selectbox("Select a Focus", ["Overall", "Demographic", "Non-Motor", "Choose-Your-Own"])
    gby = st.sidebar.selectbox("Stratifying Variable", ["None", "study_arm", "cohort", 'Phenotype'])

    if foc == "Overall":
        # columns to summarize
        columns = ['sex', 'Phenotype', 'race', 'education_level', 'primary_diagnosis', 'age_at_baseline',
                   'age_at_diagnosis', 'education_years', 'BL_durdx', 'UPDRS1', 'UPDRS2', 'UPDRS3', 'UPDRS4']

        # columns containing categorical variables
        columns = [i for i in columns if i in d.columns]

        # columns containing categorical variables
        categorical = [i for i in cats if i in columns]

        # non-normal variables
        nonnormal = [i for i in nnml if i in columns]

    if foc == "Demographic":
        # columns to summarize
        columns = ['sex', 'Phenotype', 'race', 'education_level', 'primary_diagnosis', 'age_at_baseline',
                   'age_at_diagnosis', 'education_years', 'BL_durdx']

        # columns containing categorical variables
        categorical = [i for i in cats if i in columns]

        # non-normal variables
        nonnormal = [i for i in nnml if i in columns]

    if foc == "Motor":
        # columns to summarize
        columns = ['BL_durdx', 'UPDRS3', 'UPDRS4', 'bradykinesia', 'resting_tremor',
                   'rigidity']

        # columns containing categorical variables
        categorical = [i for i in cats if i in columns]

        # non-normal variables
        nonnormal = [i for i in nnml if i in columns]


    if foc == "Non-Motor":
        columns = ["gds15_total_score", 'depress_test_name', 'depress_test_score', "MoCA",
                   'clinical_dx_depression', 'clinical_dx_insomnia',
                   'pdq39_mobility_score', 'pdq39_adl_score',
                   'smell_test_name', 'smell_test_score',
                   'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
                   'pdq39_cognition_score', 'pdq39_communication_score',
                   'pdq39_discomfort_score', ]

        # columns containing categorical variables
        categorical = [i for i in cats if i in columns]

        # non-normal variables
        nonnormal = [i for i in nnml if i in columns]

    if foc == "Choose-Your-Own":
        if st.checkbox("Include All?"):
            columns = [i for i in cols if i in d.columns[2:]]
        else:
            columns = st.multiselect("Select Variables", list(d.columns[2:]), ['sex', 'Phenotype', 'age_at_baseline'])

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
            order[col] = ["Yes", "No", "Unknown"]

    # alternative labels
    # labels={'death': 'Mortality'}

    # set decimal places for age to 0
    # decimals = {"Age": 0}

    # optionally, a categorical variable for stratification
    # groupby = ['visit_month']

    # rename the death column
    # labels={'death': 'Mortality'}

    # display minimum and maximum for listed variables
    # min_max = ['Height']

    if not columns:
        st.warning("No Variable Has Been Selected")
    elif gby == "None":
        dt = tab1.TableOne(d.drop_duplicates(subset="participant_id", keep="first"), columns=columns,
                           categorical=categorical,
                           nonnormal=nonnormal, order=order)
        st.write(dt)
    else:
        dt = tab1.TableOne(d.drop_duplicates(subset="participant_id", keep="first"), columns=columns,
                           categorical=categorical,
                           nonnormal=nonnormal, groupby=gby, order=order)
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

    strat = st.sidebar.selectbox("Pick a Stratifier", ["None"] + list(d.columns[7:]))

    if var in nnml:
        nb = st.sidebar.slider("Number of Bins", 1, 25, 50, step=1)
    else:
        nb = len(d.loc[:,var].unique())

    if strat == "None":
        fig = px.histogram(d, x=var, nbins=nb)

        if counts == True:
            fig = px.histogram(d, x=var, barnorm="percent",
                               nbins=nb)

    elif counts == True:
        d = d.dropna(subset=[var, strat])
        fig = px.histogram(d, x=var, barnorm="percent", color=strat, nbins=nb)
    else:
        d = d.dropna(subset=[strat])
        fig = px.histogram(d, x=var, color=strat, nbins=nb)

    st.write(fig)
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the new stacked bar graphs
###----------------------------------------------------------------------------------------------'''
if viz == "Bar Graph":
    var = st.sidebar.selectbox("Pick a y-axis", cats[3:])

    bl = st.sidebar.selectbox("Pick an x-axis", ["yearsDx", "yearsB", "visit_month",
                                                 "HY3", "UPDRS1", "UPDRS2", "UPDRS3",
                                                 "UPDRS4"])

    strat = st.sidebar.selectbox("Pick a Stratifier", ["None", "sex", "age_bin", "study_arm", "cohort"])
    counts = st.sidebar.selectbox("Pick a Histogram Type", ["Count", "Percent"])
    nb = st.sidebar.slider("Number of Bins", 1, 25, 50, step=1)
    d = d.dropna(subset=[var])

    if bl in ["yearsDx", "yearsB"]:
        d = d.drop_duplicates(subset=["participant_id", var])

    if var not in d.columns:
        st.write("This variable is not supported in the selected dataset")

    if counts == "Percent":
        if strat == "None":
            fig = px.histogram(d, x=bl, color=var, barnorm="percent",
                               nbins=nb)
        else:
            d = d.dropna(subset=[strat])
            fig = px.histogram(d.dropna(subset=[strat]), facet_col=strat, x=bl, color=var,
                               barnorm="percent",
                               nbins=nb)
    if counts == "Count":
        if strat == "None":
            fig = px.histogram(d, x=bl, color=var,
                               nbins=nb)
        else:
            d = d.dropna(subset=[strat])
            fig = px.histogram(d.dropna(subset=[strat]), facet_col=strat, x=bl, color=var,
                               nbins=nb)

    st.write(fig)
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the sankey diagrams
###----------------------------------------------------------------------------------------------'''
if viz == "Sankey Diagram":
    ny = st.slider("Select Number of Years Since Diagnosis", 1, 10, 5, step=1)
    voi = st.selectbox("Select a Variable of Interest", ["HY","bradykinesia","rigidity","resting_tremor","motor_fluctuations"])
    d = d[(d.yearsDx >= 0) & (d.yearsDx <= ny)]

    ydx = d.loc[:, ["participant_id", "yearsDx", voi]].dropna().drop_duplicates(keep="first",
                                                                                 subset=["participant_id", "yearsDx"])
#    ydx.loc[:, "HY"] = ydx.HY.round(0)
    ydx.loc[:, "yearsDx"] = ydx.yearsDx.round(0)
    if voi == "HY":
        ydx.loc[:, "VOI"] = ["0" if i <= 3
                             else "1" for i in ydx.loc[:,voi]]
    else:
        ydx.loc[:,"VOI"] = ydx.loc[:,voi].replace({"Yes":1,"No":0})

    SNK = ydx.pivot_table(values="VOI", index="participant_id", columns="yearsDx", aggfunc='max')
    SNK.columns = range(0, ny + 1)

    miss = st.selectbox("How are missing values handled?", ["Include Censoring", "Interpolate Values"])

    if miss == "Include Censoring":
        SNK = SNK.fillna(2)
        nl = 3
    if miss == "Interpolate Values":
        SNK = SNK.fillna(method='pad', axis=1).fillna(method='bfill', axis=1)
        nl = 2

    count = SNK.reset_index().pivot_table(values="participant_id", columns=range(0, ny + 1), aggfunc="count")

    count = count.melt()

    MC2 = pd.DataFrame()
    for i in range(0, ny):
        MC = count.iloc[:, [i, (i + 1), ny + 1]].groupby([i, (i + 1)]).agg(sum).reset_index()
        MC.columns = ["SRC", "TRGT", "VAL"]
        MC = MC.astype(int)
        MC.iloc[:, 1] = [int(x) + (nl * i) + nl for x in list(MC.iloc[:, 1])]
        if i > 0:
            MC.iloc[:, 0] = [int(x) + (nl * i) for x in list(MC.iloc[:, 0])]
        MC2 = MC2.append(MC)

    MC2 = MC2.sort_values(["SRC", "TRGT"])

    ll = ["Has"]
    if nl == 2:
        labl = ["True", "False"] * int(len(MC2) / nl)
        node_y = [0, 1] * int(len(MC2) / nl)
    elif nl == 3:
        labl = ["True", "False", "Censor"] * int(len(MC2) / nl)
        node_y = [0, 1, 2] * int(len(MC2) / nl)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labl,
            y=node_y,
            color="blue"
        ),
        link=dict(
            source=np.array(MC2.SRC),  # indices correspond to labels, eg A1, A2, A1, B1, ...
            target=np.array(MC2.TRGT),
            value=np.array(MC2.VAL)
        ))])

    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    st.write(fig)

###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the Kaplan Meier Survival plots
###----------------------------------------------------------------------------------------------'''
if viz == "Kaplan-Meier Survival":
    var = st.sidebar.selectbox("Pick a variable", ["HY", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4",
                                                   "living_status",
                                                   "depress_test_score", 'gds15_total_score',
                                                   'moca_total_score', 'MoCA', 'ledd_daily',
                                                   'rbd_summary_score', 'ess_total_score', 'pdq39_mobility_score',
                                                   'pdq39_adl_score',
                                                   'pdq39_emotional_score', 'pdq39_stigma_score', 'pdq39_social_score',
                                                   'pdq39_cognition_score', 'pdq39_communication_score',
                                                   'pdq39_discomfort_score',
                                                   "code_upd2403_time_spent_in_the_off_state",
                                                   'bradykinesia', 'resting_tremor', 'rigidity'])

    cohort_w_var = d.loc[d.loc[:,var].notna(),"study_arm"].unique()
    if len(d.study_arm.unique()) != len(cohort_w_var):
        missing_cohorts = [c for c in d.study_arm.unique() if c not in cohort_w_var]
        st.warning(f"{missing_cohorts}: Missing Variable of Interest")

    if var in ['bradykinesia', 'resting_tremor', 'rigidity']:
        d.loc[:, "KM_var"] = d.loc[:, var].map({"Yes": 1, "No": 0, "Unknown": np.nan}).astype(int)
    elif "liv" in var:
         var = f'{var}_at_censoring'
         d.loc[:, "KM_var"] = d.loc[:,var].mask(d.participant_id.duplicated(keep='last'),0)

    else:
        cutoff = st.sidebar.slider("Select a Cut-off Value", float(d.loc[:, var].min()), float(d.loc[:, var].median()),
                                   float(d.loc[:, var].max()), step=1.0)
        if var == "MoCA":
            d.loc[:, "KM_var"] = [1 if i <= cutoff
                                  else 0 for i in d.loc[:, var]]
        else:
            d.loc[:, "KM_var"] = [1 if i >= cutoff
                                  else 0 for i in d.loc[:, var]]

    d = d.dropna(subset=["daysB", var])



    kmv = d.loc[:, ["participant_id", "KM_var"]].groupby("participant_id").agg(func="max").reset_index()

    kmv.columns = ["participant_id", "KM_var_ev"]

    d = d.merge(kmv, how="outer")

    d_0 = d[d.KM_var_ev == 0].drop_duplicates(subset="participant_id", keep='last')
    d_1 = d[(d.KM_var_ev == 1) & (d.KM_var == 1)].drop_duplicates(subset="participant_id", keep="first")

    dn = d_0.append(d_1).loc[:, ["participant_id", "daysB", "KM_var"]]

    ## create a kmf object
    kmf = KaplanMeierFitter()


    if st.sidebar.checkbox("Stratify Curve?"):
        stvar = st.sidebar.selectbox("Pick a variable",
                                     ['cohort', 'Phenotype', "study_arm", 'sex', 'age_bin', 'ethnicity'])

        dn = dn.merge(d.loc[:, ["participant_id", stvar]].drop_duplicates(subset="participant_id"), how='left')
        items = dn.loc[:, stvar].unique()
        leni = len(items)

        ## fit the model for 1st cohort
        kmf.fit(dn.loc[dn.loc[:, stvar] == items[0], "daysB"], dn.loc[dn.loc[:, stvar] == items[0], "KM_var"],
                label=items[0])
        a1 = kmf.plot_survival_function()

        for i in range(1, leni):
            ## fit the model for 1st cohort
            kmf.fit(dn.loc[dn.loc[:, stvar] == items[i], "daysB"], dn.loc[dn.loc[:, stvar] == items[i], "KM_var"],
                    label=items[i])
            a1 = kmf.plot_survival_function(ax=a1, at_risk_counts=True)
    else:
        ## Fit the data into the model
        kmf.fit(dn.loc[:, "daysB"], dn.loc[:, "KM_var"], label='Kaplan Meier Estimate')

        ## Create an estimate
        kmf.plot_survival_function(ci_show=True, at_risk_counts=True)

    if st.sidebar.checkbox("Stretch y-axis?", value=True):
        plt.ylim(0, 1.05)

    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the cohort compare function
###----------------------------------------------------------------------------------------------'''

if viz == "Cohort Compare":
    type = st.selectbox("Select Type of Subset to Compare", ["study_arm", "cohort"])
    col1, col2 = st.beta_columns(2)

    with col1:
        co1 = st.selectbox("Select Cohort 1",d.loc[:,type].unique())
    with col2:
        co2 = st.selectbox("Select Cohort 2",d.loc[:,type].unique())
###----------------------------------------------------------------------------------------------###
###The following code is reponsible for producing the line plot
###----------------------------------------------------------------------------------------------'''

alt.data_transformers.enable('default', max_rows=None)

if viz == "Line Plot":

    var = st.sidebar.selectbox("Pick a Variable", ["HY", "MoCA", "UPDRS1", "UPDRS2", "UPDRS3", "All UPDRS"])
    bl = st.sidebar.selectbox("Select a Baseline", ["Study Baseline", "Diagnosis"])

    bl_name = "Baseline" if bl == "Study Baseline" else "Diagnosis"
    bl_val = "yearsB" if bl == "Study Baseline" else "yearsDx"

    # if var == "HY":
    #     # Plot for HY

    if var == "All UPDRS":
        alt.data_transformers.enable('default', max_rows=None)

        source = pd.DataFrame({
            'Years Since ' + bl_name: d[bl_val],
            'UPDRS1': d["UPDRS1"],
            'UPDRS2': d["UPDRS2"],
            'UPDRS3': d["UPDRS3"],
        })

        base = alt.Chart(source).mark_circle(opacity=0, color='white').transform_fold(
            fold=['UPDRS1', 'UPDRS2', 'UPDRS3'],
            as_=['category', 'UPDRS Score']
        ).encode(
            alt.X(('Years Since ' + bl_name + ':Q'),
                  axis=alt.Axis(grid=False)),
            alt.Y('UPDRS Score:Q',
                  axis=alt.Axis(grid=False)),
            alt.Color('category:N')
        ).properties(width=750, height=500)
        # Creating 95 CI
        band = alt.Chart(source).mark_errorband(extent='ci').transform_fold(
            fold=['UPDRS1', 'UPDRS2', 'UPDRS3'],
            as_=['category', 'UPDRS Score']
        ).encode(
            alt.X(('Years Since ' + bl_name + ':Q'),
                  axis=alt.Axis(grid=False)),
            alt.Y('UPDRS Score:Q',
                  axis=alt.Axis(grid=False)),
            alt.Color('category:N')
        )

        base = base + band + base.transform_loess(('Years Since ' + bl_name), 'UPDRS Score',
                                                  groupby=['category']).mark_line(size=4)

        st.altair_chart(base)
    else:
        alt.data_transformers.enable('default', max_rows=None)

        source = pd.DataFrame({
            ('Years Since ' + bl_name): d[bl_val],
            var: d.loc[:, var]
        })

        base = alt.Chart(source).mark_circle(opacity=0, color='white'

                                             ).encode(
            alt.X(('Years Since ' + bl_name + ':Q'),
                  axis=alt.Axis(grid=False)),
            alt.Y((var + ':Q'),
                  axis=alt.Axis(grid=False))).properties(width = 750, height = 500)

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
