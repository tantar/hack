from google.colab import drive
import streamlit as st
import pandas as pd
import numpy as np
import datetime as dt
import os
import altair as alt

drive.mount('/content/gdrive')
datafolder = f'/content/gdrive/MyDrive/Hackathon/data'
workfolder = f'/content/gdrive/MyDrive/Hackathon/work'

os.chdir(workfolder)\


st.title("Hackathon Proof of Concept")
st.sidebar.header("Variable of Interest")
dataset = st.selectbox("Pick a Dataset", ["PPMI","PDBP","All"])
var = st.sidebar.selectbox("Pick a Variable", ["HY","UPDRS1","UPDRS2","UPDRS3","UPDRS4"])

cols = ["participant_id","visit_name", "Phenotype","hoehn_and_yahr_stage","mds_updrs_part_i_summary_score",
        "mds_updrs_part_ii_summary_score","mds_updrs_part_iii_summary_score","mds_updrs_part_iv_summary_score",
        "age_at_baseline","age_at_diagnosis","date_visit","date_baseline"]

if dataset == "PPMI":
  d = pd.read_csv("data/PPMI_processed.csv", usecols = cols)
if dataset == "PDBP":
  d = pd.read_csv("data/PDBP_processed.csv", usecols = cols)
if dataset == "All":
  d1 = pd.read_csv("data/PPMI_processed.csv", usecols = cols)
  d2 = pd.read_csv("data/PDBP_processed.csv", usecols = cols)
  d = d1.append(d2)

d = d.rename({"hoehn_and_yahr_stage" : "HY","mds_updrs_part_i_summary_score" : "UPDRS1",
        "mds_updrs_part_ii_summary_score":"UPDRS2","mds_updrs_part_iii_summary_score":"UPDRS3",
        "mds_updrs_part_iv_summary_score":"UPDRS4"}, axis = 1)
d = d[d.Phenotype == "PD"]

d.loc[:,"BL_durdx"] = d.age_at_baseline - d.age_at_diagnosis
d.loc[:,"daysB"] = d.date_visit - d.date_baseline
d.loc[:,"daysDx"] = (d.daysB + d.BL_durdx * 365.25).round(2)
d.loc[:,"yearsDx"] = (d.daysDx / 365.25).round(0)
d = d[(d.yearsDx > 0)]
yearsDx = "Years Since Diagnosis"

slice_width = st.slider("Time Interval", 1.0, 7.5, 2.5, step = 0.25)
d.loc[:,"yearsDx"] = (d.yearsDx / slice_width).round(0) * slice_width
d = d.rename({"yearsDx" : "Years Since Diagnosis"},axis=1)
size = slice_width * (7.5 - slice_width/10)

if var == "HY":
  ydx = d.loc[:,["participant_id",yearsDx,"HY"]].dropna().drop_duplicates(keep="first", subset = ["participant_id",yearsDx])
  ydx.loc[:,"HY"] = ydx.HY.round()
  ydx.loc[:,"HY Stage Greater Than 3"] = ["HY<3" if i < 3 
                        else "HY>=3" for i in ydx.HY]

  c2 = alt.Chart(ydx).mark_bar(size = size).encode(
     x="Years Since Diagnosis",
     y=alt.Y("count():Q", stack="normalize"),
     color='HY Stage Greater Than 3:N'
  ).properties(width=750)
elif var == "UPDRS1":
  ydx = d.loc[:,["participant_id",yearsDx,"UPDRS1"]].dropna().drop_duplicates(keep="first", subset = ["participant_id",yearsDx])
  ydx.loc[:,"UPDRS Part 1 Score"] = [(var + "<10")  if i < 10 
                      else (var + "<20") if i < 20
                      else (var + "<30")  if i < 30
                      else (var + "<40")  if i< 40
                      else (var + "<50")  if i< 50
                      else np.nan for i in ydx.loc[:,var]]

  c2 = alt.Chart(ydx).mark_bar(size = size).encode(
     x="Years Since Diagnosis",
     y=alt.Y("count():Q", stack="normalize"),
     color='UPDRS Part 1 Score:N'
  ).properties(width=750)
elif var == "UPDRS2":
  ydx = d.loc[:,["participant_id",yearsDx,"UPDRS2"]].dropna().drop_duplicates(keep="first", subset = ["participant_id",yearsDx])
  ydx.loc[:,"UPDRS Part 2 Score"] = [(var + "<10")  if i < 10 
                      else (var + "<20") if i < 20
                      else (var + "<30")  if i < 30
                      else (var + "<40")  if i< 40
                      else (var + "<50")  if i< 50
                      else np.nan for i in ydx.loc[:,var]]
  
  c2 = alt.Chart(ydx).mark_bar(size = size).encode(
     x="Years Since Diagnosis",
     y=alt.Y("count():Q", stack="normalize"),
     color='UPDRS Part 2 Score:N'
  ).properties(width=750)
elif var == "UPDRS3":
  ydx = d.loc[:,["participant_id",yearsDx,"UPDRS3"]].dropna().drop_duplicates(keep="first", subset = ["participant_id",yearsDx])
  ydx.loc[:,"UPDRS Part 3 Score"] = [(var + "<10")  if i < 10 
                      else (var + "<20") if i < 20
                      else (var + "<30")  if i < 30
                      else (var + "<40")  if i< 40
                      else (var + "<50")  if i< 50
                      else np.nan for i in ydx.loc[:,var]]

  c2 = alt.Chart(ydx).mark_bar(size = size).encode(
     x="Years Since Diagnosis",
     y=alt.Y("count():Q", stack="normalize"),
     color='UPDRS Part 3 Score:N'
  ).properties(width=750)
elif var == "UPDRS4":
  ydx = d.loc[:,["participant_id",yearsDx,"UPDRS4"]].dropna().drop_duplicates(keep="first", subset = ["participant_id",yearsDx])
  ydx.loc[:,"UPDRS Part 4 Score"] = [(var + "<10")  if i < 10 
                      else (var + "<20") if i < 20
                      else (var + "<30")  if i < 30
                      else (var + "<40")  if i< 40
                      else (var + "<50")  if i< 50
                      else np.nan for i in ydx.loc[:,var]]

  c2 = alt.Chart(ydx).mark_bar(size = size).encode(
     x="Years Since Diagnosis",
     y=alt.Y("count():Q", stack="normalize"),
     color='UPDRS Part 4 Score:N'
  ).properties(width=750)



if st.checkbox("Show Raw Data"):
  st.subheader("Raw data")
  st.write(d)

st.altair_chart(c2)