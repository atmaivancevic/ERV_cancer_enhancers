#!/usr/bin/env python
# coding: utf-8

import plotly.express as px
import pandas as pd
import os
import plotly.graph_objects as go

### patients from Orouji et al. 2021

df = pd.read_csv('Orouji_H3K27ac_giggle_results.tab', sep='\t')
#print(df)

fig = go.Figure(data=go.Heatmap(
                   z=df,
                   x=["COAD-003","COAD-004","COAD-005","COAD-007","COAD-008","COAD-011","COAD-012","COAD-013","COAD-015","COAD-016","COAD-017","COAD-018","FAPP-01","FAPP-05"],
                   y=["LTR10F_tumor", "LTR10F_normal", "LTR10A_tumor", "LTR10A_normal"],
                   colorscale=[(0, '#F4F5FF'),(0.1, 'white'), (1, "red")],
                    ))

for template in ["plotly_white"]:    
    #fig.update_xaxes(categoryorder="category ascending")
    #fig.update_yaxes(categoryorder="category ascending")
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=1300,height=500)
    #fig.update_xaxes(tickvals=["LTR10A", "LTR10F"])

fig.show()
fig.write_image("Orouji2021_patients_H3K27ac.svg")

##################################################################

### CEMT patients

df2 = pd.read_csv('CEMT_H3K27ac_giggle_results.tab', sep='\t')
#print(df2)

fig2 = go.Figure(data=go.Heatmap(
                   z=df2,
                   x=["AKCC70", "AKCC63", "AKCC58", "AKCC55","AKCC52"],
                   y=["LTR10F_tumor", "LTR10F_normal", "LTR10A_tumor", "LTR10A_normal"],
                   colorscale=[(0, 'white'), (1, "red")],
                    ))

for template in ["plotly_white"]:    
    fig2.update_layout(template=template)
    fig2.update_layout(autosize=False,width=600,height=500)

fig2.show()
fig2.write_image("cemt_patients_H3K27ac.svg")
