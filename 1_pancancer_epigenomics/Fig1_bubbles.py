#!/usr/bin/env python
# coding: utf-8

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os

# Read table
df = pd.read_csv('SuppTable1_TCGA_giggle_results_top23TEs.tab', sep='\t')
#print(df)

# Filter table to keep positive scores only (i.e. giggle_score>0)
columns = ['giggle_score']
filter_= (df[columns] > 0).all(axis=1)
positive_df = df[filter_]
#print(positive_df)

# Make a bubble plot
fig = px.scatter(positive_df, y="repeat_family", x="cancer_type", size="giggle_score",  color="odds_ratio", size_max=18, opacity=1,
                 color_continuous_scale=[(0.0, '#FD9367'), (0.33, '#C3305D'), (0.67, '#782D65'), (1, '#432967')],
                 labels={"cancer_type":"Cancer type", "giggle_score":"Enrichment", "repeat":"Repeats"},
                 category_orders={"cancer_type":["LUAD", "COAD", "CHOL", "KIRP", "SKCM", "STAD", "CESC", 
                                            "BRCA", "ESCA", "KIRC", "BLCA", "PRAD", "LUSC", "THCA", 
                                            "LIHC", "HNSC", "GBM", "LGG", "ACC", "MESO", "PCPG"],
                                  "repeat_family": ["LTR10A", "LTR10F", "LTR5_Hs", "MER51A", "LTR2B", "MER11B", 
                                             "SVA_A", "MER11D", "MER11A", "LTR5B", "LTR10C", "PABL_A", 
                                             "Harlequin-int", "LTR9", "LTR7Y", "HERVE_a-int", "LTR13", 
                                             "MER44B", "MER41B", "LTR13A", "LTR3B", "MER39B", "LTR3A"]}
                )

# Update plot layout
for template in ["plotly_white"]:    
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=720,height=770)
    
fig.update_layout(font=dict(size=10))
fig.update_yaxes(title_standoff=12)
fig.update_xaxes(title_standoff=12)
fig.update_yaxes(tickfont_size=12, ticks="outside", tickangle=0, ticklen=4)
fig.update_xaxes(tickfont_size=12, tickangle=45, ticks="outside", ticklen=4)
fig.layout.font.family = 'Helvetica'

# Show and save fig
fig.show()
fig.write_image("bubbles.svg")
