#!/usr/bin/env python
# coding: utf-8

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os

df = pd.read_csv('tcgaMinusRoadmap_vs_TEs_23topTEs_giggle0.tab', sep='\t')
#print(df)

fig = px.scatter(df, y="repeat", x="cancer", size="gigglescore",  color="oddsratio", size_max=18, opacity=1,
                 color_continuous_scale=[(0.0, '#FD9367'), (0.33, '#C3305D'), (0.67, '#782D65'), (1, '#432967')],
                 labels={"cancer":"Cancer type", "gigglescore":"Enrichment", "repeat":"Repeats"},
                 category_orders={"cancer":["LUAD", "COAD", "CHOL", "KIRP", "SKCM", "STAD", "CESC", 
                                            "BRCA", "ESCA", "KIRC", "BLCA", "PRAD", "LUSC", "THCA", 
                                            "LIHC", "HNSC", "GBM", "LGG", "ACC", "MESO", "PCPG"],
                                  "repeat": ["LTR10A", "LTR10F", "LTR5_Hs", "MER51A", "LTR2B", "MER11B", 
                                             "SVA_A", "MER11D", "MER11A", "LTR5B", "LTR10C", "PABL_A", 
                                             "Harlequin-int", "LTR9", "LTR7Y", "HERVE_a-int", "LTR13", 
                                             "MER44B", "MER41B", "LTR13A", "LTR3B", "MER39B", "LTR3A"]}
                )

for template in ["plotly_white"]:    
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=720,height=770)
    
fig.update_layout(
    font=dict(
        size=10,
    )
)

fig.update_yaxes(title_standoff=12)
fig.update_xaxes(title_standoff=12)
fig.update_yaxes(tickfont_size=12, ticks="outside", tickangle=0, ticklen=4)
fig.update_xaxes(tickfont_size=12, tickangle=45, ticks="outside", ticklen=4)

fig.layout.font.family = 'Helvetica'

fig.show()
fig.write_image("bubbles.svg")
