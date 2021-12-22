import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os
from IPython.display import Image

df = pd.read_csv('Top_Cistrome_Repressors.tab', sep='\t')
#print(df)

fig = px.bar(df, 
                 x="combo_score", 
                 y="tf", 
                 #color="tf_role",
                 #text="name",
                 #color="cell_line",
                 #color_discrete_sequence=px.colors.sequential.Viridis,
                 #hover_data=df.columns,
                 #size="score", size_max=13,
                 #title="Top transcription factors for CRC-enriched LTR10A and LTR10F elements",
                 labels={"combo_score":"LTR10 Enrichment", "tf":"Transcription Factor", "cell":"Cell Line"}
                )

for template in ["plotly_white"]:    
    #fig.update_xaxes(categoryorder="category ascending")
    #fig.update_yaxes(categoryorder="category descending")
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=800,height=1000)
    
fig.update_yaxes(categoryorder="max ascending")

fig.add_trace(go.Scatter(
    x=[2240.834133169611+100,1422.2355803243328+100,1398.8030035511867+100,998.12299767979575+100,955.3302531309778+100,695.58492637380047+100,527.25758144153969+100,500.14928785356605+100,472.54177384823121+100,358.13381068808135+100,324.20323585960454+100,278.69483019677056+100,271.76079921372494+100],
    y=["ZNF562","TRIM28","ZEB1","CBX3","SETDB1","ZEB2","E2F6","ZNF671","SIN3A","CBX5","ZNF506","ZNF561","ZNF264"],
    mode="text",
    name="Text labels",
    text=["293T","WIBR3","H1975","HCT116","U2OS","K562","H1","293T","HCT116","GM12878","293T","293T","HEK293"],
    #xshift=30,
    textposition="middle right",
    textfont=dict(
        size=28,
        color="grey",
    )
))

fig.update_layout(xaxis_range=[-200,2700])

fig.update_layout(showlegend=False)

fig.update_layout(
    font_family="Helvetica",
    font_color="black",
    font_size=32,
)

fig.update_traces(marker=dict(size=18,
                              color='rgba(0, 0, 139, 0.6)',
                              line=dict(width=1,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))


fig.update_xaxes(title_standoff=20)
fig.update_yaxes(title_standoff=20)
fig.update_xaxes(tickfont_size=30, ticks="outside", ticklen=10, tickwidth=3)
fig.update_yaxes(tickfont_size=30, ticks="outside", ticklen=10, tickwidth=3)

fig.update_xaxes(showline=True, linewidth=4, linecolor='black', gridwidth=2)
fig.update_yaxes(showline=True, linewidth=4, linecolor='black', gridwidth=2)

fig.update_xaxes(fixedrange=True)

fig.show()
fig.write_image("cistromeTFs_repressors.svg")
