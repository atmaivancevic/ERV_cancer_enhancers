import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os
from IPython.display import Image

df = pd.read_csv('Top_Cistrome_Activators.tab', sep='\t')
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
    fig.update_layout(autosize=False,width=800,height=1200)
    
fig.update_yaxes(categoryorder="max ascending")

fig.add_trace(go.Scatter(
    x=[1058.4851545431459+100,1025.5891323526377+100,1022.2168157273652+100,982.13877172309136+100,960.91866554674981+100,874.68897475225331+100,759.38960030696413+100,748.38794959685396+100,693.13553734886261+100,565.11408622986153+100,529.63648017480485+100,498.50114417391795+100,495.78839604835667+100,488.70018997736251+100,418.86305874263667+100,269.54112156046084+100,185.3415913447491+100],
    y=["ATF3","FOSL1","USF1","POLR2A","JUND","PPARG","SRF","NR1H3","CEBPB","NR1H2","EP300","TEAD4","NIPBL","ELF1","EGR1","FOXO4","JUN"],
    mode="text",
    name="Text labels",
    text=["HCT116","HCT116","HCT116","HT29","HCT116","HT29","HCT116","HT29","HCT116","HT29","HCT116","HCT116","HCT116","HCT116","HCT116","LoVo","LoVo"],
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
#fig.write_image("cistromeTFs_activators.svg")
