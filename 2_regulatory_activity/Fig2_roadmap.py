import plotly.express as px
import pandas as pd

df = pd.read_csv('SuppTable2_LTR10AF_enriched_in_Roadmap_giggle_results.tsv', sep='\t')
#print(df)

fig = px.strip(df, x="roadmap_category", y="giggle_score",
             category_orders={"roadmap_category":["1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies"]},
             labels={"giggle_score":"LTR10 Enrichment", "roadmap_category":"Roadmap Category"})

for template in ["plotly_white"]:    
    #fig.update_xaxes(categoryorder="category ascending")
    fig.update_layout(template=template)
    fig.update_layout(autosize=False,width=800,height=600)
    
fig.update_traces(quartilemethod="linear", jitter=0.8)
    
fig.update_xaxes(title_standoff=0)
fig.update_yaxes(title_standoff=5)
fig.update_xaxes(tickfont_size=18, ticks="outside", ticklen=5, tickwidth=2)
fig.update_yaxes(tickfont_size=18, ticks="outside", ticklen=5, tickwidth=2)

fig.update_xaxes(showline=True, linewidth=2.5, linecolor='black', gridwidth=1, tickangle=45)
fig.update_yaxes(showline=True, linewidth=2.5, linecolor='black', gridwidth=1)

fig.update_layout(
    font_family="Helvetica",
    font_color="black",
    font_size=16,
)

fig.show()
fig.write_image("roadmap_enrichment.svg")
