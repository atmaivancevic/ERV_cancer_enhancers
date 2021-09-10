import plotly.express as px
import pandas as pd

df = pd.read_csv('jalview_pca_transformed_points_allLTR10_min200_alignment.csv', sep=',')
#print(df)

fig = px.scatter(df, x="PC1", y="PC2",
                color="subfamily",
                color_discrete_sequence=["red", "green", "blue", "purple", "magenta","orange","grey"],
                hover_data=df.columns,
                title="Clusterf*ck",
                marginal_y="violin",
                marginal_x="histogram",
                )

fig.for_each_trace(lambda t: t.update(name=t.name.replace("=",": ")))
            
fig.show()
#py.plot(fig, filename = 'clusterfck', auto_open=True) 
