import plotly.express as px
from sklearn.manifold import TSNE
import pandas as pd


df = pd.read_csv(r'vcfTracks-bcftoolsmerge-pbsv.txt', index_col=False)  # , header=None
tsne = TSNE(n_components=2, random_state=42)
n_rows=len(df)
n_columns=len(df.columns)

#X_tsne = tsne.fit_transform( df.iloc[ :, 13:22 ] )
#fig = px.scatter(x=X_tsne[:, 0], y=X_tsne[:, 1])

fig = px.scatter(x=df.iloc[:, 14], y=df.iloc[:, 22])

fig.show()
