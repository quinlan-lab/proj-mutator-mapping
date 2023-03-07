# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go

import pandas as pd
import numpy as np
from dash.dependencies import Input, Output

from ihd.utils import (
    compute_spectra,
    compute_genotype_similarity,
    perform_ihd_scan,
    perform_permutation_test,
)
from ihd.schema import MutationSchema, IHDResultSchema

app = dash.Dash(__name__)

markers = pd.read_csv("data/genotypes/bxd.markers")
geno = pd.read_csv("data/genotypes/bxd.geno")
mutations = pd.read_csv("data/mutations/bxd/annotated_filtered_singletons.csv")
MutationSchema.validate(mutations)

samples = mutations['sample'].unique()
# get the overlap between those and the sample names in the genotype data
samples_overlap = list(set(samples).intersection(set(geno.columns)))

cols2use = ["marker"]
cols2use.extend(samples_overlap)
geno = geno[cols2use]

mutations_filtered = mutations[mutations['sample'].isin(samples_overlap)]

# convert string genotypes to integers based on config definition
replace_dict = {"B":0, "D": 2, "H": 1, "U": np.nan}
geno_asint = geno.replace(replace_dict).replace({1: np.nan})

geno_filtered = geno_asint[samples_overlap].values
markers_filtered = geno_asint['marker'].values

app.layout = html.Div([
    dcc.Markdown("""
    ### Discovering epistasis between germline mutator alleles in mice

    Thomas A. Sasani, Aaron R. Quinlan, Kelley Harris\
    """),
    html.Div([
        html.Div(
            className="row",
            children=[
                dcc.Markdown(
                    """**Select the k-mer size to use to define the mutatin spectrum:**""",
                    style={'width': '50%'}),
                dcc.Markdown(
                    """**Select the number of permutations you want to use to define significance:**""",
                    style={'width': '50%'}),
            ],
            style={'display': 'flex'},
        ),
        html.Div(className="row",
                 children=[
                     dcc.RadioItems(
                         id="kmer-size-radio",
                         options=[
                             {
                                 'label': '1-mer',
                                 'value': 1,
                             },
                             {
                                 'label': '3-mer',
                                 'value': 3
                             },
                         ],
                         value=1,
                         style={'width': '50%'},
                     ),
                     dcc.Dropdown(id="num-perms-dropdown", options=[100, 500, 1_000, 10_000], value=100, style={'width': '50%'})
                 ],
                 style={'display': 'flex'}),

        html.Button('Submit', id='submit-val', n_clicks=0),
        dcc.Store(id='genotype-similarity-data'),
        html.Div(
            className="row",
            children=[
                dcc.Graph(
                    id='graph-with-dropdown',
                    hoverData={'points': [{
                        'customdata': 'rs31879829'
                    }]},
                    style={'width': '100%'},
                ),

            ],
            style={'display': 'flex'},
        ),
    ]),

])

@app.callback(
        Output('genotype-similarity-data', 'data'),
        Input('submit-val', 'n-clicks'),
)
def store_genotype_similarity(n_clicks):
    if n_clicks > 0:
        genotype_similarity = compute_genotype_similarity(geno_filtered)
    return genotype_similarity.to_json()

@app.callback(
    Output('graph-with-dropdown', 'figure'),
    [
        Input('kmer-size-radio', 'value'),
        Input('num-perms-dropdown', 'value'),
        Input('genotype-similarity-data', 'data')
    ],
)
def update_scan_plot(k, p):
    samples, mutation_types, spectra = compute_spectra(mutations, k=k)
    # compute the maximum cosine distance between groups of
    # haplotypes at each site in the genotype matrix
    focal_dists = perform_ihd_scan(
        spectra,
        geno_filtered,
        genotype_similarity,
    )

    res_df = pd.DataFrame({
        'marker': markers_filtered,
        'Distance': focal_dists,
        'k': k,
    })
    IHDResultSchema.validate(res_df)

    res_df = res_df.merge(markers_filtered, on="marker").sort_values(["chromosome", "Mb"])
    res_df["marker_number"] = np.arange(res_df.shape[0])

    chrom_to_max_idx = res_df.groupby('chromosome').max().to_dict()['marker_number']
    x_idxs, x_labels = [], []
    prev_idx = 0
    for chrom in np.arange(1, 20):
        max_idx = chrom_to_max_idx[str(chrom)]
        idx_to_use = (max_idx + prev_idx)  / 2
        x_idxs.append(idx_to_use)
        x_labels.append(str(chrom))
        prev_idx = max_idx

    fig = px.scatter(
        res_df,
        x="marker_number",
        y="Distance",
        #color="color",
        hover_name="marker",
        hover_data=["chromosome", "Mb", "Distance"],
        title="IHD scan results (mm10/GRCm38)",
        labels={
            'marker': "Genotype marker on mm10",
            'Distance': "Adjusted inter-haplotype " + r"$\chi^2$" + " statistic",
        },
    )

    fig.update_layout(transition_duration=500,
                      showlegend=False,
                      xaxis=dict(
                          tickmode='array',
                          tickvals=x_idxs,
                          ticktext=x_labels,
                      ))

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
