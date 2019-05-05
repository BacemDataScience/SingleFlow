import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

def plot_cell_clusters(tsne, clusters, colors):
    """Plot cell clusters on the tSNE map
    :param tsne: tSNE map
    :param clusters: Results of the determine_cell_clusters function
    """

    # Cluster colors
    n_clusters = len(set(clusters))
    cluster_colors = pd.Series(colors[:n_clusters], index = set(clusters))

    # Set up figure
    n_cols = 6
    n_rows = int(np.ceil(n_clusters / n_cols))
    fig = plt.figure(figsize=[2*n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    # Clusters
    ax = plt.subplot(gs[0:2, 2:4])
    ax.scatter(tsne['x'], tsne['y'], s=3,
               color=cluster_colors[clusters[tsne.index]])
    ax.set_axis_off()

    # Branch probabilities
    for i, cluster in enumerate(set(clusters)):
        row = int(np.floor(i/n_cols))
        ax = plt.subplot(gs[row+2, i % n_cols])
        ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3, color='lightgrey')
        cells = clusters.index[clusters == cluster]
        ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'],
                   s=3, color=cluster_colors[cluster])
        ax.set_axis_off()
        ax.set_title(cluster, fontsize=10)

def get_fig(fig=None, ax=None, figsize=[4, 4]):
    """fills in any missing axis or figure with the currently active one
    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object
    """
    if not fig:
        fig = plt.figure(figsize=figsize)
    if not ax:
        ax = plt.gca()
    return fig, ax

def highlight_cells_on_tsne(tsne, cells, fig=None, ax=None):
    """    Function to highlight specific cells on the tSNE map
    """
    fig, ax = get_fig(fig=fig, ax=ax)
    ax.scatter(tsne['x'], tsne['y'], s=3, color='lightgrey')
    ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'], s=3, color='b')
    ax.set_axis_off()
    return fig, ax


def plot_gene_expression(data, tsne, genes):
    """ Plot gene expression on tSNE maps
    :param genes: Iterable of strings to plot on tSNE
    """

    not_in_dataframe = set(genes).difference(data.columns)
    if not_in_dataframe:
        if len(not_in_dataframe) < len(genes):
            print('The following genes were either not observed in the experiment, '
                  'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
        else:
            print('None of the listed genes were observed in the experiment, or the '
                  'wrong symbols were used.')
            return 0

    # remove genes missing from experiment
    genes = set(genes).difference(not_in_dataframe)

    out_dir = '/Users/herman/dev/Bioinf/master-scrna/nextflow_scripts/results'

    # Plot
    fig = FigureGrid(len(genes), 5)
    points2 = None

    for g, ax in zip(genes, fig):
        points2 = ax.scatter(tsne['x'], tsne['y'], c=(int(len(tsne['x'])/2)*[0] + int(len(tsne['x'])/2)*[1]), s=3,
                    cmap=matplotlib.cm.Spectral_r)
        points = ax.scatter(tsne['x'], tsne['y'], c=data.loc[tsne.index,  g], s=3,
                   cmap=matplotlib.cm.Spectral_r)

        #plt.colorbar(z1_plot,cax=ax1)
        ax.set_axis_off()
        ax.set_title(g)

    t = plt.colorbar(points2, label='expression', ax=ax)
    plt.savefig(os.path.join(out_dir, 'test_expression.svg'))


class FigureGrid:
    """
    Generates a grid of axes for plotting
    axes can be iterated over or selected by number. e.g.:
    >>> # iterate over axes and plot some nonsense
    >>> fig = FigureGrid(4, max_cols=2)
    >>> for i, ax in enumerate(fig):
    >>>     plt.plot(np.arange(10) * i)
    >>> # select axis using indexing
    >>> ax3 = fig[3]
    >>> ax3.set_title("I'm axis 3")
    """

    # Figure Grid is favorable for displaying multiple graphs side by side.

    def __init__(self, n: int, max_cols=3, scale=3):
        """
        :param n: number of axes to generate
        :param max_cols: maximum number of axes in a given row
        """

        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * scale, self.nrows * scale

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def colorbar(self, points, label=None, ax=None):
        self.figure.colorbar(points, label=label, ax=ax)

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]
