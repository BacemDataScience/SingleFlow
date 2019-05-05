import os
import pandas as pd
import numpy as np

from compute_goea import goea

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from PyQt5.QtCore import QThread, pyqtSignal, Qt

from matplotlib.figure import Figure

class DEGWindow(QtWidgets.QMainWindow):
    def __init__(self, df, id_to_name, cluster1, cluster2, out_dir, goea_dir, adj_p_vals=True):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("DEG result for " + cluster1 + ' and ' + cluster2)
        layout = QtWidgets.QVBoxLayout(self._main)

        scroll = QtWidgets.QScrollArea()
        table = QtWidgets.QTableWidget()
        scroll.setWidget(table)
        layout.addWidget(table)

        self.df = df
        self.id_to_name = id_to_name
        self.cluster1 = cluster1
        self.cluster2 = cluster2
        self.out_dir = out_dir
        self.goea_dir = goea_dir

        if adj_p_vals:
            self.p_vals = 'adj-p-value'
        else:
            self.p_vals = 'p-value'

        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        from rpy2.robjects import pandas2ri
        from rpy2.robjects import r
        import rpy2.robjects as robjects
        robjects.r('''
        p <- function(ediff) {
                library(scde)
                p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
                p.values
        }
        ''')

        robjects.r('''
        p.adj <- function(ediff) {
                library(scde)
                p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
                p.values.adj
        }
        ''')

        r_p = robjects.globalenv['p']
        r_p_adj = robjects.globalenv['p.adj']
        pandas2ri.activate()
        p_vals = r_p(df)
        adj_p_vals = r_p_adj(df)

        self.df['p-value'] = pd.Series(pandas2ri.ri2py(p_vals), index = self.df.index)
        self.df['adj-p-value'] = pd.Series(pandas2ri.ri2py(adj_p_vals), index = self.df.index)

        gene_ids = self.df.index
        gene_symbols = np.array([self.id_to_name[gene_id] for gene_id in gene_ids])
        self.df['gene_symbols'] = pd.Series(gene_symbols, index = self.df.index)

        self.df_up = self.df.loc[df['Z'] > 0]
        self.df_down = self.df.loc[df['Z'] < 0]

        writer = pd.ExcelWriter(os.path.join(self.out_dir, self.cluster1 + ' deg ' + self.cluster2 + '.xlsx'))
        self.df.to_excel(writer,'Up and down')
        self.df_up.to_excel(writer,'Up')
        self.df_down.to_excel(writer,'Down')
        writer.save()

        ## Get fold change

        ## Data frame with fold change and p values

        volcano_canvas = FigureCanvas(Figure(figsize=(5, 5)))
        layout.addWidget(volcano_canvas)
        self.addToolBar(QtCore.Qt.BottomToolBarArea,
                        NavigationToolbar(volcano_canvas, self))
        self._volcano_fig = volcano_canvas.figure
        self._volcano_ax = volcano_canvas.figure.subplots()

        fig, self._volcano_ax = self.volcano_plot(df, fig=self._volcano_fig, ax=self._volcano_ax)

        self._volcano_ax.figure.canvas.draw()

        go_button = QtWidgets.QPushButton("GO analysis")
        go_button.clicked.connect(self._on_click_go)
        layout.addWidget(go_button)

        table.setColumnCount(len(df.columns))
        table.setRowCount(len(df.index))
        for i in range(len(df.index)):
            table.setVerticalHeaderItem(i,QtWidgets.QTableWidgetItem(self.id_to_name[df.index[i]]))
            for j in range(len(df.columns)):
                table.setHorizontalHeaderItem(j,QtWidgets.QTableWidgetItem(df.columns[j]))
                table.setItem(i,j,QtWidgets.QTableWidgetItem(str(df.iloc[i, j])))

    def get_fig(self, fig=None, ax=None, figsize=[4, 4]):
        """fills in any missing axis or figure with the currently active one
        :param ax: matplotlib Axis object
        :param fig: matplotlib Figure object
        """
        if not fig:
            fig = plt.figure(figsize=figsize)
        if not ax:
            ax = plt.gca()
        return fig, ax

    def volcano_plot(self, df, fig=None, ax=None):
        """Function to highlight specific cells on the tSNE map
        """
        fig, ax = self.get_fig(fig=fig, ax=ax)
        df['log-p'] = -df[self.p_vals].apply(np.log10)
        ax.scatter(df['ce'], df['log-p'], s=5, color='lightgrey')
        ax.scatter(df.loc[(df[self.p_vals] < 0.05) & (df['ce'] > 1), 'ce'],
                            df.loc[(df[self.p_vals] < 0.05) & (df['ce'] > 1), 'log-p'], s=5, color=self.cluster1)
        ax.scatter(df.loc[(df[self.p_vals] < 0.05) & (df['ce'] < -1), 'ce'],
                            df.loc[(df[self.p_vals] < 0.05) & (df['ce'] < -1), 'log-p'], s=5, color=self.cluster2)
        fig.savefig(os.path.join(self.out_dir, self.cluster1 + '_' + self.cluster2 + '_volcano.svg'))
        return fig, ax

    def _on_click_go(self):
        cluster = self.df.loc[(self.df[self.p_vals] < 0.05)].index
        print('up and down ', len(cluster), ' number of genes')
        gene_symbols = [self.id_to_name[gene_id] for gene_id in cluster]
        goea(cluster, gene_symbols, self.cluster1, self.cluster2, self.goea_dir, self.out_dir) ## list of genes represented by their ensembl id and gene symbol

        cluster = self.df.loc[(self.df[self.p_vals] < 0.05) & (self.df['ce'] > 1)].index
        print('up ', len(cluster), ' number of genes')
        gene_symbols = [self.id_to_name[gene_id] for gene_id in cluster]
        goea(cluster, gene_symbols, self.cluster1, str(self.cluster2) + 'up', self.goea_dir, self.out_dir) ## list of genes represented by their ensembl id and gene symbol

        cluster = self.df.loc[(self.df[self.p_vals] < 0.05) & (self.df['ce'] < -1)].index
        print('down ', len(cluster), ' number of genes')
        gene_symbols = [self.id_to_name[gene_id] for gene_id in cluster]
        goea(cluster, gene_symbols, self.cluster1, str(self.cluster2) + 'down', self.goea_dir, self.out_dir)## list of genes represented by their ensembl id and gene symbol
