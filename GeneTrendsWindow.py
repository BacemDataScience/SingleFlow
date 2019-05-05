import os
import palantir
import pdb
import pandas as pd
import sys
import numpy as np
from sklearn.preprocessing import StandardScaler

from compute_goea import goea

from TsneWindow import TsneWindow

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt4.QtCore import QThread, Qt
from matplotlib.figure import Figure

class GeneTrendsWindow(QtWidgets.QMainWindow):

    def __init__(self, trends, pr_res, clusters, id_to_name, out_dir, goea_dir, tsne, cell_clusters, colors):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Gene trends for ")
        layout = QtWidgets.QVBoxLayout(self._main)

        self.trends = trends
        self.clusters = clusters
        self.id_to_name = id_to_name
        self.out_dir = out_dir
        self.goea_dir = goea_dir

        self.line = False
        self.lines = list()
        self.pseudotime = 0
        self.n_lines = 0

        self.colors = colors

        ##self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

        ## drop down to choose trajectory
        self.select_cluster = QtWidgets.QComboBox()
        layout.addWidget(self.select_cluster)
        self.select_cluster.addItem('Select cluster')

        for cluster in set(self.clusters):
            self.select_cluster.addItem(str(cluster))
        self.select_cluster.currentIndexChanged.connect(self._select_cluster)

        gene_trend_canvas = FigureCanvas(Figure(figsize=(5, 5)))
        layout.addWidget(gene_trend_canvas)
        self.addToolBar(NavigationToolbar(gene_trend_canvas, self))
        self._gt_fig = gene_trend_canvas.figure
        self._gt_ax = gene_trend_canvas.figure.subplots()

        ## Add listener for mouse motion to update tsne canvas
        gene_trend_canvas.mpl_connect('motion_notify_event', self._on_mouse_move)

        go_button = QtWidgets.QPushButton("GO analysis")
        go_button.clicked.connect(self._on_click_go)
        layout.addWidget(go_button)

        self.tsne_win = TsneWindow(tsne, cell_clusters, colors, out_dir, pr_res=pr_res, trends_win=self)
        self.tsne_win.show()

    def _on_click_go(self):
        c = int(self.select_cluster.currentText())
        #cluster = self.clusters[c]

        gene_symbols = [self.id_to_name[gene_id] for gene_id in self.clusters.index[self.clusters == c]]
#        print(gene_symbols)
#        gene_symbols = [self.id_to_name[gene_id] for gene_id in cluster]
        goea(self.clusters.index[self.clusters == c], gene_symbols, 'gene trends', str(c), self.goea_dir, self.out_dir) ## list of genes represented by their ensembl id and gene symbol

    def _on_mouse_move(self, event):
        if event.xdata:
            self.pseudotime = event.xdata
            self._update_line()

    def _update_line(self):
        if self.line:#len(self._gt_ax.lines) > 3:
            self._gt_ax.lines[-1].remove()
        self.line = True
        self._gt_ax.axvline(self.pseudotime)
        self._gt_ax.figure.canvas.draw()
    #        self._update_tsne(pseudotime)

    def _update_pseudotime(self, pseudotime):
        self.pseudotime = pseudotime

    def _add_line(self):
        if self.line:#len(self._gt_ax.lines) > 3:
            self._gt_ax.lines[-1].remove()
            self._gt_ax.axvline(self.pseudotime, color=self.colors[self.n_lines])
            self.n_lines += 1
            self.line = False
            self._gt_ax.figure.canvas.draw()
            self.lines = self._gt_ax.lines[-self.n_lines:]

    def _select_cluster(self,i):
        print("Items in the list are :")
        for count in range(self.select_cluster.count()):
            print(self.select_cluster.itemText(count))
        print("Current index",i,"selection changed ",self.select_cluster.currentText())
        c = int(self.select_cluster.currentText())
        print("CLUSTER:")
        print(c)

        trends = pd.DataFrame(StandardScaler().fit_transform(self.trends.T).T,
                          index=self.trends.index, columns=self.trends.columns)
        clusters = self.clusters

        gene_symbols_in_cluster = [self.id_to_name[gene_id] for gene_id in clusters.index[clusters == c]]

        means = trends.loc[clusters.index[clusters == c], :].mean()
        std = trends.loc[clusters.index[clusters == c], :].std()

        self._gt_ax.clear()
        self.line = False

        # Plot all trends
        for g in clusters.index[clusters == c]:
            self._gt_ax.plot(means.index, np.ravel(
                    trends.loc[g, :]), linewidth=0.5, color='lightgrey')
        self._gt_ax.plot(means.index, np.ravel(means), color='#377eb8')
        self._gt_ax.plot(means.index, np.ravel(means - std), linestyle='--',
                color='#377eb8', linewidth=0.75)
        self._gt_ax.plot(means.index, np.ravel(means + std), linestyle='--',
                color='#377eb8', linewidth=0.75)
        self._gt_ax.set_title('Cluster {}'.format(c), fontsize=12)
        self._gt_ax.tick_params('both', length=2, width=1, which='major')
        self._gt_ax.tick_params(axis='both', which='major', labelsize=8, direction='in')
#        self._gt_ax.set_xticklabels([])

        for line in self.lines:
            self._gt_ax.lines.append(line)

        self._gt_ax.figure.canvas.draw()
