import pandas as pd
import palantir

from sklearn.neighbors import NearestNeighbors

import os
import numpy as np

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from PyQt5.QtCore import QThread, pyqtSignal, Qt

from matplotlib.figure import Figure

class TsneWindow(QtWidgets.QMainWindow):
    def __init__(self, tsne, cell_clusters, colors, out_dir, pr_res=None, trends_win=None):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        layout = QtWidgets.QVBoxLayout(self._main)

        self.tsne = tsne
        self.pr_res = pr_res
        self.cell_clusters = cell_clusters
        self.colors = colors
        self.out_dir = out_dir
        self.pseudotime = 0
        self.trends_win = trends_win

        self.done_button = QtWidgets.QPushButton("Done")
        self.done_button.clicked.connect(self.on_click_done)
        layout.addWidget(self.done_button)

        if pr_res:
            self.setWindowTitle("tSNE -- pseudotime")

            self.set_pseudotime = QtWidgets.QLineEdit(self)
            layout.addWidget(self.set_pseudotime)
            self.set_pseudotime_button = QtWidgets.QPushButton('set pseudotime')
            self.set_pseudotime_button.clicked.connect(self._update_pseudotime)
            layout.addWidget(self.set_pseudotime_button)

            self.slider = QtWidgets.QSlider(Qt.Horizontal)
            self.slider.setFocusPolicy(Qt.StrongFocus)
            self.slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
            self.slider.setMinimum(0)
            self.slider.setMaximum(10000)
            self.slider.setTickInterval(1)
            self.slider.setSingleStep(1)
            self.slider.valueChanged.connect(self._update_tsne)
            layout.addWidget(self.slider)

            self.pseudotime_cluster_button = QtWidgets.QPushButton('Set pseudotime cluster as cluster')
            self.pseudotime_cluster_button.clicked.connect(self._on_pseudotime_cluster_click)
            layout.addWidget(self.pseudotime_cluster_button)
        else:
            self.setWindowTitle("tSNE -- define custom clusters")

        tsne_canvas = FigureCanvas(Figure(figsize=(5, 5)))
        layout.addWidget(tsne_canvas)
        self.addToolBar(QtCore.Qt.BottomToolBarArea,
                        NavigationToolbar(tsne_canvas, self))
        self._tsne_fig = tsne_canvas.figure
        self._tsne_ax = tsne_canvas.figure.subplots()

        tsne_canvas.mpl_connect('button_release_event', self.on_mouse_click)

        fig, self._tsne_ax = palantir.plot.plot_tsne(self.tsne, fig=self._tsne_fig, ax=self._tsne_ax)
        self._tsne_ax.figure.canvas.draw()

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

    def highlight_cells_on_tsne(self, tsne, clusters, fig=None, ax=None):
        fig, ax = self.get_fig(fig=fig, ax=ax)
        ax.scatter(tsne['x'], tsne['y'], s=5, color='lightgrey')
        for i, cluster in enumerate(clusters):
            ax.scatter(tsne.loc[cluster, 'x'], tsne.loc[cluster, 'y'], s=3, color=self.colors[i])
        ax.set_axis_off()
        return fig, ax

    def highlight_cells_on_tsne2(self, tsne, cells, fig=None, ax=None):
        """    Function to highlight specific cells on the tSNE map
        """
        fig, ax = self.get_fig(fig=fig, ax=ax)
        ax.scatter(tsne['x'], tsne['y'], s=5, color='lightgrey')
        ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'], s=5)
        ax.set_axis_off()
        return fig, ax

    def on_mouse_click(self, event):
        self.n_neighbors = 100
        neigh = NearestNeighbors(n_neighbors=self.n_neighbors)
        neigh.fit(self.tsne)

        cluster = self.tsne.index[neigh.kneighbors([pd.Series([float(event.xdata), float(event.ydata)], ['x','y'])])[1][0]]
        """
        with open(os.path.join(self.pickle_input, 'exclude.pickle'), 'wb') as output_file:
            pickle.dump(cluster, output_file)
        """

        self.cell_clusters.append(cluster)
        self._tsne_ax.clear()
        fig, self._tsne_ax = self.highlight_cells_on_tsne(self.tsne, self.cell_clusters, fig=self._tsne_fig, ax=self._tsne_ax)
        #fig, self._tsne_ax = self.highlight_cells_on_tsne(self.tsne, cluster,fig=self._tsne_fig, ax=self._tsne_ax)
        self._tsne_ax.figure.canvas.draw()
        print('clicked')

    def _on_pseudotime_cluster_click(self):
        print('set pseudotime cluster as cluster')
        self.cell_clusters.append(self.pseudotime_cluster)
        self._tsne_ax.clear()
        fig, self._tsne_ax = self.highlight_cells_on_tsne(self.tsne, self.cell_clusters, fig=self._tsne_fig, ax=self._tsne_ax)
        self._tsne_ax.figure.canvas.draw()
        if self.trends_win:
            self.trends_win._add_line()

    def on_click_done(self):
        self._tsne_ax.figure.savefig(os.path.join(self.out_dir, 'custom_clusters.svg'))
        self.close()

    def _update_pseudotime(self):
        self.pseudotime = self.slider.setValue(int(float(self.set_pseudotime.text())*10000))
        self._update_tsne

    def _update_tsne(self):
        self.pseudotime = self.slider.value()/10000
        print('Pseudotime: ', self.pseudotime)
        self._tsne_ax.clear()
        self.pseudotime_cluster= self.tsne.index[(self.pr_res.pseudotime > (self.pseudotime - 0.005)) & (self.pr_res.pseudotime < (self.pseudotime + 0.005))]

        fig, self._tsne_ax = self.highlight_cells_on_tsne2(self.tsne, self.pseudotime_cluster, fig=self._tsne_fig, ax=self._tsne_ax)
        self._tsne_ax.figure.canvas.draw()

        if self.trends_win:
            self.trends_win._update_pseudotime(self.pseudotime)
            self.trends_win._update_line()
