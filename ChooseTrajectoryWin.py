import pandas as pd
import palantir

from sklearn.neighbors import NearestNeighbors

import numpy as np

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from PyQt5.QtCore import QThread, pyqtSignal, Qt

from matplotlib.figure import Figure

class ChooseTrajectoryWin(QtWidgets.QMainWindow):
    def __init__(self, pr_trajectories, trajectories):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Choose trajectory for gene trend analysis")
        layout = QtWidgets.QVBoxLayout(self._main)

        self.trajectories = trajectories

        self.done_button = QtWidgets.QPushButton("Done")
        self.done_button.clicked.connect(self.on_click_done)
        layout.addWidget(self.done_button)

        self.cb = QtWidgets.QComboBox()
        self.cb.addItems(pr_trajectories)
        layout.addWidget(self.cb)

    def on_click_done(self):
        self.trajectories.append(self.cb.currentText())
        self.close()
