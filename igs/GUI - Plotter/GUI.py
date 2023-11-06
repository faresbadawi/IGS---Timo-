import numpy as np
import sys
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QLabel, QLineEdit, QPushButton, QComboBox, QTabWidget, QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class BSpinePlotter(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("B-Spline Basisfunktionen")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        splitter = QSplitter(central_widget)
        central_widget_layout = QHBoxLayout(central_widget)
        central_widget_layout.addWidget(splitter)

        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        grad_label = QLabel("Grad (p):", left_panel)
        left_layout.addWidget(grad_label)

        self.grad_entry = QLineEdit(left_panel)
        self.grad_entry.setText("3")
        left_layout.addWidget(self.grad_entry)

        knotenvektor_label = QLabel("Knotenvektor (kommagetrennt):", left_panel)
        left_layout.addWidget(knotenvektor_label)

        self.knotenvektor_entry = QLineEdit(left_panel)
        self.knotenvektor_entry.setText("0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1")
        left_layout.addWidget(self.knotenvektor_entry)

        ableitungsordnung_label = QLabel("Ableitungsordnung (k):", left_panel)
        left_layout.addWidget(ableitungsordnung_label)

        self.ableitungsordnung_combobox = QComboBox(left_panel)
        self.ableitungsordnung_combobox.addItems(["1", "2", "3", "4", "5", "6"])
        self.ableitungsordnung_combobox.setCurrentText("1")
        left_layout.addWidget(self.ableitungsordnung_combobox)

        basis_button = QPushButton("Basisfunktionen plotten", left_panel)
        basis_button.clicked.connect(self.basisfunktionen_plot)
        left_layout.addWidget(basis_button)

        ableitungen_button = QPushButton("Ableitungen plotten", left_panel)
        ableitungen_button.clicked.connect(self.ableitungen_plot)
        left_layout.addWidget(ableitungen_button)

        save_button = QPushButton("Bild speichern", left_panel)
        save_button.clicked.connect(self.save_plot)
        left_layout.addWidget(save_button)

        left_panel.setLayout(left_layout)

        right_panel = QTabWidget()
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)

        self.tabs = right_panel

        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)

        self.tabs.currentChanged.connect(self.tab_changed)

    def tab_changed(self):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, PlotCanvas):
            current_tab.draw()

    def close_tab(self, index):
        widget = self.tabs.widget(index)
        if widget is not None:
            widget.deleteLater()
            self.tabs.removeTab(index)

    def basisfunktionen_plot(self):
        p = int(self.grad_entry.text())
        knotenvektor_str = self.knotenvektor_entry.text()
        knotenvektor = list(map(float, knotenvektor_str.split(',')))
        t_werte = np.linspace(min(knotenvektor), max(knotenvektor), 1000)
        basisfunktions_werte = np.zeros((len(t_werte), len(knotenvektor) - p - 1))

        for i in range(len(knotenvektor) - p - 1):
            for j in range(len(t_werte)):
                basisfunktions_werte[j, i] = self.basisfunktion(t_werte[j], knotenvektor, i, p)

        tab_index = self.tabs.addTab(PlotCanvas(self), f'Basisfunktionen (p={p})')
        self.tabs.setCurrentIndex(tab_index)
        self.tabs.currentWidget().plot_basiskurven(t_werte, basisfunktions_werte, p)

    def ableitungen_plot(self):
        p = int(self.grad_entry.text())
        knotenvektor_str = self.knotenvektor_entry.text()
        knotenvektor = list(map(float, knotenvektor_str.split(',')))

        n = len(knotenvektor) - p - 1
        k = int(self.ableitungsordnung_combobox.currentText()) - 1

        t_werte = np.linspace(min(knotenvektor), max(knotenvektor), 1000)
        ableitungen_werte = np.zeros((len(t_werte), n))

        for i in range(n):
            for j in range(len(t_werte)):
                ableitungen_werte[j, i] = self.ableitung_kter_ordnung(t_werte[j], knotenvektor, i, p, k)

        tab_index = self.tabs.addTab(PlotCanvas(self), f'Ableitungen (p={p}, k={k + 1})')
        self.tabs.setCurrentIndex(tab_index)
        self.tabs.currentWidget().plot_ableitungen(t_werte, ableitungen_werte, n, k)

    def basisfunktion(self, t, knotenvektor, i, p):
        if knotenvektor[-1] and i == len(knotenvektor) - p - 1:
            return 1.0
        else:
            if p == 0:
                return 1 if knotenvektor[i] <= t < knotenvektor[i + 1] else 0
            else:
                denom1 = knotenvektor[i + p] - knotenvektor[i]
                denom2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
                a = 0 if denom1 == 0 else (t - knotenvektor[i]) / denom1
                b = 0 if denom2 == 0 else (knotenvektor[i + p + 1] - t) / denom2
                return a * self.basisfunktion(t, knotenvektor, i, p - 1) + b * self.basisfunktion(t, knotenvektor, i + 1, p - 1)

    def ableitungsbasisfunktion(self, t, knotenvektor, i, p):
        if p == 0:
            return 0
        else:
            denom1 = knotenvektor[i + p] - knotenvektor[i]
            denom2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
            a = 0 if denom1 == 0 else p / denom1
            b = 0 if denom2 == 0 else -p / denom2
            return a * self.basisfunktion(t, knotenvektor, i, p - 1) + b * self.basisfunktion(t, knotenvektor, i + 1, p - 1)

    def ableitung_kter_ordnung(self, t, knotenvektor, i, p, k):
        if k == 0:
            return self.ableitungsbasisfunktion(t, knotenvektor, i, p)
        else:
            denom1 = knotenvektor[i + p] - knotenvektor[i]
            denom2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
            a = 0 if denom1 == 0 else p / denom1
            b = 0 if denom2 == 0 else -p / denom2
            return a * self.ableitung_kter_ordnung(t, knotenvektor, i, p - 1, k - 1) + b * self.ableitung_kter_ordnung(t, knotenvektor, i + 1, p - 1, k - 1)

    def save_plot(self):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, PlotCanvas):
            file_dialog = QFileDialog()
            options = QFileDialog.Options()
            options |= QFileDialog.ReadOnly
            file_path, _ = file_dialog.getSaveFileName(self, "Bild speichern", "", "PNG Files (*.png);;All Files (*)", options=options)
            if file_path:
                current_tab.figure.savefig(file_path, dpi=100)
                print(f"Das Bild wurde unter {file_path} gespeichert.")
        else:
            print("Kein gültiges Plot-Tab ausgewählt.")

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.figure = Figure(figsize=(10, 6))
        super().__init__(self.figure)
        self.setParent(parent)

    def plot_basiskurven(self, t, kurven, p):
        ax = self.figure.add_subplot(111)
        ax.clear()
        for i in range(len(kurven[0])):
            ax.plot(t, kurven[:, i], label=f'Basis {i}')
        ax.set_title(f'B-Spline Basisfunktionen (p={p})')
        ax.set_xlabel('t')
        ax.set_ylabel('Wert')
        ax.legend()
        ax.grid(False)
        ax.spines['top'].set_position('zero')
        ax.spines['right'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_tick_params(direction='out')
        ax.yaxis.set_tick_params(direction='out')
        self.draw()

    def plot_ableitungen(self, t, ableitungen, n, k):
        ax = self.figure.add_subplot(111)
        ax.clear()
        for i in range(n):
            ax.plot(t, ableitungen[:, i], label=f'Ableitung {i} der Basis')
        ax.set_title(f'Ableitungen der B-Spline Basisfunktionen (k={k + 1})')
        ax.set_xlabel('t')
        ax.set_ylabel('Wert')
        ax.legend()
        ax.grid(False)
        ax.spines['top'].set_position('zero')
        ax.spines['right'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_tick_params(direction='out')
        ax.yaxis.set_tick_params(direction='out')
        self.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = BSpinePlotter()
    window.show()
    sys.exit(app.exec_())
