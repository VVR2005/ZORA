"""Qt widgets, theming, and external tool integration for ZORA."""

import os
import importlib.util
from typing import Dict, List, Optional, Union

import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide6.QtWidgets import (
    QFrame, QLabel, QVBoxLayout, QFileDialog,
)
from PySide6.QtCore import Qt, Signal, QRect
from PySide6.QtGui import QFont, QColor, QAction, QSyntaxHighlighter, QTextCharFormat, QDragEnterEvent, QDropEvent

from zora.models import SequenceRecord, SeqUtils
from zora.project import FileParser


# ============================================================
# Sequence Syntax Highlighter
# ============================================================
class SequenceHighlighter(QSyntaxHighlighter):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._rules = []
        base_colors = {
            'A': QColor('#22C55E'),
            'T': QColor('#EF4444'),
            'C': QColor('#3B82F6'),
            'G': QColor('#EAB308'),
            'N': QColor('#6B7280'),
        }
        import re
        for base, color in base_colors.items():
            fmt = QTextCharFormat()
            fmt.setForeground(color)
            fmt.setFontWeight(QFont.Bold)
            self._rules.append((re.compile(f'[{base.lower()}{base}]'), fmt))

    def highlightBlock(self, text):
        for pattern, fmt in self._rules:
            for match in pattern.finditer(text):
                start = match.start()
                length = match.end() - match.start()
                self.setFormat(start, length, fmt)


# ============================================================
# Matplotlib Canvas
# ============================================================
class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=9, height=5.5, dpi=120):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)

    @staticmethod
    def style_ax(ax, theme=None, grid=True):
        if grid:
            ax.grid(True, alpha=0.15, linestyle='--', linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if theme:
            border = theme.get('BORDER', '#1F2A40')
            ax.spines['left'].set_color(border)
            ax.spines['bottom'].set_color(border)
            ax.tick_params(colors=theme.get('TEXT_MUTED', '#9CA3AF'), labelsize=9)
            ax.yaxis.label.set_color(theme.get('TEXT_HEADER', '#8B9DC3'))
            ax.xaxis.label.set_color(theme.get('TEXT_HEADER', '#8B9DC3'))
            ax.title.set_color(theme.get('ACCENT', '#00D1FF'))


# ============================================================
# Drag-Drop Sequence Input
# ============================================================
class DropArea(QFrame):
    sequences_dropped = Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.setMinimumHeight(100)
        self.label = QLabel("Drop FASTA/FASTQ/GenBank files here\nor click to browse", self)
        self.label.setAlignment(Qt.AlignCenter)
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        self.setObjectName("dropArea")
        self.label.setStyleSheet("color: #6B7280; font-size: 13px;")

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dragLeaveEvent(self, event):
        pass

    def dropEvent(self, event: QDropEvent):
        records = []
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path:
                records.extend(FileParser.parse_file(path))
        if records:
            self.sequences_dropped.emit(records)
        from PySide6.QtCore import QTimer
        QTimer.singleShot(2000, lambda: None)

    def mousePressEvent(self, event):
        files, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "",
            "Sequence Files (*.fasta *.fa *.fastq *.gb *.genbank *.txt *.csv)")
        if files:
            records = FileParser.parse_file(files)
            if records:
                self.sequences_dropped.emit(records)


# ============================================================
# External Tool Integration
# ============================================================
class ExternalToolLauncher:
    @staticmethod
    def find_executable(name: str) -> Optional[str]:
        import shutil
        return shutil.which(name)

    @staticmethod
    def launch_pymol(seq: Optional[str] = None, filename: Optional[Union[str, List[str]]] = None):
        exe = ExternalToolLauncher.find_executable('pymol')
        if not exe:
            return None, "PyMOL not found. Install it or add to PATH."
        args = [exe]
        if filename:
            if isinstance(filename, list):
                args.extend(filename)
            else:
                args.append(filename)
        import subprocess
        subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True, None

    @staticmethod
    def launch_chimerax(seq: Optional[str] = None, filename: Optional[str] = None):
        exe = ExternalToolLauncher.find_executable('chimerax') or ExternalToolLauncher.find_executable('ChimeraX')
        if not exe:
            return None, "ChimeraX not found."
        args = [exe]
        if filename:
            args.append(filename)
        import subprocess
        subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True, None

    @staticmethod
    def launch_snapgene(filename: Optional[str] = None):
        exe = ExternalToolLauncher.find_executable('snapgene')
        if not exe:
            return None, "SnapGene not found."
        import subprocess
        subprocess.Popen([exe, filename] if filename else [exe],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True, None

    @staticmethod
    def export_for_snapgene(records: List[SequenceRecord], path: str):
        with open(path, 'w') as f:
            for r in records:
                f.write(f">{r.name} {r.description}\n")
                for i in range(0, len(r.sequence), 80):
                    f.write(r.sequence[i:i+80] + "\n")

    @staticmethod
    def export_for_chimerax(records: List[SequenceRecord], path: str):
        ExternalToolLauncher.export_for_snapgene(records, path)

    @staticmethod
    def send_to_pymol_via_file(records: List[SequenceRecord], path: str):
        with open(path, 'w') as f:
            f.write("# PyMOL script generated by ZORA\n")
            for i, r in enumerate(records):
                f.write(f"# Sequence {i+1}: {r.name}\n")
                f.write(f"# Length: {len(r)} bp\n")
                f.write(f"# GC%: {SeqUtils.gc_content(r.sequence):.1f}\n")

    @staticmethod
    def get_available_tools() -> Dict[str, bool]:
        import shutil
        return {
            'PyMOL': shutil.which('pymol') is not None,
            'ChimeraX': shutil.which('chimerax') is not None or shutil.which('ChimeraX') is not None,
            'SnapGene': shutil.which('snapgene') is not None,
            'AutoDock Tools': shutil.which('adt') is not None or shutil.which('AutoDockTools') is not None,
        }


# ============================================================
# Plugin Manager
# ============================================================
class PluginManager:
    def __init__(self):
        self.plugins = []
        self.plugin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plugins')

    def discover_plugins(self):
        if not os.path.isdir(self.plugin_dir):
            os.makedirs(self.plugin_dir, exist_ok=True)
            return []
        self.plugins = []
        for fname in sorted(os.listdir(self.plugin_dir)):
            if fname.endswith('.py') and not fname.startswith('_'):
                try:
                    spec = importlib.util.spec_from_file_location(fname[:-3],
                        os.path.join(self.plugin_dir, fname))
                    mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(mod)
                    if hasattr(mod, 'register_plugin'):
                        info = mod.register_plugin()
                        self.plugins.append(info)
                except Exception as e:
                    print(f"Plugin load error ({fname}): {e}")
        return self.plugins

    def get_plugin_actions(self, parent) -> List[QAction]:
        actions = []
        for p in self.plugins:
            name = p.get('name', 'Unknown Plugin')
            action = QAction(name, parent)
            if 'callback' in p:
                action.triggered.connect(p['callback'])
            actions.append(action)
        return actions


# ============================================================
# Theme System
# ============================================================
class ThemeManager:
    THEMES = {
        'Deep Scientific Navy': {
            'name': 'Deep Scientific Navy',
            'BG_DARK': '#0B1020', 'BG_PANEL': '#121A2B', 'BG_MID': '#1A2740',
            'BG_HOVER': '#1F2A45', 'BG_SEL': '#24314D', 'BORDER': '#24314D',
            'BORDER_HL': '#00D1FF', 'TEXT': '#E5E7EB', 'TEXT_MUTED': '#9CA3AF',
            'TEXT_HEADER': '#8B9DC3', 'ACCENT': '#00D1FF', 'ACCENT2': '#8B5CF6',
            'SUCCESS': '#22C55E', 'WARN': '#F59E0B', 'DANGER': '#EF4444',
            'GRADIENT': 'stop:0 #00D1FF, stop:1 #8B5CF6',
            'MATPLOTLIB': 'dark_background',
        },
        'Clinical Clean White': {
            'name': 'Clinical Clean White',
            'BG_DARK': '#FFFFFF', 'BG_PANEL': '#F6F8FB', 'BG_MID': '#EEF2F7',
            'BG_HOVER': '#E3E8F0', 'BG_SEL': '#DBEAFE', 'BORDER': '#E5E7EB',
            'BORDER_HL': '#2563EB', 'TEXT': '#111827', 'TEXT_MUTED': '#4B5563',
            'TEXT_HEADER': '#1F2937', 'ACCENT': '#2563EB', 'ACCENT2': '#10B981',
            'SUCCESS': '#10B981', 'WARN': '#F59E0B', 'DANGER': '#DC2626',
            'GRADIENT': 'stop:0 #2563EB, stop:1 #10B981',
            'MATPLOTLIB': 'default',
        },
        'Molecular Graphite': {
            'name': 'Molecular Graphite',
            'BG_DARK': '#1B1B1B', 'BG_PANEL': '#2A2A2A', 'BG_MID': '#333333',
            'BG_HOVER': '#3A3A3A', 'BG_SEL': '#404040', 'BORDER': '#3F3F3F',
            'BORDER_HL': '#4FD1C5', 'TEXT': '#E5E7EB', 'TEXT_MUTED': '#A1A1AA',
            'TEXT_HEADER': '#9CA3AF', 'ACCENT': '#4FD1C5', 'ACCENT2': '#60A5FA',
            'SUCCESS': '#34D399', 'WARN': '#FBBF24', 'DANGER': '#F87171',
            'GRADIENT': 'stop:0 #4FD1C5, stop:1 #60A5FA',
            'MATPLOTLIB': 'dark_background',
        },
        'Emerald Genomics': {
            'name': 'Emerald Genomics',
            'BG_DARK': '#0F172A', 'BG_PANEL': '#1E293B', 'BG_MID': '#111827',
            'BG_HOVER': '#1E3A2E', 'BG_SEL': '#2D4A3E', 'BORDER': '#334155',
            'BORDER_HL': '#22C55E', 'TEXT': '#F1F5F9', 'TEXT_MUTED': '#94A3B8',
            'TEXT_HEADER': '#CBD5E1', 'ACCENT': '#22C55E', 'ACCENT2': '#06B6D4',
            'SUCCESS': '#22C55E', 'WARN': '#FBBF24', 'DANGER': '#EF4444',
            'GRADIENT': 'stop:0 #22C55E, stop:1 #06B6D4',
            'MATPLOTLIB': 'dark_background',
        },
        'Monochrome Paper': {
            'name': 'Monochrome Paper',
            'BG_DARK': '#FFFFFF', 'BG_PANEL': '#FAFAFA', 'BG_MID': '#F5F5F5',
            'BG_HOVER': '#EEEEEE', 'BG_SEL': '#E0E0E0', 'BORDER': '#E5E7EB',
            'BORDER_HL': '#000000', 'TEXT': '#111111', 'TEXT_MUTED': '#555555',
            'TEXT_HEADER': '#000000', 'ACCENT': '#2563EB', 'ACCENT2': '#000000',
            'SUCCESS': '#15803D', 'WARN': '#92400E', 'DANGER': '#B91C1C',
            'GRADIENT': 'stop:0 #2563EB, stop:1 #000000',
            'MATPLOTLIB': 'default',
        },
        'Professional High Contrast': {
            'name': 'Professional High Contrast',
            'BG_DARK': '#0D1117', 'BG_PANEL': '#161B22', 'BG_MID': '#21262D',
            'BG_HOVER': '#292E36', 'BG_SEL': '#30363D', 'BORDER': '#30363D',
            'BORDER_HL': '#58A6FF', 'TEXT': '#F0F6FC', 'TEXT_MUTED': '#8B949E',
            'TEXT_HEADER': '#C9D1D9', 'ACCENT': '#58A6FF', 'ACCENT2': '#F78166',
            'SUCCESS': '#3FB950', 'WARN': '#D29922', 'DANGER': '#F85149',
            'GRADIENT': 'stop:0 #58A6FF, stop:1 #F78166',
            'MATPLOTLIB': 'dark_background',
        },
    }

    # Contrast accent variants that can be applied on top of any theme
    CONTRAST_VARIANTS = {
        'Default (Theme)': {},
        'Ocean Blue':      {'BORDER_HL': '#58A6FF', 'ACCENT': '#58A6FF', 'ACCENT2': '#79C0FF', 'GRADIENT': 'stop:0 #58A6FF, stop:1 #79C0FF'},
        'Teal Cyan':       {'BORDER_HL': '#56D4DD', 'ACCENT': '#56D4DD', 'ACCENT2': '#7EE8F0', 'GRADIENT': 'stop:0 #56D4DD, stop:1 #7EE8F0'},
        'Vibrant Purple':  {'BORDER_HL': '#BC8CFF', 'ACCENT': '#BC8CFF', 'ACCENT2': '#D2A8FF', 'GRADIENT': 'stop:0 #BC8CFF, stop:1 #D2A8FF'},
        'Sunset Orange':   {'BORDER_HL': '#F78166', 'ACCENT': '#F78166', 'ACCENT2': '#FFA28B', 'GRADIENT': 'stop:0 #F78166, stop:1 #FFA28B'},
        'Lime Green':      {'BORDER_HL': '#3FB950', 'ACCENT': '#3FB950', 'ACCENT2': '#56D364', 'GRADIENT': 'stop:0 #3FB950, stop:1 #56D364'},
        'Rose Red':        {'BORDER_HL': '#F85149', 'ACCENT': '#F85149', 'ACCENT2': '#FF7B72', 'GRADIENT': 'stop:0 #F85149, stop:1 #FF7B72'},
        'Amber Gold':      {'BORDER_HL': '#D29922', 'ACCENT': '#D29922', 'ACCENT2': '#E3B341', 'GRADIENT': 'stop:0 #D29922, stop:1 #E3B341'},
    }

    STYLESHEET_TEMPLATE = """
        QMainWindow {{ background: {BG_DARK}; }}
        QWidget {{ color: {TEXT}; }}
        QTreeWidget, QTextEdit, QTableWidget, QListWidget, QPlainTextEdit {{
            background: {BG_PANEL};
            color: {TEXT};
            border: 1px solid {BORDER};
            border-radius: 4px;
            selection-background-color: {BG_SEL};
            selection-color: {ACCENT};
        }}
        QTreeWidget::item:hover, QTableWidget::item:hover {{
            background: {BG_HOVER};
        }}
        QTreeWidget::item:selected, QTableWidget::item:selected {{
            background: {BG_SEL};
            color: {ACCENT};
        }}
        QTreeWidget::item:alternate, QTableWidget::item:alternate {{
            background: {BG_DARK};
        }}
        QHeaderView::section {{
            background: {BG_MID};
            color: {TEXT_HEADER};
            border: 1px solid {BORDER};
            padding: 4px;
            font-weight: bold;
        }}
        QGroupBox {{
            font-weight: bold;
            color: {TEXT_HEADER};
            border: 1px solid {BORDER};
            border-radius: 6px;
            margin-top: 10px;
            padding-top: 10px;
        }}
        QGroupBox::title {{
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 8px;
            color: {ACCENT};
        }}
        QPushButton {{
            border: 1px solid {BORDER};
            border-radius: 4px;
            padding: 5px 12px;
            background: {BG_MID};
            color: {TEXT};
        }}
        QPushButton:hover {{
            background: {BG_SEL};
            border-color: {BORDER_HL};
        }}
        QPushButton:pressed {{
            background: {BG_PANEL};
        }}
        QComboBox {{
            padding: 3px 8px;
            border: 1px solid {BORDER};
            border-radius: 4px;
            background: {BG_PANEL};
            color: {TEXT};
        }}
        QComboBox:hover {{
            border-color: {BORDER_HL};
        }}
        QComboBox::drop-down {{
            border: none;
            background: {BG_MID};
        }}
        QComboBox QAbstractItemView {{
            background: {BG_PANEL};
            color: {TEXT};
            selection-background-color: {BG_SEL};
            border: 1px solid {BORDER};
        }}
        QSpinBox, QDoubleSpinBox {{
            padding: 3px;
            border: 1px solid {BORDER};
            border-radius: 4px;
            background: {BG_PANEL};
            color: {TEXT};
        }}
        QSpinBox:hover, QDoubleSpinBox:hover {{
            border-color: {BORDER_HL};
        }}
        QSpinBox::up-button, QDoubleSpinBox::up-button,
        QSpinBox::down-button, QDoubleSpinBox::down-button {{
            background: {BG_MID};
            border: 1px solid {BORDER};
        }}
        QCheckBox {{ color: {TEXT}; }}
        QCheckBox::indicator {{
            border: 1px solid {BORDER};
            background: {BG_PANEL};
        }}
        QCheckBox::indicator:checked {{
            background: {ACCENT};
        }}
        QRadioButton {{ color: {TEXT}; }}
        QTabWidget::pane {{
            border: 1px solid {BORDER};
            background: {BG_DARK};
        }}
        QTabBar::tab {{
            background: {BG_MID};
            color: {TEXT_MUTED};
            padding: 6px 16px;
            border: 1px solid {BORDER};
            border-bottom: none;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }}
        QTabBar::tab:selected {{
            background: {BG_DARK};
            color: {ACCENT};
            border-bottom: 1px solid {BG_DARK};
        }}
        QTabBar::tab:hover:!selected {{
            background: {BG_SEL};
            color: {TEXT};
        }}
        QScrollBar:vertical {{
            background: {BG_DARK};
            width: 10px;
            border: none;
        }}
        QScrollBar::handle:vertical {{
            background: {BORDER};
            border-radius: 5px;
            min-height: 20px;
        }}
        QScrollBar::handle:vertical:hover {{ background: {BORDER_HL}; }}
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
        QScrollBar:horizontal {{
            background: {BG_DARK};
            height: 10px;
            border: none;
        }}
        QScrollBar::handle:horizontal {{
            background: {BORDER};
            border-radius: 5px;
            min-width: 20px;
        }}
        QScrollBar::handle:horizontal:hover {{ background: {BORDER_HL}; }}
        QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{ width: 0; }}
        QMenuBar {{
            background: {BG_DARK};
            color: {TEXT_HEADER};
            border-bottom: 1px solid {BORDER};
        }}
        QMenuBar::item:selected {{
            background: {BG_SEL};
            color: {ACCENT};
        }}
        QMenu {{
            background: {BG_MID};
            color: {TEXT};
            border: 1px solid {BORDER};
        }}
        QMenu::item:selected {{
            background: {BG_SEL};
            color: {ACCENT};
        }}
        QMenu::separator {{
            height: 1px;
            background: {BORDER};
            margin: 4px 8px;
        }}
        QToolBar {{
            background: {BG_PANEL};
            border: none;
            border-bottom: 1px solid {BORDER};
            spacing: 4px;
            padding: 2px;
        }}
        QToolBar QToolButton {{
            border: 1px solid transparent;
            border-radius: 3px;
            padding: 4px 8px;
            color: {TEXT_HEADER};
        }}
        QToolBar QToolButton:hover {{
            background: {BG_SEL};
            border-color: {BORDER_HL};
            color: {ACCENT};
        }}
        QStatusBar {{
            background: {BG_PANEL};
            border-top: 1px solid {BORDER};
            color: {TEXT_MUTED};
        }}
        QSplitter::handle {{
            background: {BORDER};
            width: 1px;
        }}
        QProgressBar {{
            border: 1px solid {BORDER};
            border-radius: 4px;
            background: {BG_PANEL};
            text-align: center;
            color: {TEXT};
        }}
        QProgressBar::chunk {{
            background: qlineargradient(x1:0, y1:0, x2:1, y2:0, {GRADIENT});
            border-radius: 3px;
        }}
        QLabel {{ color: {TEXT_HEADER}; }}
    """

    @classmethod
    def generate_stylesheet(cls, theme_name: str) -> str:
        if theme_name not in cls.THEMES:
            theme_name = 'Deep Scientific Navy'
        colors = cls.THEMES[theme_name]
        return cls.STYLESHEET_TEMPLATE.format(**colors)

    @classmethod
    def get_theme_names(cls) -> List[str]:
        return list(cls.THEMES.keys())

    @staticmethod
    def _generate_from_dict(theme_dict: dict) -> str:
        return ThemeManager.STYLESHEET_TEMPLATE.format(**theme_dict)
