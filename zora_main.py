#!/usr/bin/env python3
"""
ZORA - Sequence & Mutation Analysis Workstation
A modular platform for loading sequences, analyzing mutations,
editing sequences, and visualizing results.
"""

import sys
import os
import re
import json
import csv
import hashlib
import datetime
import tempfile
import subprocess
import shutil
import importlib.util
from collections import Counter
from typing import List, Tuple, Dict, Optional, Union

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QTabWidget, QTreeWidget, QTreeWidgetItem, QTextEdit,
    QTableWidget, QTableWidgetItem, QPushButton, QLabel, QLineEdit,
    QFileDialog, QMessageBox, QMenu, QMenuBar, QStatusBar, QToolBar,
    QInputDialog, QComboBox, QCheckBox, QGroupBox, QGridLayout, QProgressBar,
    QScrollArea, QFrame, QHeaderView, QDialog, QDialogButtonBox, QFormLayout,
    QSlider, QSpinBox, QDoubleSpinBox, QRadioButton, QButtonGroup, QSizePolicy, QListWidget, QListWidgetItem
)
from PySide6.QtCore import Qt, QSize, Signal, Slot, QRect, QMimeData, QUrl, QThread, QObject, QRunnable, QThreadPool
from PySide6.QtGui import QAction, QActionGroup, QIcon, QFont, QColor, QSyntaxHighlighter, QTextCharFormat, QPainter, QPixmap, QDragEnterEvent, QDropEvent

import numpy as np
import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# Dark theme for matplotlib
plt.style.use('dark_background')
matplotlib.rcParams.update({
    'figure.facecolor': '#0B1020',
    'axes.facecolor': '#0D1528',
    'axes.edgecolor': '#1F2A40',
    'axes.labelcolor': '#8B9DC3',
    'axes.titlecolor': '#00D1FF',
    'xtick.color': '#8B9DC3',
    'ytick.color': '#8B9DC3',
    'grid.color': '#1F2A40',
    'grid.alpha': 0.2,
    'text.color': '#E5E7EB',
    'legend.facecolor': '#121A2B',
    'legend.edgecolor': '#1F2A40',
    'legend.fancybox': False,
    'legend.framealpha': 0.85,
    'font.size': 10,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
})

# ============================================================
# Constants
# ============================================================
AMINO_ACIDS = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
    'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    '*': 'Stop Codon'
}

COMPLEMENT_TABLE = str.maketrans('ATCGNatcgn', 'TAGCNTAGCN')

CODON_TABLE = {v: k for k, v in AMINO_ACIDS.items()}  # reverse lookup


# ============================================================
# Sequence Utilities
# ============================================================
class SeqUtils:
    @staticmethod
    def validate_dna(seq: str) -> bool:
        return all(c in 'ATCGNatcgn' for c in seq.strip())

    @staticmethod
    def reverse_complement(seq: str) -> str:
        return seq.translate(COMPLEMENT_TABLE)[::-1]

    @staticmethod
    def complement(seq: str) -> str:
        return seq.translate(COMPLEMENT_TABLE)

    @staticmethod
    def gc_content(seq: str) -> float:
        if not seq:
            return 0.0
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        return (gc / len(seq)) * 100.0

    @staticmethod
    def gc_skew(seq: str) -> float:
        seq = seq.upper()
        g = seq.count('G')
        c = seq.count('C')
        if g + c == 0:
            return 0.0
        return (g - c) / (g + c)

    @staticmethod
    def nucleotide_frequency(seq: str) -> Dict[str, float]:
        if not seq:
            return {}
        seq = seq.upper()
        length = len(seq)
        freqs = {}
        for base in 'ATCG':
            freqs[base] = (seq.count(base) / length) * 100.0
        return freqs

    @staticmethod
    def translate(seq: str, frame: int = 0) -> str:
        seq = seq.upper()
        seq = seq[frame:]
        protein = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            protein.append(AMINO_ACIDS.get(codon, '?'))
        return ''.join(protein)

    @staticmethod
    def find_orfs(seq: str, min_length: int = 30) -> List[Tuple[int, int, str, str]]:
        orfs = []
        seq = seq.upper()
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                codon = seq[i:i+3]
                if codon == 'ATG':
                    start = i
                    for j in range(i, len(seq) - 2, 3):
                        codon = seq[j:j+3]
                        if AMINO_ACIDS.get(codon, 'X') == '*':
                            stop = j + 3
                            length = stop - start
                            if length >= min_length:
                                protein = SeqUtils.translate(seq[start:stop], 0)
                                orfs.append((start, stop, protein, f"+{frame+1}"))
                            i = j + 3
                            break
                    else:
                        i = len(seq)
                else:
                    i += 3
        return orfs

    @staticmethod
    def find_motifs(seq: str, pattern: str, regex: bool = False) -> List[Tuple[int, str]]:
        seq = seq.upper()
        matches = []
        if regex:
            for m in re.finditer(pattern, seq, re.IGNORECASE):
                matches.append((m.start(), m.group()))
        else:
            pattern = pattern.upper()
            start = 0
            while True:
                pos = seq.find(pattern, start)
                if pos == -1:
                    break
                matches.append((pos, pattern))
                start = pos + 1
        return matches

    @staticmethod
    def codon_usage(seq: str) -> Dict[str, Dict]:
        seq = seq.upper()
        usage = {}
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) != 3:
                continue
            aa = AMINO_ACIDS.get(codon, 'X')
            if aa not in usage:
                usage[aa] = {'codons': {}, 'total': 0}
            usage[aa]['codons'][codon] = usage[aa]['codons'].get(codon, 0) + 1
            usage[aa]['total'] += 1
        for aa in usage:
            for codon in usage[aa]['codons']:
                usage[aa]['codons'][codon] = {
                    'count': usage[aa]['codons'][codon],
                    'freq': usage[aa]['codons'][codon] / usage[aa]['total'] * 100
                }
        return usage

    @staticmethod
    def molecular_weight(seq: str) -> float:
        weights = {'A': 313.21, 'T': 304.20, 'G': 329.21, 'C': 289.18}
        seq = seq.upper()
        mw = 0.0
        for base in seq:
            mw += weights.get(base, 0)
        return round(mw, 2)

    @staticmethod
    def melting_temp(seq: str) -> float:
        seq = seq.upper()
        if len(seq) < 14:
            return (seq.count('A') + seq.count('T')) * 2 + (seq.count('G') + seq.count('C')) * 4
        gc = seq.count('G') + seq.count('C')
        return 64.9 + 41.0 * (gc - 16.4) / len(seq)

    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> int:
        min_len = min(len(seq1), len(seq2))
        return sum(1 for i in range(min_len) if seq1[i].upper() != seq2[i].upper())


# ============================================================
# Data Models
# ============================================================
class SequenceRecord:
    def __init__(self, name: str = "", sequence: str = "", description: str = ""):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description
        self.id = hashlib.md5(sequence.encode()).hexdigest()[:8]

    def __len__(self):
        return len(self.sequence)

    def stats(self) -> Dict:
        return {
            'length': len(self.sequence),
            'gc_content': SeqUtils.gc_content(self.sequence),
            'gc_skew': SeqUtils.gc_skew(self.sequence),
            'molecular_weight': SeqUtils.molecular_weight(self.sequence),
            'melting_temp': SeqUtils.melting_temp(self.sequence),
            'nucleotide_freq': SeqUtils.nucleotide_frequency(self.sequence)
        }


class Mutation:
    def __init__(self, position: int, ref: str, alt: str, mut_type: str = "SNP"):
        self.position = position
        self.ref = ref
        self.alt = alt
        self.mut_type = mut_type

    def __str__(self):
        return f"{self.mut_type}\t{self.position}\t{self.ref}→{self.alt}"


class MutationResult:
    def __init__(self, position: int, ref_base: str, alt_base: str,
                 ref_codon: str = "", alt_codon: str = "",
                 ref_aa: str = "", alt_aa: str = "",
                 classification: str = ""):
        self.position = position
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.ref_aa = ref_aa
        self.alt_aa = alt_aa
        self.classification = classification

    def to_dict(self) -> Dict:
        return {
            'position': self.position,
            'ref_base': self.ref_base,
            'alt_base': self.alt_base,
            'ref_codon': self.ref_codon,
            'alt_codon': self.alt_codon,
            'ref_aa': self.ref_aa,
            'alt_aa': self.alt_aa,
            'classification': self.classification
        }


# ============================================================
# Mutation Analysis Engine
# ============================================================
class MutationEngine:
    @staticmethod
    def detect_mutations(seq1: str, seq2: str) -> List[MutationResult]:
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        results = []
        i, j = 0, 0
        len1, len2 = len(seq1), len(seq2)
        while i < len1 and j < len2:
            if seq1[i] != seq2[j]:
                if seq1[i] == '-' or seq2[j] == '-':
                    if seq1[i] == '-':
                        results.append(MutationResult(i, '-', seq2[j], '', '', '', '', 'Insertion'))
                        j += 1
                    else:
                        results.append(MutationResult(i, seq1[i], '-', '', '', '', '', 'Deletion'))
                        i += 1
                else:
                    # SNP: classify by codon impact
                    res = MutationEngine._classify_snp(seq1, seq2, i)
                    results.append(res)
                    i += 1
                    j += 1
            else:
                i += 1
                j += 1
        while i < len1:
            results.append(MutationResult(i, seq1[i], '-', '', '', '', '', 'Deletion'))
            i += 1
        while j < len2:
            results.append(MutationResult(i, '-', seq2[j], '', '', '', '', 'Insertion'))
            j += 1
        return results

    @staticmethod
    def _classify_snp(seq1: str, seq2: str, pos: int) -> MutationResult:
        ref_base = seq1[pos].upper()
        alt_base = seq2[pos].upper()
        # Find codon position for this base
        start = (pos // 3) * 3
        if start + 2 >= len(seq1) or start + 2 >= len(seq2):
            return MutationResult(pos, ref_base, alt_base, classification="Intergenic")

        ref_codon = seq1[start:start+3].upper()
        alt_codon = seq2[start:start+3].upper()
        ref_aa = AMINO_ACIDS.get(ref_codon, '?')
        alt_aa = AMINO_ACIDS.get(alt_codon, '?')

        if ref_aa == alt_aa:
            classification = 'Silent (Synonymous)'
        elif ref_aa == '*':
            classification = 'Stop Loss'
        elif alt_aa == '*':
            classification = 'Nonsense'
        else:
            classification = 'Missense (Non-synonymous)'

        return MutationResult(pos, ref_base, alt_base, ref_codon, alt_codon, ref_aa, alt_aa, classification)

    @staticmethod
    def simulate_mutation(seq: str, position: int, alt_base: str) -> str:
        seq_list = list(seq.upper())
        if position < len(seq):
            seq_list[position] = alt_base.upper()
        return ''.join(seq_list)

    @staticmethod
    def simulate_insertion(seq: str, position: int, insertion: str) -> str:
        return seq[:position] + insertion.upper() + seq[position:]

    @staticmethod
    def simulate_deletion(seq: str, start: int, end: int) -> str:
        return seq[:start] + seq[end:]

    @staticmethod
    def batch_generate_mutations(seq: str, num_mutations: int = 100,
                                 snp_rate: float = 0.7, ins_rate: float = 0.15,
                                 del_rate: float = 0.15) -> List[Tuple[str, List[MutationResult]]]:
        import random
        variants = []
        for _ in range(num_mutations):
            mutated = list(seq.upper())
            r = random.random()
            if r < snp_rate:
                pos = random.randint(0, len(seq) - 1)
                bases = ['A', 'T', 'C', 'G']
                bases.remove(mutated[pos])
                mutated[pos] = random.choice(bases)
            elif r < snp_rate + ins_rate:
                pos = random.randint(0, len(seq))
                base = random.choice(['A', 'T', 'C', 'G'])
                mutated.insert(pos, base)
            else:
                if len(mutated) > 1:
                    pos = random.randint(0, len(mutated) - 1)
                    mutated.pop(pos)
            variants.append(''.join(mutated))
        return variants


# ============================================================
# Alignment Engine
# ============================================================
class AlignmentEngine:
    @staticmethod
    def needleman_wunsch(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str, int]:
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            dp[i][0] = gap * i
        for j in range(m + 1):
            dp[0][j] = gap * j
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                up = dp[i-1][j] + gap
                left = dp[i][j-1] + gap
                dp[i][j] = max(diag, up, left)
        # Traceback
        i, j = n, m
        aln1, aln2 = '', ''
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                aln1 = seq1[i-1] + aln1
                aln2 = seq2[j-1] + aln2
                i -= 1
                j -= 1
            elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
                aln1 = seq1[i-1] + aln1
                aln2 = '-' + aln2
                i -= 1
            else:
                aln1 = '-' + aln1
                aln2 = seq2[j-1] + aln2
                j -= 1
        return aln1, aln2, dp[n][m]

    @staticmethod
    def smith_waterman(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str, int]:
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        max_score = 0
        max_pos = (0, 0)
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                up = dp[i-1][j] + gap
                left = dp[i][j-1] + gap
                dp[i][j] = max(0, diag, up, left)
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_pos = (i, j)
        i, j = max_pos
        aln1, aln2 = '', ''
        while i > 0 and j > 0 and dp[i][j] > 0:
            if dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                aln1 = seq1[i-1] + aln1
                aln2 = seq2[j-1] + aln2
                i -= 1
                j -= 1
            elif dp[i][j] == dp[i-1][j] + gap:
                aln1 = seq1[i-1] + aln1
                aln2 = '-' + aln2
                i -= 1
            else:
                aln1 = '-' + aln1
                aln2 = seq2[j-1] + aln2
                j -= 1
        return aln1, aln2, max_score


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
        self.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding))
        FigureCanvas.updateGeometry(self)

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
        self.setStyleSheet("""
            QFrame#dropArea { border: 2px dashed #1F2A40; border-radius: 8px; background: #0D1528; }
            QFrame#dropArea:hover { border-color: #00D1FF; background: #0F1A35; }
        """)
        self.setObjectName("dropArea")
        self.label.setStyleSheet("color: #6B7280; font-size: 13px;")

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            self.setStyleSheet("""
                QFrame#dropArea { border: 2px dashed #00D1FF; border-radius: 8px; background: #0F1A35; }
            """)

    def dragLeaveEvent(self, event):
        self.setStyleSheet("""
            QFrame#dropArea { border: 2px dashed #1F2A40; border-radius: 8px; background: #0D1528; }
        """)

    def dropEvent(self, event: QDropEvent):
        records = []
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path:
                records.extend(FileParser.parse_file(path))
        if records:
            self.sequences_dropped.emit(records)
        self.setStyleSheet("""
            QFrame#dropArea { border: 2px dashed #00D1FF; border-radius: 8px; background: #0F1A35; }
        """)
        # Reset after brief delay
        from PySide6.QtCore import QTimer
        QTimer.singleShot(2000, lambda: self.setStyleSheet("""
            QFrame#dropArea { border: 2px dashed #1F2A40; border-radius: 8px; background: #0D1528; }
        """))

    def mousePressEvent(self, event):
        files, _ = QFileDialog.getOpenFileName(self, "Open Sequence File",
                                                "", "Sequence Files (*.fasta *.fa *.fastq *.gb *.genbank *.txt *.csv)")
        if files:
            records = FileParser.parse_file(files)
            if records:
                self.sequences_dropped.emit(records)


# ============================================================
# File Parser
# ============================================================
class FileParser:
    @staticmethod
    def parse_file(path: str) -> List[SequenceRecord]:
        ext = os.path.splitext(path)[1].lower()
        if ext in ('.fasta', '.fa', '.fna'):
            return FileParser.parse_fasta(path)
        elif ext in ('.fastq', '.fq'):
            return FileParser.parse_fastq(path)
        elif ext in ('.gb', '.genbank'):
            return FileParser.parse_genbank(path)
        elif ext == '.csv':
            return FileParser.parse_csv(path)
        elif ext == '.txt':
            return FileParser.parse_txt(path)
        return []

    @staticmethod
    def parse_fasta(path: str) -> List[SequenceRecord]:
        records = []
        with open(path, 'r') as f:
            name, seq = '', ''
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if name and seq:
                        records.append(SequenceRecord(name, seq))
                    parts = line[1:].split(None, 1)
                    name = parts[0] if parts else ''
                    desc = parts[1] if len(parts) > 1 else ''
                    seq = ''
                else:
                    seq += line.replace(' ', '').replace('\t', '')
            if name and seq:
                records.append(SequenceRecord(name, seq))
        return records

    @staticmethod
    def parse_fastq(path: str) -> List[SequenceRecord]:
        records = []
        with open(path, 'r') as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith('@'):
                name = lines[i][1:].strip()
                seq = lines[i+1].strip() if i+1 < len(lines) else ''
                records.append(SequenceRecord(name, seq))
                i += 4
            else:
                i += 1
        return records

    @staticmethod
    def parse_genbank(path: str) -> List[SequenceRecord]:
        with open(path, 'r') as f:
            content = f.read()
        return FileParser.parse_genbank_from_text(content)

    @staticmethod
    def parse_genbank_from_text(content: str) -> List[SequenceRecord]:
        records = []
        blocks = re.split(r'^//\n?', content, flags=re.MULTILINE)
        for block in blocks:
            if not block.strip():
                continue
            name = ''
            seq = ''
            locus = re.search(r'^LOCUS\s+(\S+)', block, re.MULTILINE)
            if locus:
                name = locus.group(1)
            in_seq = False
            for line in block.split('\n'):
                if line.startswith('ORIGIN'):
                    in_seq = True
                    continue
                if in_seq:
                    seq += re.sub(r'[^a-zA-Z]', '', line)
            if name and seq:
                records.append(SequenceRecord(name, seq.upper()))
        return records

    @staticmethod
    def parse_csv(path: str) -> List[SequenceRecord]:
        records = []
        with open(path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', row.get('id', 'Unknown'))
                seq = row.get('sequence', row.get('seq', ''))
                if seq:
                    records.append(SequenceRecord(name, seq))
        return records

    @staticmethod
    def parse_fasta_from_text(text: str) -> List[SequenceRecord]:
        records = []
        name, seq = '', ''
        for line in text.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name and seq:
                    records.append(SequenceRecord(name, seq))
                parts = line[1:].split(None, 1)
                name = parts[0] if parts else ''
                seq = ''
            else:
                seq += line.replace(' ', '').replace('\t', '')
        if name and seq:
            records.append(SequenceRecord(name, seq))
        return records

    @staticmethod
    def parse_txt(path: str) -> List[SequenceRecord]:
        records = []
        with open(path, 'r') as f:
            seq = f.read().strip().replace(' ', '').replace('\n', '').replace('\r', '')
        if seq:
            name = os.path.splitext(os.path.basename(path))[0]
            records.append(SequenceRecord(name, seq))
        return records


# ============================================================
# Main Window
# ============================================================
class ZORAMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.sequences: List[SequenceRecord] = []
        self.mutation_results: List[MutationResult] = []
        self.alignment_results: List[Tuple[str, str, int]] = []
        self.current_figure = None
        self.plot_windows = []
        self.project_manager = ProjectManager()
        self.current_theme = 'Deep Scientific Navy'
        self.current_contrast = 'Default (Theme)'
        self._graph_grid_visible = True
        self._gc_cache = {}  # (seq_hash, window, step) -> (positions, gc_values)
        self.docking_results = []
        self._dock_worker = None
        self._last_docking_result = None
        self._init_ui()
        self._build_menu()
        self._build_toolbar()
        self._connect_signals()
        self.statusBar().showMessage("Ready")
        self._check_recovery()

    def closeEvent(self, event):
        if self.project_manager.is_open():
            self.project_manager.save(
                sequences=self.sequences,
                mutation_results=self.mutation_results,
                alignment_results=self.alignment_results,
                analysis_text=self.analysis_output.toPlainText(),
            )
            self._save_recovery()
        event.accept()

    def _btn(self, text: str, bg: str, fg: str = '#0B1020', bold: bool = True) -> QPushButton:
        btn = QPushButton(text)
        weight = 'bold' if bold else 'normal'
        btn.setStyleSheet(f"background: {bg}; color: {fg}; padding: 6px 12px; border: none; border-radius: 4px; font-weight: {weight};")
        return btn

    def _ensure_canvas(self):
        if not hasattr(self, 'graph_canvas') or self.graph_canvas is None:
            self.graph_canvas = MplCanvas(self)
            self.graph_toolbar = NavigationToolbar(self.graph_canvas, self)
            graphs_idx = 4
            tab = self.sequence_tabs.widget(graphs_idx)
            if tab:
                layout = tab.layout()
                # Remove placeholder
                if hasattr(self, '_canvas_placeholder') and self._canvas_placeholder is not None:
                    layout.removeWidget(self._canvas_placeholder)
                    self._canvas_placeholder.deleteLater()
                    self._canvas_placeholder = None
                layout.addWidget(self.graph_toolbar)
                layout.addWidget(self.graph_canvas)

    def _apply_mpl_rcparams(self, theme: dict, mpl_style: str):
        if mpl_style == 'default':
            matplotlib.rcParams.update({
                'figure.facecolor': '#FFFFFF', 'axes.facecolor': '#F6F8FB',
                'axes.edgecolor': '#D1D5DB', 'axes.labelcolor': '#1F2937',
                'axes.titlecolor': theme.get('ACCENT', '#2563EB'),
                'xtick.color': '#4B5563', 'ytick.color': '#4B5563',
                'grid.color': '#E5E7EB', 'text.color': '#111827',
                'legend.facecolor': '#FFFFFF', 'legend.edgecolor': '#D1D5DB',
            })
        else:
            matplotlib.rcParams.update({
                'figure.facecolor': theme.get('BG_DARK', '#0B1020'),
                'axes.facecolor': theme.get('BG_PANEL', '#0D1528'),
                'axes.edgecolor': theme.get('BORDER', '#1F2A40'),
                'axes.labelcolor': theme.get('TEXT_HEADER', '#8B9DC3'),
                'axes.titlecolor': theme.get('ACCENT', '#00D1FF'),
                'xtick.color': theme.get('TEXT_MUTED', '#9CA3AF'),
                'ytick.color': theme.get('TEXT_MUTED', '#9CA3AF'),
                'grid.color': theme.get('BORDER', '#1F2A40'),
                'text.color': theme.get('TEXT', '#E5E7EB'),
                'legend.facecolor': theme.get('BG_PANEL', '#121A2B'),
                'legend.edgecolor': theme.get('BORDER', '#1F2A40'),
            })

    def _init_ui(self):
        self.setWindowTitle("ZORA - Sequence & Mutation Analysis Workstation")
        self.setGeometry(100, 100, 1400, 900)

        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(4, 4, 4, 4)

        # Main splitter
        self.main_splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(self.main_splitter)

        # Left: Project Explorer
        self.project_tree = QTreeWidget()
        self.project_tree.setHeaderLabel("Project Explorer")
        self.project_tree.setMinimumWidth(200)
        self.project_tree.setMaximumWidth(350)
        self.project_tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.project_tree.customContextMenuRequested.connect(self._project_context_menu)

        # Center: Sequence viewer + analysis tabs
        center_panel = QSplitter(Qt.Vertical)

        # Top center: Sequence viewer
        self.sequence_tabs = QTabWidget()
        # Sequence viewer with selector
        self.seq_viewer_widget = QWidget()
        seq_viewer_layout = QVBoxLayout(self.seq_viewer_widget)
        seq_viewer_layout.setContentsMargins(0, 0, 0, 0)
        seq_selector_layout = QHBoxLayout()
        seq_selector_layout.addWidget(QLabel("Show:"))
        self.seq_viewer_combo = QComboBox()
        self.seq_viewer_combo.currentIndexChanged.connect(self._on_seq_viewer_select)
        seq_selector_layout.addWidget(self.seq_viewer_combo)
        seq_selector_layout.addStretch()
        seq_viewer_layout.addLayout(seq_selector_layout)
        self.seq_viewer = QTextEdit()
        self.seq_viewer.setReadOnly(True)
        self.seq_viewer.setFont(QFont("Courier New", 11))
        self.seq_viewer.setLineWrapMode(QTextEdit.NoWrap)
        self.seq_highlighter = SequenceHighlighter(self.seq_viewer.document())
        seq_viewer_layout.addWidget(self.seq_viewer)
        self.sequence_tabs.addTab(self.seq_viewer_widget, "Sequence Viewer")

        # Analysis tab
        self.analysis_tab = QWidget()
        self._build_analysis_tab()
        self.sequence_tabs.addTab(self.analysis_tab, "Analysis")

        # Mutation tab
        self.mutation_tab = QWidget()
        self._build_mutation_tab()
        self.sequence_tabs.addTab(self.mutation_tab, "Mutations")

        # Alignment tab
        self.alignment_tab = QWidget()
        self._build_alignment_tab()
        self.sequence_tabs.addTab(self.alignment_tab, "Alignment")

        # Graphs tab
        self.graphs_tab = QWidget()
        self._build_graphs_tab()
        self.sequence_tabs.addTab(self.graphs_tab, "Graphs")

        # Interactive viewer tab
        self.interactive_viewer = QWidget()
        self._build_interactive_viewer()
        self.sequence_tabs.addTab(self.interactive_viewer, "Interactive Viewer")

        # CRISPR tab
        self.crispr_tab = QWidget()
        self._build_crispr_tab()
        self.sequence_tabs.addTab(self.crispr_tab, "CRISPR")

        # AI tab
        self.ai_tab = QWidget()
        self._build_ai_tab()
        self.sequence_tabs.addTab(self.ai_tab, "AI Analysis")

        # Integrations tab
        self.integrations_tab = QWidget()
        self._build_integrations_tab()
        self.sequence_tabs.addTab(self.integrations_tab, "Integrations")

        # Docking tab
        self.docking_tab = QWidget()
        self._build_docking_tab()
        self.sequence_tabs.addTab(self.docking_tab, "Docking")

        center_panel.addWidget(self.sequence_tabs)

        # Drop area
        self.drop_area = DropArea()
        center_panel.addWidget(self.drop_area)
        center_panel.setSizes([600, 100])

        # Right panel: mutation table + annotation
        right_panel = QSplitter(Qt.Vertical)

        self.mutation_table = QTableWidget()
        self.mutation_table.setColumnCount(8)
        self.mutation_table.setHorizontalHeaderLabels(
            ["Pos", "Ref", "Alt", "Ref Codon", "Alt Codon", "Ref AA", "Alt AA", "Classification"]
        )
        self.mutation_table.horizontalHeader().setStretchLastSection(True)
        self.mutation_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.mutation_table.setAlternatingRowColors(True)

        # Annotation panel
        self.annotation_text = QTextEdit()
        self.annotation_text.setReadOnly(True)
        self.annotation_text.setPlaceholderText("Sequence annotations will appear here...")
        self.annotation_text.setMaximumHeight(150)

        right_panel.addWidget(self.mutation_table)
        right_panel.addWidget(self.annotation_text)

        self.main_splitter.addWidget(self.project_tree)
        self.main_splitter.addWidget(center_panel)
        self.main_splitter.addWidget(right_panel)
        self.main_splitter.setSizes([250, 700, 350])

        # Sequence info label
        self.info_label = QLabel("No sequence loaded")
        self.statusBar().addPermanentWidget(self.info_label)

    def _build_analysis_tab(self):
        layout = QVBoxLayout(self.analysis_tab)
        top = QHBoxLayout()
        top.addWidget(QLabel("Sequence:"))
        self.analysis_seq_combo = QComboBox()
        self.analysis_seq_combo.setMinimumWidth(250)
        top.addWidget(self.analysis_seq_combo)
        self.analyze_btn = self._btn("Analyze", '#00D1FF')
        top.addWidget(self.analyze_btn)
        top.addStretch()
        layout.addLayout(top)

        # Send-to-tools bar (hidden until protein seq available)
        self.send_tools_bar = QFrame()
        send_layout = QHBoxLayout(self.send_tools_bar)
        send_layout.setContentsMargins(0, 0, 0, 0)
        send_layout.addWidget(QLabel("Send translation to:"))
        self.send_pymol_btn = self._btn("PyMOL", '#8B5CF6', 'white')
        send_layout.addWidget(self.send_pymol_btn)
        self.send_chimerax_btn = self._btn("ChimeraX", '#22C55E')
        send_layout.addWidget(self.send_chimerax_btn)
        self.send_snapgene_btn = self._btn("SnapGene", '#00D1FF')
        send_layout.addWidget(self.send_snapgene_btn)
        send_layout.addStretch()
        self.send_tools_bar.setVisible(False)
        layout.addWidget(self.send_tools_bar)

        self.analysis_output = QTextEdit()
        self.analysis_output.setReadOnly(True)
        self.analysis_output.setFont(QFont("Courier New", 10))
        layout.addWidget(self.analysis_output)

    def _build_mutation_tab(self):
        layout = QVBoxLayout(self.mutation_tab)

        # Sequence selection
        sel_group = QGroupBox("Sequence Selection")
        sel_layout = QGridLayout(sel_group)
        sel_layout.addWidget(QLabel("Reference:"), 0, 0)
        self.ref_seq_combo = QComboBox()
        sel_layout.addWidget(self.ref_seq_combo, 0, 1)
        sel_layout.addWidget(QLabel("Mutated:"), 1, 0)
        self.mut_seq_combo = QComboBox()
        sel_layout.addWidget(self.mut_seq_combo, 1, 1)
        self.detect_btn = self._btn("Detect Mutations", '#00D1FF')
        sel_layout.addWidget(self.detect_btn, 2, 0, 1, 2)
        layout.addWidget(sel_group)

        # Mutation simulation
        sim_group = QGroupBox("Mutation Simulation")
        sim_layout = QGridLayout(sim_group)
        sim_layout.addWidget(QLabel("Position:"), 0, 0)
        self.sim_pos_spin = QSpinBox()
        self.sim_pos_spin.setMinimum(0)
        self.sim_pos_spin.setMaximum(100000)
        sim_layout.addWidget(self.sim_pos_spin, 0, 1)
        sim_layout.addWidget(QLabel("Change to:"), 0, 2)
        self.sim_base_combo = QComboBox()
        self.sim_base_combo.addItems(['A', 'T', 'C', 'G'])
        sim_layout.addWidget(self.sim_base_combo, 0, 3)
        self.simulate_btn = self._btn("Simulate Mutation", '#8B5CF6', 'white')
        sim_layout.addWidget(self.simulate_btn, 0, 4)

        sim_layout.addWidget(QLabel("Insert at pos:"), 1, 0)
        self.ins_pos_spin = QSpinBox()
        self.ins_pos_spin.setMinimum(0)
        self.ins_pos_spin.setMaximum(100000)
        sim_layout.addWidget(self.ins_pos_spin, 1, 1)
        self.ins_seq_edit = QLineEdit()
        self.ins_seq_edit.setPlaceholderText("Insertion seq")
        sim_layout.addWidget(self.ins_seq_edit, 1, 2)
        self.insert_btn = self._btn("Insert", '#8B5CF6', 'white')
        sim_layout.addWidget(self.insert_btn, 1, 3)

        sim_layout.addWidget(QLabel("Delete from:"), 2, 0)
        self.del_start_spin = QSpinBox()
        self.del_start_spin.setMinimum(0)
        self.del_start_spin.setMaximum(100000)
        sim_layout.addWidget(self.del_start_spin, 2, 1)
        sim_layout.addWidget(QLabel("to:"), 2, 2)
        self.del_end_spin = QSpinBox()
        self.del_end_spin.setMinimum(0)
        self.del_end_spin.setMaximum(100000)
        sim_layout.addWidget(self.del_end_spin, 2, 3)
        self.delete_btn = self._btn("Delete", '#EF4444', 'white')
        sim_layout.addWidget(self.delete_btn, 2, 4)
        layout.addWidget(sim_group)

        # Batch mutation
        batch_group = QGroupBox("Batch Mutation Generator")
        batch_layout = QHBoxLayout(batch_group)
        batch_layout.addWidget(QLabel("Count:"))
        self.batch_count_spin = QSpinBox()
        self.batch_count_spin.setRange(1, 100000)
        self.batch_count_spin.setValue(100)
        batch_layout.addWidget(self.batch_count_spin)
        self.batch_btn = self._btn("Generate Variants", '#00D1FF')
        batch_layout.addWidget(self.batch_btn)
        self.batch_target_combo = QComboBox()
        batch_layout.addWidget(QLabel("Target:"))
        batch_layout.addWidget(self.batch_target_combo)
        batch_layout.addStretch()
        layout.addWidget(batch_group)

        # Mutation output
        mut_view_layout = QHBoxLayout()
        self.view_mut_in_seq_btn = QPushButton("View Mutations in Sequence Viewer")
        self.view_mut_in_seq_btn.setStyleSheet("background: #22C55E; color: #0B1020; font-weight: bold; padding: 4px; border: none; border-radius: 4px;")
        mut_view_layout.addWidget(self.view_mut_in_seq_btn)
        mut_view_layout.addStretch()
        layout.addLayout(mut_view_layout)

        self.mutation_output = QTextEdit()
        self.mutation_output.setReadOnly(True)
        self.mutation_output.setFont(QFont("Courier New", 10))
        layout.addWidget(self.mutation_output)

    def _build_alignment_tab(self):
        layout = QVBoxLayout(self.alignment_tab)
        sel_group = QGroupBox("Select Sequences")
        sel_layout = QGridLayout(sel_group)
        sel_layout.addWidget(QLabel("Sequence 1:"), 0, 0)
        self.aln_seq1_combo = QComboBox()
        sel_layout.addWidget(self.aln_seq1_combo, 0, 1)
        sel_layout.addWidget(QLabel("Sequence 2:"), 1, 0)
        self.aln_seq2_combo = QComboBox()
        sel_layout.addWidget(self.aln_seq2_combo, 1, 1)

        self.aln_nw_btn = QPushButton("Needleman-Wunsch (Global)")
        self.aln_nw_btn.setStyleSheet("background: #00D1FF; color: #0B1020; padding: 6px; border: none; border-radius: 4px;")
        sel_layout.addWidget(self.aln_nw_btn, 0, 2)
        self.aln_sw_btn = QPushButton("Smith-Waterman (Local)")
        self.aln_sw_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        sel_layout.addWidget(self.aln_sw_btn, 1, 2)
        layout.addWidget(sel_group)

        self.aln_output = QTextEdit()
        self.aln_output.setReadOnly(True)
        self.aln_output.setFont(QFont("Courier New", 10))
        layout.addWidget(self.aln_output)

    def _build_graphs_tab(self):
        layout = QVBoxLayout(self.graphs_tab)
        layout.setSpacing(6)

        # Sequence selection for comparison
        sel_layout = QHBoxLayout()
        sel_layout.addWidget(QLabel("Select sequences to compare:"))
        sel_layout.addStretch()
        layout.addLayout(sel_layout)
        self.graph_seq_list = QListWidget()
        self.graph_seq_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.graph_seq_list.setMaximumHeight(80)
        self.graph_seq_list.setAlternatingRowColors(True)
        layout.addWidget(self.graph_seq_list)

        btn_layout = QHBoxLayout()

        self.gc_plot_btn = self._btn("GC Content", '#22C55E')
        btn_layout.addWidget(self.gc_plot_btn)

        self.mut_dist_btn = self._btn("Mutation Dist.", '#EF4444', 'white')
        btn_layout.addWidget(self.mut_dist_btn)

        self.nuc_freq_btn = self._btn("Nucleotide Freq.", '#00D1FF')
        btn_layout.addWidget(self.nuc_freq_btn)

        self.codon_usage_btn = self._btn("Codon Usage", '#8B5CF6', 'white')
        btn_layout.addWidget(self.codon_usage_btn)

        self.mut_pie_btn = self._btn("Mutation Pie", '#FFB74D')
        btn_layout.addWidget(self.mut_pie_btn)

        self.heatmap_btn = self._btn("SNP Heatmap", '#F472B6')
        btn_layout.addWidget(self.heatmap_btn)

        self.conservation_btn = self._btn("Conservation", '#22C55E')
        btn_layout.addWidget(self.conservation_btn)

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        # Controls row: GC window + export + grid toggle
        controls_row = QHBoxLayout()
        controls_row.addWidget(QLabel("GC Window:"))
        self.gc_window_spin = QSpinBox()
        self.gc_window_spin.setRange(10, 5000)
        self.gc_window_spin.setValue(100)
        self.gc_window_spin.setSingleStep(10)
        controls_row.addWidget(self.gc_window_spin)
        controls_row.addWidget(QLabel("Step:"))
        self.gc_step_spin = QSpinBox()
        self.gc_step_spin.setRange(1, 1000)
        self.gc_step_spin.setValue(10)
        controls_row.addWidget(self.gc_step_spin)
        controls_row.addStretch()

        self.toggle_grid_btn = self._btn("Toggle Grid", '#6B7280', 'white')
        controls_row.addWidget(self.toggle_grid_btn)

        self.export_svg_btn = self._btn("Export SVG", '#6B7280', 'white')
        controls_row.addWidget(self.export_svg_btn)

        self.export_png_btn = self._btn("Export PNG", '#6B7280', 'white')
        controls_row.addWidget(self.export_png_btn)
        layout.addLayout(controls_row)

        # Placeholder for lazy MplCanvas
        self.graph_canvas = None
        self.graph_toolbar = None
        self._canvas_placeholder = QLabel("Graphs will appear here after clicking a plot button")
        self._canvas_placeholder.setAlignment(Qt.AlignCenter)
        self._canvas_placeholder.setStyleSheet("color: #6B7280; padding: 40px; font-size: 14px;")
        layout.addWidget(self._canvas_placeholder)

    def _build_interactive_viewer(self):
        layout = QVBoxLayout(self.interactive_viewer)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Position:"))
        self.viewer_pos_spin = QSpinBox()
        self.viewer_pos_spin.setMinimum(0)
        self.viewer_pos_spin.setMaximum(1000000)
        controls.addWidget(self.viewer_pos_spin)
        self.viewer_goto_btn = QPushButton("Go")
        controls.addWidget(self.viewer_goto_btn)
        controls.addStretch()
        layout.addLayout(controls)

        self.viewer_seq_display = QTextEdit()
        self.viewer_seq_display.setReadOnly(True)
        self.viewer_seq_display.setFont(QFont("Courier New", 14))
        self.viewer_seq_display.setLineWrapMode(QTextEdit.NoWrap)
        self.viewer_seq_display.setMinimumHeight(100)
        # Highlighter
        self.viewer_highlighter = SequenceHighlighter(self.viewer_seq_display.document())
        layout.addWidget(self.viewer_seq_display)

        # Info panel for selected base
        self.viewer_info = QTextEdit()
        self.viewer_info.setReadOnly(True)
        self.viewer_info.setMaximumHeight(120)
        self.viewer_info.setPlaceholderText("Click on a base above to see codon, amino acid, and mutation info")
        layout.addWidget(self.viewer_info)

    def _build_crispr_tab(self):
        layout = QVBoxLayout(self.crispr_tab)
        layout.addWidget(QLabel("CRISPR Guide RNA Design"))
        layout.addWidget(QLabel("PAM: NGG (SpCas9)"))
        self.crispr_seq_combo = QComboBox()
        layout.addWidget(QLabel("Target sequence:"))
        layout.addWidget(self.crispr_seq_combo)
        self.crispr_design_btn = QPushButton("Design sgRNAs")
        self.crispr_design_btn.setStyleSheet("background: #00D1FF; color: #0B1020; font-weight: bold; padding: 6px; border: none; border-radius: 4px;")
        layout.addWidget(self.crispr_design_btn)
        self.crispr_output = QTextEdit()
        self.crispr_output.setReadOnly(True)
        self.crispr_output.setFont(QFont("Courier New", 10))
        layout.addWidget(self.crispr_output)

    def _build_ai_tab(self):
        layout = QVBoxLayout(self.ai_tab)
        self.ai_target_combo = QComboBox()
        layout.addWidget(QLabel("Select target sequence:"))
        layout.addWidget(self.ai_target_combo)
        self.ai_analyze_btn = QPushButton("Run AI Pathogenicity Prediction (Placeholder)")
        self.ai_analyze_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        layout.addWidget(self.ai_analyze_btn)
        self.ai_output = QTextEdit()
        self.ai_output.setReadOnly(True)
        self.ai_output.setPlaceholderText("AI analysis results will appear here (requires trained model)")
        layout.addWidget(self.ai_output)

    def _build_integrations_tab(self):
        layout = QVBoxLayout(self.integrations_tab)
        layout.setSpacing(8)

        # Database Fetch Section
        db_group = QGroupBox("Database Fetch")
        db_layout = QVBoxLayout(db_group)
        db_form = QHBoxLayout()
        db_form.addWidget(QLabel("Accession/ID:"))
        self.db_accession_edit = QLineEdit()
        self.db_accession_edit.setPlaceholderText("e.g. NM_007294, NP_009225, BRCA1")
        db_form.addWidget(self.db_accession_edit)
        self.db_fetch_btn = QPushButton("Fetch from NCBI")
        self.db_fetch_btn.setStyleSheet("background: #00D1FF; color: #0B1020; font-weight: bold; padding: 6px; border: none; border-radius: 4px;")
        db_form.addWidget(self.db_fetch_btn)
        self.db_search_btn = QPushButton("Search")
        self.db_search_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        db_form.addWidget(self.db_search_btn)
        db_layout.addLayout(db_form)
        self.db_output = QTextEdit()
        self.db_output.setReadOnly(True)
        self.db_output.setMaximumHeight(120)
        self.db_output.setPlaceholderText("Database fetch results...")
        db_layout.addWidget(self.db_output)
        layout.addWidget(db_group)

        # PDB Import
        pdb_group = QGroupBox("PDB Import")
        pdb_layout = QHBoxLayout(pdb_group)
        self.pdb_import_btn = QPushButton("Import PDB File")
        self.pdb_import_btn.setStyleSheet("background: #22C55E; color: #0B1020; font-weight: bold; padding: 6px; border: none; border-radius: 4px;")
        pdb_layout.addWidget(self.pdb_import_btn)
        self.pdb_fetch_btn = QPushButton("Fetch PDB from RCSB")
        self.pdb_fetch_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        pdb_layout.addWidget(self.pdb_fetch_btn)
        pdb_layout.addStretch()
        layout.addWidget(pdb_group)

        # PDB Browser (local project files)
        pdb_browse_group = QGroupBox("Project PDB Files")
        pdb_browse_layout = QVBoxLayout(pdb_browse_group)
        pdb_top = QHBoxLayout()
        self.pdb_refresh_btn = QPushButton("Refresh")
        self.pdb_refresh_btn.setStyleSheet("background: #00D1FF; color: #0B1020; padding: 4px 10px; border: none; border-radius: 4px;")
        pdb_top.addWidget(self.pdb_refresh_btn)
        self.pdb_open_pymol_btn = QPushButton("Open in PyMOL")
        self.pdb_open_pymol_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 4px 10px; border: none; border-radius: 4px;")
        pdb_top.addWidget(self.pdb_open_pymol_btn)
        self.pdb_open_chimerax_btn = QPushButton("Open in ChimeraX")
        self.pdb_open_chimerax_btn.setStyleSheet("background: #22C55E; color: #0B1020; padding: 4px 10px; border: none; border-radius: 4px;")
        pdb_top.addWidget(self.pdb_open_chimerax_btn)
        self.pdb_delete_btn = QPushButton("Delete")
        self.pdb_delete_btn.setStyleSheet("background: #EF5350; color: white; padding: 4px 10px; border: none; border-radius: 4px;")
        pdb_top.addWidget(self.pdb_delete_btn)
        self.pdb_load_seq_btn = QPushButton("Load as Sequence")
        self.pdb_load_seq_btn.setStyleSheet("background: #FFB74D; color: #0B1020; padding: 4px 10px; border: none; border-radius: 4px;")
        pdb_top.addWidget(self.pdb_load_seq_btn)
        pdb_top.addStretch()
        pdb_browse_layout.addLayout(pdb_top)
        self.pdb_list = QTreeWidget()
        self.pdb_list.setHeaderLabels(["PDB File", "Size", "Modified"])
        self.pdb_list.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.pdb_list.setAlternatingRowColors(True)
        self.pdb_list.setMinimumHeight(120)
        pdb_browse_layout.addWidget(self.pdb_list)
        self.pdb_project_label = QLabel("No project open")
        pdb_browse_layout.addWidget(self.pdb_project_label)
        layout.addWidget(pdb_browse_group)

        # External Tools Section
        tools_group = QGroupBox("External Tools")
        tools_layout = QVBoxLayout(tools_group)

        # Tool status
        self.tool_status_label = QLabel()
        tools_layout.addWidget(self.tool_status_label)
        self._refresh_tool_status()

        tool_buttons = QHBoxLayout()
        self.launch_pymol_btn = QPushButton("Launch PyMOL")
        self.launch_pymol_btn.setStyleSheet("background: #00D1FF; color: #0B1020; padding: 6px; border: none; border-radius: 4px;")
        tool_buttons.addWidget(self.launch_pymol_btn)
        self.launch_chimerax_btn = QPushButton("Launch ChimeraX")
        self.launch_chimerax_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        tool_buttons.addWidget(self.launch_chimerax_btn)
        self.launch_snapgene_btn = QPushButton("Launch SnapGene")
        self.launch_snapgene_btn.setStyleSheet("background: #22C55E; color: #0B1020; padding: 6px; border: none; border-radius: 4px;")
        tool_buttons.addWidget(self.launch_snapgene_btn)
        self.export_for_tools_btn = QPushButton("Export FASTA for Tools")
        self.export_for_tools_btn.setStyleSheet("background: #8B5CF6; color: white; padding: 6px; border: none; border-radius: 4px;")
        tool_buttons.addWidget(self.export_for_tools_btn)
        tool_buttons.addStretch()
        tools_layout.addLayout(tool_buttons)
        layout.addWidget(tools_group)

        # Plugins Section
        plugin_group = QGroupBox("Plugins")
        plugin_layout = QVBoxLayout(plugin_group)
        plugin_top = QHBoxLayout()
        self.discover_plugins_btn = QPushButton("Discover Plugins")
        self.discover_plugins_btn.setStyleSheet("background: #00D1FF; color: #0B1020; padding: 6px; border: none; border-radius: 4px;")
        plugin_top.addWidget(self.discover_plugins_btn)
        self.plugin_status_label = QLabel("No plugins discovered")
        plugin_top.addWidget(self.plugin_status_label)
        plugin_top.addStretch()
        plugin_layout.addLayout(plugin_top)
        self.plugin_output = QTextEdit()
        self.plugin_output.setReadOnly(True)
        self.plugin_output.setMaximumHeight(100)
        self.plugin_output.setPlaceholderText("Plugin output...")
        plugin_layout.addWidget(self.plugin_output)
        layout.addWidget(plugin_group)

        # Initialise plugin manager
        self.plugin_manager = PluginManager()
        layout.addStretch()

    def _build_docking_tab(self):
        layout = QVBoxLayout(self.docking_tab)
        layout.setSpacing(6)

        # Status banner
        dock_status = DockingEngine.check_available()
        status_text = ' | '.join(f"{'✓' if v else '✗'} {k}" for k, v in dock_status.items())
        status_bar = QLabel(f"Docking tools: {status_text}")
        status_bar.setStyleSheet("color: #8B9DC3; padding: 4px;")
        layout.addWidget(status_bar)

        # Receptor section
        rec_group = QGroupBox("Receptor")
        rec_layout = QVBoxLayout(rec_group)
        rec_row = QHBoxLayout()
        rec_row.addWidget(QLabel("PDB file:"))
        self.dock_receptor_path = QLineEdit()
        self.dock_receptor_path.setPlaceholderText("Path to receptor PDB file...")
        rec_row.addWidget(self.dock_receptor_path)
        self.dock_rec_browse_btn = self._btn("Browse...", '#6B7280', 'white')
        rec_row.addWidget(self.dock_rec_browse_btn)
        self.dock_rec_auto_btn = self._btn("Auto-detect box", '#00D1FF')
        rec_row.addWidget(self.dock_rec_auto_btn)
        rec_layout.addLayout(rec_row)
        fetch_row = QHBoxLayout()
        fetch_row.addWidget(QLabel("Fetch PDB ID:"))
        self.dock_pdb_id = QLineEdit()
        self.dock_pdb_id.setPlaceholderText("e.g. 4HHB, 1CRN...")
        self.dock_pdb_id.setMaximumWidth(120)
        fetch_row.addWidget(self.dock_pdb_id)
        self.dock_fetch_pdb_btn = self._btn("Download", '#00D1FF')
        fetch_row.addWidget(self.dock_fetch_pdb_btn)
        fetch_row.addStretch()
        rec_layout.addLayout(fetch_row)
        layout.addWidget(rec_group)

        # Ligand section
        lig_group = QGroupBox("Ligand")
        lig_layout = QVBoxLayout(lig_group)
        lig_row = QHBoxLayout()
        lig_row.addWidget(QLabel("Source:"))
        self.dock_lig_source = QComboBox()
        self.dock_lig_source.addItems(["SMILES", "File (PDB/SDF/MOL2)"])
        lig_row.addWidget(self.dock_lig_source)
        lig_row.addStretch()
        lig_layout.addLayout(lig_row)
        self.dock_lig_smiles = QLineEdit()
        self.dock_lig_smiles.setPlaceholderText("e.g. CCO (ethanol), CC1=CC=C(C=C1)C(=O)O (ibuprofen)")
        self.dock_lig_smiles.setToolTip("Enter a SMILES string for the ligand")
        lig_layout.addWidget(self.dock_lig_smiles)
        lig_file_row = QHBoxLayout()
        self.dock_lig_path = QLineEdit()
        self.dock_lig_path.setPlaceholderText("Path to ligand file...")
        self.dock_lig_path.setVisible(False)
        lig_file_row.addWidget(self.dock_lig_path)
        self.dock_lig_browse_btn = self._btn("Browse...", '#6B7280', 'white')
        self.dock_lig_browse_btn.setVisible(False)
        lig_file_row.addWidget(self.dock_lig_browse_btn)
        lig_file_row.addStretch()
        lig_layout.addLayout(lig_file_row)
        layout.addWidget(lig_group)

        # Binding box section
        box_group = QGroupBox("Binding Box")
        box_layout = QGridLayout(box_group)
        box_layout.addWidget(QLabel("Center X:"), 0, 0)
        self.dock_cx = QDoubleSpinBox()
        self.dock_cx.setRange(-999, 999)
        self.dock_cx.setDecimals(2)
        box_layout.addWidget(self.dock_cx, 0, 1)
        box_layout.addWidget(QLabel("Y:"), 0, 2)
        self.dock_cy = QDoubleSpinBox()
        self.dock_cy.setRange(-999, 999)
        self.dock_cy.setDecimals(2)
        box_layout.addWidget(self.dock_cy, 0, 3)
        box_layout.addWidget(QLabel("Z:"), 0, 4)
        self.dock_cz = QDoubleSpinBox()
        self.dock_cz.setRange(-999, 999)
        self.dock_cz.setDecimals(2)
        box_layout.addWidget(self.dock_cz, 0, 5)

        box_layout.addWidget(QLabel("Size X:"), 1, 0)
        self.dock_sx = QDoubleSpinBox()
        self.dock_sx.setRange(1, 500)
        self.dock_sx.setValue(30)
        box_layout.addWidget(self.dock_sx, 1, 1)
        box_layout.addWidget(QLabel("Y:"), 1, 2)
        self.dock_sy = QDoubleSpinBox()
        self.dock_sy.setRange(1, 500)
        self.dock_sy.setValue(30)
        box_layout.addWidget(self.dock_sy, 1, 3)
        box_layout.addWidget(QLabel("Z:"), 1, 4)
        self.dock_sz = QDoubleSpinBox()
        self.dock_sz.setRange(1, 500)
        self.dock_sz.setValue(30)
        box_layout.addWidget(self.dock_sz, 1, 5)

        box_layout.addWidget(QLabel("Exhaustiveness:"), 2, 0)
        self.dock_exhaust = QSpinBox()
        self.dock_exhaust.setRange(1, 64)
        self.dock_exhaust.setValue(8)
        box_layout.addWidget(self.dock_exhaust, 2, 1)
        box_layout.addWidget(QLabel("Modes:"), 2, 2)
        self.dock_modes = QSpinBox()
        self.dock_modes.setRange(1, 50)
        self.dock_modes.setValue(9)
        box_layout.addWidget(self.dock_modes, 2, 3)
        layout.addWidget(box_group)

        # Action buttons
        action_row = QHBoxLayout()
        self.dock_run_btn = self._btn("Run Docking", '#22C55E', bold=True)
        action_row.addWidget(self.dock_run_btn)
        self.dock_cancel_btn = self._btn("Cancel", '#EF4444', 'white')
        self.dock_cancel_btn.setEnabled(False)
        action_row.addWidget(self.dock_cancel_btn)
        action_row.addStretch()
        layout.addLayout(action_row)

        # Output log
        self.dock_output = QTextEdit()
        self.dock_output.setReadOnly(True)
        self.dock_output.setFont(QFont("Courier New", 10))
        self.dock_output.setMaximumHeight(150)
        self.dock_output.setPlaceholderText("Docking output will appear here...")
        layout.addWidget(self.dock_output)

        # Results table
        res_group = QGroupBox("Docking Results")
        res_layout = QVBoxLayout(res_group)
        self.dock_results_table = QTableWidget()
        self.dock_results_table.setColumnCount(5)
        self.dock_results_table.setHorizontalHeaderLabels(
            ["Mode", "Affinity (kcal/mol)", "RMSD LB", "RMSD UB", "Rank"]
        )
        self.dock_results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.dock_results_table.setAlternatingRowColors(True)
        res_layout.addWidget(self.dock_results_table)
        res_buttons = QHBoxLayout()
        self.dock_vis_btn = self._btn("Visualize (chart)", '#00D1FF')
        res_buttons.addWidget(self.dock_vis_btn)
        self.dock_open_pymol_btn = self._btn("Open both in PyMOL", '#8B5CF6', 'white', bold=True)
        self.dock_open_pymol_btn.setToolTip("Opens receptor + best docking pose together in one PyMOL session")
        res_buttons.addWidget(self.dock_open_pymol_btn)
        self.dock_poses_pymol_btn = self._btn("All poses in PyMOL", '#A78BFA', 'white')
        self.dock_poses_pymol_btn.setToolTip("Opens all docking poses as overlaid structures in PyMOL")
        res_buttons.addWidget(self.dock_poses_pymol_btn)
        self.dock_export_csv_btn = self._btn("Export CSV", '#6B7280', 'white')
        res_buttons.addWidget(self.dock_export_csv_btn)
        res_buttons.addStretch()
        res_layout.addLayout(res_buttons)
        layout.addWidget(res_group)

    def _refresh_tool_status(self):
        tools = ExternalToolLauncher.get_available_tools()
        statuses = []
        for name, available in tools.items():
            icon = "✓" if available else "✗"
            statuses.append(f"{icon} {name}")
        self.tool_status_label.setText(" | ".join(statuses))
        self.tool_status_label.setStyleSheet("color: #8B9DC3; padding: 4px;")

    # -------- Menu ----------
    def _build_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        open_action = QAction("&Open Sequence...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_file)
        file_menu.addAction(open_action)

        import_menu = file_menu.addMenu("&Import")
        import_fasta = QAction("FASTA", self)
        import_fasta.triggered.connect(lambda: self._open_file())
        import_menu.addAction(import_fasta)
        import_fastq = QAction("FASTQ", self)
        import_fastq.triggered.connect(lambda: self._open_file())
        import_menu.addAction(import_fastq)
        import_genbank = QAction("GenBank", self)
        import_genbank.triggered.connect(lambda: self._open_file())
        import_menu.addAction(import_genbank)
        import_csv = QAction("CSV", self)
        import_csv.triggered.connect(lambda: self._open_file())
        import_menu.addAction(import_csv)

        file_menu.addSeparator()

        # Project actions
        new_project_act = QAction("&New Project...", self)
        new_project_act.setShortcut("Ctrl+Shift+N")
        new_project_act.triggered.connect(self._new_project)
        file_menu.addAction(new_project_act)

        open_project_act = QAction("&Open Project...", self)
        open_project_act.setShortcut("Ctrl+Shift+O")
        open_project_act.triggered.connect(self._open_project)
        file_menu.addAction(open_project_act)

        self.save_project_act = QAction("&Save Project", self)
        self.save_project_act.setShortcut("Ctrl+S")
        self.save_project_act.triggered.connect(self._save_project)
        file_menu.addAction(self.save_project_act)

        self.save_project_as_act = QAction("Save Project &As...", self)
        self.save_project_as_act.setShortcut("Ctrl+Shift+S")
        self.save_project_as_act.triggered.connect(self._save_project_as)
        file_menu.addAction(self.save_project_as_act)

        self.close_project_act = QAction("&Close Project", self)
        self.close_project_act.triggered.connect(self._close_project)
        file_menu.addAction(self.close_project_act)

        file_menu.addSeparator()

        paste_action = QAction("&Paste Sequence", self)
        paste_action.setShortcut("Ctrl+V")
        paste_action.triggered.connect(self._paste_sequence)
        file_menu.addAction(paste_action)

        file_menu.addSeparator()

        self.export_csv_action = QAction("Export Mutations as CSV", self)
        self.export_csv_action.triggered.connect(self._export_csv)
        file_menu.addAction(self.export_csv_action)

        self.export_json_action = QAction("Export as JSON", self)
        self.export_json_action.triggered.connect(self._export_json)
        file_menu.addAction(self.export_json_action)

        self.export_fasta_action = QAction("Export Sequence as FASTA", self)
        self.export_fasta_action.triggered.connect(self._export_fasta)
        file_menu.addAction(self.export_fasta_action)

        self.export_report_action = QAction("Export Report (Text)", self)
        self.export_report_action.triggered.connect(self._export_report)
        file_menu.addAction(self.export_report_action)

        file_menu.addSeparator()
        export_svg_action = QAction("Export Graph as SVG", self)
        export_svg_action.triggered.connect(self._export_svg)
        file_menu.addAction(export_svg_action)
        export_png_action = QAction("Export Graph as PNG", self)
        export_png_action.triggered.connect(self._export_png)
        file_menu.addAction(export_png_action)

        file_menu.addSeparator()
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        edit_menu = menubar.addMenu("&Edit")
        clear_action = QAction("Clear All Sequences", self)
        clear_action.triggered.connect(self._clear_all)
        edit_menu.addAction(clear_action)

        tools_menu = menubar.addMenu("&Tools")
        analyze_action = QAction("Analyze Sequence", self)
        analyze_action.triggered.connect(self._run_analysis)
        tools_menu.addAction(analyze_action)
        detect_action = QAction("Detect Mutations", self)
        detect_action.triggered.connect(self._detect_mutations)
        tools_menu.addAction(detect_action)
        tools_menu.addSeparator()
        find_motif_action = QAction("Find Motif / Pattern...", self)
        find_motif_action.setShortcut("Ctrl+F")
        find_motif_action.triggered.connect(self._find_motif_dialog)
        tools_menu.addAction(find_motif_action)
        revcomp_action = QAction("Reverse Complement", self)
        revcomp_action.setShortcut("Ctrl+R")
        revcomp_action.triggered.connect(self._show_revcomp)
        tools_menu.addAction(revcomp_action)
        orf_action = QAction("Show Open Reading Frames", self)
        orf_action.setShortcut("Ctrl+G")
        orf_action.triggered.connect(self._show_orfs)
        tools_menu.addAction(orf_action)
        tools_menu.addSeparator()
        annot_action = QAction("Lookup Annotation (ClinVar/dbSNP)...", self)
        annot_action.triggered.connect(self._annotation_lookup)
        tools_menu.addAction(annot_action)
        tools_menu.addSeparator()
        ncbi_act = QAction("Fetch from NCBI...", self)
        ncbi_act.triggered.connect(self._fetch_from_ncbi)
        tools_menu.addAction(ncbi_act)
        pdb_act = QAction("Import PDB File...", self)
        pdb_act.triggered.connect(self._import_pdb)
        tools_menu.addAction(pdb_act)
        pdb_fetch_act = QAction("Fetch PDB from RCSB...", self)
        pdb_fetch_act.triggered.connect(self._fetch_pdb_rcsb)
        tools_menu.addAction(pdb_fetch_act)
        tools_menu.addSeparator()
        ext_menu = tools_menu.addMenu("External Tools")
        pm_act = QAction("Launch PyMOL", self)
        pm_act.triggered.connect(self._launch_pymol)
        ext_menu.addAction(pm_act)
        cx_act = QAction("Launch ChimeraX", self)
        cx_act.triggered.connect(self._launch_chimerax)
        ext_menu.addAction(cx_act)
        sg_act = QAction("Launch SnapGene", self)
        sg_act.triggered.connect(self._launch_snapgene)
        ext_menu.addAction(sg_act)
        exp_act = QAction("Export FASTA for External Tools", self)
        exp_act.triggered.connect(self._export_for_tools)
        ext_menu.addAction(exp_act)
        tools_menu.addSeparator()
        disc_act = QAction("Discover Plugins", self)
        disc_act.triggered.connect(self._discover_plugins)
        tools_menu.addAction(disc_act)

        view_menu = menubar.addMenu("&View")
        self.toggle_project = QAction("Toggle Project Explorer", self, checkable=True)
        self.toggle_project.setChecked(True)
        self.toggle_project.triggered.connect(lambda: self.project_tree.setVisible(self.toggle_project.isChecked()))
        view_menu.addAction(self.toggle_project)
        view_menu.addSeparator()

        theme_menu = view_menu.addMenu("Theme")
        self.theme_actions = {}
        theme_group = QActionGroup(self)
        theme_group.setExclusive(True)
        for tname in ThemeManager.get_theme_names():
            act = QAction(tname, self, checkable=True)
            act.setChecked(tname == 'Deep Scientific Navy')
            act.triggered.connect(lambda checked, n=tname: self._switch_theme(n))
            theme_group.addAction(act)
            theme_menu.addAction(act)
            self.theme_actions[tname] = act

        theme_menu.addSeparator()
        contrast_menu = theme_menu.addMenu("Contrast Variation")
        self.contrast_actions = {}
        contrast_group = QActionGroup(self)
        contrast_group.setExclusive(True)
        self.current_contrast = 'Default (Theme)'
        for vname in ThemeManager.CONTRAST_VARIANTS:
            act = QAction(vname, self, checkable=True)
            act.setChecked(vname == 'Default (Theme)')
            act.triggered.connect(lambda checked, v=vname: self._apply_contrast_variant(v))
            contrast_group.addAction(act)
            contrast_menu.addAction(act)
            self.contrast_actions[vname] = act

        help_menu = menubar.addMenu("&Help")
        about_action = QAction("About ZORA", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    def _build_toolbar(self):
        toolbar = QToolBar("Main Toolbar")
        toolbar.setIconSize(QSize(20, 20))
        self.addToolBar(toolbar)

        open_btn = QAction("Open", self)
        open_btn.setToolTip("Open sequence file")
        open_btn.triggered.connect(self._open_file)
        toolbar.addAction(open_btn)

        paste_btn = QAction("Paste", self)
        paste_btn.setToolTip("Paste sequence from clipboard")
        paste_btn.triggered.connect(self._paste_sequence)
        toolbar.addAction(paste_btn)

        toolbar.addSeparator()

        analyze_btn = QAction("Analyze", self)
        analyze_btn.setToolTip("Analyze current sequence")
        analyze_btn.triggered.connect(self._run_analysis)
        toolbar.addAction(analyze_btn)

        detect_btn = QAction("Detect SNPs", self)
        detect_btn.setToolTip("Detect mutations between two sequences")
        detect_btn.triggered.connect(self._detect_mutations)
        toolbar.addAction(detect_btn)

        find_btn = QAction("Find", self)
        find_btn.setToolTip("Find motif/pattern in sequence")
        find_btn.triggered.connect(self._find_motif_dialog)
        toolbar.addAction(find_btn)

        revcomp_btn = QAction("RevComp", self)
        revcomp_btn.setToolTip("Show reverse complement")
        revcomp_btn.triggered.connect(self._show_revcomp)
        toolbar.addAction(revcomp_btn)

        toolbar.addSeparator()

        gc_btn = QAction("GC Plot", self)
        gc_btn.setToolTip("GC Content Plot")
        gc_btn.triggered.connect(self._plot_gc_content)
        toolbar.addAction(gc_btn)

        export_btn = QAction("Export CSV", self)
        export_btn.setToolTip("Export mutations as CSV")
        export_btn.triggered.connect(self._export_csv)
        toolbar.addAction(export_btn)

    def _connect_signals(self):
        self.drop_area.sequences_dropped.connect(self._add_sequences)
        self.analyze_btn.clicked.connect(self._run_analysis)
        self.detect_btn.clicked.connect(self._detect_mutations)
        self.simulate_btn.clicked.connect(self._simulate_mutation)
        self.insert_btn.clicked.connect(self._simulate_insertion)
        self.delete_btn.clicked.connect(self._simulate_deletion)
        self.batch_btn.clicked.connect(self._batch_mutations)
        self.gc_plot_btn.clicked.connect(self._plot_gc_content)
        self.mut_dist_btn.clicked.connect(self._plot_mutation_distribution)
        self.nuc_freq_btn.clicked.connect(self._plot_nucleotide_frequency)
        self.codon_usage_btn.clicked.connect(self._plot_codon_usage)
        self.mut_pie_btn.clicked.connect(self._plot_mutation_pie)
        self.heatmap_btn.clicked.connect(self._plot_snp_heatmap)
        self.conservation_btn.clicked.connect(self._plot_conservation)
        self.toggle_grid_btn.clicked.connect(self._toggle_graph_grid)
        self.export_svg_btn.clicked.connect(self._export_svg)
        self.export_png_btn.clicked.connect(self._export_png)
        self.aln_nw_btn.clicked.connect(lambda: self._run_alignment('nw'))
        self.aln_sw_btn.clicked.connect(lambda: self._run_alignment('sw'))
        self.viewer_goto_btn.clicked.connect(self._viewer_goto)
        self.viewer_seq_display.cursorPositionChanged.connect(self._viewer_cursor_changed)
        self.crispr_design_btn.clicked.connect(self._design_crispr)
        self.ai_analyze_btn.clicked.connect(self._run_ai_analysis)
        self.mutation_table.itemSelectionChanged.connect(self._mutation_selected)
        self.view_mut_in_seq_btn.clicked.connect(self._view_mutations_in_sequence)
        # Integrations signals
        self.db_fetch_btn.clicked.connect(self._fetch_from_ncbi)
        self.db_search_btn.clicked.connect(self._search_ncbi_dialog)
        self.pdb_import_btn.clicked.connect(self._import_pdb)
        self.pdb_fetch_btn.clicked.connect(self._fetch_pdb_rcsb)
        self.launch_pymol_btn.clicked.connect(self._launch_pymol)
        self.launch_chimerax_btn.clicked.connect(self._launch_chimerax)
        self.send_pymol_btn.clicked.connect(self._send_seq_to_pymol)
        self.send_chimerax_btn.clicked.connect(self._send_seq_to_chimerax)
        self.send_snapgene_btn.clicked.connect(self._send_seq_to_snapgene)
        self.launch_snapgene_btn.clicked.connect(self._launch_snapgene)
        self.export_for_tools_btn.clicked.connect(self._export_for_tools)
        self.discover_plugins_btn.clicked.connect(self._discover_plugins)
        # PDB browser signals
        self.pdb_refresh_btn.clicked.connect(self._refresh_pdb_list)
        self.pdb_open_pymol_btn.clicked.connect(self._open_selected_pdb_in_pymol)
        self.pdb_open_chimerax_btn.clicked.connect(self._open_selected_pdb_in_chimerax)
        self.pdb_delete_btn.clicked.connect(self._delete_selected_pdb)
        self.pdb_load_seq_btn.clicked.connect(self._load_selected_pdb_as_sequence)

        # Docking signals
        self.dock_rec_browse_btn.clicked.connect(self._dock_browse_receptor)
        self.dock_lig_browse_btn.clicked.connect(self._dock_browse_ligand)
        self.dock_lig_source.currentTextChanged.connect(self._dock_lig_source_changed)
        self.dock_rec_auto_btn.clicked.connect(self._dock_auto_detect)
        self.dock_fetch_pdb_btn.clicked.connect(self._dock_fetch_pdb)
        self.dock_pdb_id.returnPressed.connect(self._dock_fetch_pdb)
        self.dock_run_btn.clicked.connect(self._dock_run)
        self.dock_cancel_btn.clicked.connect(self._dock_cancel)
        self.dock_vis_btn.clicked.connect(self._dock_visualize)
        self.dock_open_pymol_btn.clicked.connect(self._dock_open_pymol)
        self.dock_poses_pymol_btn.clicked.connect(self._dock_open_all_poses)
        self.dock_export_csv_btn.clicked.connect(self._dock_export_csv)

    # -------- Project Methods ----------
    def _check_recovery(self):
        rec_dir = os.path.join(os.path.expanduser("~"), ".zora", "recovery")
        if os.path.isdir(rec_dir):
            meta = os.path.join(rec_dir, "project.json")
            if os.path.exists(meta):
                reply = QMessageBox.question(self, "Recovery Found",
                    "A recovery project was found from a previous session.\n"
                    "Would you like to restore it?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.Yes:
                    self.project_manager.open(rec_dir)
                    self._restore_from_project()
                    self.statusBar().showMessage("Recovery project restored", 5000)
                else:
                    import shutil
                    shutil.rmtree(rec_dir, ignore_errors=True)

    def _new_project(self):
        path = QFileDialog.getExistingDirectory(self, "Choose Project Directory")
        if not path:
            return
        name = os.path.basename(path)
        ok = self.project_manager.create(path, name)
        if ok:
            self._update_project_title()
            self._clear_all()
            self._refresh_pdb_list()
            self._update_project_tree()
            self.statusBar().showMessage(f"Created project: {name}", 5000)
        else:
            QMessageBox.warning(self, "Error", "Failed to create project.")

    def _open_project(self):
        path = QFileDialog.getExistingDirectory(self, "Open Project Directory",
                                                self.project_manager.project_dir or "")
        if not path:
            return
        ok = self.project_manager.open(path)
        if ok:
            self._restore_from_project()
            self._update_project_title()
            self._refresh_pdb_list()
            self.statusBar().showMessage(f"Opened project: {self.project_manager.metadata.get('name', path)}", 5000)
        else:
            QMessageBox.warning(self, "Error", "Failed to open project.")

    def _save_project(self):
        if not self.project_manager.is_open():
            self._save_project_as()
            return
        self.project_manager.save(
            sequences=self.sequences,
            mutation_results=self.mutation_results,
            alignment_results=self.alignment_results,
            analysis_text=self.analysis_output.toPlainText(),
        )
        self._save_recovery()
        self.statusBar().showMessage("Project saved", 3000)

    def _save_project_as(self):
        path = QFileDialog.getExistingDirectory(self, "Save Project As",
                                                self.project_manager.project_dir or "")
        if not path:
            return
        name = os.path.basename(path)
        if not self.project_manager.is_open() or self.project_manager.project_dir != path:
            self.project_manager.create(path, name)
        self._save_project()

    def _close_project(self):
        self.project_manager.save(
            sequences=self.sequences,
            mutation_results=self.mutation_results,
        )
        self.project_manager.close()
        self._clear_all()
        self._update_project_tree()
        self._refresh_pdb_list()
        self._update_project_title()
        self._remove_recovery()
        self.statusBar().showMessage("Project closed", 3000)

    def _restore_from_project(self):
        self._clear_all()
        records = self.project_manager.load_sequences()
        if records:
            self._add_sequences(records)
        muts = self.project_manager.load_mutations()
        if muts:
            for m in muts:
                self.mutation_results.append(MutationResult(**m))
        analysis_path = os.path.join(self.project_manager.project_dir, 'analysis', 'analysis.txt')
        if os.path.exists(analysis_path):
            with open(analysis_path) as f:
                self.analysis_output.setPlainText(f.read())
        self._update_project_title()

    def _save_recovery(self):
        if not self.project_manager.is_open():
            return
        rec_dir = os.path.join(os.path.expanduser("~"), ".zora", "recovery")
        src = self.project_manager.project_dir
        if not src or not os.path.isdir(src):
            return
        # Avoid copy onto itself when already in recovery mode
        if os.path.realpath(src) == os.path.realpath(rec_dir):
            return
        import shutil
        if os.path.exists(rec_dir):
            shutil.rmtree(rec_dir, ignore_errors=True)
        shutil.copytree(src, rec_dir, dirs_exist_ok=True)

    def _remove_recovery(self):
        rec_dir = os.path.join(os.path.expanduser("~"), ".zora", "recovery")
        import shutil
        if os.path.exists(rec_dir):
            shutil.rmtree(rec_dir, ignore_errors=True)

    def _update_project_title(self):
        name = self.project_manager.metadata.get('name', 'Untitled')
        if self.project_manager.is_open():
            self.setWindowTitle(f"ZORA - {name}")
        else:
            self.setWindowTitle("ZORA - Sequence & Mutation Analysis Workstation")

    # -------- PDB Browser Methods ----------
    def _refresh_pdb_list(self):
        self.pdb_list.clear()
        if not self.project_manager.is_open():
            self.pdb_project_label.setText("No project open — PDBs not saved locally")
            self.pdb_open_pymol_btn.setEnabled(False)
            self.pdb_open_chimerax_btn.setEnabled(False)
            self.pdb_delete_btn.setEnabled(False)
            self.pdb_load_seq_btn.setEnabled(False)
            return
        self.pdb_open_pymol_btn.setEnabled(True)
        self.pdb_open_chimerax_btn.setEnabled(True)
        self.pdb_delete_btn.setEnabled(True)
        self.pdb_load_seq_btn.setEnabled(True)
        proj_name = self.project_manager.metadata.get('name', 'Untitled')
        self.pdb_project_label.setText(f"Project: {proj_name}  —  {self.project_manager.project_dir}")
        files = self.project_manager.get_pdb_files()
        if not files:
            item = QTreeWidgetItem(["<no PDB files in project>"])
            item.setFlags(Qt.NoItemFlags)
            self.pdb_list.addTopLevelItem(item)
            return
        for f in files:
            size_str = f"{f['size']:,} B"
            modified = f.get('modified', '')[:19]
            item = QTreeWidgetItem([f['name'], size_str, modified])
            item.setData(0, Qt.UserRole, f['path'])
            self.pdb_list.addTopLevelItem(item)
        self.pdb_list.resizeColumnToContents(0)

    def _get_selected_pdb_paths(self) -> list:
        items = self.pdb_list.selectedItems()
        paths = []
        for item in items:
            path = item.data(0, Qt.UserRole)
            if path and os.path.exists(path):
                paths.append(path)
        return paths

    def _open_selected_pdb_in_pymol(self):
        paths = self._get_selected_pdb_paths()
        if not paths:
            QMessageBox.information(self, "No Selection", "Select PDB files from the list first.")
            return
        for p in paths:
            ExternalToolLauncher.launch_pymol(filename=p)
        self.statusBar().showMessage(f"Opened {len(paths)} PDB(s) in PyMOL", 3000)

    def _open_selected_pdb_in_chimerax(self):
        paths = self._get_selected_pdb_paths()
        if not paths:
            QMessageBox.information(self, "No Selection", "Select PDB files from the list first.")
            return
        for p in paths:
            ExternalToolLauncher.launch_chimerax(filename=p)
        self.statusBar().showMessage(f"Opened {len(paths)} PDB(s) in ChimeraX", 3000)

    def _delete_selected_pdb(self):
        items = self.pdb_list.selectedItems()
        if not items:
            QMessageBox.information(self, "No Selection", "Select PDB files to delete.")
            return
        names = [item.text(0) for item in items]
        reply = QMessageBox.warning(self, "Confirm Delete",
            f"Delete {len(names)} PDB file(s)?\n\n" + "\n".join(names),
            QMessageBox.Yes | QMessageBox.No)
        if reply != QMessageBox.Yes:
            return
        for item in items:
            self.project_manager.delete_pdb(item.text(0))
        self._refresh_pdb_list()
        self.statusBar().showMessage(f"Deleted {len(names)} PDB file(s)", 3000)

    def _load_selected_pdb_as_sequence(self):
        paths = self._get_selected_pdb_paths()
        if not paths:
            QMessageBox.information(self, "No Selection", "Select a PDB file first.")
            return
        records = []
        for p in paths:
            records.extend(PDBParser.parse(p))
        if records:
            self._add_sequences(records)
            self.statusBar().showMessage(f"Loaded {len(records)} chain(s) from PDB", 3000)
        else:
            QMessageBox.warning(self, "No Data", "Could not extract sequences from selected PDB(s).")

    # -------- Docking Methods ----------
    def _dock_browse_receptor(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Receptor", "", "PDB/PDBQT (*.pdb *.pdbqt);;All files (*)")
        if path:
            self.dock_receptor_path.setText(path)

    def _dock_browse_ligand(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Ligand", "", "Molecule files (*.pdb *.pdbqt *.sdf *.mol2);;All files (*)")
        if path:
            self.dock_lig_path.setText(path)

    def _dock_lig_source_changed(self, source):
        is_smiles = source == "SMILES"
        self.dock_lig_smiles.setVisible(is_smiles)
        self.dock_lig_path.setVisible(not is_smiles)
        self.dock_lig_browse_btn.setVisible(not is_smiles)

    def _dock_auto_detect(self):
        pdb_path = self.dock_receptor_path.text()
        if not pdb_path or not os.path.exists(pdb_path):
            QMessageBox.warning(self, "No Receptor", "Set a receptor PDB file first.")
            return
        box = DockingEngine.get_pdb_center(pdb_path)
        if box is None:
            QMessageBox.warning(self, "Parse Error", "Could not parse PDB coordinates.")
            return
        self.dock_cx.setValue(box['center_x'])
        self.dock_cy.setValue(box['center_y'])
        self.dock_cz.setValue(box['center_z'])
        self.dock_sx.setValue(box['size_x'])
        self.dock_sy.setValue(box['size_y'])
        self.dock_sz.setValue(box['size_z'])
        self.dock_output.append("✓ Binding box auto-detected from receptor coordinates.")

    def _dock_run(self):
        if not DockingEngine.check_available().get('vina', False):
            QMessageBox.critical(self, "Vina Not Found",
                                "Vina executable not found. Install it at /usr/local/bin/vina")
            return

        receptor_path = self.dock_receptor_path.text()
        if not receptor_path or not os.path.exists(receptor_path):
            QMessageBox.warning(self, "Receptor Required", "Select a valid receptor PDB file.")
            return

        # Prepare output directory
        output_dir = os.path.join(
            self.project_manager.project_dir if self.project_manager.is_open() else os.path.expanduser("~"),
            ".zora_docking")
        os.makedirs(output_dir, exist_ok=True)

        self.dock_output.clear()
        self.dock_output.append("Preparing receptor...")
        QApplication.processEvents()

        if receptor_path.lower().endswith('.pdbqt'):
            receptor_pdbqt = DockingEngine._sanitize_receptor_pdbqt(receptor_path, output_dir)
            self.dock_output.append(f"✓ Receptor PDBQT: {receptor_pdbqt}")
        else:
            prepared = DockingEngine.prepare_receptor(receptor_path, output_dir)
            if not prepared:
                self.dock_output.append("✗ Failed to prepare receptor PDBQT.")
                return
            receptor_pdbqt = DockingEngine._sanitize_receptor_pdbqt(prepared, output_dir)
            self.dock_output.append(f"✓ Receptor: {receptor_pdbqt}")

        # Prepare ligand
        self.dock_output.append("Preparing ligand...")
        QApplication.processEvents()

        ligand_pdbqt = None
        if self.dock_lig_source.currentText() == "SMILES":
            smiles = self.dock_lig_smiles.text().strip()
            if smiles:
                ligand_pdbqt = DockingEngine.ligand_from_smiles(smiles, output_dir)
                if ligand_pdbqt:
                    self.dock_output.append(f"✓ Ligand from SMILES: {ligand_pdbqt}")
                else:
                    self.dock_output.append("✗ Failed to generate ligand from SMILES.")
                    return
        else:
            lig_path = self.dock_lig_path.text()
            if lig_path and os.path.exists(lig_path):
                if lig_path.lower().endswith('.pdbqt'):
                    ligand_pdbqt = lig_path
                    self.dock_output.append(f"✓ Ligand PDBQT: {ligand_pdbqt}")
                else:
                    ligand_pdbqt = DockingEngine.prepare_ligand(lig_path, output_dir)
                    if not ligand_pdbqt:
                        self.dock_output.append("✗ Failed to prepare ligand PDBQT.")
                        return
                    self.dock_output.append(f"✓ Ligand: {ligand_pdbqt}")

        if not ligand_pdbqt:
            QMessageBox.warning(self, "Ligand Required", "Provide a ligand (SMILES or file).")
            return

        # Disable run button
        self.dock_run_btn.setEnabled(False)
        self.dock_cancel_btn.setEnabled(True)
        self.dock_output.append("Running docking (this may take a while)...")

        box_params = {
            'center_x': self.dock_cx.value(),
            'center_y': self.dock_cy.value(),
            'center_z': self.dock_cz.value(),
            'size_x': self.dock_sx.value(),
            'size_y': self.dock_sy.value(),
            'size_z': self.dock_sz.value(),
            'exhaustiveness': self.dock_exhaust.value(),
            'num_modes': self.dock_modes.value(),
        }

        self._dock_worker = DockingWorker(receptor_pdbqt, ligand_pdbqt, output_dir, box_params)
        self._dock_worker.signals.finished.connect(self._dock_finished)
        self._dock_worker.signals.error.connect(self._dock_error)
        QThreadPool.globalInstance().start(self._dock_worker)

    def _dock_finished(self, result, modes):
        self._last_docking_result = result
        self.dock_run_btn.setEnabled(True)
        self.dock_cancel_btn.setEnabled(False)
        self.docking_results = modes

        if result['returncode'] != 0:
            self.dock_output.append(f"✗ Vina returned error code {result['returncode']}")
            if result['stderr']:
                self.dock_output.append(f"STDERR: {result['stderr']}")
            return

        self.dock_output.append("✓ Docking complete!")
        if result['out_path']:
            self.dock_output.append(f"Results: {result['out_path']}")
            # Save to project and split poses
            if self.project_manager.is_open():
                with open(result['out_path'], 'rb') as f:
                    content = f.read()
                saved = self.project_manager.save_docking_pdbqt('docking_result', content)
                if saved:
                    self.dock_output.append(f"✓ Saved to project: {saved}")
                    poses = self.project_manager.split_docking_poses(saved)
                    if poses:
                        self.dock_output.append(f"✓ Split into {len(poses)} individual poses")
                    self._update_project_tree()

        # Populate table
        self.dock_results_table.setRowCount(len(modes))
        for i, m in enumerate(modes):
            self.dock_results_table.setItem(i, 0, QTableWidgetItem(str(m['mode'])))
            aff = QTableWidgetItem(f"{m['affinity']:.2f}")
            aff.setForeground(QColor('#22C55E') if m['affinity'] < -7 else QColor('#F59E0B'))
            self.dock_results_table.setItem(i, 1, aff)
            self.dock_results_table.setItem(i, 2, QTableWidgetItem(f"{m['rmsd_lb']:.3f}"))
            self.dock_results_table.setItem(i, 3, QTableWidgetItem(f"{m['rmsd_ub']:.3f}"))
            self.dock_results_table.setItem(i, 4, QTableWidgetItem(f"{i + 1}"))
        self.dock_results_table.resizeColumnsToContents()

    def _dock_error(self, err_msg):
        self.dock_run_btn.setEnabled(True)
        self.dock_cancel_btn.setEnabled(False)
        self.dock_output.append(f"✗ Error: {err_msg}")

    def _dock_cancel(self):
        if hasattr(self, '_dock_worker'):
            self.dock_output.append("Cancelling...")
            self.dock_run_btn.setEnabled(True)
            self.dock_cancel_btn.setEnabled(False)

    def _dock_visualize(self):
        if not self.docking_results:
            QMessageBox.information(self, "No Results", "Run a docking calculation first.")
            return
        self._ensure_canvas()
        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        MplCanvas.style_ax(ax)

        modes = self.docking_results
        affinities = [m['affinity'] for m in modes]
        labels = [str(m['mode']) for m in modes]
        colors = ['#22C55E' if a < -7 else '#F59E0B' if a < -5 else '#EF4444' for a in affinities]

        bars = ax.bar(labels, affinities, color=colors, edgecolor='#1F2A40', linewidth=0.6)
        ax.axhline(y=-7, color='#22C55E', linestyle='--', alpha=0.5, label='Good (-7)')
        ax.axhline(y=-5, color='#F59E0B', linestyle='--', alpha=0.5, label='Moderate (-5)')

        for bar, aff in zip(bars, affinities):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f'{aff:.1f}', ha='center', va='bottom' if aff >= 0 else 'top',
                    fontsize=9, color='#E5E7EB')

        ax.set_xlabel('Mode', color='#8B9DC3')
        ax.set_ylabel('Affinity (kcal/mol)', color='#8B9DC3')
        ax.set_title('Docking Binding Affinities', color='#00D1FF', fontweight='bold')
        ax.legend(loc='upper left', facecolor='#121A2B', edgecolor='#1F2A40')
        ax.set_facecolor('#0D1528')

        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _dock_open_pymol(self):
        files = []
        receptor = self.dock_receptor_path.text()
        if receptor and os.path.exists(receptor):
            files.append(receptor)
        if self._last_docking_result and self._last_docking_result.get('out_path'):
            out = self._last_docking_result['out_path']
            if os.path.exists(out):
                files.append(out)
        if not files:
            QMessageBox.information(self, "Nothing to Open", "No receptor or docking results to show.")
            return
        ExternalToolLauncher.launch_pymol(filename=files)
        self.statusBar().showMessage(f"Opened {len(files)} file(s) in PyMOL", 3000)

    def _dock_open_all_poses(self):
        files = []
        receptor = self.dock_receptor_path.text()
        if receptor and os.path.exists(receptor):
            files.append(receptor)
        if self._last_docking_result and self._last_docking_result.get('out_path'):
            out = self._last_docking_result['out_path']
            if out and os.path.exists(out):
                if self.project_manager.is_open():
                    poses = self.project_manager.split_docking_poses(out)
                    files.extend(poses)
                else:
                    files.append(out)
        if not files:
            QMessageBox.information(self, "Nothing to Open", "No docking results to show.")
            return
        ExternalToolLauncher.launch_pymol(filename=files)
        self.statusBar().showMessage(f"Opened {len(files)} file(s) in PyMOL", 3000)

    def _dock_fetch_pdb(self):
        pdb_id = self.dock_pdb_id.text().strip().upper()
        if not pdb_id:
            QMessageBox.warning(self, "PDB ID Required", "Enter a PDB ID (e.g. 4HHB).")
            return
        self.dock_output.append(f"Fetching PDB {pdb_id}...")
        self.dock_fetch_pdb_btn.setEnabled(False)
        self.dock_fetch_pdb_btn.setText("Fetching...")
        QApplication.processEvents()

        self._dock_fetch_pdb_id = pdb_id
        self._dock_pdb_fetch_thread = PdbFetchThread(pdb_id)
        self._dock_pdb_fetch_thread.pdb_fetched.connect(self._dock_pdb_fetched)
        self._dock_pdb_fetch_thread.start()

    def _dock_pdb_fetched(self, pdb_id, data, error):
        self.dock_fetch_pdb_btn.setEnabled(True)
        self.dock_fetch_pdb_btn.setText("Download")
        if error:
            self.dock_output.append(f"✗ Failed to fetch PDB {pdb_id}: {error}")
            return
        dest = os.path.join(
            self.project_manager.project_dir if self.project_manager.is_open() else os.path.expanduser("~"),
            ".zora_docking", f"{pdb_id}.pdb")
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(dest, 'w') as f:
            f.write(data)
        self.dock_receptor_path.setText(dest)
        self.dock_output.append(f"✓ Downloaded {pdb_id} → {dest}")
        self.statusBar().showMessage(f"Loaded PDB {pdb_id}", 3000)

    def _dock_export_csv(self):
        if not self.docking_results:
            QMessageBox.information(self, "No Results", "Run a docking calculation first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", "docking_results.csv", "CSV (*.csv)")
        if not path:
            return
        import csv
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Mode", "Affinity (kcal/mol)", "RMSD LB", "RMSD UB"])
            for m in self.docking_results:
                writer.writerow([m['mode'], m['affinity'], m['rmsd_lb'], m['rmsd_ub']])
        self.statusBar().showMessage(f"Exported docking results to {path}", 3000)

    # -------- Public Methods ----------
    def _add_sequences(self, records: List[SequenceRecord]):
        existing_ids = {rec.id for rec in self.sequences}
        added = 0
        for rec in records:
            if rec.sequence and rec.id not in existing_ids:
                self.sequences.append(rec)
                existing_ids.add(rec.id)
                self._show_sequence(rec)
                self._update_info()
                self.statusBar().showMessage(f"Loaded: {rec.name} ({len(rec)} bp)", 5000)
                added += 1
        if added:
            self._update_project_tree()
            self._update_combos()
            if self.project_manager.is_open():
                self.project_manager.save(sequences=self.sequences)

    def _update_project_tree(self):
        self.project_tree.clear()
        # Sequence items
        seq_root = QTreeWidgetItem(["Sequences"])
        seq_root.setExpanded(True)
        for rec in self.sequences:
            item = QTreeWidgetItem([f"{rec.name} ({len(rec)} bp)"])
            item.setData(0, Qt.UserRole, rec.id)
            item.setData(0, Qt.UserRole + 1, 'seq')
            seq_root.addChild(item)
        self.project_tree.addTopLevelItem(seq_root)
        # PDB items (from project)
        pdb_root = QTreeWidgetItem(["Structures (PDB)"])
        pdb_root.setExpanded(True)
        if self.project_manager.is_open():
            for pf in self.project_manager.get_pdb_files():
                item = QTreeWidgetItem([pf['name']])
                item.setData(0, Qt.UserRole, pf['path'])
                item.setData(0, Qt.UserRole + 1, 'pdb')
                pdb_root.addChild(item)
        else:
            pdb_root.setToolTip(0, "Open a project to save and browse PDB files")
        self.project_tree.addTopLevelItem(pdb_root)
        # Docking results (from project)
        dock_root = QTreeWidgetItem(["Docking Results"])
        dock_root.setExpanded(True)
        if self.project_manager.is_open():
            for df in self.project_manager.get_docking_files():
                item = QTreeWidgetItem([df['name']])
                item.setData(0, Qt.UserRole, df['path'])
                item.setData(0, Qt.UserRole + 1, 'docking')
                dock_root.addChild(item)
            # Also list individual pose PDBs
            for df in self.project_manager.get_docking_files():
                base = df['name'].replace('.pdbqt', '')
                d_dir = os.path.dirname(df['path'])
                for fname in sorted(os.listdir(d_dir)):
                    if fname.startswith(base + '_pose_') and fname.endswith('.pdb'):
                        item = QTreeWidgetItem([f"  {fname}"])
                        item.setData(0, Qt.UserRole, os.path.join(d_dir, fname))
                        item.setData(0, Qt.UserRole + 1, 'pdb')
                        dock_root.addChild(item)
        else:
            dock_root.setToolTip(0, "Open a project to save docking results")
        self.project_tree.addTopLevelItem(dock_root)

    def _update_combos(self):
        names = [f"{r.name} ({len(r)} bp)" for r in self.sequences]
        for combo in [self.ref_seq_combo, self.mut_seq_combo, self.aln_seq1_combo,
                       self.aln_seq2_combo, self.batch_target_combo, self.crispr_seq_combo,
                       self.ai_target_combo, self.seq_viewer_combo, self.analysis_seq_combo]:
            combo.clear()
            combo.addItems(names)
        self.graph_seq_list.clear()
        for r in self.sequences:
            item = QListWidgetItem(f"{r.name} ({len(r)} bp)")
            item.setData(Qt.UserRole, self.sequences.index(r))
            item.setSelected(True)
            self.graph_seq_list.addItem(item)

    def _update_info(self):
        if self.sequences:
            rec = self.sequences[0]
            s = rec.stats()
            info = f"{rec.name} | {s['length']} bp | GC: {s['gc_content']:.1f}% | Tm: {s['melting_temp']:.1f}°C | MW: {s['molecular_weight']:.0f} Da"
            self.info_label.setText(info)
        else:
            self.info_label.setText("No sequence loaded")

    def _show_sequence(self, rec: SequenceRecord):
        self.seq_viewer.clear()
        formatted = ""
        seq = rec.sequence
        for i in range(0, len(seq), 60):
            line = seq[i:i+60]
            num = f"{i+1:>8} "
            formatted += num + line + "\n"
        self.seq_viewer.setPlainText(formatted)
        self.sequence_tabs.setCurrentIndex(0)

    def _on_seq_viewer_select(self, idx: int):
        if 0 <= idx < len(self.sequences):
            self._show_sequence(self.sequences[idx])

    def _view_mutations_in_sequence(self):
        if not self.mutation_results or not self.sequences:
            return
        ref_idx = self.ref_seq_combo.currentIndex()
        if ref_idx >= 0 and ref_idx < len(self.sequences):
            rec = self.sequences[ref_idx]
            self.seq_viewer.clear()
            seq = rec.sequence
            mut_positions = {r.position for r in self.mutation_results}
            formatted = ""
            for i in range(0, len(seq), 60):
                line = seq[i:i+60]
                num = f"{i+1:>8} "
                # Mark mutated bases with brackets
                marked = ""
                for j, base in enumerate(line):
                    pos = i + j
                    if pos in mut_positions:
                        marked += f"[{base}]"
                    else:
                        marked += f" {base} "
                formatted += num + marked + "\n"
            self.seq_viewer.setPlainText(formatted)
            self.annotation_text.setHtml(
                f"<b>Showing mutations in {rec.name}</b><br>"
                f"Found {len(self.mutation_results)} mutation(s). "
                f"Bases in [brackets] = mutated positions."
            )
            self.sequence_tabs.setCurrentIndex(0)

    # -------- Slots ----------
    def _open_file(self):
        files, _ = QFileDialog.getOpenFileName(self, "Open Sequence File",
                                                "", "Sequence Files (*.fasta *.fa *.fastq *.fq *.gb *.genbank *.txt *.csv);;All Files (*)")
        if files:
            records = FileParser.parse_file(files)
            if records:
                self._add_sequences(records)
            else:
                QMessageBox.warning(self, "Error", "Could not parse file. Unsupported format or empty.")

    def _paste_sequence(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Paste Sequence")
        dialog.setMinimumSize(600, 400)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel("Enter sequence (FASTA or raw):"))
        editor = QTextEdit()
        editor.setFont(QFont("Courier New", 11))
        layout.addWidget(editor)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        layout.addWidget(buttons)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)

        if dialog.exec():
            text = editor.toPlainText().strip()
            if not text:
                return
            if text.startswith('>'):
                records = FileParser.parse_fasta_from_text(text)
            else:
                seq = ''.join(c for c in text if c.isalpha()).upper()
                if seq:
                    rec = SequenceRecord("Pasted_Sequence", seq, "Manually pasted")
                    records = [rec]
                else:
                    records = []
            if records:
                self._add_sequences(records)

    def _clear_all(self):
        self.sequences.clear()
        self.mutation_results.clear()
        self.project_tree.clear()
        self.seq_viewer.clear()
        self.analysis_output.clear()
        self.mutation_output.clear()
        self.mutation_table.setRowCount(0)
        self.aln_output.clear()
        self.annotation_text.clear()
        self.viewer_seq_display.clear()
        self.viewer_info.clear()
        self.info_label.setText("No sequence loaded")
        self.statusBar().showMessage("All data cleared")

    def _run_analysis(self):
        if not self.sequences:
            QMessageBox.warning(self, "No Sequence", "Load a sequence first.")
            return
        idx = self.analysis_seq_combo.currentIndex()
        if idx < 0 or idx >= len(self.sequences):
            idx = 0
        rec = self.sequences[idx]
        seq = rec.sequence
        stats = rec.stats()
        orfs = SeqUtils.find_orfs(seq)
        rev_comp = SeqUtils.reverse_complement(seq)
        translation = SeqUtils.translate(seq)

        output = []
        output.append("=" * 60)
        output.append(f"SEQUENCE ANALYSIS REPORT: {rec.name}")
        output.append("=" * 60)
        output.append(f"Length: {stats['length']} bp")
        output.append(f"GC Content: {stats['gc_content']:.2f}%")
        output.append(f"GC Skew: {stats['gc_skew']:.4f}")
        output.append(f"Molecular Weight: {stats['molecular_weight']:.2f} Da")
        output.append(f"Melting Temperature: {stats['melting_temp']:.2f} °C")
        output.append("")
        output.append("Nucleotide Frequency:")
        for base, freq in stats['nucleotide_freq'].items():
            bar = '█' * int(freq / 2)
            output.append(f"  {base}: {freq:.2f}% {bar}")
        output.append("")
        output.append(f"Reverse Complement:")
        output.append(f"  {rev_comp[:100]}..." if len(rev_comp) > 100 else f"  {rev_comp}")
        output.append("")
        output.append("Translation (Frame 1):")
        prot = translation[:80]
        output.append(f"  {prot}")
        output.append("")
        output.append(f"Open Reading Frames (min 30bp): {len(orfs)}")
        for start, stop, protein, frame in orfs[:20]:
            output.append(f"  Frame {frame}: {start}-{stop} ({len(protein)} aa) : {protein[:50]}")
        if len(orfs) > 20:
            output.append(f"  ... and {len(orfs) - 20} more ORFs")
        output.append("")
        output.append("Codon Usage (top codons per amino acid):")
        cu = SeqUtils.codon_usage(seq)
        for aa in sorted(cu.keys()):
            if aa == '?':
                continue
            codons = cu[aa]['codons']
            if codons:
                best_codon = max(codons.items(), key=lambda x: x[1]['freq'])
                output.append(f"  {aa} ({AMINO_ACID_NAMES.get(aa, aa)}): best codon {best_codon[0]} ({best_codon[1]['freq']:.1f}%)")

        self.analysis_output.setPlainText('\n'.join(output))
        # Show send-to-tools bar if protein translation exists
        has_protein = bool(translation.strip('*?'))
        self.send_tools_bar.setVisible(has_protein)
        self.sequence_tabs.setCurrentWidget(self.analysis_tab)
        self.statusBar().showMessage("Analysis complete", 3000)

    def _detect_mutations(self):
        if len(self.sequences) < 2:
            QMessageBox.warning(self, "Need Two Sequences",
                                "Load at least two sequences (reference and mutated) to detect mutations.")
            return
        ref_idx = self.ref_seq_combo.currentIndex()
        mut_idx = self.mut_seq_combo.currentIndex()
        if ref_idx < 0 or mut_idx < 0:
            return
        ref_seq = self.sequences[ref_idx].sequence
        mut_seq = self.sequences[mut_idx].sequence

        results = MutationEngine.detect_mutations(ref_seq, mut_seq)
        self.mutation_results = results

        # Fill table
        self.mutation_table.setRowCount(len(results))
        for i, r in enumerate(results):
            self.mutation_table.setItem(i, 0, QTableWidgetItem(str(r.position)))
            self.mutation_table.setItem(i, 1, QTableWidgetItem(r.ref_base))
            self.mutation_table.setItem(i, 2, QTableWidgetItem(r.alt_base))
            self.mutation_table.setItem(i, 3, QTableWidgetItem(r.ref_codon))
            self.mutation_table.setItem(i, 4, QTableWidgetItem(r.alt_codon))
            self.mutation_table.setItem(i, 5, QTableWidgetItem(r.ref_aa))
            self.mutation_table.setItem(i, 6, QTableWidgetItem(r.alt_aa))
            self.mutation_table.setItem(i, 7, QTableWidgetItem(r.classification))

        # Summary
        summary = []
        summary.append(f"MUTATION REPORT")
        summary.append(f"Reference: {self.sequences[ref_idx].name}")
        summary.append(f"Mutated:   {self.sequences[mut_idx].name}")
        summary.append(f"Total differences: {len(results)}")
        types = Counter(r.classification for r in results)
        for t, count in types.most_common():
            summary.append(f"  {t}: {count}")
        if results:
            summary.append("")
            summary.append("Detailed Results (first 50):")
            for r in results[:50]:
                aa_info = ""
                if r.ref_aa and r.alt_aa:
                    aa_info = f" | {r.ref_codon}→{r.alt_codon} {r.ref_aa}→{r.alt_aa}"
                summary.append(f"  Pos {r.position}: {r.ref_base}→{r.alt_base} [{r.classification}]{aa_info}")
            if len(results) > 50:
                summary.append(f"  ... and {len(results) - 50} more")

        self.mutation_output.setPlainText('\n'.join(summary))
        self.sequence_tabs.setCurrentWidget(self.mutation_tab)
        self.statusBar().showMessage(f"Detected {len(results)} mutations", 3000)
        if self.project_manager.is_open():
            self.project_manager.save(sequences=self.sequences, mutation_results=self.mutation_results)

    def _mutation_selected(self):
        rows = self.mutation_table.selectedItems()
        if not rows:
            return
        row = rows[0].row()
        if row < len(self.mutation_results):
            r = self.mutation_results[row]
            aa_info = ""
            if r.ref_aa and r.alt_aa:
                ref_name = AMINO_ACID_NAMES.get(r.ref_aa, r.ref_aa)
                alt_name = AMINO_ACID_NAMES.get(r.alt_aa, r.alt_aa)
                aa_info = (f"\n\n<b>Amino Acid Change:</b>\n"
                           f"  Codon: {r.ref_codon} → {r.alt_codon}\n"
                           f"  Residue: {r.ref_aa} ({ref_name}) → {r.alt_aa} ({alt_name})")
            self.annotation_text.setHtml(
                f"<h3>Mutation at Position {r.position + 1}</h3>"
                f"<b>Type:</b> {r.classification}<br>"
                f"<b>Change:</b> {r.ref_base} → {r.alt_base}"
                f"{aa_info}"
            )

    def _simulate_mutation(self):
        if not self.sequences:
            return
        ref_idx = self.ref_seq_combo.currentIndex()
        if ref_idx < 0:
            ref_idx = 0
        seq = self.sequences[0].sequence
        pos = self.sim_pos_spin.value()
        alt = self.sim_base_combo.currentText()
        if pos < len(seq):
            mutated = MutationEngine.simulate_mutation(seq, pos, alt)
            rec = SequenceRecord(f"{self.sequences[0].name}_mut_{pos}{alt}", mutated,
                                 f"SNP at {pos}: {seq[pos]}→{alt}")
            self.sequences.append(rec)
            self._update_project_tree()
            self._update_combos()
            self._show_sequence(rec)
            self._update_info()
            self.statusBar().showMessage(f"Simulated mutation at position {pos}: {seq[pos]} → {alt}")

    def _simulate_insertion(self):
        if not self.sequences:
            return
        seq = self.sequences[0].sequence
        pos = self.ins_pos_spin.value()
        ins = self.ins_seq_edit.text().strip()
        if ins:
            mutated = MutationEngine.simulate_insertion(seq, pos, ins)
            rec = SequenceRecord(f"{self.sequences[0].name}_ins_{pos}", mutated,
                                 f"Insertion at {pos}: +{ins}")
            self.sequences.append(rec)
            self._update_project_tree()
            self._update_combos()
            self._show_sequence(rec)
            self._update_info()

    def _simulate_deletion(self):
        if not self.sequences:
            return
        seq = self.sequences[0].sequence
        start = self.del_start_spin.value()
        end = self.del_end_spin.value()
        if start < end <= len(seq):
            mutated = MutationEngine.simulate_deletion(seq, start, end)
            rec = SequenceRecord(f"{self.sequences[0].name}_del_{start}-{end}", mutated,
                                 f"Deletion {start}-{end} ({end-start} bp)")
            self.sequences.append(rec)
            self._update_project_tree()
            self._update_combos()
            self._show_sequence(rec)
            self._update_info()

    def _batch_mutations(self):
        if not self.sequences:
            return
        idx = self.batch_target_combo.currentIndex()
        if idx < 0:
            idx = 0
        seq = self.sequences[idx].sequence
        count = self.batch_count_spin.value()

        variants = MutationEngine.batch_generate_mutations(seq, count)
        added = 0
        for v in variants:
            if v != seq:
                rec = SequenceRecord(f"{self.sequences[idx].name}_variant_{added+1}", v)
                self.sequences.append(rec)
                added += 1
                if added > 500:
                    break
        self._update_project_tree()
        self._update_combos()
        self.statusBar().showMessage(f"Generated {added} variant sequences")

    def _find_motif_dialog(self):
        if not self.sequences:
            QMessageBox.warning(self, "No Sequence", "Load a sequence first.")
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Find Motif / Pattern")
        dialog.setMinimumWidth(500)
        layout = QVBoxLayout(dialog)

        form = QFormLayout()
        pattern_edit = QLineEdit()
        pattern_edit.setPlaceholderText("e.g. ATG, TATA, or regex pattern")
        form.addRow("Pattern:", pattern_edit)

        seq_combo = QComboBox()
        for r in self.sequences:
            seq_combo.addItem(f"{r.name} ({len(r)} bp)")
        form.addRow("In sequence:", seq_combo)

        regex_cb = QCheckBox("Use regex")
        form.addRow("", regex_cb)

        case_cb = QCheckBox("Case sensitive")
        form.addRow("", case_cb)

        layout.addLayout(form)

        result_text = QTextEdit()
        result_text.setReadOnly(True)
        result_text.setFont(QFont("Courier New", 10))
        layout.addWidget(result_text)

        buttons = QHBoxLayout()
        search_btn = QPushButton("Search")
        search_btn.setStyleSheet("background: #00D1FF; color: #0B1020; font-weight: bold; padding: 6px; border: none; border-radius: 4px;")
        buttons.addWidget(search_btn)
        buttons.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.accept)
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)

        def do_search():
            pattern = pattern_edit.text()
            if not pattern:
                return
            idx = seq_combo.currentIndex()
            rec = self.sequences[idx]
            seq = rec.sequence if case_cb.isChecked() else rec.sequence.upper()
            pat = pattern if case_cb.isChecked() else pattern.upper()
            matches = SeqUtils.find_motifs(seq, pat, regex_cb.isChecked())
            if not matches:
                result_text.setPlainText(f"No matches found for '{pattern}' in {rec.name}")
                return
            lines = [f"Found {len(matches)} match(es) for '{pattern}' in {rec.name}:",
                     f"{'Pos':>8}  {'Match':<20}  Context"]
            lines.append("-" * 60)
            for pos, m in matches[:200]:
                ctx_start = max(0, pos - 5)
                ctx_end = min(len(seq), pos + len(m) + 5)
                context = seq[ctx_start:pos] + '[' + seq[pos:pos+len(m)] + ']' + seq[pos+len(m):ctx_end]
                lines.append(f"{pos+1:>8}  {m:<20}  ...{context}...")
            if len(matches) > 200:
                lines.append(f"\n... and {len(matches) - 200} more matches")
            result_text.setPlainText('\n'.join(lines))

        search_btn.clicked.connect(do_search)
        dialog.exec()

    def _show_revcomp(self):
        if not self.sequences:
            return
        rec = self.sequences[0]
        rc = SeqUtils.reverse_complement(rec.sequence)
        comp = SeqUtils.complement(rec.sequence)
        text = (f"Original ({rec.name}):\n  {rec.sequence[:200]}{'...' if len(rec.sequence) > 200 else ''}\n\n"
                f"Complement:\n  {comp[:200]}{'...' if len(comp) > 200 else ''}\n\n"
                f"Reverse Complement:\n  {rc[:200]}{'...' if len(rc) > 200 else ''}")
        # Create new record for revcomp
        rc_rec = SequenceRecord(f"{rec.name}_revcomp", rc, "Reverse complement")
        self.sequences.append(rc_rec)
        self._update_project_tree()
        self._update_combos()
        self._show_sequence(rc_rec)
        self._update_info()
        self.statusBar().showMessage("Reverse complement added to project", 3000)

    def _show_orfs(self):
        if not self.sequences:
            return
        rec = self.sequences[0]
        orfs = SeqUtils.find_orfs(rec.sequence)
        if not orfs:
            QMessageBox.information(self, "ORFs", "No open reading frames found (min 30bp).")
            return

        text = (f"Open Reading Frames in {rec.name}\n"
                f"{'Frame':>8} {'Start':>8} {'Stop':>8} {'Length(aa)':>12} {'Sequence'}")
        text += "\n" + "-" * 70
        for start, stop, protein, frame in sorted(orfs, key=lambda x: len(x[2]), reverse=True)[:50]:
            text += f"\n{frame:>8} {start+1:>8} {stop:>8} {len(protein):>8} aa  {protein[:40]}"
        if len(orfs) > 50:
            text += f"\n\n... and {len(orfs) - 50} more ORFs"

        self.analysis_output.setPlainText(text)
        self.sequence_tabs.setCurrentWidget(self.analysis_tab)

    def _annotation_lookup(self):
        if not self.mutation_results:
            QMessageBox.warning(self, "No Mutations", "Detect mutations first to look up annotations.")
            return
        QMessageBox.information(self, "Annotation Lookup",
            "<h3>Annotation Sources</h3>"
            "<p>This feature integrates with:</p>"
            "<ul>"
            "<li><b>ClinVar</b> - https://www.ncbi.nlm.nih.gov/clinvar/</li>"
            "<li><b>dbSNP</b> - https://www.ncbi.nlm.nih.gov/snp/</li>"
            "<li><b>Ensembl</b> - https://www.ensembl.org/</li>"
            "</ul>"
            "<p><b>How to use:</b></p>"
            "<ol>"
            "<li>Select a mutation in the mutation table</li>"
            "<li>Look up the position in dbSNP/ClinVar using the links above</li>"
            "<li>OR select 'Open in Browser' to search automatically</li>"
            "</ol>"
            "<p><i>Full API integration coming in ZORA v2.0</i></p>"
        )

    def _run_alignment(self, method: str = 'nw'):
        if len(self.sequences) < 2:
            QMessageBox.warning(self, "Need Two Sequences", "Load at least two sequences for alignment.")
            return
        idx1 = self.aln_seq1_combo.currentIndex()
        idx2 = self.aln_seq2_combo.currentIndex()
        if idx1 < 0 or idx2 < 0:
            return
        seq1 = self.sequences[idx1].sequence[:2000]
        seq2 = self.sequences[idx2].sequence[:2000]
        name1 = self.sequences[idx1].name
        name2 = self.sequences[idx2].name

        self.statusBar().showMessage(f"Running {method} alignment...")
        worker = AlignWorker(seq1, seq2, method, name1, name2)
        worker.signals.finished.connect(lambda a1, a2, sc, mn, n1, n2: self._alignment_done(a1, a2, sc, mn, n1, n2))
        worker.signals.error.connect(lambda err: self.statusBar().showMessage(f"Alignment error: {err}", 5000))
        QThreadPool.globalInstance().start(worker)

    def _alignment_done(self, aln1, aln2, score, method_name, name1, name2):
        output = []
        output.append(f"ALIGNMENT: {method_name}")
        output.append(f"Score: {score}")
        output.append(f"Sequence 1: {name1}")
        output.append(f"Sequence 2: {name2}")
        output.append(f"Alignment length: {len(aln1)}")
        output.append("")
        match_str = ''
        for a, b in zip(aln1, aln2):
            if a == b:
                match_str += '|'
            elif a == '-' or b == '-':
                match_str += ' '
            else:
                match_str += '.'
        for i in range(0, len(aln1), 80):
            end = min(i + 80, len(aln1))
            output.append(f"Q: {aln1[i:end]}")
            output.append(f"   {match_str[i:end]}")
            output.append(f"R: {aln2[i:end]}")
            output.append("")
        self.aln_output.setPlainText('\n'.join(output))
        self.sequence_tabs.setCurrentWidget(self.alignment_tab)
        self.statusBar().showMessage(f"Alignment score: {score}", 3000)

    def _get_graph_seqs(self) -> list:
        idxs = [item.data(Qt.UserRole) for item in self.graph_seq_list.selectedItems()]
        return [self.sequences[i] for i in idxs if 0 <= i < len(self.sequences)]

    def _plot_gc_content(self):
        self._ensure_canvas()
        seqs = self._get_graph_seqs()
        if not seqs:
            QMessageBox.information(self, "No Selection", "Select sequences in the list above.")
            return
        window = self.gc_window_spin.value()
        step = self.gc_step_spin.value()
        theme = ThemeManager.THEMES.get(self.current_theme, {})
        colors = ['#22C55E', '#EF4444', '#00D1FF', '#8B5CF6', '#F59E0B', '#F472B6']

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        for idx, rec in enumerate(seqs):
            seq = rec.sequence
            key = (id(seq), window, step)
            cached = self._gc_cache.get(key)
            if cached:
                positions, gc_values = cached
            else:
                gc_values = []
                positions = []
                for i in range(0, max(len(seq) - window + 1, 1), step):
                    frag = seq[i:i+window]
                    gc = SeqUtils.gc_content(frag)
                    gc_values.append(gc)
                    positions.append(i + window // 2)
                if len(self._gc_cache) > 100:
                    self._gc_cache.pop(next(iter(self._gc_cache)))
                self._gc_cache[key] = (positions, gc_values)
            color = colors[idx % len(colors)]
            ax.plot(positions, gc_values, color=color, linewidth=1.8, label=rec.name, alpha=0.85)
            ax.fill_between(positions, 0, gc_values, color=color, alpha=0.08)

        ax.axhline(y=50, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('GC Content (%)')
        ax.set_title(f'GC Content — sliding window ({window} bp, step {step})')
        ax.set_ylim(0, 100)
        ax.legend(frameon=True, fontsize=9, loc='upper right')
        MplCanvas.style_ax(ax, theme, grid=self._graph_grid_visible)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_mutation_distribution(self):
        self._ensure_canvas()
        if not self.mutation_results:
            QMessageBox.warning(self, "No Mutations", "Run mutation detection first.")
            return
        positions = [r.position for r in self.mutation_results]
        theme = ThemeManager.THEMES.get(self.current_theme, {})

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        n, bins, patches = ax.hist(positions, bins=50, color='#EF4444', alpha=0.6,
                                    edgecolor='white', linewidth=0.5)
        for patch in patches:
            patch.set_facecolor(plt.cm.Reds(patch.get_height() / max(n + 1e-9)))
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Mutation Count')
        ax.set_title(f'Mutation Distribution ({len(positions)} total)')
        # Stats annotation
        mean_pos = np.mean(positions)
        median_pos = np.median(positions)
        ax.axvline(mean_pos, color='#00D1FF', linestyle='--', linewidth=1.2, alpha=0.7, label=f'Mean: {mean_pos:.0f}')
        ax.axvline(median_pos, color='#22C55E', linestyle=':', linewidth=1.2, alpha=0.7, label=f'Median: {median_pos:.0f}')
        ax.legend(fontsize=8, frameon=True)
        # Stats box
        stats_text = f'Total: {len(positions)}\nMean: {mean_pos:.1f}\nMedian: {median_pos:.0f}\nRange: {min(positions)}–{max(positions)}'
        ax.text(0.02, 0.97, stats_text, transform=ax.transAxes, fontsize=7,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor=theme.get('BG_PANEL', '#121A2B'), alpha=0.8))
        MplCanvas.style_ax(ax, theme, grid=self._graph_grid_visible)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_nucleotide_frequency(self):
        self._ensure_canvas()
        seqs = self._get_graph_seqs()
        if not seqs:
            QMessageBox.information(self, "No Selection", "Select sequences in the list above.")
            return
        theme = ThemeManager.THEMES.get(self.current_theme, {})
        colors = ['#22C55E', '#EF4444', '#2563EB', '#F59E0B', '#8B5CF6', '#00BCD4']

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        bases = ['A', 'T', 'G', 'C']
        n_seq = len(seqs)
        bar_width = 0.7 / n_seq
        for idx, rec in enumerate(seqs):
            freqs = SeqUtils.nucleotide_frequency(rec.sequence)
            values = [freqs.get(b, 0) for b in bases]
            x = [i + idx * bar_width for i in range(len(bases))]
            bars = ax.bar(x, values, width=bar_width, color=colors[idx % len(colors)],
                          edgecolor='white', linewidth=0.5, label=rec.name)
            for bar, v in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        f'{v:.1f}%', ha='center', fontsize=7, fontweight='bold')
        ax.set_xticks([i + bar_width * (n_seq - 1) / 2 for i in range(len(bases))])
        ax.set_xticklabels(bases)
        ax.set_ylabel('Frequency (%)')
        ax.set_title('Nucleotide Frequency Comparison')
        ax.legend(fontsize=8, frameon=True)
        ax.set_ylim(0, 70)
        MplCanvas.style_ax(ax, theme, grid=self._graph_grid_visible)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_codon_usage(self):
        self._ensure_canvas()
        seqs = self._get_graph_seqs()
        if not seqs:
            QMessageBox.information(self, "No Selection", "Select sequences in the list above.")
            return
        theme = ThemeManager.THEMES.get(self.current_theme, {})
        colors = ['#22C55E', '#EF4444', '#00D1FF', '#8B5CF6', '#F59E0B', '#F472B6']

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        n_seq = len(seqs)
        bar_height = 0.8 / n_seq
        all_labels = []
        for idx, rec in enumerate(seqs):
            seq = rec.sequence
            codons = Counter()
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if len(codon) == 3 and all(c in 'ATCGatcg' for c in codon):
                    codons[codon.upper()] += 1
            total = sum(codons.values()) or 1
            items = codons.most_common(15)
            labels = [f"{c[0]} ({AMINO_ACIDS.get(c[0], '?')})" for c in items]
            values = [c[1] / total * 100 for c in items]
            if idx == 0:
                all_labels = labels
            y = [i + idx * bar_height for i in range(len(items))]
            bars = ax.barh(y, values, height=bar_height, color=colors[idx % len(colors)],
                           edgecolor='white', linewidth=0.3, label=rec.name, alpha=0.85)
            for v, yi in zip(values, y):
                if v > 1:
                    ax.text(v + 0.3, yi, f'{v:.1f}%', va='center', fontsize=6.5, fontweight='bold')
        ax.set_yticks([i + bar_height * (n_seq - 1) / 2 for i in range(len(all_labels))])
        ax.set_yticklabels(all_labels, fontsize=7)
        ax.set_xlabel('Frequency (%)')
        ax.set_title('Codon Usage Comparison (top 15)')
        ax.legend(fontsize=7, frameon=True, loc='lower right')
        MplCanvas.style_ax(ax, theme, grid=self._graph_grid_visible)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_mutation_pie(self):
        self._ensure_canvas()
        if not self.mutation_results:
            QMessageBox.warning(self, "No Mutations", "Run mutation detection first.")
            return
        types = Counter(r.classification for r in self.mutation_results)
        theme = ThemeManager.THEMES.get(self.current_theme, {})

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        labels = list(types.keys())
        sizes = list(types.values())
        pie_colors = ['#EF4444', '#22C55E', '#00D1FF', '#8B5CF6', '#F59E0B', '#F472B6', '#6B7280']
        explode = [0.05] * len(labels)
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                                            colors=pie_colors[:len(labels)], startangle=90,
                                            explode=explode, pctdistance=0.78,
                                            wedgeprops=dict(width=0.35, edgecolor=theme.get('BG_DARK', '#0B1020'), linewidth=2))
        for t in autotexts:
            t.set_fontsize(9)
            t.set_fontweight('bold')
        for t in texts:
            t.set_fontsize(9)
        ax.set_title('Mutation Type Distribution (donut)', pad=15)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_snp_heatmap(self):
        self._ensure_canvas()
        if not self.mutation_results:
            QMessageBox.warning(self, "No Mutations", "Run mutation detection first.")
            return
        positions = [r.position for r in self.mutation_results]
        if not positions:
            return
        theme = ThemeManager.THEMES.get(self.current_theme, {})

        # Proper 1D density heatmap
        if len(set(positions)) > 1:
            bins = min(100, len(set(positions)))
            hist, edges = np.histogram(positions, bins=bins)
            data_2d = hist.reshape(1, -1)
            extent = [edges[0], edges[-1], 0, 1]
        else:
            data_2d = np.array([[1]])
            extent = [positions[0] - 1, positions[0] + 1, 0, 1]

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)
        cmap = 'hot_r' if self.current_theme == 'Clinical Clean White' else 'hot'
        im = ax.imshow(data_2d, cmap=cmap, aspect='auto', interpolation='bilinear', extent=extent)
        cb = self.graph_canvas.fig.colorbar(im, ax=ax, label='SNP density', pad=0.02, shrink=0.8)
        cb.ax.tick_params(labelsize=8)
        ax.set_title('SNP Density Heatmap')
        ax.set_xlabel('Position (bp)')
        ax.set_yticks([])
        MplCanvas.style_ax(ax, theme, grid=False)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _plot_conservation(self):
        self._ensure_canvas()
        seqs = self._get_graph_seqs()
        if len(seqs) < 2:
            QMessageBox.warning(self, "Need Two Sequences", "Select at least two sequences for conservation plot.")
            return
        theme = ThemeManager.THEMES.get(self.current_theme, {})

        names = [s.name for s in seqs]
        raw_seqs = [s.sequence.upper()[:5000] for s in seqs]
        min_len = min(len(s) for s in raw_seqs)
        n_seq = len(raw_seqs)

        window = 30
        positions = list(range(0, min_len - window + 1))
        pos_center = [p + window // 2 for p in positions]

        # Pairwise conservation matrix
        n_pairs = n_seq * (n_seq - 1) // 2
        conservation_matrix = np.zeros((n_seq, len(positions)))

        for si in range(n_seq):
            for pi, p in enumerate(positions):
                identical = sum(1 for j in range(n_seq) if si != j and
                                raw_seqs[si][p:p+window] == raw_seqs[j][p:p+window])
                conservation_matrix[si, pi] = identical / (n_seq - 1) * 100

        self.graph_canvas.fig.clear()
        ax = self.graph_canvas.fig.add_subplot(111)

        mean_conservation = np.mean(conservation_matrix, axis=0)
        std_conservation = np.std(conservation_matrix, axis=0)
        ax.fill_between(pos_center, mean_conservation - std_conservation,
                         mean_conservation + std_conservation, alpha=0.15, color='#00D1FF')
        ax.plot(pos_center, mean_conservation, color='#00D1FF', linewidth=2, label='Mean conservation')
        for si in range(n_seq):
            color = ['#22C55E', '#EF4444', '#8B5CF6', '#F59E0B'][si % 4]
            ax.plot(pos_center, conservation_matrix[si], color=color, linewidth=0.8, alpha=0.4, label=f'{names[si]}')
        ax.axhline(y=50, color='gray', linestyle='--', alpha=0.3, linewidth=0.8)
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Conservation (%)')
        ax.set_title(f'Multi-Sequence Conservation (window={window})')
        ax.set_ylim(-5, 105)
        ax.legend(fontsize=7, frameon=True, loc='lower left', ncol=2)
        MplCanvas.style_ax(ax, theme, grid=self._graph_grid_visible)
        self.graph_canvas.draw()
        self.sequence_tabs.setCurrentWidget(self.graphs_tab)

    def _toggle_graph_grid(self):
        if not hasattr(self, 'graph_canvas') or self.graph_canvas is None:
            return
        self._graph_grid_visible = not self._graph_grid_visible
        for ax in self.graph_canvas.fig.axes:
            if self._graph_grid_visible:
                ax.grid(True, alpha=0.15, linestyle='--', linewidth=0.5)
            else:
                ax.grid(False)
        self.graph_canvas.draw()

    def _export_svg(self):
        if not hasattr(self, 'graph_canvas') or self.graph_canvas is None:
            QMessageBox.information(self, "No Graph", "Generate a graph first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save SVG", "", "SVG Files (*.svg)")
        if not path:
            return
        self.graph_canvas.fig.savefig(path, format='svg', dpi=300, bbox_inches='tight')
        self.statusBar().showMessage(f"Graph saved as SVG: {path}", 3000)

    def _export_png(self):
        if not hasattr(self, 'graph_canvas') or self.graph_canvas is None:
            QMessageBox.information(self, "No Graph", "Generate a graph first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save PNG", "", "PNG Files (*.png)")
        if not path:
            return
        self.graph_canvas.fig.savefig(path, format='png', dpi=300, bbox_inches='tight')
        self.statusBar().showMessage(f"Graph saved as PNG: {path}", 3000)

    def _viewer_goto(self):
        pos = self.viewer_pos_spin.value()
        if not self.sequences:
            return
        seq = self.sequences[0].sequence
        if pos >= len(seq):
            return
        start = max(0, pos - 75)
        end = min(len(seq), pos + 75)
        fragment = seq[start:end]

        display = []
        for i, base in enumerate(fragment):
            if (start + i) % 60 == 0:
                display.append(f"\n{start + i + 1:>8} ")
            # Mark mutations with brackets
            is_mut = any(r.position == start + i for r in self.mutation_results)
            if is_mut:
                display.append(f"[{base}]")
            else:
                display.append(f" {base} ")
        self.viewer_seq_display.setPlainText(''.join(display))

        pos_info = []
        for offset in range(-2, 3):
            check_pos = pos + offset
            if 0 <= check_pos < len(seq):
                codon_start = (check_pos // 3) * 3
                if codon_start + 2 < len(seq):
                    codon = seq[codon_start:codon_start+3]
                    aa = AMINO_ACIDS.get(codon, '?')
                    aa_name = AMINO_ACID_NAMES.get(aa, 'Unknown')
                    mark = " ←" if check_pos == pos else ""
                    pos_info.append(f"  Base {check_pos + 1}: {seq[check_pos]} → Codon {codon_start+1}-{codon_start+3}: {codon} → {aa} ({aa_name}){mark}")
        self.viewer_info.setPlainText('\n'.join(pos_info) if pos_info else "Position info not available")

    def _viewer_cursor_changed(self):
        cursor = self.viewer_seq_display.textCursor()
        text = self.viewer_seq_display.toPlainText()
        if not text or not self.sequences:
            return
        char_pos = cursor.position()
        lines = text.split('\n')
        cum_len = 0
        for line in lines:
            line_clean = ''.join(c for c in line if c.isalpha())
            if char_pos <= cum_len + len(line_clean):
                local_pos = char_pos - cum_len
                if 0 <= local_pos < len(line_clean):
                    global_pos = local_pos + cum_len
                    if global_pos < len(self.sequences[0].sequence):
                        self.viewer_pos_spin.setValue(global_pos)
                        self._viewer_goto()
                break
            cum_len += len(line_clean) + 1

    def _design_crispr(self):
        if not self.sequences:
            return
        idx = self.crispr_seq_combo.currentIndex()
        if idx < 0:
            idx = 0
        seq = self.sequences[idx].sequence

        guides = []
        pam = 'NGG'
        # Find PAM sites (NGG)
        for i in range(len(seq) - 22):
            if seq[i+20:i+22].upper() == 'GG':
                guide = seq[i:i+20]
                if set(guide).issubset('ATCGatcg'):
                    guides.append({
                        'seq': guide,
                        'pam_pos': i + 20,
                        'strand': '+',
                        'gc': SeqUtils.gc_content(guide),
                        'offtargets': 0
                    })
            # Reverse strand
            rev = SeqUtils.reverse_complement(seq)
            if rev[i+20:i+22].upper() == 'GG':
                guide = rev[i:i+20]
                if set(guide).issubset('ATCGatcg'):
                    guides.append({
                        'seq': guide,
                        'pam_pos': len(seq) - i - 20,
                        'strand': '-',
                        'gc': SeqUtils.gc_content(guide),
                        'offtargets': 0
                    })

        if not guides:
            self.crispr_output.setPlainText("No SpCas9 (NGG) PAM sites found in this sequence.")
            return

        # Rank by GC content (best: 40-60%)
        def score(g):
            gc = g['gc']
            if 40 <= gc <= 60:
                return 100 - abs(gc - 50)
            return 0
        guides.sort(key=score, reverse=True)

        output = []
        output.append(f"CRISPR sgRNA Design for {self.sequences[idx].name}")
        output.append(f"PAM: {pam} (SpCas9)")
        output.append(f"Total guides found: {len(guides)}")
        output.append("")
        output.append("Top 50 Guides (ranked):")
        output.append(f"{'#':>4} {'Guide (20bp)':<22} {'PAM Pos':>8} {'Str':>4} {'GC%':>6} {'Score':>6}")
        output.append("-" * 60)
        for i, g in enumerate(guides[:50]):
            output.append(f"{i+1:>4} {g['seq']:<22} {g['pam_pos']:>8} {g['strand']:>4} {g['gc']:>5.1f}% {score(g):>6.1f}")

        self.crispr_output.setPlainText('\n'.join(output))
        self.sequence_tabs.setCurrentWidget(self.crispr_tab)

    def _run_ai_analysis(self):
        if not self.sequences:
            return
        self.ai_output.setPlainText(
            "AI Pathogenicity Prediction Module\n"
            "===============================\n\n"
            "This module will predict pathogenicity of mutations using ML models.\n\n"
            "STATUS: Placeholder\n\n"
            "To enable this feature:\n"
            "  1. Train a model (e.g., using ClinVar data)\n"
            "  2. Load the model weights\n"
            "  3. Run predictions on detected mutations\n\n"
            "Planned features:\n"
            "  - Pathogenicity scoring (benign/pathogenic)\n"
            "  - Mutation impact score\n"
            "  - SHAP explainability\n"
            "  - Sequence classification (cancer vs normal)\n"
            "  - Mutation hotspot clustering\n"
        )
        self.sequence_tabs.setCurrentWidget(self.ai_tab)

    # -------- Integrations ----------
    def _fetch_from_ncbi(self):
        accession = self.db_accession_edit.text().strip()
        if not accession:
            QMessageBox.warning(self, "No Accession", "Enter an NCBI accession number (e.g. NM_007294).")
            return
        self.db_output.setPlainText(f"Fetching {accession} from NCBI...")
        self.db_fetch_btn.setEnabled(False)
        self.db_fetch_btn.setText("Fetching...")
        self._ncbi_thread = NcbiFetchThread(accession)
        self._ncbi_thread.ncbi_fetched.connect(self._on_ncbi_fetched)
        self._ncbi_thread.start()

    def _on_ncbi_fetched(self, accession: str, rec_json: str, error: str):
        self.db_fetch_btn.setEnabled(True)
        self.db_fetch_btn.setText("Fetch from NCBI")
        if error:
            self.db_output.setPlainText(f"✗ Failed to fetch {accession}: {error}")
        else:
            import json
            rec_data = json.loads(rec_json)
            rec = SequenceRecord(rec_data['name'], rec_data['sequence'], rec_data['description'])
            self._add_sequences([rec])
            self.db_output.setPlainText(f"✓ Successfully fetched {accession}\n  {rec.name} | {len(rec)} bp")
        self._ncbi_thread = None

    def _search_ncbi_dialog(self):
        query = self.db_accession_edit.text().strip()
        if not query:
            QMessageBox.warning(self, "No Query", "Enter a search term (e.g. BRCA1, TP53, insulin).")
            return
        self.db_output.setPlainText(f"Searching NCBI for '{query}'...")
        QApplication.processEvents()
        results = NCBIFetcher.search_ncbi(query)
        if results:
            lines = [f"Found {len(results)} results for '{query}':"]
            for r in results[:20]:
                lines.append(f"  • {r['id']} ({r.get('accession', 'N/A')})")
            self.db_output.setPlainText('\n'.join(lines))
        else:
            self.db_output.setPlainText(f"No results found for '{query}'.")

    def _import_pdb(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB Files (*.pdb *.ent)")
        if not path:
            return
        fname = os.path.splitext(os.path.basename(path))[0]
        if self.project_manager.is_open():
            import shutil
            dest = self.project_manager.get_path('pdb', f"{fname}.pdb")
            shutil.copy2(path, dest)
            self._refresh_pdb_list()
            self._update_project_tree()
            self.statusBar().showMessage(f"Copied PDB to project: {fname}.pdb", 3000)
        else:
            # No project open — offer to start one or load as sequence
            reply = QMessageBox.question(self, "No Project Open",
                "PDB files are saved to a project.\n"
                "Open a project first (File → New Project) to store PDB files locally.\n\n"
                "Load sequences from this PDB into ZORA instead?",
                QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                records = PDBParser.parse(path)
                if records:
                    self._add_sequences(records)
                    self.statusBar().showMessage(f"Imported {len(records)} chains from PDB", 3000)
                else:
                    QMessageBox.warning(self, "No Data", "Could not extract sequences from PDB file.")

    def _fetch_pdb_rcsb(self):
        accession = self.db_accession_edit.text().strip().upper()
        if not accession:
            QMessageBox.warning(self, "No PDB ID", "Enter a PDB ID (e.g. 1CRN, 4HHB).")
            return
        if getattr(self, '_fetch_thread', None) is not None and self._fetch_thread.isRunning():
            QMessageBox.information(self, "Already Fetching", "A fetch is already in progress.")
            return
        self.db_output.setPlainText(f"Fetching {accession}...")
        self.pdb_fetch_btn.setEnabled(False)
        self.pdb_fetch_btn.setText("Fetching...")
        # Run in background thread so GUI doesn't freeze
        self._fetch_thread = PdbFetchThread(accession)
        self._fetch_thread.pdb_fetched.connect(self._on_pdb_fetched)
        self._fetch_thread.start()

    def _on_pdb_fetched(self, accession: str, data: str, error: str):
        self.pdb_fetch_btn.setEnabled(True)
        self.pdb_fetch_btn.setText("Fetch PDB from RCSB")
        if error:
            self.db_output.setPlainText(f"✗ Failed to fetch {accession}: {error}")
            self._fetch_thread = None
            return
        if self.project_manager.is_open():
            self.project_manager.save_pdb(accession, data)
            self._refresh_pdb_list()
            self._update_project_tree()
            self.db_output.setPlainText(f"✓ Saved PDB {accession} to project ({len(data):,} bytes)")
        else:
            self.db_output.setPlainText(
                f"✓ Fetched PDB {accession} ({len(data):,} bytes)\n"
                "  Open a project first to save PDB files locally.")
        self._fetch_thread = None

    def _launch_pymol(self):
        success, err = ExternalToolLauncher.launch_pymol()
        if err:
            QMessageBox.warning(self, "Launch Failed", err)
        else:
            self.statusBar().showMessage("PyMOL launched", 3000)

    def _launch_chimerax(self):
        success, err = ExternalToolLauncher.launch_chimerax()
        if err:
            QMessageBox.warning(self, "Launch Failed", err)
        else:
            self.statusBar().showMessage("ChimeraX launched", 3000)

    def _launch_snapgene(self):
        success, err = ExternalToolLauncher.launch_snapgene()
        if err:
            QMessageBox.warning(self, "Launch Failed", err)
        else:
            self.statusBar().showMessage("SnapGene launched", 3000)

    def _clean_seq_for_tool(self, rec: SequenceRecord) -> Optional[str]:
        STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')
        VALID_DNA = set('ATCGatcg')
        seq = rec.sequence.upper()

        # Check if it's DNA
        if set(seq).issubset(VALID_DNA):
            reply = QMessageBox.question(self, "DNA Sequence Detected",
                f"'{rec.name}' appears to be DNA ({len(seq)} bp, contains only A/T/C/G).\n\n"
                "Would you like to translate it to protein before sending?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return None
            if reply == QMessageBox.Yes:
                seq = SeqUtils.translate(seq)
                # Clean trailing stop codons for structure building
                seq = seq.replace('*', '')
                QMessageBox.information(self, "Translated",
                    f"Translated to protein ({len(seq)} aa):\n{seq[:60]}{'...' if len(seq) > 60 else ''}")

        # Check for non-standard AA
        bad = set(seq) - STANDARD_AA - set('*X?')
        if bad:
            AA_MAP = {
                'U': 'C',  # selenocysteine → cysteine
                'O': 'K',  # pyrrolysine → lysine
                'B': 'N',  # ambiguous Asx → asparagine
                'Z': 'Q',  # ambiguous Glx → glutamine
                'J': 'L',  # ambiguous Xle → leucine
            }
            mapped = {b: AA_MAP.get(b, 'A') for b in bad}
            map_desc = ', '.join(f"{b}→{m}" for b, m in mapped.items())
            reply = QMessageBox.warning(self, "Non-standard Amino Acids",
                f"Sequence contains non-standard residue(s): {', '.join(sorted(bad))}\n\n"
                "These will cause errors in PyMOL/ChimeraX structure building.\n\n"
                f"Recommended mapping: {map_desc}\n"
                "Apply mapping and send?",
                QMessageBox.Yes | QMessageBox.Cancel)
            if reply != QMessageBox.Yes:
                return None
            cleaned = []
            for c in seq:
                if c in AA_MAP:
                    cleaned.append(AA_MAP[c])
                elif c in STANDARD_AA:
                    cleaned.append(c)
                elif c == 'X':
                    cleaned.append('A')
            seq = ''.join(cleaned)
        # Replace X→A and remove stops/unknowns regardless of cleaning
        seq = seq.replace('X', 'A').replace('*', '').replace('?', '')
        return seq

    def _send_seq_to_pymol(self):
        if not self.sequences:
            return
        rec = self.sequences[0]
        seq = self._clean_seq_for_tool(rec)
        if seq is None:
            return
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        tmp.write(f">{rec.name}\n{seq}\n")
        tmp.close()
        success, err = ExternalToolLauncher.launch_pymol(filename=tmp.name)
        if err:
            QMessageBox.warning(self, "Send Failed", err)
        else:
            self.statusBar().showMessage(f"Sent {rec.name} ({len(seq)} aa) to PyMOL", 3000)

    def _send_seq_to_chimerax(self):
        if not self.sequences:
            return
        rec = self.sequences[0]
        seq = self._clean_seq_for_tool(rec)
        if seq is None:
            return
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        tmp.write(f">{rec.name}\n{seq}\n")
        tmp.close()
        success, err = ExternalToolLauncher.launch_chimerax(filename=tmp.name)
        if err:
            QMessageBox.warning(self, "Send Failed", err)
        else:
            self.statusBar().showMessage(f"Sent {rec.name} ({len(seq)} aa) to ChimeraX", 3000)

    def _send_seq_to_snapgene(self):
        if not self.sequences:
            return
        rec = self.sequences[0]
        seq = self._clean_seq_for_tool(rec)
        if seq is None:
            return
        import tempfile
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        tmp.write(f">{rec.name}\n{seq}\n")
        tmp.close()
        success, err = ExternalToolLauncher.launch_snapgene(filename=tmp.name)
        if err:
            QMessageBox.warning(self, "Send Failed", err)
        else:
            self.statusBar().showMessage(f"Sent {rec.name} ({len(seq)} aa) to SnapGene", 3000)

    def _export_for_tools(self):
        if not self.sequences:
            QMessageBox.warning(self, "No Sequences", "Load sequences first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export FASTA", "", "FASTA Files (*.fasta)")
        if not path:
            return
        ExternalToolLauncher.export_for_chimerax(self.sequences, path)
        self.statusBar().showMessage(f"Exported {len(self.sequences)} sequences to {path}", 3000)

    def _discover_plugins(self):
        plugins = self.plugin_manager.discover_plugins()
        if plugins:
            lines = [f"Found {len(plugins)} plugin(s):"]
            for p in plugins:
                lines.append(f"  • {p.get('name', 'Unknown')} v{p.get('version', '?')} - {p.get('description', '')}")
            self.plugin_output.setPlainText('\n'.join(lines))
            self.plugin_status_label.setText(f"{len(plugins)} plugin(s) loaded")
        else:
            self.plugin_output.setPlainText("No plugins found. Add .py files to the 'plugins/' directory.\n\nEach plugin must have a 'register_plugin()' function returning:\n  {'name': str, 'version': str, 'description': str, 'callback': callable}")
            self.plugin_status_label.setText("No plugins discovered")

    def _export_csv(self):
        if not self.mutation_results:
            QMessageBox.warning(self, "No Data", "Run mutation detection first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if not path:
            return
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'position', 'ref_base', 'alt_base', 'ref_codon', 'alt_codon',
                'ref_aa', 'alt_aa', 'classification'
            ])
            writer.writeheader()
            for r in self.mutation_results:
                writer.writerow(r.to_dict())
        self.statusBar().showMessage(f"Exported {len(self.mutation_results)} mutations to {path}", 3000)

    def _export_json(self):
        data = {
            'sequences': [{'name': r.name, 'description': r.description, 'length': len(r), 'md5': r.id}
                          for r in self.sequences],
            'mutations': [r.to_dict() for r in self.mutation_results]
        }
        path, _ = QFileDialog.getSaveFileName(self, "Save JSON", "", "JSON Files (*.json)")
        if not path:
            return
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        self.statusBar().showMessage(f"Exported to {path}", 3000)

    def _export_fasta(self):
        if not self.sequences:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save FASTA", "", "FASTA Files (*.fasta)")
        if not path:
            return
        with open(path, 'w') as f:
            for r in self.sequences:
                f.write(f">{r.name} {r.description}\n")
                for i in range(0, len(r.sequence), 80):
                    f.write(r.sequence[i:i+80] + "\n")
        self.statusBar().showMessage(f"Exported {len(self.sequences)} sequences to {path}", 3000)

    def _export_report(self):
        if not self.sequences:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Report", "", "Text Files (*.txt)")
        if not path:
            return

        lines = []
        lines.append("=" * 70)
        lines.append("ZORA ANALYSIS REPORT")
        lines.append(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("=" * 70)
        lines.append("")

        for rec in self.sequences:
            stats = rec.stats()
            lines.append(f"Sequence: {rec.name}")
            lines.append(f"  Length: {stats['length']} bp")
            lines.append(f"  GC: {stats['gc_content']:.2f}%")
            lines.append(f"  MW: {stats['molecular_weight']:.2f} Da")
            lines.append(f"  Tm: {stats['melting_temp']:.2f} °C")
            lines.append("")

        if self.mutation_results:
            lines.append("Mutations:")
            lines.append(f"  Total: {len(self.mutation_results)}")
            types = Counter(r.classification for r in self.mutation_results)
            for t, c in types.most_common():
                lines.append(f"  {t}: {c}")
            lines.append("")

        with open(path, 'w') as f:
            f.write('\n'.join(lines))
        self.statusBar().showMessage(f"Report saved to {path}", 3000)

    def _export_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save SVG", "", "SVG Files (*.svg)")
        if not path:
            return
        self.graph_canvas.fig.savefig(path, format='svg', dpi=300, bbox_inches='tight')
        self.statusBar().showMessage(f"Graph saved as SVG: {path}", 3000)

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save PNG", "", "PNG Files (*.png)")
        if not path:
            return
        self.graph_canvas.fig.savefig(path, format='png', dpi=300, bbox_inches='tight')
        self.statusBar().showMessage(f"Graph saved as PNG: {path}", 3000)

    def _show_about(self):
        QMessageBox.about(self, "About ZORA",
                          "<h2>ZORA</h2>"
                          "<p><i>Where DNA meets Python, biology meets algorithms, "
                          "and researchers meet fewer headaches.</i></p>"
                          "<hr>"
                          "<p>Built for scientists, students, and curious minds, "
                          "ZORA provides tools for sequence analysis, mutation detection, "
                          "CRISPR guide design, molecular docking, visualization, "
                          "and much more—all in one place.</p>"
                          "<p>Developed by <b>Veera Rahul</b>, a biotechnology graduate "
                          "who apparently thought writing thousands of lines of Python "
                          "was easier than opening five different bioinformatics tools.</p>"
                          "<p><b>Warning:</b> ZORA may cause excessive curiosity, "
                          "spontaneous genome analysis, and an irresistible urge "
                          "to automate everything.</p>"
                          "<hr>"
                          "<p><b>Version:</b> 1.0.0</p>"
                          "<p><b>Powered by:</b> Python, curiosity, "
                          "and just enough caffeine.</p>"
                          )

    def _switch_theme(self, theme_name: str):
        app = QApplication.instance()
        if app:
            app.setStyleSheet(ThemeManager.generate_stylesheet(theme_name))
        self.current_theme = theme_name
        theme = ThemeManager.THEMES.get(theme_name, {})
        mpl_style = theme.get('MATPLOTLIB', 'dark_background')
        try:
            plt.style.use(mpl_style)
        except Exception:
            pass
        self._apply_mpl_rcparams(theme, mpl_style)
        if hasattr(self, 'graph_canvas') and self.graph_canvas is not None:
            self.graph_canvas.draw()
        self.statusBar().showMessage(f"Theme switched to {theme_name}", 3000)
        for name, act in self.theme_actions.items():
            act.setChecked(name == theme_name)
        self._apply_contrast_variant(self.current_contrast)

    def _apply_contrast_variant(self, variant_name: str):
        self.current_contrast = variant_name
        app = QApplication.instance()
        if not app:
            return
        overrides = ThemeManager.CONTRAST_VARIANTS.get(variant_name, {})
        theme = dict(ThemeManager.THEMES.get(self.current_theme, {}))
        if overrides:
            theme.update(overrides)
            app.setStyleSheet(ThemeManager._generate_from_dict(theme))
        else:
            app.setStyleSheet(ThemeManager.generate_stylesheet(self.current_theme))
        mpl_style = theme.get('MATPLOTLIB', 'dark_background')
        try:
            plt.style.use(mpl_style)
        except Exception:
            pass
        self._apply_mpl_rcparams(theme, mpl_style)
        if hasattr(self, 'graph_canvas') and self.graph_canvas is not None:
            self.graph_canvas.draw()
        self.statusBar().showMessage(f"Contrast: {variant_name}", 2000)
        for name, act in self.contrast_actions.items():
            act.setChecked(name == variant_name)

    def _project_context_menu(self, pos):
        item = self.project_tree.itemAt(pos)
        if not item:
            return
        item_type = item.data(0, Qt.UserRole + 1)
        item_key = item.data(0, Qt.UserRole)
        menu = QMenu()

        if item_type == 'seq':
            show_action = menu.addAction("Show Sequence")
            menu.addSeparator()
            send_menu = menu.addMenu("Send to")
            send_pymol = send_menu.addAction("PyMOL")
            send_chimerax = send_menu.addAction("ChimeraX")
            send_snapgene = send_menu.addAction("SnapGene")
            menu.addSeparator()
            export_menu = menu.addMenu("Export")
            export_fasta = export_menu.addAction("FASTA")
            export_csv = export_menu.addAction("CSV")
            menu.addSeparator()
            remove_action = menu.addAction("Remove")
        elif item_type == 'pdb':
            open_pymol = menu.addAction("Open in PyMOL")
            open_chimerax = menu.addAction("Open in ChimeraX")
            menu.addSeparator()
            load_seq = menu.addAction("Load as Sequence")
            menu.addSeparator()
            delete_pdb = menu.addAction("Delete")
        else:
            return

        action = menu.exec(self.project_tree.viewport().mapToGlobal(pos))

        if item_type == 'seq':
            if action == show_action:
                for r in self.sequences:
                    if r.id == item_key:
                        self._show_sequence(r)
                        break
            elif action in (send_pymol, send_chimerax, send_snapgene):
                rec = next((r for r in self.sequences if r.id == item_key), None)
                if rec:
                    seq = self._clean_seq_for_tool(rec)
                    if seq is None:
                        return
                    import tempfile
                    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
                    tmp.write(f">{rec.name}\n{seq}\n")
                    tmp.close()
                    if action == send_pymol:
                        ExternalToolLauncher.launch_pymol(filename=tmp.name)
                    elif action == send_chimerax:
                        ExternalToolLauncher.launch_chimerax(filename=tmp.name)
                    elif action == send_snapgene:
                        ExternalToolLauncher.launch_snapgene(filename=tmp.name)
                    self.statusBar().showMessage(f"Sent {rec.name} to {action.text()}", 3000)
            elif action in (export_fasta, export_csv):
                rec = next((r for r in self.sequences if r.id == item_key), None)
                if rec:
                    path, _ = QFileDialog.getSaveFileName(self, f"Export {rec.name}", "", "FASTA Files (*.fasta)" if action == export_fasta else "CSV Files (*.csv)")
                    if path:
                        if action == export_fasta:
                            with open(path, 'w') as f:
                                f.write(f">{rec.name} {rec.description}\n")
                                for i in range(0, len(rec.sequence), 80):
                                    f.write(rec.sequence[i:i+80] + "\n")
                        else:
                            import csv
                            with open(path, 'w', newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow(['position', 'base'])
                                for i, b in enumerate(rec.sequence):
                                    writer.writerow([i+1, b])
                        self.statusBar().showMessage(f"Exported {rec.name}", 3000)
            elif action == remove_action:
                self.sequences = [r for r in self.sequences if r.id != item_key]
                self._update_project_tree()
                self._update_combos()
                self.statusBar().showMessage("Sequence removed", 3000)
        elif item_type == 'pdb':
            if action == open_pymol:
                ExternalToolLauncher.launch_pymol(filename=item_key)
                self.statusBar().showMessage(f"Opened {os.path.basename(item_key)} in PyMOL", 3000)
            elif action == open_chimerax:
                ExternalToolLauncher.launch_chimerax(filename=item_key)
                self.statusBar().showMessage(f"Opened {os.path.basename(item_key)} in ChimeraX", 3000)
            elif action == load_seq:
                records = PDBParser.parse(item_key)
                if records:
                    self._add_sequences(records)
                else:
                    QMessageBox.warning(self, "No Data", "Could not extract sequences from PDB.")
            elif action == delete_pdb:
                fname = os.path.basename(item_key)
                reply = QMessageBox.warning(self, "Confirm Delete", f"Delete {fname}?", QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.Yes:
                    self.project_manager.delete_pdb(fname)
                    self._update_project_tree()
                    self._refresh_pdb_list()
                    self.statusBar().showMessage(f"Deleted {fname}", 3000)


# ============================================================
# Alignment Worker (Threaded)
# ============================================================
class AlignWorker(QRunnable):
    def __init__(self, seq1: str, seq2: str, method: str, name1: str, name2: str):
        super().__init__()
        self.seq1 = seq1
        self.seq2 = seq2
        self.method = method
        self.name1 = name1
        self.name2 = name2
        self.signals = _AlignSignals()

    def run(self):
        try:
            if self.method == 'nw':
                aln1, aln2, score = AlignmentEngine.needleman_wunsch(self.seq1, self.seq2)
                method_name = "Needleman-Wunsch (Global)"
            else:
                aln1, aln2, score = AlignmentEngine.smith_waterman(self.seq1, self.seq2)
                method_name = "Smith-Waterman (Local)"
            self.signals.finished.emit(aln1, aln2, score, method_name, self.name1, self.name2)
        except Exception as e:
            self.signals.error.emit(str(e))


class _AlignSignals(QObject):
    finished = Signal(str, str, int, str, str, str)
    error = Signal(str)


# ============================================================
# Docking Worker (Threaded)
# ============================================================
class DockingWorker(QRunnable):
    def __init__(self, receptor_pdbqt, ligand_pdbqt, output_dir, box_params):
        super().__init__()
        self.receptor_pdbqt = receptor_pdbqt
        self.ligand_pdbqt = ligand_pdbqt
        self.output_dir = output_dir
        self.box_params = box_params
        self.signals = _DockingSignals()

    def run(self):
        try:
            result = DockingEngine.run_vina(
                self.receptor_pdbqt, self.ligand_pdbqt, self.output_dir,
                **self.box_params
            )
            modes = DockingEngine.parse_vina_log(result['log_path'])
            self.signals.finished.emit(result, modes)
        except Exception as e:
            self.signals.error.emit(str(e))


class _DockingSignals(QObject):
    finished = Signal(dict, list)
    error = Signal(str)


# ============================================================
# Fetch Threads
# ============================================================
class PdbFetchThread(QThread):
    pdb_fetched = Signal(str, str, str)  # accession, data, error

    def __init__(self, accession: str):
        super().__init__()
        self.accession = accession

    def run(self):
        try:
            import urllib.request
            import ssl
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            url = f"https://files.rcsb.org/download/{self.accession}.pdb"
            req = urllib.request.Request(url, headers={'User-Agent': 'ZORA/1.0'})
            resp = urllib.request.urlopen(req, timeout=30, context=ctx)
            data = resp.read().decode('utf-8')
            self.pdb_fetched.emit(self.accession, data, "")
        except Exception as e:
            self.pdb_fetched.emit(self.accession, "", str(e))


class NcbiFetchThread(QThread):
    ncbi_fetched = Signal(str, str, str)  # accession, rec_json, error

    def __init__(self, accession: str):
        super().__init__()
        self.accession = accession

    def run(self):
        try:
            rec = NCBIFetcher.fetch_fasta(self.accession)
            if rec:
                import json
                self.ncbi_fetched.emit(self.accession, json.dumps({'name': rec.name, 'sequence': rec.sequence, 'description': rec.description}), "")
            else:
                self.ncbi_fetched.emit(self.accession, "", "No record returned")
        except Exception as e:
            self.ncbi_fetched.emit(self.accession, "", str(e))


# ============================================================
# PDB Parser
# ============================================================
class PDBParser:
    @staticmethod
    def parse(path: str) -> List[SequenceRecord]:
        with open(path, 'r') as f:
            content = f.read()
        name_base = os.path.splitext(os.path.basename(path))[0]
        return PDBParser._parse_content(name_base, content, path)

    @staticmethod
    def parse_from_text(accession: str, content: str) -> List[SequenceRecord]:
        return PDBParser._parse_content(accession, content, f"PDB:{accession}")

    @staticmethod
    def _parse_content(name_base: str, content: str, source: str) -> List[SequenceRecord]:
        records = []
        chains = {}
        for line in content.split('\n'):
            if line.startswith('SEQRES'):
                if len(line) > 18:
                    chain = line[11]
                    seq = ''.join(line[19:].split())
                    if chain not in chains:
                        chains[chain] = ''
                    chains[chain] += seq
        for chain, seq in chains.items():
            name = f"{name_base}_chain{chain}"
            records.append(SequenceRecord(name, seq, f"From PDB: {source}"))
        if not chains:
            atoms = {}
            for line in content.split('\n'):
                if line.startswith('ATOM') and len(line) > 22:
                    chain = line[21]
                    resname = line[17:20].strip()
                    if chain not in atoms:
                        atoms[chain] = []
                    if resname not in atoms[chain]:
                        atoms[chain].append(resname)
            three_to_one = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            }
            for chain, residues in atoms.items():
                seq = ''.join(three_to_one.get(r, 'X') for r in residues)
                if seq:
                    name = f"{name_base}_chain{chain}"
                    records.append(SequenceRecord(name, seq, f"From PDB ATOM: {source}"))
        return records


# ============================================================
# NCBI / Database Fetcher
# ============================================================
class NCBIFetcher:
    @staticmethod
    def fetch_genbank(accession: str) -> Optional[SequenceRecord]:
        import urllib.request
        import urllib.parse
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=gb&retmode=text"
            req = urllib.request.Request(url, headers={'User-Agent': 'ZORA/1.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read().decode('utf-8')
            records = FileParser.parse_genbank_from_text(data)
            return records[0] if records else None
        except Exception as e:
            print(f"NCBI fetch error: {e}")
            return None

    @staticmethod
    def fetch_fasta(accession: str) -> Optional[SequenceRecord]:
        import urllib.request
        import urllib.parse
        try:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
            req = urllib.request.Request(url, headers={'User-Agent': 'ZORA/1.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read().decode('utf-8')
            records = FileParser.parse_fasta_from_text(data)
            return records[0] if records else None
        except Exception as e:
            print(f"NCBI fetch error: {e}")
            return None

    @staticmethod
    def search_ncbi(query: str, db: str = 'nucleotide', max_results: int = 10) -> List[Dict]:
        import urllib.request
        import urllib.parse
        import xml.etree.ElementTree as ET
        try:
            params = urllib.parse.urlencode({'db': db, 'term': query, 'retmax': max_results, 'retmode': 'json'})
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
            req = urllib.request.Request(url, headers={'User-Agent': 'ZORA/1.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
            ids = data.get('esearchresult', {}).get('idlist', [])
            results = []
            for rid in ids:
                results.append({'id': rid, 'accession': rid})
            return results
        except Exception as e:
            print(f"NCBI search error: {e}")
            return []


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
        exe = ExternalToolLauncher.find_executable('chimerax')
        if not exe:
            exe = ExternalToolLauncher.find_executable('ChimeraX')
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
        args = [exe]
        if filename:
            args.append(filename)
        import subprocess
        subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
        """Export as FASTA for ChimeraX sequence loading"""
        ExternalToolLauncher.export_for_snapgene(records, path)

    @staticmethod
    def send_to_pymol_via_file(records: List[SequenceRecord], path: str):
        """Write a PyMOL script that loads sequence data"""
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
# Docking Engine
# ============================================================
class DockingEngine:
    VINA_PATH = shutil.which('vina') or os.environ.get('ZORA_VINA_PATH', '/usr/local/bin/vina')
    _obabel_which = shutil.which('obabel')
    OBABEL_PATH = _obabel_which or os.environ.get('ZORA_OBABEL_PATH', None)
    _OBABEL_ENV = None

    @staticmethod
    def _obabel_env():
        if DockingEngine._OBABEL_ENV is None:
            env = os.environ.copy()
            mgltools_root = os.environ.get('ZORA_MGLTOOLS_ROOT', '')
            lib_dir = os.path.join(mgltools_root, 'lib') if mgltools_root else ''
            existing = env.get('LD_LIBRARY_PATH', '')
            if lib_dir not in existing:
                env['LD_LIBRARY_PATH'] = f'{lib_dir}:{existing}' if existing else lib_dir
            DockingEngine._OBABEL_ENV = env
        return DockingEngine._OBABEL_ENV

    @staticmethod
    def check_available() -> Dict[str, bool]:
        return {
            'vina': os.path.exists(DockingEngine.VINA_PATH) if DockingEngine.VINA_PATH else False,
            'obabel': DockingEngine.OBABEL_PATH is not None,
        }

    @staticmethod
    def _sanitize_receptor_pdbqt(pdbqt_path: str, output_dir: str) -> str:
        """Convert any PDB/PDBQT to a clean rigid PDBQT with valid AD4 atom types."""
        out_path = os.path.join(output_dir, 'receptor_rigid.pdbqt')
        _ELEM_TO_AD = {
            'H': 'H', 'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P',
            'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I', 'SI': 'Si',
            'FE': 'Fe', 'CA': 'Ca', 'MG': 'Mg', 'ZN': 'Zn', 'MN': 'Mn',
            'CU': 'Cu', 'NA': 'Na', 'K': 'K', 'SE': 'Se',
        }
        _VALID_AD = set(_ELEM_TO_AD.values())
        try:
            has_root = False
            is_valid_pdbqt = False
            with open(pdbqt_path) as f:
                for line in f:
                    if line.startswith('ROOT'):
                        has_root = True
                    if line.startswith(('ATOM', 'HETATM')) and len(line) >= 80:
                        charge_col = line[66:76]
                        if charge_col.strip() == '':
                            continue
                        ad = line[77:79].rstrip()
                        if ad in _VALID_AD:
                            is_valid_pdbqt = True
            if not has_root and is_valid_pdbqt:
                return pdbqt_path

            with open(pdbqt_path) as f_in, open(out_path, 'w') as f_out:
                for line in f_in:
                    if line.startswith(('ATOM', 'HETATM')):
                        if is_valid_pdbqt:
                            f_out.write(line)
                        else:
                            if len(line) < 54:
                                continue
                            try:
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                            except ValueError:
                                continue
                            elem = line[76:78].strip().upper() if len(line) >= 78 else ''
                            if not elem or elem == '**':
                                raw = line[12:16].strip().upper()
                                elem = raw.rstrip('0123456789+-')[:2]
                            ad = _ELEM_TO_AD.get(elem, 'C')
                            prefix = line[:66].rstrip()
                            f_out.write(f"{prefix:<66}{'0.000':>10} {ad:<2}\n")
                    elif line.startswith(('REMARK', 'END', 'TER')):
                        f_out.write(line)
            return out_path if os.path.exists(out_path) else pdbqt_path
        except Exception:
            return pdbqt_path

    @staticmethod
    def prepare_receptor(pdb_path: str, output_dir: str) -> Optional[str]:
        pdbqt_path = os.path.join(output_dir, 'receptor.pdbqt')
        # 1) Try obabel
        if DockingEngine.OBABEL_PATH:
            cmd = [DockingEngine.OBABEL_PATH, pdb_path, '-o', 'pdbqt', '-O', pdbqt_path, '-xr']
            subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=DockingEngine._obabel_env())
            if os.path.exists(pdbqt_path):
                return pdbqt_path
        # 2) Try MGLTools prepare_receptor4.py via pythonsh
        obabel_dir = os.path.dirname(DockingEngine.OBABEL_PATH) if DockingEngine.OBABEL_PATH else ''
        pythonsh = os.path.join(obabel_dir, 'pythonsh')
        prep4 = None
        if os.path.exists(pythonsh):
            mgl_root = os.path.dirname(obabel_dir)
            prep4 = os.path.join(mgl_root, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        if prep4 and os.path.exists(prep4):
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    shutil.copy2(pdb_path, tmpdir)
                    base = os.path.basename(pdb_path)
                    out_base = os.path.basename(pdbqt_path)
                    subprocess.run(
                        [pythonsh, prep4, '-r', os.path.join(tmpdir, base), '-o', os.path.join(tmpdir, out_base)],
                        capture_output=True, text=True, timeout=120)
                    tmp_out = os.path.join(tmpdir, out_base)
                    if os.path.exists(tmp_out):
                        shutil.copy2(tmp_out, pdbqt_path)
                        return pdbqt_path
            except Exception:
                pass
        # 3) Fallback: write minimal PDBQT (will be sanitized by _sanitize_receptor_pdbqt)
        try:
            with open(pdb_path) as f_in, open(pdbqt_path, 'w') as f_out:
                for line in f_in:
                    if line.startswith(('ATOM', 'HETATM')) and len(line) > 54:
                        f_out.write(line.rstrip() + '\n')
            return pdbqt_path if os.path.exists(pdbqt_path) else None
        except Exception:
            return None

    @staticmethod
    def prepare_ligand(input_path: str, output_dir: str) -> Optional[str]:
        if not DockingEngine.OBABEL_PATH:
            return None
        pdbqt_path = os.path.join(output_dir, 'ligand.pdbqt')
        cmd = [DockingEngine.OBABEL_PATH, input_path, '-o', 'pdbqt', '-O', pdbqt_path]
        subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=DockingEngine._obabel_env())
        if not os.path.exists(pdbqt_path):
            cmd.append('--gen3d')
            subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=DockingEngine._obabel_env())
        if os.path.exists(pdbqt_path):
            return pdbqt_path
        # Fallback: meeko via RDKit
        try:
            from rdkit import Chem
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            mol = Chem.MolFromMolFile(input_path) or Chem.MolFromPDBFile(input_path)
            if mol is None:
                return None
            mol = Chem.AddHs(mol)
            prep = MoleculePreparation()
            setups = prep.prepare(mol)
            if not setups:
                return None
            pdbqt_string, is_ok, log = PDBQTWriterLegacy.write_string(setups[0])
            if not is_ok:
                return None
            with open(pdbqt_path, 'w') as f:
                f.write(pdbqt_string)
            return pdbqt_path if os.path.exists(pdbqt_path) else None
        except Exception:
            return None

    @staticmethod
    def ligand_from_smiles(smiles: str, output_dir: str, name: str = 'ligand') -> Optional[str]:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mol = Chem.AddHs(mol)
            ret = AllChem.EmbedMolecule(mol, randomSeed=42)
            if ret != 0:
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
            AllChem.MMFFOptimizeMolecule(mol)
            prep = MoleculePreparation()
            setups = prep.prepare(mol)
            if not setups:
                return None
            pdbqt_string, is_ok, log = PDBQTWriterLegacy.write_string(setups[0])
            if not is_ok:
                return None
            pdbqt_path = os.path.join(output_dir, f'{name}.pdbqt')
            with open(pdbqt_path, 'w') as f:
                f.write(pdbqt_string)
            return pdbqt_path if os.path.exists(pdbqt_path) else None
        except Exception:
            return None

    @staticmethod
    def run_vina(receptor_pdbqt: str, ligand_pdbqt: str, output_dir: str,
                 center_x: float = 0, center_y: float = 0, center_z: float = 0,
                 size_x: float = 30, size_y: float = 30, size_z: float = 30,
                 exhaustiveness: int = 8, num_modes: int = 9) -> dict:
        out_path = os.path.join(output_dir, 'docking_results.pdbqt')
        log_path = os.path.join(output_dir, 'docking_log.txt')
        cmd = [
            DockingEngine.VINA_PATH,
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', out_path,
            '--center_x', str(center_x),
            '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(size_x),
            '--size_y', str(size_y),
            '--size_z', str(size_z),
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        with open(log_path, 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write('\nSTDERR:\n' + result.stderr)
        return {
            'out_path': out_path if os.path.exists(out_path) else None,
            'log_path': log_path,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'returncode': result.returncode,
            'cmd': ' '.join(cmd),
        }

    @staticmethod
    def parse_vina_log(log_path: str) -> List[dict]:
        import re
        results = []
        try:
            with open(log_path) as f:
                text = f.read()
        except Exception:
            return results
        mode_start = False
        for line in text.split('\n'):
            if 'mode |   affinity | dist from best' in line:
                mode_start = True
                continue
            if mode_start and line.strip() and not line.startswith('---'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        results.append({
                            'mode': int(parts[0]),
                            'affinity': float(parts[1]),
                            'rmsd_lb': float(parts[2]),
                            'rmsd_ub': float(parts[3]),
                        })
                    except ValueError:
                        pass
        return results

    @staticmethod
    def get_pdb_center(pdb_path: str) -> Optional[Dict[str, float]]:
        try:
            xs, ys, zs = [], [], []
            with open(pdb_path) as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')) and len(line) > 54:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
            if xs:
                return {
                    'center_x': (max(xs) + min(xs)) / 2,
                    'center_y': (max(ys) + min(ys)) / 2,
                    'center_z': (max(zs) + min(zs)) / 2,
                    'size_x': max(xs) - min(xs) + 10,
                    'size_y': max(ys) - min(ys) + 10,
                    'size_z': max(zs) - min(zs) + 10,
                }
        except Exception:
            pass
        return None


# ============================================================
# Plugin Manager
# ============================================================
class PluginManager:
    def __init__(self):
        self.plugins = []
        default = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plugins')
        self.plugin_dir = os.environ.get('ZORA_PLUGIN_DIR', default)

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
# Project Manager
# ============================================================
class ProjectManager:
    def __init__(self):
        self.project_dir: Optional[str] = None
        self.metadata: Dict = {
            'name': 'Untitled',
            'version': '1.0',
            'created': '',
            'modified': '',
            'sequences': [],
            'mutations': [],
        }

    def is_open(self) -> bool:
        return self.project_dir is not None and os.path.isdir(self.project_dir)

    def create(self, path: str, name: str = "Untitled") -> bool:
        import datetime
        self.project_dir = path
        os.makedirs(path, exist_ok=True)
        os.makedirs(os.path.join(path, 'sequences'), exist_ok=True)
        os.makedirs(os.path.join(path, 'pdb'), exist_ok=True)
        os.makedirs(os.path.join(path, 'analysis'), exist_ok=True)
        self.metadata = {
            'name': name,
            'version': '1.0',
            'created': datetime.datetime.now().isoformat(),
            'modified': datetime.datetime.now().isoformat(),
            'sequences': [],
            'mutations': [],
        }
        return self._save_metadata()

    def open(self, path: str) -> bool:
        if not os.path.isdir(path):
            return False
        self.project_dir = path
        os.makedirs(os.path.join(path, 'sequences'), exist_ok=True)
        os.makedirs(os.path.join(path, 'pdb'), exist_ok=True)
        os.makedirs(os.path.join(path, 'analysis'), exist_ok=True)
        meta_path = os.path.join(path, 'project.json')
        if os.path.exists(meta_path):
            with open(meta_path) as f:
                self.metadata = json.load(f)
        return True

    def close(self):
        self.project_dir = None
        self.metadata = {'name': 'Untitled', 'version': '1.0', 'created': '', 'modified': '', 'sequences': [], 'mutations': []}

    def save(self, sequences=None, mutation_results=None, alignment_results=None, analysis_text=""):
        if not self.is_open():
            return False
        import datetime
        self.metadata['modified'] = datetime.datetime.now().isoformat()
        if sequences is not None:
            self.metadata['sequence_count'] = len(sequences)
            seq_dir = os.path.join(self.project_dir, 'sequences')
            seq_data = []
            for rec in sequences:
                fname = f"{rec.id}.fasta"
                fpath = os.path.join(seq_dir, fname)
                with open(fpath, 'w') as f:
                    f.write(f">{rec.name} {rec.description}\n")
                    for i in range(0, len(rec.sequence), 80):
                        f.write(rec.sequence[i:i+80] + "\n")
                seq_data.append({
                    'id': rec.id, 'name': rec.name, 'description': rec.description,
                    'length': len(rec), 'file': fname,
                })
            self.metadata['sequences'] = seq_data
        if mutation_results is not None:
            self.metadata['mutations'] = [r.to_dict() for r in mutation_results]
            analysis_dir = os.path.join(self.project_dir, 'analysis')
            with open(os.path.join(analysis_dir, 'mutations.json'), 'w') as f:
                json.dump(self.metadata['mutations'], f, indent=2)
        if alignment_results is not None:
            aln_data = [{'seq1': a[0], 'seq2': a[1], 'score': a[2]} for a in alignment_results]
            analysis_dir = os.path.join(self.project_dir, 'analysis')
            with open(os.path.join(analysis_dir, 'alignments.json'), 'w') as f:
                json.dump(aln_data, f, indent=2)
        if analysis_text:
            analysis_dir = os.path.join(self.project_dir, 'analysis')
            with open(os.path.join(analysis_dir, 'analysis.txt'), 'w') as f:
                f.write(analysis_text)
        return self._save_metadata()

    def load_sequences(self) -> list:
        if not self.is_open():
            return []
        records = []
        for sd in self.metadata.get('sequences', []):
            fpath = os.path.join(self.project_dir, 'sequences', sd['file'])
            if os.path.exists(fpath):
                import Bio.SeqIO  # fallback: use FileParser
                seq = open(fpath).read()
                recs = FileParser.parse_fasta_from_text(seq)
                if recs:
                    records.extend(recs)
        return records

    def load_mutations(self) -> list:
        if not self.is_open():
            return []
        analysis_dir = os.path.join(self.project_dir, 'analysis')
        mpath = os.path.join(analysis_dir, 'mutations.json')
        if os.path.exists(mpath):
            with open(mpath) as f:
                return json.load(f)
        return []

    def get_pdb_files(self) -> list:
        if not self.is_open():
            return []
        pdb_dir = os.path.join(self.project_dir, 'pdb')
        if not os.path.isdir(pdb_dir):
            return []
        files = []
        for fname in sorted(os.listdir(pdb_dir)):
            if fname.lower().endswith('.pdb'):
                fpath = os.path.join(pdb_dir, fname)
                files.append({
                    'name': fname,
                    'path': fpath,
                    'size': os.path.getsize(fpath),
                    'modified': datetime.datetime.fromtimestamp(os.path.getmtime(fpath)).isoformat(),
                })
        return files

    def save_pdb(self, pdb_id: str, content: str) -> Optional[str]:
        if not self.is_open():
            return None
        pdb_dir = os.path.join(self.project_dir, 'pdb')
        os.makedirs(pdb_dir, exist_ok=True)
        fname = f"{pdb_id}.pdb"
        fpath = os.path.join(pdb_dir, fname)
        with open(fpath, 'w') as f:
            f.write(content)
        return fpath

    def delete_pdb(self, filename: str) -> bool:
        if not self.is_open():
            return False
        fpath = os.path.join(self.project_dir, 'pdb', filename)
        if os.path.exists(fpath):
            os.remove(fpath)
            return True
        return False

    def get_path(self, *parts: str) -> str:
        return os.path.join(self.project_dir, *parts)

    def save_docking_pdbqt(self, name: str, content: bytes) -> Optional[str]:
        if not self.is_open():
            return None
        d_dir = self.get_path('docking')
        os.makedirs(d_dir, exist_ok=True)
        fpath = os.path.join(d_dir, f"{name}.pdbqt")
        with open(fpath, 'wb') as f:
            f.write(content)
        return fpath

    def get_docking_files(self) -> list:
        if not self.is_open():
            return []
        d_dir = os.path.join(self.project_dir, 'docking')
        if not os.path.isdir(d_dir):
            return []
        files = []
        for fname in sorted(os.listdir(d_dir)):
            if fname.lower().endswith('.pdbqt'):
                fpath = os.path.join(d_dir, fname)
                files.append({
                    'name': fname,
                    'path': fpath,
                    'size': os.path.getsize(fpath),
                    'modified': datetime.datetime.fromtimestamp(os.path.getmtime(fpath)).isoformat(),
                })
        return files

    def split_docking_poses(self, pdbqt_path: str) -> List[str]:
        """Split multi-model PDBQT into individual pose PDB files."""
        poses = []
        try:
            d_dir = os.path.dirname(pdbqt_path)
            with open(pdbqt_path) as f:
                content = f.read()
            models = content.split('ENDMDL\n')
            for i, model in enumerate(models):
                model = model.strip()
                if not model:
                    continue
                model_name = f"pose_{i+1}"
                # Build PDB from PDBQT lines (strip charge/AD type columns)
                pdb_lines = []
                for line in model.split('\n'):
                    if line.startswith(('ATOM', 'HETATM')):
                        pdb_line = line[:66].rstrip() + '\n'
                        if len(pdb_line) > 54:
                            pdb_lines.append(pdb_line)
                    elif line == 'END':
                        pdb_lines.append('END\n')
                if pdb_lines:
                    pdb_path = os.path.join(d_dir, f"{model_name}.pdb")
                    with open(pdb_path, 'w') as f:
                        f.writelines(pdb_lines)
                    poses.append(pdb_path)
        except Exception:
            pass
        return poses

    def _save_metadata(self) -> bool:
        if not self.is_open():
            return False
        with open(os.path.join(self.project_dir, 'project.json'), 'w') as f:
            json.dump(self.metadata, f, indent=2)
        return True


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


# ============================================================
# Application Entry Point
# ============================================================
def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')

    app.setStyleSheet(ThemeManager.generate_stylesheet('Deep Scientific Navy'))

    window = ZORAMainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
