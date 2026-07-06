"""Thread workers for ZORA (QRunnable, QThread)."""

from PySide6.QtCore import QThread, QObject, QRunnable, Signal

from zora.models import AlignmentEngine, SequenceRecord
from zora.docking import DockingEngine
from zora.project import NCBIFetcher, FileParser


# ============================================================
# Alignment Signals
# ============================================================
class _AlignSignals(QObject):
    finished = Signal(str, str, int, str, str, str)
    error = Signal(str)


# ============================================================
# Alignment Worker (Threaded via QRunnable)
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


# ============================================================
# Docking Signals
# ============================================================
class _DockingSignals(QObject):
    finished = Signal(dict, list)
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


# ============================================================
# PDB Fetch Thread
# ============================================================
class PdbFetchThread(QThread):
    pdb_fetched = Signal(str, str, str)

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


# ============================================================
# NCBI Fetch Thread
# ============================================================
class NcbiFetchThread(QThread):
    ncbi_fetched = Signal(str, str, str)

    def __init__(self, accession: str):
        super().__init__()
        self.accession = accession

    def run(self):
        try:
            rec = NCBIFetcher.fetch_fasta(self.accession)
            if rec:
                import json
                self.ncbi_fetched.emit(
                    self.accession,
                    json.dumps({'name': rec.name, 'sequence': rec.sequence, 'description': rec.description}),
                    ""
                )
            else:
                self.ncbi_fetched.emit(self.accession, "", "No record returned")
        except Exception as e:
            self.ncbi_fetched.emit(self.accession, "", str(e))
