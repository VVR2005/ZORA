"""Project management, file I/O, and remote fetch for ZORA."""

import os
import re
import json
import csv
import datetime
from typing import Dict, List, Optional

from zora.models import SequenceRecord, SeqUtils


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
        try:
            params = urllib.parse.urlencode({'db': db, 'term': query, 'retmax': max_results, 'retmode': 'json'})
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
            req = urllib.request.Request(url, headers={'User-Agent': 'ZORA/1.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
            ids = data.get('esearchresult', {}).get('idlist', [])
            return [{'id': rid, 'accession': rid} for rid in ids]
        except Exception as e:
            print(f"NCBI search error: {e}")
            return []


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
        self.project_dir = path
        os.makedirs(path, exist_ok=True)
        for sub in ('sequences', 'pdb', 'analysis'):
            os.makedirs(os.path.join(path, sub), exist_ok=True)
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
        for sub in ('sequences', 'pdb', 'analysis'):
            os.makedirs(os.path.join(path, sub), exist_ok=True)
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
