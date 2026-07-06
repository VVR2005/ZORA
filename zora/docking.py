"""AutoDock Vina docking engine for ZORA."""

import os
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional


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
            existing = env.get('LD_LIBRARY_PATH', '')
            if mgltools_root:
                lib_dir = os.path.join(mgltools_root, 'lib')
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
        if DockingEngine.OBABEL_PATH:
            cmd = [DockingEngine.OBABEL_PATH, pdb_path, '-o', 'pdbqt', '-O', pdbqt_path, '-xr']
            subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=DockingEngine._obabel_env())
            if os.path.exists(pdbqt_path):
                return pdbqt_path
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
