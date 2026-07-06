"""Pure data structures and algorithms for ZORA (no Qt dependencies)."""

import re
import hashlib
from collections import Counter
from typing import List, Tuple, Dict, Optional

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
    '*': 'Stop Codon',
}

COMPLEMENT_TABLE = str.maketrans('ATCGNatcgn', 'TAGCNTAGCN')

MUTATION_COLORS = {
    'Missense (Non-synonymous)': '#FF6B6B',
    'Nonsense': '#FF4444',
    'Silent (Synonymous)': '#51CF66',
    'Stop Loss': '#FF922B',
    'Intergenic': '#868E96',
    'Insertion': '#339AF0',
    'Deletion': '#F06595',
}

MUTATION_DESCRIPTIONS = {
    'Missense': 'Single nucleotide change that results in a different amino acid',
    'Nonsense': 'Single nucleotide change causing premature stop codon',
    'Silent': 'Single nucleotide change with no amino acid change',
    'Frameshift': 'Insertion/deletion that shifts the reading frame',
}


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
        return {base: (seq.count(base) / length) * 100.0 for base in 'ATCG'}

    @staticmethod
    def translate(seq: str, frame: int = 0) -> str:
        seq = seq.upper()[frame:]
        return ''.join(AMINO_ACIDS.get(seq[i:i+3], '?') for i in range(0, len(seq) - 2, 3))

    @staticmethod
    def find_orfs(seq: str, min_length: int = 30) -> List[Tuple[int, int, str, str]]:
        orfs = []
        seq = seq.upper()
        for frame in range(3):
            i = frame
            while i < len(seq) - 2:
                if seq[i:i+3] == 'ATG':
                    start = i
                    for j in range(i, len(seq) - 2, 3):
                        if AMINO_ACIDS.get(seq[j:j+3], 'X') == '*':
                            length = j + 3 - start
                            if length >= min_length:
                                protein = SeqUtils.translate(seq[start:j+3], 0)
                                orfs.append((start, j + 3, protein, f"+{frame+1}"))
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
        if regex:
            return [(m.start(), m.group()) for m in re.finditer(pattern, seq, re.IGNORECASE)]
        pattern = pattern.upper()
        matches = []
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
                    'freq': usage[aa]['codons'][codon] / usage[aa]['total'] * 100,
                }
        return usage

    @staticmethod
    def molecular_weight(seq: str) -> float:
        weights = {'A': 313.21, 'T': 304.20, 'G': 329.21, 'C': 289.18}
        return round(sum(weights.get(b, 0) for b in seq.upper()), 2)

    @staticmethod
    def melting_temp(seq: str) -> float:
        seq = seq.upper()
        if len(seq) < 14:
            return seq.count('A') * 2 + seq.count('T') * 2 + seq.count('G') * 4 + seq.count('C') * 4
        gc = seq.count('G') + seq.count('C')
        return 64.9 + 41.0 * (gc - 16.4) / len(seq)

    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> int:
        return sum(1 for i in range(min(len(seq1), len(seq2))) if seq1[i].upper() != seq2[i].upper())


# ============================================================
# Data Models
# ============================================================
class SequenceRecord:
    def __init__(self, name: str = "", sequence: str = "", description: str = ""):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description
        self.id = hashlib.md5(sequence.encode()).hexdigest()[:8]
        self._stats_cache: Optional[dict] = None

    def __len__(self):
        return len(self.sequence)

    def clear_cache(self):
        self._stats_cache = None

    def stats(self) -> Dict:
        if self._stats_cache is not None:
            return self._stats_cache
        self._stats_cache = {
            'length': len(self.sequence),
            'gc_content': SeqUtils.gc_content(self.sequence),
            'gc_skew': SeqUtils.gc_skew(self.sequence),
            'molecular_weight': SeqUtils.molecular_weight(self.sequence),
            'melting_temp': SeqUtils.melting_temp(self.sequence),
            'nucleotide_freq': SeqUtils.nucleotide_frequency(self.sequence),
        }
        return self._stats_cache


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
            'classification': self.classification,
        }


# ============================================================
# Mutation Analysis Engine
# ============================================================
class MutationEngine:
    @staticmethod
    def detect_mutations(seq1: str, seq2: str) -> List[MutationResult]:
        seq1, seq2 = seq1.upper(), seq2.upper()
        results = []
        i = j = 0
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
                    results.append(MutationEngine._classify_snp(seq1, seq2, i))
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
        ref_base, alt_base = seq1[pos].upper(), seq2[pos].upper()
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
        i, j = n, m
        aln1 = aln2 = ''
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
        aln1 = aln2 = ''
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
# PDB Utilities
# ============================================================
NONSTANDARD_AA_MAP = {
    'U': 'C', 'O': 'K', 'B': 'N', 'Z': 'Q', 'J': 'L', 'X': 'A',
}


def map_nonstandard_aa(seq: str) -> str:
    """Map non-standard amino acids to standard ones for external tools."""
    return ''.join(NONSTANDARD_AA_MAP.get(c, c) for c in seq)
