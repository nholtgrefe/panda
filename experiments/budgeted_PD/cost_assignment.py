"""
Assign random integer taxon costs (lognormal) for networks listed in ``nonbinary_nets.csv``.

Networks are parsed with PhyloZoo's extended-Newick reader
(:meth:`phylozoo.DirectedPhyNetwork.from_string`; there is no separate
``phylozoo.load_from_string`` symbol—this class method is the supported I/O entry
point).  Taxon names come from :attr:`DirectedPhyNetwork.taxa`.

Writes ``costs.csv`` in the same directory: one line per network, space-separated
``net_id`` and a concatenated list of ``(taxon,cost)`` pairs.

This script was run with phylozoo version 0.1.2.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from phylozoo import DirectedPhyNetwork

# Lognormal parameters: underlying normal has mean ``mu`` and std ``sigma``.
_MU = math.log(5.0)
_SIGMA = 0.5


def sample_taxon_costs(taxa: list[str], *, rng: np.random.Generator | None = None) -> list[tuple[str, int]]:
    """
    Sample an integer cost per taxon from a lognormal distribution.

    Draws are i.i.d. ``LogNormal(mu=ln(5), sigma=0.5)``, then rounded to the
    nearest integer (zero is allowed).

    Parameters
    ----------
    taxa : list[str]
        Taxon names in the order costs should be paired with.
    rng : numpy.random.Generator | None, optional
        Random generator. If ``None``, uses ``numpy.random.default_rng()``.

    Returns
    -------
    list[tuple[str, int]]
        ``(taxon, cost)`` pairs aligned with ``taxa``.
    """
    if rng is None:
        rng = np.random.default_rng()
    n = len(taxa)
    if n == 0:
        return []
    raw = rng.lognormal(mean=_MU, sigma=_SIGMA, size=n)
    costs = np.rint(raw).astype(int)
    return list(zip(taxa, costs.tolist()))


def _read_network_rows(csv_path: Path) -> list[tuple[str, str]]:
    """Return ``(net_id, newick)`` for each data row after the header."""
    rows: list[tuple[str, str]] = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=" ")
        header = next(reader, None)
        if header is None:
            return rows
        for parts in reader:
            if not parts or not parts[0].strip():
                continue
            if len(parts) < 4:
                continue
            net_id = parts[0].strip()
            newick = " ".join(parts[3:]).strip()
            rows.append((net_id, newick))
    return rows


def main() -> None:
    """Build ``costs.csv`` for every network in ``nonbinary_nets.csv``."""
    here = Path(__file__).resolve().parent
    nets_path = here / "nonbinary_nets.csv"
    out_path = here / "costs.csv"

    rng = np.random.default_rng()
    lines_out: list[str] = []

    for net_id, newick in _read_network_rows(nets_path):
        network = DirectedPhyNetwork.from_string(newick.strip())
        taxa = sorted(network.taxa)
        pairs = sample_taxon_costs(taxa, rng=rng)
        pair_str = ",".join(f"({tax},{cost})" for tax, cost in pairs)
        lines_out.append(f"{net_id} {pair_str}")

    out_path.write_text("\n".join(lines_out) + ("\n" if lines_out else ""), encoding="utf-8")
    print(f"Wrote {len(lines_out)} rows to {out_path}")


if __name__ == "__main__":
    main()
