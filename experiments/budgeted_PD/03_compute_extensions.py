"""
XP node-scanwidth on PhyloZoo networks and batch export of optimal extensions.

CSV eNewick cells are parsed with :meth:`~phylozoo.utils.io.IOMixin.from_string`
(default format ``enewick``). Output ``extensions_nodewidth.csv`` uses space-separated columns
(like ``networks.csv``); the ``extension`` field is one column whose value is
the vertex order as comma-separated PhyloZoo labels (taxa and internal names).
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Any, Tuple

from phylozoo import DirectedPhyNetwork
from scanwidth import DAG, Extension, node_scanwidth
from scanwidth.node_scanwidth.reduction.config import ReducerConfig


def _node_label_or_id(network: DirectedPhyNetwork, node: Any) -> str:
    """Taxon / internal label from PhyloZoo, falling back to the graph node id."""
    label = network.get_label(node)
    return str(node) if label is None else str(label)


def optimize_nsw(
    network: DirectedPhyNetwork,
) -> Tuple[int, Extension, float, list[str]]:
    """
    Compute optimal node-scanwidth via the XP solver (reduction on, no s-block pools).

    Timing covers the scanwidth call only (after the network object exists).

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input DAG (no parallel edges required by scanwidth's underlying simple DAG view).

    Returns
    -------
    nsw : int
        Optimal node-scanwidth value.
    extension : Extension
        Corresponding scanwidth :class:`~scanwidth.extension.Extension`.
    comp_time_xp : float
        Wall time in seconds for ``node_scanwidth`` (XP + reducer).
    extension_labels : list of str
        ``extension.ordering`` with each vertex replaced by
        :meth:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork.get_label`
        when present, else the node id as string.

    Notes
    -----
    Uses ``algorithm='xp'``, ``reduce=True``, and
    ``ReducerConfig(parallel_sblocks=False)`` so s-blocks are not solved in parallel.
    """
    dag = DAG(network._graph._graph)
    t0 = time.perf_counter()
    nsw, extension = node_scanwidth(
        dag,
        algorithm="xp",
        reduce=True,
        reducer_config=ReducerConfig(parallel_sblocks=False),
    )
    elapsed = time.perf_counter() - t0
    extension_labels = [
        _node_label_or_id(network, v) for v in extension.ordering
    ]
    return nsw, extension, elapsed, extension_labels


def _run_nonbinary_extensions_csv() -> None:
    """Read ``networks.csv`` and write ``extensions_nodewidth.csv`` with XP NSW results."""
    base = Path(__file__).resolve().parent
    src = base / "networks.csv"
    dst = base / "extensions_nodewidth.csv"

    with src.open(encoding="utf-8") as f_in, dst.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        writer = csv.writer(f_out, delimiter=" ")
        writer.writerow(["id", "nsw", "comp_time_xp", "extension"])
        for raw in f_in:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(maxsplit=3)
            if parts[0] == "id" or len(parts) < 4:
                continue
            net_id, _ntips, _level, newick = parts
            try:
                net = DirectedPhyNetwork.from_string(newick.strip())
                nsw, _ext, t_xp, ext_labels = optimize_nsw(net)
            except Exception as exc:
                print(f"{net_id}: skipped ({exc})")
                continue
            writer.writerow(
                [net_id, nsw, f"{t_xp:.9f}", ",".join(ext_labels)]
            )


if __name__ == "__main__":
    _run_nonbinary_extensions_csv()
