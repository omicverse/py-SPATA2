from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from anndata import AnnData

import spata2py as spata


def make_fixture(n_side: int = 80, n_genes: int = 100, seed: int = 42) -> AnnData:
    rng = np.random.default_rng(seed)
    grid_x, grid_y = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.column_stack([grid_x.ravel(), grid_y.ravel()]).astype(float)
    outliers = np.array([[n_side * 4.0, n_side * 4.0], [n_side * 4.2, n_side * 4.1]])
    coords = np.vstack([coords, outliers])
    counts = rng.poisson(1.2, size=(coords.shape[0], n_genes)).astype(float)
    obs = pd.DataFrame(
        {"region": np.where(coords[:, 0] < n_side / 2, "left", "right")},
        index=[f"spot_{idx}" for idx in range(coords.shape[0])],
    )
    var = pd.DataFrame(index=[f"Gene{idx}" for idx in range(n_genes)])
    adata = AnnData(X=counts, obs=obs, var=var)
    adata.obsm["spatial"] = coords
    return adata


def timed(name: str, fn):
    start = time.perf_counter()
    result = fn()
    elapsed = time.perf_counter() - start
    return name, elapsed, result


def main() -> None:
    adata = make_fixture()
    runs = []

    for name, fn in [
        ("getCoordsDf", lambda: spata.getCoordsDf(adata)),
        ("getFeatureDf", lambda: spata.getFeatureDf(adata, ["region", "Gene0", "Gene1"])),
        ("identifyTissueOutline", lambda: spata.identifyTissueOutline(adata)),
        (
            "identifySpatialOutliers",
            lambda: spata.identifySpatialOutliers(adata, radius=2.0, min_neighbors=2),
        ),
        ("removeSpatialOutliers", lambda: spata.removeSpatialOutliers(adata)),
    ]:
        label, elapsed, result = timed(name, fn)
        runs.append(
            {
                "operation": label,
                "seconds": elapsed,
                "result_shape": list(getattr(result, "shape", (len(result),))),
            }
        )

    result = {
        "dataset": {
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "coordinate_key": "obsm['spatial']",
        },
        "runs": runs,
    }

    out_dir = Path("benchmarks/results")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "spata2py_benchmark.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# py-SPATA2 Benchmark",
        "",
        f"Dataset: {adata.n_obs} observations x {adata.n_vars} genes.",
        "",
        "| Operation | Seconds | Result shape |",
        "|---|---:|---|",
    ]
    for row in runs:
        lines.append(
            f"| `{row['operation']}` | {row['seconds']:.6f} | `{tuple(row['result_shape'])}` |"
        )
    (out_dir / "spata2py_benchmark.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
