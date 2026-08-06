#!/usr/bin/env -S uv run --python 3.12 --with jump_portrait --with 'zarr>=3' --with pandas --with numpy -s
"""Rebuild jump_rab30_mini/ -- the demo bundle that demo.sql queries.

Pulls the RAB30 ORF and CRISPR sites used by the JUMP hub's
`14_display_perturbation_images` notebook, center-crops each to 384x384 at
native resolution (keeps real pixel statistics, including saturation), and
writes a Zarr store plus two CSV lookup tables.

    ./make_bundle.py      # -> jump_rab30_mini/
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from jump_portrait.fetch import get_item_location_metadata, get_jump_image

CHANNELS = ["AGP", "DNA", "ER", "Mito", "RNA"]
KEYS = ("Source", "Batch", "Plate", "Well", "Site")
CROP = 384
OUT = Path("jump_rab30_mini")

# Two real RAB30 sites: ORF is JCP2022_9*, CRISPR is JCP2022_8*.
# Sorted because get_item_location_metadata row order is not stable.
info = get_item_location_metadata("RAB30").to_pandas().sort_values(
    [f"Metadata_{k}" for k in KEYS]
)
sites = [
    {"idx": i, "modality": m, **{k: r[f"Metadata_{k}"] for k in KEYS}}
    for i, (m, p) in enumerate((("ORF", "JCP2022_9"), ("CRISPR", "JCP2022_8")))
    for r in [info[info.Metadata_JCP2022.str.startswith(p)].iloc[0]]
]

ctr = lambda a: a[
    (a.shape[0] - CROP) // 2 : (a.shape[0] + CROP) // 2,
    (a.shape[1] - CROP) // 2 : (a.shape[1] + CROP) // 2,
]
stack = np.stack([
    np.stack([ctr(get_jump_image(*[s[k] for k in KEYS[:4]], c, s["Site"])) for c in CHANNELS])
    for s in sites
])
print("stack", stack.shape, stack.dtype)

shutil.rmtree(OUT, ignore_errors=True)
OUT.mkdir()
# gzip, not zarr-python's default zstd: duckdb_zarr reads raw + gzip only.
zarr.open_group(OUT / "images.zarr", mode="w").create_array(
    "images", shape=stack.shape, chunks=(1, 1, CROP, CROP), dtype="u2",
    compressors=zarr.codecs.GzipCodec(),
)[:] = stack

pd.DataFrame(sites).to_csv(OUT / "sites.csv", index=False)
pd.DataFrame(enumerate(CHANNELS), columns=["idx", "channel"]).to_csv(OUT / "channels.csv", index=False)
print("wrote", OUT)
