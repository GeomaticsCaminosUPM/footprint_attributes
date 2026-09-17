#!/usr/bin/env python3
"""Build the interactive pilot-region map embedded in docs/examples.rst,
using this package's own footprint_attributes.visualization.build_map
against the three example pilot-region datasets under examples/data/.

See examples/generate_interactive_map.py for the same thing set up as a
standalone worked example.

Run with this package's own venv (needs the ``vis`` extra):
    .venv/bin/python docs/generate_map_data.py
"""

from __future__ import annotations

import os

import geopandas as gpd

from footprint_attributes.visualization import build_map

EXAMPLES_DATA = os.path.join(os.path.dirname(__file__), "..", "examples", "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "_static", "maps")

DATASETS = {
    "guatemala": {
        "file": "guatemala_pilot_region.gpkg",
        "label": "Guatemala — Zona 10",
    },
    "san_jose": {
        "file": "san_jose_pilot_region.gpkg",
        "label": "San José — Mata Redonda",
    },
    "santo_domingo": {
        "file": "santo_domingo_pilot_region.gpkg",
        "label": "Santo Domingo — Ensanche Quisquella",
    },
}


def main() -> None:
    datasets = {
        name: gpd.read_file(os.path.join(EXAMPLES_DATA, info["file"]))
        for name, info in DATASETS.items()
    }
    labels = {name: info["label"] for name, info in DATASETS.items()}

    build_map(
        datasets,
        OUTPUT_DIR,
        labels=labels,
        default_dataset="guatemala",
        title="Pilot regions — footprint_attributes",
    )
    print(f"wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
