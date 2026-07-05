"""
Building position within the urban block.

This module classifies each building footprint into one of five relative-position
categories based on a "contact force" analogy: forces proportional to the shared
wall area are computed for every touching pair of footprints, and their resultant
magnitude, confinement ratio, and weighted angle are used to assign the class.

The public ``position`` object is both callable (runs the full pipeline) and
exposes sub-functions as attributes:

Typical usage
-------------
>>> import geopandas as gpd
>>> from footprint_attributes import position
>>>
>>> footprints = gpd.read_file("footprints.gpkg")
>>>
>>> # Full pipeline — computes forces and classifies
>>> result = position(footprints)
>>> result[['relativePosition', 'contact_force', 'contact_confinementRatio']]
>>>
>>> # Classification only — reuses pre-computed force columns when present
>>> labels = position.relative_position(result)          # fast: columns exist
>>> labels = position.relative_position(footprints)      # slow: recomputes forces
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd

from .geometry import (
    ensure_projected,
    to_gdf,
    validate_geodataframe,
    fill_holes,
    calc_inertia_z,
    select_touching_edges,
    explode_edges,
    edge_normal,
    edge_momentum,
)
from .config import POSITION_DEFAULTS


# ─────────────────────────────────────────────────────────────────────────────
# Contact force helpers
# ─────────────────────────────────────────────────────────────────────────────


def resultant_angle(
    gdf: gpd.GeoDataFrame,
    vector_column: str,
    id_column: str,
) -> gpd.GeoDataFrame:
    """Compute angle between each force vector and the resultant per building.

    Args:
        gdf: GeoDataFrame with vector_column and id_column.
        vector_column: Name of column containing 2-D force vectors.
        id_column: Name of column with building ID.

    Returns:
        GeoDataFrame with added 'angle' column.
    """
    gdf = gdf.copy()

    # Compute resultant per building
    resultants = gdf.groupby(id_column)[vector_column].apply(
        lambda vecs: np.sum(np.array(vecs.tolist()), axis=0)
    )

    # Compute angle between each force and resultant
    def angle_to_resultant(row):
        res = resultants[row[id_column]]
        force = row[vector_column]
        if np.linalg.norm(res) < 1e-12 or np.linalg.norm(force) < 1e-12:
            return 0.0
        cos_angle = np.dot(force, res) / (
            np.linalg.norm(force) * np.linalg.norm(res) + 1e-12
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.arccos(cos_angle)

    gdf["angle"] = gdf.apply(angle_to_resultant, axis=1)
    return gdf


# ─────────────────────────────────────────────────────────────────────────────
# Main contact forces function
# ─────────────────────────────────────────────────────────────────────────────


def contact_forces_df(
    geoms: gpd.GeoDataFrame,
    buffer: float = 0,
    height_column: str | None = None,
    minRadius: float = 0,
) -> pd.DataFrame:
    """Compute direction-aware contact-force metrics for each building footprint.

    For every touching pair of buildings a virtual unit pressure (1 Pa) is
    applied perpendicularly to each shared wall segment.  The force magnitude
    is proportional to the wall area (segment length × building height).
    The resulting per-building metrics are:

    - **force** – magnitude of the resultant contact force, normalised by the
      square root of the footprint area.
    - **confinementRatio** – fraction of the total force that is
      "cancelled out" by opposing forces (0 = fully lateral, 1 = fully
      enclosed).
    - **angularAcc** – angular acceleration proxy (momentum × area / I_z).
    - **angle** – force-weighted mean angle between individual forces and the
      resultant (in radians).
    - **height** – building height used in the calculation.

    Args:
        geoms: GeoDataFrame of building footprints (any CRS).
        buffer: Contact detection buffer in metres.  Use the pixel size of the
            aerial image used to digitise the footprints (e.g., 0.15 for
            15 cm resolution imagery).
        height_column: Name of a column in *geoms* containing building heights
            in metres.  If ``None`` all buildings get height = 1.
        minRadius: Fraction of the equivalent circle radius below which a
            force's momentum is only counted when it *reduces* the net torque
            (conservative filter for very close contacts).

    Returns:
        DataFrame indexed like *geoms* with columns:
        ``height``, ``force``, ``confinementRatio``, ``angularAcc``,
        ``angle``.  Buildings with no neighbours have zeros in all metrics
        except *height*.
    """
    geoms = to_gdf(geoms)
    geoms_orig = geoms[[geoms.geometry.name]].copy()

    # ------------------------------------------------------------------
    # Preparation
    # ------------------------------------------------------------------
    geoms["geom_id"] = geoms.index.copy()
    geoms.geometry = geoms.geometry.force_2d()
    geoms = ensure_projected(geoms)
    validate_geodataframe(geoms, context="contact_forces_df")

    if height_column is None:
        geoms["height"] = 1.0
    else:
        geoms["height"] = geoms[height_column].astype(float)
        nan_mask = geoms["height"].isna()
        if nan_mask.any():
            bad = list(geoms.index[nan_mask][:10])
            raise ValueError(
                f"[contact_forces_df] Found {nan_mask.sum()} null value(s) in "
                f"height_column '{height_column}' at index positions {bad}"
                f"{'…' if nan_mask.sum() > 10 else ''}. A missing height silently "
                "zeroes that building's contact force/angularAcc, making it look "
                "'isolated' regardless of real neighbours. Fill or drop null "
                "heights before calling contact_forces_df()/position()."
            )

    geoms["inertia"] = calc_inertia_z(geoms.geometry)
    geoms["area"] = geoms.geometry.area
    geoms["centroid"] = geoms.geometry.centroid
    orig_height = geoms["height"].copy()

    # ------------------------------------------------------------------
    # Select touching edge segments and compute forces
    # ------------------------------------------------------------------
    # Interior holes (courtyards) are irrelevant to neighbour contact and
    # would otherwise show up as spurious "touching" edges (their ring lies
    # inside the eroded union just like a real shared wall).  Fill them
    # before detecting touching edges; area/centroid/inertia above were
    # already computed from the real (unfilled) footprint.
    geoms.geometry = fill_holes(geoms.geometry)
    geoms = select_touching_edges(geoms, buffer=buffer)
    geoms = explode_edges(geoms, min_length=buffer)

    if len(geoms) == 0:
        # No building in the batch touches any other (e.g. a single isolated
        # footprint, or a set of footprints with no shared walls at all).
        result = geoms_orig.copy()
        result["height"] = orig_height
        for col in ("force", "confinementRatio", "angularAcc", "angle"):
            result[col] = 0.0
        return result[
            ["height", "force", "confinementRatio", "angularAcc", "angle"]
        ].astype(float)

    # Compute edge normals (force vectors). Use the exploded two-point
    # "edges" column explicitly, not r.geometry -- explode_edges() drops the
    # stale pre-split geometry column precisely to prevent this.
    normal_results = geoms.apply(
        lambda r: edge_normal(r["edges"], scale=r["height"]),
        axis=1,
        result_type="expand",
    )
    geoms["edge_center"] = normal_results[0]
    geoms["normal_vector"] = normal_results[1]

    # Compute momentum for each force
    geoms["momentum"] = geoms.apply(
        lambda r: edge_momentum(
            r["edge_center"],
            r["normal_vector"],
            r["centroid"],
            min_dist=minRadius * np.sqrt(r["area"] / np.pi),
        ),
        axis=1,
    )

    geoms["abs_force"] = geoms["normal_vector"].apply(np.linalg.norm)

    # Angle between each force and the resultant
    geoms = resultant_angle(geoms, vector_column="normal_vector", id_column="geom_id")
    geoms["angle"] = geoms["angle"] * 2 * geoms["abs_force"]  # weighted contribution

    # ------------------------------------------------------------------
    # Aggregate per building
    # ------------------------------------------------------------------
    agg = geoms.groupby("geom_id").agg(
        height=("height", "first"),
        normal_vector=("normal_vector", "sum"),
        abs_force=("abs_force", "sum"),
        angle=("angle", "sum"),
        momentum=("momentum", "sum"),
        inertia=("inertia", "first"),
        area=("area", "first"),
    )

    res_force = agg["normal_vector"].apply(np.linalg.norm)
    agg["force"] = res_force / np.sqrt(agg["area"] + 1e-12)
    agg["confinementRatio"] = (agg["abs_force"] - res_force) / (
        agg["abs_force"] + 1e-12
    )
    agg["momentum_mag"] = agg["momentum"].apply(
        lambda m: np.abs(m).min() if isinstance(m, np.ndarray) else abs(m)
    )
    agg["angle"] = agg["angle"] / (agg["abs_force"] + 1e-12)
    agg["angularAcc"] = agg["momentum_mag"] / (agg["inertia"] + 1e-12) * agg["area"]

    # ------------------------------------------------------------------
    # Merge back and fill missing (isolated buildings)
    # ------------------------------------------------------------------
    result = geoms_orig.merge(
        agg[["height", "force", "confinementRatio", "angularAcc", "angle"]],
        left_index=True,
        right_index=True,
        how="left",
    )
    for col in ("force", "confinementRatio", "angularAcc", "angle"):
        result[col] = result[col].fillna(0.0).astype(float)

    return result[
        ["height", "force", "confinementRatio", "angularAcc", "angle"]
    ].astype(float)


# ─────────────────────────────────────────────────────────────────────────────
# Callable class — exposes position() and position.relative_position()
# ─────────────────────────────────────────────────────────────────────────────

# Mapping from the prefixed column names stored by position() to the plain
# column names expected by the classification logic.
_FORCE_COL_MAP: dict[str, str] = {
    "contact_force": "force",
    "contact_confinementRatio": "confinementRatio",
    "contact_angularAcc": "angularAcc",
    "contact_angle": "angle",
}

# The plain column names required for classification.
_REQUIRED_PLAIN: frozenset[str] = frozenset(_FORCE_COL_MAP.values())


class _Position:
    """Callable namespace for building position attributes.

    Call the object itself for the full pipeline, or use the
    :meth:`relative_position` method to (re-)classify an existing GeoDataFrame
    without recomputing forces.

    Examples
    --------
    >>> result = position(footprints)                        # full pipeline
    >>> labels = position.relative_position(result)          # reuses columns
    >>> labels = position.relative_position(footprints)      # computes forces
    """

    # ------------------------------------------------------------------
    # __call__ — full pipeline
    # ------------------------------------------------------------------

    def __call__(
        self,
        footprints_gdf: gpd.GeoDataFrame,
        columns: list[str] | None = None,
        buffer: float = POSITION_DEFAULTS["buffer"],
        height_column: str | None = None,
        minRadius: float = POSITION_DEFAULTS["minRadius"],
        minForce: float = POSITION_DEFAULTS["minForce"],
        minAngle: float = POSITION_DEFAULTS["minAngle"],
        minConfinement: float = POSITION_DEFAULTS["minConfinement"],
        minAngularAcc: float = POSITION_DEFAULTS["minAngularAcc"],
    ) -> gpd.GeoDataFrame:
        """Compute all position attributes for buildings.

        Args:
            footprints_gdf: GeoDataFrame of building footprints.
            columns: Column names to keep in output; ``None`` returns all.
            buffer: Contact detection buffer (metres).
            height_column: Column with building heights in metres.
            minRadius: Minimum radius fraction for momentum filtering.
            minForce: Force threshold for the *lateral* class.
            minAngle: Weighted-angle threshold (rad) for the *corner* class.
            minConfinement: Confinement threshold for the *confined* class.
            minAngularAcc: Angular-acceleration threshold for the *torque* class.

        Returns:
            GeoDataFrame with position columns added:
            ``contact_force``, ``contact_confinementRatio``,
            ``contact_angularAcc``, ``contact_angle``, ``relativePosition``.
        """
        gdf = ensure_projected(to_gdf(footprints_gdf))

        # Compute contact forces
        forces = contact_forces_df(
            gdf,
            buffer=buffer,
            height_column=height_column,
            minRadius=minRadius,
        )

        # Store force columns with the prefixed API names. `contact_height`
        # is kept too so a later `relative_position()` call reusing these
        # prefixed columns can still undo the height scaling of
        # contact_force/contact_angularAcc when classifying (see _classify).
        gdf["contact_force"] = forces["force"]
        gdf["contact_confinementRatio"] = forces["confinementRatio"]
        gdf["contact_angularAcc"] = forces["angularAcc"]
        gdf["contact_angle"] = forces["angle"]
        gdf["contact_height"] = forces["height"]

        # Classify (forces already has the plain column names)
        gdf["relativePosition"] = self._classify(
            forces,
            minAngularAcc=minAngularAcc,
            minConfinement=minConfinement,
            minAngle=minAngle,
            minForce=minForce,
        )

        if columns is not None:
            available = [c for c in columns if c in gdf.columns or c == "geometry"]
            gdf = gdf[available]

        return gdf

    # ------------------------------------------------------------------
    # relative_position — classification only, reuses columns when present
    # ------------------------------------------------------------------

    def relative_position(
        self,
        footprints_gdf: gpd.GeoDataFrame | pd.DataFrame,
        minAngularAcc: float = POSITION_DEFAULTS["minAngularAcc"],
        minConfinement: float = POSITION_DEFAULTS["minConfinement"],
        minAngle: float = POSITION_DEFAULTS["minAngle"],
        minForce: float = POSITION_DEFAULTS["minForce"],
        buffer: float = POSITION_DEFAULTS["buffer"],
        height_column: str | None = None,
        minRadius: float = POSITION_DEFAULTS["minRadius"],
    ) -> list[str]:
        """Classify each building into a relative-position category.

        Pre-computed force columns are reused when available, avoiding an
        expensive re-computation.  Two column naming conventions are accepted:

        * **Plain** (output of :func:`contact_forces_df`):
          ``force``, ``confinementRatio``, ``angularAcc``, ``angle``.
        * **Prefixed** (output of :meth:`__call__`):
          ``contact_force``, ``contact_confinementRatio``,
          ``contact_angularAcc``, ``contact_angle``.

        If neither set is present the contact forces are computed from the
        footprint geometries before classification.

        Categories (in decreasing priority):

        1. **"torque"** – Confined or corner building with high angular
           acceleration, indicating potentially damaging torsional response.
        2. **"confined"** – Touches neighbours on both lateral sides
           (high confinement ratio).
        3. **"corner"** – Touches at two perpendicular sides (high resultant
           force AND high weighted angle).
        4. **"lateral"** – Touches on one side (high resultant force only).
        5. **"isolated"** – No touching neighbours.

        Args:
            footprints_gdf: GeoDataFrame of footprints **or** a forces
                DataFrame.  When force columns are already present (plain or
                prefixed) they are reused directly.
            minAngularAcc: Angular-acceleration threshold for the *torque*
                class.  Default: ``POSITION_DEFAULTS["minAngularAcc"]``.
            minConfinement: Confinement-ratio threshold for the *confined*
                class.  Default: ``POSITION_DEFAULTS["minConfinement"]``.
            minAngle: Weighted-angle threshold (rad) for the *corner* class.
                Default: ``POSITION_DEFAULTS["minAngle"]``.
            minForce: Resultant-force threshold for the *lateral* class.
                Default: ``POSITION_DEFAULTS["minForce"]``.
            buffer: Forwarded to :func:`contact_forces_df` when forces must be
                computed.
            height_column: Forwarded to :func:`contact_forces_df`.
            minRadius: Forwarded to :func:`contact_forces_df`.

        Returns:
            List of category strings aligned with the input DataFrame index.
        """
        df = footprints_gdf.copy()

        if _REQUIRED_PLAIN.issubset(df.columns):
            # Plain columns already present — use them directly. `height` is
            # carried along too when available, so _classify can undo the
            # height scaling of force/angularAcc instead of assuming height=1.
            cols = list(_REQUIRED_PLAIN)
            if "height" in df.columns:
                cols.append("height")
            forces = df[cols]

        elif _FORCE_COL_MAP.keys() <= set(df.columns):
            # Prefixed columns present (output of __call__) — rename to plain.
            cols = list(_FORCE_COL_MAP)
            rename = dict(_FORCE_COL_MAP)
            if "contact_height" in df.columns:
                cols.append("contact_height")
                rename["contact_height"] = "height"
            forces = df[cols].rename(columns=rename)

        else:
            # No force columns found — compute from footprint geometries.
            forces = contact_forces_df(
                df,
                buffer=buffer,
                height_column=height_column,
                minRadius=minRadius,
            )

        return self._classify(
            forces,
            minAngularAcc=minAngularAcc,
            minConfinement=minConfinement,
            minAngle=minAngle,
            minForce=minForce,
        )

    # ------------------------------------------------------------------
    # Internal: pure classification logic (operates on plain column names)
    # ------------------------------------------------------------------

    @staticmethod
    def _classify(
        forces: pd.DataFrame,
        *,
        minAngularAcc: float,
        minConfinement: float,
        minAngle: float,
        minForce: float,
    ) -> list[str]:
        """Apply threshold rules to a DataFrame with plain force columns.

        Args:
            forces: DataFrame with columns ``force``, ``angle``,
                ``confinementRatio``, ``angularAcc``, and optionally
                ``height`` (used to undo the height scaling below).
            minAngularAcc: Threshold for the *torque* class.
            minConfinement: Threshold for the *confined* class.
            minAngle: Threshold for the *corner* class.
            minForce: Threshold for the *lateral* class.

        Returns:
            List of category strings aligned with *forces*.
        """
        out = forces.copy()
        out["relativePosition"] = "isolated"

        # `force`/`angularAcc` scale linearly with building height (a taller
        # shared wall really does carry more contact force), but minForce/
        # minAngularAcc are calibrated for height=1 (see POSITION_DEFAULTS).
        # Compare against the height=1-equivalent value so classification
        # doesn't depend on whether/what height_column was supplied;
        # confinementRatio/angle are already height-invariant ratios and
        # need no such correction.
        if "height" in out.columns:
            height = out["height"].replace(0, np.nan).fillna(1.0)
        else:
            height = 1.0
        force_for_classification = out["force"] / height
        angular_acc_for_classification = out["angularAcc"] / height

        out.loc[force_for_classification > minForce, "relativePosition"] = "lateral"

        out.loc[
            (out["angle"] > minAngle) & (out["relativePosition"] == "lateral"),
            "relativePosition",
        ] = "corner"

        out.loc[out["confinementRatio"] > minConfinement, "relativePosition"] = (
            "confined"
        )

        out.loc[
            out["relativePosition"].isin(["corner", "confined"])
            & (angular_acc_for_classification > minAngularAcc),
            "relativePosition",
        ] = "torque"

        return list(out["relativePosition"])

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:  # pragma: no cover
        return (
            "position  (callable)\n"
            "  position(footprints_gdf, ...)           → full pipeline\n"
            "  position.relative_position(gdf, ...)    → classification only"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Module-level singleton — this is what ``from footprint_attributes import
# position`` gives you.
# ─────────────────────────────────────────────────────────────────────────────

#: The public ``position`` object.  Call it or use its sub-methods.
position = _Position()
