"""Build the recurring "one FancyFolium map, one dataset per switcher entry,
same list of colour-by layers" comparison map every example notebook was
hand-rolling (a ~30-80 line per-city loop of
``background_layer``/``vector_layer`` calls + ``merge_maps``).

This is the 2D, static-column FancyFolium map (folium/Leaflet) used
alongside :func:`.maps.build_map`'s 3D MapLibre + deck.gl viewer -- the two
are complementary, not alternatives: this one shows every layer as a flat
choropleth with FancyFolium's own layer-control panel and popups; that one
extrudes buildings by height and adds the auto-playing showcase tour.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import geopandas as gpd

if TYPE_CHECKING:
    import folium

_DEFAULT_STYLE = {"stroke_color": "black", "weight": 1, "fillOpacity": 0.5}


def _default_label(dataset_id: str) -> str:
    return dataset_id.replace("_", " ").title()


def build_comparison_map(
    datasets: dict[str, gpd.GeoDataFrame],
    layers: list[dict],
    *,
    popup_columns: list[str] | None = None,
    labels: dict[str, str] | None = None,
    style: dict | None = None,
    backgrounds: tuple[str, str] = ("google hybrid", "cartodb light"),
    extra_layers: Callable[[str, gpd.GeoDataFrame, "folium.Map"], "folium.Map"]
    | None = None,
) -> "folium.Map":
    """Build one FancyFolium map with a dataset switcher, each entry getting
    the same list of colour-by layers.

    Args:
        datasets: ``{dataset_id: GeoDataFrame}`` -- one switcher entry each.
        layers: One dict per "Color by" layer, each with:

            - ``column``: Column to colour by.
            - ``layer_name``: Display name in the layer control.
            - ``categorical``: ``True`` for a discrete-category column
              (FancyFolium picks its own swatch colours); otherwise a
              continuous ramp, using ``cmap``/``vmin``/``vmax`` (all
              optional -- FancyFolium falls back to its own defaults/data
              range if omitted).
            - ``active``: Whether this layer starts visible (default
              ``False`` -- exactly one layer should normally have
              ``active=True``, since these are ``overlay=False`` /
              radio-style).
            - ``style``: Per-layer style override (default: *style* /
              :data:`_DEFAULT_STYLE`).
        popup_columns: Columns shown in each feature's popup/tooltip
            (default: every column named in *layers*).
        labels: ``{dataset_id: switcher label}`` override (default:
            title-cased dataset id).
        style: Default per-feature style for layers that don't set their
            own (default: black 1px outline, 50% fill opacity).
        backgrounds: ``(default tile, alternate tile)`` -- both added as a
            radio choice, default first.
        extra_layers: Optional ``(dataset_id, gdf, m) -> m`` callback run
            after the standard layers for each dataset, for anything
            dataset-specific the caller wants to draw on top (e.g. force
            arrows) that doesn't fit the flat column-list model above.

    Returns:
        The merged map (a single :func:`FancyFolium.merge_maps` result if
        *datasets* has more than one entry; otherwise that one map as-is).
    """
    import FancyFolium

    style = style or _DEFAULT_STYLE
    popup_columns = popup_columns or [layer["column"] for layer in layers]

    maps = []
    names = []
    for dataset_id, gdf in datasets.items():
        m = FancyFolium.background_layer(backgrounds[0])
        m = FancyFolium.background_layer(
            backgrounds[1], m=m, overlay=False, active=False
        )

        for layer in layers:
            kwargs = {
                "gdf": gdf,
                "layer_name": layer["layer_name"],
                "column": layer["column"],
                "overlay": False,
                "popup": popup_columns,
                "active": layer.get("active", False),
                "style": layer.get("style", style),
                "m": m,
            }
            if layer.get("categorical", False):
                kwargs["categorical"] = True
            else:
                for key in ("cmap", "vmin", "vmax"):
                    if key in layer:
                        kwargs[key] = layer[key]
            m = FancyFolium.vector_layer(**kwargs)

        if extra_layers is not None:
            m = extra_layers(dataset_id, gdf, m)

        maps.append(m)
        names.append((labels or {}).get(dataset_id, _default_label(dataset_id)))

    return FancyFolium.merge_maps(maps, names) if len(maps) > 1 else maps[0]
