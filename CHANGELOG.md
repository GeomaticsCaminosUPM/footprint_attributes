# Changelog

## Unreleased

- License changed from GPL-3.0-or-later to the European Union Public Licence
  v1.2 (EUPL-1.2), the license of the other repositories of the project.

## 1.1.0

Changes the values of `contact_confinementRatio` and some `blockPosition`
classes. Recompute stored position features when upgrading.

- **Breaking:** the position class column is renamed from `relativePosition`
  to `blockPosition` (and `position.relative_position()` to
  `position.block_position()`).
- **Breaking:** new `contact_confinementRatio`. Each wall force is reduced by
  the resultant projected on its own direction, clipped to [0, |f|], and the
  clipped values are summed and divided by the total force. Only walls
  pushed against by walls on the other side count: perpendicular walls (a
  corner) and several segments on the same side give 0, two opposite walls
  of lengths a and b give 2 min(a, b) / (a + b). The old ratio,
  (sum |f| - |sum f|) / sum |f|, gave 0.29 to two equal perpendicular walls
  and let an unopposed back wall hide a pair of opposite walls, so buildings
  touched on three sides could be classified as corners.
- **Breaking:** `minConfinement` default changes from 0.5 to 0.2.
- **Breaking:** the default contact `buffer` changes from 0 to 0.1 m, the value the
  documentation already described. With 0, touching footprints whose edges are not
  exactly collinear could lose most of a shared wall.
- Docs, interactive maps and figure file names say "block position" everywhere.
- `run(..., n_jobs=...)` computes shape indices in chunks and position in
  groups of touching buildings on several CPU cores (about 6x faster with
  8 workers). Shape indices and classes are identical to a serial run.
- Contact detection no longer fails with a GEOS "side location conflict" on nearly
  coincident edges (robust snap-rounding overlay as a fallback).
- EC8 radius ratio: the label direction of the compliance check is fixed.
- Fix of a circle-centre bug in `geometry.py`.

## 1.0.0

First release.
