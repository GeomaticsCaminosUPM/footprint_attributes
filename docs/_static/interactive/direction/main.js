// footprint_attributes docs -- "building direction" map. Every building's
// principal axes drawn two ways: from the minimum-rotated bounding box
// (direction.bbox -- dotted shaft, solid tip) and from the principal
// moments of inertia (direction.inertia -- solid throughout), both the
// same thickness so only the line pattern tells them apart. Buildings are
// colored by bearing (a full-circle hue wheel, one full turn per 180
// degrees since an axis has no front/back) -- inertia-based bearing by
// default, switching to bbox-based bearing every few seconds in video mode.
const { GeoJsonLayer } = deck;
const { MapboxOverlay } = deck;
const { PathStyleExtension } = deck;

const DATASETS = {
  guatemala: { label: "Guatemala City (Zona 10)", dir: "../data/guatemala" },
  san_jose: { label: "San José (Mata Redonda)", dir: "../data/san_jose" },
  santo_domingo: { label: "Santo Domingo (Ensanche Quisquella)", dir: "../data/santo_domingo" },
};
const DEFAULT_DATASET = "san_jose";
const DATASET_CYCLE_MS = 10_000;
const PINNED_DATASET = null;
const ACTIVE_DEFAULT_DATASET = PINNED_DATASET || DEFAULT_DATASET;
const BASE_TITLE = document.title;
function updateTitle(datasetId) {
  document.title = `${DATASETS[datasetId].label} – ${BASE_TITLE}`;
}

// Video mode toggles which method's bearing drives the building color,
// independently of (and faster than) the per-city cycle above.
const COLOR_CYCLE_MS = 5_000;

const ARROW_COLOR_BBOX = "#111111";
const ARROW_COLOR_INERTIA = "#111111";
const ARROW_WIDTH = 4;
const ARROW_DASH = [1.4, 1.2];

const DEFAULT_BUILDING_HEIGHT_M = 8;
const SHOWCASE_IDLE_RESUME_MS = 30_000;

function hexToRgb(hex) {
  const v = parseInt(hex.replace("#", ""), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}
// Bearing (0-180 degrees, an axis has no direction) -> a full color wheel,
// so every orientation gets a visually distinct hue.
function hslToRgb(h, s, l) {
  h = ((h % 360) + 360) % 360;
  const c = (1 - Math.abs(2 * l - 1)) * s;
  const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
  const m = l - c / 2;
  let [r, g, b] = [0, 0, 0];
  if (h < 60) [r, g, b] = [c, x, 0];
  else if (h < 120) [r, g, b] = [x, c, 0];
  else if (h < 180) [r, g, b] = [0, c, x];
  else if (h < 240) [r, g, b] = [0, x, c];
  else if (h < 300) [r, g, b] = [x, 0, c];
  else [r, g, b] = [c, 0, x];
  return [Math.round((r + m) * 255), Math.round((g + m) * 255), Math.round((b + m) * 255)];
}
function colorForBearing(bearingDeg) {
  if (typeof bearingDeg !== "number" || Number.isNaN(bearingDeg)) return [140, 140, 140];
  const hue = ((bearingDeg % 180) / 180) * 360;
  return hslToRgb(hue, 0.62, 0.55);
}

const state = {
  datasetId: ACTIVE_DEFAULT_DATASET,
  data: { type: "FeatureCollection", features: [] },
  datasetCenter: null,
  colorMethod: "inertia", // "inertia" | "bbox"
  selectedBuildingId: null,
  showcaseActive: false,
  bboxData: null,
  inertiaData: null,
};

const URL_PARAMS = new URLSearchParams(location.search);
const DARK_THEME = URL_PARAMS.get("theme") === "dark";
if (DARK_THEME) document.body.classList.add("theme-dark");
const MAP_STYLE = DARK_THEME
  ? "https://tiles.openfreemap.org/styles/dark"
  : "https://tiles.openfreemap.org/styles/positron";
const map = new maplibregl.Map({
  container: "map",
  style: MAP_STYLE,
  center: [-69.94, 18.46],
  zoom: 15,
  pitch: 0,
  bearing: 0,
  attributionControl: false,
});

map.addControl(new maplibregl.AttributionControl({ compact: true }));
(() => {
  const attribEl = document.querySelector(".maplibregl-ctrl-attrib");
  if (!attribEl) return;
  attribEl.removeAttribute("open");
  const observer = new MutationObserver(() => attribEl.removeAttribute("open"));
  observer.observe(attribEl, { attributes: true, attributeFilter: ["open"] });
  setTimeout(() => observer.disconnect(), 5000);
})();

const overlay = new MapboxOverlay({ layers: [] });
map.addControl(overlay);
map.setMaxPitch(85);

// This map is always 2D -- no middle-drag pitch handler needed.

// Zoom in past the bounds-fit level so the direction arrows are legible,
// even if that crops some of the AOI off-screen (same reasoning as the
// building_sizes map).
const EXTRA_ZOOM = 1.6;
function fitBoundsCloser(bounds, { duration = 0 } = {}) {
  const camera = map.cameraForBounds(bounds, { padding: 60 });
  if (!camera) {
    map.fitBounds(bounds, { padding: 60, duration });
    return;
  }
  map.easeTo({ center: camera.center, zoom: camera.zoom + EXTRA_ZOOM, duration });
}

function computeBbox(collection) {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  const visit = (coords) => {
    if (Array.isArray(coords) && typeof coords[0] === "number") {
      const [x, y] = coords;
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    } else if (Array.isArray(coords)) {
      coords.forEach(visit);
    }
  };
  for (const f of collection.features) {
    if (f.geometry && "coordinates" in f.geometry) visit(f.geometry.coordinates);
  }
  return [[minX, minY], [maxX, maxY]];
}

function isSelected(feature) {
  return feature.properties?.building_uid === state.selectedBuildingId;
}

function bearingOf(feature) {
  const p = feature.properties ?? {};
  return state.colorMethod === "bbox" ? p.bearing_bbox : p.bearing_inertia;
}
function getFillColor(feature) {
  return [...colorForBearing(bearingOf(feature)), 200];
}
function getLineColor(feature) {
  return isSelected(feature) ? [...hexToRgb("#e2a33f"), 255] : [45, 55, 72, 190];
}
function getLineWidth(feature) {
  return isSelected(feature) ? 3 : 1;
}

let onBuildingClick = () => {};

function buildingTooltip({ object, layer }) {
  if (!object || !layer?.id?.startsWith("buildings-")) return null;
  const p = object.properties ?? {};
  const bearing = typeof bearingOf(object) === "number" ? `${bearingOf(object).toFixed(1)}°` : "n/a";
  return {
    html: `<div><strong>Bearing (${state.colorMethod}):</strong> ${bearing}</div><div class="hint" style="margin-top:4px">Click for full details</div>`,
    className: "deck-tooltip",
  };
}

function renderLayer() {
  const buildings = new GeoJsonLayer({
    id: `buildings-${state.datasetId}`,
    data: state.data,
    filled: true,
    stroked: true,
    pickable: true,
    extruded: false,
    getFillColor,
    getLineColor,
    getLineWidth,
    lineWidthUnits: "pixels",
    lineWidthMinPixels: 1,
    updateTriggers: {
      getFillColor: [state.colorMethod],
      getLineColor: [state.selectedBuildingId],
      getLineWidth: [state.selectedBuildingId],
    },
    onClick: (info) => {
      const id = info.object?.properties?.building_uid;
      if (id !== undefined && id !== null) onBuildingClick(String(id), info.object.properties);
    },
  });

  const layers = [buildings];
  if (state.bboxData) {
    const shaft = state.bboxData.features.filter((f) => f.properties?.part === "shaft");
    const head = state.bboxData.features.filter((f) => f.properties?.part === "head");
    layers.push(
      // bbox axis: dotted shaft...
      new GeoJsonLayer({
        id: `overlay-bbox-shaft-${state.datasetId}`,
        data: { type: "FeatureCollection", features: shaft },
        filled: false,
        stroked: true,
        pickable: false,
        getLineColor: [...hexToRgb(ARROW_COLOR_BBOX), 235],
        getLineWidth: ARROW_WIDTH,
        lineWidthUnits: "pixels",
        lineWidthMinPixels: ARROW_WIDTH,
        extensions: [new PathStyleExtension({ dash: true })],
        getDashArray: ARROW_DASH,
        dashJustified: true,
      }),
      // ...but a solid filled-triangle tip, so it always reads as an arrow.
      new GeoJsonLayer({
        id: `overlay-bbox-head-${state.datasetId}`,
        data: { type: "FeatureCollection", features: head },
        filled: true,
        stroked: true,
        pickable: false,
        getFillColor: [...hexToRgb(ARROW_COLOR_BBOX), 235],
        getLineColor: [...hexToRgb(ARROW_COLOR_BBOX), 235],
        lineWidthUnits: "pixels",
        lineWidthMinPixels: 1,
        lineWidthMaxPixels: 1.5,
      }),
    );
  }
  if (state.inertiaData) {
    layers.push(
      // inertia axis: solid throughout ("normal").
      new GeoJsonLayer({
        id: `overlay-inertia-${state.datasetId}`,
        data: state.inertiaData,
        filled: false,
        stroked: true,
        pickable: false,
        getLineColor: [...hexToRgb(ARROW_COLOR_INERTIA), 235],
        getLineWidth: ARROW_WIDTH,
        lineWidthUnits: "pixels",
        lineWidthMinPixels: ARROW_WIDTH,
      }),
    );
  }
  overlay.setProps({ layers, getTooltip: buildingTooltip });
}

// ---------------------------------------------------------------------------
// Legend.
const LEGEND_ICON_WIDTH = 44;
const LEGEND_ICON_HEIGHT = 18;
function legendArrowSvg(color, dashed) {
  const y = LEGEND_ICON_HEIGHT / 2;
  const x1 = 3;
  const x2 = LEGEND_ICON_WIDTH - 3;
  const headLen = 8;
  const headHalf = 4;
  const dashAttr = dashed ? ` stroke-dasharray="${ARROW_DASH.join(",")}"` : "";
  return `<svg viewBox="0 0 ${LEGEND_ICON_WIDTH} ${LEGEND_ICON_HEIGHT}" width="${LEGEND_ICON_WIDTH}" height="${LEGEND_ICON_HEIGHT}" class="legend-line-icon">
    <line x1="${x1}" y1="${y}" x2="${x2 - headLen}" y2="${y}" stroke="${color}" stroke-width="3"${dashAttr}></line>
    <polygon points="${x2},${y} ${x2 - headLen},${y - headHalf} ${x2 - headLen},${y + headHalf}" fill="${color}"></polygon>
  </svg>`;
}
function legendColorWheelSvg() {
  const size = 60;
  const r = size / 2 - 2;
  const cx = size / 2;
  const cy = size / 2;
  let stops = "";
  for (let deg = 0; deg <= 360; deg += 20) {
    const [r_, g_, b_] = hslToRgb(deg, 0.62, 0.55);
    stops += `<stop offset="${(deg / 360) * 100}%" stop-color="rgb(${r_},${g_},${b_})"></stop>`;
  }
  const gradId = "bearingWheelGrad";
  return `<svg viewBox="0 0 ${size} ${size}" width="${size}" height="${size}">
    <defs><linearGradient id="${gradId}" x1="0" y1="0" x2="1" y2="0">${stops}</linearGradient></defs>
    <circle cx="${cx}" cy="${cy}" r="${r}" fill="url(#${gradId})" stroke="rgba(255,255,255,0.3)"></circle>
  </svg>`;
}
function renderLegend() {
  const container = document.getElementById("legend");
  container.innerHTML = "";
  const heading = document.createElement("h2");
  heading.textContent = "Principal axes";
  container.appendChild(heading);

  const list = document.createElement("ul");
  list.className = "legend-list";
  const items = [
    { color: ARROW_COLOR_BBOX, dashed: true, label: "Bounding box (bbox)" },
    { color: ARROW_COLOR_INERTIA, dashed: false, label: "Principal inertia" },
  ];
  for (const { color, dashed, label } of items) {
    const item = document.createElement("li");
    const icon = document.createElement("span");
    icon.className = "legend-line";
    icon.innerHTML = legendArrowSvg(color, dashed);
    const text = document.createElement("span");
    text.textContent = label;
    item.append(icon, text);
    list.appendChild(item);
  }
  container.appendChild(list);

  const colorHeading = document.createElement("h2");
  colorHeading.style.marginTop = "14px";
  colorHeading.textContent = `Color: bearing (${state.colorMethod})`;
  container.appendChild(colorHeading);
  const wheelRow = document.createElement("div");
  wheelRow.style.display = "flex";
  wheelRow.style.alignItems = "center";
  wheelRow.style.gap = "10px";
  const wheel = document.createElement("span");
  wheel.innerHTML = legendColorWheelSvg();
  const hint = document.createElement("span");
  hint.className = "hint";
  hint.textContent = "Hue = orientation (0-180°)";
  wheelRow.append(wheel, hint);
  container.appendChild(wheelRow);
}

// ---------------------------------------------------------------------------
// Building popup.
function renderBuildingPanel(id, properties) {
  const panel = document.getElementById("building-panel");
  const content = document.getElementById("building-panel-content");
  content.innerHTML = "";
  const heading = document.createElement("h2");
  heading.textContent = `Building ${id}`;
  content.appendChild(heading);

  const table = document.createElement("table");
  table.className = "building-summary";
  const tbody = document.createElement("tbody");
  const rows = [
    ["Height", typeof properties.height === "number" ? `${properties.height.toFixed(1)} m` : "n/a"],
    ["Bearing (bbox)", typeof properties.bearing_bbox === "number" ? `${properties.bearing_bbox.toFixed(1)}°` : "n/a"],
    ["Bearing (inertia)", typeof properties.bearing_inertia === "number" ? `${properties.bearing_inertia.toFixed(1)}°` : "n/a"],
  ];
  for (const [label, value] of rows) {
    const row = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = label;
    const td = document.createElement("td");
    td.textContent = value;
    row.append(th, td);
    tbody.appendChild(row);
  }
  table.appendChild(tbody);
  content.appendChild(table);
  panel.classList.remove("hidden");
}
function closeBuildingPanel() {
  document.getElementById("building-panel").classList.add("hidden");
  state.selectedBuildingId = null;
  renderLayer();
}
function selectBuilding(id, properties) {
  state.selectedBuildingId = id;
  renderLayer();
  renderBuildingPanel(id, properties);
}
onBuildingClick = selectBuilding;

// ---------------------------------------------------------------------------
// Dataset switching.
let datasetDropdown = null;

async function loadDataset(datasetId) {
  const dir = DATASETS[datasetId].dir;
  const buildings = await fetch(`${dir}/buildings.geojson`).then((r) => r.json());
  state.datasetId = datasetId;
  state.data = buildings;

  try {
    state.bboxData = await fetch(`${dir}/overlays/bbox_axis.geojson`).then((r) => r.json());
  } catch {
    state.bboxData = null;
  }
  try {
    state.inertiaData = await fetch(`${dir}/overlays/inertia_axis.geojson`).then((r) => r.json());
  } catch {
    state.inertiaData = null;
  }

  const [[minX, minY], [maxX, maxY]] = computeBbox(buildings);
  state.datasetCenter = { lng: (minX + maxX) / 2, lat: (minY + maxY) / 2 };

  renderLegend();
  renderLayer();
}

async function setDataset(datasetId, { fromShowcase = false } = {}) {
  await loadDataset(datasetId);
  const h1 = document.querySelector(".subtitle");
  if (h1) h1.textContent = DATASETS[datasetId].label;
  updateTitle(datasetId);
  datasetDropdown?.setValue(datasetId);
  const [[minX, minY], [maxX, maxY]] = computeBbox(state.data);
  fitBoundsCloser([[minX, minY], [maxX, maxY]], { duration: fromShowcase ? 0 : 500 });
  if (!fromShowcase) stopShowcase({ resumeAfterIdle: true });
}

function setColorMethod(method, { fromShowcase = false } = {}) {
  state.colorMethod = method;
  const btn = document.getElementById("color-method-toggle");
  if (btn) btn.textContent = method === "bbox" ? "Color: bbox" : "Color: inertia";
  renderLegend();
  renderLayer();
  if (!fromShowcase) stopShowcase({ resumeAfterIdle: true });
}

function createDropdown(container, options, onChange) {
  container.innerHTML = "";
  container.classList.add("dropdown");
  const toggle = document.createElement("button");
  toggle.type = "button";
  toggle.className = "dropdown-toggle";
  const valueEl = document.createElement("span");
  valueEl.className = "dropdown-value";
  toggle.appendChild(valueEl);
  const caret = document.createElement("span");
  caret.className = "dropdown-caret";
  caret.textContent = "▾";
  toggle.appendChild(caret);

  const list = document.createElement("ul");
  list.className = "dropdown-list";
  list.hidden = true;

  function renderOptions(current) {
    list.innerHTML = "";
    for (const option of options) {
      const item = document.createElement("li");
      item.textContent = option.label;
      if (option.value === current) item.classList.add("selected");
      item.addEventListener("click", () => {
        onChange(option.value);
        close();
      });
      list.appendChild(item);
    }
  }
  function close() {
    list.hidden = true;
    container.removeAttribute("data-open");
  }
  function open() {
    list.hidden = false;
    container.setAttribute("data-open", "");
  }
  toggle.addEventListener("click", () => (list.hidden ? open() : close()));
  document.addEventListener("click", (event) => {
    if (!container.contains(event.target)) close();
  });

  container.append(toggle, list);
  renderOptions(options[0]?.value);

  return {
    setValue(value) {
      const option = options.find((o) => o.value === value);
      if (option) valueEl.textContent = option.label;
      renderOptions(value);
    },
  };
}

// ---------------------------------------------------------------------------
// Showcase (video) mode: every COLOR_CYCLE_MS toggle which method's bearing
// drives the building color; every DATASET_CYCLE_MS move to the next city.
// Always 2D -- there's no orbit here, the arrows read best top-down.
let showcaseColorTimer = null;
let showcaseDatasetTimer = null;
let showcaseIdleTimer = null;

function scheduleColorCycle() {
  showcaseColorTimer = setTimeout(() => {
    if (!state.showcaseActive) return;
    setColorMethod(state.colorMethod === "inertia" ? "bbox" : "inertia", { fromShowcase: true });
    scheduleColorCycle();
  }, COLOR_CYCLE_MS);
}
function stopColorCycle() {
  if (showcaseColorTimer !== null) clearTimeout(showcaseColorTimer);
  showcaseColorTimer = null;
}

function scheduleDatasetCycle() {
  if (PINNED_DATASET) return;
  const ids = Object.keys(DATASETS);
  if (ids.length < 2) return;
  showcaseDatasetTimer = setTimeout(() => {
    if (!state.showcaseActive) return;
    const next = ids[(ids.indexOf(state.datasetId) + 1) % ids.length];
    setDataset(next, { fromShowcase: true }).then(scheduleDatasetCycle);
  }, DATASET_CYCLE_MS);
}
function stopDatasetCycle() {
  if (showcaseDatasetTimer !== null) clearTimeout(showcaseDatasetTimer);
  showcaseDatasetTimer = null;
}

function startShowcase() {
  if (state.showcaseActive) return;
  state.showcaseActive = true;
  scheduleColorCycle();
  scheduleDatasetCycle();
}
function stopShowcase({ resumeAfterIdle = true } = {}) {
  stopColorCycle();
  stopDatasetCycle();
  state.showcaseActive = false;
  if (showcaseIdleTimer !== null) clearTimeout(showcaseIdleTimer);
  if (resumeAfterIdle) showcaseIdleTimer = setTimeout(startShowcase, SHOWCASE_IDLE_RESUME_MS);
}
function registerUserInteraction() {
  stopShowcase({ resumeAfterIdle: true });
}
["dragstart", "zoomstart", "rotatestart", "pitchstart"].forEach((event) => {
  map.on(event, (e) => {
    if (e.originalEvent) registerUserInteraction();
  });
});
map.on("click", () => registerUserInteraction());

function resetOrientation() {
  map.easeTo({ bearing: 0, pitch: 0, duration: 500 });
}

async function bootstrap() {
  datasetDropdown = createDropdown(
    document.getElementById("dataset-select"),
    Object.entries(DATASETS).map(([value, { label }]) => ({ value, label })),
    (value) => setDataset(value),
  );

  await loadDataset(ACTIVE_DEFAULT_DATASET);
  const h1 = document.querySelector(".subtitle");
  if (h1) h1.textContent = DATASETS[ACTIVE_DEFAULT_DATASET].label;
  updateTitle(ACTIVE_DEFAULT_DATASET);
  datasetDropdown.setValue(ACTIVE_DEFAULT_DATASET);
  const [[minX, minY], [maxX, maxY]] = computeBbox(state.data);
  fitBoundsCloser([[minX, minY], [maxX, maxY]], { duration: 0 });

  document.getElementById("controls-toggle").addEventListener("click", (event) => {
    document.getElementById("controls-fields").classList.toggle("hidden");
    event.currentTarget.classList.toggle("collapsed");
  });
  document.getElementById("legend-toggle").addEventListener("click", (event) => {
    document.getElementById("legend").classList.toggle("hidden");
    event.currentTarget.classList.toggle("active");
  });
  document.getElementById("reorient-toggle").addEventListener("click", () => {
    resetOrientation();
    registerUserInteraction();
  });
  document.getElementById("resume-showcase-toggle").addEventListener("click", () => {
    if (state.showcaseActive) stopShowcase({ resumeAfterIdle: false });
    else startShowcase();
  });
  document.getElementById("settings-toggle").addEventListener("click", () => {
    document.getElementById("settings-panel").classList.toggle("hidden");
  });
  document.getElementById("settings-panel-close").addEventListener("click", () => {
    document.getElementById("settings-panel").classList.add("hidden");
  });
  document.getElementById("building-panel-close").addEventListener("click", closeBuildingPanel);
  document.getElementById("color-method-toggle")?.addEventListener("click", () => {
    setColorMethod(state.colorMethod === "inertia" ? "bbox" : "inertia");
  });

  startShowcase();
}

map.on("load", () => {
  bootstrap().catch((error) => {
    console.error(error);
    const controls = document.getElementById("controls-body");
    const message = document.createElement("p");
    message.className = "error";
    message.textContent = `Failed to load: ${error.message}. Serve this folder over HTTP so fetch() can read the data files.`;
    controls.appendChild(message);
  });
});
