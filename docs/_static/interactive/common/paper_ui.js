// Map chrome of footprint_attributes' MapLibre + deck.gl pages, in the look
// of the structural-system paper's web map. Load it after maplibre-gl.js and
// before the page's own main.js: it wraps maplibregl.Map so that the map the
// page creates gets
//   - the basemap chosen here (light, dark or satellite; ?basemap= or ?theme=dark),
//   - zoom/compass buttons and one row of toggles at the top right: dark mode,
//     3D, presentation mode, plots, auto tour and reset orientation (the
//     page's own buttons are moved into that row, so their handlers still work),
//   - a basemap selector in the left panel,
//   - presentation mode (P key or ?presentation): only the legend and plots.
(() => {
  const BASEMAPS = {
    light: { label: "Light", tone: "light", style: "https://tiles.openfreemap.org/styles/positron" },
    dark: { label: "Dark", tone: "dark", style: "https://tiles.openfreemap.org/styles/dark" },
    satellite: {
      label: "Satellite", tone: "dark",
      style: {
        version: 8,
        sources: { sat: { type: "raster", tileSize: 256, attribution: "&copy; Esri, Maxar, Earthstar Geographics",
          tiles: ["https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"] } },
        layers: [{ id: "sat", type: "raster", source: "sat" }],
      },
    },
  };
  const ICONS = {
    present: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="3" y="4" width="18" height="12" rx="1"/><path d="M12 16v4M8 20h8"/></svg>',
    moon: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>',
    sun: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>',
  };
  // The page's own buttons, in the order they appear in the row.
  const PAGE_BUTTONS = ["view-3d-toggle", "color-method-toggle", "charts-toggle", "resume-showcase-toggle", "reorient-toggle"];

  const params = new URLSearchParams(window.location.search);
  const $ = (id) => document.getElementById(id);
  const ui = {
    map: null,
    basemap: BASEMAPS[params.get("basemap")] ? params.get("basemap") : params.get("theme") === "dark" ? "dark" : "light",
  };
  window.paperUI = ui;
  // The tab title before the page appends the dataset name to it.
  const pageTitle = document.querySelector('meta[name="map-title"]')?.content
    || document.title.replace(/\s*[–—-]\s*footprint_attributes\s*$/, "").replace(/^footprint_attributes\s*[–—-]\s*/, "");

  function isDark() { return document.body.classList.contains("dark"); }

  // The basemap's own extruded buildings would hide the page's buildings.
  function hideBasemapExtrusions(map) {
    for (const layer of map.getStyle().layers || []) {
      if (layer.type === "fill-extrusion") map.setLayoutProperty(layer.id, "visibility", "none");
    }
  }

  function setBasemap(id) {
    if (!BASEMAPS[id] || ui.basemap === id) return;
    ui.basemap = id;
    const select = $("basemap");
    if (select) select.value = id;
    ui.map.setStyle(BASEMAPS[id].style, { diff: false });
    ui.map.once("style.load", () => hideBasemapExtrusions(ui.map));
  }

  // Dark or light panels. The basemap follows (light <-> dark), except the satellite view.
  function setDark(on) {
    document.body.classList.toggle("dark", on);
    document.body.classList.toggle("theme-dark", on);
    const b = $("btnTheme");
    if (b) { b.innerHTML = on ? ICONS.sun : ICONS.moon; b.title = on ? "Light mode" : "Dark mode"; }
    if (ui.basemap !== "satellite") setBasemap(on ? "dark" : "light");
  }

  function setPresentation(on) {
    document.body.classList.toggle("presentation", on);
    $("btnPresent")?.classList.toggle("on", on);
  }

  // The plots and building panels sit right under the map buttons.
  function placePanels() {
    const ctrl = document.querySelector(".maplibregl-ctrl-top-right");
    if (!ctrl) return;
    const top = `${ctrl.getBoundingClientRect().bottom + 8}px`;
    for (const id of ["charts-panel", "building-panel"]) if ($(id)) $(id).style.top = top;
  }

  class ToggleControl {
    onAdd() {
      this.el = document.createElement("div");
      this.el.className = "maplibregl-ctrl maplibregl-ctrl-group";
      const add = (id, title, html, onClick) => {
        const b = document.createElement("button");
        b.id = id; b.type = "button"; b.title = title; b.className = "toggle"; b.innerHTML = html;
        b.addEventListener("click", onClick);
        this.el.appendChild(b);
      };
      add("btnTheme", "Dark mode", ICONS.moon, () => setDark(!isDark()));
      const addPresent = () => add("btnPresent", "Presentation mode (P)", ICONS.present,
        () => setPresentation(!document.body.classList.contains("presentation")));
      let presentAdded = false;
      for (const id of PAGE_BUTTONS) {
        const b = $(id);
        if (!b) continue;
        if (id === "charts-toggle" && !$("charts-panel")) { b.remove(); continue; }
        b.classList.remove("icon-button");
        b.classList.add("toggle");
        this.el.appendChild(b);
        if (id === "view-3d-toggle") { addPresent(); presentAdded = true; }
      }
      if (!presentAdded) addPresent();
      return this.el;
    }
    onRemove() { this.el.remove(); }
  }

  function buildPanel() {
    const body = $("controls-body");
    if (!body) return;
    // Controls always visible, as in the paper's map (no collapse arrow).
    $("controls-fields")?.classList.remove("hidden");
    $("controls")?.classList.remove("collapsed");
    const titleRow = body.querySelector(".panel-title-row");
    if (titleRow && !titleRow.querySelector(".map-title")) {
      const h1 = document.createElement("div");
      h1.className = "map-title";
      h1.textContent = pageTitle;
      titleRow.prepend(h1);
    }
    const fields = $("controls-fields");
    if (fields && !$("basemap")) {
      const label = document.createElement("label");
      label.className = "ctl";
      label.textContent = "Basemap";
      const select = document.createElement("select");
      select.id = "basemap";
      select.innerHTML = Object.entries(BASEMAPS).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join("");
      select.value = ui.basemap;
      select.addEventListener("change", (e) => setBasemap(e.target.value));
      label.appendChild(select);
      fields.appendChild(label);
    }
  }

  function setup(map) {
    ui.map = map;
    map.addControl(new maplibregl.NavigationControl({ visualizePitch: true }), "top-right");
    map.addControl(new ToggleControl(), "top-right");
    buildPanel();
    if (BASEMAPS[ui.basemap].tone === "dark" && ui.basemap !== "satellite") {
      document.body.classList.add("dark", "theme-dark");
      $("btnTheme").innerHTML = ICONS.sun;
      $("btnTheme").title = "Light mode";
    }
    map.once("style.load", () => hideBasemapExtrusions(map));
    if (params.has("presentation")) setPresentation(true);
    document.addEventListener("keydown", (e) => {
      const typing = e.target instanceof HTMLSelectElement || e.target instanceof HTMLInputElement;
      if ((e.key === "p" || e.key === "P") && !typing) setPresentation(!document.body.classList.contains("presentation"));
    });
    window.addEventListener("resize", placePanels);
    requestAnimationFrame(placePanels);
    map.once("load", placePanels);
  }

  // The page creates its map after this script: give it this basemap and the chrome.
  const BaseMap = maplibregl.Map;
  class PaperMap extends BaseMap {
    constructor(options) {
      super({ ...options, style: BASEMAPS[ui.basemap].style });
      setup(this);
    }
  }
  maplibregl.Map = PaperMap;
  ui.setBasemap = setBasemap;
  ui.setDark = setDark;
  ui.setPresentation = setPresentation;
})();
