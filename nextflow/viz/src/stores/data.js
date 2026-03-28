/**
 * Reactive data stores for the microscape amplicon pipeline visualization.
 * Uses Svelte 5 runes ($state) for fine-grained reactivity.
 */

// ── Core data ───────────────────────────────────────────────────────────────

/** Sample coordinates and metadata: { id, x, y, reads, richness, ... } */
export let samples = $state([]);

/** ASV objects: { id, x, y, taxonomy, group, total_reads, prevalence, ... } */
export let asvs = $state([]);

/** Sparse count matrix: [{ sample_idx, asv_idx, count }, ...] */
export let counts = $state([]);

/** Network edges: [{ source, target, weight }, ...] */
export let network = $state([]);

/** Per-database taxonomy assignments: { silva: [...], unite: [...], ... } */
export let taxonomy = $state({});

// ── Selection state ─────────────────────────────────────────────────────────

export let selectedSample = $state(null);
export let selectedAsv = $state(null);

// ── UI state ────────────────────────────────────────────────────────────────

export let loading = $state(true);
export let error = $state(null);

// ── Data loading ────────────────────────────────────────────────────────────

async function fetchJson(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status}`);
  return res.json();
}

/**
 * Load all data files from /data/*.json.
 * Missing files are tolerated (logged, not thrown).
 */
export async function loadData() {
  loading = true;
  error = null;

  const files = [
    { key: 'samples', url: '/data/samples.json' },
    { key: 'asvs', url: '/data/asvs.json' },
    { key: 'counts', url: '/data/counts.json' },
    { key: 'network', url: '/data/network.json' },
    { key: 'taxonomy', url: '/data/taxonomy.json' },
  ];

  const results = {};

  await Promise.all(
    files.map(async ({ key, url }) => {
      try {
        results[key] = await fetchJson(url);
      } catch (e) {
        console.warn(`[microscape-viz] Could not load ${url}:`, e.message);
        results[key] = key === 'taxonomy' ? {} : [];
      }
    })
  );

  samples = results.samples;
  asvs = results.asvs;
  counts = results.counts;
  network = results.network;
  taxonomy = results.taxonomy;

  loading = false;
}

// ── Derived helpers ─────────────────────────────────────────────────────────

/** Build a Map from sample_idx -> [{asv_idx, count}, ...] */
export function countsBySample() {
  const map = new Map();
  for (const c of counts) {
    if (!map.has(c.sample_idx)) map.set(c.sample_idx, []);
    map.get(c.sample_idx).push({ asv_idx: c.asv_idx, count: c.count });
  }
  return map;
}

/** Build a Map from asv_idx -> [{sample_idx, count}, ...] */
export function countsByAsv() {
  const map = new Map();
  for (const c of counts) {
    if (!map.has(c.asv_idx)) map.set(c.asv_idx, []);
    map.get(c.asv_idx).push({ sample_idx: c.sample_idx, count: c.count });
  }
  return map;
}

/** Group colors for the four ASV groups */
export const GROUP_COLORS = {
  prokaryote: [0.3, 0.5, 1.0, 0.8],    // blue
  eukaryote: [1.0, 0.3, 0.3, 0.8],      // red
  chloroplast: [0.2, 0.85, 0.4, 0.8],   // green
  mitochondria: [0.2, 0.9, 0.9, 0.8],   // cyan
};

/** Group colors as hex for UI elements */
export const GROUP_HEX = {
  prokaryote: '#4d80ff',
  eukaryote: '#ff4d4d',
  chloroplast: '#33d966',
  mitochondria: '#33e6e6',
};
