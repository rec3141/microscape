/**
 * Reactive data stores for the microscape amplicon pipeline visualization.
 * Uses Svelte 5 runes ($state) with setter functions for cross-module access.
 */

// ── Core data ───────────────────────────────────────────────────────────────

export let samples = $state([]);
export let asvs = $state([]);
export let counts = $state([]);
export let network = $state([]);
export let taxonomy = $state({});

// ── Selection state ─────────────────────────────────────────────────────────

export let selectedSample = $state(null);
export let selectedAsv = $state(null);

/** Setter functions for cross-module mutation */
export function setSelectedSample(v) { selectedSample = v; }
export function setSelectedAsv(v) { selectedAsv = v; }

// ── UI state ────────────────────────────────────────────────────────────────

export let loading = $state(true);
export let error = $state(null);

// ── Data loading ────────────────────────────────────────────────────────────

async function fetchJson(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status}`);
  return res.json();
}

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
        // Try gzipped first
        let gzUrl = url + '.gz';
        let res = await fetch(gzUrl);
        if (res.ok) {
          // Check if it's actually gzipped
          const buf = await res.arrayBuffer();
          const bytes = new Uint8Array(buf);
          if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
            // Decompress
            const ds = new DecompressionStream('gzip');
            const reader = new Blob([buf]).stream().pipeThrough(ds).getReader();
            const chunks = [];
            while (true) {
              const { done, value } = await reader.read();
              if (done) break;
              chunks.push(value);
            }
            const text = new TextDecoder().decode(
              new Uint8Array(chunks.reduce((acc, c) => [...acc, ...c], []))
            );
            results[key] = JSON.parse(text);
          } else {
            // Not actually gzipped, parse as text
            results[key] = JSON.parse(new TextDecoder().decode(buf));
          }
        } else {
          // Fall back to plain JSON
          results[key] = await fetchJson(url);
        }
      } catch (e) {
        console.warn(`[microscape-viz] Could not load ${url}:`, e.message);
        results[key] = key === 'taxonomy' ? {} : [];
      }
    })
  );

  samples = results.samples || [];
  asvs = results.asvs || [];
  counts = results.counts || [];
  network = results.network || [];
  taxonomy = results.taxonomy || {};

  loading = false;
}

// ── Derived helpers ─────────────────────────────────────────────────────────

export function countsBySample() {
  const map = new Map();
  if (!counts || !counts.data) return map;
  for (const [si, ai, count, prop] of counts.data) {
    if (!map.has(si)) map.set(si, []);
    map.get(si).push({ asv_idx: ai, count, proportion: prop });
  }
  return map;
}

export function countsByAsv() {
  const map = new Map();
  if (!counts || !counts.data) return map;
  for (const [si, ai, count, prop] of counts.data) {
    if (!map.has(ai)) map.set(ai, []);
    map.get(ai).push({ sample_idx: si, count, proportion: prop });
  }
  return map;
}

/** Group colors as RGBA arrays for regl-scatterplot */
export const GROUP_COLORS = {
  prokaryote: [0.3, 0.5, 1.0, 0.8],
  eukaryote: [1.0, 0.3, 0.3, 0.8],
  chloroplast: [0.2, 0.85, 0.4, 0.8],
  mitochondria: [0.2, 0.9, 0.9, 0.8],
  unknown: [0.6, 0.6, 0.6, 0.5],
};

/** Group colors as hex for UI elements */
export const GROUP_HEX = {
  prokaryote: '#4d80ff',
  eukaryote: '#ff4d4d',
  chloroplast: '#33d966',
  mitochondria: '#33e6e6',
  unknown: '#999999',
};
