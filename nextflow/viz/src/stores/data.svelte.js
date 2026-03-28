/**
 * Reactive data stores for microscape visualization.
 *
 * Svelte 5 module-level $state can't be exported AND reassigned.
 * Solution: wrap all state in a single exported object whose
 * properties are mutated (not reassigned).
 */

export const store = $state({
  // Core data
  samples: [],
  asvs: [],
  counts: { data: [], samples: [], asvs: [] },
  network: { edges: [] },
  taxonomy: {},

  // Selection
  selectedSample: null,
  selectedAsv: null,

  // UI
  loading: true,
  error: null,
});

/** Group colors as RGBA for regl-scatterplot */
export const GROUP_COLORS = {
  prokaryote: [0.3, 0.5, 1.0, 0.8],
  eukaryote: [1.0, 0.3, 0.3, 0.8],
  chloroplast: [0.2, 0.85, 0.4, 0.8],
  mitochondria: [0.2, 0.9, 0.9, 0.8],
  unknown: [0.6, 0.6, 0.6, 0.5],
};

/** Group colors as hex for UI */
export const GROUP_HEX = {
  prokaryote: '#4d80ff',
  eukaryote: '#ff4d4d',
  chloroplast: '#33d966',
  mitochondria: '#33e6e6',
  unknown: '#999999',
};

// ── Data loading ────────────────────────────────────────────────────────────

async function fetchJson(url) {
  // Try .gz first
  try {
    const gzRes = await fetch(url + '.gz');
    if (gzRes.ok) {
      const buf = await gzRes.arrayBuffer();
      const bytes = new Uint8Array(buf);
      if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
        const ds = new DecompressionStream('gzip');
        const reader = new Blob([buf]).stream().pipeThrough(ds).getReader();
        const chunks = [];
        while (true) {
          const { done, value } = await reader.read();
          if (done) break;
          chunks.push(value);
        }
        const combined = new Uint8Array(chunks.reduce((a, c) => a + c.length, 0));
        let offset = 0;
        for (const c of chunks) { combined.set(c, offset); offset += c.length; }
        return JSON.parse(new TextDecoder().decode(combined));
      }
      return JSON.parse(new TextDecoder().decode(buf));
    }
  } catch (_) { /* fall through */ }

  const res = await fetch(url);
  if (!res.ok) throw new Error(`${res.status} ${url}`);
  return res.json();
}

export async function loadData() {
  store.loading = true;
  store.error = null;

  try {
    const [samples, asvs, counts, network, taxonomy] = await Promise.all([
      fetchJson('/data/samples.json').catch(() => []),
      fetchJson('/data/asvs.json').catch(() => []),
      fetchJson('/data/counts.json').catch(() => ({ data: [], samples: [], asvs: [] })),
      fetchJson('/data/network.json').catch(() => ({ edges: [] })),
      fetchJson('/data/taxonomy.json').catch(() => ({})),
    ]);

    store.samples = samples;
    store.asvs = asvs;
    store.counts = counts;
    store.network = network;
    store.taxonomy = taxonomy;
  } catch (e) {
    store.error = e.message;
    console.error('[microscape] Data load failed:', e);
  }

  store.loading = false;
}

// ── Helpers ─────────────────────────────────────────────────────────────────

export function countsBySample() {
  const map = new Map();
  const data = store.counts?.data;
  if (!data) return map;
  for (const row of data) {
    const [si, ai, count, prop] = row;
    if (!map.has(si)) map.set(si, []);
    map.get(si).push({ asv_idx: ai, count, proportion: prop });
  }
  return map;
}

export function countsByAsv() {
  const map = new Map();
  const data = store.counts?.data;
  if (!data) return map;
  for (const row of data) {
    const [si, ai, count, prop] = row;
    if (!map.has(ai)) map.set(ai, []);
    map.get(ai).push({ sample_idx: si, count, proportion: prop });
  }
  return map;
}
