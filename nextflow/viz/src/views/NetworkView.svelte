<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import {
    store, countsByAsv,
    GROUP_COLORS, GROUP_HEX,
  } from '../stores/data.svelte.js';

  // ── Controls ──────────────────────────────────────────────────────────────
  let minPrevalence = $state(0);
  let corrThreshold = $state(0.3);
  let taxonFilter = $state('');
  let showEdges = $state(true);

  // ── Canvas ────────────────────────────────────────────────────────────────
  let canvasContainer = $state(null);
  let canvas = $state(null);
  let scatterplot = $state(null);
  let svgOverlay = $state(null);
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });

  // ── Derived ───────────────────────────────────────────────────────────────

  let taxonRe = $derived(() => {
    try { return taxonFilter ? new RegExp(taxonFilter, 'i') : null; }
    catch { return null; }
  });

  /** Filtered ASVs */
  let filteredAsvs = $derived.by(() => {
    const re = taxonRe();
    return store.asvs.filter(a => {
      if ((a.prevalence ?? 0) < minPrevalence) return false;
      if (re && !(re.test(a.taxonomy ?? '') || re.test(a.id ?? ''))) return false;
      return true;
    });
  });

  /** Index mapping: original ASV index -> filtered index */
  let idxMap = $derived.by(() => {
    const m = new Map();
    filteredAsvs.forEach((a, fi) => {
      const oi = store.asvs.indexOf(a);
      m.set(oi, fi);
    });
    return m;
  });

  /** Filtered edges */
  let filteredEdges = $derived.by(() => {
    if (!showEdges) return [];
    return (store.network?.edges ?? store.network ?? []).filter(e => {
      if (Math.abs(e.weight ?? 0) < corrThreshold) return false;
      return idxMap.has(e.source) && idxMap.has(e.target);
    });
  });

  // ── Selected ASV detail ───────────────────────────────────────────────────
  let selectedAsvObj = $derived(
    store.selectedAsv != null ? store.asvs[store.selectedAsv] : null
  );

  // ── Scatterplot lifecycle ─────────────────────────────────────────────────
  onMount(() => {
    return () => {
      if (scatterplot) scatterplot.destroy();
    };
  });

  $effect(() => {
    if (canvas && !scatterplot) {
      const rect = canvasContainer.getBoundingClientRect();
      const sp = createScatterplot({
        canvas,
        width: rect.width,
        height: rect.height,
        pointSize: 4,
        opacity: 0.85,
        lassoOnLongPress: true,
        backgroundColor: [0.02, 0.06, 0.1, 1],
      });

      sp.subscribe('pointover', (idx) => {
        const a = filteredAsvs[idx];
        if (a) {
          tooltip = {
            show: true,
            x: 0, y: 0,
            text: `${a.id ?? 'ASV'} | ${a.taxonomy ?? ''} | ${(a.total_reads ?? 0).toLocaleString()} reads`,
          };
        }
      });

      sp.subscribe('pointout', () => {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      });

      sp.subscribe('select', ({ points: indices }) => {
        if (indices.length > 0) {
          const oi = store.asvs.indexOf(filteredAsvs[indices[0]]);
          store.selectedAsv = oi >= 0 ? oi : null;
        }
      });

      scatterplot = sp;
    }
  });

  // ── Redraw points ─────────────────────────────────────────────────────────
  $effect(() => {
    if (!scatterplot || filteredAsvs.length === 0) return;

    const xArr = filteredAsvs.map(a => a.x ?? 0);
    const yArr = filteredAsvs.map(a => a.y ?? 0);
    const sizes = filteredAsvs.map(a => Math.max(2, Math.log2((a.total_reads ?? 1) + 1) * 0.8));
    const colors = filteredAsvs.map(a => GROUP_COLORS[a.group ?? 'prokaryote'] ?? GROUP_COLORS.prokaryote);

    scatterplot.draw({ x: xArr, y: yArr, size: sizes, color: colors });
  });

  // ── Draw edges on SVG overlay ─────────────────────────────────────────────
  let edgeLines = $derived.by(() => {
    if (!showEdges || filteredEdges.length === 0 || filteredAsvs.length === 0) return [];

    // We need the scatterplot view transform to map data coords to screen coords.
    // For a simple fallback, normalize x/y to canvas dimensions.
    // In practice the scatterplot's camera handles this, so we approximate.
    const xs = filteredAsvs.map(a => a.x ?? 0);
    const ys = filteredAsvs.map(a => a.y ?? 0);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const rangeX = maxX - minX || 1;
    const rangeY = maxY - minY || 1;

    const rect = canvasContainer?.getBoundingClientRect();
    const w = rect?.width ?? 800;
    const h = rect?.height ?? 600;
    const pad = 40;

    function toScreen(ax, ay) {
      return {
        sx: pad + ((ax - minX) / rangeX) * (w - 2 * pad),
        sy: pad + ((ay - minY) / rangeY) * (h - 2 * pad),
      };
    }

    return filteredEdges.slice(0, 500).map(e => {
      const srcAsv = filteredAsvs[idxMap.get(e.source)];
      const tgtAsv = filteredAsvs[idxMap.get(e.target)];
      if (!srcAsv || !tgtAsv) return null;

      const s = toScreen(srcAsv.x ?? 0, srcAsv.y ?? 0);
      const t = toScreen(tgtAsv.x ?? 0, tgtAsv.y ?? 0);

      return {
        x1: s.sx, y1: s.sy,
        x2: t.sx, y2: t.sy,
        opacity: Math.min(1, Math.abs(e.weight ?? 0)),
        color: (e.weight ?? 0) > 0 ? 'rgba(100,160,255,0.25)' : 'rgba(255,100,100,0.25)',
      };
    }).filter(Boolean);
  });

  function handleResize() {
    if (scatterplot && canvasContainer) {
      const rect = canvasContainer.getBoundingClientRect();
      scatterplot.set({ width: rect.width, height: rect.height });
    }
  }
</script>

<svelte:window onresize={handleResize} />

<div class="flex h-full">
  <!-- Sidebar -->
  <aside class="flex w-64 flex-col gap-4 overflow-y-auto border-r border-slate-800 bg-slate-900/60 p-4">
    <h2 class="text-xs font-semibold uppercase tracking-wider text-slate-400">Network Controls</h2>

    <label class="block">
      <span class="text-xs text-slate-400">Min prevalence: {minPrevalence}</span>
      <input
        type="range" min="0" max="100" step="1"
        bind:value={minPrevalence}
        class="mt-1 w-full accent-blue-500"
      />
    </label>

    <label class="block">
      <span class="text-xs text-slate-400">Correlation threshold: {corrThreshold.toFixed(2)}</span>
      <input
        type="range" min="0" max="1" step="0.01"
        bind:value={corrThreshold}
        class="mt-1 w-full accent-blue-500"
      />
    </label>

    <label class="block">
      <span class="text-xs text-slate-400">Taxonomy filter (regex)</span>
      <input
        type="text"
        bind:value={taxonFilter}
        placeholder="e.g. Bacill"
        class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200 placeholder-slate-500 focus:border-blue-500 focus:outline-none"
      />
    </label>

    <label class="flex items-center gap-2 text-sm">
      <input type="checkbox" bind:checked={showEdges} class="accent-blue-500" />
      Show edges (max 500)
    </label>

    <!-- Legend -->
    <div class="mt-auto space-y-1 border-t border-slate-800 pt-3">
      {#each Object.entries(GROUP_HEX) as [group, hex]}
        <div class="flex items-center gap-2 text-xs capitalize text-slate-400">
          <span class="inline-block h-2.5 w-2.5 rounded-full" style="background:{hex}"></span>
          {group}
        </div>
      {/each}
      <p class="text-xs text-slate-500 mt-2">{filteredAsvs.length} / {store.asvs.length} ASVs shown</p>
      <p class="text-xs text-slate-500">{filteredEdges.length} edges above threshold</p>
    </div>
  </aside>

  <!-- Main area -->
  <div class="flex flex-1 flex-col overflow-hidden">
    <div class="relative flex-1" bind:this={canvasContainer}>
      <canvas bind:this={canvas} class="absolute inset-0 h-full w-full"></canvas>

      <!-- SVG edge overlay -->
      {#if edgeLines.length > 0}
        <svg class="pointer-events-none absolute inset-0 h-full w-full" bind:this={svgOverlay}>
          {#each edgeLines as line}
            <line
              x1={line.x1} y1={line.y1}
              x2={line.x2} y2={line.y2}
              stroke={line.color}
              stroke-width="1"
            />
          {/each}
        </svg>
      {/if}

      {#if tooltip.show}
        <div class="pointer-events-none absolute left-4 top-4 rounded bg-slate-800/90 px-3 py-1.5 text-xs text-slate-200 shadow-lg">
          {tooltip.text}
        </div>
      {/if}
    </div>

    <!-- Selected ASV detail -->
    {#if selectedAsvObj}
      <div class="border-t border-slate-800 bg-slate-900/80 p-4">
        <div class="flex items-center justify-between">
          <h3 class="text-sm font-semibold text-slate-200">
            {selectedAsvObj.id ?? 'ASV'} &mdash; {selectedAsvObj.group ?? ''}
          </h3>
          <button
            class="text-xs text-slate-500 hover:text-slate-300"
            onclick={() => store.selectedAsv = null}
          >Close</button>
        </div>
        <p class="mt-1 text-xs text-slate-400">{selectedAsvObj.taxonomy ?? 'No taxonomy'}</p>
        <p class="text-xs text-slate-500">
          {(selectedAsvObj.total_reads ?? 0).toLocaleString()} total reads |
          Prevalence: {selectedAsvObj.prevalence ?? 0}
        </p>
      </div>
    {/if}
  </div>
</div>
