<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import {
    store, countsBySample,
    GROUP_COLORS, GROUP_HEX,
  } from '../stores/data.svelte.js';

  // ── Sidebar controls ──────────────────────────────────────────────────────
  let minReads = $state(0);
  let taxonFilter = $state('');
  let groupFlags = $state({
    prokaryote: true,
    eukaryote: true,
    chloroplast: true,
    mitochondria: true,
  });
  let showOverlay = $state(true);

  // ── Canvas + scatterplot ──────────────────────────────────────────────────
  let canvasContainer = $state(null);
  let canvas = $state(null);
  let scatterplot = $state(null);
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });

  // ── Derived data ──────────────────────────────────────────────────────────

  /** Filtered store.samples by min read depth */
  let filteredSamples = $derived(
    store.samples.filter(s => (s.reads ?? 0) >= minReads)
  );

  /** Precompute store.counts-by-sample map */
  let cMap = $derived(countsBySample());

  /** Regex for taxonomy filtering (case-insensitive) */
  let taxonRe = $derived(() => {
    try { return taxonFilter ? new RegExp(taxonFilter, 'i') : null; }
    catch { return null; }
  });

  /** Build overlay points: per-sample taxa scatter */
  let overlayPoints = $derived.by(() => {
    if (!showOverlay || filteredSamples.length === 0 || store.asvs.length === 0) return [];

    const re = taxonRe();
    const pts = [];

    for (const sample of filteredSamples) {
      const sIdx = store.samples.indexOf(sample);
      const entries = cMap.get(sIdx) ?? [];
      const totalCount = entries.reduce((s, e) => s + e.count, 0) || 1;

      for (const { asv_idx, count } of entries) {
        const asv = store.asvs[asv_idx];
        if (!asv) continue;

        const group = asv.group ?? 'prokaryote';
        if (!groupFlags[group]) continue;
        if (re && !(re.test(asv.taxonomy ?? '') || re.test(asv.id ?? ''))) continue;

        const proportion = count / totalCount;
        const color = GROUP_COLORS[group] ?? GROUP_COLORS.prokaryote;

        pts.push({
          x: sample.x + (Math.random() - 0.5) * 0.3,
          y: sample.y + (Math.random() - 0.5) * 0.3,
          size: Math.max(1, proportion * 20),
          color,
          sampleIdx: sIdx,
          asvIdx: asv_idx,
        });
      }
    }
    return pts;
  });

  // ── Top taxa for selected sample ──────────────────────────────────────────
  let topTaxa = $derived.by(() => {
    if (store.selectedSample == null) return [];
    const entries = cMap.get(store.selectedSample) ?? [];
    const total = entries.reduce((s, e) => s + e.count, 0) || 1;
    return entries
      .map(e => ({
        asv: store.asvs[e.asv_idx],
        count: e.count,
        pct: ((e.count / total) * 100).toFixed(1),
      }))
      .filter(e => e.asv)
      .sort((a, b) => b.count - a.count)
      .slice(0, 20);
  });

  let selectedSampleObj = $derived(
    store.selectedSample != null ? store.samples[store.selectedSample] : null
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
        pointSize: 5,
        opacity: 0.8,
        lassoOnLongPress: true,
        backgroundColor: [0.02, 0.06, 0.1, 1],
      });

      sp.subscribe('pointover', (idx) => {
        const s = filteredSamples[idx];
        if (s) {
          tooltip = {
            show: true,
            x: 0, y: 0,
            text: `${s.id ?? 'Sample ' + idx} | ${(s.reads ?? 0).toLocaleString()} reads`,
          };
        }
      });

      sp.subscribe('pointout', () => {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      });

      sp.subscribe('select', ({ points: indices }) => {
        if (indices.length > 0) {
          const sIdx = store.samples.indexOf(filteredSamples[indices[0]]);
          store.selectedSample = sIdx >= 0 ? sIdx : null;
        }
      });

      scatterplot = sp;
    }
  });

  // ── Redraw base layer when filtered store.samples change ────────────────────────
  $effect(() => {
    if (!scatterplot || filteredSamples.length === 0) return;

    const xArr = filteredSamples.map(s => s.x);
    const yArr = filteredSamples.map(s => s.y);
    const sizes = filteredSamples.map(s => Math.max(2, Math.log2((s.reads ?? 1) + 1)));
    const colors = filteredSamples.map(() => [0.5, 0.55, 0.65, 0.5]);

    scatterplot.draw({
      x: xArr,
      y: yArr,
      size: sizes,
      color: colors,
    });
  });

  // ── Handle resize ─────────────────────────────────────────────────────────
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
    <h2 class="text-xs font-semibold uppercase tracking-wider text-slate-400">Controls</h2>

    <!-- Min reads slider -->
    <label class="block">
      <span class="text-xs text-slate-400">Min reads: {minReads.toLocaleString()}</span>
      <input
        type="range"
        min="0" max="100000" step="100"
        bind:value={minReads}
        class="mt-1 w-full accent-blue-500"
      />
    </label>

    <!-- Show overlay -->
    <label class="flex items-center gap-2 text-sm">
      <input type="checkbox" bind:checked={showOverlay} class="accent-blue-500" />
      Show taxa overlay
    </label>

    <!-- Taxonomy filter -->
    <label class="block">
      <span class="text-xs text-slate-400">Taxonomy filter (regex)</span>
      <input
        type="text"
        bind:value={taxonFilter}
        placeholder="e.g. Bacill"
        class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200 placeholder-slate-500 focus:border-blue-500 focus:outline-none"
      />
    </label>

    <!-- Group checkboxes -->
    <fieldset class="space-y-1">
      <legend class="text-xs text-slate-400">Groups</legend>
      {#each Object.keys(groupFlags) as group}
        <label class="flex items-center gap-2 text-sm capitalize">
          <input type="checkbox" bind:checked={groupFlags[group]} class="accent-blue-500" />
          <span class="inline-block h-2.5 w-2.5 rounded-full" style="background:{GROUP_HEX[group]}"></span>
          {group}
        </label>
      {/each}
    </fieldset>

    <!-- Legend -->
    <div class="mt-auto border-t border-slate-800 pt-3">
      <p class="text-xs text-slate-500">Base points sized by log(read depth).</p>
      <p class="text-xs text-slate-500">Overlay points sized by ASV proportion.</p>
      <p class="text-xs text-slate-500 mt-1">{filteredSamples.length} / {store.samples.length} samples shown</p>
    </div>
  </aside>

  <!-- Main content -->
  <div class="flex flex-1 flex-col overflow-hidden">
    <!-- Canvas area -->
    <div class="relative flex-1" bind:this={canvasContainer}>
      <canvas bind:this={canvas} class="absolute inset-0 h-full w-full"></canvas>

      {#if tooltip.show}
        <div class="pointer-events-none absolute left-4 top-4 rounded bg-slate-800/90 px-3 py-1.5 text-xs text-slate-200 shadow-lg">
          {tooltip.text}
        </div>
      {/if}
    </div>

    <!-- Selected sample detail -->
    {#if selectedSampleObj}
      <div class="border-t border-slate-800 bg-slate-900/80 p-4">
        <div class="mb-2 flex items-center justify-between">
          <h3 class="text-sm font-semibold text-slate-200">
            {selectedSampleObj.id ?? 'Sample'} &mdash; {(selectedSampleObj.reads ?? 0).toLocaleString()} reads
          </h3>
          <button
            class="text-xs text-slate-500 hover:text-slate-300"
            onclick={() => store.selectedSample = null}
          >Close</button>
        </div>

        {#if topTaxa.length > 0}
          <div class="max-h-48 overflow-y-auto">
            <table class="w-full text-xs">
              <thead class="sticky top-0 bg-slate-900 text-left text-slate-400">
                <tr>
                  <th class="py-1 pr-4">ASV</th>
                  <th class="py-1 pr-4">Taxonomy</th>
                  <th class="py-1 pr-4">Group</th>
                  <th class="py-1 pr-4 text-right">Reads</th>
                  <th class="py-1 text-right">%</th>
                </tr>
              </thead>
              <tbody class="text-slate-300">
                {#each topTaxa as row}
                  <tr class="border-t border-slate-800/50 hover:bg-slate-800/30">
                    <td class="py-1 pr-4 font-mono">{row.asv.id ?? ''}</td>
                    <td class="py-1 pr-4 max-w-xs truncate">{row.asv.taxonomy ?? ''}</td>
                    <td class="py-1 pr-4">
                      <span class="inline-block h-2 w-2 rounded-full mr-1" style="background:{GROUP_HEX[row.asv.group] ?? GROUP_HEX.prokaryote}"></span>
                      {row.asv.group ?? ''}
                    </td>
                    <td class="py-1 pr-4 text-right font-mono">{row.count.toLocaleString()}</td>
                    <td class="py-1 text-right font-mono">{row.pct}</td>
                  </tr>
                {/each}
              </tbody>
            </table>
          </div>
        {:else}
          <p class="text-xs text-slate-500">No taxa for this sample.</p>
        {/if}
      </div>
    {/if}
  </div>
</div>
