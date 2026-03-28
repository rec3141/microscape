<script>
  import { store, GROUP_HEX } from '../stores/data.svelte.js';

  // ── Tab state ─────────────────────────────────────────────────────────────
  let activeTable = $state('asvs'); // 'asvs' | 'samples'
  let search = $state('');
  let sortCol = $state(null);
  let sortAsc = $state(true);

  // ── ASV table ─────────────────────────────────────────────────────────────
  const asvCols = [
    { key: 'id', label: 'ASV ID' },
    { key: 'taxonomy', label: 'Taxonomy' },
    { key: 'group', label: 'Group' },
    { key: 'total_reads', label: 'Total Reads', numeric: true },
    { key: 'prevalence', label: 'Prevalence', numeric: true },
  ];

  const sampleCols = [
    { key: 'id', label: 'Sample ID' },
    { key: 'reads', label: 'Reads', numeric: true },
    { key: 'richness', label: 'Richness', numeric: true },
  ];

  let cols = $derived(activeTable === 'asvs' ? asvCols : sampleCols);
  let rawRows = $derived(activeTable === 'asvs' ? store.asvs : store.samples);

  let filteredRows = $derived.by(() => {
    let rows = rawRows;

    // Search filter
    if (search) {
      const q = search.toLowerCase();
      rows = rows.filter(r =>
        Object.values(r).some(v =>
          v != null && String(v).toLowerCase().includes(q)
        )
      );
    }

    // Sort
    if (sortCol) {
      const col = cols.find(c => c.key === sortCol);
      const numeric = col?.numeric ?? false;
      rows = [...rows].sort((a, b) => {
        let va = a[sortCol] ?? '';
        let vb = b[sortCol] ?? '';
        if (numeric) {
          va = Number(va) || 0;
          vb = Number(vb) || 0;
          return sortAsc ? va - vb : vb - va;
        }
        va = String(va).toLowerCase();
        vb = String(vb).toLowerCase();
        return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
    }

    return rows;
  });

  // ── Pagination ────────────────────────────────────────────────────────────
  let page = $state(0);
  const perPage = 100;

  let pageRows = $derived(
    filteredRows.slice(page * perPage, (page + 1) * perPage)
  );

  let totalPages = $derived(Math.ceil(filteredRows.length / perPage));

  // Reset page when search changes
  $effect(() => { search; page = 0; });

  // ── Sort handler ──────────────────────────────────────────────────────────
  function toggleSort(key) {
    if (sortCol === key) {
      sortAsc = !sortAsc;
    } else {
      sortCol = key;
      sortAsc = true;
    }
  }

  // ── CSV export ────────────────────────────────────────────────────────────
  function exportCsv() {
    const header = cols.map(c => c.label).join(',');
    const body = filteredRows.map(r =>
      cols.map(c => {
        const v = r[c.key] ?? '';
        // Quote strings that contain commas
        return String(v).includes(',') ? `"${v}"` : v;
      }).join(',')
    ).join('\n');

    const blob = new Blob([header + '\n' + body], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${activeTable}.csv`;
    a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="flex h-full flex-col overflow-hidden p-4">
  <!-- Header bar -->
  <div class="mb-4 flex flex-wrap items-center gap-3">
    <!-- Table selector -->
    <div class="flex rounded-lg border border-slate-700 bg-slate-800">
      <button
        class="px-4 py-1.5 text-sm font-medium transition-colors rounded-l-lg
          {activeTable === 'asvs' ? 'bg-blue-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
        onclick={() => { activeTable = 'asvs'; sortCol = null; page = 0; }}
      >ASVs</button>
      <button
        class="px-4 py-1.5 text-sm font-medium transition-colors rounded-r-lg
          {activeTable === 'samples' ? 'bg-blue-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
        onclick={() => { activeTable = 'samples'; sortCol = null; page = 0; }}
      >Samples</button>
    </div>

    <!-- Search -->
    <input
      type="text"
      bind:value={search}
      placeholder="Search..."
      class="rounded border border-slate-700 bg-slate-800 px-3 py-1.5 text-sm text-slate-200 placeholder-slate-500 focus:border-blue-500 focus:outline-none w-64"
    />

    <span class="text-xs text-slate-500">
      {filteredRows.length} rows
    </span>

    <button
      class="ml-auto rounded border border-slate-700 bg-slate-800 px-3 py-1.5 text-xs font-medium text-slate-300 hover:bg-slate-700 hover:text-slate-100 transition-colors"
      onclick={exportCsv}
    >
      Export CSV
    </button>
  </div>

  <!-- Table -->
  <div class="flex-1 overflow-auto rounded-lg border border-slate-800">
    <table class="w-full text-sm">
      <thead class="sticky top-0 z-10 bg-slate-900 text-left">
        <tr>
          {#each cols as col}
            <th
              class="cursor-pointer select-none border-b border-slate-800 px-4 py-2.5 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200 transition-colors"
              onclick={() => toggleSort(col.key)}
            >
              <span class="inline-flex items-center gap-1">
                {col.label}
                {#if sortCol === col.key}
                  <span class="text-blue-400">{sortAsc ? '▲' : '▼'}</span>
                {/if}
              </span>
            </th>
          {/each}
        </tr>
      </thead>
      <tbody>
        {#each pageRows as row, i}
          <tr class="border-t border-slate-800/50 hover:bg-slate-800/30 transition-colors">
            {#each cols as col}
              <td class="px-4 py-2 text-slate-300 {col.numeric ? 'text-right font-mono' : ''}">
                {#if col.key === 'group' && row[col.key]}
                  <span class="inline-flex items-center gap-1.5 capitalize">
                    <span class="inline-block h-2 w-2 rounded-full" style="background:{GROUP_HEX[row[col.key]] ?? '#888'}"></span>
                    {row[col.key]}
                  </span>
                {:else if col.numeric}
                  {(row[col.key] ?? 0).toLocaleString()}
                {:else}
                  <span class="max-w-md truncate block">{row[col.key] ?? ''}</span>
                {/if}
              </td>
            {/each}
          </tr>
        {/each}

        {#if pageRows.length === 0}
          <tr>
            <td colspan={cols.length} class="px-4 py-8 text-center text-slate-500">
              No data available.
            </td>
          </tr>
        {/if}
      </tbody>
    </table>
  </div>

  <!-- Pagination -->
  {#if totalPages > 1}
    <div class="mt-3 flex items-center justify-center gap-2">
      <button
        class="rounded border border-slate-700 bg-slate-800 px-3 py-1 text-xs text-slate-400 hover:text-slate-200 disabled:opacity-30"
        disabled={page === 0}
        onclick={() => page--}
      >Prev</button>
      <span class="text-xs text-slate-500">
        Page {page + 1} / {totalPages}
      </span>
      <button
        class="rounded border border-slate-700 bg-slate-800 px-3 py-1 text-xs text-slate-400 hover:text-slate-200 disabled:opacity-30"
        disabled={page >= totalPages - 1}
        onclick={() => page++}
      >Next</button>
    </div>
  {/if}
</div>
