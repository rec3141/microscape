<script>
  import { onMount } from 'svelte';
  import NavBar from './components/NavBar.svelte';
  import SampleView from './views/SampleView.svelte';
  import NetworkView from './views/NetworkView.svelte';
  import TablesView from './views/TablesView.svelte';
  import { store, loadData } from './stores/data.svelte.js';

  let activeTab = $state('samples');

  function updateTab() {
    const hash = window.location.hash.replace('#', '') || 'samples';
    if (['samples', 'network', 'tables'].includes(hash)) {
      activeTab = hash;
    }
  }

  onMount(() => {
    updateTab();
    window.addEventListener('hashchange', updateTab);
    loadData();
    return () => window.removeEventListener('hashchange', updateTab);
  });
</script>

<div class="flex h-screen flex-col">
  <NavBar {activeTab} />

  <main class="flex-1 overflow-hidden">
    {#if store.loading}
      <div class="flex h-full items-center justify-center">
        <div class="text-center">
          <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
          <p class="text-sm text-slate-400">Loading data...</p>
        </div>
      </div>
    {:else if store.error}
      <div class="flex h-full items-center justify-center">
        <div class="rounded-lg border border-red-800 bg-red-950/50 p-6 text-center">
          <p class="text-red-400">{store.error}</p>
          <button
            class="mt-3 rounded bg-red-800 px-4 py-1.5 text-sm text-red-100 hover:bg-red-700"
            onclick={() => loadData()}
          >Retry</button>
        </div>
      </div>
    {:else}
      {#if activeTab === 'samples'}
        <SampleView />
      {:else if activeTab === 'network'}
        <NetworkView />
      {:else if activeTab === 'tables'}
        <TablesView />
      {/if}
    {/if}
  </main>
</div>
