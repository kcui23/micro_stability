<script>
  import { createEventDispatcher } from 'svelte';
  import { scale } from 'svelte/transition';
  import InteractiveValcano from './InteractiveValcano.svelte';
  import { selectedPoints } from '../store.js';

  export let selectedMethod;
  
  const dispatch = createEventDispatcher();

  function closePage() {
    dispatch('close');
  }

  let searchTerm = '';

  function removeASV(asvName) {
    selectedPoints.update(points => points.filter(p => p.name !== asvName));
  }

  function clearAllASVs() {
    selectedPoints.set([]);
  }

  function downloadASVs() {
    const asvList = $selectedPoints.map(p => p.name).join('\n');
    const blob = new Blob([asvList], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'selected_asvs.tsv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  $: filteredASVs = $selectedPoints.filter(p => 
    p.name.toLowerCase().includes(searchTerm.toLowerCase())
  );
</script>

<div 
  class="overlay"
  transition:scale={{
    duration: 300,
    start: 1,
    opacity: 0
  }}
  style="backdrop-filter: blur(5px);"
>
  <div 
    class="asv-page"
    transition:scale={{
      duration: 300,
      delay: 150,
      start: 0.98
    }}
  >
    <div 
      class="close-button" 
      on:click={closePage}
      on:keydown={e => e.key === 'Enter' && closePage()}
      role="button"
      tabindex="0"
      aria-label="Close page">×</div>
    <div class="page-content">
      <div class="left-section">
        <div class="graph-placeholder">
            <InteractiveValcano {selectedMethod} specific_interact=true />
        </div>
      </div>
      <div class="right-section">
        <div class="selected-asvs">
          <div class="header">
            <h3>Selected ASVs ({$selectedPoints.length})</h3>
            {#if $selectedPoints.length > 0}
              <button class="clear-all" on:click={clearAllASVs}>Clear All</button>
            {/if}
          </div>
          <input
            type="text"
            placeholder="Search ASVs..."
            bind:value={searchTerm}
            class="search-input"
          />
          <div class="asv-list">
            {#each filteredASVs as asv}
              <div class="asv-item">
                <span>{asv.name}</span>
                <button class="remove-btn" on:click={() => removeASV(asv.name)}>×</button>
              </div>
            {/each}
          </div>
        </div>
        <button 
          class="download-button" 
          on:click={downloadASVs}
          disabled={$selectedPoints.length === 0}
        >
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round">
            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" stroke-width="3"></path>
            <polyline points="7 10 12 15 17 10" stroke-width="2" color="#a3a3a3"></polyline>
            <line x1="12" y1="15" x2="12" y2="3" stroke-width="2" color="#a3a3a3"></line>
          </svg>
          Download selected data
        </button>
      </div>
    </div>
  </div>
</div>

<style>
  .overlay {
    position: fixed;
    top: 0;
    left: 0;
    width: 100vw;
    height: 100vh;
    background: rgba(0, 0, 0, 0.15);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 1000;
    -webkit-backdrop-filter: blur(5px);
    backdrop-filter: blur(5px);
  }

  .asv-page {
    width: calc(100vw - 200px);
    height: calc(100vh - 200px);
    background: white;
    padding: 40px;
    border-radius: 15px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
    position: relative;
  }

  .close-button {
    position: absolute;
    top: -20px;
    right: -20px;
    width: 40px;
    height: 40px;
    background: #ececec;
    color: black;
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    font-size: 24px;
    z-index: 1001;
    transition: transform 0.2s ease, background-color 0.2s ease;
  }

  .close-button:hover {
    transform: scale(1.1);
    background: #bdbdbd;
  }

  .close-button:active {
    transform: scale(0.9);
    background: #bdbdbd;
  }

  .page-content {
    display: flex;
    gap: 20px;
    height: 100%;
    width: 100%;
  }

  .left-section {
    flex: 2;
    display: flex;
    flex-direction: column;
  }

  .right-section {
    flex: 1;
    display: flex;
    flex-direction: column;
  }

  .graph-placeholder {
    flex: 1;
    background: #f5f5f5;
    border-radius: 10px;
    aspect-ratio: 1;
  }

  .content-placeholder {
    flex: 1;
    background: #f5f5f5;
    border-radius: 10px;
    margin-bottom: 20px;
  }

  .download-button {
    font-size: 14px;
    letter-spacing: 0.1px;
    gap: 5px;
  }
  .download-button svg {
    width: 16px;
    height: 16px;
  }

  .selected-asvs {
    flex: 1;
    background: #f5f5f5;
    border-radius: 10px;
    padding: 20px;
    display: flex;
    flex-direction: column;
    gap: 15px;
    margin-bottom: 20px;
    min-height: 0;
    overflow: hidden;
  }

  .header {
    display: flex;
    justify-content: space-between;
    align-items: center;
  }

  .header h3 {
    margin: 0;
  }

  .clear-all {
    background: none;
    border: none;
    color: #ff4444;
    cursor: pointer;
    font-size: 14px;
  }

  .search-input {
    padding: 8px;
    border: 1px solid #ddd;
    border-radius: 5px;
    width: 100%;
  }

  .asv-list {
    flex: 1;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    gap: 8px;
    padding-right: 5px;
    &::-webkit-scrollbar {
      width: 8px;
    }
    &::-webkit-scrollbar-track {
      background: #f1f1f1;
      border-radius: 4px;
    }
    &::-webkit-scrollbar-thumb {
      background: #888;
      border-radius: 4px;
    }
    &::-webkit-scrollbar-thumb:hover {
      background: #555;
    }
  }

  .asv-item {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 8px;
    background: white;
    border-radius: 5px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
  }

  .remove-btn {
    background: none;
    border: none;
    color: #666;
    cursor: pointer;
    font-size: 18px;
    padding: 0 5px;
  }

  .remove-btn:hover {
    color: #ff4444;
  }

  .download-button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
  }
</style> 