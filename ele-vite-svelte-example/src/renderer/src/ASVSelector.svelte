<script>
    import { onMount, tick } from 'svelte';
  
    let asvList = [];
    let selectedASVs = [];
    let searchTerm = '';
    let combinedResults = null;
    let significantMethods = [];
    let dropdownOpen = false;
  
    onMount(async () => {
      // Fetch the list of ASVs and combined results
      const response = await fetch('http://localhost:8000/get_asv_list');
      asvList = await response.json();
  
      const resultsResponse = await fetch('http://localhost:8000/get_combined_results');
      combinedResults = await resultsResponse.json();
    });
  
    $: filteredASVs = asvList.filter(asv => asv.toLowerCase().includes(searchTerm.toLowerCase()));
  
    function toggleASV(asv) {
      const index = selectedASVs.indexOf(asv);
      if (index === -1) {
        selectedASVs = [...selectedASVs, asv];
      } else {
        selectedASVs = selectedASVs.filter(a => a !== asv);
      }
      updateSignificantMethods();
    }

    function clickOutside(node) {
        const handleClick = event => {
            if (node && !node.contains(event.target) && !event.defaultPrevented) {
            node.dispatchEvent(
                new CustomEvent('clickoutside', node)
            )
            }
        }

        document.addEventListener('click', handleClick, true);

        return {
            destroy() {
            document.removeEventListener('click', handleClick, true);
            }
        }
    }

    function removeASV(asv) {
      selectedASVs = selectedASVs.filter(a => a !== asv);
      updateSignificantMethods();
    }
  
    function updateSignificantMethods() {
      if (!combinedResults) return;
  
      const methods = ['deseq2', 'aldex2', 'edger', 'maaslin2'];
      significantMethods = methods.filter(method => 
        selectedASVs.every(asv => 
          combinedResults.find(result => result.asv_name === asv)[`${method}_significant`]
        )
      );
    }
  
    function toggleDropdown() {
      dropdownOpen = !dropdownOpen;
    }
  
    function closeDropdown() {
      dropdownOpen = false;
    }
  </script>
  
  <div class="asv-selector">
    <div class="asv-list">
      <div class="dropdown" use:clickOutside on:clickoutside={closeDropdown}>
        <div class="input-wrapper">
          <input 
            type="text" 
            bind:value={searchTerm} 
            placeholder="Search ASVs..." 
            on:focus={() => dropdownOpen = true}
          >
          <button class="dropdown-toggle" on:click={toggleDropdown}>
            {dropdownOpen ? '▲' : '▼'}
          </button>
        </div>
        {#if dropdownOpen}
          <div class="dropdown-content">
            {#each filteredASVs as asv}
              <label class="asv-item">
                <input type="checkbox" checked={selectedASVs.includes(asv)} on:change={() => toggleASV(asv)}>
                <span>{asv}</span>
              </label>
            {/each}
          </div>
        {/if}
      </div>
    </div>
    
    <div class="selected-asvs">
        <h3>Selected ASVs</h3>
        <ul>
          {#each selectedASVs as asv}
            <li class="selected-asv-item">
              <span 
                role="button" 
                tabindex="0" 
                on:click={() => removeASV(asv)}
                on:keypress={(e) => e.key === 'Enter' && removeASV(asv)}
                aria-label="Remove {asv}"
              >
                {asv}
              </span>
              <button on:click={() => removeASV(asv)} aria-label="Remove {asv}">×</button>
            </li>
          {/each}
        </ul>
    </div>
  </div>
  
  <div class="significant-methods">
    <h3>Methods with all selected ASVs significant:</h3>
    <ul>
      {#each significantMethods as method}
        <li>{method}</li>
      {/each}
    </ul>
  </div>
  
  <style>
    .asv-selector {
      display: flex;
      justify-content: space-between;
      top: 10px;
    }
    .asv-list, .selected-asvs {
      width: 45%;
    }
    .dropdown {
      position: relative;
    }
    .input-wrapper {
      display: flex;
      align-items: center;
    }
    input[type="text"] {
      width: 100%;
      padding: 5px;
    }
    .dropdown-toggle {
      background: none;
      border: none;
      cursor: pointer;
      padding: 5px;
    }
    .dropdown-content {
      position: absolute;
      top: 100%;
      left: 0;
      right: 0;
      background-color: white;
      border: 1px solid #ccc;
      max-height: 200px;
      overflow-y: auto;
      z-index: 1;
    }
    .asv-item, .selected-asv-item {
      display: flex;
      align-items: center;
      /* padding: 5px; */
      cursor: pointer;
    }
    .asv-item:hover, .selected-asv-item:hover {
      background-color: #f0f0f0;
    }
    .asv-item input[type="checkbox"] {
      margin-right: 5px;
    }
    .selected-asvs ul {
      list-style-type: none;
      padding: 0;
    }
    .selected-asv-item button {
      margin-left: auto;
      background: none;
      border: none;
      cursor: pointer;
    }

    .selected-asv-item span {
    cursor: pointer;
    flex-grow: 1;
    padding: 5px;
  }

    /* .selected-asv-item span:hover,
    .selected-asv-item span:focus, */
    .selected-asv-item button:hover,
    .selected-asv-item button:focus {
    color: #ff3e00;
    outline: 2px solid #ff3e00
  }

  </style>