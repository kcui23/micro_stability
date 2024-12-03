<script>
  import { onMount, onDestroy } from 'svelte';
  import FileUploader from './FileUploader.svelte';
  import { quintOut } from 'svelte/easing';
  import { slide } from 'svelte/transition';
  import { fileUploaded, StabilityMetricCalFinished } from '../store.js';

  export let handleFileChange;
  export let handleGroupingsChange;
  export let updatePreVars;
  export let asvFiles;
  export let groupingsFile;
  export let startApp;
  export let missingMethods;
  export let showStartPage;
  let startSelectedMethods = ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'];
  let DataPerturbationMethods = ['DESeq2', 'EdgeR', 'Maaslin2', 'Aldex2', 'MetagenomeSeq'];
  let raw_total_data = null;
  // const test_json_files = async () => {
  //   const res = await fetch(`http://localhost:8000/test_json_files`)
  //   raw_total_data = await res.json()
  // }
  // let result_data = []
  // const testdeseq2 = async () => {
  //   const res = await fetch(`http://localhost:8000/test-deseq2`)
  //   result_data = await res.json()
  // }

  const sections = [
    {
      title: 'Filtering',
      items: ['Low Abundance', 'Prevalence', 'Variance', 'No Filtering']
    },
    {
      title: 'Zero-Handling',
      items: ['Pseudocount Addition', 'k-NN Imputation', 'No Zero-Handling']
    },
    {
      title: 'Normalization',
      items: ['TSS', 'CSS', 'TMM', 'CLR', 'No Normalization']
    },
    {
      title: 'Transformation',
      items: ['Log', 'Logit', 'AST', 'No Transformation']
    },
    {
      title: 'Model Perturbation',
      items: ['DESeq2', 'EdgeR', 'Maaslin2', 'Aldex2', 'MetagenomeSeq']
    }
  ];

  let openSections = {};
  sections.forEach(section => {
    openSections[section.title] = false;
  });

  function toggleSection(title) {
    openSections[title] = !openSections[title];
    openSections = openSections; // Trigger reactivity
  }

  let isSubmitting = false;
  let isSubmitted = false;

  async function handleSubmit() {
    isSubmitting = true;
    missingMethods = startSelectedMethods;
    startApp();
    setTimeout(() => {
      isSubmitting = false;
      isSubmitted = true;
    }, 2000);
  }

  $: progress = missingMethods ? 
    ((startSelectedMethods.length - missingMethods.length ) / startSelectedMethods.length) * 100 : 
    0;

  let isOpen = false;

  function handleClickOutside(event) {
    const dropdown = document.querySelector('.dropdown');
    if (dropdown && !dropdown.contains(event.target)) {
      isOpen = false;
    }
  }

  let loaded_web_json_files = false;
  async function loadWebJsonFiles() {
    try {
      const response = await fetch('http://localhost:8000/load_web_json_files', {
        method: 'POST'
      });

      if (!response.ok) {
        throw new Error('Failed to load web JSON files');
      }
      const result = await response.json();
      console.log('Web JSON files loaded successfully:', result);
      loaded_web_json_files = true;
      return true;
    } catch (error) {
      console.error('Error loading web JSON files:', error);
      return false;
    }
  }

  onMount(async () => {
    // First load web JSON files
    // const jsonFilesLoaded = await loadWebJsonFiles();
    
    if (!jsonFilesLoaded) {
      console.error('Failed to load initial JSON files');
      return;
    }
    document.addEventListener('click', handleClickOutside);
  });

  onDestroy(() => {
    document.removeEventListener('click', handleClickOutside);
  });
</script>

<style>
  .start-page {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background: #ffffff;
    z-index: 1000;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
  }

  .start-page-content {
    background: #ffffff;
    padding: 50px 30px;
    width: 80%;
  }

  .start-page-content h2 {
    margin-bottom: 15px;
    font-size: 2.2rem;
    color: #333;
    text-align: center;
  }

  .columns-container {
    display: flex;
    justify-content: space-between;
    gap: 20px;
    margin-top: 20px;
  }

  .left-column {
    flex: 0 0 70%;
  }

  .right-column {
    flex: 0 0 30%;
    max-height: 400px;
    overflow-y: hidden;
    margin-left: 10px;
  }

  .start-page-content p {
    text-align: left;
    margin: 12px 0;
    font-size: 1.3rem;
    color: #555;
  }

  .submit-button {
    position: relative;
    overflow: hidden;
    background-color: #f0f0f0;
    width: 120px;
    height: 40px;
    font-size: 1.1rem;
    margin: 25px 8px 0;
    border: none;
    border-radius: 4px;
    transition: all 0.3s ease;
  }

  .submit-button.submitting {
    background-color: #4CAF50;
    color: white;
    cursor: wait;
  }

  .submit-button.submitted {
    background-color: #45a049;
    color: white;
    cursor: default;
  }

  .submit-button.submitting::after {
    content: '';
    position: absolute;
    top: 0;
    left: -100%;
    width: 100%;
    height: 100%;
    background: linear-gradient(
      90deg,
      transparent,
      rgba(255, 255, 255, 0.2),
      transparent
    );
    animation: loading 1.5s infinite;
  }

  .submit-button:not(:disabled):hover {
    background-color: #e2e8ef;
  }

  .submit-button:not(:disabled):active {
    background-color: #b6c2ce;
  }

  @keyframes loading {
    100% {
      left: 100%;
    }
  }

  .go-to-main-page-button {
    position: relative;
    overflow: hidden;
    background-color: #f0f0f0;
    width: 155px;
    height: 40px;
    font-size: 1.1rem;
    margin-top: 25px;
    border: none;
    border-radius: 4px;
    transition: all 0.3s ease;
  }

  .go-to-main-page-button:hover {
    background-color: #e2e8ef;
  }

  .go-to-main-page-button:active {
    background-color: #b6c2ce;
  }

  .step-1-upload-section {
    position: relative;
    padding: 15px;
    margin-left: 15px;
  }

  .step-1-upload-section::before {
    content: '';
    position: absolute;
    border-radius: 8px;
    top: 0;
    left: 0;
    height: 100%;
    width: calc(100% - 50px);
    background-color: #636363;
    z-index: -1;
  }

  .step-2 {
    display: flex;
    align-items: center;
    gap: 8px;
    flex-wrap: wrap;
    margin-top: 15px;
  }

  .step-2 p {
    margin: 0;
    font-size: 1.3rem;
    color: #555;
    white-space: nowrap;
  }

.multiverse-column {
    background: #f8f8f8;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    height: 100%;
    overflow-y: auto;
  }

  .multiverse-column h3 {
    color: #333;
    font-size: 1.6rem;
    margin-bottom: 10px;
  }

  .multiverse-column ul {
    list-style: none;
    padding-left: 10px;
    margin: 8px 0;
  }

  .multiverse-column li {
    color: #666;
    padding: 4px 0;
    font-size: 1.1rem;
    transition: color 0.2s ease;
  }

  .multiverse-column li:hover {
    color: #444;
  }

  .details-wrapper {
    margin-bottom: 8px;
  }

  .summary-button {
    width: 100%;
    text-align: left;
    padding: 8px 12px;
    background: none;
    border: none;
    font-size: 1.2rem;
    font-weight: 500;
    color: #555;
    cursor: pointer;
    display: flex;
    justify-content: space-between;
    align-items: center;
    transition: color 0.2s ease;
  }

  .summary-button:hover {
    color: #333;
  }

  .arrow {
    transition: transform 0.3s ease;
    font-size: 0.8em;
  }

  .arrow.open {
    transform: rotate(180deg);
  }

  .content {
    padding: 0 12px;
    overflow: hidden;
  }

  .content ul {
    list-style: none;
    padding-left: 10px;
    margin: 8px 0;
  }

  .content li {
    color: #666;
    padding: 4px 0;
    font-size: 1.1rem;
    transition: color 0.2s ease;
  }

  .content li:hover {
    color: #444;
  }

  .submit-section {
    display: flex;
    align-items: center;
    gap: 20px;
    margin-top: 25px;
  }

  .progress-section {
    display: flex;
    align-items: center;
    gap: 10px;
    margin-top: 25px;
    margin-left: 15px;
  }

  .progress-container {
    width: 400px;
    height: 20px;
    background-color: #f0f0f0;
    border-radius: 10px;
    overflow: hidden;
  }

  .progress-bar {
    height: 100%;
    background-color: #4CAF50;
    transition: width 0.3s ease;
  }

  .progress-text {
    font-size: 1.1rem;
    color: #666;
    min-width: 48px;
  }

  .checkbox-group {
    display: flex;
    flex-direction: column;
    gap: 10px;
    margin: 10px 0;
  }

  .checkbox-label {
    display: flex;
    align-items: center;
    gap: 8px;
    cursor: pointer;
    user-select: none;
  }

  .checkbox-label input[type="checkbox"] {
    width: 16px;
    height: 16px;
    cursor: pointer;
  }

  .dropdown {
    position: relative;
    display: inline-block;
  }

  .dropdown-button {
    height: 30px;
    margin-top: 5px;
    background-color: #f0f0f0;
    border: 1px solid #ddd;
    border-radius: 4px;
    cursor: pointer;
    font-size: 1rem;
    min-width: 200px;
    text-align: left;
  }

  .dropdown-button:hover {
    background-color: #e2e8ef;
  }

  .dropdown-content {
    position: absolute;
    top: 100%;
    left: 0;
    background-color: white;
    border: 1px solid #ddd;
    border-radius: 4px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    z-index: 1000;
    margin-top: 4px;
    min-width: 200px;
    max-height: 300px;
    overflow-y: auto;
  }

  .checkbox-group {
    padding: 8px;
  }

  .checkbox-label {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 6px 8px;
    cursor: pointer;
    user-select: none;
  }

  .checkbox-label:hover {
    background-color: #f5f5f5;
  }

  .checkbox-label input[type="checkbox"] {
    width: 16px;
    height: 16px;
    cursor: pointer;
  }
</style>

<div class="start-page">
  <div class="start-page-content">
    <h2>Welcome to the Micro Stability App</h2>
    <div class="columns-container">
      <div class="left-column">
        <p><strong>Step 1:</strong> Upload the file, update the ASV file, and grouping file.</p>
        <div class="step-1-upload-section">
          <FileUploader 
            {handleFileChange} 
            {handleGroupingsChange} 
            {updatePreVars}
            bind:asvFiles
            bind:groupingsFile
          />
        </div>

        <div class="step-2">
          <p><strong>Step 2:</strong> Choose methods you like to start: </p>
          <div class="dropdown">
            <button class="dropdown-button" on:click|preventDefault={() => isOpen = !isOpen}>
              Select Methods ({startSelectedMethods.length} selected) ▼
            </button>
            
            {#if isOpen}
              <div class="dropdown-content" transition:slide|local={{ duration: 200 }}>
                <div class="checkbox-group">
                  {#each DataPerturbationMethods as method}
                    <label class="checkbox-label">
                      <input
                        type="checkbox"
                        value={method}
                        checked={startSelectedMethods.includes(method.toLowerCase())}
                        on:change={(e) => {
                          if (e.target.checked) {
                            let tmp_method = method.toLowerCase();
                            startSelectedMethods = [...startSelectedMethods, tmp_method];
                          } else {
                            startSelectedMethods = startSelectedMethods.filter(m => m !== method.toLowerCase());
                          }
                        }}
                      />
                      {method}
                    </label>
                  {/each}
                </div>
              </div>
            {/if}
          </div>
        </div>
        <p><strong>Step 3:</strong> Start the app by submitting your job. It may take hours to complete.</p>
        <div class="submit-section">
          <!-- <button on:click={() => showStartPage = false}>Skip</button> -->
          <button 
            on:click={handleSubmit} 
            class="submit-button" 
            class:submitting={isSubmitting}
            class:submitted={isSubmitted}
            disabled={!$fileUploaded || isSubmitting || isSubmitted}
          >
            {#if isSubmitting}
              Submitting...
            {:else if isSubmitted}
              Submitted!
            {:else}
              Submit Job
            {/if}
          </button>
          
          {#if isSubmitted}
            <div class="progress-section">
              <div class="progress-container">
                <div class="progress-bar" style="width: {progress}%"></div>
              </div>
              <span class="progress-text">{Math.round(progress)}%</span>
            </div>
          {/if}
          
          {#if $StabilityMetricCalFinished}
            <button 
              class="go-to-main-page-button"
              on:click={() => showStartPage = false}
            >
              Go to Main Page
            </button>
          {/if}
        </div>
        <!-- <p>missingMethods: {missingMethods}</p>
        <p>startSelectedMethods: {startSelectedMethods}</p>
        <button on:click={testdeseq2}>test deseq2</button>
        <div style="max-height: 200px; overflow-y: auto; border: 1px solid #ccc; padding: 8px;">
          <p>result_data: {JSON.stringify(result_data)}</p>
        </div> -->
      </div>
      <div class="right-column">
        <div class="multiverse-column">
          <h3>Explore the Multiverse</h3>
          
          {#each sections as section}
            <div class="details-wrapper">
              <button 
                class="summary-button" 
                class:open={openSections[section.title]}
                on:click={() => toggleSection(section.title)}
              >
                {section.title}
                <span class="arrow" class:open={openSections[section.title]}>▼</span>
              </button>
              
              {#if openSections[section.title]}
                <div
                  class="content"
                  transition:slide|local={{ duration: 200, easing: quintOut }}
                >
                  <ul>
                    {#each section.items as item}
                      <li>{item}</li>
                    {/each}
                  </ul>
                </div>
              {/if}
            </div>
          {/each}
          <!--
          <div>
            <button 
              on:click={loadWebJsonFiles}
              style="background-color: {loaded_web_json_files ? '#4CAF50' : '#f44336'}; color: white;"
            >
              load web json files
            </button>
          </div>
          <button on:click={test_json_files}>test_json_files</button>
          <div style="max-height: 200px; overflow-y: auto; border: 1px solid #ccc; padding: 8px;">
            <p>raw_total_data: {JSON.stringify(raw_total_data)}</p>
          </div> -->
        </div>
      
      
      </div>
    </div>
  </div>
</div>
