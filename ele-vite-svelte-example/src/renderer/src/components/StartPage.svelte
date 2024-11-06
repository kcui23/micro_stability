<script>
  import FileUploader from './FileUploader.svelte';
  import { quintOut } from 'svelte/easing';
  import { slide } from 'svelte/transition';
  import { fileUploaded } from '../store.js';

  export let handleFileChange;
  export let handleGroupingsChange;
  export let updatePreVars;
  export let asvFiles;
  export let groupingsFile;
  export let startMethod;
  export let startApp;
  export let DataPerturbationMethods;

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
</script>

<style>
  /* Start Page Styles */
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
  }

  .start-page-content p {
    text-align: left;
    margin: 12px 0;
    font-size: 1.3rem;
    color: #555;
  }

  .submit-button {
    background-color: #f0f0f0;
    width: 120px;
    height: 40px;
    font-size: 1.1rem;
    margin: 25px 8px 0;
    border: none;
    border-radius: 4px;
    transition: background-color 0.1s ease;
  }

  .submit-button:hover {
    background-color: #e2e8ef;
  }
  .submit-button:not([disabled]):active {
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

  .step-2 select {
    flex: 0 0 auto;
    width: 100px;
    max-width: 100px;
      margin-left: 10px;
      padding: 5px;
    font-size: 1rem;
    border-radius: 4px;
    overflow: hidden;
    text-overflow: ellipsis;
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
          <p><strong>Step 2:</strong> Choose a method you like to start: </p>
          <select bind:value={startMethod}>
            {#each DataPerturbationMethods as method}
              <option value={method}>{method}</option>
            {/each}
          </select>
        </div>
        <p><strong>Step 3:</strong> Start the app by submitting your job. It may take hours to complete.</p>
        <button on:click={startApp} class="submit-button" data-content="Submit Job" disabled={!$fileUploaded}>
          Submit Job
        </button>
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
        </div>
      
      
      </div>
    </div>
  </div>
</div>
