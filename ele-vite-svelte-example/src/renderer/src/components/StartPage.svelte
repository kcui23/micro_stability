<script>
  import FileUploader from './FileUploader.svelte';

  export let handleFileChange;
  export let handleGroupingsChange;
  export let updatePreVars;
  export let asvFiles;
  export let groupingsFile;
  export let startMethod;
  export let startApp;
  export let startAppSkip;
  export let DataPerturbationMethods;
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
    border: 1px solid #c668d9;
  }

  .start-page-content {
    background: #ffffff;
    padding: 50px 30px;
    text-align: center;
    max-width: 900px;
    width: 95%;
    border: 1px solid #68d9d1;
  }

  .start-page-content h2 {
    margin-bottom: 15px;
    font-size: 2.2rem;
    color: #333;
  }

  .start-page-content p {
    text-align: left;
    margin: 12px 0;
    font-size: 1.3rem;
    color: #555;
  }

  .start-page-content button {
    width: 120px;
    height: 40px;
    font-size: 1.1rem;
    margin: 25px 8px 0;
    border: none;
    border-radius: 4px;
    cursor: pointer;
    transition: background-color 0.3s ease;
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
    border: 1px solid #d1d968;
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
    border: 1px solid #ccc;
    border-radius: 4px;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  .start-page-content button:first-of-type {
    background-color: #b7f3b7;
  }

  .start-page-content button:first-of-type:hover {
    background-color: #90ee90;
  }

  .start-page-content button:last-of-type {
    background-color: #ffb6b6;
  }

  .start-page-content button:last-of-type:hover {
    background-color: #ff9999;
  }
</style>

<div class="start-page">
  <div class="start-page-content">
    <h2>Welcome to the Micro Stability App</h2>
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
    <button on:click={startApp}>Start</button>
    <button on:click={startAppSkip}>Skip</button>
  </div>
</div>