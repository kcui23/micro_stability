<script>
  import { onMount } from 'svelte';
  import { autoLoaded } from '../store.js';

  export let handleFileChange;
  export let handleGroupingsChange;
  export let updatePreVars;
  export let asvFiles = [];
  export let groupingsFile = null;

  let asvContent = '';
  let groupingsContent = '';

  function getFileName(files) {
    return files && files.length > 0 ? `Selected file: ${files[0].name}` : '';
  }

  function triggerFileInput(id) {
    document.getElementById(id).click();
  }

  async function handleASVFileChange(event) {
    const file = event.target.files[0];
    if (file) {
      asvFiles = [file];
      asvContent = await file.text();
      handleFileChange(event);
    }
  }

  async function handleGroupingsFileChange(event) {
    const file = event.target.files[0];
    if (file) {
      groupingsFile = file;
      groupingsContent = await file.text();
      handleGroupingsChange(event);
    }
  }

  async function uploadFiles() {
    if (asvContent && groupingsContent) {
      try {
        const response = await fetch('http://localhost:8000/store_files', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            asv: asvContent,
            groupings: groupingsContent,
          }),
        });

        if (response.ok) {
          const result = await response.json();
          console.log('Files stored successfully:', result);
          updatePreVars()
          return true;
        } else {
          console.error('Failed to store files');
          return false;
        }
      } catch (error) {
        console.error('Error uploading files:', error);
        return false;
      }
    } else {
      console.warn('Please select both ASV and groupings files before uploading');
      return false;
    }
  }

  async function autoLoadFiles() {
    if ($autoLoaded) {
      return;
    }
    const asvPath = '../../datasets/Blueberry/Blueberry_ASVs_table.tsv';
    const groupingPath = '../../datasets/Blueberry/Blueberry_metadata.tsv';

    try {
      const asvResponse = await fetch(asvPath);
      const groupingResponse = await fetch(groupingPath);

      if (asvResponse.ok && groupingResponse.ok) {
        const asvFile = new File([await asvResponse.blob()], 'Blueberry_ASVs_table.tsv');
        const groupingFile = new File([await groupingResponse.blob()], 'Blueberry_metadata.tsv');

        asvFiles = [asvFile];
        groupingsFile = groupingFile;

        asvContent = await asvFile.text();
        groupingsContent = await groupingFile.text();

        handleFileChange({ target: { files: [asvFile] } });
        handleGroupingsChange({ target: { files: [groupingFile] } });

        const uploadSuccess = await uploadFiles();
        if (uploadSuccess) {
          console.log('Files auto-loaded and uploaded successfully');
          autoLoaded.set(true);
          updatePreVars();
        } else {
          console.error('Failed to upload auto-loaded files');
        }
      }
    } catch (error) {
      console.error('Error auto-loading files:', error);
    }
  }

  onMount(() => {
    autoLoadFiles();
  });
</script>

<div class="upload-container">
  <!-- ASV file upload section -->
  <div class="upload-section">
    <span class="file-label">Upload ASV File:</span>
    <div class="button-group">
      <button type="button" on:click={() => triggerFileInput('fileInput1')}>Choose File</button>
      <button type="button" class="clear-button" on:click={() => { asvFiles = []; asvContent = ''; }}>Clear</button>
    </div>
    <input id="fileInput1" type="file" accept=".tsv" on:change={handleASVFileChange} style="display: none;" />
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName1" class="file-name">{getFileName(asvFiles)}</div>
  </div>
  
  <!-- Groupings file upload section -->
  <div class="upload-section">
    <span class="file-label">Upload Groupings File:</span>
    <div class="button-group">
      <button type="button" on:click={() => triggerFileInput('fileInput2')}>Choose File</button>
      <button type="button" class="clear-button" on:click={() => { groupingsFile = null; groupingsContent = ''; }}>Clear</button>
    </div>
    <input id="fileInput2" type="file" accept=".tsv" on:change={handleGroupingsFileChange} style="display: none;" />
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName2" class="file-name">{getFileName(groupingsFile ? [groupingsFile] : [])}</div>
  </div>
</div>

<button on:click={uploadFiles} class="upload-button">Upload Files</button>

<style>
  .upload-container {
    display: flex;
    justify-content: space-between;
    margin-bottom: 20px;
  }

  .upload-section {
    flex: 1;
    margin-right: 20px;
  }

  .button-group {
    display: flex;
    align-items: center;
  }

  .upload-button {
    margin-left: 40%;
    margin-right: auto;
  }

  button {
    background-color: #f0f0f0;
    color: black;
    padding: 10px 16px;
    font-size: 12px;
    font-weight: bold;
    text-align: center;
    border-radius: 4px;
    border: none;
    cursor: pointer;
    transition: background-color 0.1s ease;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-right: 10px;
  }

  button:hover {
    background-color: #e2e8ef;
  }

  button:active {
    background-color: #b6c2ce;
  }

  .clear-button {
    background-color: #f8d7da;
    color: #721c24;
  }

  .clear-button:hover {
    background-color: #f5c6cb;
  }

  .clear-button:active {
    background-color: #f1b0b7;
  }

  .file-label {
    text-align: left;
    display: block;
    margin-bottom: 5px;
    font-weight: bold;
  }

  .note {
    text-align: left;
    font-size: 0.9em;
    color: gray;
    margin-top: 5px;
  }

  .file-name {
    text-align: left;
    font-size: 0.9em;
    color: green;
    margin-top: 5px;
  }
</style>
