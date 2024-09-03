<script>
  export let handleFileChange;
  export let handleGroupingsChange;
  export let asvFiles = [];
  export let groupingsFile = null;

  function getFileName(files) {
    return files && files.length > 0 ? `Selected file: ${files[0].name}` : '';
  }

  function triggerFileInput(id) {
    document.getElementById(id).click();
  }
</script>

<div class="upload-container">
  <div class="upload-section">
    <span class="file-label">Upload ASV File:</span>
    <button type="button" on:click={() => triggerFileInput('fileInput1')}>Choose File</button>
    <input id="fileInput1" type="file" accept=".tsv" on:change={handleFileChange} style="display: none;" />
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName1" class="file-name">{getFileName(asvFiles)}</div>
  </div>

  <div class="upload-section">
    <span class="file-label">Upload Groupings File:</span>
    <button type="button" on:click={() => triggerFileInput('fileInput2')}>Choose File</button>
    <input id="fileInput2" type="file" accept=".tsv" on:change={handleGroupingsChange} style="display: none;" />
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName2" class="file-name">{getFileName(groupingsFile ? [groupingsFile] : [])}</div>
  </div>
</div>

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

  button {
    background-color: #f0f0f0;
    color: black;
    padding: 10px 16px;
    font-size: 12px;
    font-weight: bold;
    text-align: center;
    border-radius: 4px;
    border: 1px solid #ccc;
    cursor: pointer;
    box-shadow: 0px 2px 6px rgba(0, 0, 0, 0.1);
    transition: background-color 0.1s ease;

    display: flex;
    align-items: center;
    justify-content: center;
  }

  button:hover {
    background-color: #e2e8ef;
  }

  button:active {
    background-color: #b6c2ce;
  }

  .file-label {
    display: block;
    margin-bottom: 5px;
    font-weight: bold;
  }

  .note {
    font-size: 0.9em;
    color: gray;
    margin-top: 5px;
  }

  .file-name {
    font-size: 0.9em;
    color: green;
    margin-top: 5px;
  }
</style>
