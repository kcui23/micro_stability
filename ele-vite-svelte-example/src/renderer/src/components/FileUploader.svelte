<script>
  export let handleFileChange
  export let handleGroupingsChange
  export let asvFiles = []
  export let groupingsFile = null

  function getFileName(files) {
    return files && files.length > 0 ? `Selected file: ${files[0].name}` : ''
  }
</script>

<div class="upload-container">
  <div class="upload-section">
    <span class="file-label">Upload ASV File:</span>
    <div class="custom-file-input">
      <label for="fileInput1">Choose File</label>
      <input id="fileInput1" type="file" accept=".tsv" on:change={handleFileChange} />
    </div>
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName1" class="file-name">{getFileName(asvFiles)}</div>
  </div>

  <div class="upload-section">
    <span class="file-label">Upload Groupings File:</span>
    <div class="custom-file-input">
      <label for="fileInput2">Choose File</label>
      <input id="fileInput2" type="file" accept=".tsv" on:change={handleGroupingsChange} />
    </div>
    <div class="note">Note: Please upload a .tsv file</div>
    <div id="fileName2" class="file-name">{getFileName(groupingsFile ? [groupingsFile] : [])}</div>
  </div>
</div>

<style>
    .upload-container {
      display: flex;
      flex-direction: column;
    }
  
    .upload-section {
      margin-bottom: 20px;
    }
  
    .custom-file-input {
      position: relative;
      overflow: hidden;
      display: inline-block;
    }
  
    .custom-file-input input[type="file"] {
      position: absolute;
      width: 0;
      height: 0;
      opacity: 0;
      right: 0;
      top: 0;
    }
  
    .custom-file-input label {
      background-color: rgb(235, 235, 235);
      color: black;
      padding: 0 7px;
      font-size: 14px;
      text-align: center;
      height: 25px;
      line-height: 25px;
      border-radius: 2px;
      cursor: pointer;
      display: inline-block;
      border: 1px solid black;
    }
  
    .custom-file-input label:hover {
      background-color: rgb(223, 223, 223);
    }
  </style>