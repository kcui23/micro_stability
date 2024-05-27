<script>
  let files = [];
  let visualizations = { plot1: '', plot2: '', plot3: '' };
  let selectedMethod = '';
  let groupingsFile = null;
  let currentStep = 'Raw data';
  let previewContent = [];
  let groupingsContentPreview = [];
  let fileDimensions = { rows: 0, columns: 0 };
  let groupingsDimensions = { rows: 0, columns: 0 };
  let threshold = 0.0;
  let filteredContent = [];
  let filteredDimensions = { rows: 0, columns: 0 };

  const handleFileChange = (event) => {
    files = event.target.files;
    if (files.length > 0) {
      const reader = new FileReader();
      reader.onload = () => {
        const content = reader.result;
        previewContent = previewFileContent(content);
        fileDimensions = getFileDimensions(content);
        filteredContent = previewContent;
        filteredDimensions = fileDimensions;
      };
      reader.readAsText(files[0]);
    }
  };

  const handleGroupingsChange = (event) => {
    groupingsFile = event.target.files[0];
    if (groupingsFile) {
      const reader = new FileReader();
      reader.onload = () => {
        const content = reader.result;
        groupingsContentPreview = previewFileContent(content);
        groupingsDimensions = getFileDimensions(content);
      };
      reader.readAsText(groupingsFile);
    }
  };

  const handleMethodChange = (event) => {
    selectedMethod = event.target.value;
  };

  const handleSubmit = async () => {
    const file = files[0];
    const groupings = groupingsFile;

    const asvReader = new FileReader();
    const groupingsReader = new FileReader();

    asvReader.onload = () => {
      const asvContent = asvReader.result;

      groupingsReader.onload = async () => {
        const groupingsContent = groupingsReader.result;

        const response = await fetch(`http://localhost:8000/process?method=${selectedMethod}`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            asv: asvContent,
            groupings: groupingsContent,
            threshold: threshold
          })
        });

        const result = await response.json();
        visualizations = {
          plot1: `data:image/png;base64,${result.plot1}`,
          plot2: `data:image/png;base64,${result.plot2}`,
          plot3: `data:image/png;base64,${result.plot3}`
        };
      };

      groupingsReader.readAsText(groupings);
    };

    asvReader.readAsText(file);
  };

  const handleFilter = async () => {
    const file = files[0];

    const asvReader = new FileReader();
    asvReader.onload = async () => {
      const asvContent = asvReader.result;

      const response = await fetch(`http://localhost:8000/filter`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          asv: asvContent,
          threshold: threshold
        })
      });

      const result = await response.json();
      filteredContent = previewFileContent(result.filteredAsv);
      filteredDimensions = getFileDimensions(result.filteredAsv);
    };

    asvReader.readAsText(file);
  };

  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Prediction Evaluation Metric', 'Stability Metric'];

  const goToStep = (step) => {
    if (currentStep === 'Raw data' && (files.length === 0 || !groupingsFile)) {
      return; // Do not allow navigation if files are not uploaded
    }
    currentStep = step;
  };

  const previewFileContent = (fileContent) => {
    const rows = fileContent.split('\n').slice(0, 5); // Get first 5 rows
    return rows.map(row => row.split('\t').slice(0, 5)); // Get first 5 columns of each row
  };

  const getFileDimensions = (fileContent) => {
    const rows = fileContent.split('\n');
    const columns = rows[0].split('\t').length;
    return { rows: rows.length, columns };
  };
</script>

<style>
  .steps {
    display: flex;
    justify-content: space-around;
    margin-bottom: 20px;
  }
  .step {
    cursor: pointer;
    padding: 10px;
    border: 1px solid #ccc;
    border-radius: 5px;
  }
  .step.active {
    background-color: #007bff;
    color: white;
  }
  .preview {
    margin: 20px 0;
  }
  .filters {
    margin: 20px 0;
  }
  .methods {
    display: flex;
    flex-direction: column;
    margin: 20px 0;
  }
  .methods button {
    margin: 5px 0;
  }
  .navigation {
    display: flex;
    justify-content: space-between;
    margin: 20px 0;
  }
  table {
    width: 100%;
    border-collapse: collapse;
  }
  th, td {
    border: 1px solid #ccc;
    padding: 8px;
    text-align: left;
  }
  th {
    background-color: #f2f2f2;
  }
</style>

<div>
  <!-- Logo Section -->
  <div class="logo">
    <h1>Logo</h1>
  </div>

  <!-- Steps Navigation -->
  <div class="steps">
    {#each steps as step}
      <div class:active={step === currentStep} class="step" on:click={() => goToStep(step)}>
        {step}
      </div>
    {/each}
  </div>

  <!-- Content for Current Step -->
  <div class="upload-section" hidden={currentStep !== 'Raw data'}>
    <h1>Upload ASV Dataset</h1>
    <input type="file" accept=".tsv" on:change={handleFileChange} />
    <input type="file" accept=".tsv" on:change={handleGroupingsChange} />

    {#if files.length > 0}
      <div class="preview">
        <h2>Preview of TSV File</h2>
        <p>Dimensions: {fileDimensions.rows} rows, {fileDimensions.columns} columns</p>
        <table>
          {#each previewContent as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
        </table>
      </div>
    {/if}

    {#if groupingsFile}
      <div class="preview">
        <h2>Preview of Groupings File</h2>
        <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns} columns</p>
        <table>
          {#each groupingsContentPreview as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
        </table>
      </div>
    {/if}
  </div>

  <!-- Data Perturbation Step Preview and Filters -->
  {#if currentStep === 'Data Perturbation'}
    {#if files.length > 0}
      <div class="preview">
        <h2>Preview of TSV File</h2>
        <p>Dimensions: {filteredDimensions.rows} rows, {filteredDimensions.columns} columns</p>
        <table>
          {#each filteredContent as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
        </table>
      </div>
    {/if}

    {#if groupingsFile}
      <div class="preview">
        <h2>Preview of Groupings File</h2>
        <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns} columns</p>
        <table>
          {#each groupingsContentPreview as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
        </table>
      </div>
    {/if}

    <!-- Placeholder for Filtering -->
    <div class="filters">
      <label for="threshold">Threshold for Rare Genes:</label>
      <input type="range" id="threshold" bind:value={threshold} min="0" max="1" step="0.01" />
      <span>{threshold}</span>
      <button on:click={handleFilter}>Apply Threshold</button>
    </div>

    <!-- Placeholders for additional options -->
    <div class="filters">
      <label>Additional Option 1:</label>
      <input type="text" />
    </div>

    <div class="filters">
      <label>Additional Option 2:</label>
      <input type="text" />
    </div>
  {/if}

  <!-- Method Selection Section -->
  {#if currentStep === 'Model Perturbation'}
    <div class="methods">
      <button on:click={() => handleMethodChange({ target: { value: 'deseq2' } })}>Method 1 (DESeq2)</button>
      <button on:click={() => handleMethodChange({ target: { value: 'aldex2' } })}>Method 2 (ALDEx2)</button>
      <button on:click={() => handleMethodChange({ target: { value: 'method3' } })}>Method 3</button>
      <button on:click={() => handleMethodChange({ target: { value: 'method4' } })}>Method 4</button>
      <button on:click={() => handleMethodChange({ target: { value: 'method5' } })}>Method 5</button>
    </div>
  {/if}

  <!-- Visualizations Section -->
  {#if visualizations.plot1}
    <h2>Visualizations</h2>
    <img src={visualizations.plot1} alt="Plot 1" style="width: 300px; height: auto;" />
    <img src={visualizations.plot2} alt="Plot 2" style="width: 300px; height: auto;" />
    <img src={visualizations.plot3} alt="Plot 3" style="width: 300px; height: auto;" />
  {/if}

  <!-- Navigation Buttons -->
  <div class="navigation">
    <button on:click={() => goToStep(steps[Math.max(0, steps.indexOf(currentStep) - 1)])}>Previous</button>
    <button on:click={() => goToStep(steps[Math.min(steps.length - 1, steps.indexOf(currentStep) + 1)])}>Next</button>
  </div>

  <!-- Submit Button -->
  {#if files.length > 0 && groupingsFile && selectedMethod}
    <button on:click={handleSubmit}>Submit</button>
  {/if}
</div>
