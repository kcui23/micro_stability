<script>
  let files = [];
  let visualizations = { plot1: '', plot2: '', plot3: '' };
  let selectedMethod = '';
  let groupingsFile = null;
  let currentStep = 'Raw data';

  const handleFileChange = (event) => {
    files = event.target.files;
  };

  const handleGroupingsChange = (event) => {
    groupingsFile = event.target.files[0];
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
            groupings: groupingsContent
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

  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Prediction Evaluation Metric', 'Stability Metric'];

  const goToStep = (step) => {
    currentStep = step;
  };

  const previewFileContent = (fileContent) => {
    // Function to display a preview of the TSV file content
    const rows = fileContent.split('\n').slice(0, 5); // Get first 5 rows
    return rows.map(row => row.split('\t').slice(0, 5)); // Get first 5 columns of each row
  };

  let previewContent = [];
  if (files.length > 0) {
    const reader = new FileReader();
    reader.onload = () => {
      previewContent = previewFileContent(reader.result);
    };
    reader.readAsText(files[0]);
  }
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

  <!-- File Upload and Preview Section -->
  <div class="upload-section">
    <h1>Upload ASV Dataset</h1>
    <input type="file" accept=".tsv" on:change={handleFileChange} />
    <input type="file" accept=".tsv" on:change={handleGroupingsChange} />
    
    {#if files.length > 0 && groupingsFile}
      <div class="preview">
        <h2>Preview of TSV File</h2>
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
  </div>

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
