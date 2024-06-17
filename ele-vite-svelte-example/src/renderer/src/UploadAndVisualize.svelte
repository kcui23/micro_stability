<script>
  import { onMount } from 'svelte';
  import ADGPlot from './ADGPlot.svelte';
  let files = [];
  let visualizations = {
    deseq2_plot1: '',
    deseq2_plot2: '',
    deseq2_plot3: '',
    aldex2_plot1: '',
    aldex2_plot2: '',
    aldex2_plot3: '',
    overlap_volcano: '',
    overlap_pvalue_distribution: ''
  };
  let selectedMethod = '';
  let groupingsFile = null;
  let currentStep = 'Raw data';
  let lastStep = 'Raw data';
  let previewContent = [];
  let groupingsContentPreview = [];
  let fileDimensions = { rows: 0, columns: 0 };
  let groupingsDimensions = { rows: 0, columns: 0 };
  let threshold = 0.0;
  let filteredContent = [];
  let filteredDimensions = { rows: 0, columns: 0 };
  let filteredAsvContent = null;
  let showAllPlots = false;
  let showDetailedPlots = false;
  let isCalculating = false;
  let randomSeed = 1234;
  let edgeThicknesses = [1,2,3,4];

  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Prediction Evaluation Metric', 'Stability Metric'];

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
        filteredAsvContent = content;
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

  const handleMethodChange = (method) => {
    selectedMethod = method;
    showAllPlots = false;
    showDetailedPlots = false;
  };

  const handleQuickExplore = async () => {
    showAllPlots = true;
    isCalculating = true;
    const asvContent = filteredAsvContent || files[0];
    const groupings = groupingsFile;

    if (typeof asvContent === 'string') {
      await processQuickExplore(asvContent, groupings);
    } else {
      const asvReader = new FileReader();
      asvReader.onload = async () => {
        const asvContentText = asvReader.result;
        await processQuickExplore(asvContentText, groupings);
      };
      asvReader.readAsText(asvContent);
    }
  };

  const handleSubmit = async () => {
    isCalculating = true;
    const asvContent = filteredAsvContent || files[0];
    const groupings = groupingsFile;

    if (typeof asvContent === 'string') {
      await processSubmit(asvContent, groupings);
    } else {
      const asvReader = new FileReader();
      asvReader.onload = async () => {
        const asvContentText = asvReader.result;
        await processSubmit(asvContentText, groupings);
      };
      asvReader.readAsText(asvContent);
    }
  };

  const processSubmit = async (asvContentText, groupings) => {
    const groupingsReader = new FileReader();
    groupingsReader.onload = async () => {
      const groupingsContent = groupingsReader.result;

      try {
        const response = await fetch(`http://localhost:8000/process?method=${selectedMethod}`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            asv: asvContentText,
            groupings: groupingsContent,
            threshold: threshold,
            seed: randomSeed
          })
        });

        if (!response.ok) {
          throw new Error('Network response was not ok');
        }

        const result = await response.json();
        visualizations = {
          deseq2_plot1: selectedMethod === 'deseq2' ? `data:image/png;base64,${result.plot1}` : visualizations.deseq2_plot1,
          deseq2_plot2: selectedMethod === 'deseq2' ? `data:image/png;base64,${result.plot2}` : visualizations.deseq2_plot2,
          deseq2_plot3: selectedMethod === 'deseq2' ? `data:image/png;base64,${result.plot3}` : visualizations.deseq2_plot3,
          aldex2_plot1: selectedMethod === 'aldex2' ? `data:image/png;base64,${result.plot1}` : visualizations.aldex2_plot1,
          aldex2_plot2: selectedMethod === 'aldex2' ? `data:image/png;base64,${result.plot2}` : visualizations.aldex2_plot2,
          aldex2_plot3: selectedMethod === 'aldex2' ? `data:image/png;base64,${result.plot3}` : visualizations.aldex2_plot3,
          overlap_volcano: visualizations.overlap_volcano,
          overlap_pvalue_distribution: visualizations.overlap_pvalue_distribution
        };
        isCalculating = false;
      } catch (error) {
        console.error('Fetch error:', error);
        isCalculating = false;
      }
    };

    groupingsReader.readAsText(groupings);
  };

  const processQuickExplore = async (asvContentText, groupings) => {
    const groupingsReader = new FileReader();
    groupingsReader.onload = async () => {
      const groupingsContent = groupingsReader.result;

      try {
        const response = await fetch(`http://localhost:8000/quick_explore`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            asv: asvContentText,
            groupings: groupingsContent,
            method: selectedMethod,
            seed: randomSeed
          })
        });

        if (!response.ok) {
          throw new Error('Network response was not ok');
        }

        const result = await response.json();
        visualizations = {
          deseq2_plot1: `data:image/png;base64,${result.deseq2_plot1}`,
          deseq2_plot2: `data:image/png;base64,${result.deseq2_plot2}`,
          deseq2_plot3: `data:image/png;base64,${result.deseq2_plot3}`,
          aldex2_plot1: `data:image/png;base64,${result.aldex2_plot1}`,
          aldex2_plot2: `data:image/png;base64,${result.aldex2_plot2}`,
          aldex2_plot3: `data:image/png;base64,${result.aldex2_plot3}`,
          overlap_volcano: `data:image/png;base64,${result.overlap_volcano}`,
          overlap_pvalue_distribution: `data:image/png;base64,${result.overlap_pvalue_distribution}`
        };
        isCalculating = false;
      } catch (error) {
        console.error('Fetch error:', error);
        isCalculating = false;
      }
    };

    groupingsReader.readAsText(groupings);
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
      const filteredAsv = result.filteredAsv;

      if (typeof filteredAsv === 'string') {
        filteredContent = previewFileContent(filteredAsv);
        filteredDimensions = getFileDimensions(filteredAsv);
        filteredAsvContent = filteredAsv;
      } else {
        console.error('Filtered ASV content is not a string:', filteredAsv);
        if (Array.isArray(filteredAsv)) {
          const filteredAsvString = filteredAsv.join('\n');
          filteredContent = previewFileContent(filteredAsvString);
          filteredDimensions = getFileDimensions(filteredAsvString);
          filteredAsvContent = filteredAsvString;
        }
      }
    };

    asvReader.readAsText(file);
  };

  const goToStep = (step) => {
    if (currentStep === 'Raw data' && (files.length === 0 || !groupingsFile)) {
      return;
    }
    lastStep = currentStep;
    currentStep = step;

    const contentContainer = document.querySelector('.content');
    if (contentContainer) {
      contentContainer.scrollTo(0, 0);
    }
  };

  const previewFileContent = (fileContent) => {
    const rows = fileContent.split('\n').slice(0, 5);
    return rows.map(row => row.split('\t').slice(0, 5));
  };

  const getFileDimensions = (fileContent) => {
    const rows = fileContent.split('\n');
    const columns = rows[0].split('\t').length;
    return { rows: rows.length, columns };
  };

  const showMoreDetails = () => {
    showDetailedPlots = true;
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
  .methods button.selected {
    background-color: #007bff;
    color: white;
  }
  .navigation {
    display: flex;
    justify-content: space-between;
    margin: 20px 40px 20px 0;
  }
</style>

<div id="app" class="container">
  <!-- ADG Sidebar -->
  <div class="sidebar">
    <ADGPlot steps={steps} currentStep={currentStep} setCurrentStep={goToStep} edgeThicknesses={edgeThicknesses} lastStep={lastStep} />
  </div>

  <!-- Main Content Area -->
  <div class="content">
    <div class="logo">
      <h1>Logo</h1>
    </div>

    <!-- <div class="steps">
      {#each steps as step}
        <button class:active={step === currentStep} class="step" on:click={() => goToStep(step)}>
          {step}
        </button>
      {/each}
    </div> -->

    <div class="upload-section" hidden={currentStep !== 'Raw data'}>
      <h1>Upload ASV Dataset</h1>
      <input type="file" accept=".tsv" on:change={handleFileChange} />
      <input type="file" accept=".tsv" on:change={handleGroupingsChange} />
      <div>
        <label for="randomSeed">Set Random Seed:</label>
        <input type="number" id="randomSeed" bind:value={randomSeed} min="1" />
      </div>

      {#if files.length > 0}
        <div class="preview">
          <h2>Preview of TSV File</h2>
          <p>Dimensions: {fileDimensions.rows} rows, {fileDimensions.columns} columns</p>
          <div class="table-container">
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
        </div>
      {/if}

      {#if groupingsFile}
        <div class="preview">
          <h2>Preview of Groupings File</h2>
          <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns}</p>
          <div class="table-container">
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
        </div>
      {/if}

      <div class="tooltip">
        <button on:click={handleQuickExplore} disabled={files.length === 0 || !groupingsFile}>Quick Explore</button>
        <span class="tooltiptext">Upload files to continue</span>
      </div>
    </div>

    <div class="preview-section" hidden={currentStep !== 'Data Perturbation'}>
      {#if files.length > 0}
        <div class="preview">
          <h2>Preview of TSV File</h2>
          <p>Dimensions: {filteredDimensions.rows} rows, {filteredDimensions.columns} columns</p>
          <div class="table-container">
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
        </div>
      {/if}

      {#if groupingsFile}
        <div class="preview">
          <h2>Preview of Groupings File</h2>
          <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns}</p>
          <div class="table-container">
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
        </div>
      {/if}

      <div class="filters">
        <label for="threshold">Threshold for Rare Genes:</label>
        <input type="range" id="threshold" bind:value={threshold} min="0" max="1" step="0.01" />
        <span>{threshold}</span>
        <button on:click={handleFilter}>Apply Threshold</button>
      </div>

      <div class="filters">
        <label for="option1">Additional Option 1:</label>
        <input type="text" id="option1" />
      </div>

      <div class="filters">
        <label for="option2">Additional Option 2:</label>
        <input type="text" id="option2" />
      </div>
    </div>

    <div class="methods-section" hidden={currentStep !== 'Model Perturbation'}>
      <div class="methods">
        <button on:click={() => handleMethodChange('deseq2')} class:selected={selectedMethod === 'deseq2'}>Method 1 (DESeq2)</button>
        <button on:click={() => handleMethodChange('aldex2')} class:selected={selectedMethod === 'aldex2'}>Method 2 (ALDEx2)</button>
        <button on:click={() => handleMethodChange('method3')} class:selected={selectedMethod === 'method3'}>Method 3</button>
        <button on:click={() => handleMethodChange('method4')} class:selected={selectedMethod === 'method4'}>Method 4</button>
        <button on:click={() => handleMethodChange('method5')} class:selected={selectedMethod === 'method5'}>Method 5</button>
      </div>
    </div>

    <div class="visualizations-section" hidden={!showAllPlots && !selectedMethod}>
      <h2>Visualizations</h2>
      {#if isCalculating}
        <div class="loader">
          <p>Loading...</p>
        </div>
      {:else}
        {#if showAllPlots}
          <h3>Overlap Visualizations</h3>
          <img class="large" src={visualizations.overlap_volcano} alt="Overlap Volcano Plot" />
          <img class="large" src={visualizations.overlap_pvalue_distribution} alt="Overlap P-value Distribution" />
          {#if !showDetailedPlots}
            <button on:click={showMoreDetails}>Show More Details</button>
          {/if}
          {#if showDetailedPlots}
            <h3>DESeq2 Plots</h3>
            <img src={visualizations.deseq2_plot1} alt="DESeq2 Plot 1" style="width: 300px; height: auto;" />
            <img src={visualizations.deseq2_plot2} alt="DESeq2 Plot 2" style="width: 300px; height: auto;" />
            <img src={visualizations.deseq2_plot3} alt="DESeq2 Plot 3" style="width: 300px; height: auto;" />
            <h3>ALDEx2 Plots</h3>
            <img src={visualizations.aldex2_plot1} alt="ALDEx2 Plot 1" style="width: 300px; height: auto;" />
            <img src={visualizations.aldex2_plot2} alt="ALDEx2 Plot 2" style="width: 300px; height: auto;" />
            <img src={visualizations.aldex2_plot3} alt="ALDEx2 Plot 3" style="width: 300px; height: auto;" />
          {/if}
        {:else if selectedMethod === 'deseq2'}
          <h3>DESeq2 Plots</h3>
          <img src={visualizations.deseq2_plot1} alt="DESeq2 Plot 1" style="width: 300px; height: auto;" />
          <img src={visualizations.deseq2_plot2} alt="DESeq2 Plot 2" style="width: 300px; height: auto;" />
          <img src={visualizations.deseq2_plot3} alt="DESeq2 Plot 3" style="width: 300px; height: auto;" />
        {:else if selectedMethod === 'aldex2'}
          <h3>ALDEx2 Plots</h3>
          <img src={visualizations.aldex2_plot1} alt="ALDEx2 Plot 1" style="width: 300px; height: auto;" />
          <img src={visualizations.aldex2_plot2} alt="ALDEx2 Plot 2" style="width: 300px; height: auto;" />
          <img src={visualizations.aldex2_plot3} alt="ALDEx2 Plot 3" style="width: 300px; height: auto;" />
        {/if}
      {/if}
    </div>

    <div class="navigation">
      <button on:click={() => goToStep(steps[Math.max(0, steps.indexOf(currentStep) - 1)])}>Previous</button>
      <div class="tooltip">
        <button on:click={() => goToStep(steps[Math.min(steps.length - 1, steps.indexOf(currentStep) + 1)])} disabled={files.length === 0 || !groupingsFile}>Next</button>
        <span class="tooltiptext">Upload files to continue</span>
      </div>
    </div>

    {#if currentStep === 'Model Perturbation' && files.length > 0 && groupingsFile && selectedMethod}
      <button on:click={handleSubmit}>Submit</button>
    {/if}
  </div>
</div>
