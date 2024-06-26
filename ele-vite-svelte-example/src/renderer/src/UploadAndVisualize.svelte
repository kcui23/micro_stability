<script>
  import { onMount } from 'svelte';
  import ADGPlot from './ADGPlot.svelte';

  let asvFiles = [];
  let groupingsFile = null;
  let visualizations = {
    deseq2_plot1: '', deseq2_plot2: '', deseq2_plot3: '',
    aldex2_plot1: '', aldex2_plot2: '', aldex2_plot3: '',
    edger_plot1: '', edger_plot2: '', edger_plot3: '',
    maaslin2_plot1: '', maaslin2_plot2: '', maaslin2_plot3: '',
    overlap_volcano: '', overlap_pvalue_distribution: ''
  };
  let methodFileStatus = {
    deseq2: false,
    aldex2: false,
    edger: false,
    maaslin2: false
  };
  let previousMethodFileStatus = { ...methodFileStatus };
  let dataChanged = false;
  let selectedMethod = '';
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
  let isCalculatingMissing = false;
  let randomSeed = 1234;
  let edgeThicknesses = [1,2,3,4];
  let isSubmitted = false;

  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Prediction Evaluation Metric', 'Stability Metric'];

  const showNotification = () => {
    const notification = document.getElementById('notification');
    notification.classList.add('show');
    setTimeout(() => {
      notification.classList.add('hide');
      setTimeout(() => {
        notification.classList.remove('show', 'hide');
      }, 500); // Match the duration of the hide transition
    }, 3000); // Duration to show the notification
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

  // Event handlers
  const handleFileChange = (event) => {
    const fileInput = event.target;
    const fileNameDisplay = document.getElementById('fileName1');
    asvFiles = fileInput.files;
    if (fileInput.files.length > 0) {
      const reader = new FileReader();
      reader.onload = () => {
        const content = reader.result;
        previewContent = previewFileContent(content);
        fileDimensions = getFileDimensions(content);
        filteredContent = previewContent;
        filteredDimensions = fileDimensions;
        filteredAsvContent = content;
      };
      fileNameDisplay.textContent = `Selected file: ${fileInput.files[0].name}`;
      reader.readAsText(fileInput.files[0]);
    } else {
      fileNameDisplay.textContent = '';
    }
    resetMethodStatus();
  };

  const handleGroupingsChange = (event) => {
    const fileInput = event.target;
    const fileNameDisplay = document.getElementById('fileName2');
    groupingsFile = fileInput.files[0];
    if (groupingsFile) {
      const reader = new FileReader();
      reader.onload = () => {
        const content = reader.result;
        groupingsContentPreview = previewFileContent(content);
        groupingsDimensions = getFileDimensions(content);
      };
      fileNameDisplay.textContent = `Selected file: ${fileInput.files[0].name}`;
      reader.readAsText(groupingsFile);
    } else {
      fileNameDisplay.textContent = '';
    }
    resetMethodStatus();
  };

  const handleMethodChange = (method) => {
    selectedMethod = method;
    showAllPlots = false;
    showDetailedPlots = false;
    isSubmitted = false;
  };

  const handleQuickExplore = async () => {
    showAllPlots = true;
    isCalculating = true;
    const asvContent = filteredAsvContent || asvFiles[0];

    if (typeof asvContent === 'string') {
      await processQuickExplore(asvContent, groupingsFile);
    } else {
      const asvReader = new FileReader();
      asvReader.onload = async () => {
        const asvContentText = asvReader.result;
        await processQuickExplore(asvContentText, groupingsFile);
      };
      asvReader.readAsText(asvContent);
    }
  };

  const handleSubmit = async () => {
    isCalculating = true;
    isSubmitted = true;
    const asvContent = filteredAsvContent || asvFiles[0];

    if (typeof asvContent === 'string') {
      await processSubmit(asvContent, groupingsFile);
      await checkMethodFileStatus();
     } else {
      const asvReader = new FileReader();
      asvReader.onload = async () => {
        const asvContentText = asvReader.result;
        await processSubmit(asvContentText, groupingsFile);
      };
      asvReader.readAsText(asvContent);
    }
  };

  // a reactive statement to check method file status when entering the Stability Metric step
  $: if (currentStep === 'Stability Metric')  {
    checkMethodFileStatus();
  }

  const handleFilter = async () => {
    const file = asvFiles[0];

    const asvReader = new FileReader();
    asvReader.onload = async () => {
      const asvContent = asvReader.result;

      try {
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

        // Parse the JSON response
        const result = await response.json();
        let filteredAsv = result.filteredAsv;

        // Ensure filteredAsv is a string
        if (typeof filteredAsv !== 'string') {
          filteredAsv = String(filteredAsv);
        }

        filteredContent = previewFileContent(filteredAsv);
        filteredDimensions = getFileDimensions(filteredAsv);
        filteredAsvContent = filteredAsv;
      } catch (error) {
        console.error('Error in filtering:', error);
      }
    };

    asvReader.readAsText(file);
    resetMethodStatus();
  };

  const resetMethodStatus = () => {
    previousMethodFileStatus = { ...methodFileStatus };
    methodFileStatus = {
      deseq2: [false],
      aldex2: [false],
      edger: [false],
      maaslin2: [false]
    };
    dataChanged = true;
  };

  const goToStep = (step) => {
    if (currentStep === 'Raw data' && (asvFiles.length === 0 || !groupingsFile)) {
      showNotification();
      return;
    }
    lastStep = currentStep;
    currentStep = step;

    const contentContainer = document.querySelector('.content');
    if (contentContainer) {
      contentContainer.scrollTo(0, 0);
    }
  };

  const calculateMissingMethods = async () => {
    isCalculatingMissing = true;
    const methodsToCalculate = dataChanged ? 
      Object.keys(methodFileStatus) : 
      Object.entries(methodFileStatus).filter(([_, status]) => !status[0]).map(([method, _]) => method);

    for (const method of methodsToCalculate) {
      try {
        const response = await fetch(`http://localhost:8000/process?method=${method}`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            asv: filteredAsvContent,
            groupings: await (new Response(groupingsFile)).text(),
            threshold: threshold,
            seed: randomSeed
          })
        });

        if (!response.ok) {
          throw new Error(`Failed to process ${method}`);
        }

        methodFileStatus[method] = [true];
      } catch (error) {
        console.error(`Error processing ${method}:`, error);
      }
    }

    isCalculatingMissing = false;
    dataChanged = false;
    await checkMethodFileStatus();
  };

  // API calls
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
          edger_plot1: `data:image/png;base64,${result.edgeR_plot1}`,
          edger_plot2: `data:image/png;base64,${result.edgeR_plot2}`,
          edger_plot3: `data:image/png;base64,${result.edgeR_plot3}`,
          maaslin2_plot1: `data:image/png;base64,${result.maaslin2_plot1}`,
          maaslin2_plot2: `data:image/png;base64,${result.maaslin2_plot2}`,
          maaslin2_plot3: `data:image/png;base64,${result.maaslin2_plot3}`,
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
          ...visualizations,
          [`${selectedMethod}_plot1`]: `data:image/png;base64,${result.plot1}`,
          [`${selectedMethod}_plot2`]: `data:image/png;base64,${result.plot2}`,
          [`${selectedMethod}_plot3`]: `data:image/png;base64,${result.plot3}`
        };
        isCalculating = false;
      } catch (error) {
        console.error('Fetch error:', error);
        isCalculating = false;
      }
    };

    groupingsReader.readAsText(groupings);
  };

  const handleDownload = async () => {
    try {
      const response = await fetch(`http://localhost:8000/download?method=${selectedMethod}`, {
        method: 'GET'
      });

      if (!response.ok) {
        throw new Error('Network response was not ok');
      }

      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = `${selectedMethod}_results.tsv`;
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error('Fetch error:', error);
    }
  };

  async function autoLoadFiles() {
    const asvPath = '../../Blueberry/Blueberry_ASVs_table.tsv';
    const groupingPath = '../../Blueberry/Blueberry_metadata.tsv';

    try {
      const asvResponse = await fetch(asvPath);
      const groupingResponse = await fetch(groupingPath);

      if (asvResponse.ok && groupingResponse.ok) {
        const asvFile = new File([await asvResponse.blob()], 'Blueberry_ASVs_table.tsv');
        const groupingFile = new File([await groupingResponse.blob()], 'Blueberry_metadata.tsv');

        handleFileChange({ target: { files: [asvFile] } });
        handleGroupingsChange({ target: { files: [groupingFile] } });

        console.log('Debug: Files auto-loaded successfully');
      } else {
        console.log('Debug: One or both files not found');
      }
    } catch (error) {
      console.error('Debug: Error auto-loading files:', error);
    }
  }

  const checkMethodFileStatus = async () => {
    try {
      const response = await fetch('http://localhost:8000/check_method_files');
      if (response.ok) {
        methodFileStatus = await response.json();
        console.log("Updated method file status:", methodFileStatus);  // Debug log
      } else {
        console.error('Failed to check method file status');
      }
    } catch (error) {
      console.error('Error checking method file status:', error);
    }
  };

  const downloadMethodFile = async (method) => {
    if (!methodFileStatus[method]) {
      console.error(`File for ${method} is not ready for download`);
      return;
    }
    try {
      const response = await fetch(`http://localhost:8000/download_method_file?method=${method}`, {
        method: 'GET'
      });

      if (response.ok) {
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.style.display = 'none';
        a.href = url;
        a.download = `${method}_results.tsv`;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
      } else {
        console.error(`Failed to download ${method} file`);
      }
    } catch (error) {
      console.error(`Error downloading ${method} file:`, error);
    }
  };

  onMount(() => {
    autoLoadFiles();
  });
</script>

<style>
  /* .steps {
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
  } */
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

  .navigation .previous-button {
    margin-right: auto;
  }

  .navigation .next-button {
    margin-left: auto;
  }

  .method-file-status {
    margin-top: 20px;
  }

  .method-status {
    display: flex;
    align-items: center;
    margin-bottom: 10px;
  }

  .method-status span {
    margin-right: 10px;
  }

  .method-status button {
    padding: 5px 10px;
    background-color: #4CAF50;
    color: white;
    border: none;
    cursor: pointer;
  }

  .method-status button:disabled {
    background-color: #cccccc;
    cursor: not-allowed;
  }

  .method-status button:hover {
    background-color: #45a049;
  }

  .method-status button:disabled:hover {
    background-color: #919191;
  }

  .stability-vis button {
    margin-top: 10px;
    padding: 5px 10px;
    background-color: #4CAF50;
    color: white;
    border: none;
    cursor: pointer;
  }

  .stability-vis button:disabled {
    background-color: #cccccc;
    cursor: not-allowed;
  }
</style>

<div id="app" class="container">
  <!-- Non-modal notification -->
  <div id="notification" class="notification">
    <span class="notification-icon">⚠️</span>
    <span>Please upload the required files first.</span>
  </div>

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

    <div hidden={currentStep !== 'Raw data'}>
      <h1>Raw Data</h1>
      <div class="upload-container">
        <div class="upload-section">
          <span class="file-label">Upload ASV File:</span>
          <div class="custom-file-input">
            <label for="fileInput1">Choose File</label>
            <input id="fileInput1" type="file" accept=".tsv" on:change={handleFileChange} />
          </div>
          <div class="note">Note: Please upload a .tsv file</div>
          <div id="fileName1" class="file-name"></div>
        </div>
  
        <div class="upload-section">
          <span class="file-label">Upload Groupings File:</span>
          <div class="custom-file-input">
            <label for="fileInput2">Choose File</label>
            <input id="fileInput2" type="file" accept=".tsv" on:change={handleGroupingsChange} />
          </div>
          <div class="note">Note: Please upload a .tsv file</div>
          <div id="fileName2" class="file-name"></div>
        </div>
      </div>

      
      <div>
        <label for="randomSeed">Set Random Seed:</label>
        <input type="number" id="randomSeed" bind:value={randomSeed} min="1" />
      </div>
    </div>

    <div class="preview-section" hidden={currentStep !== 'Raw data' && currentStep !== 'Data Perturbation'}>
      {#if asvFiles.length > 0}
        <div class="preview">
          <h2>Preview of ASV File</h2>
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
    </div>

    <div class="quick-explore-section" hidden={currentStep !== 'Raw data'}>
      <div class="tooltip">
        <button on:click={handleQuickExplore} disabled={asvFiles.length === 0 || !groupingsFile}>Quick Explore</button>
        <span class="tooltiptext">Upload files to continue</span>
      </div>
    </div>

    <div class="filter-section" hidden={currentStep !== 'Data Perturbation'}>
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
        <button on:click={() => handleMethodChange('edger')} class:selected={selectedMethod === 'edger'}>Method 3 (edgeR)</button>
        <button on:click={() => handleMethodChange('maaslin2')} class:selected={selectedMethod === 'maaslin2'}>Method 4 (Maaslin2)</button>
        <button on:click={() => handleMethodChange('method5')} class:selected={selectedMethod === 'method5'}>Method 5</button>
      </div>
    </div>

    <div class="stability-vis" hidden={currentStep !== 'Stability Metric'}>
      <h2>Stability Metric</h2>
      <div class="method-file-status">
        <h3>Method File Status:</h3>
        {#each Object.entries(methodFileStatus) as [method, status]}
          <div class="method-status">
            <span>{method}: {status[0] ? '✅ Ready' : '❌ Not ready'}</span>
            {#if status[0]}
              <button on:click={() => downloadMethodFile(method)}>Download {method} file</button>
            {:else}
              <button disabled>Download {method} file</button>
            {/if}
            {#if dataChanged && previousMethodFileStatus[method][0]}
              <span class="warning">⚠️ Data changed, recalculation needed</span>
            {/if}
          </div>
        {/each}
      </div>
      
      <button on:click={calculateMissingMethods} disabled={isCalculatingMissing}>
        {#if isCalculatingMissing}
          Calculating...
        {:else if dataChanged}
          Recalculate All Methods
        {:else if Object.values(methodFileStatus).some(status => !status[0])}
          Calculate Missing Methods
        {:else}
          All Methods Calculated
        {/if}
      </button>
    </div>    

    <div class="visualizations-section" hidden={!showAllPlots && !isSubmitted}>
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
            <button on:click={() => showDetailedPlots = true}>Show More Details</button>
          {:else}
            <button on:click={() => showDetailedPlots = false}>Collapse</button>
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
            <h3>edgeR Plots</h3>
            <img src={visualizations.edger_plot1} alt="edgeR Plot 1" style="width: 300px; height: auto;" />
            <img src={visualizations.edger_plot2} alt="edgeR Plot 2" style="width: 300px; height: auto;" />
            <img src={visualizations.edger_plot3} alt="edgeR Plot 3" style="width: 300px; height: auto;" />
            <h3>Maaslin2 Plots</h3>
            <img src={visualizations.maaslin2_plot1} alt="Maaslin2 Plot 1" style="width: 300px; height: auto;" />
            <img src={visualizations.maaslin2_plot2} alt="Maaslin2 Plot 2" style="width: 300px; height: auto;" />
            <img src={visualizations.maaslin2_plot3} alt="Maaslin2 Plot 3" style="width: 300px; height: auto;" />
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
        {:else if selectedMethod === 'edger'}
          <h3>edgeR Plots</h3>
          <img src={visualizations.edger_plot1} alt="edgeR Plot 1" style="width: 300px; height: auto;" />
          <img src={visualizations.edger_plot2} alt="edgeR Plot 2" style="width: 300px; height: auto;" />
          <img src={visualizations.edger_plot3} alt="edgeR Plot 3" style="width: 300px; height: auto;" />
        {:else if selectedMethod === 'maaslin2'}
          <h3>Maaslin2 Plots</h3>
          <img src={visualizations.maaslin2_plot1} alt="Maaslin2 Plot 1" style="width: 300px; height: auto;" />
          <img src={visualizations.maaslin2_plot2} alt="Maaslin2 Plot 2" style="width: 300px; height: auto;" />
          <img src={visualizations.maaslin2_plot3} alt="Maaslin2 Plot 3" style="width: 300px; height: auto;" />
        {/if}
      {/if}
    </div>
    
    {#if currentStep === 'Model Perturbation' && asvFiles.length > 0 && groupingsFile && selectedMethod}
      <button on:click={handleSubmit}>Submit</button>
      <button on:click={handleDownload} disabled={!selectedMethod || !isSubmitted}>Download</button>
    {/if}

    <div class="navigation">
      {#if currentStep != 'Raw data'}
        <div class="previous-button">
          <button on:click={() => goToStep(steps[Math.max(0, steps.indexOf(currentStep) - 1)])}>Previous</button>
        </div>
      {/if}
      {#if currentStep != 'Stability Metric'}
        <div class="tooltip next-button">
          <button on:click={() => goToStep(steps[Math.min(steps.length - 1, steps.indexOf(currentStep) + 1)])} disabled={asvFiles.length === 0 || !groupingsFile}>Next</button>
          <span class="tooltiptext">Upload files to continue</span>
        </div>
      {/if}
    </div>

  </div>
</div>
