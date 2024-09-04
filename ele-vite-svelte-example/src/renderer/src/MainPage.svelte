<script>
  import { fade } from 'svelte/transition';
  import { onMount, onDestroy } from 'svelte';
  import { selectedPoints, 
    currentPath, 
    stepStatus, 
    selectedOperations, 
    singleSelectOperations 
  } from './store.js';

  import D3Tree from './components/D3Tree.svelte';
  import FileUploader from './components/FileUploader.svelte';
  import PreviewSection from './components/PreviewSection.svelte';
  import VisualizationSection from './components/VisualizationSection.svelte';
  import SidebarComponent from './components/SidebarComponent.svelte';
  import ASVSelector from './components/ASVSelector.svelte';
  import ScatterPlot from './components/ScatterPlot.svelte';


  let DataPerturbationMethods = ['deseq2', 'edger', 'maaslin2', 'aldex2', 'method5']
  let treeData;
  let d3TreeContainer;
  let d3TreeComponent;
  let data_points_updated_counter = 0;
  let selectedPointsList = [];
  $: selectedPoints.subscribe(value => {
    selectedPointsList = value;
    console.log("selectedPointsList updated:", selectedPointsList);
  });

  let asvFiles = [];
  let groupingsFile = null;
  let isStatic = true;
  let zoomedImage = null;
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
  let selectedMethod = 'deseq2';
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
  let isCombiningResults = false;
  let combinedResultsReady = false;
  let randomSeed = 1234;
  let isSubmitted = false;
  let showASVSelector = false;
  let stabilityPlot = '';
  let shuffledAnalysisProgress = 0;
  let shuffledAnalysisTotal = 10;
  let shuffledAnalysisPlot = '';
  let isShuffledAnalysisRunning = false;
  let ws;
  let ws_id;

  const interactWithJson = async () => {
    console.log("Interacting with JSON in main page");

    try {
      const response = await fetch(`http://localhost:8000/update_leaf_data`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
      });

      if (response.ok) {
        data_points_updated_counter += 1;
        console.log("Data points updated successfully in main page.");
        const jsonResponse = await fetch('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json');
        if (jsonResponse.ok) {
          const jsonData = await jsonResponse.json();
          console.log("Local data points JSON file content:", jsonData);
        } else {
          console.error("Not able to fetch local data points JSON file");
        }
      } else {
        console.error("Failed to update data points JSON file");
      }
    } catch (error) {
      console.error("Error during fetch:", error);
    }
  }

  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Stability Metric'];

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
      reader.readAsText(fileInput.files[0]);
    }
    resetMethodStatus();
  };

  const handleGroupingsChange = (event) => {
    const fileInput = event.target;
    groupingsFile = fileInput.files[0];
    if (groupingsFile) {
      const reader = new FileReader();
      reader.onload = () => {
        const content = reader.result;
        groupingsContentPreview = previewFileContent(content);
        groupingsDimensions = getFileDimensions(content);
      };
      reader.readAsText(groupingsFile);
    }
    resetMethodStatus();
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
    showAllPlots = false;
    showDetailedPlots = false;
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
    const asvPath = '../../datasets/Blueberry/Blueberry_ASVs_table.tsv';
    const groupingPath = '../../datasets/Blueberry/Blueberry_metadata.tsv';

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

  const generateCombinedResults = async () => {
  isCombiningResults = true;
  try {
    const response = await fetch('http://localhost:8000/generate_combined_results', {
      method: 'POST'
    });

    if (response.ok) {
      const result = await response.json();
      console.log(result.message);
      combinedResultsReady = true;
      await fetchStabilityPlot();  // Fetch the stability plot after generating combined results
    } else {
      console.error('Failed to generate combined results');
    }
  } catch (error) {
    console.error('Error generating combined results:', error);
  }
  isCombiningResults = false;
};

const downloadCombinedResults = async () => {
  try {
    const response = await fetch('http://localhost:8000/download_combined_results', {
      method: 'GET'
    });

    if (response.ok) {
      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = 'combined_results.tsv';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } else {
      console.error('Failed to download combined results');
    }
  } catch (error) {
    console.error('Error downloading combined results:', error);
  }
};

const fetchStabilityPlot = async () => {
  try {
    const response = await fetch('http://localhost:8000/generate_stability_plot');
    if (response.ok) {
      const blob = await response.blob();
      stabilityPlot = URL.createObjectURL(blob);
    } else {
      console.error('Failed to fetch stability plot');
    }
  } catch (error) {
    console.error('Error fetching stability plot:', error);
  }
};

const runShuffledAnalysis = async () => {
    if (!ws_id) {
      console.error('WebSocket not connected');
      return;
    }

    isShuffledAnalysisRunning = true;
    shuffledAnalysisProgress = 0;

    try {
      const response = await fetch('http://localhost:8000/shuffled_analysis', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          asv: filteredAsvContent,
          groupings: await (new Response(groupingsFile)).text(),
          seed: randomSeed,
          iterations: shuffledAnalysisTotal,
          ws_id: ws_id
        })
      });

      if (!response.ok) {
        throw new Error('Failed to run shuffled analysis');
      }

      const result = await response.json();
      shuffledAnalysisPlot = `data:image/png;base64,${result.plot}`;
    } catch (error) {
      console.error('Error running shuffled analysis:', error);
    } finally {
      isShuffledAnalysisRunning = false;
    }
  };

  function handleScatterPointClick(event) {
    const { path } = event.detail;

    // Update d3 tree
    currentPath.set(path);
    if (d3TreeComponent) {
      d3TreeComponent.highlightPath(path);
    } else {
      console.error("D3Tree component not found");
    }

    // Update sidebar
    stepStatus.update(status => {
      return Object.fromEntries(
        Object.keys(status).map(key => [key, 'Enabled'])
      );
    });
    
    selectedOperations.update(operations => {
      operations[path[1]] = [path[2]];
      operations['Model Perturbation'] = ['Select Method'];
      return operations;
    });

    // deal with single select and multi-select in 'Stability Metric' step
    selectedOperations.update((selections) => {
			const isSingleSelect = $singleSelectOperations[path[5]] && $singleSelectOperations[path[5]].includes(path[6]);
			if (isSingleSelect) {
				const existingSingleSelect = $singleSelectOperations[path[5]].find(op => selections[path[5]].includes(op));
				console.log('existingSingleSelect:', existingSingleSelect);
				if (existingSingleSelect) {
					selections[path[5]] = selections[path[5]].filter(op => op !== existingSingleSelect);
				}
				selections[path[5]].push(path[6]);
			} else {
				if (selections[path[5]].includes(path[6])) {
					selections[path[5]] = selections[path[5]].filter(op => op !== path[6]);
				} else {
					selections[path[5]] = [...selections[path[5]], path[6]];
				}
			}

		return selections;
	});

    // Update selected methods
    selectedMethod = path[4];
  }

  function toggleView(view) {
    isStatic = view === 'static';
  }

  function zoomImage(image) {
    zoomedImage = zoomedImage === image ? null : image;
  }

  async function downloadImage(image, method) {
    try {
      const response = await fetch(`http://localhost:8000/download_image?method=${method}&plot=${image}`);
      if (response.ok) {
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.style.display = 'none';
        a.href = url;
        a.download = `${method}_${image}.png`;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
      } else {
        const errorText = await response.text();
        console.error('Failed to download image:', errorText);
      }
    } catch (error) {
      console.error('Error downloading image:', error);
    }
  }

  function handleKeydown(event) {
    if (event.key === 'Escape' && zoomedImage) {
      zoomImage(null);
    }
  }

  $: {
    console.log('Current step:', currentStep);
    console.log('Selected operations:', JSON.stringify($selectedOperations, null, 2));
    console.log("Tree data in main:", treeData);
  }

  onMount(() => {
    autoLoadFiles();
    window.addEventListener('keydown', handleKeydown);

    ws = new WebSocket('ws://localhost:8000/ws');
    ws.onopen = () => {
      console.log('WebSocket connected');
      ws_id = ws.url.split('/').pop();
    };
    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      if (data.progress && data.total) {
        shuffledAnalysisProgress = data.progress;
        shuffledAnalysisTotal = data.total;
      }
    };

    fetch('https://raw.githubusercontent.com/kcui23/micro_stability/main/ele-vite-svelte-example/src/renderer/src/public/data.json')
      .then(response => response.json())
      .then(data => {
        console.log("Data loaded successfully:", data);
        treeData = data
      })
      .catch(error => console.error('Error loading or parsing data.json:', error));


  });

  onDestroy(() => {
    if (ws) ws.close();
    window.removeEventListener('keydown', handleKeydown);

    d3.select(d3TreeContainer).select("svg").remove(); // Clean up on component destruction
  });
</script>

<style>
  .sidebar {
    width: 300px;
    padding: 20px;
    background-color: #f5f5f5;
    height: 100vh;
    overflow-y: auto;
  }

  .filters {
    margin: 20px 0;
  }
  .methods {
    display: flex;
    flex-direction: column;
    margin: 20px 0;
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

#asv-selector-main {
  margin-top: 10px;
}
.progress-bar {
    width: 100%;
    height: 20px;
    background-color: #f0f0f0;
    border-radius: 10px;
    overflow: hidden;
  }
  .progress-bar-fill {
    height: 100%;
    background-color: #4CAF50;
    transition: width 0.5s ease-in-out;
  }

  .input-group {
    display: flex;
    align-items: center;
    margin-bottom: 10px;
  }
  .input-group label {
    margin-right: 10px;
  }
  .input-group input {
    width: 60px;
  }

  img {
    max-width: 100%;
    height: auto;
    margin-bottom: 1rem;
  }

  .tree-container-wrapper {
    display: flex;
    width: 100%;
    height: 200px;
    margin-bottom: 20px;
    padding-right: 10px;
  }
  .scatter-plot-container {
    flex: 0 0 20%;
    height: 100%;
    background-color: #fff;
    border: 1px solid #ddd;
    border-radius: 4px;
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    box-sizing: border-box;
    border-radius: 8px;
    margin-left: 10px;
  }

</style>

<div id="app" class="container">
  <!-- Non-modal notification -->
  <div id="notification" class="notification">
    <span class="notification-icon">⚠️</span>
    <span>Please upload the required files first.</span>
  </div>

  <!-- Sidebar -->
  <div class="sidebar">
    <SidebarComponent 
      steps={steps} 
      currentStep={currentStep} 
      setCurrentStep={goToStep}
    />
  </div>

  <!-- Main Content Area -->
  <div class="content">
    <div class="tree-container-wrapper">
      <!-- D3 Tree Container -->
      <D3Tree 
        {treeData}
        setCurrentStep={goToStep}
        bind:this={d3TreeComponent} 
      />
      <!-- Scatter Plot -->
      <div class="scatter-plot-container">
        <ScatterPlot on:pointClick={handleScatterPointClick} {data_points_updated_counter}/>
      </div> 
    </div>

    <div class="logo">
      <h1>Logo</h1>
    </div>


    {#if currentStep === 'Raw data'}
      <div key='raw-data' in:fade class="step-content" class:active={currentStep === 'Raw data'}>
        <h2>Raw Data</h2>
        <!-- ASV File Upload UI -->
        <FileUploader 
          {handleFileChange} 
          {handleGroupingsChange} 
          bind:asvFiles
          bind:groupingsFile
        />
        
        {#if $selectedOperations['Raw data']?.includes('Preview')}
          <PreviewSection
            {filteredContent}
            {groupingsContentPreview}
            {filteredDimensions}
            {groupingsDimensions}
            {asvFiles}
            {groupingsFile}
          />
        {/if}

        {#if $selectedOperations['Raw data']?.includes('Set Random Seed')}
          <!-- Random Seed Input UI -->
          <div>
            <label for="randomSeed">Set Random Seed:</label>
            <input type="number" id="randomSeed" bind:value={randomSeed} min="1" />
          </div>
        {/if}

        {#if $selectedOperations['Raw data']?.includes('Quick Explore')}
          <div class="tooltip">
            <button on:click={handleQuickExplore} disabled={asvFiles.length === 0 || !groupingsFile}>Quick Explore</button>
            <span class="tooltiptext">Upload files to continue</span>
          </div>
        {/if}

      </div>
    {:else if currentStep === 'Data Perturbation'}
      <div key='data-perturbation' in:fade class="step-content" class:active={currentStep === 'Data Perturbation'}>
        <h2>Data Perturbation</h2>
        <p>Dimensions: {filteredDimensions.rows} rows, {filteredDimensions.columns} columns</p>
        {#if $selectedOperations['Data Perturbation']?.includes('Filter')}
          <div class="filters">
            <label for="filter">Filter:</label>
            <select id="filter">
              <option value="genefilter">genefilter</option>
              <option value="koverAfilter">koverAfilter</option>
            </select>
            <button>Apply Filter</button>
          </div>
        {/if}
        {#if $selectedOperations['Data Perturbation']?.includes('Threshold')}
          <!-- Threshold Application UI -->
            <div class="filters">
              <label for="threshold">Threshold for Rare Genes:</label>
              <input type="range" id="threshold" bind:value={threshold} min="0" max="1" step="0.01" />
              <span>{threshold}</span>
              <button on:click={handleFilter}>Apply Threshold</button>
            </div>
        {/if}
        {#if $selectedOperations['Data Perturbation']?.includes('Transformation')}
          <!-- Transformation Application UI -->
          <div class="filters">
            <label for="transformation">Transformation:</label>
            <select id="transformation">
              <option value="log">tidybult</option>
              <option value="sqrt">log</option>
              <option value="boxcox">normalization</option>
            </select>
            <button>Apply Transformation</button>
          </div>
        {/if}
        {#if $selectedOperations['Data Perturbation']?.includes('R/A Abundance')}
          <!-- R/A Abundance Application UI -->
          <div class="filters">
            <label for="raAbundance">R/A Abundance:</label>
            <select id="raAbundance">
              <option value="Compositional data analysis (CoDa)">Compositional data analysis (CoDa)</option>
              <option value="Additive Log-ratio Transformation">Additive Log-ratio Transformation</option>
            </select>
            <button>Apply R/A Abundance</button>
          </div>
        {/if}
        {#if $selectedOperations['Data Perturbation']?.includes('Data Splitting')}
          <!-- Data Splitting Application UI -->
          <div class="filters">
            <label for="dataSplitting">Data Splitting:</label>
            <select id="dataSplitting">
              <option value="splitByGroup">K-Fold</option>
            </select>
            <button>Apply Data Splitting</button>
          </div>
        {/if}
        {#if $selectedOperations['Data Perturbation']?.includes('Batch Effect Removal')}
          <!-- Batch Effect Removal Application UI -->
          <div class="filters">
            <label for="batchEffectRemoval">Batch Effect Removal:</label>
            <select id="batchEffectRemoval">
              <option value="ComBat">ComBat</option>
              <option value="RUV">RUV</option>
              <option value="CQN">CQN</option>
            </select>
            <button>Apply Batch Effect Removal</button>
          </div>
        {/if}
      </div>
    {:else if currentStep === 'Model Perturbation'}
      <div key='model-perturbation' in:fade class="step-content" class:active={currentStep === 'Model Perturbation'}>
        <h2>Model Perturbation</h2>
        {#if $selectedOperations['Model Perturbation']?.includes('Select Method')}
          <!-- Method Selection UI -->
          <div class="methods">
            <label for="method-select">Select Method:</label>
            <select id="method-select" bind:value={selectedMethod}>
              {#each DataPerturbationMethods as method}
                <option value={method}>{method}</option>
              {/each}
            </select>
          </div>
        {/if}

        {#if $selectedOperations['Model Perturbation']?.includes('Select Method') && asvFiles.length > 0 && groupingsFile && selectedMethod}
          <button on:click={handleSubmit}>Submit</button>
          <button on:click={handleDownload} disabled={!selectedMethod || !isSubmitted}>Download</button>
        {/if}

      </div>
    {:else if currentStep === 'Stability Metric'}
      <div key='stability-metric' in:fade class="step-content" class:active={currentStep === 'Stability Metric'}>
        <h2>Stability Metric</h2>

        {#if $selectedOperations['Stability Metric']?.includes('Differences in ASVs')}
          <div>
            <button>Calculate Differences in ASVs</button>
          </div>
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('AUROC')}
          <div>
            <button>Calculate AUROC</button>
          </div>
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('FDR')}
          <div>
            <button>Calculate FDR</button>
          </div>
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('All methods calculation')}
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
    
          <div class="combined-results-actions">
            <button on:click={generateCombinedResults} disabled={isCombiningResults || !Object.values(methodFileStatus).every(status => status[0])}>
              {isCombiningResults ? 'Generating...' : 'Generate Combined Results'}
            </button>
            <button on:click={downloadCombinedResults} disabled={!combinedResultsReady}>
              Download Combined Results
            </button>
          </div>
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('View Stability Plot')}
          <!-- Stability Plot Viewing UI -->
          {#if stabilityPlot}
            <div class="stability-plot">
              <h3>Stability Plot</h3>
              <img src={stabilityPlot} alt="Stability Plot" style="width: 100%; max-width: 800px; height: auto;" />
            </div>
          {/if}
        {/if}
        {#if $selectedOperations['Stability Metric']?.includes('Run Shuffled Analysis')}
          <!-- Shuffled Analysis UI -->
          <h3>Shuffled Analysis</h3>

          <div class="input-group">
            <label for="iterations">Number of Iterations:</label>
            <input 
              type="number" 
              id="iterations" 
              bind:value={shuffledAnalysisTotal} 
              min="1" 
              max="1000"
              disabled={isShuffledAnalysisRunning}
            >
          </div>
  
          <button on:click={runShuffledAnalysis} disabled={isShuffledAnalysisRunning}>
            {isShuffledAnalysisRunning ? 'Running...' : 'Run Shuffled Analysis'}
          </button>
      
          {#if isShuffledAnalysisRunning}
            <div class="progress-bar">
              <div class="progress-bar-fill" style="width: {(shuffledAnalysisProgress / shuffledAnalysisTotal) * 100}%"></div>
            </div>
            <p>Progress: {shuffledAnalysisProgress} / {shuffledAnalysisTotal}</p>
          {/if}
      
          {#if shuffledAnalysisPlot}
            <div class="shuffled-analysis-plot">
              <h4>Stability of Significant ASVs Across Shuffled Analyses</h4>
              <img src={shuffledAnalysisPlot} alt="Shuffled Analysis Plot" style="width: 100%; max-width: 800px; height: auto;" />
            </div>
          {/if}
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('ASV Selector')}
          <button on:click={() => showASVSelector = !showASVSelector} disabled={!combinedResultsReady}>
            {showASVSelector ? 'Hide' : 'Show'} ASV Selector
          </button>
          <div id="asv-selector-main">
            {#if showASVSelector}
              <ASVSelector />
            {/if}
          </div>
        {/if}

        {#if $selectedOperations['Stability Metric']?.includes('json interaction')}
          <button on:click={() => interactWithJson()}>Interact with JSON</button>
        {/if}
        
      </div>
    {/if}

    <div class="stability-vis" hidden={currentStep !== 'Stability Metric'}>

      <div class="shuffled-analysis">
        
      </div>
    </div>

    <VisualizationSection
      {visualizations} 
      {zoomImage}
      {isCalculating}
      {showAllPlots}
      bind:showDetailedPlots
      {isStatic}
      {toggleView}
      {selectedPointsList}
      {selectedMethod}
      {isSubmitted}
      {handleDownload}
      {zoomedImage}
    />

  </div>
</div>
