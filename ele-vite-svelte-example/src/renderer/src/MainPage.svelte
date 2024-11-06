<script>
  import { fade } from 'svelte/transition';
  import { onMount, onDestroy } from 'svelte';
  import { selectedPoints, 
    currentPath, 
    stepStatus, 
    selectedOperations,
    singleSelectOperations,
    params,
    StabilityMetricCalFinished
  } from './store.js';
  $: console.log("startMethod in MainPage.svelte:", startMethod);
  import D3Tree from './components/D3Tree.svelte';
  import FileUploader from './components/FileUploader.svelte';
  import PreviewSection from './components/PreviewSection.svelte';
  import VisualizationSection from './components/VisualizationSection.svelte';
  import SidebarComponent from './components/SidebarComponent.svelte';
  import ASVSelector from './components/ASVSelector.svelte';
  import ScatterPlot from './components/ScatterPlot.svelte';
  import ControlPanel from './components/ControlPanel.svelte';
  import StartPage from './components/StartPage.svelte';

  let startMethod = 'deseq2';
  let specific_interact = false;
  let scatterPlotClicked = false;
  let missingMethods = ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'];
  let DataPerturbationMethods = ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'];
  let abundance_threshold = 5;
  let prevalence_threshold = 10;
  let variance_threshold = 0.0;
  let pseudocount = 1;
  let knn = 5;
  let knn_bound = 50;
  let treeData;
  let d3TreeContainer;
  let d3TreeComponent;
  let data_points_updated_counter = 0;
  let preview_update_counter = 0;
  let highlight_point_path = [];
  let selectedPointsList = [];
  $: selectedPoints.subscribe(value => {
    selectedPointsList = value;
  });

  $: if (preview_update_counter) {
    console.log("preview_update_counter:", preview_update_counter);
    updatePreVars();
  }

  let asvPreview = [];
  let groupingsPreview = [];
  let asvDimensions = [];
  let groupingsDimensions = [];
  let oldAsvDimensions = [];
  let oldGroupingsDimensions = [];

  async function updatePreVars() {
    try {
      const response = await fetch('http://localhost:8000/preview_data');
      if (!response.ok) {
        throw new Error('Network response was not ok');
      }
      const previewData = await response.json();
      
      oldAsvDimensions = asvDimensions;
      oldGroupingsDimensions = groupingsDimensions;
      asvPreview = previewData.asv_preview;
      asvDimensions = previewData.asv_dimensions;
      groupingsPreview = previewData.groupings_preview;
      groupingsDimensions = previewData.groupings_dimensions;
      // console.log("asvPreview:", asvPreview);
      // console.log("asvDimensions:", asvDimensions);
      // console.log("groupingsPreview:", groupingsPreview);
      // console.log("groupingsDimensions:", groupingsDimensions);
    } catch (err) {
      console.error('Error fetching preview data:', err);
    }
  }

  $: console.log("selectedOperations:", $selectedOperations);

  const calculateStabilityMetric = async (method, missing_methods, destroy=false) => {
    try {
      const response = await fetch('http://localhost:8000/calculate_stability_metric', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          method: method,
          missing_methods: missing_methods,
          destroy: destroy
        })
      });

      if (response.ok) {
        data_points_updated_counter += 1;
        console.log("Data points updated counter in calculateStabilityMetric:", data_points_updated_counter);

        missingMethods = missingMethods.filter(m => m !== method);
        console.log("missingMethods:", missingMethods);
        console.log("length of missingMethods:", missingMethods.length);
        if (missingMethods.length === 0) {
          StabilityMetricCalFinished.set(true);
          console.log("StabilityMetricCalFinished set to true in calculateStabilityMetric");
        }
        console.log("StabilityMetricCalFinished:", StabilityMetricCalFinished);

        const result = await response.json();
        return result.message; // Return the success message from the API
      } else {
        const errorMessage = await response.json();
        console.error('Failed to calculate stability metric:', errorMessage.error);
        throw new Error(errorMessage.error);
      }
    } catch (error) {
      console.error('Error calculating stability metric:', error);
      throw error;
    }
  };

  let asvFiles = [];
  let groupingsFile = null;
  let isStatic = true;
  let zoomedImage = null;
  let zero_distribution_plot = null;
  let visualizations = {
    deseq2_plot1: '', deseq2_plot2: '', deseq2_plot3: '',
    aldex2_plot1: '', aldex2_plot2: '', aldex2_plot3: '',
    edger_plot1: '', edger_plot2: '', edger_plot3: '',
    maaslin2_plot1: '', maaslin2_plot2: '', maaslin2_plot3: '',
    metagenomeseq_plot1: '', metagenomeseq_plot2: '', metagenomeseq_plot3: '',
    overlap_volcano: '', overlap_pvalue_distribution: ''
  };
  let methodFileStatus = {
    deseq2: false,
    aldex2: false,
    edger: false,
    maaslin2: false,
    metagenomeseq: false
  };
  let previousMethodFileStatus = { ...methodFileStatus };
  let dataChanged = false;
  let selectedMethod = 'deseq2';
  $: {
    if ($selectedOperations['Model Perturbation'] && $selectedOperations['Model Perturbation'].length > 0) {
      selectedMethod = $selectedOperations['Model Perturbation'][0].toLowerCase();
      console.log("selectedMethod updated by /$selectedOperations['Model Perturbation'][0]:", selectedMethod);
    }
  }

  let currentStep = 'Raw data';
  let lastStep = 'Raw data';
  let previewContent = [];
  let groupingsContentPreview = [];
  let fileDimensions = { rows: 0, columns: 0 };
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

    try {
      const response = await fetch(`http://localhost:8000/update_leaf_data`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
      });

      if (response.ok) {
        data_points_updated_counter += 1;
        const jsonResponse = await fetch('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json');
        if (jsonResponse.ok) {
          const jsonData = await jsonResponse.json();
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

  const steps = ['Raw data', 'Filtering', 'Zero-Handling', 'Normalization', 'Transformation', 'Model Perturbation', 'Stability Metric'];

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
    specific_interact = false;
    await processQuickExplore();
  };

  const handleSubmit = async () => {
    isCalculating = true;
    isSubmitted = true;
    specific_interact = true;
    showAllPlots = false;
    showDetailedPlots = false;
    await processSubmit();
    await checkMethodFileStatus();
  };

  $: if (currentStep === 'Stability Metric')  {
    checkMethodFileStatus();
  }

  const handleFilter = async () => {
    let filter_method = $selectedOperations['Filtering'][0];
    let tmp_threshold;
    if (filter_method === 'Low Abundance Filtering') {
      tmp_threshold = abundance_threshold;
    } else if (filter_method === 'Prevalence Filtering') {
      tmp_threshold = prevalence_threshold;
    } else if (filter_method === 'Variance Filtering') {
      tmp_threshold = variance_threshold;
    } else {
      tmp_threshold = threshold;
    }
    try {
      const response = await fetch(`http://localhost:8000/filter`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          filter_method: filter_method,
          threshold: tmp_threshold
        })
      });

      const result = await response.json();
      if (result.success) {
        console.log('Filtering applied successfully');
        preview_update_counter += 1;
      } else {
        console.error('Filtering failed');
      }
    } catch (error) {
      console.error('Error in filtering:', error);
    }
  };

  $: if (currentStep === 'Zero-Handling') {
    fetchZeroDistributionPlot();
  }

  async function fetchZeroDistributionPlot() {
    try {
      const response = await fetch('http://localhost:8000/zero_distribution_plot');
      if (!response.ok) {
        throw new Error('Network response was not ok');
      }
      const blob = await response.blob();
      zero_distribution_plot = URL.createObjectURL(blob);
    } catch (error) {
      console.error('Error fetching zero distribution plot:', error);
    }
  }

  const handleZeroHandling = async () => {
    let zero_handling_method = $selectedOperations['Zero-Handling'][0];
    let zero_value;
    if (zero_handling_method === 'Pseudocount Addition') {
      zero_value = pseudocount;
    } else if (zero_handling_method === 'k-NN Imputation') {
      zero_value = [knn, knn_bound];
    }
    try {
      const response = await fetch(`http://localhost:8000/zero_handling`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        }, 
        body: JSON.stringify({
          zero_handling_method: zero_handling_method,
          param_value: zero_value
        })
      });

      const result = await response.json();
      if (result.success) {
        console.log('Zero-handling applied successfully'); 
        preview_update_counter += 1;
        fetchZeroDistributionPlot();
      } else {
        console.error('Zero-handling failed');
      }
    } catch (error) {
      console.error('Error in zero-handling:', error);
    }
  };
  
  const handleNormalization = async () => {
    let norm_method = $selectedOperations['Normalization'][0];
    try {
      const response = await fetch(`http://localhost:8000/normalization`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        }, 
        body: JSON.stringify({
          norm_method: norm_method
        })
      });

      const result = await response.json();
      if (result.success) {
        console.log('Normalization applied successfully'); 
        preview_update_counter += 1;
        fetchZeroDistributionPlot();
      } else {
        console.error('Normalization failed');
      }
    } catch (error) {
      console.error('Error in normalization:', error);
    }
  };

  const handleTransformation = async () => {
    let trans_method = $selectedOperations['Transformation'][0];
    try {
      const response = await fetch(`http://localhost:8000/transformation`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        }, 
        body: JSON.stringify({
          trans_method: trans_method
        })
      });

      const result = await response.json();
      if (result.success) {
        console.log('Transformation applied successfully'); 
        preview_update_counter += 1;
        fetchZeroDistributionPlot();
      } else {
        console.error('Transformation failed');
      }
    } catch (error) {
      console.error('Error in transformation:', error);
    }
  };

  const resetMethodStatus = () => {
    previousMethodFileStatus = { ...methodFileStatus };
    methodFileStatus = {
      deseq2: [false],
      aldex2: [false],
      edger: [false],
      maaslin2: [false],
      metagenomeseq: [false]
    };
    dataChanged = true;
  };

  const goToStep = (step) => {
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
  const processQuickExplore = async () => {
      try {
        const response = await fetch(`http://localhost:8000/quick_explore`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            seed: randomSeed
          })
        });

        if (!response.ok) {
          throw new Error('Network response was not ok in processQuickExplore()');
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
          metagenomeseq_plot1: `data:image/png;base64,${result.metagenomeseq_plot1}`,
          metagenomeseq_plot2: `data:image/png;base64,${result.metagenomeseq_plot2}`,
          metagenomeseq_plot3: `data:image/png;base64,${result.metagenomeseq_plot3}`,
          overlap_volcano: `data:image/png;base64,${result.overlap_volcano}`,
          overlap_pvalue_distribution: `data:image/png;base64,${result.overlap_pvalue_distribution}`
        };
        isCalculating = false;
      } catch (error) {
        console.error('Fetch error:', error);
        isCalculating = false;
      }
  };

  const processSubmit = async () => {
    try {
      const response = await fetch(`http://localhost:8000/process?method=${selectedMethod}`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          seed: randomSeed
        })
      });

      if (!response.ok) {
        throw new Error('Network response was not ok in processSubmit()');
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

  const checkMethodFileStatus = async () => {
    try {
      const response = await fetch('http://localhost:8000/check_method_files');
      if (response.ok) {
        methodFileStatus = await response.json();
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
    scatterPlotClicked = true;
    showAllPlots = false;
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
      operations['Filtering'] = [path[1]];
      operations['Zero-Handling'] = [path[2]];
      operations['Normalization'] = [path[3]];
      operations['Transformation'] = [path[4]];
      operations['Model Perturbation'] = [path[5]];
      return operations;
    });

    // deal with single select and multi-select in 'Stability Metric' step
    selectedOperations.update((selections) => {
      // these 6 and 7 are the 'Stability Metric' step
      const isSingleSelect = $singleSelectOperations['Stability Metric'] && $singleSelectOperations['Stability Metric'].includes(path[6]);
      if (isSingleSelect) {
        const existingSingleSelect = $singleSelectOperations['Stability Metric'].find(op => selections['Stability Metric'].includes(op));
        if (existingSingleSelect) {
          selections['Stability Metric'] = selections['Stability Metric'].filter(op => op !== existingSingleSelect);
        }
        selections['Stability Metric'].push(path[6]);
      } else {
        if (selections['Stability Metric'].includes(path[6])) {
          selections['Stability Metric'] = selections['Stability Metric'].filter(op => op !== path[6]);
        } else {
          selections['Stability Metric'] = [...selections['Stability Metric'], path[6]];
        }
      }

      return selections;
    });
  }

  function highlightPoint(path, parent_func='') {
    console.log("Called from " + parent_func + ": highlightPoint:", path);
    highlight_point_path = path;
  }

  function updateTreeandScatterplot() {
    currentPath.update(path => {
      path[5] = selectedMethod;
      return path;
    });
    handlePathChangeFromSidebar($currentPath, true);
    highlightPoint($currentPath, "updateTreeandScatterplot");
  }

  function handlePathChangeFromSidebar(event, direct=false) {
    let data;
    let existingSingleSelect;
    if (direct) { data = event } else { 
      data = event.detail 
      // steps[5] is the method perturbation, deseq2, ...
      existingSingleSelect = $singleSelectOperations[steps[5]].find(op => data[steps[5]].includes(op));
    }
    let allEnabled = Object.values($stepStatus).every(value => value === 'Enabled');;
    
    if (allEnabled) {
      let path = new Array(7).fill('');
      path[0] = steps[0]
      path[1] = data[steps[1]][0]
      // path[2] = data[steps[1]]&&data[steps[1]].length > 0 ? data[steps[1]][0] : data[2]
      path[2] = data[steps[2]][0]
      path[3] = data[steps[3]][0]
      path[4] = data[steps[4]][0]
      path[5] = data[steps[5]][0]
      path[6] = data[steps[6]][0]
      // Update d3 tree
      currentPath.set(path);
      if (d3TreeComponent) {
        d3TreeComponent.highlightPath(path);
      } else {
        console.error("D3Tree component not found");
      }
      highlightPoint(path, "handlePathChangeFromSidebar");
    }
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

  $: params.update(currentParams => {
    const updatedParams = {
      ...currentParams,
      'Low Abundance Filtering': abundance_threshold,
      'Prevalence Filtering': prevalence_threshold,
      'Variance Filtering': variance_threshold,
      'Pseudocount Addition': pseudocount,
      'knn': knn,
      'knn_bound': knn_bound,
      'Random Seed': randomSeed
    };
    console.log('Updated params:', updatedParams);
    return updatedParams;
  })

  async function start_generate_all_codes() {
    try {
      const response = await fetch('http://localhost:8000/generate_all_r_code', {
        method: 'POST'
      });

      if (!response.ok) {
        throw new Error('Failed to generate all codes');
      }
      const result = await response.json();
      console.log('All codes generated successfully:', result);
    } catch (error) {
      console.error('Error generating all codes:', error);
    }
  }

  onMount(() => {
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

    fetch('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json')
      .then(response => response.json())
      .then(data => {
        treeData = data
      })
      .catch(error => console.error('Error loading or parsing data.json:', error));

    start_generate_all_codes();


  });

  onDestroy(() => {
    if (ws) ws.close();
    window.removeEventListener('keydown', handleKeydown);

    d3.select(d3TreeContainer).select("svg").remove(); // Clean up on component destruction
  });

  // Added State for Start Page
  let showStartPage = true;
  $: if (!showStartPage) {
    selectedOperations.update((operations) => {
      operations['Raw data'].push('Preview');
      return operations;
    });
  }
  let tmp_missingMethods = [];
  // Function to Hide Start Page
  function startApp() {
    calculateStabilityMetric(1,1,true); // set all data points to 0,0
    tmp_missingMethods = missingMethods.filter(m => m !== startMethod);
    calculateStabilityMetric(startMethod, missingMethods); // start with the first method selected from start page
    // commented for test
    for (let method of tmp_missingMethods) {
      calculateStabilityMetric(method, missingMethods);
    }
  }
  // function startAppSkip() {
  //   showStartPage = false;
  // }
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
      width: 100%;
      height: 200px;
      margin-top: 10px;
      margin-bottom: 20px;
  }

    .control-panel {
      width: 100%;
      height: 200px;
      padding: 10px;
  }

    .scatter-plot-ctn {
      width: 100%;
      height: 400px;
      padding: 10px;
      margin-bottom: 20px;
      overflow: hidden;
    }

  .normalization-description {
    font-size: 1.1rem;
    color: #4a4a4a;
    line-height: 1.6;
    margin-top: 10px;
    margin-bottom: 10px;
    background-color: #f0f4f8;
    padding: 15px 20px;
    border-left: 4px solid #007BFF;
    border-radius: 4px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    transition: background-color 0.3s ease, border-left-color 0.3s ease;
  }

  .normalization-description:hover {
    background-color: #e6f0fa;
    border-left-color: #0056b3;
  }

  strong {
    font-weight: bold;
  }

  .fix-height-img-container {
    margin-bottom: 0px;
    padding-bottom: 0px;
  }

  .fix-height-img-container img {
  display: block;
  margin-bottom: 0;
}

.apply-button {
  width: 200px;
  height: 35px;
  margin-top: 20px;
  font-size: 1rem;
  font-weight: bold;
  display: block;
  margin-left: auto;
  margin-right: auto;
  }

  div[key='model-perturbation'] p {
  margin: 5px;
    font-size: 1.1rem;
  }

  div[key='model-perturbation'] button {
    margin-left: 5px;
  }

</style>

<div id="app" class="container">
  <!-- Start Page Overlay -->
  {#if showStartPage}
    <StartPage 
      {handleFileChange} 
      {handleGroupingsChange} 
      {updatePreVars}
      {DataPerturbationMethods}
      {startApp}
      bind:asvFiles
      bind:groupingsFile
      bind:showStartPage
      bind:missingMethods
    />
  {/if}

  <!-- Sidebar -->
  <div class="sidebar">
    <SidebarComponent 
      steps={steps} 
      currentStep={currentStep} 
      setCurrentStep={goToStep}
      on:pathChange={handlePathChangeFromSidebar}
    />
  </div>

  <!-- Main Content Area -->
  <div class="content-left">
    <div class="tree-container-wrapper">
      <!-- D3 Tree Container -->
      <D3Tree 
        {treeData}
        setCurrentStep={goToStep}
        bind:this={d3TreeComponent} 
      />
    </div>

    <div class="logo">
      <h1>Logo</h1>
    </div>


    {#if currentStep === 'Raw data'}
      <div key='raw-data' in:fade class="step-content" class:active={currentStep === 'Raw data'}>
        <h2>Raw Data</h2>
        
        {#if $selectedOperations['Raw data']?.includes('Preview')}
          <PreviewSection
            {asvPreview}
            {groupingsPreview}
            {asvDimensions}
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
    {:else if currentStep === 'Filtering'}
      <div key='filter-rare-units' in:fade class="step-content" class:active={currentStep === 'Filtering'}>
        <h2>Filtering</h2>
        <div class="tidy-table-container">
          <table class="tidy-table">
            <thead>
              <tr>
                <th>Dimension</th>
                <th>Rows</th>
                <th>Columns</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>Original</td>
                <td>{oldAsvDimensions[0]}</td>
                <td>{oldAsvDimensions[1]}</td>
              </tr>
              <tr>
                <td>Filtered</td>
                <td>{asvDimensions[0]}</td>
                <td>{asvDimensions[1]}</td>
              </tr>
            </tbody>
          </table>
        </div>
        {#if $selectedOperations['Filtering']?.includes('Low Abundance Filtering')}
          <!-- Threshold Application UI -->
            <div class="filters">
              <p class="normalization-description">
                <strong>Low Abundance Filtering:</strong> Set ASV abundance to zero if it is below the specified threshold.
              </p>
              <label for="abundance_threshold">Low Abundance Filtering Threshold:</label>
              <input type="range" id="abundance_threshold" bind:value={abundance_threshold} min="0" max="100" step="1" />
              <span>{abundance_threshold}</span>
              <button on:click={handleFilter}>Apply Threshold</button>
            </div>
        {/if}
        {#if $selectedOperations['Filtering']?.includes('Prevalence Filtering')}
          <!-- Threshold Application UI -->
            <div class="filters">
              <p class="normalization-description">
                <strong>Prevalence Filtering:</strong> Remove ASVs with a prevalence below the specified threshold.
              </p>
              <label for="prevalence_threshold">Prevalence Filtering Threshold (%):</label>
              <input type="range" id="prevalence_threshold" bind:value={prevalence_threshold} min="0" max="100" step="1" />
              <span>{prevalence_threshold}</span>
              <button on:click={handleFilter}>Apply Threshold</button>
            </div>
        {/if}
        {#if $selectedOperations['Filtering']?.includes('Variance Filtering')}
          <!-- Threshold Application UI -->
            <div class="filters">
              <p class="normalization-description">
                <strong>Variance Filtering:</strong> Remove ASVs with a variance below the specified threshold.
              </p>
              <label for="variance_threshold">Variance Filtering Threshold:</label>
              <input type="range" id="variance_threshold" bind:value={variance_threshold} min="0" max="1" step="0.01" />
              <span>{variance_threshold}</span>
              <button on:click={handleFilter}>Apply Threshold</button>
            </div>
        {/if}
      </div>

      {:else if currentStep === 'Zero-Handling'}
      <div key='zero-handling' in:fade class="step-content" class:active={currentStep === 'Zero-Handling'}>
        <h2>Zero-Handling</h2>
        {#if zero_distribution_plot}
          <div class="fix-height-img-container">
            <img src={zero_distribution_plot} alt="Zero Distribution Plot" />
          </div>
        {/if}
        {#if $selectedOperations['Zero-Handling']?.includes('Pseudocount Addition')}
          <div class="filters">
            <p class="normalization-description">
              <strong>Pseudocount Addition:</strong> Add a constant value to all values (including zeros).
            </p>
            <label for="pseudocount">Pseudocount Addition:</label>
            <input type="range" id="pseudocount" bind:value={pseudocount} min="0" max="10" step="0.1" />
            <span>{pseudocount}</span>
            <button on:click={handleZeroHandling}>Apply Pseudocount Addition</button>
          </div>
        {/if}
        {#if $selectedOperations['Zero-Handling']?.includes('k-NN Imputation')}
          <div class="filters">
            <p class="normalization-description">
              <strong>k-NN Imputation:</strong> Impute missing values using k-nearest neighbors. <br>
              <strong>Note:</strong> Only applied to ASVs with zero percentage less than {Math.min(knn_bound, 50)}%.
              {#if knn_bound > 50}
                For rows with more than 50% missing values, mean imputation is applied.
              {/if}
            </p>
            <label for="knn">k-NN Imputation number of neighbors:</label>
            <input type="range" id="knn" bind:value={knn} min="0" max="10" step="1" />
            <span>{knn}</span> <br>
            <label for="knn-bound">k-NN Imputation bound of missing values (%):</label>
            <input type="range" id="knn-bound" bind:value={knn_bound} min="0" max="100" step="1" />
            <span>{knn_bound}</span> <br>
            <button on:click={handleZeroHandling}>Apply k-NN Imputation</button>
          </div>
        {/if}
      </div>

    {:else if currentStep === 'Normalization'}
      <div key='normalization' in:fade class="step-content" class:active={currentStep === 'Normalization'}>
        <h2>Normalization</h2>
          {#if $selectedOperations['Normalization']?.includes('TSS')}
            <p class="normalization-description">
              <strong>TSS:</strong> Total Sum Scaling<br>
              <strong>Description:</strong> Each value is divided by the total sum of all values in its group, converting data into relative proportions. <br>
              <strong>Usually used when:</strong> Adjusting data to account for differences in total sums across groups or samples.
            </p>
            <button on:click={handleNormalization} class="apply-button">Apply TSS</button>
          {/if}
          {#if $selectedOperations['Normalization']?.includes('CSS')}
            <p class="normalization-description">
              <strong>CSS:</strong> Cumulative Sum Scaling<br>
              <strong>Description:</strong> Values are scaled using a factor derived from cumulative distributions, reducing the influence of large values. <br>
              <strong>Usually used when:</strong> Normalizing data with skewed distributions to lessen the impact of extreme values.
            </p>
            <button on:click={handleNormalization} class="apply-button">Apply CSS</button>
          {/if}
          {#if $selectedOperations['Normalization']?.includes('TMM')}
            <p class="normalization-description">
              <strong>TMM:</strong> Trimmed Mean of M-values<br>
              <strong>Description:</strong> Computes scaling factors by averaging log ratios after trimming extreme values, adjusting for sample differences. <br>
              <strong>Usually used when:</strong> Correcting data for composition differences while minimizing the effect of outliers.
            </p>
            <button on:click={handleNormalization} class="apply-button">Apply TMM</button>
          {/if}
          {#if $selectedOperations['Normalization']?.includes('CLR')}
            <p class="normalization-description">
              <strong>CLR:</strong> Centered Log-Ratio<br>
              <strong>Description:</strong> Transforms data by taking the logarithm of each value divided by the group's geometric mean, emphasizing relative differences. <br>
              <strong>Usually used when:</strong> Analyzing data where only relative differences are meaningful.
            </p>
            <button on:click={handleNormalization} class="apply-button">Apply CLR</button>
          {/if}
          {#if $selectedOperations['Normalization']?.includes('No Normalization')}
            <p class="normalization-description">
              <strong>No Normalization:</strong> Do not perform normalization.
            </p>
          {/if}
      </div>

    {:else if currentStep === 'Transformation'}
      <div key='transformation' in:fade class="step-content" class:active={currentStep === 'Transformation'}>
        <h2>Transformation</h2>
        
          {#if $selectedOperations['Transformation']?.includes('Log')}
            <p class="normalization-description">
              <strong>Log</strong> <br>
              <strong>Description:</strong> Applies a logarithmic scale to data to reduce skewness and stabilize variance. <br>
              <strong>Usually used when:</strong> Dealing with data that spans multiple orders of magnitude and has a skewed distribution.
            </p>
            <button on:click={handleTransformation} class="apply-button">Apply Log</button>
          {/if}
          {#if $selectedOperations['Transformation']?.includes('Logit')}
            <p class="normalization-description"> 
              <strong>Logit</strong> <br>
              <strong>Description:</strong> Converts proportion data to a log-odds scale, linearizing relationships. <br>
              <strong>Usually used when:</strong> Analyzing proportions or probabilities that are bounded between 0 and 1.
            </p>
            <button on:click={handleTransformation} class="apply-button">Apply Logit</button>
          {/if}
          {#if $selectedOperations['Transformation']?.includes('AST')}
            <p class="normalization-description">
              <strong>AST:</strong> Arcsine Transformation<br>
              <strong>Description:</strong> Applies the arcsine square root transformation to each proportion, stabilizing variance for proportional data. <br>
              <strong>Usually used when:</strong> Dealing with proportional data, especially when proportions are near 0 or 1.
            </p>
            <button on:click={handleTransformation} class="apply-button">Apply AST</button>
          {/if}
          {#if $selectedOperations['Transformation']?.includes('No Transformation')}
            <strong>No Transformation:</strong> Do not perform transformation.
          {/if}
      </div>

    {:else if currentStep === 'Model Perturbation'}
      <div key='model-perturbation' in:fade class="step-content" class:active={currentStep === 'Model Perturbation'}>
        <h2>Model Perturbation</h2>
        {#if $selectedOperations['Model Perturbation']?.includes('deseq2')}
          <p>DESeq2</p>
        {/if}
        {#if $selectedOperations['Model Perturbation']?.includes('edger')}
          <p>edgeR</p>
        {/if}
        {#if $selectedOperations['Model Perturbation']?.includes('maaslin2')}
          <p>MaAsLin2</p>
        {/if}
        {#if $selectedOperations['Model Perturbation']?.includes('aldex2')}
          <p>ALDEx2</p>
        {/if}
        {#if $selectedOperations['Model Perturbation']?.includes('metagenomeseq')}
          <p>metagenomeSeq</p>
        {/if}

        <!-- Method Selection UI -->
        <!-- <div class="methods">
          <label for="method-select">Select Method:</label>
          <select id="method-select" bind:value={selectedMethod} on:change={updateTreeandScatterplot}>
            {#each DataPerturbationMethods as method}
              <option value={method}>{method}</option>
            {/each}
          </select>
        </div> -->

        {#if $selectedOperations['Model Perturbation'].length > 0 && asvFiles.length > 0 && groupingsFile && selectedMethod}
          <button on:click={handleSubmit}>Submit</button>
          <button on:click={handleDownload} disabled={!selectedMethod || !isSubmitted || isCalculating}>Download</button>
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
      {downloadImage}
      {zoomedImage}
      {specific_interact}
      {scatterPlotClicked}
    />

  </div>

  <div class="content-right">
    <!-- control panel -->
    <div class="control-panel">
      <ControlPanel />
    </div>
    <!-- Scatter Plot -->
    <div class="scatter-plot-ctn">
      <ScatterPlot 
        on:pointClick={handleScatterPointClick} 
        {data_points_updated_counter}
        {highlight_point_path}
      />
    </div> 
  </div>
</div>
