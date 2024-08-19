<script>
  import { fade, scale } from 'svelte/transition';
  import { onMount, onDestroy } from 'svelte';
  import SidebarComponent from './SidebarComponent.svelte';
  import ASVSelector from './ASVSelector.svelte';
  import InteractiveValcano from './InteractiveValcano.svelte';
  import ScatterPlot from './ScatterPlot.svelte';
  import {selectedPoints} from './store.js';
  import * as d3 from 'd3';


  let treeData;
  let treeRoot;
  let d3TreeContainer;
  let selectedPointsList = [];
  $: selectedPoints.subscribe(value => {
    selectedPointsList = value;
    console.log("selectedPointsList updated:", selectedPointsList);
  });
  let isCardExpanded = true;

  const toggleCard = () => {
    isCardExpanded = !isCardExpanded;
  };

  const downloadSelectedPoints = async () => {
    const data = selectedPointsList.map(point => ({
      name: point.name,
      x: point.x,
      y: point.y
    }));

    const response = await fetch('http://localhost:8000/download_selected_points', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ points: data })
    });

    if (response.ok) {
      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = 'selected_points.tsv';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } else {
      console.error('Error downloading selected points:', response.statusText);
    }
  };
  let asvFiles = [];
  let groupingsFile = null;
  let selectedOperations = {};
  let isInteractive = false; // Whether the user is interacting with the visualization
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
  let isCombiningResults = false;
  let combinedResultsReady = false;
  let randomSeed = 1234;
  let edgeThicknesses = [1,2,3,4];
  let isSubmitted = false;
  let showASVSelector = false;
  let stabilityPlot = '';
  let shuffledAnalysisProgress = 0;
  let shuffledAnalysisTotal = 10;
  let shuffledAnalysisPlot = '';
  let isShuffledAnalysisRunning = false;
  let ws;
  let ws_id;


  const steps = ['Raw data', 'Data Perturbation', 'Model Perturbation', 'Stability Metric'];

  const removePoint = (pointToRemove) => {
    selectedPoints.update(points => points.filter(point => point !== pointToRemove));
  };
  
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

  function renderD3Tree(container, data) {
    const width = container.offsetWidth || 400;
    const marginTop = 10;
    const marginRight = 10;
    const marginBottom = 10;
    const marginLeft = 40;

    const root = d3.hierarchy(data);
    const dx = 25;
    const dy = (width - marginRight - marginLeft) / (1.5 + root.height);

    const tree = d3.tree().nodeSize([dx, dy]);
    const diagonal = d3.linkHorizontal().x(d => d.y).y(d => d.x);

    const svg = d3.select(container).append("svg")
      .attr("width", width)
      .attr("height", 300)
      .attr("viewBox", [-marginLeft, -marginTop, width, 300])
      .attr("style", "max-width: 100%; height: auto; font: 12px sans-serif; user-select: none;");

    const gLink = svg.append("g")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-opacity", 0.4)
      .attr("stroke-width", 1.5);

    const gNode = svg.append("g")
      .attr("cursor", "pointer")
      .attr("pointer-events", "all");

    function update(event, source) {
      const duration = event?.altKey ? 2500 : 250;
      const nodes = root.descendants().reverse();
      const links = root.links();

      tree(root);

      let left = root;
      let right = root;
      root.eachBefore(node => {
        if (node.x < left.x) left = node;
        if (node.x > right.x) right = node;
      });

      const height = right.x - left.x + marginTop + marginBottom + 20;

      const transition = svg.transition()
        .duration(duration)
        .attr("height", height)
        .attr("viewBox", [-marginLeft, left.x - marginTop, width, height]);

      const node = gNode.selectAll("g")
        .data(nodes, d => d.id);

      const nodeEnter = node.enter().append("g")
        .attr("transform", d => `translate(${source.y0},${source.x0})`)
        .attr("fill-opacity", 0)
        .attr("stroke-opacity", 0)
        .on("click", (event, d) => {
          d.children = d.children ? null : d._children;
          update(event, d);
        });

      nodeEnter.append("circle")
        .attr("r", 2.5)
        .attr("fill", d => d._children ? "#555" : "#999")
        .attr("stroke-width", 10);

      nodeEnter.append("text")
        .attr("dy", "1em")
        .attr("x", d => d._children ? -6 : 6)
        .attr("text-anchor", d => d._children ? "end" : "start")
        .text(d => d.data.name)
        .clone(true).lower()
        .attr("stroke-linejoin", "round")
        .attr("stroke-width", 3)
        .attr("stroke", "white");

      const nodeUpdate = node.merge(nodeEnter).transition(transition)
        .attr("transform", d => `translate(${d.y},${d.x})`)
        .attr("fill-opacity", 1)
        .attr("stroke-opacity", 1);

      const nodeExit = node.exit().transition(transition).remove()
        .attr("transform", d => `translate(${source.y},${source.x})`)
        .attr("fill-opacity", 0)
        .attr("stroke-opacity", 0);

      const link = gLink.selectAll("path")
        .data(links, d => d.target.id);

      const linkEnter = link.enter().append("path")
        .attr("d", d => {
          const o = { x: source.x0, y: source.y0 };
          return diagonal({ source: o, target: o });
        });

      link.merge(linkEnter).transition(transition)
        .attr("d", diagonal);

      link.exit().transition(transition).remove()
        .attr("d", d => {
          const o = { x: source.x, y: source.y };
          return diagonal({ source: o, target: o });
        });

      root.eachBefore(d => {
        d.x0 = d.x;
        d.y0 = d.y;
      });
    }

    root.x0 = dy / 2;
    root.y0 = 0;
    root.descendants().forEach((d, i) => {
      d.id = i;
      d._children = d.children;
      if (d.depth && d.data.name.length !== 7) d.children = null;
    });

    treeRoot = root;
    update(null, root);
  }

  function highlightPath(path) {
    // Remove previous highlights
    d3.select(d3TreeContainer).selectAll('.node').classed('highlighted', false);
    d3.select(d3TreeContainer).selectAll('.link').classed('highlighted', false);

    // Highlight the path
    let currentNode = treeRoot;
    path.forEach(nodeName => {
      const node = currentNode.children.find(child => child.data.name === nodeName);
      if (node) {
        d3.select(d3TreeContainer).select(`[data-name="${node.data.name}"]`).classed('highlighted', true);
        d3.select(d3TreeContainer).select(`[data-target="${node.data.name}"]`).classed('highlighted', true);
        currentNode = node;
      }
    });

    // Expand the path if necessary
    expandPath(treeRoot, path);
  }

  function expandPath(node, path) {
    if (path.length === 0) return;
    
    const childName = path[0];
    const child = node.children?.find(c => c.data.name === childName);
    
    if (child) {
      child.children = child._children;
      child._children = null;
      update(null, child);
      expandPath(child, path.slice(1));
    }
  }

  function handleScatterPointClick(event) {
    const { path } = event.detail;
    highlightPath(path);
  }

  function handleStepSelected(event) {
    const { step } = event.detail;
    currentStep = step;
  }

  function handleOperationsChanged(event) {
    const { step, operations } = event.detail;
    selectedOperations[step] = operations;
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
    console.log('Selected operations:', selectedOperations);
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

    fetch('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json')
      .then(response => response.json())
      .then(data => {
        console.log("Data loaded successfully:", data);
        treeData = data
        console.log("treeData after assignment:", treeData);
        if (d3TreeContainer) {
          renderD3Tree(d3TreeContainer, data);
        }
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

  .sidebar {
    width: 300px;
    padding: 20px;
    background-color: #f5f5f5;
    height: 100vh;
    overflow-y: auto;
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

  .stability-vis button:hover {
  background-color: #45a049;
}

.stability-vis button:disabled:hover {
  background-color: #919191;
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

  .view-toggle {
    display: flex;
    justify-content: center;
    margin-bottom: 1rem;
  }

  .view-toggle button {
    padding: 0.5rem 1rem;
    margin: 0 0.5rem;
    border: none;
    background-color: #f0f0f0;
    cursor: pointer;
    transition: background-color 0.3s;
  }

  .view-toggle button.active {
    background-color: #007bff;
    color: white;
  }

  .visualizations-section {
    margin-top: 2rem;
  }

  .visualization-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 1rem;
  }

  .visualization-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 1rem;
  }

  .card {
    border: 1px solid #ddd;
    border-radius: 4px;
    overflow: hidden;
    background-color: #fff;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }

  .card-header {
    background-color: #f5f5f5;
    padding: 1rem;
    border-bottom: 1px solid #ddd;
  }

  .card-header h3 {
    margin: 0;
    font-size: 1.2rem;
  }

  .card-content {
    padding: 1rem;
  }

  .interactive-placeholder {
    height: 300px;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: #f0f0f0;
    border-radius: 4px;
    font-style: italic;
    color: #666;
  }

  img {
    max-width: 100%;
    height: auto;
    margin-bottom: 1rem;
  }

  .card-container {
    display: flex;
    justify-content: space-between;
  }

  .card {
    flex: 1;
    border: 1px solid #ddd;
    border-radius: 4px;
    overflow: hidden;
    background-color: #fff;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin-right: 20px;
  }

  .card-header {
    background-color: #f5f5f5;
    padding: 1rem;
    border-bottom: 1px solid #ddd;
  }

  .card-header h3 {
    margin: 0;
    font-size: 1.2rem;
  }

  .card-content {
    padding: 1rem;
  }

  .floating-card {
    flex: 0.3;
    border: 1px solid #ddd;
    border-radius: 4px;
    background-color: #fff;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    overflow: hidden;
  }

  .floating-card .card-header {
    background-color: #f5f5f5;
    padding: 1rem;
    border-bottom: 1px solid #ddd;
  }

  .floating-card .card-content {
    padding: 1rem;
    max-height: calc(100% - 4rem);
    overflow-y: auto;
  }

  .floating-card ul {
    list-style-type: none;
    padding: 0;
    margin: 0;
  }

  .floating-card li {
    display: flex;
    justify-content: space-between;
    margin-bottom: 0.5rem;
  }

  .floating-card li button {
    background-color: #e74c3c;
    color: white;
    border: none;
    cursor: pointer;
    padding: 0.2rem 0.5rem;
  }

  .download-button {
    margin-top: 1rem;
    padding: 0.5rem 1rem;
    background-color: #007bff;
    color: white;
    border: none;
    cursor: pointer;
    transition: background-color 0.3s;
  }

  .download-button:hover {
    background-color: #0056b3;
  }

  .zoomed-image-container {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background-color: rgba(0, 0, 0, 0.8);
    display: flex;
    justify-content: center;
    align-items: center;
    z-index: 1000;
  }

  .zoomed-image-content {
    position: relative;
    max-width: 90%;
    max-height: 90%;
    background-color: white;
    padding: 20px;
    box-shadow: 0 0 20px rgba(0, 0, 0, 0.3);
  }

  .zoomed-image-content img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
  }

  .download-button {
    position: absolute;
    bottom: 1rem;
    right: 1rem;
    padding: 0.5rem 1rem;
    background-color: #007bff;
    color: white;
    border: none;
    cursor: pointer;
    transition: background-color 0.3s;
  }

  .download-button:hover {
    background-color: #0056b3;
  }

  .expand-collapse-button {
    position: sticky;
    bottom: 1rem;
    left: 50%;
    transform: translateX(-50%);
    z-index: 10;
    padding: 0.5rem 1rem;
    background-color: #007bff;
    color: white;
    border: none;
    cursor: pointer;
    transition: background-color 0.3s;
  }

  .expand-collapse-button:hover {
    background-color: #0056b3;
  }


  .tree-container-wrapper {
    display: flex;
    width: 100%;
    height: 200px;
    margin-bottom: 20px;
  }

  .d3-tree-container {
    flex: 0 0 80%; /* Take up 80% of the available space */
    height: 100%; /* Full height of the parent */
    border: 1px solid #ddd;
    background-color: #fff;
    overflow: auto; /* Allow scrolling if the tree content overflows */
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
    padding: 10px;
    border-radius: 8px;
    box-sizing: border-box;
  }

  .scatter-plot-container {
    flex: 0 0 20%; /* Take up 20% of the available space */
    height: 100%; /* Full height of the parent */
    background-color: #fff;
    border: 1px solid #ddd;
    border-radius: 4px;
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    box-sizing: border-box;
    border-radius: 8px;
    margin-left: 10px; /* Add some space between the containers */
  }

  .node.highlighted circle {
    fill: #ff0000;
  }

  .link.highlighted {
    stroke: #ff0000;
    stroke-width: 2px;
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
    <SidebarComponent 
      steps={steps} 
      currentStep={currentStep} 
      setCurrentStep={goToStep}
      on:stepSelected={handleStepSelected}
      on:operationsChanged={handleOperationsChanged}
    />
  </div>

  <!-- Main Content Area -->
  <div class="content">
    <div class="tree-container-wrapper">
      <!-- D3 Tree Container -->
      <div class="d3-tree-container" bind:this={d3TreeContainer}></div>
      <!-- Scatter Plot -->
      <div class="scatter-plot-container">
        <ScatterPlot on:pointClick={handleScatterPointClick} />
      </div> 
    </div>

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

    {#if currentStep === 'Raw data'}
      <div>
        <h1>Raw Data</h1>
        <!-- ASV File Upload UI -->
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
        
        {#if selectedOperations['Raw data']?.includes('Preview')}
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
        {/if}

        {#if selectedOperations['Raw data']?.includes('Set Random Seed')}
          <!-- Random Seed Input UI -->
          <div>
            <label for="randomSeed">Set Random Seed:</label>
            <input type="number" id="randomSeed" bind:value={randomSeed} min="1" />
          </div>
        {/if}

        {#if selectedOperations['Raw data']?.includes('Quick Explore')}
          <div class="tooltip">
            <button on:click={handleQuickExplore} disabled={asvFiles.length === 0 || !groupingsFile}>Quick Explore</button>
            <span class="tooltiptext">Upload files to continue</span>
          </div>
        {/if}

      </div>
    {:else if currentStep === 'Data Perturbation'}
      <div>
        <p>Dimensions: {filteredDimensions.rows} rows, {filteredDimensions.columns} columns</p>
        {#if selectedOperations['Data Perturbation']?.includes('Threshold')}
          <!-- Threshold Application UI -->
            <div class="filters">
              <label for="threshold">Threshold for Rare Genes:</label>
              <input type="range" id="threshold" bind:value={threshold} min="0" max="1" step="0.01" />
              <span>{threshold}</span>
              <button on:click={handleFilter}>Apply Threshold</button>
            </div>
        {/if}
        {#if selectedOperations['Data Perturbation']?.includes('Additional Option 1')}
          <!-- Additional Option 1 UI -->
          <div class="filters">
            <label for="option1">Additional Option 1:</label>
            <input type="text" id="option1" />
          </div>
        {/if}
        {#if selectedOperations['Data Perturbation']?.includes('Additional Option 2')}
          <!-- Additional Option 2 UI -->
          <div class="filters">
            <label for="option2">Additional Option 2:</label>
            <input type="text" id="option2" />
          </div>
        {/if}
      </div>
    {:else if currentStep === 'Model Perturbation'}
      <div>
        {#if selectedOperations['Model Perturbation']?.includes('Select Method')}
          <!-- Method Selection UI -->
          <div class="methods">
            <button on:click={() => handleMethodChange('deseq2')} class:selected={selectedMethod === 'deseq2'}>Method 1 (DESeq2)</button>
            <button on:click={() => handleMethodChange('aldex2')} class:selected={selectedMethod === 'aldex2'}>Method 2 (ALDEx2)</button>
            <button on:click={() => handleMethodChange('edger')} class:selected={selectedMethod === 'edger'}>Method 3 (edgeR)</button>
            <button on:click={() => handleMethodChange('maaslin2')} class:selected={selectedMethod === 'maaslin2'}>Method 4 (Maaslin2)</button>
            <button on:click={() => handleMethodChange('method5')} class:selected={selectedMethod === 'method5'}>Method 5</button>
          </div>
        {/if}
      </div>
    {:else if currentStep === 'Stability Metric'}
      <div>
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
  
        <div class="combined-results-actions">
          <button on:click={generateCombinedResults} disabled={isCombiningResults || !Object.values(methodFileStatus).every(status => status[0])}>
            {isCombiningResults ? 'Generating...' : 'Generate Combined Results'}
          </button>
          <button on:click={downloadCombinedResults} disabled={!combinedResultsReady}>
            Download Combined Results
          </button>
        </div>

        {#if selectedOperations['Stability Metric']?.includes('View Stability Plot')}
          <!-- Stability Plot Viewing UI -->
          {#if stabilityPlot}
            <div class="stability-plot">
              <h3>Stability Plot</h3>
              <img src={stabilityPlot} alt="Stability Plot" style="width: 100%; max-width: 800px; height: auto;" />
            </div>
          {/if}
        {/if}
        {#if selectedOperations['Stability Metric']?.includes('Run Shuffled Analysis')}
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

        {#if selectedOperations['Stability Metric']?.includes('ASV Selector')}
          <button on:click={() => showASVSelector = !showASVSelector} disabled={!combinedResultsReady}>
            {showASVSelector ? 'Hide' : 'Show'} ASV Selector
          </button>
          <div id="asv-selector-main">
            {#if showASVSelector}
              <ASVSelector />
            {/if}
          </div>
        {/if}
        
      </div>
    {/if}

    <div class="stability-vis" hidden={currentStep !== 'Stability Metric'}>

      <div class="shuffled-analysis">
        
      </div>
    </div>    

    <div class="visualizations-section" hidden={!showAllPlots && !isSubmitted}>
      <div class="visualization-header">
        <h2>Visualizations</h2>
        <div class="view-toggle">
          <button class:active={isStatic} on:click={() => toggleView('static')}>Static</button>
          <button class:active={!isStatic} on:click={() => toggleView('interactive')}>Interactive</button>
        </div>
      </div>
    
      {#if isCalculating}
        <div class="loader">
          <p>Loading...</p>
        </div>
      {:else}
        {#if showAllPlots}
          <div class="card-container">
            <div class="card">
              <div class="card-header">
                <h3>Overlap Visualizations</h3>
              </div>
              <div class="card-content">
                {#if isStatic}
                  <img src={visualizations.overlap_volcano} alt="Overlap Volcano Plot" on:click={() => zoomImage('overlap_volcano')} />
                  <img src={visualizations.overlap_pvalue_distribution} alt="Overlap P-value Distribution" on:click={() => zoomImage('overlap_pvalue_distribution')} />
                {:else}
                  <InteractiveValcano />
                {/if}
              </div>
            </div>

            <!-- floating card -->
            {#if selectedPointsList.length > 0}
              <div class="floating-card">
                <div class="card-header">
                  <h3>Selected Points</h3>
                </div>
                <div class="card-content">
                  <ul>
                    {#each selectedPointsList as point (point.name)}
                      <li>
                        <span>{point.name}</span>
                        <button on:click={() => removePoint(point)}>Remove</button>
                      </li>
                    {/each}
                  </ul>
                  <button class="download-button" on:click={downloadSelectedPoints}>Download Points</button>
                </div>
              </div>
            {/if}
    
            {#if showDetailedPlots}
              {#each ['deseq2', 'aldex2', 'edger', 'maaslin2'] as method}
                <div class="card">
                  <div class="card-header">
                    <h3>{method.toUpperCase()} Plots</h3>
                  </div>
                  <div class="card-content">
                    {#if isStatic}
                      {#each [1, 2, 3] as plotNumber}
                        <img 
                          src={visualizations[`${method}_plot${plotNumber}`]} 
                          alt={`${method.toUpperCase()} Plot ${plotNumber}`} 
                          on:click={() => zoomImage(`${method}_plot${plotNumber}`)}
                        />
                      {/each}
                    {:else}
                      <div class="interactive-placeholder">
                        Interactive {method.toUpperCase()} Plots coming soon...
                      </div>
                    {/if}
                  </div>
                </div>
              {/each}
            {/if}
          </div>
    
          <button class="expand-collapse-button" on:click={() => showDetailedPlots = !showDetailedPlots}>
            {showDetailedPlots ? 'Collapse Details' : 'Show More Details'}
          </button>
        {:else if selectedMethod}
          <div class="card">
            <div class="card-header">
              <h3>{selectedMethod.toUpperCase()} Plots</h3>
            </div>
            <div class="card-content">
              {#if isStatic}
                {#each [1, 2, 3] as plotNumber}
                  <img 
                    src={visualizations[`${selectedMethod}_plot${plotNumber}`]} 
                    alt={`${selectedMethod.toUpperCase()} Plot ${plotNumber}`} 
                    on:click={() => zoomImage(`${selectedMethod}_plot${plotNumber}`)}
                  />
                {/each}
              {:else}
                <div class="interactive-placeholder">
                  Interactive {selectedMethod.toUpperCase()} Plots coming soon...
                </div>
              {/if}
            </div>
          </div>
        {/if}
      {/if}
    </div>

    
    {#if zoomedImage}
      <div class="zoomed-image-container" on:click={() => zoomImage(null)} transition:fade>
        <div class="zoomed-image-content" on:click|stopPropagation transition:scale>
          <img src={visualizations[zoomedImage]} alt="Zoomed Image" />
          <button class="download-button" on:click|stopPropagation={() => downloadImage(zoomedImage, zoomedImage.split('_')[0])}>
            Download
          </button>
        </div>
      </div>
    {/if}
    
    
    {#if currentStep === 'Model Perturbation' && asvFiles.length > 0 && groupingsFile && selectedMethod}
      <button on:click={handleSubmit}>Submit</button>
      <button on:click={handleDownload} disabled={!selectedMethod || !isSubmitted}>Download</button>
    {/if}

    <!-- <div class="navigation">
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
    </div> -->

  </div>
</div>
