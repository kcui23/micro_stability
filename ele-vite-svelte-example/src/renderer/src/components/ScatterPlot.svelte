<script>
  import { onMount, createEventDispatcher } from 'svelte';
  import * as d3 from 'd3';
  import { singleSelectOperations, 
    selectedColorStep, 
    scatterPlotColors,
    colorStatus } from '../store.js';
  import { fade } from 'svelte/transition';

  const dispatch = createEventDispatcher();
  let svg;
  let data;
  let leafIdDataPoints;
  let selectedPointId = null;
  let gradientId;
  let overlayVisible = true;

  export let data_points_updated_counter;
  export let highlight_point_path;

  const pointRadius = 3;
  const selectedPointRadius = 7;

  $: if (data_points_updated_counter) {
    main_function();
  }

  $: if (highlight_point_path) {
    highlightPoint(highlight_point_path);
  }

  $: if ($selectedColorStep || $colorStatus) {
    updateScatterPlot();
  }

  async function fetchData(url) {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error('Failed to fetch data from ' + url);
    }
    return await response.json();
  }

  onMount(() => {
    main_function();
  });

  function updatePointsOrder() {
    d3.select(svg).selectAll('circle')
      .sort((a, b) => {
        if (a.id === selectedPointId) return 1;
        if (b.id === selectedPointId) return -1;
        return 0;
      });
  }

  async function main_function() {
    try {
      const response = await fetch('http://localhost:8000/get_plot_data');
      if (!response.ok) {
        throw new Error('Failed to fetch plot data');
      }
      const plotData = await response.json();
      data = plotData.tree_data;
      leafIdDataPoints = plotData.leaf_data;

      initializeSVG();
      updateScatterPlot();
    } catch (error) {
      console.error('Error fetching or processing data:', error);
    }
  }

  function initializeSVG() {
    const width = 300;
    const height = 300;
    const margin = { top: 20, right: 20, bottom: 40, left: 40 };

    d3.select(svg).selectAll('*').remove();

    const svgElement = d3.select(svg)
      .attr('width', '100%')
      .attr('height', '100%')
      .attr('viewBox', `0 0 ${width} ${height}`);

    // Initialize gradient
    const defs = svgElement.append('defs');
    gradientId = `gradient-${Math.random().toString(36).substr(2, 9)}`;
    const gradient = defs.append('linearGradient')
      .attr('id', gradientId)
      .attr('x1', '0%')
      .attr('y1', '0%')
      .attr('x2', '100%')
      .attr('y2', '100%');
    
    gradient.append('stop')
      .attr('offset', '0%')
      .attr('stop-color', 'CornflowerBlue');

    gradient.append('stop')
      .attr('offset', '65%')
      .attr('stop-color', 'DarkSalmon');
    
    gradient.append('stop')
      .attr('offset', '100%')
      .attr('stop-color', 'MistyRose');

    // Initialize axes
    svgElement.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${height - margin.bottom})`);

    svgElement.append('g')
      .attr('class', 'y-axis')
      .attr('transform', `translate(${margin.left},0)`);

    // Initialize axis labels
    svgElement.append('text')
      .attr('class', 'x-axis-label')
      .attr('x', width / 2 + 10)
      .attr('y', height - margin.bottom + 25)
      .attr('fill', 'black')
      .attr('text-anchor', 'middle')
      .attr('font-size', '12px')
      .text('UMAP Dimension 1');

    svgElement.append('text')
      .attr('class', 'y-axis-label')
      .attr('transform', `rotate(-90)`)
      .attr('x', -height / 2)
      .attr('y', margin.left - 20)
      .attr('fill', 'black')
      .attr('text-anchor', 'middle')
      .attr('font-size', '12px')
      .text('UMAP Dimension 2');
  }

  function updateScatterPlot() {
    if (!data) return;

    const width = 300;
    const height = 300;
    const margin = { top: 20, right: 20, bottom: 40, left: 40 };

    const x = d3.scaleLinear().range([margin.left, width - margin.right]);
    const y = d3.scaleLinear().range([height - margin.bottom, margin.top]);

    const svgElement = d3.select(svg);

    const points = extractDataPoints(data).filter(d => isPointVisible(d));

    // Remove existing "No data available" text
    svgElement.select('.no-data-text').remove();

    const currentStepStatus = $colorStatus[$selectedColorStep];
    if (points.length > 0 && currentStepStatus && currentStepStatus.length > 0) {
      x.domain(d3.extent(points, d => d.x));
      y.domain(d3.extent(points, d => d.y));

      svgElement.selectAll('circle')
        .data(points, d => d.id)
        .join(
          enter => enter.append('circle')
            .attr('cx', d => x(d.x))
            .attr('cy', d => y(d.y))
            .attr('r', pointRadius)
            .attr('fill', d => getColorForPoint(d))
            .on('click', (event, d) => handlePointClick(d)),
          update => update
            .attr('cx', d => x(d.x))
            .attr('cy', d => y(d.y))
            .attr('fill', d => d.id === selectedPointId ? `url(#${gradientId})` : getColorForPoint(d)),
          exit => exit.remove()
        )
        .attr('r', d => d.id === selectedPointId ? selectedPointRadius : pointRadius)
        .attr('stroke', d => d.id === selectedPointId ? 'black' : 'none')
        .attr('stroke-width', d => d.id === selectedPointId ? 2 : 0);

      updatePointsOrder();

      // Update axes
      svgElement.select('.x-axis').call(d3.axisBottom(x));
      svgElement.select('.y-axis').call(d3.axisLeft(y));
    } else {
      svgElement.selectAll('circle').remove();
      
      // Add the "No data available" text if there are no points to display
      if (currentStepStatus && currentStepStatus.length === 0) {
        svgElement.append('text')
          .attr('class', 'no-data-text')
          .attr('x', width / 2)
          .attr('y', height / 2)
          .attr('text-anchor', 'middle')
          .text('No color tag selected');
      } else {
        svgElement.append('text')
          .attr('class', 'no-data-text')
          .attr('x', width / 2)
          .attr('y', height / 2)
          .attr('text-anchor', 'middle')
          .text('No data points available');
      }
    }
  }

  function getColorForPoint(d) {
    const stepColors = $scatterPlotColors[$selectedColorStep];
    if (!stepColors) return 'pink';

    const legal_ops = $singleSelectOperations[$selectedColorStep];
    const stepName = d.path.find(op => legal_ops.includes(op));
    const index = legal_ops.indexOf(stepName);
    return stepColors[index] || 'black';
  }

  function isPointVisible(d) {
    const currentStepStatus = $colorStatus[$selectedColorStep];
    if (!currentStepStatus || currentStepStatus.length === 0) return false;
    return d.path.some(step => currentStepStatus.includes(step));
  }

  function handlePointClick(d) {
    overlayVisible = false;
    selectedPointId = d.id;
    dispatch('pointClick', { path: d.path });
    highlightPoint(d.path);
    updatePointsOrder();
  }

  function highlightPoint(path) {
    if (!path) return;
    
    const points = extractDataPoints(data);
    const selectedPoint = points.find(p => p.path.toString() === path.toString());
    
    if (selectedPoint) {
      selectedPointId = selectedPoint.id;
      d3.select(svg).selectAll('circle')
        .transition()
        .duration(100)
        .attr('stroke', d => d.id === selectedPointId ? 'black' : 'none')
        .attr('stroke-width', d => d.id === selectedPointId ? 2 : 0)
        .attr('fill', d => d.id === selectedPointId ? `url(#${gradientId})` : getColorForPoint(d))
        .attr('r', d => d.id === selectedPointId ? selectedPointRadius : pointRadius);
      updatePointsOrder();
    }
  }

  function extractDataPoints(node, path = []) {
    let points = [];
    if (node && node.id && leafIdDataPoints[node.id]) {
      const dataPoint = leafIdDataPoints[node.id].data_point;
      if (dataPoint) {
        const flatPath = [...path, node.name].flat();
        points.push({ 
          id: node.id, 
          x: dataPoint[0], 
          y: dataPoint[1], 
          path: flatPath 
        });
      }
    }
    if (node && node.children) {
      for (const child of node.children) {
        points = points.concat(extractDataPoints(child, [...path, node.name].flat()));
      }
    }
    return points;
  }

  function handleOverlayClick() {
    overlayVisible = false;
  }
</script>

<div class="scatter-plot-container">
  <svg bind:this={svg}></svg>
  {#if overlayVisible}
    <div 
      class="overlay"
      transition:fade={{ duration: 300 }}
      role="button" 
      tabindex="0"
      on:click={handleOverlayClick}
      on:keydown={e => e.key === 'Enter' && handleOverlayClick()}>
      <span>Click any point to start</span>
    </div>
  {/if}
</div>

<style>
  .scatter-plot-container {
    position: relative;
    flex: 0 0 80%;
    height: 100%;
    background-color: #fff;
    overflow: hidden;
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
    border-radius: 8px;
    box-sizing: border-box;
    z-index: 998;
  }

  .overlay {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background-color: rgba(200, 200, 200, 0.7);
    display: flex;
    justify-content: center;
    align-items: center;
    cursor: default;
    backdrop-filter: blur(2px);
    z-index: 999;
  }

  .overlay span {
    color: #333;
    font-size: 1.2em;
    font-weight: bold;
    text-shadow: 0 0 10px white;
  }
</style>