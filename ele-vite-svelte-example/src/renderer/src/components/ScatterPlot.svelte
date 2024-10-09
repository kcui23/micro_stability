<script>
  import { onMount, createEventDispatcher } from 'svelte';
  import * as d3 from 'd3';
  import { singleSelectOperations, selectedColorStep, scatterPlotColors } from '../store.js'; // Assume scatterPlotColors is an object with color schemes

  const dispatch = createEventDispatcher();
  let svg;
  let data;
  let leafIdDataPoints;
  let selectedPointId = null;
  
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

  // Watch for selectedColorStep updates and refresh the scatter plot colors
  $: if (selectedColorStep) {
    updatePointColors();
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
      const [dataResponse, leafIdDataPointsResponse] = await Promise.all([
        fetchData('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json'),
        fetchData('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json')
      ]);
      data = dataResponse;
      leafIdDataPoints = leafIdDataPointsResponse;

      const width = 300;
      const height = 300;
      const margin = { top: 20, right: 20, bottom: 40, left: 40 };

      const x = d3.scaleLinear().range([margin.left, width - margin.right]);
      const y = d3.scaleLinear().range([height - margin.bottom, margin.top]);

      d3.select(svg).selectAll('*').remove();

      const svgElement = d3.select(svg)
        .attr('width', '100%')
        .attr('height', '100%')
        .attr('viewBox', `0 0 ${width} ${height}`);

      if (svgElement.select('.x-axis').empty()) {
        svgElement.append('g')
          .attr('class', 'x-axis')
          .attr('transform', `translate(0,${height - margin.bottom})`);
      }
      if (svgElement.select('.y-axis').empty()) {
        svgElement.append('g')
          .attr('class', 'y-axis')
          .attr('transform', `translate(${margin.left},0)`);
      }

      const points = extractDataPoints(data);

      if (points.length > 0) {
        x.domain(d3.extent(points, d => d.x));
        y.domain(d3.extent(points, d => d.y));

        svgElement.selectAll('circle')
          .data(points, d => d.id)
          .join(
            enter => enter.append('circle')
              .attr('cx', d => x(d.x))
              .attr('cy', d => y(d.y))
              .attr('r', pointRadius)
              .attr('fill', d => getColorForPoint(d)) // Update to use the color function
              .on('click', (event, d) => handlePointClick(d)),
            update => update
              .attr('cx', d => x(d.x))
              .attr('cy', d => y(d.y))
          )
          .attr('r', d => d.id === selectedPointId ? selectedPointRadius : pointRadius)
          .attr('fill', d => d.id === selectedPointId ? 'orange' : getColorForPoint(d)) // Use color based on step
          .attr('stroke', d => d.id === selectedPointId ? 'black' : 'none')
          .attr('stroke-width', d => d.id === selectedPointId ? 2 : 0);

        updatePointsOrder();

        // Update x-axis
        svgElement.select('.x-axis')
          .call(d3.axisBottom(x));

        // Update y-axis
        svgElement.select('.y-axis')
          .call(d3.axisLeft(y));

        // Add or update x-axis label
        let xLabel = svgElement.select('.x-axis-label');
        if (xLabel.empty()) {
          xLabel = svgElement.append('text')
            .attr('class', 'x-axis-label')
            .attr('x', width / 2 + 10)
            .attr('y', height - margin.bottom + 25)
            .attr('fill', 'black')
            .attr('text-anchor', 'middle')
            .attr('font-size', '12px')
            .text('Stability Metric 1');
        }

        // Add or update y-axis label
        let yLabel = svgElement.select('.y-axis-label');
        if (yLabel.empty()) {
          yLabel = svgElement.append('text')
            .attr('class', 'y-axis-label')
            .attr('transform', `rotate(-90)`)
            .attr('x', -height / 2)
            .attr('y', margin.left - 20)
            .attr('fill', 'black')
            .attr('text-anchor', 'middle')
            .attr('font-size', '12px')
            .text('Stability Metric 2');
        }
      } else {
        svgElement.append('text')
          .attr('x', width / 2)
          .attr('y', height / 2)
          .attr('text-anchor', 'middle')
          .text('No data points available');
      }
    } catch (error) {
      console.error('Error fetching or processing data:', error);
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

  function updatePointColors() {
    d3.select(svg).selectAll('circle')
      .attr('fill', d => d.id === selectedPointId ? 'orange' : getColorForPoint(d));
  }

  function handlePointClick(d) {
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
        .attr('fill', d => d.id === selectedPointId ? 'orange' : getColorForPoint(d))
        .attr('r', d => d.id === selectedPointId ? selectedPointRadius : pointRadius);
      updatePointsOrder();
    }
  }

  function extractDataPoints(node, path = []) {
    let points = [];
    if (node && node.id && leafIdDataPoints[node.id]) {
      const dataPoint = leafIdDataPoints[node.id].data_point;
      if (dataPoint) {
        points.push({ 
          id: node.id, 
          x: dataPoint[0], 
          y: dataPoint[1], 
          path: [...path, node.name] 
        });
      }
    }
    if (node && node.children) {
      for (const child of node.children) {
        points = points.concat(extractDataPoints(child, [...path, node.name]));
      }
    }
    return points;
  }
</script>

<svg bind:this={svg}></svg>
