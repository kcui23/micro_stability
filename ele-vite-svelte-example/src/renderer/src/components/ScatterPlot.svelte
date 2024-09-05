<script>
  import { onMount, createEventDispatcher } from 'svelte';
  import * as d3 from 'd3';

  const dispatch = createEventDispatcher();
  let svg;
  let data;
  let leafIdDataPoints;
  let selectedPointId = null;
  
  export let data_points_updated_counter;
  export let highlight_point_path;

  $: if (data_points_updated_counter) {
    main_function();
  }

  $: if (highlight_point_path) {
    highlightPoint(highlight_point_path);
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

  async function main_function() {
    try {
      const [dataResponse, leafIdDataPointsResponse] = await Promise.all([
        fetchData('https://raw.githubusercontent.com/kcui23/micro_stability/main/ele-vite-svelte-example/src/renderer/src/public/data.json'),
        fetchData('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/leaf_id_data_points.json')
      ]);
      data = dataResponse;
      leafIdDataPoints = leafIdDataPointsResponse;

      const width = 300;
      const height = 300;
      const margin = { top: 20, right: 20, bottom: 30, left: 40 };

      const x = d3.scaleLinear().range([margin.left, width - margin.right]);
      const y = d3.scaleLinear().range([height - margin.bottom, margin.top]);

      d3.select(svg).selectAll('*').remove();

      const svgElement = d3.select(svg)
        .attr('width', '100%')
        .attr('height', '100%')
        .attr('viewBox', `0 0 ${width} ${height}`);

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
              .attr('r', 5)
              .attr('fill', 'steelblue')
              .on('click', (event, d) => handlePointClick(d)),
            update => update
              .attr('cx', d => x(d.x))
              .attr('cy', d => y(d.y))
          )
          .attr('r', d => d.id === selectedPointId ? 10 : 5)
          .attr('fill', d => d.id === selectedPointId ? 'orange' : 'steelblue');

        svgElement.append('g')
          .attr('transform', `translate(0,${height - margin.bottom})`)
          .call(d3.axisBottom(x));

        svgElement.append('g')
          .attr('transform', `translate(${margin.left},0)`)
          .call(d3.axisLeft(y));
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

  function handlePointClick(d) {
    selectedPointId = d.id;
    dispatch('pointClick', { path: d.path });
    highlightPoint(d.path);
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
        .attr('fill', d => d.id === selectedPointId ? 'orange' : 'steelblue')
        .attr('r', d => d.id === selectedPointId ? 10 : 5);
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