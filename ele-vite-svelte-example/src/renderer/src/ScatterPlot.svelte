<script>
    import { onMount, createEventDispatcher } from 'svelte';
    import * as d3 from 'd3';
  
    const dispatch = createEventDispatcher();
    let svg;
    let data;
  
    onMount(async () => {
      try {
        const response = await fetch('/Users/kai/Desktop/MSDS/micro_stability/ele-vite-svelte-example/src/renderer/src/public/data.json');
        if (!response.ok) {
          throw new Error('Failed to fetch data');
        }
        data = await response.json();
  
        const width = 300;
        const height = 300;
        const margin = { top: 20, right: 20, bottom: 30, left: 40 };
  
        const x = d3.scaleLinear().range([margin.left, width - margin.right]);
        const y = d3.scaleLinear().range([height - margin.bottom, margin.top]);
  
        const svgElement = d3.select(svg)
          .attr('width', '100%')
          .attr('height', '100%')
          .attr('viewBox', `0 0 ${width} ${height}`);
  
        const points = extractDataPoints(data);
  
        if (points.length > 0) {
          x.domain(d3.extent(points, d => d.x));
          y.domain(d3.extent(points, d => d.y));
  
          svgElement.selectAll('circle')
            .data(points)
            .enter()
            .append('circle')
            .attr('cx', d => x(d.x))
            .attr('cy', d => y(d.y))
            .attr('r', 5)
            .attr('fill', 'steelblue')
            .on('click', (event, d) => {
              dispatch('pointClick', { path: d.path });
            });
  
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
    });
  
    function extractDataPoints(node, path = []) {
      let points = [];
      if (node && node.data_point) {
        points.push({ ...node.data_point, path: [...path, node.name] });
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
