<script>
    import { onMount } from 'svelte';
    import { select } from 'd3-selection';
    import { graphlib, layout } from 'dagre';
    export let steps;
    export let currentStep;
    export let setCurrentStep;
  
    const drawADG = () => {
      const g = new graphlib.Graph().setGraph({}).setDefaultEdgeLabel(() => ({}));
  
      // Nodes
      steps.forEach((step, index) => {
        g.setNode(index, { label: step });
      });
  
      // Edges
      for (let i = 0; i < steps.length - 1; i++) {
        g.setEdge(i, i + 1);
      }
  
      layout(g);
  
      const svg = select('#adg-container');
      svg.selectAll('*').remove(); // Clear any previous content
      const inner = svg.append('g');
  
      g.nodes().forEach((v) => {
        const node = g.node(v);
        inner.append('circle')
          .attr('r', 10)
          .attr('cx', node.x)
          .attr('cy', node.y)
          .attr('class', currentStep === steps[v] ? 'active' : '')
          .style('fill', currentStep === steps[v] ? '#007bff' : '#ccc')
          .style('cursor', 'pointer')
          .on('click', () => {
            setCurrentStep(steps[v]);
          });
  
        inner.append('text')
          .attr('x', node.x)
          .attr('y', node.y - 15)
          .attr('dy', '0.35em')
          .attr('text-anchor', 'middle')
          .text(node.label)
          .style('cursor', 'pointer')
          .on('click', () => {
            setCurrentStep(steps[v]);
          });
      });
  
      g.edges().forEach((e) => {
        const edge = g.edge(e);
        inner.append('path')
          .attr('d', `M${edge.points[0].x},${edge.points[0].y}L${edge.points[1].x},${edge.points[1].y}`)
          .attr('marker-end', 'url(#arrowhead)');
      });
  
      // Define arrowhead marker
      svg.append('defs').append('marker')
        .attr('id', 'arrowhead')
        .attr('viewBox', '-0 -5 10 10')
        .attr('refX', 13)
        .attr('refY', 0)
        .attr('orient', 'auto')
        .attr('markerWidth', 6)
        .attr('markerHeight', 6)
        .append('path')
        .attr('d', 'M0,-5L10,0L0,5')
        .attr('stroke', '#000')
        .attr('fill', '#000');
    };
  
    onMount(() => {
      drawADG();
    });
  
    $: if (steps.length > 0) {
      drawADG();
    }
  </script>
  
  <svg id="adg-container" width="200" height="600"></svg>
  
  <style>
    svg {
    overflow: visible;
  }
    circle {
      cursor: pointer;
    }
    circle.active {
      fill: #007bff;
    }
    text {
      font-size: 12px;
      fill: #333;
    }
  </style>
  