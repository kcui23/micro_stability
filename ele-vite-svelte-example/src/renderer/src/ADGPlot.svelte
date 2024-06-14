<script>
    import { onMount } from 'svelte';
    import { select } from 'd3-selection';
    import { graphlib, layout } from 'dagre';
  
    export let steps;
    export let currentStep;
    export let setCurrentStep;
    export let edgeThicknesses = [];
    const basicThickness = 2;
  
    const drawADG = () => {
        const g = new graphlib.Graph().setGraph({}).setDefaultEdgeLabel(() => ({}));
  
        // Nodes
        steps.forEach((step, index) => {
            g.setNode(index, { label: step, x: 0, y: index * 50 }); // Initial guess for positions
        });
  
        // Edges
        for (let i = 0; i < steps.length - 1; i++) {
            g.setEdge(i, i + 1, { thickness: edgeThicknesses[i] * basicThickness || 2 });
        }
  
        layout(g);
  
        const svg = select('#adg-container');
        svg.selectAll('*').remove(); // Clear any previous content
        const inner = svg.append('g');
  
        g.nodes().forEach((v) => {
            const node = g.node(v);
            if (node && !isNaN(node.x) && !isNaN(node.y)) {
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
                    .attr('x', node.x + 15) // Position text to the right of the node
                    .attr('y', node.y)
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'start') // Anchor text to the start (left side)
                    .text(node.label)
                    .style('cursor', 'pointer')
                    .on('click', () => {
                        setCurrentStep(steps[v]);
                    });
            } else {
                console.error('Invalid node coordinates:', node);
            }
        });
  
        g.edges().forEach((e) => {
            const edge = g.edge(e);
            const points = edge.points;
            if (points && points.length > 1) {
                points.forEach(point => {
                    if (isNaN(point.y)) {
                        point.y = g.node(e.w).y;
                    }
                });
  
                inner.append('path')
                    .attr('d', `M${points[0].x},${points[0].y-40} ${points.slice(1).map(p => `L${p.x},${p.y-10}`).join(' ')}`)
                    .attr('stroke-width', edge.thickness)
                    .attr('stroke', 'darkgray')
            } else {
                console.error('Invalid edge points:', points);
            }
        });
    };
  
    onMount(() => {
        drawADG();
    });
  
    $: if (steps.length > 0) {
        drawADG();
    }
  
    $: if (currentStep) {
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
  