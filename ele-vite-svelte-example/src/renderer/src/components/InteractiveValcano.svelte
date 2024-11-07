<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { selectedPoints, currentHighlightedPath } from '../store.js';

  export let specific_interact;
  export let selectedMethod;
  let plotDiv;

  onMount(async () => {
    try {
      const response = await fetch(`http://localhost:8000/get_overlap_combined_results?specific_interact=${specific_interact}&selectedMethod=${selectedMethod}`,
        {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({ path: $currentHighlightedPath })
        }
      );
      if (!response.ok) {
        throw new Error('Failed to fetch data');
      }
      const data = await response.json();
      
      if (data.error) {
        throw new Error(data.error);
      }

      // Extract arrays from the JSON data
      const asvNames = data.asv_name;
      const log2FoldChanges = data.log2FoldChange.map(val => parseFloat(val));
      const pvalues = data.pvalue.map(val => parseFloat(val));
      const negLog10Pvalues = pvalues.map(p => -Math.log10(p));
      const significantValues = pvalues.map(p => p < 0.05 ? 'Significant' : 'Not Significant');
      const methods = data.method || Array(asvNames.length).fill(selectedMethod);

      const trace = {
        x: log2FoldChanges,
        y: negLog10Pvalues,
        text: asvNames,
        mode: 'markers',
        type: 'scatter',
        marker: {
          size: 8,
          color: significantValues.map(val => val === 'Significant' ? '#0000FF' : '#808080'),
          symbol: methods.map(method => {
            switch (method) {
              case 'DESeq2': return 'circle';
              case 'ALDEx2': return 'triangle-up';
              case 'edgeR': return 'square';
              case 'MaAsLin2': return 'triangle-down';
              case 'metagenomeSeq': return 'diamond';
              default: return 'circle';
            }
          }),
          opacity: 0.5
        },
        hoverinfo: 'text',
        hovertext: asvNames.map((name, i) => 
          `Coordinates: (${log2FoldChanges[i].toFixed(2)}, ${negLog10Pvalues[i].toFixed(2)})<br>` +
          `Method: ${methods[i]}<br>` +
          `ASV: ${name}`
        )
      };

      const layout = {
        title: specific_interact ? 'Click to select ASV' : 'Volcano Plot Overlap',
        xaxis: { title: 'log2FoldChange' },
        yaxis: { title: '-log10(pvalue)' },
        showlegend: false,
        hovermode: 'closest'
      };

      // Draw the plot
      Plotly.newPlot(plotDiv, [trace], layout, { responsive: true });

      // Add click event listener
      plotDiv.on('plotly_click', function(data) {
        const point = data.points[0];
        selectedPoints.update(points => {
          // Check if point already exists
          const existingPointIndex = points.findIndex(p => p.name === point.text);
          
          if (existingPointIndex !== -1) {
            // Remove point if already selected
            return points.filter((_, index) => index !== existingPointIndex);
          } else {
            // Add new point
            return [...points, { 
              name: point.text,
              x: point.x,
              y: point.y,
              method: methods[asvNames.indexOf(point.text)]
            }];
          }
        });
      });

    } catch (error) {
      console.error('Error fetching or plotting data:', error);
    }
  });
</script>

<div bind:this={plotDiv} style="width: 100%; height: 600px;"></div>

<style>
  div {
    margin: 0 auto;
  }
</style>