<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { selectedPoints } from '../store.js';

  export let specific_interact
  export let selectedMethod
  let plotDiv;

  onMount(async () => {
    try {
      const response = await fetch(`http://localhost:8000/get_overlap_combined_results?specific_interact=${specific_interact}&selectedMethod=${selectedMethod}`);
      if (!response.ok) {
        throw new Error('Failed to fetch data');
      }
      const tsvData = await response.text();

      // Parse TSV data
      const rows = tsvData.trim().split('\n').map(row => row.split('\t'));
      const headers = rows[0];
      const values = rows.slice(1);

      const asvNames = values.map(row => row[headers.indexOf('asv_name')]);
      const log2FoldChanges = values.map(row => parseFloat(row[headers.indexOf('log2FoldChange')]));
      const negLog10Pvalues = values.map(row => -Math.log10(parseFloat(row[headers.indexOf('pvalue')])));
      const significantValues = values.map(row => parseFloat(row[headers.indexOf('pvalue')]) < 0.05 ? 'Significant' : 'Not Significant');
      const methods = values.map(row => row[headers.indexOf('method')]);

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
        title: 'Volcano plot overlap',
        xaxis: { title: 'log2FoldChange' },
        yaxis: { title: '-log10(pvalue)' },
        showlegend: false,
        hovermode: 'closest'
      };

      // draw the plot
      Plotly.newPlot(plotDiv, [trace], layout, { responsive: true });

      // add click event listener
      plotDiv.on('plotly_click', function(data) {
        // get the point that was clicked
        const point = data.points[0];
        console.log('Clicked point:', point);
        console.log('x:', point.x);
        console.log('y:', point.y);
        console.log('ASV Name:', point.text);
      });
    } catch (error) {
      console.error('Error fetching or plotting data:', error);
    }
    
    plotDiv.on('plotly_click', function(data) {
      const point = data.points[0];
      selectedPoints.update(points => {
        const newPoints = [...points, { name: point.text }];
        console.log("Updated selectedPoints:", newPoints)
        return newPoints;
      });
    });
  });
</script>

<div bind:this={plotDiv} style="width: 100%; height: 600px;"></div>

<style>
  div {
    margin: 0 auto;
  }
</style>