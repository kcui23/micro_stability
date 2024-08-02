<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let plotDiv;

  onMount(async () => {
    try {
      const response = await fetch('/Users/kai/Downloads/overlap_combined_results.tsv');
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
              case 'ALDEx2': return 'diamond';
              case 'edgeR': return 'square';
              case 'MaAsLin2': return 'triangle-up';
              default: return 'circle';
            }
          }),
          opacity: 0.5
        }
      };

      const layout = {
        title: 'Volcano plot overlap',
        xaxis: { title: 'log2FoldChange' },
        yaxis: { title: '-log10(pvalue)' },
        showlegend: false,
        hovermode: 'closest'
      };

      Plotly.newPlot(plotDiv, [trace], layout, { responsive: true });
    } catch (error) {
      console.error('Error fetching or plotting data:', error);
      // You might want to display an error message to the user here
    }
  });
</script>

<div bind:this={plotDiv} style="width: 100%; height: 600px;"></div>

<style>
  div {
    margin: 0 auto;
  }
</style>