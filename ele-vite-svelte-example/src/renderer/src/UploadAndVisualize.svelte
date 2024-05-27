<script>
  let files = [];
  let visualizations = { plot1: '', plot2: '', plot3: '' };
  let selectedMethod = '';

  const handleFileChange = async (event) => {
    files = event.target.files;
  };

  const handleMethodChange = (event) => {
    selectedMethod = event.target.value;
  };

  const handleSubmit = async () => {
    const file = files[0];
    const reader = new FileReader();
    reader.onload = async () => {
      const fileContent = reader.result;

      const response = await fetch(`http://localhost:8000/process?method=${selectedMethod}`, {
        method: 'POST',
        headers: {
          'Content-Type': 'text/plain'
        },
        body: fileContent
      });

      const result = await response.json();
      visualizations = {
        plot1: `data:image/png;base64,${result.plot1}`,
        plot2: `data:image/png;base64,${result.plot2}`,
        plot3: `data:image/png;base64,${result.plot3}`
      };
    };
    reader.readAsText(file);
  };
</script>

<div>
  <h1>Upload ASV Dataset</h1>
  <input type="file" accept=".tsv" on:change={handleFileChange} />

  {#if files.length > 0}
    <p>File: {files[0].name}</p>
    
    <label>
      Select Method:
      <select on:change={handleMethodChange}>
        <option value="" disabled selected>Select your option</option>
        <option value="deseq2">DESeq2</option>
        <option value="aldex2">Aldex2</option>
      </select>
    </label>
    
    <button on:click={handleSubmit}>Submit</button>
  {/if}

  {#if visualizations.plot1}
    <h2>Visualizations</h2>
    <img src={visualizations.plot1} alt="Plot 1" style="width: 300px; height: auto;" />
    <img src={visualizations.plot2} alt="Plot 2" style="width: 300px; height: auto;" />
    <img src={visualizations.plot3} alt="Plot 3" style="width: 300px; height: auto;" />
  {/if}
</div>
