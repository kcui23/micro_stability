<script>
  let files = [];
  let visualizations = { plot1: '', plot2: '', plot3: '' };
  let selectedMethod = '';
  let groupingsFile = null;

  const handleFileChange = (event) => {
    files = event.target.files;
  };

  const handleGroupingsChange = (event) => {
    groupingsFile = event.target.files[0];
  };

  const handleMethodChange = (event) => {
    selectedMethod = event.target.value;
  };

  const handleSubmit = async () => {
    const file = files[0];
    const groupings = groupingsFile;

    const asvReader = new FileReader();
    const groupingsReader = new FileReader();

    asvReader.onload = () => {
      const asvContent = asvReader.result;

      groupingsReader.onload = async () => {
        const groupingsContent = groupingsReader.result;

        const response = await fetch(`http://localhost:8000/process?method=${selectedMethod}`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify({
            asv: asvContent,
            groupings: groupingsContent
          })
        });

        const result = await response.json();
        visualizations = {
          plot1: `data:image/png;base64,${result.plot1}`,
          plot2: `data:image/png;base64,${result.plot2}`,
          plot3: `data:image/png;base64,${result.plot3}`
        };
      };

      groupingsReader.readAsText(groupings);
    };

    asvReader.readAsText(file);
  };
</script>

<div>
  <h1>Upload ASV Dataset</h1>
  <input type="file" accept=".tsv" on:change={handleFileChange} />
  <input type="file" accept=".tsv" on:change={handleGroupingsChange} />

  {#if files.length > 0 && groupingsFile}
    <p>File: {files[0].name}</p>
    <p>Groupings File: {groupingsFile.name}</p>
    
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
