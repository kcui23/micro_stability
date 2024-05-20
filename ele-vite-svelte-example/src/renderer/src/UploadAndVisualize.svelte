<script>
  let files = [];
  let visualizations = { plot1: '', plot2: '', plot3: '' };

  const handleFileChange = async (event) => {
    files = event.target.files;
    const file = files[0];
    
    const reader = new FileReader();
    reader.onload = async () => {
      const csvContent = reader.result;

      const response = await fetch('http://localhost:8000/upload', {
        method: 'POST',
        headers: {
          'Content-Type': 'text/plain'
        },
        body: csvContent
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
  <h1>Upload Iris Dataset</h1>
  <input type="file" accept=".csv" on:change={handleFileChange} />

  {#if files.length > 0}
    <p>File: {files[0].name}</p>
  {/if}

  {#if visualizations.plot1}
    <h2>Visualizations</h2>
    <img src={visualizations.plot1} alt="Sepal Length by Species" />
    <img src={visualizations.plot2} alt="Petal Length Density by Species" />
    <img src={visualizations.plot3} alt="Petal Length vs Petal Width" />
  {/if}
</div>
