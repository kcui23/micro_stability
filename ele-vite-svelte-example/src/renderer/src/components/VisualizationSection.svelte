<script>
  import InteractiveValcano from './InteractiveValcano.svelte'
  import { fade, scale } from 'svelte/transition'

  export let visualizations
  export let zoomImage
  export let isCalculating
  export let showAllPlots
  export let showDetailedPlots
  export let isStatic
  export let toggleView
  export let selectedPointsList
  export let selectedMethod
  export let isSubmitted
  export let zoomedImage

  const saveImage = (imageData, name) => {
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    const img = new Image();
    img.onload = () => {
      canvas.width = img.width;
      canvas.height = img.height;
      ctx.fillStyle = 'white';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0);
      const link = document.createElement('a');
      link.href = canvas.toDataURL('image/png');
      link.download = name + '.png';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    };
    img.src = imageData;
  }

  const removePoint = (pointToRemove) => {
    selectedPoints.update(points => points.filter(point => point !== pointToRemove));
  };

  const downloadSelectedPoints = async () => {
    const data = selectedPointsList.map(point => ({
      name: point.name,
      x: point.x,
      y: point.y
    }));

    const response = await fetch('http://localhost:8000/download_selected_points', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ points: data })
    });

    if (response.ok) {
      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = 'selected_points.tsv';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } else {
      console.error('Error downloading selected points:', response.statusText);
    }
  };
</script>

<div class="visualizations-section" hidden={!showAllPlots && !isSubmitted}>
  <div class="visualization-header">
    <h2>Visualizations</h2>
    <div class="view-toggle">
      <button class:activeToggle={isStatic} on:click={() => toggleView('static')}>
        Static
      </button>
      <button class:activeToggle={!isStatic} on:click={() => toggleView('interactive')}>
        Interactive
      </button>
    </div>
  </div>

  {#if isCalculating}
    <div class="loader">
      <p>Loading...</p>
    </div>
  {:else if showAllPlots}
    <div class="card-container">
      <div class="card">
        <div class="card-header">
          <h3>Overlap Visualizations</h3>
        </div>
        <div class="card-content">
          {#if isStatic}
            <img
              src={visualizations.overlap_volcano}
              alt="Overlap Volcano Plot"
              on:click={() => zoomImage('overlap_volcano')}
            />
            <img
              src={visualizations.overlap_pvalue_distribution}
              alt="Overlap P-value Distribution"
              on:click={() => zoomImage('overlap_pvalue_distribution')}
            />
          {:else}
            <InteractiveValcano />
          {/if}
        </div>
      </div>

      {#if selectedPointsList.length > 0}
        <div class="floating-card">
          <div class="card-header">
            <h3>Selected Points</h3>
          </div>
          <div class="card-content">
            <ul>
              {#each selectedPointsList as point (point.name)}
                <li>
                  <span>{point.name}</span>
                  <button on:click={() => removePoint(point)}>Remove</button>
                </li>
              {/each}
            </ul>
            <button class="download-button" on:click={downloadSelectedPoints}
              >Download Points</button
            >
          </div>
        </div>
      {/if}

      {#if showDetailedPlots}
        {#each ['deseq2', 'aldex2', 'edger', 'maaslin2', 'metagenomeseq'] as method}
          <div class="card">
            <div class="card-header">
              <h3>{method.toUpperCase()} Plots</h3>
            </div>
            <div class="card-content">
              {#if isStatic}
                {#each [1, 2] as plotNumber}
                  <img
                    src={visualizations[`${method}_plot${plotNumber}`]}
                    alt={`${method.toUpperCase()} Plot ${plotNumber}`}
                    on:click={() => zoomImage(`${method}_plot${plotNumber}`)}
                  />
                {/each}
              {:else}
                <div class="interactive-placeholder">
                  Interactive {method.toUpperCase()} Plots coming soon...
                </div>
              {/if}
            </div>
          </div>
        {/each}
      {/if}
    </div>

    <button
      class="expand-collapse-button"
      on:click={() => showDetailedPlots = !showDetailedPlots}
    >
      {showDetailedPlots ? 'Collapse Details' : 'Show More Details'}
    </button>
  {:else if selectedMethod}
    <div class="card">
      <div class="card-header">
        <h3>{selectedMethod.toUpperCase()} Plots</h3>
      </div>
      <div class="card-content">
        {#if isStatic}
          {#each [1, 2, 3] as plotNumber}
            <img
              src={visualizations[`${selectedMethod}_plot${plotNumber}`]}
              alt={`${selectedMethod.toUpperCase()} Plot ${plotNumber}`}
              on:click={() => zoomImage(`${selectedMethod}_plot${plotNumber}`)}
            />
          {/each}
        {:else}
          <div class="interactive-placeholder">
            Interactive {selectedMethod.toUpperCase()} Plots coming soon...
          </div>
        {/if}
      </div>
    </div>
  {/if}

  {#if zoomedImage}
    <div class="zoomed-image-container" on:click={() => zoomImage(null)} transition:fade>
      <div class="zoomed-image-content" on:click|stopPropagation transition:scale>
        <img src={visualizations[zoomedImage]} alt="Zoomed Image" class="bounded-image" />
        <button
          class="download-button"
          on:click|stopPropagation={() => saveImage(visualizations[zoomedImage], zoomedImage)}
        >
          Save Image
        </button>
      </div>
    </div>
  {/if}
</div>

<style>
  .visualizations-section {
    margin-top: 2rem;
  }

  .visualization-header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 1rem;
  }

  .view-toggle {
    display: flex;
  }

  .view-toggle button {
    background-color: #f0f0f0;
    border: 1px solid #ddd;
    border-radius: 4px;
    cursor: pointer;
    width: 100px;
    transition: background-color 0.3s ease, color 0.3s ease;
  }

  .view-toggle button:hover {
    background-color: #e0e0e0;
  }

  .view-toggle button.activeToggle {
    background-color: #007bff;
    color: white;
  }

  .view-toggle button.activeToggle:hover {
    background-color: #0073ee;
  }

  .card-container {
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
  }

  .card {
    flex: 1 1 calc(50% - 0.5rem);
    border: 1px solid #ddd;
    border-radius: 4px;
    overflow: hidden;
    background-color: #fff;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    position: relative;
  }

  .card-header {
    background-color: #f5f5f5;
    padding: 1rem;
    border-bottom: 1px solid #ddd;
  }

  .card-header h3 {
    margin: 0;
    font-size: 1.2rem;
  }

  .card-content {
    padding: 1rem;
    display: flex;
    flex-wrap: wrap;
    gap: 1rem;
  }

  .card-content img {
    max-width: calc(50% - 0.5rem);
    height: auto;
    cursor: pointer;
  }

  .floating-card {
    position: fixed;
    right: 1rem;
    top: 1rem;
    width: 250px;
    border: 1px solid #ddd;
    border-radius: 4px;
    background-color: #fff;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    z-index: 100;
  }

  .floating-card .card-content {
    max-height: 300px;
    overflow-y: auto;
  }

  .floating-card ul {
    list-style-type: none;
    padding: 0;
    margin: 0;
  }

  .floating-card li {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0.5rem 0;
  }

  .floating-card li button {
    background-color: #e74c3c;
    color: white;
    border: none;
    cursor: pointer;
    padding: 0.2rem 0.5rem;
  }

  .expand-collapse-button {
    display: block;
    margin: 1rem auto;
    padding: 0.5rem 1rem;
    background-color: #007bff;
    color: white;
    border: none;
    cursor: pointer;
    transition: background-color 0.3s;
  }

  .expand-collapse-button:hover {
    background-color: #0056b3;
  }

  .interactive-placeholder {
    width: 100%;
    height: 300px;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: #f0f0f0;
    border-radius: 4px;
    font-style: italic;
    color: #666;
  }

  .zoomed-image-container {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background-color: rgba(0, 0, 0, 0.8);
    display: flex;
    justify-content: center;
    align-items: center;
    z-index: 1000;
  }

  .zoomed-image-content {
    position: relative;
    max-width: 90vw;
    max-height: 90vh;
    background-color: white;
    padding: 20px 20px 45px 20px;
    box-shadow: 0 0 20px rgba(0, 0, 0, 0.3);
    border-radius: 8px;
    display: flex;
    flex-direction: column;
  }

  .zoomed-image-content img {
    max-width: 100vw;
    max-height: 75vh;
    object-fit: contain;
  }

  .download-button {
    position: absolute;
    width: 20%;
    bottom: 10px;
    left: 50%;
    transform: translateX(-50%);
  }
  .loader {
    display: flex;
    justify-content: center;
    align-items: center;
    height: 200px;
  }
</style>