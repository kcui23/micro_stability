<script>
  export let asvPreview
  export let groupingsPreview
  export let asvDimensions
  export let groupingsDimensions
  export let asvFiles
  export let groupingsFile
</script>

<div class="preview-section">
  {#if asvFiles.length > 0}
    <div class="preview">
      <h2>Preview of ASV File</h2>
      <p>Dimensions: {asvDimensions[0]} rows, {asvDimensions[1]} columns</p>
      <div class="tidy-table-container">
        <table class="tidy-table">
          <thead>
            <tr>
              {#each Object.keys(asvPreview[0]) as header}
                <th>{header}</th>
              {/each}
            </tr>
          </thead>
          <tbody>
            {#each asvPreview as row}
              <tr>
                {#each Object.values(row) as cell}
                  <td>{cell}</td>
                {/each}
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    </div>
  {/if}

  {#if groupingsFile}
    <div class="preview">
      <h2>Preview of Groupings File</h2>
      <p>Dimensions: {groupingsDimensions[0]} rows, {groupingsDimensions[1]} columns</p>
      <div class="tidy-table-container">
        <table class="tidy-table">
          <thead>
            <tr>
              {#each Object.keys(groupingsPreview[0]) as header}
                <th>{header}</th>
              {/each}
            </tr>
          </thead>
          <tbody>
            {#each groupingsPreview as row}
              <tr>
                {#each Object.values(row) as cell}
                  <td>{cell}</td>
                {/each}
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    </div>
  {/if}
</div>

<style>
  .preview-section {
    margin-top: 20px;
  }

  .preview {
    margin-top: 20px;
  }

</style>
