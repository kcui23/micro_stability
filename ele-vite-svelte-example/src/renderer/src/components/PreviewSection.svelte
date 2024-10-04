<script>
  export let filteredContent
  export let groupingsContentPreview
  export let filteredDimensions
  export let groupingsDimensions
  export let asvFiles
  export let groupingsFile
</script>

<div class="preview-section">
  {#if asvFiles.length > 0}
    <div class="preview">
      <h2>Preview of ASV File</h2>
      <p>Dimensions: {filteredDimensions.rows} rows, {filteredDimensions.columns} columns</p>
      <div class="tidy-table-container">
        <table class="tidy-table">
          <thead>
            <tr>
              {#each filteredContent[0] as header}
                <th>{header}</th>
              {/each}
            </tr>
          </thead>
          <tbody>
            {#each filteredContent.slice(1) as row}
              <tr>
                {#each row as cell}
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
      <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns} columns</p>
      <div class="tidy-table-container">
        <table class="tidy-table">
          <thead>
            <tr>
              {#each groupingsContentPreview[0] as header}
                <th>{header}</th>
              {/each}
            </tr>
          </thead>
          <tbody>
            {#each groupingsContentPreview.slice(1) as row}
              <tr>
                {#each row as cell}
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
