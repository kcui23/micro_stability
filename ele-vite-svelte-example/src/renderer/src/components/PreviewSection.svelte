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
      <div class="table-container">
        <table>
          {#each filteredContent as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
        </table>
      </div>
    </div>
  {/if}

  {#if groupingsFile}
    <div class="preview">
      <h2>Preview of Groupings File</h2>
      <p>Dimensions: {groupingsDimensions.rows} rows, {groupingsDimensions.columns}</p>
      <div class="table-container">
        <table>
          {#each groupingsContentPreview as row}
            <tr>
              {#each row as cell}
                <td>{cell}</td>
              {/each}
            </tr>
          {/each}
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

  .table-container {
    margin-top: 10px;
    background-color: #ffffff;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 10px;
    box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
    overflow-x: auto;
  }

  table {
    width: 100%;
    border-collapse: collapse;
  }

  th,
  td {
    padding: 8px;
    text-align: left;
    border-bottom: 1px solid #ddd;
  }

  td {
    font-size: 0.9rem;
  }
</style>
