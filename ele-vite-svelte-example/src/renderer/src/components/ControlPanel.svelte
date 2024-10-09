<script>
    import { singleSelectOperations, selectedColorStep, scatterPlotColors } from '../store.js';
    
    $: operations = $singleSelectOperations || {};

    function isLightColor(color) {
        const hex = color.replace('#', '');
        const r = parseInt(hex.substr(0, 2), 16);
        const g = parseInt(hex.substr(2, 2), 16);
        const b = parseInt(hex.substr(4, 2), 16);
        const brightness = (r * 299 + g * 587 + b * 114) / 1000;
        return brightness > 128;
    }

</script>

<h2>Color Panel</h2>
<select bind:value={$selectedColorStep}>
    {#each Object.keys(operations) as step}
        <option value={step}>{step}</option>
    {/each}
</select>

<div class="color-options">
    {#each operations[$selectedColorStep] as operation, index}
        {@const bgColor = $scatterPlotColors[$selectedColorStep][index]}
        <div class="color-item" 
             style="background-color: {bgColor}; color: {isLightColor(bgColor) ? 'black' : 'white'}">
            {operation}
        </div>
    {/each}
</div>

<style>
    h2 {
        margin-top: 0;
    }

    select {
        width: 100%;
        margin-bottom: 1rem;
    }

    .color-options {
        display: flex;
        flex-wrap: wrap;
        gap: 0.5rem;
    }

    .color-item {
        display: inline-block;
        padding: 0.25rem 0.5rem;
        border-radius: 4px;
        color: white;
        font-size: 0.9rem;
        text-align: center;
        margin-bottom: 0.25rem;
    }
</style>