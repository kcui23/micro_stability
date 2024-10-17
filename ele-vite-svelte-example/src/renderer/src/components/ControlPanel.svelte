<script>
    import { singleSelectOperations, 
        selectedColorStep, 
        scatterPlotColors,
        colorStatus } from '../store.js';
    
    $: operations = $singleSelectOperations || {};

    function isLightColor(color) {
        const hex = color.replace('#', '');
        const r = parseInt(hex.substr(0, 2), 16);
        const g = parseInt(hex.substr(2, 2), 16);
        const b = parseInt(hex.substr(4, 2), 16);
        const brightness = (r * 299 + g * 587 + b * 114) / 1000;
        return brightness > 128;
    }

    function handleButtonClick(operation) {
        let tmp = $colorStatus
        if (tmp[$selectedColorStep].includes(operation)) {
            colorStatus.update(d => {
                d[$selectedColorStep] = d[$selectedColorStep].filter(op => op !== operation);
                return d;
            })
        }
        else {
            colorStatus.update(d => {
                d[$selectedColorStep] = [...d[$selectedColorStep], operation];
                return d;
            })
        }
        console.log("updated $colorStatus:", $colorStatus);
    }
</script>

<div class="panel-container">
    <h2>Color Panel</h2>
    <select bind:value={$selectedColorStep}>
        {#each Object.keys(operations) as step}
            <option value={step}>{step}</option>
        {/each}
    </select>

    <div class="color-options">
        {#each operations[$selectedColorStep] as operation, index}
            {@const bgColor = $scatterPlotColors[$selectedColorStep][index]}
            <button class="color-item" 
                    style="background-color: {bgColor}; 
                    color: {isLightColor(bgColor) ? 'black' : 'white'};
                    opacity: {$colorStatus[$selectedColorStep].includes(operation) ? 1 : 0.5}"
                    on:click={() => handleButtonClick(operation)}>
                {operation}
            </button>
        {/each}
    </div>
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
        height: 2rem;
        display: inline-block;
        padding: 0.25rem 0.5rem;
        border-radius: 4px;
        color: white;
        font-size: 0.9rem;
        text-align: center;
        margin-bottom: 0.25rem;
        border: none;
    }

    .color-item:hover {
        opacity: 0.8;
    }

    .panel-container {
        flex: 0 0 80%;
        height: 100%;
        background-color: #fff;
        overflow: auto; /* Allow scrolling if the tree content overflows */
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        padding: 10px;
        border-radius: 8px;
        box-sizing: border-box;
    }
</style>
