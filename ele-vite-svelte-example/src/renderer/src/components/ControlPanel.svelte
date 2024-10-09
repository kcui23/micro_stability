<script>
    import { singleSelectOperations, selectedStep } from '../store.js';
    
    let colors = {
        'Filtering': ['#FF5733', '#33FF57', '#3357FF', '#FF33F1'],
        'Zero-Handling': ['#FFC300', '#DAF7A6', '#FF5733'],
        'Normalization': ['#C70039', '#900C3F', '#581845', '#FFC300', '#DAF7A6'],
        'Transformation': ['#FF5733', '#C70039', '#900C3F', '#581845'],
        'Model Perturbation': ['#FFC300', '#DAF7A6', '#FF5733', '#C70039', '#900C3F'],
        'Stability Metric': ['#581845', '#FFC300', '#DAF7A6']
    };
    $: operations = $singleSelectOperations || {};
    $: if ($selectedStep) {
        console.log($selectedStep);
    }

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
<select bind:value={$selectedStep}>
    {#each Object.keys(operations) as step}
        <option value={step}>{step}</option>
    {/each}
</select>

<div class="color-options">
    {#each operations[$selectedStep] as operation, index}
        {@const bgColor = colors[$selectedStep][index]}
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