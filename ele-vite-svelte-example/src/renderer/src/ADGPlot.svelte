<script>
    import { onMount, createEventDispatcher } from 'svelte';
    import { writable } from 'svelte/store';
    
    export let steps;
    export let currentStep;
    export let setCurrentStep;
    
    const dispatch = createEventDispatcher();
    
    let stepStatus = writable(steps.reduce((acc, step) => {
        acc[step] = step === 'Raw data' ? 'Enabled' : 'Disabled';
        return acc;
    }, {}));
    
    let subOperations = writable({
        'Raw data': ['Set Random Seed'],
        'Data Perturbation': ['Apply Threshold', 'Additional Option 1', 'Additional Option 2'],
        'Model Perturbation': ['Select Method'],
        'Prediction Evaluation Metric': ['View Results'],
        'Stability Metric': ['View Stability Plot', 'Run Shuffled Analysis']
    });
    
    let selectedOperations = writable(Object.fromEntries(steps.map(step => [step, []])));
    

    function toggleStepStatus(step) {
        if (step === 'Raw data') return; // Raw data is always enabled
        stepStatus.update(status => {
            status[step] = status[step] === 'Enabled' ? 'Disabled' : 'Enabled';
            return status;
        });
    }

    function selectStep(step) {
        if ($stepStatus[step] === 'Enabled') {
            setCurrentStep(step);
            dispatch('stepSelected', { step });
        }
    }

    function toggleOperation(step, operation) {
        selectedOperations.update(selections => {
            if (selections[step].includes(operation)) {
                selections[step] = selections[step].filter(op => op !== operation);
            } else {
                selections[step] = [...selections[step], operation];
            }
            return selections;
        });
        dispatch('operationsChanged', { step, operations: $selectedOperations[step] });
    }

    onMount(() => {
        // Any initialization logic if needed
    });
</script>

<div class="adg-container">
    <h2>Steps</h2>
    <div class="steps-list">
        {#each steps as step}
            <div class="step-item" class:disabled={$stepStatus[step] === 'Disabled'}>
                <button 
                    class="step-button" 
                    on:click={() => selectStep(step)}
                    disabled={$stepStatus[step] === 'Disabled'}
                >
                    {step}
                </button>
                {#if step !== 'Raw data'}
                    <div class="status-toggles">
                        <button 
                            class="status-toggle disabled" 
                            class:inactive={$stepStatus[step] === 'Enabled'}
                            on:click={() => $stepStatus[step] === 'Enabled' ? toggleStepStatus(step) : null}
                        >
                            Disabled
                        </button>
                        <button 
                            class="status-toggle enabled" 
                            class:inactive={$stepStatus[step] === 'Disabled'}
                            on:click={() => $stepStatus[step] === 'Disabled' ? toggleStepStatus(step) : null}
                        >
                            Enabled
                        </button>
                    </div>
                {/if}
            </div>
            {#if step === currentStep}
                <div class="sub-operations">
                    {#each $subOperations[step] as operation}
                        <label class="operation-checkbox">
                            <input 
                                type="checkbox" 
                                checked={$selectedOperations[step].includes(operation)}
                                on:change={() => toggleOperation(step, operation)}
                            />
                            {operation}
                        </label>
                    {/each}
                </div>
            {/if}
        {/each}
    </div>
</div>

<style>
    .adg-container {
        width: 100%;
        max-width: 300px;
        padding: 20px;
        background-color: #f5f5f5;
        border-radius: 8px;
    }

    h2 {
        margin-bottom: 20px;
        color: #333;
        font-size: 1.5em;
    }

    .steps-list {
        display: flex;
        flex-direction: column;
        gap: 10px;
    }

    .step-item {
        display: flex;
        flex-direction: column;
        transition: all 0.3s ease;
    }

    .step-item.disabled {
        opacity: 0.6;
        transform: scale(0.95);
    }

    .step-button {
        padding: 10px;
        background-color: #fff;
        border: 1px solid #ddd;
        border-radius: 4px;
        cursor: pointer;
        transition: background-color 0.3s ease;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
    }

    .step-button:hover:not(:disabled) {
        background-color: #e9e9e9;
    }

    .step-button:disabled {
        cursor: not-allowed;
    }

    .status-toggles {
        display: flex;
        justify-content: space-between;
        margin-top: 5px;
    }

    .status-toggle {
        flex: 1;
        padding: 5px;
        border: none;
        cursor: pointer;
        transition: all 0.3s ease;
        font-size: 0.85em; 
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }


    .status-toggle.disabled {
        background-color: #ffcccc;
        color: #cc0000;
    }

    .status-toggle.enabled {
        background-color: #ccffcc;
        color: #006600;
    }

    .status-toggle.inactive {
        background-color: #e0e0e0;
        color: #7d7d7d;
    }

    .status-toggle.active {
        font-weight: bold;
    }

    .status-toggle.disabled:hover {
        background-color: #ff9999;
    }

    .status-toggle.enabled:hover {
        background-color: #99ff99;
    }

    .sub-operations {
        margin-top: 10px;
        padding-left: 20px;
    }

    .operation-checkbox {
        display: block;
        margin-bottom: 5px;
    }

    .operation-checkbox input {
        margin-right: 5px;
}
</style>
