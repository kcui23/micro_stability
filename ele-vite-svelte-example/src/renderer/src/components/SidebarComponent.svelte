<script>
  import { onMount, createEventDispatcher } from 'svelte';
  import {
    stepStatus,
    params,
    subOperations,
    selectedOperations,
    openMenus,
    singleSelectOperations,
    crossStepMutuallyExclusiveOptions
  } from '../store.js';
  import { fade, slide } from 'svelte/transition';

  export let steps;
  export let currentStep;
  export let setCurrentStep;

  const dispatch = createEventDispatcher();

  // Add operation name mapping
  const operationDisplayNames = {
    'Model Perturbation': {
      'deseq2': 'DESeq2',
      'edger': 'edgeR',
      'maaslin2': 'MaAsLin2',
      'aldex2': 'ALDEx2',
      'metagenomeseq': 'metagenomeSeq'
    }
  };

  function getDisplayName(step, operation) {
    if (operationDisplayNames[step] && operationDisplayNames[step][operation]) {
      return operationDisplayNames[step][operation];
    }
    return operation;
  }

  async function downloadCodes() {
    try {
      const response = await fetch('http://localhost:8000/download_code', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          selectedOperations: $selectedOperations,
          params: $params
        })
      });
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = 'r_code.R';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error('Error downloading the code:', error);
    }
  }

  stepStatus.set(
    steps.reduce((acc, step) => {
      acc[step] = step === 'Raw data' ? 'Enabled' : 'Disabled';
      return acc;
    }, {})
  );

  stepStatus.subscribe(value => {
    console.log("stepStatus:", value);
  });

  selectedOperations.set(Object.fromEntries(steps.map((step) => [step, []])));
  console.log("selectedOperations:", $selectedOperations);
  openMenus.set(Object.fromEntries(steps.map((step) => [step, false])));

  function toggleStepStatus(step) {
    if (step === 'Raw data') return; // Raw data is always enabled
    stepStatus.update((status) => {
      const newStatus = status[step] === 'Enabled' ? 'Disabled' : 'Enabled';
      status[step] = newStatus;

      // If the new status is 'Disabled', close the menu and update currentStep
      if (newStatus === 'Disabled') {
        openMenus.update(menus => {
          menus[step] = false;
          return menus;
        });
        if (currentStep === step) {
          setCurrentStep(null);
        }
      }

      // If the new status is 'Enabled', open the menu and update currentStep
      if (newStatus === 'Enabled') {
        openMenus.update((menus) => {
          menus[step] = true;
          // set all the other menus[step] to false
          Object.keys(menus).forEach((key) => {
            if (key !== step) menus[key] = false;
          });
          return menus;
        });
        setCurrentStep(step);
      }
      
      return status;
    });
  }

  function selectStep(step) {
    if ($stepStatus[step] === 'Enabled') {
      setCurrentStep(step);
      toggleMenu(step);
    }
  }

  function toggleOperation(step, operation) {
    selectedOperations.update((selections) => {
      const isSingleSelect = $singleSelectOperations[step] && $singleSelectOperations[step].includes(operation);
      if (isSingleSelect) {
        const existingSingleSelect = $singleSelectOperations[step].find(op => selections[step].includes(op));
        if (existingSingleSelect) {
          selections[step] = selections[step].filter(op => op !== existingSingleSelect);
        }
        selections[step].push(operation);
      } else {
        if (selections[step].includes(operation)) {
          selections[step] = selections[step].filter(op => op !== operation);
        } else {
          selections[step] = [...selections[step], operation];
        }
      }

      if (crossStepMutuallyExclusiveOptions[step] && crossStepMutuallyExclusiveOptions[step][operation]) {
        Object.entries(crossStepMutuallyExclusiveOptions[step][operation]).forEach(([mutexStep, mutexOps]) => {
          selections[mutexStep] = selections[mutexStep].filter(op => !mutexOps.includes(op));
        });
      }

      dispatch('pathChange', selections);
      return selections;
    });
  }

  function isOptionDisabled(step, operation) {
  for (const [otherStep, options] of Object.entries(crossStepMutuallyExclusiveOptions)) {
    for (const [otherOp, mutexInfo] of Object.entries(options)) {
      if ($selectedOperations[otherStep].includes(otherOp)) {
        for (const [mutexStep, mutexOps] of Object.entries(mutexInfo)) {
          if (step === mutexStep && mutexOps.includes(operation)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

  function toggleMenu(step) {
    openMenus.update((menus) => {
      Object.keys(menus).forEach((key) => {
        if (key !== step) menus[key] = false;
      });
      menus[step] = !menus[step];
      return menus;
    });
  }

  onMount(() => {
    openMenus.update((menus) => {
      menus[currentStep] = true;
      return menus;
    });
  });
</script>

<div class="adg-container">
  <h2>Analysis Steps</h2>
  <div class="steps-list">
    {#each steps as step, index}
      <div
        class="step-item"
        class:disabled={$stepStatus[step] === 'Disabled'}
        class:active={step === currentStep}
      >
        <button
          class="step-button"
          on:click={() => selectStep(step)}
          disabled={$stepStatus[step] === 'Disabled'}
        >
          <span class="step-number" class:disabled={$stepStatus[step] === 'Disabled'}>
            {index + 1}
          </span>
          <span class="step-text">{step}</span>
          <span class="dropdown-indicator" class:open={$openMenus[step]}>▼</span>
        </button>
        <!-- {#if step !== 'Raw data'}
          <div class="status-toggles">
            <button
              class="status-toggle disabled"
              class:inactive={$stepStatus[step] === 'Enabled'}
              on:click={() => ($stepStatus[step] === 'Enabled' ? toggleStepStatus(step) : null)}
            >
              {$stepStatus[step] === 'Enabled' ? 'Disable' : 'Disabled'}
            </button>
            <button
              class="status-toggle enabled"
              class:inactive={$stepStatus[step] === 'Disabled'}
              on:click={() => ($stepStatus[step] === 'Disabled' ? toggleStepStatus(step) : null)}
            >
              {$stepStatus[step] === 'Disabled' ? 'Enable' : 'Enabled'}
            </button>
          </div>
        {/if} -->
        {#if $openMenus[step]}
          <div class="sub-operations" transition:slide={{ duration: 300 }}>
            {#if $singleSelectOperations[step]}
              <!-- Single Select (Radio) Operations -->
              {#each $singleSelectOperations[step] as operation}
                <label class="operation-radio" transition:fade={{ duration: 200 }}>
                  <input
                    type="radio"
                    name={step + "-singleSelect"}
                    checked={$selectedOperations[step].includes(operation)}
                    on:change={() => toggleOperation(step, operation)}
                    disabled={isOptionDisabled(step, operation)}
                  />
                  {getDisplayName(step, operation)}
                </label>
              {/each}
              {#if $subOperations[step].length > $singleSelectOperations[step].length}
                <hr /> <!-- Divider between radio and checkbox selections -->
              {/if}
            {/if}
            <!-- Multiple Select (Checkbox) Operations -->
            
            {#each $subOperations[step] as operation}
              {#if !$singleSelectOperations[step] || !$singleSelectOperations[step].includes(operation)}
                <label class="operation-checkbox" transition:fade={{ duration: 200 }}>
                  <input
                    type="checkbox"
                    checked={$selectedOperations[step].includes(operation)}
                    on:change={() => toggleOperation(step, operation)}
                    disabled={isOptionDisabled(step, operation)}
                  />
                  {getDisplayName(step, operation)}
                </label>
              {/if}
            {/each}
          </div>
        {/if}
      </div>
    {/each}
  </div>
  <div class="button-group">
    <button class="sidebar-button" on:click={() => downloadCodes()}>
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round">
        <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" stroke-width="3"></path>
        <polyline points="7 10 12 15 17 10" stroke-width="2" color="#a3a3a3"></polyline>
        <line x1="12" y1="15" x2="12" y2="3" stroke-width="2" color="#a3a3a3"></line>
      </svg>
      Download codes
    </button>
    <button class="sidebar-button">
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round">
        <path d="M22 3H2l8 9.46V19l4 2v-8.54L22 3z" stroke-width="3"/>
      </svg>
      Select ASVs
    </button>
  </div>
</div>

<style>
  .button-group {
    display: flex;
    flex-direction: column;
    gap: 10px;
    margin-top: 20px;
  }

  .sidebar-button {
    font-size: 14px;
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 5px;
  }

  .sidebar-button svg {
    width: 16px;
    height: 16px;
  }

  .adg-container {
    width: 100%;
    max-width: 300px;
    padding: 20px;
    background-color: #f5f5f5;
    border-radius: 8px;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
  }

  h2 {
    margin-bottom: 20px;
    color: #333;
    font-size: 1.5em;
    text-align: center;
  }

  .steps-list {
    display: flex;
    flex-direction: column;
    gap: 15px;
  }

  .step-item {
    display: flex;
    flex-direction: column;
    transition: all 0.3s ease;
    background-color: #fff;
    border-radius: 8px;
    overflow: hidden;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  }

  .step-item.disabled {
    opacity: 0.6;
  }

  .step-item.active {
    box-shadow: 0 0 0 2px #007bff;
    filter: drop-shadow(0 0 3px #bdddff);
  }

  .step-button {
    padding: 10px;
    background-color: #fff;
    border: none;
    cursor: pointer;
    transition: background-color 0.3s ease;
    display: flex;
    align-items: center;
    justify-content: space-between;
    text-align: left;
    width: 100%;
    font-size: 1em;
    min-height: 60px;
  }

  .step-number {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    flex: 0 0 24px;
    width: 24px;
    height: 24px;
    background-color: #007bff;
    color: #fff;
    border-radius: 50%;
    font-size: 0.8em;
    font-weight: bold;
    transition: background-color 0.3s ease;
  }

  .step-number.disabled {
    background-color: #999999;
  }

  .step-text {
    flex: 1;
    margin: 0 10px;
    word-break: break-word;
  }

  .dropdown-indicator {
    flex: 0 0 20px;
    text-align: center;
  }

  .step-button:hover:not(:disabled) {
    background-color: #e9e9e9;
  }

  .step-button:disabled {
    cursor: not-allowed;
  }

  /* .status-toggles {
    display: flex;
    justify-content: space-between;
    padding: 5px 10px;
    background-color: #f0f0f0;
  }

  .status-toggle {
    flex: 1;
    padding: 5px;
    border: none;
    cursor: pointer;
    transition: all 0.3s ease;
    font-size: 0.8em;
    border-radius: 4px;
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

  .step-item:not(.disabled) .status-toggle.disabled:hover {
    background-color: #fdbbbb;
    color: #cc0000;
  }

  .step-item.disabled .status-toggle.enabled:hover {
    background-color: #66ff66;
    color: #006600;
  } */

  .sub-operations {
    padding: 10px;
    background-color: #f9f9f9;
    border-top: 1px solid #e0e0e0;
  }

  hr {
    border: none;
    border-top: 1px solid #e0e0e0;
    margin: 10px 0;
  }

  .operation-checkbox,
  .operation-radio {
    display: flex;
    align-items: center;
    margin-bottom: 5px;
    font-size: 0.9em;
  }

  .operation-checkbox input,
  .operation-radio input {
    margin-right: 5px;
  }

  .operation-checkbox input:disabled,
  .operation-radio input:disabled {
    opacity: 0.5;
    cursor: not-allowed;
  }

  .operation-checkbox:has(input:disabled),
  .operation-radio:has(input:disabled) {
    color: #999;
  }
  
  .dropdown-indicator {
    flex: 0 0 20px;
    text-align: center;
    transition: transform 0.3s ease;
  }

  .dropdown-indicator.open {
    transform: rotate(180deg);
  }
</style>
