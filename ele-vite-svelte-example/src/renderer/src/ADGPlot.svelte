<script>
	import { onMount, createEventDispatcher } from 'svelte';
	import { writable } from 'svelte/store';
	import { fade, slide } from 'svelte/transition';

	export let steps;
	export let currentStep;
	export let setCurrentStep;

	const dispatch = createEventDispatcher();

	let stepStatus = writable(
		steps.reduce((acc, step) => {
			acc[step] = step === 'Raw data' ? 'Enabled' : 'Disabled';
			return acc;
		}, {})
	);

	let subOperations = writable({
		'Raw data': ['Set Random Seed'],
		'Data Perturbation': ['Apply Threshold', 'Additional Option 1', 'Additional Option 2'],
		'Model Perturbation': ['Select Method'],
		'Prediction Evaluation Metric': ['View Results'],
		'Stability Metric': ['View Stability Plot', 'Run Shuffled Analysis']
	});

	let selectedOperations = writable(Object.fromEntries(steps.map((step) => [step, []])));
	let openMenus = writable(Object.fromEntries(steps.map((step) => [step, false])));

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
				openMenus.update(menus => {
					menus[step] = true;
					// set all the other menus[step] to false
					Object.keys(menus).forEach(key => {
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
			dispatch('stepSelected', { step });
			toggleMenu(step);
		}
	}

	function toggleOperation(step, operation) {
		selectedOperations.update((selections) => {
			if (selections[step].includes(operation)) {
				selections[step] = selections[step].filter((op) => op !== operation);
			} else {
				selections[step] = [...selections[step], operation];
			}
			return selections;
		});
		dispatch('operationsChanged', { step, operations: $selectedOperations[step] });
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
					<!-- <span class="step-number">{index + 1}</span> -->
					<span class="step-number" class:disabled={$stepStatus[step] === 'Disabled'}>{index + 1}</span>
					<span class="step-text">{step}</span>
					<span class="dropdown-indicator">{$openMenus[step] ? '▲' : '▼'}</span>
				</button>
				{#if step !== 'Raw data'}
					<div class="status-toggles">
						<button
							class="status-toggle disabled"
							class:inactive={$stepStatus[step] === 'Enabled'}
							on:click={() => ($stepStatus[step] === 'Enabled' ? toggleStepStatus(step) : null)}
						>
							Disabled
						</button>
						<button
							class="status-toggle enabled"
							class:inactive={$stepStatus[step] === 'Disabled'}
							on:click={() => ($stepStatus[step] === 'Disabled' ? toggleStepStatus(step) : null)}
						>
							Enabled
						</button>
					</div>
				{/if}
				{#if $openMenus[step]}
					<div class="sub-operations" transition:slide={{ duration: 300 }}>
						{#each $subOperations[step] as operation}
							<label class="operation-checkbox" transition:fade={{ duration: 200 }}>
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
			</div>
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

	.status-toggles {
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
	}

	.sub-operations {
		padding: 10px;
		background-color: #f9f9f9;
		border-top: 1px solid #e0e0e0;
	}

	.operation-checkbox {
		display: flex;
		align-items: center;
		margin-bottom: 5px;
		font-size: 0.9em;
	}

	.operation-checkbox input {
		margin-right: 5px;
	}
</style>