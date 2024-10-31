// store.js
import { readable, writable } from 'svelte/store';

export const selectedPoints = writable([]);
export const currentPath = writable(['Raw data', 'Low Abundance Filtering', 'TSS', 'Log', 'deseq2', 'Differences in ASVs']);
export const stepStatus = writable({});
export const subOperations = readable({
    'Raw data': ['Set Random Seed', 'Preview', 'Quick Explore'],
    // 'Data Perturbation': ['Filter', 'Threshold', 'Transformation', 'R/A Abundance', 'Data Splitting', 'Batch Effect Removal'],
    'Filtering': ['Low Abundance Filtering', 'Prevalence Filtering', 'Variance Filtering', 'No Filtering'],
    'Zero-Handling': ['Pseudocount Addition', 'k-NN Imputation', 'No Zero-Handling'],
    'Normalization': ['TSS', 'CSS', 'TMM', 'CLR', 'No Normalization'],
    'Transformation': ['Log', 'Logit', 'AST', 'No Transformation'],
    'Model Perturbation': ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR', 'All methods calculation', 'View Stability Plot', 'Run Shuffled Analysis', 'ASV Selector', 'json interaction']
});
export const singleSelectOperations = readable({
    // 'Data Perturbation': ['Filter', 'Threshold', 'Transformation', 'R/A Abundance', 'Data Splitting', 'Batch Effect Removal'],
    'Filtering': ['Low Abundance Filtering', 'Prevalence Filtering', 'Variance Filtering', 'No Filtering'],
    'Zero-Handling': ['Pseudocount Addition', 'k-NN Imputation', 'No Zero-Handling'],
    'Normalization': ['TSS', 'CSS', 'TMM', 'CLR', 'No Normalization'],
    'Transformation': ['Log', 'Logit', 'AST', 'No Transformation'],
    'Model Perturbation': ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR']
})
export const selectedOperations = writable({});
export const openMenus = writable({});
export const crossStepMutuallyExclusiveOptions = {
    'Normalization': {
        'CLR': { 'Transformation': ['Log', 'Logit'], 'Zero-Handling': ['No Zero-Handling'] },
        'TMM': { 'Transformation': ['Logit'] },
        'CSS': { 'Transformation': ['Logit'], 'Model Perturbation': ['metagenomeseq'] },
        'No Normalization': { 'Transformation': ['AST'] }
    },
    'Transformation': {
        'Log': { 'Normalization': ['CLR'], 'Zero-Handling': ['No Zero-Handling'] },
        'Logit': { 'Normalization': ['CLR', 'TMM', 'CSS'], 'Zero-Handling': ['No Zero-Handling'] },
        'AST': { 'Normalization': ['No Normalization'] }
    },
    'Zero-Handling': {
        'No Zero-Handling': { 'Transformation': ['Log', 'Logit'], 'Normalization': ['CLR'] }
    },
    'Model Perturbation': {
        'metagenomeseq': { 'Normalization': ['CSS'] }
    }
};
export const selectedColorStep = writable('Filtering');
export const scatterPlotColors = readable({
    'Filtering': ['#FF5733', '#33FF57', '#3357FF', '#FF33F1'],
    'Zero-Handling': ['#FFC300', '#DAF7A6', '#FF5733'],
    'Normalization': ['#C70039', '#900C3F', '#581845', '#FFC300', '#DAF7A6'],
    'Transformation': ['#FF5733', '#C70039', '#900C3F', '#581845'],
    'Model Perturbation': ['#FFC300', '#DAF7A6', '#FF5733', '#C70039', '#900C3F'],
    'Stability Metric': ['#581845', '#FFC300', '#DAF7A6']
});
export const colorStatus = writable({
    'Filtering': ['Low Abundance Filtering', 'Prevalence Filtering', 'Variance Filtering', 'No Filtering'],
    'Zero-Handling': ['Pseudocount Addition', 'k-NN Imputation', 'No Zero-Handling'],
    'Normalization': ['TSS', 'CSS', 'TMM', 'CLR', 'No Normalization'],
    'Transformation': ['Log', 'Logit', 'AST', 'No Transformation'],
    'Model Perturbation': ['deseq2', 'edger', 'maaslin2', 'aldex2', 'metagenomeseq'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR']
})
export const autoLoaded = writable(false);
export const params = writable({});
export const currentHighlightedPath = writable([]);
export const showStartPage = writable(true);