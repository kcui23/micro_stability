// store.js
import { readable, writable } from 'svelte/store';

export const selectedPoints = writable([]);
export const currentPath = writable(['Raw data', 'Data Perturbation', 'Filter', 'Model Perturbation', 'deseq2', 'Stability Metric', 'Differences in ASVs']);
export const stepStatus = writable({});
export const subOperations = readable({
    'Raw data': ['Set Random Seed', 'Preview', 'Quick Explore'],
    // 'Data Perturbation': ['Filter', 'Threshold', 'Transformation', 'R/A Abundance', 'Data Splitting', 'Batch Effect Removal'],
    'Filtering': ['Low Abundance Filtering', 'Prevalence Filtering', 'Variance Filtering', 'None'],
    'Zero-Handling': ['Pseudocount Addition', 'k-NN Imputation', 'None'],
    'Normalization': ['TSS', 'CSS', 'TMM', 'CLR', 'None'],
    'Transformation': ['Log', 'Logit', 'AST', 'None'],
    'Model Perturbation': ['Select Method'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR', 'All methods calculation', 'View Stability Plot', 'Run Shuffled Analysis', 'ASV Selector', 'json interaction']
});
export const singleSelectOperations = readable({
    // 'Data Perturbation': ['Filter', 'Threshold', 'Transformation', 'R/A Abundance', 'Data Splitting', 'Batch Effect Removal'],
    'Filtering': ['Low Abundance Filtering', 'Prevalence Filtering', 'Variance Filtering', 'None'],
    'Zero-Handling': ['Pseudocount Addition', 'k-NN Imputation', 'None'],
    'Normalization': ['TSS', 'CSS', 'TMM', 'CLR', 'None'],
    'Transformation': ['Log', 'Logit', 'AST', 'None'],
    'Model Perturbation': ['Select Method'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR']
})
export const selectedOperations = writable({});
export const openMenus = writable({});