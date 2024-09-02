// store.js
import { readable, writable } from 'svelte/store';

export const selectedPoints = writable([]);
export const currentPath = writable([]);
export const stepStatus = writable({});
export const subOperations = readable({
    'Raw data': ['Set Random Seed', 'Preview', 'Quick Explore'],
    'Data Perturbation': ['Filter', 'Threshold', 'Transformation', 'R/A Abundance', 'Data Splitting', 'Batch Effect Removal'],
    'Model Perturbation': ['Select Method'],
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR', 'All methods calculation', 'View Stability Plot', 'Run Shuffled Analysis', 'ASV Selector', 'json interaction']
});
export const singleSelectOperations = readable({
    'Stability Metric': ['Differences in ASVs', 'AUROC', 'FDR']
})
export const selectedOperations = writable({});
export const openMenus = writable({});
