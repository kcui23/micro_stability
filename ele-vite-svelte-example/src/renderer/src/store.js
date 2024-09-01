// store.js
import { writable } from 'svelte/store';

export const selectedPoints = writable([]);
export const currentPath = writable([]);
export const stepStatus = writable({});
export const subOperations = writable({});
export const selectedOperations = writable({});
export const openMenus = writable({});
