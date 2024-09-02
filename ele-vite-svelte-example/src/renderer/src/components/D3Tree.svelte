<script>
    import { onDestroy } from 'svelte';
    import * as d3 from 'd3';

    export let treeData;

    let d3TreeContainer;
    let treeRoot;
    let updateTree;
    
    // Track the previously expanded and highlighted path
    let previousPath = [];

    $: if (treeData && d3TreeContainer) {
        console.log("Tree data updated, rendering tree...");
        d3.select(d3TreeContainer).select("svg").remove(); // Clean up any existing SVG
        renderD3Tree(d3TreeContainer, treeData);
    }

    onDestroy(() => {
        d3.select(d3TreeContainer).select("svg").remove(); // Clean up on component destruction
    });

    function renderD3Tree(container, data) {
        const width = container.offsetWidth || 400;
        const marginTop = 10;
        const marginRight = 10;
        const marginBottom = 10;
        const marginLeft = 40;

        const root = d3.hierarchy(data);
        const dx = 25;
        const dy = (width - marginRight - marginLeft) / (1.5 + root.height);

        const tree = d3.tree().nodeSize([dx, dy]);
        const diagonal = d3.linkHorizontal().x(d => d.y).y(d => d.x);

        const svg = d3.select(container).append("svg")
            .attr("width", width)
            .attr("height", 300)
            .attr("viewBox", [-marginLeft, -marginTop, width, 300])
            .attr("style", "max-width: 100%; height: auto; font: 12px sans-serif; user-select: none;");

            const gLink = svg.append("g")
            .attr("fill", "none")
            .attr("stroke", "#555")
            .attr("stroke-opacity", 0.4)
            .attr("stroke-width", 1.5);

        const gNode = svg.append("g")
            .attr("cursor", "pointer")
            .attr("pointer-events", "all");

        function update(event, source) {
            return new Promise((resolve) => {
                const duration = event?.altKey ? 2500 : 250;
                const nodes = root.descendants().reverse();
                const links = root.links();

                tree(root);

                // Reset all paths and nodes to default state
                d3.selectAll('path')
                    .attr('stroke-opacity', 0.4)
                    .attr('stroke-width', 1.5)
                    .attr('stroke', '#555');

                d3.selectAll('g').selectAll("circle")
                    .attr('stroke', null)
                    .attr('stroke-width', null);

                d3.selectAll('g').selectAll("text")
                    .attr('font-weight', null)
                    .attr('fill', null);

                let left = root;
                let right = root;
                root.eachBefore(node => {
                    if (node.x < left.x) left = node;
                    if (node.x > right.x) right = node;
                });

                const height = right.x - left.x + marginTop + marginBottom + 20;
                const transition = svg.transition()
                    .duration(duration)
                    .attr("height", height)
                    .attr("viewBox", [-marginLeft, left.x - marginTop, width, height])
                    .on('end', resolve);  // Resolve the promise when transition ends

                const node = gNode.selectAll("g")
                    .data(nodes, d => d.id);

                const nodeEnter = node.enter().append("g")
                    .attr("transform", d => `translate(${source.y0},${source.x0})`)
                    .attr("fill-opacity", 0)
                    .attr("stroke-opacity", 0)
                    .attr("data-name", d => d.data.name)
                    .on("click", (event, d) => {
                        console.log("Clicked node:", d);
                    });

                nodeEnter.append("circle")
                    .attr("r", 2.5)
                    .attr("fill", d => d._children ? "#555" : "#999")
                    .attr("stroke-width", 10);

                nodeEnter.append("text")
                    .attr("dy", "1em")
                    .attr("x", d => d._children ? -6 : 6)
                    .attr("text-anchor", d => d._children ? "end" : "start")
                    .text(d => d.data.name)
                    .clone(true).lower()
                    .attr("stroke-linejoin", "round")
                    .attr("stroke-width", 3)
                    .attr("stroke", "white");

                const nodeUpdate = node.merge(nodeEnter).transition(transition)
                    .attr("transform", d => `translate(${d.y},${d.x})`)
                    .attr("fill-opacity", 1)
                    .attr("stroke-opacity", 1);

                const nodeExit = node.exit().transition(transition).remove()
                    .attr("transform", (d) => `translate(${source.y},${source.x})`)
                    .attr("fill-opacity", 0)
                    .attr("stroke-opacity", 0);

                const link = gLink.selectAll("path")
                    .data(links, d => d.target.id);

                const linkEnter = link.enter().append("path")
                    .attr("d", (d) => {
                        const o = { x: source.x0, y: source.y0 };
                        return diagonal({ source: o, target: o });
                    })
                    .attr("data-source", d => d.source.data.name)
                    .attr("data-target", d => d.target.data.name);

                link.merge(linkEnter).transition(transition)
                    .attr("d", diagonal);

                link.exit().transition(transition).remove()
                    .attr("d", (d) => {
                        const o = { x: source.x, y: source.y };
                        return diagonal({ source: o, target: o });
                    });

                root.eachBefore(d => {
                    d.x0 = d.x;
                    d.y0 = d.y;
                });
            });
        }

        root.x0 = dy / 2;
        root.y0 = 0;
        root.descendants().forEach((d, i) => {
            d.id = i;
            d._children = d.children;
            if (d.depth > 0) d.children = null;
        });

        treeRoot = root;
        updateTree = update;
        update(null, root);
    }

    export async function highlightPath(path) {
        console.log("=====Begin of *highlightPath()*=====");
        console.log("Highlighting path in D3tree:", path);

        if (!treeRoot || !updateTree) {
            console.error("Tree root or update function not available");
            return;
        }

        if (previousPath.length) {
            console.log("Dehighlighting previous path:", previousPath);
            dehighlightPathElements(previousPath);
            collapsePath(treeRoot);
        }

        previousPath = path.slice();
        expandPath(treeRoot, path);
        await updateTree(null, treeRoot);

        highlightPathElements(path);

        console.log("=====End of *highlightPath()*=====");
    }

    function expandPath(node, path) {
        if (path.length === 0) return;

        const childName = path[0];
        const child = findNodeByName(node, childName);

        if (child) {
            if (child._children) {
                child.children = child._children;
                child._children = null;
            }
            expandPath(child, path.slice(1));
        }
    }

    function collapsePath(node) {
        if (!node) return;

        if (node.children) {
            node._children = node.children;
            node.children = null;
        }

        if (node._children) {
            node._children.forEach(collapsePath);
        }
    }

    function findNodeByName(node, name) {
        if (node.data.name === name) return node;

        if (node.children) {
            for (let child of node.children) {
                const found = findNodeByName(child, name);
                if (found) return found;
            }
        }

        if (node._children) {
            for (let child of node._children) {
                const found = findNodeByName(child, name);
                if (found) return found;
            }
        }

        return null;
    }

    function highlightPathElements(path) {
        path.forEach((nodeName, index) => {
            d3.select(d3TreeContainer).select(`g[data-name="${nodeName}"]`).classed('highlighted', true);
            
            if (index > 0) {
                const previousNodeName = path[index - 1];
                const pathSelection = d3.select(d3TreeContainer).select(`path[data-source="${previousNodeName}"][data-target="${nodeName}"]`);

                pathSelection.classed('highlighted', true)
                    .transition()
                    .duration(250)
                    .attr('stroke-opacity', 1.0)
                    .attr('stroke-width', 2)
                    .attr('stroke', '#f00');
            }
        });
    }

    function dehighlightPathElements(path) {
        path.forEach((nodeName, index) => {
            d3.select(d3TreeContainer).select(`g[data-name="${nodeName}"]`).classed('highlighted', false);

            if (index > 0) {
                const previousNodeName = path[index - 1];
                const pathSelection = d3.select(d3TreeContainer).select(`path[data-source="${previousNodeName}"][data-target="${nodeName}"]`);

                pathSelection.classed('highlighted', false)
                    .transition()
                    .duration(250)
                    .attr('stroke-opacity', 0.4)
                    .attr('stroke-width', 1.5)
                    .attr('stroke', '#555');
            }
        });

        d3.select(d3TreeContainer).selectAll('path:not(.highlighted)')
            .transition()
            .duration(250)
            .attr("fill", "none")
            .attr('stroke-opacity', 0.4)
            .attr('stroke-width', 1.5)
            .attr('stroke', '#555');
    }
</script>

<style>
    :global(.highlighted circle) {
        stroke: #f00;
        stroke-width: 4px;
        transition: stroke 0.25s, stroke-width 0.25s;
    }

    :global(.highlighted text) {
        font-weight: bold;
        fill: #f00;
        transition: fill 0.25s;
    }

    .d3-tree-container {
        flex: 0 0 80%;
        height: 100%;
        border: 1px solid #ddd;
        background-color: #fff;
        overflow: auto; /* Allow scrolling if the tree content overflows */
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        padding: 10px;
        border-radius: 8px;
        box-sizing: border-box;
    }
</style>

<div class="d3-tree-container" bind:this={d3TreeContainer}></div>