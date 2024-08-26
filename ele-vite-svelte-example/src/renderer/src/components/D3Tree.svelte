<script>
    import { onDestroy } from 'svelte';
    import * as d3 from 'd3';
    import { tick } from 'svelte';

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
            const duration = event?.altKey ? 2500 : 250;
            const nodes = root.descendants().reverse();
            const links = root.links();
  
            tree(root);
  
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
                .attr("viewBox", [-marginLeft, left.x - marginTop, width, height]);
  
            const node = gNode.selectAll("g")
                .data(nodes, d => d.id);
  
            const nodeEnter = node.enter().append("g")
                .attr("transform", d => `translate(${source.y0},${source.x0})`)
                .attr("fill-opacity", 0)
                .attr("stroke-opacity", 0)
                .attr("data-name", d => d.data.name)
                .on("click", (event, d) => {
                    d.children = d.children ? null : d._children;
                    update(event, d);
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
  
            // Assign node's x0 and y0 values
            root.eachBefore(d => {
                d.x0 = d.x;
                d.y0 = d.y;
            });
        }
  
        root.x0 = dy / 2;
        root.y0 = 0;
        root.descendants().forEach((d, i) => {
            d.id = i;
            d._children = d.children;
            if (d.depth && d.data.name.length !== 7) d.children = null;
        });
  
        treeRoot = root;
        updateTree = update;
        update(null, root);
    }
  
    export async function highlightPath(path) {
        console.log("Highlighting path in D3tree:", path);
  
        if (!treeRoot || !updateTree) {
            console.error("Tree root or update function not available");
            return;
        }

        // De-highlight the previous path before it's updated
        if (previousPath.length) {
            dehighlightPath(previousPath);
            collapsePath(treeRoot, previousPath);
        }
  
        // Store the new path as the current path
        previousPath = path.slice();
        
        // Expand the new path
        expandPath(treeRoot, path);
        
        // Ensure the DOM is updated before applying highlighting
        updateTree(null, treeRoot);
        await tick(); // Wait for DOM updates to complete

        // Highlight the expanded new path
        highlightPathElements(path);
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

    function collapsePath(node, path) {
        if (path.length === 0) return;

        const childName = path[0];
        const child = findNodeByName(node, childName);

        if (child && child.children) {
            child._children = child.children;
            child.children = null;
            collapsePath(child, path.slice(1));
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
            // Highlight the node in the path
            d3.select(d3TreeContainer).select(`g[data-name="${nodeName}"]`).classed('highlighted', true);
            // Highlight the edge leading to this node
            if (index > 0) {
                d3.select(d3TreeContainer).select(`path[data-target="${nodeName}"]`).classed('highlighted', true);
            }
        });
    }

    function dehighlightPath(path) {
        path.forEach((nodeName) => {
            // Remove highlight from the node
            d3.select(d3TreeContainer).select(`g[data-name="${nodeName}"]`).classed('highlighted', false);
            // Remove highlight from the edge leading to this node
            d3.select(d3TreeContainer).select(`path[data-target="${nodeName}"]`).classed('highlighted', false);
        });
    }
</script>

<style>
  :global(.highlighted circle) {
      stroke: #f00; /* Red color stroke for the highlighted circle */
      stroke-width: 5px;
  }

  :global(.highlighted text) {
      font-weight: bold;
      fill: #f00; /* Red color text for highlighted node */
  }

  :global(.highlighted path) {
      stroke: #f00; /* Red color stroke for the highlighted path */
      stroke-width: 2px;
  }

  .d3-tree-container {
      flex: 0 0 80%; /* Take up 80% of the available space */
      height: 100%; /* Full height of the parent */
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