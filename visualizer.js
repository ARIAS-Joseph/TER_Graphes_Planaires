// ── PALETTE FOR CYCLES ──────────────────────────────────────────────────────
const PALETTE = [
    '#4f7cff', '#ff6b6b', '#51cf66', '#ffd43b',
    '#cc5de8', '#ff922b', '#22d3ee', '#f06595',
    '#74c0fc', '#a9e34b', '#ff8787', '#63e6be'
];

// ── STATE ───────────────────────────────────────────────────────────────────
let graphData = null;
let activeBasis = null;
let activeCycle = null;
let selectedVertex = null;  // Vertex ID when clicked
let selectedEdge = null;    // Edge ID when clicked
let scale = 1;
let panX = 0, panY = 0;
let isPanning = false, panStart = null;

// ── FILE LOADING ─────────────────────────────────────────────────────────────
document.getElementById('fileInput').addEventListener('change', e => {
    const f = e.target.files[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = ev => {
        try {
            graphData = parseGraphFile(ev.target.result);
            renderAll();
        } catch(err) {
            alert('Error parsing file:\n' + err.message);
        }
    };
    reader.readAsText(f);
});

// ── PARSER ───────────────────────────────────────────────────────────────────
function parseGraphFile(text) {
    const lines = text.trim().split('\n').map(l => l.trim()).filter(l => l.length > 0);
    let idx = 0;

    // Line 1: V E M D
    const [V, E, M, D] = lines[idx++].split(/\s+/).map(Number);

    // Vertices
    const vertices = [];
    for (let i = 0; i < V; i++) {
        const [id, x, y] = lines[idx++].split(/\s+/).map(Number);
        vertices.push({ id, x, y });
    }

    // Edges
    const edges = [];
    for (let i = 0; i < E; i++) {
        const parts = lines[idx++].split(/\s+/).map(Number);
        edges.push({ id: parts[0], u: parts[1], v: parts[2], weight: parts[3] });
    }

    // Bases (M bases, each with D cycles)
    const bases = [];
    for (let b = 0; b < M; b++) {
        const cycles = [];
        for (let c = 0; c < D; c++) {
            const edgeIndices = lines[idx++].split(/\s+/).map(Number).filter(n => !isNaN(n));
            cycles.push(edgeIndices);
        }
        bases.push(cycles);
    }
    return { V, E, M, D, vertices, edges, bases };
}

// ── RECONSTRUCT CYCLE PATH ──────────────────────────────────────────────────
/**
 * Reconstruit le chemin ordonné d'un cycle à partir de ses arêtes.
 * @param {number[]} edgeIndices - Liste des IDs d'arêtes du cycle
 * @param {Array} edges - Liste de toutes les arêtes du graphe
 * @returns {number[]} - Chemin ordonné des sommets
 */
function reconstructCyclePath(edgeIndices, edges) {
    if (edgeIndices.length === 0) return [];

    // Créer une map des arêtes du cycle
    const cycleEdges = edgeIndices.map(eId => edges.find(e => e.id === eId));

    // Commencer par la première arête
    const path = [cycleEdges[0].u, cycleEdges[0].v];
    const used = new Set([0]);

    // Construire le chemin en trouvant les arêtes adjacentes
    while (used.size < cycleEdges.length) {
        const current = path[path.length - 1];
        let found = false;

        for (let i = 0; i < cycleEdges.length; i++) {
            if (used.has(i)) continue;

            const edge = cycleEdges[i];
            if (edge.u === current) {
                path.push(edge.v);
                used.add(i);
                found = true;
                break;
            } else if (edge.v === current) {
                path.push(edge.u);
                used.add(i);
                found = true;
                break;
            }
        }

        if (!found) break; // Cycle incomplet ou malformé
    }

    return path;
}

// ── RENDER ALL ───────────────────────────────────────────────────────────────
function renderAll() {
    const d = graphData;
    activeBasis = null;
    activeCycle = null;
    selectedVertex = null;
    selectedEdge = null;

    // Update header meta
    document.getElementById('headerMeta').style.display = 'flex';
    document.getElementById('metaV').textContent = d.V;
    document.getElementById('metaE').textContent = d.E;
    document.getElementById('metaM').textContent = d.M;
    document.getElementById('metaD').textContent = d.D;

    // Hide empty state, show SVG
    document.getElementById('emptyState').style.display = 'none';
    document.getElementById('graphSvg').style.display = 'block';
    document.getElementById('controls').style.display = 'flex';

    buildSidebar();
    fitGraph();
    renderGraph();
}

// ── SIDEBAR ──────────────────────────────────────────────────────────────────
function buildSidebar() {
    const d = graphData;
    const container = document.getElementById('basesList');
    container.innerHTML = '';

    d.bases.forEach((cycles, bIdx) => {
        const block = document.createElement('div');
        block.className = 'basis-block';

        const title = document.createElement('div');
        title.className = 'basis-title';
        title.innerHTML = `<span>Basis ${bIdx + 1}</span><span class="basis-badge">${d.D} cycles</span>`;
        title.addEventListener('click', () => toggleBasis(bIdx, title, cyclesList));
        block.appendChild(title);

        const cyclesList = document.createElement('div');
        cyclesList.className = 'cycles-list';

        cycles.forEach((edgeIndices, cIdx) => {
            const item = document.createElement('div');
            item.className = 'cycle-item';
            const color = PALETTE[(bIdx * d.D + cIdx) % PALETTE.length];

            // Reconstruire le chemin du cycle
            const path = reconstructCyclePath(edgeIndices, d.edges);

            // Construire l'affichage
            let cyclePath = '';
            if (path.length > 0) {
                cyclePath = path.join(' → ');
            }

            // Construire la liste des arêtes
            const edgesList = edgeIndices.map(eId => {
                const e = d.edges.find(edge => edge.id === eId);
                return `(${e.u},${e.v})`;
            }).join(', ');

            item.innerHTML = `
        <div class="cycle-swatch" style="background:${color}"></div>
        <div class="cycle-info">
          <div class="cycle-name">Cycle ${cIdx + 1}</div>
          <div class="cycle-detail">${edgeIndices.length} edges: ${edgesList}</div>
          <div class="cycle-detail">Path: ${cyclePath}</div>
        </div>`;
            item.addEventListener('click', (e) => {
                e.stopPropagation();
                selectCycle(bIdx, cIdx, item);
            });
            cyclesList.appendChild(item);
        });

        block.appendChild(cyclesList);
        container.appendChild(block);
    });
}

function toggleBasis(bIdx, titleEl, cyclesListEl) {
    const isOpen = cyclesListEl.classList.contains('open');

    // Close all
    document.querySelectorAll('.cycles-list').forEach(el => el.classList.remove('open'));
    document.querySelectorAll('.basis-title').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.cycle-item').forEach(el => el.classList.remove('active'));

    // Reset element selection
    selectedVertex = null;
    selectedEdge = null;

    if (!isOpen || activeBasis !== bIdx) {
        cyclesListEl.classList.add('open');
        titleEl.classList.add('active');
        activeBasis = bIdx;
        activeCycle = null;
    } else {
        activeBasis = null;
        activeCycle = null;
    }
    renderGraph();
}

function selectCycle(bIdx, cIdx, itemEl) {
    document.querySelectorAll('.cycle-item').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.basis-title').forEach(el => el.classList.remove('active'));

    // Reset element selection
    selectedVertex = null;
    selectedEdge = null;

    const alreadyActive = activeBasis === bIdx && activeCycle === cIdx;
    if (alreadyActive) {
        activeCycle = null;
    } else {
        activeBasis = bIdx;
        activeCycle = cIdx;
        itemEl.classList.add('active');
        // Also highlight the parent basis title
        const titles = document.querySelectorAll('.basis-title');
        titles[bIdx].classList.add('active');
    }
    renderGraph();
}

// ── ELEMENT SELECTION ────────────────────────────────────────────────────────
function selectVertexElement(vertexId) {

    if (activeBasis === null) return;

    // Reset sidebar selection
    document.querySelectorAll('.cycle-item').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.basis-title').forEach(el => el.classList.remove('active'));
    activeCycle = null;

    const alreadyActive = selectedVertex === vertexId;
    if (alreadyActive) {
        selectedVertex = null;
    } else {
        selectedVertex = vertexId;
    }
    renderGraph();
}

function selectEdgeElement(edgeId) {

    if (activeBasis === null) return;

    // Reset sidebar selection
    document.querySelectorAll('.cycle-item').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.basis-title').forEach(el => el.classList.remove('active'));
    activeCycle = null;

    const alreadyActive = selectedEdge === edgeId;
    if (alreadyActive) {
        selectedEdge = null;
    } else {
        selectedEdge = edgeId;
        selectedVertex = null;
    }
    renderGraph();
}

// ── GRAPH RENDERING ──────────────────────────────────────────────────────────
function fitGraph() {
    const d = graphData;
    const area = document.getElementById('canvasArea');
    const W = area.clientWidth;
    const H = area.clientHeight;
    const PAD = 60;

    const xs = d.vertices.map(v => v.x);
    const ys = d.vertices.map(v => v.y);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const rangeX = maxX - minX || 1;
    const rangeY = maxY - minY || 1;

    scale = Math.min((W - PAD*2) / rangeX, (H - PAD*2) / rangeY);
    panX = (W - rangeX * scale) / 2 - minX * scale;
    panY = (H - rangeY * scale) / 2 - minY * scale;
}

function project(x, y) {
    return { px: x * scale + panX, py: y * scale + panY };
}

function renderGraph() {
    const d = graphData;
    const edgesG = document.getElementById('edgesGroup');
    const edgeLabG = document.getElementById('edgeLabelsGroup');
    const nodesG = document.getElementById('nodesGroup');
    edgesG.innerHTML = '';
    edgeLabG.innerHTML = '';
    nodesG.innerHTML = '';

    // Determine which edges & nodes are highlighted
    let highlightedEdges = new Set();
    let highlightedNodes = new Set();
    let cycleColors = {}; // edgeId → color

    if (activeBasis !== null && activeCycle !== null) {
        // MODE 1: Single cycle highlighted (clicked in sidebar)
        const edgeIds = d.bases[activeBasis][activeCycle];
        const color = PALETTE[(activeBasis * d.D + activeCycle) % PALETTE.length];
        edgeIds.forEach(eId => {
            highlightedEdges.add(eId);
            cycleColors[eId] = color;
            const edge = d.edges.find(e => e.id === eId);
            highlightedNodes.add(edge.u);
            highlightedNodes.add(edge.v);
        });
    } else if (selectedVertex !== null && activeBasis !== null) {
        // MODE 2: Vertex selected → highlight all cycles containing this vertex
        d.bases[activeBasis].forEach((edgeIds, cIdx) => {
            const path = reconstructCyclePath(edgeIds, d.edges);
            if (path.includes(selectedVertex)) {
                const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
                edgeIds.forEach(eId => {
                    highlightedEdges.add(eId);
                    cycleColors[eId] = color;
                    const edge = d.edges.find(e => e.id === eId);
                    highlightedNodes.add(edge.u);
                    highlightedNodes.add(edge.v);
                });
            }
        });
    } else if (selectedEdge !== null && activeBasis !== null) {
        // MODE 3: Edge selected → highlight all cycles containing this edge
        d.bases[activeBasis].forEach((edgeIds, cIdx) => {
            if (edgeIds.includes(selectedEdge)) {
                const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
                edgeIds.forEach(eId => {
                    highlightedEdges.add(eId);
                    cycleColors[eId] = color;
                    const edge = d.edges.find(e => e.id === eId);
                    highlightedNodes.add(edge.u);
                    highlightedNodes.add(edge.v);
                });
            }
        });
    } else if (activeBasis !== null) {
        // MODE 4: Whole basis highlighted (clicked basis title)
        d.bases[activeBasis].forEach((edgeIds, cIdx) => {
            const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
            edgeIds.forEach(eId => {
                highlightedEdges.add(eId);
                cycleColors[eId] = color;
                const edge = d.edges.find(e => e.id === eId);
                highlightedNodes.add(edge.u);
                highlightedNodes.add(edge.v);
            });
        });
    }

    const hasHighlight = highlightedEdges.size > 0;

    // Draw edges
    d.edges.forEach(edge => {
        const v1 = d.vertices[edge.u];
        const v2 = d.vertices[edge.v];
        const { px: x1, py: y1 } = project(v1.x, v1.y);
        const { px: x2, py: y2 } = project(v2.x, v2.y);

        const isHL = highlightedEdges.has(edge.id);
        const color = isHL ? cycleColors[edge.id] : null;

        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', x1); line.setAttribute('y1', y1);
        line.setAttribute('x2', x2); line.setAttribute('y2', y2);
        line.classList.add('edge');
        if (isHL) {
            line.classList.add('highlighted');
            line.style.stroke = color;
        } else if (hasHighlight) {
            line.style.opacity = '0.2';
        }
        edgesG.appendChild(line);

        // Hover tooltip
        const hitLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        hitLine.setAttribute('x1', x1); hitLine.setAttribute('y1', y1);
        hitLine.setAttribute('x2', x2); hitLine.setAttribute('y2', y2);
        hitLine.style.stroke = 'transparent';
        hitLine.style.strokeWidth = '12';
        hitLine.style.cursor = 'pointer';
        hitLine.addEventListener('mouseenter', (e) => showTooltip(e, `Edge ${edge.id}  (${edge.u} → ${edge.v})`));
        hitLine.addEventListener('mousemove', moveTooltip);
        hitLine.addEventListener('mouseleave', hideTooltip);
        hitLine.addEventListener('click', (e) => {
            e.stopPropagation();
            selectEdgeElement(edge.id);
        });
        edgesG.appendChild(hitLine);
    });

    // Draw nodes
    const NODE_R = Math.max(10, Math.min(18, scale * 0.35));
    d.vertices.forEach(v => {
        const { px, py } = project(v.x, v.y);
        const isHL = highlightedNodes.has(v.id);

        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');

        const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        circle.setAttribute('cx', px); circle.setAttribute('cy', py);
        circle.setAttribute('r', NODE_R);
        circle.classList.add('node-circle');
        if (isHL) circle.classList.add('highlighted');
        else if (hasHighlight) circle.style.opacity = '0.4';

        const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        label.setAttribute('x', px); label.setAttribute('y', py);
        label.textContent = v.id;
        label.classList.add('node-label');
        if (isHL) label.classList.add('highlighted');
        else if (hasHighlight) label.style.opacity = '0.3';

        g.appendChild(circle);
        g.appendChild(label);
        g.style.cursor = 'pointer';
        g.addEventListener('mouseenter', (e) => showTooltip(e, `Vertex ${v.id}  (${v.x.toFixed(2)}, ${v.y.toFixed(2)})`));
        g.addEventListener('mousemove', moveTooltip);
        g.addEventListener('mouseleave', hideTooltip);
        g.addEventListener('click', (e) => {
            e.stopPropagation();
            selectVertexElement(v.id);
        });
        g.addEventListener('mousedown', (e)=>{
            const rect = canvas.getBoundingClientRect();
            mouseX = e.clientX - rect.left;
            mouseY = e.clientY - rect.top;
            v.x = mouseX;
            v.y = mouseY;
            renderGraph()
        })
        nodesG.appendChild(g);
    });
}

// ── TOOLTIP ──────────────────────────────────────────────────────────────────
const tooltip = document.getElementById('tooltip');
function showTooltip(e, text) {
    tooltip.textContent = text;
    tooltip.style.display = 'block';
    moveTooltip(e);
}
function moveTooltip(e) {
    const area = document.getElementById('canvasArea');
    const rect = area.getBoundingClientRect();
    tooltip.style.left = (e.clientX - rect.left + 12) + 'px';
    tooltip.style.top  = (e.clientY - rect.top  - 28) + 'px';
}
function hideTooltip() { tooltip.style.display = 'none'; }

// ── ZOOM & PAN ────────────────────────────────────────────────────────────────
const svg = document.getElementById('graphSvg');

// Click on background to deselect
svg.addEventListener('click', (e) => {
    if (e.target === svg || e.target.tagName === 'g') {
        selectedVertex = null;
        selectedEdge = null;
        renderGraph();
    }
});

svg.addEventListener('wheel', e => {
    if (!graphData) return;
    e.preventDefault();
    const area = document.getElementById('canvasArea');
    const rect = area.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const factor = e.deltaY < 0 ? 1.12 : 1/1.12;
    panX = mx - (mx - panX) * factor;
    panY = my - (my - panY) * factor;
    scale *= factor;
    renderGraph();
}, { passive: false });

svg.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    isPanning = true;
    panStart = { x: e.clientX - panX, y: e.clientY - panY };
    svg.style.cursor = 'grabbing';
});

window.addEventListener('mousemove', e => {
    if (!isPanning) return;
    panX = e.clientX - panStart.x;
    panY = e.clientY - panStart.y;
    renderGraph();
});

window.addEventListener('mouseup', () => {
    isPanning = false;
    svg.style.cursor = 'default';
});

document.getElementById('zoomIn').addEventListener('click', () => {
    if (!graphData) return;
    const area = document.getElementById('canvasArea');
    const cx = area.clientWidth / 2, cy = area.clientHeight / 2;
    panX = cx - (cx - panX) * 1.2;
    panY = cy - (cy - panY) * 1.2;
    scale *= 1.2;
    renderGraph();
});

document.getElementById('zoomOut').addEventListener('click', () => {
    if (!graphData) return;
    const area = document.getElementById('canvasArea');
    const cx = area.clientWidth / 2, cy = area.clientHeight / 2;
    panX = cx - (cx - panX) / 1.2;
    panY = cy - (cy - panY) / 1.2;
    scale /= 1.2;
    renderGraph();
});

document.getElementById('zoomFit').addEventListener('click', () => {
    if (!graphData) return;
    fitGraph();
    renderGraph();
});


// ── RESIZE ───────────────────────────────────────────────────────────────────
window.addEventListener('resize', () => {
    if (graphData) { fitGraph(); renderGraph(); }
});