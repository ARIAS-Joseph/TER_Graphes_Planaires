// ── PALETTE FOR CYCLES ──────────────────────────────────────────────────────
const PALETTE = [
    '#4f7cff', '#ff6b6b', '#51cf66', '#ffd43b',
    '#cc5de8', '#ff922b', '#22d3ee', '#f06595',
    '#74c0fc', '#a9e34b', '#ff8787', '#63e6be'
];

// ── STATE ────────────────────────────────────────────────────────────────────
let graphData= null;
let activeBasis= null;
let activeCycle= null;
let selectedVertex= null;
let selectedEdge= null;
let scale= 1, panX = 0, panY = 0;
let isPanning= false, panStart = null;
let draggingVertex= null, dragMoved = false;

// Modes d'edition
const MODE = { NORMAL: 'normal', ADD_VERTEX: 'add_vertex', ADD_EDGE: 'add_edge' };
let editMode     = MODE.NORMAL;
let addEdgeFirst = null;

let nbIterations = 10; // variable globale Horton

// ── SERVER ───────────────────────────────────────────────────────────────────
const SERVER = 'http://localhost:8080';

function setStatus(msg, type) {
    const el = document.getElementById('statusMsg');
    el.textContent = msg;
    el.className = 'status-msg' + (type ? ' ' + type : '');
}

function setLoading(on) {
    document.getElementById('spinner').style.display = on ? 'block' : 'none';
}

async function callC(action, args) {
    if (!args) args = [];
    if (action !== 'generate_graph' && !graphData) return;
    setLoading(true);
    setStatus('Execution : ' + action + '...');
    try {
        const body = {
            action: action,
            args: args,
            graph: (action !== 'generate_graph') ? buildGraphText() : null
        };
        const res = await fetch(SERVER, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(body)
        });
        if (!res.ok) throw new Error('HTTP ' + res.status);
        const data = await res.json();
        if (data.error) throw new Error(data.error);
        graphData = parseGraphFile(data.graph);
        renderAll();
        setStatus('OK', 'ok');
    } catch (err) {
        setStatus('Erreur : ' + err.message, 'err');
        console.error(err);
    } finally {
        setLoading(false);
    }
}

function runHorton() {
    callC('run_horton', [nbIterations]);
}

function autoRunHorton() {
    if (document.getElementById('autoHorton').checked) {
        runHorton();
    }
}

// ── FILE LOADING ─────────────────────────────────────────────────────────────
document.getElementById('fileInput').addEventListener('change', function(e) {
    const f = e.target.files[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = function(ev) {
        try {
            graphData = parseGraphFile(ev.target.result);
            renderAll();
            setStatus('Fichier charge', 'ok');
        } catch (err) {
            setStatus('Erreur parsing : ' + err.message, 'err');
        }
    };
    reader.readAsText(f);
});

// ── PARSER ───────────────────────────────────────────────────────────────────
function parseGraphFile(text) {
    const lines = text.trim().split('\n').map(function(l) { return l.trim(); }).filter(function(l) { return l.length > 0; });
    let idx = 0;
    const header = lines[idx++].split(/\s+/).map(Number);
    const V = header[0], E = header[1], M = header[2], D = header[3];

    const vertices = [];
    for (let i = 0; i < V; i++) {
        const p = lines[idx++].split(/\s+/).map(Number);
        vertices.push({ id: p[0], x: p[1], y: p[2] });
    }

    const edges = [];
    for (let i = 0; i < E; i++) {
        const p = lines[idx++].split(/\s+/).map(Number);
        edges.push({ id: p[0], u: p[1], v: p[2], label: p[3] !== undefined ? p[3] : p[0] });
    }

    const bases = [];
    for (let b = 0; b < M; b++) {
        const cycles = [];
        for (let c = 0; c < D; c++) {
            const ei = lines[idx++].split(/\s+/).map(Number).filter(function(n) { return !isNaN(n); });
            cycles.push(ei);
        }
        bases.push(cycles);
    }

    const vertexMap = {};
    vertices.forEach(v => { vertexMap[v.id] = v; });

    return { V, E, M, D, vertices, edges, bases, vertexMap };
}

// ── SERIALIZE ────────────────────────────────────────────────────────────────
function buildGraphText() {
    const d = graphData;
    const lines = [];
    lines.push(d.vertices.length + ' ' + d.edges.length + ' ' + d.M + ' ' + d.D);
    d.vertices.forEach(function(v) { lines.push(v.id + ' ' + v.x.toFixed(6) + ' ' + v.y.toFixed(6)); });
    d.edges.forEach(function(e) { lines.push(e.id + ' ' + e.u + ' ' + e.v); }); // 3 champs: id u v
    d.bases.forEach(function(cycles) { cycles.forEach(function(ei) { lines.push(ei.join(' ')); }); });
    return lines.join('\n') + '\n';
}

// ── SAVE ─────────────────────────────────────────────────────────────────────
async function saveGraph() {
    if (!graphData) return;
    const text = buildGraphText();
    if (window.showSaveFilePicker) {
        try {
            const handle = await window.showSaveFilePicker({
                suggestedName: 'graph.txt',
                types: [{ description: 'Fichier texte', accept: { 'text/plain': ['.txt'] } }]
            });
            const w = await handle.createWritable();
            await w.write(text);
            await w.close();
            setStatus('Sauvegarde OK', 'ok');
        } catch (err) {
            if (err.name !== 'AbortError') { setStatus('Erreur sauvegarde', 'err'); }
        }
    } else {
        const blob = new Blob([text], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url; a.download = 'graph.txt'; a.click();
        URL.revokeObjectURL(url);
    }
}
document.getElementById('saveBtn').addEventListener('click', saveGraph);

// ── EDIT MODES ────────────────────────────────────────────────────────────────
document.getElementById('addVertexBtn').addEventListener('click', function() {
    if (!graphData) return;
    setEditMode(editMode === MODE.ADD_VERTEX ? MODE.NORMAL : MODE.ADD_VERTEX);
});

document.getElementById('addEdgeBtn').addEventListener('click', function() {
    if (!graphData) return;
    setEditMode(editMode === MODE.ADD_EDGE ? MODE.NORMAL : MODE.ADD_EDGE);
});

function setEditMode(mode) {
    editMode = mode;
    addEdgeFirst = null;
    const svgEl = document.getElementById('graphSvg');
    svgEl.className = '';
    if (mode === MODE.ADD_VERTEX) svgEl.classList.add('mode-add-vertex');
    if (mode === MODE.ADD_EDGE)   svgEl.classList.add('mode-add-edge');

    document.getElementById('addVertexBtn').classList.toggle('active', mode === MODE.ADD_VERTEX);
    document.getElementById('addVertexHint').style.display = mode === MODE.ADD_VERTEX ? 'block' : 'none';
    document.getElementById('addEdgeBtn').classList.toggle('active', mode === MODE.ADD_EDGE);
    document.getElementById('addEdgeHint').style.display = mode === MODE.ADD_EDGE ? 'block' : 'none';
    if (mode === MODE.ADD_EDGE) {
        document.getElementById('addEdgeHint').textContent = 'Cliquez sur le 1er sommet...';
    }
    updateActionButtons();
    renderGraph();
}

// ── DELETE ────────────────────────────────────────────────────────────────────
document.getElementById('deleteVertexBtn').addEventListener('click', function() {
    if (selectedVertex === null) return;
    const vid = selectedVertex;
    selectedVertex = null;
    updateActionButtons();
    callC('delete_vertex', [vid, nbIterations]);
});

document.getElementById('deleteEdgeBtn').addEventListener('click', function() {
    if (selectedEdge === null) return;
    const eid = selectedEdge;
    selectedEdge = null;
    updateActionButtons();
    callC('delete_edge', [eid, nbIterations]);
});

// ── SPLIT MODAL ───────────────────────────────────────────────────────────────
document.getElementById('splitEdgeBtn').addEventListener('click', function() {
    if (!graphData) return;
    const list = document.getElementById('splitEdgeList');
    list.innerHTML = '';
    graphData.edges.forEach(function(e) {
        const lbl = document.createElement('label');
        lbl.className = 'edge-checkbox';
        const chk = document.createElement('input');
        chk.type = 'checkbox';
        chk.value = e.id;
        if (e.id === selectedEdge) chk.checked = true;
        lbl.appendChild(chk);
        lbl.appendChild(document.createTextNode(' Arete ' + e.id + ' (' + e.u + ' -> ' + e.v + ')'));
        list.appendChild(lbl);
    });
    document.getElementById('splitModal').classList.add('open');
});

document.getElementById('splitCancel').addEventListener('click', function() {
    document.getElementById('splitModal').classList.remove('open');
});

document.getElementById('splitConfirm').addEventListener('click', function() {
    const k = parseInt(document.getElementById('splitK').value) || 1;
    const checked = Array.from(document.querySelectorAll('#splitEdgeList input:checked'));
    if (checked.length === 0) { setStatus('Selectionnez au moins une arete', 'err'); return; }
    const edgeIds = checked.map(function(cb) { return parseInt(cb.value); });
    document.getElementById('splitModal').classList.remove('open');
    callC('split_edges', [k, nbIterations].concat(edgeIds));
});

// ── HORTON BUTTON ─────────────────────────────────────────────────────────────
document.getElementById('runHortonBtn').addEventListener('click', runHorton);

document.getElementById('hortonN').addEventListener('input', function() {
    nbIterations = parseInt(this.value) || 10;
});
document.getElementById('hortonN').addEventListener('change', function() {
    nbIterations = parseInt(this.value) || 10;
});

document.getElementById('generateBtn').addEventListener('click', function() {
    document.getElementById('generateModal').classList.add('open');
});
document.getElementById('generateCancel').addEventListener('click', function() {
    document.getElementById('generateModal').classList.remove('open');
});
document.getElementById('generateConfirm').addEventListener('click', function() {
    var type = document.getElementById('genType').value;
    var nbV  = parseInt(document.getElementById('genVertices').value) || 10;
    var nbE  = parseInt(document.getElementById('genEdges').value) || 20;
    document.getElementById('generateModal').classList.remove('open');
    callC('generate_graph', [type, nbV, nbE, nbIterations]);
});

// ── HELPERS ───────────────────────────────────────────────────────────────────
function updateActionButtons() {
    document.getElementById('deleteVertexBtn').disabled = selectedVertex === null;
    document.getElementById('deleteEdgeBtn').disabled   = selectedEdge === null;
    document.getElementById('splitEdgeBtn').disabled    = !graphData;
    document.getElementById('runHortonBtn').disabled    = !graphData;
}

// ── RECONSTRUCT CYCLE PATH ────────────────────────────────────────────────────
function reconstructCyclePath(edgeIndices, edges) {
    if (edgeIndices.length === 0) return [];
    const cycleEdges = edgeIndices.map(function(eId) { return edges.find(function(e) { return e.id === eId; }); });
    const path = [cycleEdges[0].u, cycleEdges[0].v];
    const used = new Set([0]);
    while (used.size < cycleEdges.length) {
        const current = path[path.length - 1];
        let found = false;
        for (let i = 0; i < cycleEdges.length; i++) {
            if (used.has(i)) continue;
            const edge = cycleEdges[i];
            if (edge.u === current)      { path.push(edge.v); used.add(i); found = true; break; }
            else if (edge.v === current) { path.push(edge.u); used.add(i); found = true; break; }
        }
        if (!found) break;
    }
    return path;
}

// ── RENDER ALL ────────────────────────────────────────────────────────────────
function renderAll() {
    const d = graphData;
    activeBasis = null; activeCycle = null;
    selectedVertex = null; selectedEdge = null;

    document.getElementById('headerMeta').style.display = 'flex';
    document.getElementById('metaV').textContent = d.V;
    document.getElementById('metaE').textContent = d.E;
    document.getElementById('metaM').textContent = d.M;
    document.getElementById('metaD').textContent = d.D;
    document.getElementById('emptyState').style.display = 'none';
    document.getElementById('graphSvg').style.display = 'block';
    document.getElementById('controls').style.display = 'flex';

    updateActionButtons();
    buildSidebar();
    fitGraph();
    renderGraph();
}

// ── SIDEBAR ───────────────────────────────────────────────────────────────────
function buildSidebar() {
    const d = graphData;
    const container = document.getElementById('basesList');
    container.innerHTML = '';
    d.bases.forEach(function(cycles, bIdx) {
        const block = document.createElement('div');
        block.className = 'basis-block';

        const title = document.createElement('div');
        title.className = 'basis-title';
        title.innerHTML = '<span>Basis ' + (bIdx+1) + '</span><span class="basis-badge">' + d.D + ' cycles</span>';
        const cyclesList = document.createElement('div');
        cyclesList.className = 'cycles-list';
        title.addEventListener('click', function() { toggleBasis(bIdx, title, cyclesList); });
        block.appendChild(title);

        cycles.forEach(function(edgeIndices, cIdx) {
            const item = document.createElement('div');
            item.className = 'cycle-item';
            const color = PALETTE[(bIdx * d.D + cIdx) % PALETTE.length];
            const path = reconstructCyclePath(edgeIndices, d.edges);
            const edgesList = edgeIndices.map(function(eId) {
                const e = d.edges.find(function(edge) { return edge.id === eId; });
                return '(' + e.u + ',' + e.v + ')';
            }).join(', ');
            item.innerHTML = '<div class="cycle-swatch" style="background:' + color + '"></div>' +
                '<div class="cycle-info">' +
                '<div class="cycle-name">Cycle ' + (cIdx+1) + '</div>' +
                '<div class="cycle-detail">' + edgeIndices.length + ' edges: ' + edgesList + '</div>' +
                '<div class="cycle-detail">Path: ' + path.join(' -> ') + '</div>' +
                '</div>';
            item.addEventListener('click', function(e) { e.stopPropagation(); selectCycle(bIdx, cIdx, item); });
            cyclesList.appendChild(item);
        });

        block.appendChild(cyclesList);
        container.appendChild(block);
    });
}

function toggleBasis(bIdx, titleEl, cyclesListEl) {
    const isOpen = cyclesListEl.classList.contains('open');
    document.querySelectorAll('.cycles-list').forEach(function(el) { el.classList.remove('open'); });
    document.querySelectorAll('.basis-title').forEach(function(el) { el.classList.remove('active'); });
    document.querySelectorAll('.cycle-item').forEach(function(el) { el.classList.remove('active'); });
    selectedVertex = null; selectedEdge = null;
    if (!isOpen || activeBasis !== bIdx) {
        cyclesListEl.classList.add('open');
        titleEl.classList.add('active');
        activeBasis = bIdx; activeCycle = null;
    } else {
        activeBasis = null; activeCycle = null;
    }
    updateActionButtons();
    renderGraph();
}

function selectCycle(bIdx, cIdx, itemEl) {
    document.querySelectorAll('.cycle-item').forEach(function(el) { el.classList.remove('active'); });
    document.querySelectorAll('.basis-title').forEach(function(el) { el.classList.remove('active'); });
    selectedVertex = null; selectedEdge = null;
    const alreadyActive = activeBasis === bIdx && activeCycle === cIdx;
    if (alreadyActive) {
        activeCycle = null;
    } else {
        activeBasis = bIdx; activeCycle = cIdx;
        itemEl.classList.add('active');
        document.querySelectorAll('.basis-title')[bIdx].classList.add('active');
    }
    updateActionButtons();
    renderGraph();
}

// ── ELEMENT SELECTION ─────────────────────────────────────────────────────────
function selectVertexElement(vid) {
    if (editMode === MODE.ADD_EDGE) {
        if (addEdgeFirst === null) {
            addEdgeFirst = vid;
            document.getElementById('addEdgeHint').textContent = 'Sommet ' + vid + ' selectionne, cliquez sur le 2e...';
            renderGraph();
        } else {
            if (addEdgeFirst === vid) { setStatus('Meme sommet !', 'err'); return; }
            const u = addEdgeFirst;
            addEdgeFirst = null;
            setEditMode(MODE.NORMAL);
            callC('add_edge', [u, vid, nbIterations]);
        }
        return;
    }
    document.querySelectorAll('.cycle-item').forEach(function(el) { el.classList.remove('active'); });
    document.querySelectorAll('.basis-title').forEach(function(el) { el.classList.remove('active'); });
    activeCycle = null;
    selectedVertex = selectedVertex === vid ? null : vid;
    selectedEdge = null;
    updateActionButtons();
    renderGraph();
}

function selectEdgeElement(eid) {
    document.querySelectorAll('.cycle-item').forEach(function(el) { el.classList.remove('active'); });
    document.querySelectorAll('.basis-title').forEach(function(el) { el.classList.remove('active'); });
    activeCycle = null;
    selectedEdge = selectedEdge === eid ? null : eid;
    selectedVertex = null;
    updateActionButtons();
    renderGraph();
}

// ── FIT / PROJECT / UNPROJECT ────────────────────────────────────────────────
function fitGraph() {
    const d = graphData;
    const area = document.getElementById('canvasArea');
    const W = area.clientWidth, H = area.clientHeight, PAD = 60;
    const xs = d.vertices.map(function(v) { return v.x; });
    const ys = d.vertices.map(function(v) { return v.y; });
    const minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs);
    const minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
    const rangeX = maxX - minX || 1, rangeY = maxY - minY || 1;
    scale = Math.min((W - PAD*2) / rangeX, (H - PAD*2) / rangeY);
    panX = (W - rangeX * scale) / 2 - minX * scale;
    panY = (H - rangeY * scale) / 2 - minY * scale;
}

function project(x, y)   { return { px: x * scale + panX, py: y * scale + panY }; }
function unproject(px, py){ return { x: (px - panX) / scale, y: (py - panY) / scale }; }

// ── INTERSECTIONS ─────────────────────────────────────────────────────────────
function segmentsIntersect(p1, p2, p3, p4) {
    const d1x = p2.x-p1.x, d1y = p2.y-p1.y, d2x = p4.x-p3.x, d2y = p4.y-p3.y;
    const cross = d1x*d2y - d1y*d2x;
    if (Math.abs(cross) < 1e-10) return null;
    const dx = p3.x-p1.x, dy = p3.y-p1.y;
    const t = (dx*d2y - dy*d2x) / cross, u = (dx*d1y - dy*d1x) / cross;
    const EPS = 1e-9;
    if (t > EPS && t < 1-EPS && u > EPS && u < 1-EPS)
        return { x: p1.x + t*d1x, y: p1.y + t*d1y };
    return null;
}

function computeIntersections() {
    const d = graphData;
    const crossingEdges = new Set(), points = [];
    for (let i = 0; i < d.edges.length; i++) {
        const ei = d.edges[i];
        const a = d.vertexMap[ei.u], b = d.vertexMap[ei.v];
        if (!a || !b) continue; // sommet supprimé
        for (let j = i+1; j < d.edges.length; j++) {
            const ej = d.edges[j];
            if (ei.u===ej.u||ei.u===ej.v||ei.v===ej.u||ei.v===ej.v) continue;
            const c = d.vertexMap[ej.u], dd = d.vertexMap[ej.v];
            if (!c || !dd) continue; // sommet supprimé
            const pt = segmentsIntersect(a, b, c, dd);
            if (pt) { crossingEdges.add(ei.id); crossingEdges.add(ej.id); points.push(pt); }
        }
    }
    return { crossingEdges: crossingEdges, points: points };
}

// ── RENDER GRAPH ──────────────────────────────────────────────────────────────
const svg = document.getElementById('graphSvg');

function renderGraph() {
    const d = graphData;
    if (!d) return;
    document.getElementById('edgesGroup').innerHTML = '';
    document.getElementById('edgeLabelsGroup').innerHTML = '';
    document.getElementById('nodesGroup').innerHTML = '';
    document.getElementById('intersectionsGroup').innerHTML = '';

    const inter = computeIntersections();
    const crossingEdges = inter.crossingEdges, intersectionPoints = inter.points;
    const nbInter = intersectionPoints.length;

    // Texte en bas
    let infoText = document.getElementById('intersectionInfo');
    if (!infoText) {
        infoText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        infoText.setAttribute('id', 'intersectionInfo');
        infoText.setAttribute('text-anchor', 'middle');
        infoText.setAttribute('dominant-baseline', 'auto');
        infoText.style.fontSize = '13px';
        infoText.style.pointerEvents = 'none';
        svg.appendChild(infoText);
    }
    const svgR = svg.getBoundingClientRect();
    infoText.setAttribute('x', svgR.width / 2);
    infoText.setAttribute('y', svgR.height - 10);
    if (nbInter === 0) {
        infoText.textContent = 'Aucune intersection';
        infoText.style.fill = '#51cf66';
    } else {
        infoText.textContent = nbInter + ' intersection' + (nbInter > 1 ? 's' : '');
        infoText.style.fill = '#ff6b6b';
    }

    // Highlight logic
    const highlightedEdges = new Set(), highlightedNodes = new Set(), cycleColors = {};
    if (activeBasis !== null && activeCycle !== null) {
        const eids = d.bases[activeBasis][activeCycle];
        const color = PALETTE[(activeBasis * d.D + activeCycle) % PALETTE.length];
        eids.forEach(function(eId) {
            highlightedEdges.add(eId); cycleColors[eId] = color;
            const edge = d.edges.find(function(e) { return e.id === eId; });
            highlightedNodes.add(edge.u); highlightedNodes.add(edge.v);
        });
    } else if (selectedVertex !== null && activeBasis !== null) {
        d.bases[activeBasis].forEach(function(edgeIds, cIdx) {
            const path = reconstructCyclePath(edgeIds, d.edges);
            if (path.indexOf(selectedVertex) >= 0) {
                const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
                edgeIds.forEach(function(eId) {
                    highlightedEdges.add(eId); cycleColors[eId] = color;
                    const edge = d.edges.find(function(e) { return e.id === eId; });
                    highlightedNodes.add(edge.u); highlightedNodes.add(edge.v);
                });
            }
        });
    } else if (selectedEdge !== null && activeBasis !== null) {
        d.bases[activeBasis].forEach(function(edgeIds, cIdx) {
            if (edgeIds.indexOf(selectedEdge) >= 0) {
                const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
                edgeIds.forEach(function(eId) {
                    highlightedEdges.add(eId); cycleColors[eId] = color;
                    const edge = d.edges.find(function(e) { return e.id === eId; });
                    highlightedNodes.add(edge.u); highlightedNodes.add(edge.v);
                });
            }
        });
    } else if (activeBasis !== null) {
        d.bases[activeBasis].forEach(function(edgeIds, cIdx) {
            const color = PALETTE[(activeBasis * d.D + cIdx) % PALETTE.length];
            edgeIds.forEach(function(eId) {
                highlightedEdges.add(eId); cycleColors[eId] = color;
                const edge = d.edges.find(function(e) { return e.id === eId; });
                highlightedNodes.add(edge.u); highlightedNodes.add(edge.v);
            });
        });
    }

    const hasHL = highlightedEdges.size > 0;
    const edgesG  = document.getElementById('edgesGroup');
    const intersG = document.getElementById('intersectionsGroup');
    const nodesG  = document.getElementById('nodesGroup');

    // Draw edges
    d.edges.forEach(function(edge) {
        const v1 = d.vertexMap[edge.u];
        const v2 = d.vertexMap[edge.v];
        if (!v1 || !v2) return;
        const p1 = project(v1.x, v1.y), p2 = project(v2.x, v2.y);
        const isHL = highlightedEdges.has(edge.id);
        const isCross = crossingEdges.has(edge.id);
        const isSel = edge.id === selectedEdge;

        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', p1.px); line.setAttribute('y1', p1.py);
        line.setAttribute('x2', p2.px); line.setAttribute('y2', p2.py);
        line.classList.add('edge');
        if (isSel)       { line.classList.add('selected'); }
        else if (isHL)   { line.classList.add('highlighted'); line.style.stroke = cycleColors[edge.id]; }
        else if (isCross && !hasHL) { line.classList.add('crossing'); }
        edgesG.appendChild(line);

        const hit = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        hit.setAttribute('x1', p1.px); hit.setAttribute('y1', p1.py);
        hit.setAttribute('x2', p2.px); hit.setAttribute('y2', p2.py);
        hit.style.stroke = 'transparent'; hit.style.strokeWidth = '12'; hit.style.cursor = 'pointer';
        hit.addEventListener('mouseenter', function(e) { showTooltip(e, 'Arete ' + edge.id + '  (' + edge.u + ' -> ' + edge.v + ')'); });
        hit.addEventListener('mousemove', moveTooltip);
        hit.addEventListener('mouseleave', hideTooltip);
        hit.addEventListener('click', function(e) { e.stopPropagation(); selectEdgeElement(edge.id); });
        edgesG.appendChild(hit);
    });

    // Intersection markers
    intersectionPoints.forEach(function(pt) {
        const p = project(pt.x, pt.y); const R = 5;
        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const l1 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        l1.setAttribute('x1', p.px-R); l1.setAttribute('y1', p.py-R);
        l1.setAttribute('x2', p.px+R); l1.setAttribute('y2', p.py+R);
        l1.classList.add('intersection-marker');
        const l2 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        l2.setAttribute('x1', p.px+R); l2.setAttribute('y1', p.py-R);
        l2.setAttribute('x2', p.px-R); l2.setAttribute('y2', p.py+R);
        l2.classList.add('intersection-marker');
        g.appendChild(l1); g.appendChild(l2);
        intersG.appendChild(g);
    });

    // Draw nodes
    const NODE_R = Math.max(10, Math.min(18, scale * 0.35));
    d.vertices.forEach(function(v) {
        const p = project(v.x, v.y);
        const isHL   = highlightedNodes.has(v.id);
        const isSel  = v.id === selectedVertex;
        const isFirst = v.id === addEdgeFirst;

        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        circle.setAttribute('cx', p.px); circle.setAttribute('cy', p.py);
        circle.setAttribute('r', isFirst ? NODE_R * 1.3 : NODE_R);
        circle.classList.add('node-circle');
        if (isSel)       circle.classList.add('selected');
        else if (isHL)   circle.classList.add('highlighted');
        if (isFirst) circle.style.stroke = '#ffd43b';

        const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        label.setAttribute('x', p.px); label.setAttribute('y', p.py);
        label.textContent = v.id;
        label.classList.add('node-label');
        if (isHL) label.classList.add('highlighted');


        g.appendChild(circle); g.appendChild(label);
        g.style.cursor = 'pointer';
        g.addEventListener('mouseenter', function(e) { showTooltip(e, 'Sommet ' + v.id + '  (' + v.x.toFixed(1) + ', ' + v.y.toFixed(1) + ')'); });
        g.addEventListener('mousemove', moveTooltip);
        g.addEventListener('mouseleave', hideTooltip);
        g.addEventListener('click', function(e) {
            if (dragMoved) { dragMoved = false; return; }
            e.stopPropagation();
            selectVertexElement(v.id);
        });
        g.addEventListener('mousedown', function(e) {
            if (editMode !== MODE.NORMAL) return;
            e.stopPropagation();
            draggingVertex = v; dragMoved = false;
            svg.style.cursor = 'grabbing';
        });
        nodesG.appendChild(g);
    });
}

// ── TOOLTIP ───────────────────────────────────────────────────────────────────
const tooltip = document.getElementById('tooltip');
function showTooltip(e, text) { tooltip.textContent = text; tooltip.style.display = 'block'; moveTooltip(e); }
function moveTooltip(e) {
    const rect = document.getElementById('canvasArea').getBoundingClientRect();
    tooltip.style.left = (e.clientX - rect.left + 12) + 'px';
    tooltip.style.top  = (e.clientY - rect.top  - 28) + 'px';
}
function hideTooltip() { tooltip.style.display = 'none'; }

// ── CANVAS CLICK : add vertex ─────────────────────────────────────────────────
document.getElementById('canvasArea').addEventListener('click', function(e) {
    if (editMode !== MODE.ADD_VERTEX) return;
    const rect = document.getElementById('canvasArea').getBoundingClientRect();
    const pos = unproject(e.clientX - rect.left, e.clientY - rect.top);
    setEditMode(MODE.NORMAL);
    callC('add_vertex', [pos.x.toFixed(3), pos.y.toFixed(3), nbIterations]);
});

// ── ZOOM & PAN ────────────────────────────────────────────────────────────────
svg.addEventListener('click', function(e) {
    if (editMode !== MODE.NORMAL) return;
    if (e.target === svg) {
        selectedVertex = null; selectedEdge = null;
        updateActionButtons(); renderGraph();
    }
});

svg.addEventListener('wheel', function(e) {
    if (!graphData) return;
    e.preventDefault();
    const rect = document.getElementById('canvasArea').getBoundingClientRect();
    const mx = e.clientX - rect.left, my = e.clientY - rect.top;
    const factor = e.deltaY < 0 ? 1.12 : 1/1.12;
    panX = mx - (mx - panX) * factor;
    panY = my - (my - panY) * factor;
    scale *= factor;
    renderGraph();
}, { passive: false });

svg.addEventListener('mousedown', function(e) {
    if (e.button !== 0 || editMode !== MODE.NORMAL) return;
    isPanning = true;
    panStart = { x: e.clientX - panX, y: e.clientY - panY };
    svg.style.cursor = 'grabbing';
});

window.addEventListener('mousemove', function(e) {
    if (draggingVertex !== null) {
        const rect = document.getElementById('canvasArea').getBoundingClientRect();
        const pos = unproject(e.clientX - rect.left, e.clientY - rect.top);
        draggingVertex.x = pos.x; draggingVertex.y = pos.y;
        dragMoved = true; renderGraph(); return;
    }
    if (!isPanning) return;
    panX = e.clientX - panStart.x;
    panY = e.clientY - panStart.y;
    renderGraph();
});

window.addEventListener('mouseup', function() {
    draggingVertex = null; isPanning = false; svg.style.cursor = '';
});

document.getElementById('zoomIn').addEventListener('click', function() {
    if (!graphData) return;
    const a = document.getElementById('canvasArea');
    const cx = a.clientWidth/2, cy = a.clientHeight/2;
    panX = cx-(cx-panX)*1.2; panY = cy-(cy-panY)*1.2; scale *= 1.2; renderGraph();
});
document.getElementById('zoomOut').addEventListener('click', function() {
    if (!graphData) return;
    const a = document.getElementById('canvasArea');
    const cx = a.clientWidth/2, cy = a.clientHeight/2;
    panX = cx-(cx-panX)/1.2; panY = cy-(cy-panY)/1.2; scale /= 1.2; renderGraph();
});
document.getElementById('zoomFit').addEventListener('click', function() {
    if (!graphData) return; fitGraph(); renderGraph();
});
window.addEventListener('resize', function() { if (graphData) { fitGraph(); renderGraph(); } });

// ── KEYBOARD SHORTCUTS ────────────────────────────────────────────────────────
window.addEventListener('keydown', function(e) {
    if (e.key === 'Escape') {
        setEditMode(MODE.NORMAL);
        document.getElementById('splitModal').classList.remove('open');
    }
    if ((e.key === 'Delete' || e.key === 'Backspace') && document.activeElement.tagName !== 'INPUT') {
        if (selectedEdge !== null)   document.getElementById('deleteEdgeBtn').click();
        else if (selectedVertex !== null) document.getElementById('deleteVertexBtn').click();
    }
});