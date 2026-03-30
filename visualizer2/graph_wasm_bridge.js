const GraphWasm = (() => {

    /** @type {EmscriptenModule|null} */
    let M = null;

    function readGraph() {
        const nbVTotal = M._wasm_nb_vertices();
        const nbETotal = M._wasm_nb_edges();
        const nbBases = M._wasm_nb_bases();
        const D = M._wasm_basis_dimension();
        const faceBasis = M._wasm_face_basis();

        /* Read the full array of outer-face-basis indices */
        const nbFaceBasisOuter = typeof M._wasm_nb_face_basis_outer === 'function'
            ? M._wasm_nb_face_basis_outer() : 0;
        const faceBasisOuter = [];
        if (typeof M._wasm_face_basis_outer_at === 'function') {
            for (let i = 0; i < nbFaceBasisOuter; i++)
                faceBasisOuter.push(M._wasm_face_basis_outer_at(i));
        }


        const vertices  = [];
        const vertexMap = {};
        for (let i = 0; i < nbVTotal; i++) {
            if (M._wasm_vertex_deleted(i)) continue;
            const v = { id: i, x: M._wasm_vertex_x(i), y: M._wasm_vertex_y(i) };
            vertices.push(v);
            vertexMap[i] = v;
        }


        const edges = [];
        for (let i = 0; i < nbETotal; i++) {
            if (M._wasm_edge_deleted(i)) continue;
            edges.push({
                id: i,
                u: M._wasm_edge_u(i),
                v: M._wasm_edge_v(i),
                label: M._wasm_edge_label(i)
            });
        }

        const bases = [];
        for (let b = 0; b < nbBases; b++) {
            const cycles = [];
            for (let c = 0; c < D; c++) {
                const edgeIds = [];
                /*
                 * Use wasm_cycle_edge_at() instead of reading HEAPU32 directly.
                 *
                 * With ALLOW_MEMORY_GROWTH=1, loading a large graph (many bases)
                 * triggers multiple WASM memory growths.  Each growth replaces the
                 * underlying ArrayBuffer, making any previously obtained HEAPU32
                 * reference stale → out-of-bounds reads.
                 * Routing through this C accessor always reads current WASM memory.
                 */
                for (let i = 0; i < nbETotal; i++) {
                    if (M._wasm_cycle_edge_at(b, c, i)) edgeIds.push(i);
                }
                cycles.push(edgeIds);
            }
            bases.push(cycles);
        }

        return {
            V: vertices.length,
            E: edges.length,
            M: nbBases,
            D,
            vertices,
            edges,
            bases,
            vertexMap,
            faceBasis,
            faceBasisOuter        /* now an array of 1-based indices, e.g. [] or [2, 5] */
        };
    }


    return {

        /**
         * Load and initialise the WASM module.  Must be awaited once before
         * any other call.  Mirrors the signature used in visualizer2.js:
         *   await GraphWasm.init();
         */
        async init() {
            M = await GraphModule();    // GraphModule comes from graph_tool.js
            M._wasm_init();
        },

        /**
         * Load a graph from a text string (save_graph format).
         *
         * For large files (many bases), ccall() with a 'string' argument fails
         * because Emscripten copies the string via allocateUTF8OnStack(), which
         * is limited to the WASM stack size (~64 KB).  A file with 1000 bases
         * can be several MB, causing "index out of bounds".
         *
         * Fix: write the text directly to MEMFS with Module.FS.writeFile() and
         * call wasm_load_from_file() with no arguments.
         *
         * @param {string} text
         */
        loadFromText(text) {
            M.FS.writeFile('/tmp/wasm_graph.txt', text);
            M._wasm_load_from_file();
        },

        /**
         * Read the current Graph* and return a plain JS object compatible with
         * parseGraphFile()'s output format.  Call this after any action to
         * refresh graphData in visualizer2.js.
         *
         * @returns {{ V, E, M, D, vertices, edges, bases, vertexMap, faceBasis }}
         */
        readGraph,

        /**
         * Dispatch an action to the WASM module.
         *
         * Accepts the same { action, args } object that the old HTTP bridge
         * used to receive, minus the `graph` text field (which is no longer
         * needed: state lives in WASM memory).
         *
         * Returns a resolved Promise so that the `await` in callC() still works.
         *
         * @param {{ action: string, args: (string|number)[] }} param0
         * @returns {Promise<void>}
         */
        async call({ action, args = [] }) {

            switch (action) {

                /* ── Graph generation ─────────────────────────────── */

                case 'generate_graph': {
                    // args = [type, nb_vertices, nb_edges, n_horton]
                    const [type, nbV, nbE, nH = 10] = args;
                    if (type === 'outer_planar')
                        M._wasm_generate_outer_planar(+nbV, +nbE, +nH);
                    else
                        M._wasm_generate_planar(+nbV, +nbE, +nH);
                    break;
                }

                case 'run_horton': {
                    // args = [n]
                    M._wasm_run_horton(+(args[0] ?? 10));
                    break;
                }

                case 'add_vertex': {
                    // args = [x, y, n_horton]
                    const [x, y, nH = 10] = args;
                    M._wasm_add_vertex(+x, +y, +nH);
                    break;
                }

                case 'delete_vertex': {
                    // args = [vid, n_horton]
                    const [vid, nH = 10] = args;
                    M._wasm_delete_vertex(+vid, +nH);
                    break;
                }

                case 'move_vertex': {
                    // args = [vid, x, y]
                    const [vid, x, y] = args;
                    M._wasm_move_vertex(+vid, +x, +y);
                    break;
                }

                case 'add_edge': {
                    // args = [u, v, n_horton]
                    const [u, v, nH = 10] = args;
                    M._wasm_add_edge(+u, +v, +nH);
                    break;
                }

                case 'delete_edge': {
                    // args = [eid, n_horton]
                    const [eid, nH = 10] = args;
                    M._wasm_delete_edge(+eid, +nH);
                    break;
                }

                case 'split_edges': {
                    const k    = +args[0];
                    const nH   = +(args[1] ?? 10);
                    const eids = args.slice(2).map(Number);
                    const nEids = eids.length;

                    const ptr = M._malloc(nEids * 4);
                    /*
                     * Use wasm_write_int() instead of writing to HEAP32 directly.
                     * HEAP32 can become stale if _malloc itself triggered a memory
                     * growth and replaced the ArrayBuffer.
                     */
                    for (let i = 0; i < nEids; i++)
                        M._wasm_write_int(ptr, i, eids[i]);

                    M._wasm_split_edges(ptr, nEids, k, nH);
                    M._free(ptr);
                    break;
                }

                default:
                    throw new Error('GraphWasm: unknown action "' + action + '"');
            }
            // No return value: caller reads state via GraphWasm.readGraph()
        }
    };

})();