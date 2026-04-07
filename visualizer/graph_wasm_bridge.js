const GraphWasm = (() => {

    /** @type {EmscriptenModule|null} */
    let M = null;

    function readGraph() {
        const nbVertices = M._wasm_nb_vertices();
        const nbEdges = M._wasm_nb_edges();
        const nbBases = M._wasm_nb_bases();
        const Dimension = M._wasm_basis_dimension();
        const faceBasis = M._wasm_face_basis();

        const nbFaceBasisOuter = typeof M._wasm_nb_face_basis_outer === 'function'
            ? M._wasm_nb_face_basis_outer() : 0;
        const faceBasisOuter = [];
        if (typeof M._wasm_face_basis_outer_at === 'function') {
            for (let i = 0; i < nbFaceBasisOuter; i++)
                faceBasisOuter.push(M._wasm_face_basis_outer_at(i));
        }

        const vertices  = [];
        const vertexMap = {};
        for (let i = 0; i < nbVertices; i++) {
            if (M._wasm_vertex_deleted(i)) continue;
            const v = { id: i, x: M._wasm_vertex_x(i), y: M._wasm_vertex_y(i) };
            vertices.push(v);
            vertexMap[i] = v;
        }

        const edges = [];
        for (let i = 0; i < nbEdges; i++) {
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
            for (let c = 0; c < Dimension; c++) {
                const edgeIds = [];
                for (let i = 0; i < nbEdges; i++) {
                    if (M._wasm_cycle_edge_at(b, c, i)) edgeIds.push(i);
                }
                cycles.push(edgeIds);
            }
            bases.push(cycles);
        }

        return {
            V: vertices.length,
            E: edges.length,
            B: nbBases,
            D: Dimension,
            vertices,
            edges,
            bases,
            vertexMap,
            faceBasis,
            faceBasisOuter
        };
    }

    return {

        /**
         * Load and initialise the WASM module. Must be awaited once before any other call.
         */
        async init() {
            M = await GraphModule();
            M._wasm_init();
        },

        /**
         * Load a graph from a text string (save_graph format).
         *
         * @param {string} text
         */
        loadFromText(text) {
            M.FS.writeFile('/tmp/wasm_graph.txt', text);
            M._wasm_load_from_file();
        },

        /**
         * Read the current Graph* and return a plain JS object. Call
         * this after any action to refresh graphData in visualizer.js.
         *
         * @returns {{ V, E, B, D, vertices, edges, bases, vertexMap, faceBasis }}
         */
        readGraph,

        /**
         * Dispatch an action to the WASM module. The action is a string that determines which operation to perform,
         * and args is an array of arguments for that operation. The exact meaning of the arguments depends on the
         * action.
         *
         * @param {{ action: string, args: (string|number)[] }} param0
         * @returns {Promise<void>}
         */
        async call({ action, args = [] }) {

            switch (action) {

                case 'generate_graph': {
                    // args = [type, nb_vertices, nb_edges, n_horton]
                    const [type,  nbVertices,  nbEdges,  nbHorton = 1000] = args;
                    if (type === 'outer_planar')
                        M._wasm_generate_outer_planar(+ nbVertices, + nbEdges, + nbHorton);
                    else
                        M._wasm_generate_planar(+ nbVertices, + nbEdges, + nbHorton);
                    break;
                }

                case 'run_horton': {
                    // args = [n]
                    M._wasm_run_horton(+(args[0] ?? 1000));
                    break;
                }

                case 'add_vertex': {
                    // args = [x, y, n_horton]
                    const [x, y,  nbHorton = 1000] = args;
                    M._wasm_add_vertex(+x, +y, + nbHorton);
                    break;
                }

                case 'delete_vertex': {
                    // args = [vertex_id, n_horton]
                    const [vertex_id,  nbHorton = 1000] = args;
                    M._wasm_delete_vertex(+vertex_id, + nbHorton);
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
                    const [u, v,  nbHorton = 1000] = args;
                    M._wasm_add_edge(+u, +v, + nbHorton);
                    break;
                }

                case 'delete_edge': {
                    // args = [eid, n_horton]
                    const [eid,  nbHorton = 1000] = args;
                    M._wasm_delete_edge(+eid, + nbHorton);
                    break;
                }

                case 'split_edges': {
                    const k    = +args[0];
                    const  nbHorton   = +(args[1] ?? 1000);
                    const eids = args.slice(2).map(Number);
                    const nEids = eids.length;

                    const ptr = M._malloc(nEids * 4);
                    for (let i = 0; i < nEids; i++)
                        M._wasm_write_int(ptr, i, eids[i]);

                    M._wasm_split_edges(ptr, nEids, k,  nbHorton);
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