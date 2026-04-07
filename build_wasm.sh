#!/usr/bin/env bash
# Compile WebAssembly with Emscripten.
# emsdk need to be activated before running this script.
set -e

OUT_DIR="./visualizer"
mkdir -p "$OUT_DIR"

emcc graph.c permutations.c planar_graph_creator.c graph_wasm_api.c \
  -o "$OUT_DIR/graph_visualizer_wasm.js" \
  -s WASM=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME="GraphModule" \
  -s EXPORTED_FUNCTIONS='["_wasm_init","_wasm_generate_planar","_wasm_generate_outer_planar","_wasm_load_from_text","_wasm_load_from_file","_wasm_run_horton","_wasm_add_vertex","_wasm_delete_vertex","_wasm_add_edge","_wasm_delete_edge","_wasm_move_vertex","_wasm_split_edges","_wasm_nb_vertices","_wasm_nb_edges","_wasm_nb_bases","_wasm_basis_dimension","_wasm_face_basis","_wasm_nb_face_basis_outer","_wasm_face_basis_outer_at","_wasm_vertex_deleted","_wasm_vertex_x","_wasm_vertex_y","_wasm_edge_deleted","_wasm_edge_u","_wasm_edge_v","_wasm_edge_label","_wasm_cycle_edge_at","_wasm_write_int","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","FS","HEAPU8","HEAPU32","HEAP32"]' \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s INITIAL_MEMORY=33554432 \
  -s FORCE_FILESYSTEM=1 \
  -sFILESYSTEM=1 \
  -lm -O2

echo "Build WASM OK → $OUT_DIR/graph_visualizer_wasm.{js,wasm}"