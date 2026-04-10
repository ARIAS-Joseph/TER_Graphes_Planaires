# Bases minimales de cycles dans les graphes planaires

Projet TER - Les faces d'un graphe planaire forment-elles toujours une base minimale de cycles ? Etudes des bases de cycles minimales dans les graphes planaires.
---

Ce projet a été réalisé dans le cadre du TER du Master 1 Algorithmique et Modélisation à l'Interface des Sciences (AMIS) de l'Université de Versailles Saint-Quentin-en-Yvelines (UVSQ) en 2026, sous la supervision de Dominique BARTH et Chloé GODET. Le rapport issu de ce TER est présent dans le dossier `Rapport LateX` et est rédigé en français. Il détaille le contexte, les motivations, méthodes, résultats et conclusions de ce projet.

---

## Table des matières

1. [Visualisateur interactif](#visualisateur-interactif)
2. [Lancer le visualisateur](#lancer-le-visualisateur)
3. [Compiler depuis les sources](#compiler-depuis-les-sources)
4. [Étude de masse](#étude-statistique-studyc)
5. [Format de fichier graphe](#format-de-fichier-graphe)

---

## Visualisateur interactif

Le visualisateur est une application web utilisant le code C compilé en **WebAssembly** via [Emscripten](https://emscripten.org/).

### Fonctionnalités

| Fonctionnalité | Description                                                                                            |
|---|--------------------------------------------------------------------------------------------------------|
|Génération de graphes | Graphes planaires ou extérieurement planaires aléatoires (nombre de sommets et d'arêtes paramétrables) |
| Chargement / Sauvegarde | Import/export au format texte                                                                          |
| Visualisation des bases | Toutes les bases minimales trouvées sont listées                                                       |
| Inspection des cycles | Clic sur une base ou un cycle pour le mettre en évidence sur le graphe en couleur                      |
| Édition interactive | Ajout/suppression de sommets et d'arêtes, déplacement par glisser-déposer                              |
| Subdivision d'arêtes| Insertion de $k$ sommets sur une ou plusieurs arêtes                                                   |
| Détection d'intersections | Les croisements d'arêtes sont détectés et affichés                                                     |

---

## Lancer le visualisateur

### Option 1 - Fichiers précompilés (recommandé)

Les fichiers WebAssembly compilés (`graph_visualizer_wasm.js` et `graph_visualizer_wasm.wasm`) sont inclus dans le dossier `visualizer/`. Il suffit de servir ce dossier avec un serveur HTTP local, car les navigateurs bloquent le chargement de fichiers `.wasm` depuis le système de fichiers en local (`file://`).

**Avec Python :**
```bash
cd visualizer
python3 -m http.server 8080
# Ouvrir http://localhost:8080/visualizer.html
```

**Avec Node.js (`npx serve`) :**
```bash
cd visualizer
npx http-server -p 8080
# Ouvrir http://localhost:8080/visualizer.html
```

**Avec VS Code :** installer l'extension *Live Server* et cliquer sur *Go Live* depuis `visualizer.html`.

### Option 2 - Compiler soi-même

Voir la section [Compiler depuis les sources](#compiler-depuis-les-sources) ci-dessous.

---

## Compiler depuis les sources

### Prérequis

- [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html) (`emcc` dans le PATH)
- `make` ou un shell Bash

### Compilation

La commande complète est dans `script_emcc.txt` :

```bash
emcc graph.c permutations.c planar_graph_creator.c graph_wasm_api.c \
  -o ./visualizer/graph_visualizer_wasm.js \
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
```

Cela génère `visualizer/graph_visualizer_wasm.js` et `visualizer/graph_visualizer_wasm.wasm`.

---

## Étude de masse

`study.c` est un programme en ligne de commande (Windows, utilise `_mkdir`) qui génère des milliers de graphes planaires aléatoires et enregistre les résultats dans `result.csv` :

```
Type, nb_vertices, nb_edges, nb_minimal_basis, has_faces, has_faces_outer
```

Actuellement, seuls les graphes dont aucune base minimale trouvée ne peut être l'ensemble des faces sont sauvegardés dans `all_graphs/` pour une inspection manuelle ultérieure.

**Compilation (MSVC ou MinGW) :**
```bash
gcc study.c graph.c permutations.c planar_graph_creator.c -o study -lm
./study
```

Le programme teste des graphes de 10 à 20 sommets, 1000 instances par taille, et sauvegarde les graphes « intéressants » (pour lesquels aucune base trouvée n'est la base des faces) dans `all_graphs/`.

---

## Format de fichier graphe

Si vous souhaitez créer manuellement des graphes à tester via un .txt, le format de fichier texte est le suivant:

```
<nb_sommets> <nb_arêtes> 0 0 -1
<id_sommet> <x> <y>
...
<id_arête> <sommet_u> <sommet_v>
...
```

Exemple (triangle) :
```
3 3 0 0 -1
0 0.0 0.0
1 1.0 0.0
2 0.5 0.5
0 0 1
1 1 2
2 0 2
```

Les graphes enregistrés par le visualisateur suivent ce format, avec les coordonnées des sommets et les étiquettes d'arêtes. Les paramètres 0 0 -1 correspondent respectivement au nombre de bases minimales trouvées, à la dimension des bases et à l'indice d'une base de faces intérieures (ou -1 si inexistante/introuvée). Dans ce cas, après les lignes de sommets et d'arêtes, il y a une section supplémentaire listant les bases minimales trouvées, avec leurs cycles et les arêtes correspondantes.