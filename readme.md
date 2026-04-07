# Bases minimales de cycles dans les graphes planaires

> Projet TER - Exploration de la conjecture "*les faces d'un graphe planaire forment-elles toujours une base minimale de cycles ?*" et études des bases minimales de cycles dans les graphes planaires.

---

## Table des matières

1. [Contexte mathématique](#contexte-mathématique)
2. [Structure du projet](#structure-du-projet)
3. [Algorithme de Horton](#algorithme-de-horton)
4. [Visualisateur interactif](#visualisateur-interactif)
5. [Lancer le visualisateur](#lancer-le-visualisateur)
6. [Compiler depuis les sources](#compiler-depuis-les-sources)
7. [Étude statistique](#étude-statistique-studyc)
8. [Format de fichier graphe](#format-de-fichier-graphe)
9. [Architecture du code](#architecture-du-code)

---

## Contexte mathématique

Soit $G = (V, E)$ un graphe connexe non orienté. L'**espace des cycles** de $G$ est un espace vectoriel sur $\mathbb{F}_2$ de dimension $\mu = |E| - |V| + 1$ (le **nombre cyclomatique**). Une **base minimale de cycles** (MCB - *Minimum Cycle Basis*) est une base de cet espace dont la somme des longueurs de cycles est minimale.

Pour un graphe planaire, la décomposition en faces fournit naturellement $\mu$ cycles (les faces intérieures). La question centrale de ce projet est :

> **L'ensemble des faces intérieures d'un graphe planaire constitue-t-il toujours une base minimale de cycles ?**

Pour répondre à cette question, le projet :

1. **Génère** des graphes planaires et extérieurement planaires aléatoires.
2. **Calcule** toutes les bases minimales distinctes trouvées via l'algorithme de Horton avec permutations des étiquettes d'arêtes.
3. **Vérifie** si l'une de ces bases correspond exactement aux faces du graphe (face intérieure complète, ou face extérieure + $\mu-1$ faces intérieures).
4. **Visualise** les résultats de façon interactive.

---

## Structure du projet

```
.
├── graph.h / graph.c               # Structure de données du graphe + algorithme de Horton
├── permutations.h / permutations.c # Tables d'inversion, permutations, Fisher-Yates
├── planar_graph_creator.h / .c     # Génération de graphes planaires/ext-planaires aléatoires
├── graph_wasm_api.c                # API C exposée à WebAssembly via Emscripten
├── study.c                         # Programme CLI pour étude statistique en masse
├── script_emcc.txt                 # Commande de compilation
├── buid_wasm.sh                    # Script de compilation Emscripten
└── visualizer/
    ├── visualizer.html             # Interface principale
    ├── visualizer.css              # Styles
    ├── visualizer.js               # Logique de rendu SVG et interactions
    ├── graph_wasm_bridge.js        # Pont JS ↔ WASM (GraphWasm)
    ├── graph_visualizer_wasm.js    # [généré] Glue JS Emscripten
    └── graph_visualizer_wasm.wasm  # [généré] Module WebAssembly
```

---

## Algorithme de Horton

L'algorithme de Horton (1984) calcule une base minimale de cycles en deux étapes :

1. **Génération des cycles candidats** : pour chaque sommet $v$ et chaque arête $(u, w)$, on construit le cycle formé du plus court chemin $v \to u$, de l'arête $(u, w)$, et du plus court chemin $w \to v$.
2. **Sélection par élimination gaussienne** sur $\mathbb{F}_2$ : on trie les cycles par longueur croissante et on sélectionne les cycles linéairement indépendants jusqu'à obtenir une base de dimension $\mu$.

### Permutations pour trouver plusieurs bases

Pour explorer **plusieurs bases minimales distinctes**, les étiquettes des arêtes sont permutées avant chaque exécution de Horton, ce qui influence l'ordre de sélection à longueur égale. Le projet parcourt des tables d'inversion successives (ou aléatoires) pour accumuler des bases distinctes.

Pour chaque base trouvée, le programme vérifie si elle correspond :
- **exactement aux faces intérieures**;
- ou à **la face extérieure + $\mu-1$ faces intérieures**.

---

## Visualisateur interactif

Le visualisateur est une application web utilisant le code C compilé en **WebAssembly** via [Emscripten](https://emscripten.org/).

### Fonctionnalités

| Fonctionnalité | Description                                                                                               |
|---|-----------------------------------------------------------------------------------------------------------|
| **Génération de graphes** | Graphes planaires ou extérieurement planaires aléatoires (nombre de sommets et d'arêtes paramétrables)    |
| **Chargement / Sauvegarde** | Import/export au format texte natif                                                                       |
| **Visualisation des bases** | Toutes les bases minimales trouvées sont listées                                                          |
| **Inspection des cycles** | Clic sur une base ou un cycle pour le mettre en évidence sur le graphe en couleur                         |
| **Édition interactive** | Ajout/suppression de sommets et d'arêtes, déplacement par glisser-déposer                                 |
| **Subdivision d'arêtes** | Insertion de $k$ sommets sur une ou plusieurs arêtes                                                      |
| **Détection d'intersections** | Les croisements d'arêtes sont détectés et affichés |

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

## Étude statistique

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

Les graphes enregistrés par le visualisateur suivent ce format, avec les coordonnées des sommets et les étiquettes d'arêtes. Les gens 0 0 -1 correspondent alors au nombre de bases minimales trouvées, au nombre de bases trouvées, à la dimension des bases, à la présence d'une base de faces intérieures (son indice), et à la présence d'une ou plusieurs bases composées de la face extérieure et les faces intérieures sauf une. Dans ce cas, après les lignes de sommets et d'arêtes, il y a une section supplémentaire listant les bases minimales trouvées, avec leurs cycles et les arêtes correspondantes.

---

## Architecture du code

```
graph.c
 ├── Gestion du graphe (ajout et suppression de sommets/arêtes, déplacement, subdivision)
 ├── Calcul des plus courts chemins
 ├── Recherche des faces pour le plongement donné du graphe
 ├── Algorithme de Horton (avec ou sans itération sur les permutations d'étiquettes)
 └── Comparaison entre les bases minimales trouvées et les faces

planar_graph_creator.c
 └── Création de graphes planaires et planaires extérieurs aléatoires

permutations.c
 ├── Passer d'une table d'inversion à une permutation
 ├── Générer une table d'inversion aléatoire ou succesive à une autre
 └── Mélange de Fisher-Yates

graph_wasm_api.c - wrappers EMSCRIPTEN_KEEPALIVE exposant l'API C au JS

graph_wasm_bridge.js - GraphWasm

visualizer.js - rendu html, interactions et sauvegarde
```