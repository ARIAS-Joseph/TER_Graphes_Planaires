cmake_minimum_required(VERSION 4.1)
project(TER_Graphes_Planaires C)

set(CMAKE_C_STANDARD 17)

add_executable(TER_Graphes_Planaires graph.c
        permutations.c
        permutations.h
        planar_graph_creator.c
        planar_graph_creator.h
        graph.h
        main.c)
