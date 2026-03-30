cmake_minimum_required(VERSION 4.1)
project(TER_Graphes_Planaires C)

set(CMAKE_C_STANDARD 17)

add_executable(TER_Graphes_Planaires
        permutations.c
        permutations.h
        planar_graph_creator.c
        planar_graph_creator.h
        graph.h
        graph.c
        study.c)

target_link_libraries(TER_Graphes_Planaires m)

add_executable(Test_Graphes_Planaires
        permutations.c
        permutations.h
        planar_graph_creator.c
        planar_graph_creator.h
        graph.h
        test.c
        graph.c)

target_link_libraries(Test_Graphes_Planaires m)

add_executable(graph_tool
        planar_graph_creator.c
        planar_graph_creator.h
        permutations.c
        permutations.h
        graph.h
        graph.c
        study.c)

target_link_libraries(graph_tool m)

set(VISUALISEUR_DIR "C:\\Users\\arias\\CLionProjects\\TER_Graphes_Planaires\\visualizer2")

add_custom_command(TARGET graph_tool POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:graph_tool>
        ${VISUALISEUR_DIR}/graph_tool.exe
)