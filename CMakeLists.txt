cmake_minimum_required(VERSION 3.15)
project(TER_Graphes_Planaires C)

set(CMAKE_C_STANDARD 17)

set(COMMON_SOURCES
        graph.c
        permutations.c
        planar_graph_creator.c
        cycle_graph.c
        planar_map_generator.c
        dif_faces_vertices.c
)

add_executable(TER_Graphes_Planaires
        ${COMMON_SOURCES}
        study.c
)
target_link_libraries(TER_Graphes_Planaires m)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test.c")
    add_executable(Test_Graphes_Planaires
            ${COMMON_SOURCES}
            test.c
    )
    target_link_libraries(Test_Graphes_Planaires m)
endif()

add_executable(graph_tool
        ${COMMON_SOURCES}
        study.c
        cycle_graph.c
        cycle_graph.h
        planar_map_generator.c
        planar_generator_test.c
        dif_faces_vertices.c
)
target_link_libraries(graph_tool m)

add_executable(test_dfs
        ${COMMON_SOURCES}
        planar_generator_test.c
)
target_link_libraries(test_dfs m)

add_executable(dif_faces
        ${COMMON_SOURCES}
        dif_faces_vertices.c
)
target_link_libraries(dif_faces m)

set(VISUALISEUR_DIR "${CMAKE_SOURCE_DIR}/visualizer")

add_custom_command(TARGET graph_tool POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:graph_tool>
        "${VISUALISEUR_DIR}/$<TARGET_FILE_NAME:graph_tool>"
        COMMENT "Copie graph_tool → ${VISUALISEUR_DIR}"
)