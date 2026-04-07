cmake_minimum_required(VERSION 3.15)
project(TER_Graphes_Planaires C)

set(CMAKE_C_STANDARD 17)

set(COMMON_SOURCES
        graph.c
        graph.h
        permutations.c
        permutations.h
        planar_graph_creator.c
        planar_graph_creator.h
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
)
target_link_libraries(graph_tool m)

set(VISUALISEUR_DIR "${CMAKE_SOURCE_DIR}/visualizer")

add_custom_command(TARGET graph_tool POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:graph_tool>
        "${VISUALISEUR_DIR}/$<TARGET_FILE_NAME:graph_tool>"
        COMMENT "Copie graph_tool → ${VISUALISEUR_DIR}"
)