cmake_minimum_required(VERSION 3.22)
project(benedek)

set(
        INCLUDE_DIRECTORIES
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\include"
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\concert\\include"
)
set(
        LIB_DIRECTORIES
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\lib\\x64_windows_msvc14\\stat_mda"
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\concert\\lib\\x64_windows_msvc14\\stat_mda"
)

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LIB_DIRECTORIES})

set(
        SOURCE_FILES
        source/gen_game.cpp
        source/gen_game.h
        source/common.cpp
        source/common.h
        source/BNFTest.cpp
        source/BNF.h
)

add_executable (benedek ${SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(benedek INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(benedek cplex2211.lib concert.lib ilocplex.lib)