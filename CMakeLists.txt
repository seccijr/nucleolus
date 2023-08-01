cmake_minimum_required(VERSION 3.22)
project(nucleolus)

set(
        INCLUDE_DIRECTORIES
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\include"
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\concert\\include"
        "C:\\Workspace\\glpk-5.0\\src"
)
set(
        LIB_DIRECTORIES
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\lib\\x64_windows_msvc14\\stat_mda"
        "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\concert\\lib\\x64_windows_msvc14\\stat_mda"
        "C:\\Workspace\\glpk-5.0\\w64"
)

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LIB_DIRECTORIES})

set(
        SECCI_SOURCE_FILES
        source/gen_game.cpp
        source/gen_game.h
        source/common.cpp
        source/common.h
        source/secci.cpp
        source/secci.h
)

add_executable (secci ${SECCI_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(secci INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(secci glpk)

set_target_properties(
  secci
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

set(
        BNF_SOURCE_FILES
        source/gen_game.cpp
        source/gen_game.h
        source/common.cpp
        source/common.h
        source/BNF.cpp
        source/BNF.h
)

add_executable (BNF ${BNF_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(BNF INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(BNF cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  BNF
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)