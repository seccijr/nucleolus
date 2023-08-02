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
        COMMON_SOURCE_FILES
        source/gen_game.cpp
        source/gen_game.h
        source/common.cpp
        source/common.h
)

#####################
#        BNF        #
#####################

set(
        BNF_SOURCE_FILES
        source/BNF.cpp
        source/BNF.h
)

add_executable (BNF ${COMMON_SOURCE_FILES} ${BNF_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(BNF INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(BNF cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  BNF
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#        DK         #
#####################

set(
        DK_SOURCE_FILES
        source/DK.cpp
        source/DK.h
)

add_executable (DK ${COMMON_SOURCE_FILES} ${DK_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(DK INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(DK cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  DK
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#        DP         #
#####################

set(
        DP_SOURCE_FILES
        source/DP.cpp
        source/DP.h
)

add_executable (DP ${COMMON_SOURCE_FILES} ${DP_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(DP INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(DP cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  DP
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#        PRA        #
#####################

set(
        PRA_SOURCE_FILES
        source/PRA.cpp
        source/PRA.h
)

add_executable (PRA ${COMMON_SOURCE_FILES} ${PRA_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(PRA INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(PRA cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  PRA
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#        SD         #
#####################

set(
        SD_SOURCE_FILES
        source/SD.cpp
        source/SD.h
)

add_executable (SD ${COMMON_SOURCE_FILES} ${SD_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(SD INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(SD cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  SD
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#        SP         #
#####################

set(
        SP_SOURCE_FILES
        source/SP.cpp
        source/SP.h
)

add_executable (SP ${COMMON_SOURCE_FILES} ${SP_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(SP INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(SP cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  SP
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#       GLPK        #
#####################

set(
        GLPK_SOURCE_FILES
        source/GLPK.cpp
        source/GLPK.h
)

add_executable (GLPK ${COMMON_SOURCE_FILES} ${GLPK_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(GLPK INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(GLPK glpk)

set_target_properties(
  GLPK
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)

#####################
#    OUTPUT_GAME    #
#####################

set(
        OUTPUT_GAME_SOURCE_FILES
        source/output_game.cpp
)

add_executable (OUTPUT_GAME ${COMMON_SOURCE_FILES} ${OUTPUT_GAME_SOURCE_FILES})

add_compile_definitions(WIN32 CONSOLE IL_STD)

target_include_directories(OUTPUT_GAME INTERFACE ${INCLUDE_DIRECTORIES})

target_link_libraries(OUTPUT_GAME cplex2211.lib concert.lib ilocplex.lib)

set_target_properties(
  OUTPUT_GAME
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/executable"
)