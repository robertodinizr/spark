include(cmake/cpm.cmake)

CPMAddPackage("gh:luihabl/parallel-hashmap#master")

CPMAddPackage(
        NAME hypre
        GITHUB_REPOSITORY "hypre-space/hypre"
        VERSION 2.32.0
        SOURCE_SUBDIR src
)