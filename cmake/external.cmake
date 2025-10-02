include(cmake/cpm.cmake)

CPMAddPackage("gh:greg7mdp/parallel-hashmap#master")

CPMAddPackage(
        NAME hypre
        GITHUB_REPOSITORY "hypre-space/hypre"
        VERSION 2.32.0
        SOURCE_SUBDIR src
        OPTIONS
        "HYPRE_WITH_MPI Off"
        "HYPRE_ENABLE_SHARED Off"
)