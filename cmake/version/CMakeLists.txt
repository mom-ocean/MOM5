set(VERSION_DIR "${SRC_DIR}/version")

# common files for both builds
file(GLOB_RECURSE VERSION_SOURCE LIST_DIRECTORIES false FOLLOW_SYMLINKS ${VERSION_DIR}/*.[fF]90 )

add_library(version ${VERSION_SOURCE})