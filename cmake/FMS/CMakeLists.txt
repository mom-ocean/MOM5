set(FMS_DIR "${SRC_DIR}/shared")

# FMS static library
file(GLOB_RECURSE FMS_SOURCE LIST_DIRECTORIES false
  ${FMS_DIR}/*.[fF]90 ${FMS_DIR}/*.c)

add_library(FMS ${FMS_SOURCE})

target_include_directories(FMS PRIVATE
  ${FMS_DIR}/include 
  ${FMS_DIR}/mosaic 
  ${FMS_DIR}/drifters 
  ${FMS_DIR}/fms 
  ${FMS_DIR}/mpp/include
)