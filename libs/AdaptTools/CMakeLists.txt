cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

macro(DATE RESULT)
  execute_process(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${RESULT}
                                             OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(TIME RESULT)
  execute_process(COMMAND "date" "+%H:%M:%S" OUTPUT_VARIABLE ${RESULT}
                                             OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

# MSHMET EXECUTABLE
project(mshmet LANGUAGES C)
set(CMAKE_C_FLAGS "-fcommon")

# Add the source files to the mshmet executable
add_library(mshmet STATIC sources/mshmet/boulep.c
                          sources/mshmet/defmet.c
                          sources/commonFiles/eigenv.c
                          sources/mshmet/gradls.c
                          sources/commonFiles/_chrono.c
                          sources/mshmet/hash.c
                          sources/mshmet/hessls.c
                          sources/mshmet/inout.c
                          sources/mshmet/lissag.c
                          sources/mshmet/norpoi.c
                          sources/commonFiles/libmesh5.c
                          sources/commonFiles/memory.c
                          sources/mshmet/mshme1.c
                          sources/mshmet/mshmet.c
                          sources/mshmet/redsim.c
                          sources/mshmet/scalem.c)

# Find the math library
find_library(M_LIB m)

# Set all targets and properties of the mshmet library
target_link_libraries(mshmet PRIVATE ${M_LIB})
target_include_directories(mshmet PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/mshmet>
                                  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/mshmet/include>
                                  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/commonFiles>)
set_property(TARGET mshmet PROPERTY C_STANDARD 99)

# Set public headers
set_target_properties(mshmet PROPERTIES PUBLIC_HEADER $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/mshmet/include>)

## Set the installation properties of the mshmet executable (only the non-already copied files)
install(TARGETS mshmet 
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}")