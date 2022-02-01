PROJECT(ASLI_MANUAL NONE)
cmake_minimum_required(VERSION 2.8)

INCLUDE(UseLATEX.cmake)

# 
SET(ASLI_MANUAL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# 
file(GLOB_RECURSE contents_files RELATIVE ${PROJECT_SOURCE_DIR} *.tex)

# 
ADD_LATEX_DOCUMENT(manual.tex
                   INPUTS ${contents_files}
                   IMAGE_DIRS figures)