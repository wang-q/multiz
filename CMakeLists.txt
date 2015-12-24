cmake_minimum_required(VERSION 3.3)
project(multiz)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(ZLIB REQUIRED)

set(SOURCE_FILES
        multiz.c
        maf.h
        maf.c
        util.h
        util.c
        multi_util.h
        multi_util.c
        mz_preyama.h
        mz_preyama.c
        mz_scores.h
        mz_scores.c
        mz_yama.h
        mz_yama.c
)

include_directories(${ZLIB_INCLUDE_DIRS})

add_executable(multiz ${SOURCE_FILES})
target_link_libraries(multiz ${ZLIB_LIBRARIES})
