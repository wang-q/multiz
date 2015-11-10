cmake_minimum_required(VERSION 3.3)
project(multiz_tba)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

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

add_executable(multiz ${SOURCE_FILES})
