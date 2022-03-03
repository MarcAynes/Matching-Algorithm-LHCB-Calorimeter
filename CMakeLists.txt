cmake_minimum_required(VERSION 3.21)
project(MatchingAlgorithm)

set(CMAKE_CXX_STANDARD 14)

add_executable(Matching-Algorithm main.cpp dataAcquisition.h hit.h RTree.h)
