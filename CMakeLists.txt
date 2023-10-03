cmake_minimum_required(VERSION 3.12)
project(nuc_tmech23)

find_package(pybind11 REQUIRED)

pybind11_add_module(nuc_tmech23 SHARED nuc.cpp)
