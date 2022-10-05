#!/bin/bash

path=../$1
echo "The first script parameter is the folder path: src/"$1

# CppCheck
cppcheck --std=c99 --enable=all --force ${path}/*.c
cppcheck --std=c++20 --enable=all --force ${path}/*.cpp
