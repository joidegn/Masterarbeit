#!/bin/bash

# replaces comma with space for all files in the directory found in the first argument

for file in "$1"/*
do
    sed 's/ //g' < $file > $file.new && mv $file.new $file  # first remove whitespace
    sed 's/,/ /g' < $file > $file.new && mv $file.new $file # then replace comma
done
