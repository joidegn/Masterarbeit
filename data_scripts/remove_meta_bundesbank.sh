#!/bin/bash

# removes the metadata (first 5 lines)

for file in "$1"/*
do
    tail -n +6 < $file > $file.new && mv $file.new $file
done
