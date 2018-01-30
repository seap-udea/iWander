#!/bin/bash
echo "Installing templates..."
cp templates/* .

echo "Unpacking large files..."
make unpack

echo "Testing compilation..."
make

