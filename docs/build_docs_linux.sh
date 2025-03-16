#!/bin/sh
./clean.sh
make html
./copy_ipynb2example.sh
