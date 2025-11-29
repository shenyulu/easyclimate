#!/bin/sh
./clean_all.sh
make html
./copy_ipynb2example.sh
