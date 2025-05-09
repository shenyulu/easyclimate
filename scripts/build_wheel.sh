#!/bin/sh
python -m build --wheel --no-isolation --outdir dist/
rm -rf ./build
