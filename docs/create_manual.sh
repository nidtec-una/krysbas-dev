#!/bin/bash

sphinx-build -M latexpdf source build

cp build/latex/fixme.pdf FIXME.pdf
