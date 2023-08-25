#!/bin/bash

pwd

/usr/bin/time -v jupyter nbconvert --execute time_estimation_SF11082.ipynb --to html>info.txt 2>&1