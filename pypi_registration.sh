#!/usr/bin/env bash

cp README.md README
git tag 0.3.0
git push --tags
git add .
git commit -m '0.3.0 mlst functionality, meningo update'
git push origin -u master
python3 setup.py sdist upload -r pypi
