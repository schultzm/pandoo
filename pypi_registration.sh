#!/usr/bin/env bash

cp README.md README
git tag 0.2.4
git push --tags
git add .
git commit -m 'Fixed MLST bug'
git push origin -u master
python3 setup.py sdist upload -r pypi
