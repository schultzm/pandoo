#!/usr/bin/env bash

cp README.md README
git tag 0.3.1
git push --tags
git add .
git commit -m 'Fixed tag and version'
git push origin -u master
python3 setup.py sdist upload -r pypi
