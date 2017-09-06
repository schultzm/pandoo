#!/usr/bin/env bash

cp README.md README
git tag 0.2.9
git push --tags
git add .
git commit -m 'Changed alerts header to move it to first col'
git push origin -u master
python3 setup.py sdist upload -r pypi
