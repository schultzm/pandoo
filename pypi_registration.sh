#!/usr/bin/env bash

cp README.md README
git tag 0.2.11
git push --tags
git add .
git commit -m 'Gene classification in place'
git push origin -u master
python3 setup.py sdist upload -r pypi
