#!/usr/bin/env bash

cp README.md README
git tag 0.2.7
git push --tags
git add .
git commit -m 'Rolled back ruffus to 2.6.2 bug in 2.6.3'
git push origin -u master
python3 setup.py sdist upload -r pypi
