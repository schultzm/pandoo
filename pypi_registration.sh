#!/usr/bin/env bash

git tag 0.2.0
git push --tags
git add .
git commit -m 'Fixed MLST bug'
git push origin -u master
cp README.md README
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
