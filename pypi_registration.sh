#!/usr/bin/env bash

git tag 0.1.81
git push --tags
git add .
git commit -m 'Fixed prints to stderr'
git push origin -u mashtree
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
