#!/usr/bin/env bash

git tag 0.1.83
git push --tags
git add .
git commit -m 'Fixing README, install instructions, and minor bugs.'
git push origin -u mashtree
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
