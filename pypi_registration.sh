#!/usr/bin/env bash

git tag 0.1.88
git push --tags
git add .
git commit -m 'sistr,meningo,ngmaster, only run if sp detected, and accepts input table from stdin'
git push origin -u master
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
