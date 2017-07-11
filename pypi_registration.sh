#!/usr/bin/env bash

git tag 0.1.90
git push --tags
git add .
git commit -m 'MLST was updated, added the legacy option'
git push origin -u master
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
