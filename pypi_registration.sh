#!/usr/bin/env bash

git tag 0.1.87
git push --tags
git add .
git commit -m 'runs sistr, meningotype, ngmaster'
git push origin -u master
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
