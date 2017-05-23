#!/usr/bin/env bash

git tag 0.1.89
git push --tags
git add .
git commit -m 'working to modularise the pipeline'
git push origin -u master
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
