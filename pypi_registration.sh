#!/usr/bin/env bash

cp README.md README
git tag 0.2.6
git push --tags
git add .
git commit -m 'replaced Popen calls with check_output'
git push origin -u master
python3 setup.py sdist upload -r pypi
