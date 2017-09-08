#!/usr/bin/env bash

cp README.md README
git tag 0.2.10
git push --tags
git add .
git commit -m 'Enhancing gene classification scheme'
git push origin -u master
python3 setup.py sdist upload -r pypi
