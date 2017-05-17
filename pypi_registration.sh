#!/usr/bin/env bash

git tag 0.1.82
git push --tags
git add .
git commit -m 'Allowed manual filter of abricate dbs, added qacGenes db'
git push origin -u mashtree
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
