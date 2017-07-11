#!/usr/bin/env bash

git tag 0.1.91
git push --tags
git add .
git commit -m 'Fixed the out bug on line 461 in pandoo_tasks.py'
git push origin -u master
python3 setup.py register -r pypitest
python3 setup.py sdist upload -r pypitest
python3 setup.py register -r pypi
python3 setup.py sdist upload -r pypi
