#!/bin/bash
# -*- ENCODING: UTF-8 -*-

rm -r dist
rm -r thomsonpy.egg-info
python setup.py sdist
pip3 install dist/thomsonpy*
exit