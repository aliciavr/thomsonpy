#!/bin/bash
# -*- ENCODING: UTF-8 -*-

rm -r dist
rm -r thomsonpy.egg-info
python setup.py sdist
source ~/miniconda3/etc/profile.d/conda.sh
conda activate iaa
pip3 install dist/thomsonpy*
exit