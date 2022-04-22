RMDIR dist
RMDIR thomsonpy.egg-info
python setup.py sdist
pip3 install dist/thomsonpy-1.0.0b0.tar.gz