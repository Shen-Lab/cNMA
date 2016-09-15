#!/bin/bash
cp compare.py ~/anaconda2/lib/python2.7/site-packages/prody/proteins/
cp anm.py ~/anaconda2/lib/python2.7/site-packages/prody/dynamics/
cp nmdfile.py ~/anaconda2/lib/python2.7/site-packages/prody/dynamics/

rm -f ~/anaconda2/lib/python2.7/site-packages/prody/proteins/compare.pyc
rm -f ~/anaconda2/lib/python2.7/site-packages/prody/dynamics/anm.pyc
rm -f ~/anaconda2/lib/python2.7/site-packages/prody/dynamics/nmdfile.pyc
