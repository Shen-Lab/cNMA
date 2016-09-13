'''
Created on Apr 7, 2014

@author: oliwa
'''

# The methods in scriptutils are useful for IO and other tasks used by other scripts

import os as os
import shutil as shutil
import errno as errno
import numpy as np
import traceback
import subprocess
from itertools import izip_longest
import sys

def customCopytree(src, dst, symlinks=False, ignore=None):
    """ http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth 
    
        Args:
            src: the source to be copied
            dst: the destination to copy to
            symlinks: 3rd argument symlinks for shutil.copytree
            ignore: 4rd argument ignore for shutil.copytree 
    """
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def mkdir_p(path):
    """ Reimplement mkdir -p, 
    from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python 
    
    Args:
        path: the path to the folder to be created
    """
    try:
        os.makedirs(path)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def makeStringEndWith(inputString, symbol):
    """ If the inputString does not end with symbol, append it to inputString and return the result.
        Else, return the inputString. 
        
        Args:
            inputString: the string to be examined to make the output
            symbol: the symbol that the output will end with
            
        Returns:
            inputString that ends with symbol
    """
    if inputString.endswith(symbol):
        return inputString
    else:
        return inputString+symbol
    
def makeStringNotEndWith(inputString, symbol):
    """ If the inputString does end with symbol, remove it from the end of the inputString and return the shortend string.
        Else, return the inputString. 
        
        Args:
            inputString: the string to be examined to make the output
            symbol: the symbol that the output will not end with
            
        Returns:
            inputString that does not end with symbol
    """        
        
    if inputString.endswith(symbol):
        return inputString[:-len(symbol)]
    else:
        return inputString
    
def getFullPathOfURI(uri):
    """ Get the full path of a file, adapted from http://stamat.wordpress.com/2013/06/30/dynamic-module-import-in-python/"""
    return os.path.normpath(os.path.join(os.path.dirname(__file__), uri))

def getParentDirOfURI(uri):
    """ Get the parent directory path of a file, from: http://stackoverflow.com/questions/9856683/using-pythons-os-path-how-do-i-go-up-one-directory"""
    return os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir, uri))

def tryToLoadData(folderPath, dataName, dataFormat=None):
    """ Try to load an results array. """
    try:
        if dataFormat == None:
            return np.loadtxt(folderPath+dataName)
        else:
            return np.loadtxt(folderPath+dataName, dtype=dataFormat)
    except IOError, err:
        print "IOError occurred with tryToLoadData, probably there is no such file at the path: ", err
        print traceback.format_exc()
        ###
        sys.exit()
        ###
        return None
    except Exception, err:
        print "Exception occurred with tryToLoadData when trying to load data, not an IOError: ", err
        print traceback.format_exc()
        ###
        sys.exit()
        ###
        return None
    
def percentageDecrease(x, y):
    """Return the percentage decrease from x towards y
    
        Args:
            x: original value
            y: changed value
        
        Returns: percentage decrease from x towards y.
            
    """
    return (float(x)-float(y))/float(x)

def non_decreasing(L):
    """ Assert non increasing list 
    
        Args:
            L: list
            
        Returns: true if list is non increasing
    """
    threshold = 0.0000001
    result = all(x<= (y+threshold) for x, y in zip(L, L[1:]))
    assert result == True
    return result

def startModeIndexFromOne(arr1):
    """ If the modeindex starts from zero, the values in arr1 correspond to array 
    indices and not mode numbers, so add 1 to each element and return it, else 
    return arr1 unchanged. 
    
    Args:
        arr1: array with stepsize information, to be plotted on the x-axis
    Returns: 
        the array asserted to start with 1
    """
    assert isinstance(arr1, np.ndarray)
    if isinstance(arr1, np.ndarray):
        if arr1[0] == 0:
            return arr1+1
        else:
            return arr1
    elif arr1 is None:
        return arr1
    else:
        raise Exception("array is neither np.ndarray nor None.")
    
def getAllAfterLastChar(string1, lastChar):
    """ Return a substring of string1, where the content of the substring is everything after the last occurence of lastChar in the input string.
    
        Args:
            string1: the input string
            lastChar: the character after which last occurence the subset will be returned
            
        Returns: he substring  containing everything after the last occurence of lastChar in the input string
     """
    return string1[string1.rfind(lastChar)+len(lastChar):]

def runBashCommand(bashCommand):
    """ Run the command in the bash shell and return the result. 
    From: http://stackoverflow.com/questions/4256107/running-bash-commands-in-python
    
    Args:
        command: the full string to execute in the bash shell
        
    Returns: a copy of the output of the command
     """
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    return output[:]

def grouper(n, iterable, fillvalue=None):
    """ from http://docs.python.org/2/library/itertools.html#recipes"""
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)