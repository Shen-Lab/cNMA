'''
Created on Apr 2, 2014

@author: oliwa
'''

import argparse
import glob as glob
import sys as sys
from scriptutils import mkdir_p, customCopytree, getFullPathOfURI, makeStringEndWith

def main():
    """ Run the ResultsCombiner and output the absolute path of the created combined results folder. 
    
        Args:
            customArgs: optional args object from an argparser, when this program is not called directly in the shell
    
        Returns:
            Absolute absolute path of the created combined results folder.
    """
    parser = argparse.ArgumentParser(description='Script to copy the content of multiple folders with a common prefix \
    (such as a experimental setup description) into one folder. This script relies on the fact that the string given by --prefix and \
    stored in the variable "afterUniquePrefix" appears after the unique part of the experiment setup, and will not be part of the output name anymore',
    epilog="Returns: Absolute absolute path of the created combined results folder.\n Example: python ResultsCombinter.py /home/oliwa/workspace/TNMA1/output/2ballAtomcomplexk025A6_U1 \n -> all files from all folders /home/oliwa/workspace/TNMA1/output/2ballAtomcomplexk025A6_U1* are copied in the new created folder /home/oliwa/workspace/TNMA1/output/2ballAtomcomplexk025A6_U1")
    parser.add_argument('resultsPath', help='Unique prefix string of the folders whose content should be copied into a new folder')
    parser.add_argument('--prefix', help='String common among folders after unique prefix, which will not be part of the name anymore, default is \"NumberOf\"')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # afterUniquePrefix denotes the string after the unique prefix, this afterUniquePrefix string and all that follows will not be part of the output name
    if args.prefix:
        afterUniquePrefix = args.prefix
    else:
        afterUniquePrefix = "NumberOf"
    
    foldersToCombine = glob.glob(args.resultsPath+"*")
    
    outputfolder = foldersToCombine[0][0:foldersToCombine[0].find(afterUniquePrefix)]
    mkdir_p(outputfolder)
    for resultPath in foldersToCombine:
        customCopytree(resultPath, outputfolder)

    outputfolder = makeStringEndWith(getFullPathOfURI(outputfolder), "/")
    print outputfolder
    return outputfolder
if __name__ == '__main__':
    main()