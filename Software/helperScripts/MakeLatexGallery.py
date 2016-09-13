'''
Created on Apr 17, 2014

@author: oliwa
'''

from scriptutils import makeStringEndWith, getFullPathOfURI, grouper, runBashCommand
import argparse
import sys
import glob
import os
import traceback
import ntpath

def outputLatex(path, experimentDescriptor, shortDescription):
    """ Write the latex document with the plot files as figures 
     
    Args:
        path: path with the pdfs and also the output path
        experimentDescriptor: string describing the experiment
    """
    path = makeStringEndWith(path, "/")
    pictures = sorted(glob.glob(path+"*.pdf"))
     
    # the header
    LaTeXHeader = """
\documentclass[12pt]{article}
\usepackage[latin1]{inputenc}
\usepackage[hmargin=1in,vmargin=1in]{geometry}
\usepackage[pdftex]{graphicx}
\usepackage{subfigure}
\usepackage{amsmath}
\usepackage{textcomp}
\usepackage{float}
\usepackage{txfonts}
\usepackage{graphicx}
\usepackage{fancyhdr}
\usepackage{subfigure}
\usepackage{textpos}
 
\\renewcommand{\headrulewidth}{1pt} 
 
\\begin{document}
 
\\title{Protein-Protein Docking Benchmark 4.0 Experimental Results}
 
\\author{Tomasz Oliwa, oliwa@ttic.edu \\\ TTIC \\\ 
{\small\em \  Date:  \\today }}
\date{}
\maketitle
 
\\begin{abstract}
Benchmark 4.0 plots are presented in this document based on the following plot file:
%s 
\end{abstract}
 
\\fancyhead[R]{\small{Benchmark 4.0 RMSD reductions} }
\\fancyfoot[C]{\\thepage}
 
""" % (experimentDescriptor)
 
# end document
    LaTeXEnd = """\n \end{document} """
 
#one latex page with results
    LaTeXBeginPage = """
\\newpage
 
\\begin{figure}
\\vspace{-4em}
\\begin{textblock}{10.2}[0,0](-0.2,-0.2)\hspace{0.2cm} %s \\ \end{textblock}\n
""" % (shortDescription)
 
    LaTeXEndPage = """
\end{figure}
"""
    pages = ""
    for a,b,c,d,e,f,g,h in grouper(8, pictures):
        singlePage = LaTeXBeginPage
        if a is not None:
            a = ntpath.basename(a)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(a)+"}}\n"
        if b is not None:
            b = ntpath.basename(b)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(b)+"}}\n"
        if c is not None:
            c = ntpath.basename(c)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(c)+"}}\n"                
        if d is not None:
            d = ntpath.basename(d)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(d)+"}}\n"
        if e is not None:
            e = ntpath.basename(e)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(e)+"}}\n"
        if f is not None:
            f = ntpath.basename(f)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(f)+"}}\n"
        if g is not None:
            g = ntpath.basename(g)
            singlePage += "\subfigure{\hspace{-2em}\includegraphics[width=7.7cm]{"+str(g)+"}}\n"
        if h is not None:
            h = ntpath.basename(h)
            singlePage += "\subfigure{\hspace{4.4em}\includegraphics[width=7.7cm]{"+str(h)+"}}\n"
        pages += singlePage
        pages += LaTeXEndPage
    outputLaTeX = LaTeXHeader + str(pages) + LaTeXEnd
    outputName = "latexCombined.tex"
    text_file = open(path+outputName, "w")
    text_file.write(outputLaTeX)
    text_file.close()
     
    # run latex
    try:
        os.chdir(path)
        shellCommand = "pdflatex "+ outputName
        os.system(shellCommand)
        # change the working path back to the directory of the program
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
    except Exception, err:
        print "Exception occurred when combining direct result files: ", err
        print traceback.format_exc()

def prepareExperimentDescriptor(fileWithPathsToAnalysisResults):
    """ Read the input plotfile, and output a subset of it and its first line."""
    # read input plot file
    output = "\n"
    firstLine = None
    with open(fileWithPathsToAnalysisResults, 'r') as inputPlotFile:
        for line in inputPlotFile:
            firstSlashIndex = line.rfind('/')
            secondSlashIndex = line[:firstSlashIndex].rfind('/')
            output = output + line[secondSlashIndex:] + "\n"
            if firstLine == None:
                firstLine = line[secondSlashIndex+1:firstSlashIndex]
    output = output.replace("_", "-")
    firstLine = firstLine.replace("_", "-")
    return output, firstLine

def main():
    """
    usage: MakeLatexGallery.py [-h] filePath
    
    Makes a latex file containing all pdf images in the provided path and runs
    pdflatex on it
    
    positional arguments:
      filePath    Absolute path to the pdf files
    
    optional arguments:
      -h, --help  show this help message and exit
    
    Result: latexCombined.* files in the provided path
    """
    parser = argparse.ArgumentParser(description='Makes a latex file containing all pdf images in the provided path and runs pdflatex on it',
                                     epilog='Result: latexCombined.* files in the provided path')
    parser.add_argument('filePath', help='Absolute path to the pdf files')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    plotFilePath = getFullPathOfURI(glob.glob(args.filePath+"plotFile.txt")[0])
    experimentDescriptor, shortExperimentDescriptor = prepareExperimentDescriptor(plotFilePath)
    outputLatex(args.filePath, experimentDescriptor, shortExperimentDescriptor)

if __name__ == '__main__':
    main()