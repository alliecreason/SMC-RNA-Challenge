#! /usr/bin/env python

import sys
import os
import math
import getopt
import warnings
import numpy
from collections import Counter 
from collections import defaultdict
from scipy import stats

def usage():
    print """
    quantificationEvaluator -t <truth-tsv>  -i <input-tsv>
    
    Requested Parameters:
        -t/--truth-tsv    [string:    path to truth tsv]
        -i/--input-tsv         [string:    path to input tsv]

    Version:                    1.0.0
          """


#parameters
inFile = ''               #input tsv
inTruthFile = ''          #input truth tsv

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"ht:i:",["help",
                                                 "truth-tsv=",
                                                 "input-tsv="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-i", "--input-tsv"):
            global inFile
            inFile= arg
        elif opt in ("-t","--truth-tsv"):
            global inTruthFile
            inTruthFile = arg



input_values_dic = defaultdict(lambda: 0.0)

def getInputDic():
    infile = "%s" % inFile
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        elif not line.startswith("ENST"):
            continue
        else:
            tmp=line.split("\t")
            name=tmp[0]
            value=float(tmp[1][0:len(tmp[1])-1])
            input_values_dic[name]=value

truth_values = []
input_values = []
expressed_truth_values = []
nonexpressed_truth_values = []
expressed_input_values = []
nonexpressed_input_values = []

def getBothValues():
    infile = "%s" % inTruthFile
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        elif not line.startswith("ENST"):
            continue
        else:
            tmp=line.split("\t")
            name=tmp[0]
            value=float(tmp[1][0:len(tmp[1])-1])
            truth_values.append(value)
            
            if input_values_dic.get(name) is not None:
                input_values.append(input_values_dic[name])
            else:
                input_values.append(0)
            if value >= 1:
                expressed_truth_values.append(value)
                if input_values_dic.get(name) is not None:
                    expressed_input_values.append(input_values_dic[name])
                else:
                    expressed_input_values.append(0)
            else: 
                nonexpressed_truth_values.append(value)
                if input_values_dic.get(name) is not None:
                    nonexpressed_input_values.append(input_values_dic[name])
                else:
                    nonexpressed_input_values.append(0)

def calculateCor():
    cor,p_value=stats.spearmanr(truth_values,input_values)
    exp_spearman,exp_spear_pvalue=stats.spearmanr(expressed_truth_values,expressed_input_values)
    pearson,pearson_pvalue=stats.pearsonr(truth_values,input_values)
    exp_pearson,exp_pear_pvalue=stats.pearsonr(expressed_truth_values,expressed_input_values)
    log_pearson,log_pearson_pvalue=stats.pearsonr(numpy.log(numpy.add(truth_values,0.01)),numpy.log(numpy.add(input_values,0.01)))
    exp_log,exp_log_pvalue=stats.pearsonr(numpy.log(numpy.add(expressed_truth_values,0.01)),numpy.log(numpy.add(expressed_input_values,0.01)))

    total_exp=len(expressed_truth_values)
    total_nonexp=len(nonexpressed_truth_values)

    fp=numpy.count_nonzero(nonexpressed_input_values)
    tn=len(nonexpressed_truth_values) - fp
    spec=float(tn)/(tn+fp)
    
    tp=numpy.count_nonzero(expressed_input_values)
    fn=len(expressed_truth_values) - tp
    sen=float(tp)/(tp+fn)

    final = "isoforms\tspearman\tpearson\tlog_pearson\nAll\t%s\t%s\t%s\nExpressed\t%s\t%s\t%s\n\nSpecificity:\t%s\nSensitivity:\t%s" % (cor,pearson,log_pearson,exp_spearman,exp_pearson,exp_log,spec,sen) 

    print(final)
    return(final)

def main(argv):
    getParameters(argv[1:])
    if inFile=='' or inTruthFile=='':
        usage()
        return 1
    getInputDic()
    getBothValues()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final = calculateCor()
        with open("result.out",'w') as results:
            results.write(final)
            results.close()
    return(0)
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
