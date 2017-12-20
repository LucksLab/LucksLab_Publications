#!/usr/bin/python
#Based on an R script written by Krishna Choudhary, 2017
#Adapted and generalized by Kyle Watters, 2017
#
#Original deriviation in Supplemental Data Analysis of Watters, et al. Nucleic Acids Research, 2017
#"Probing of RNA structures in a positive sense RNA virus reveal selection pressures for structural elements" 
#

from __future__ import division
import sys, getopt, itertools, glob

name = "overlapping_reactivities.py"
help_message = '''

overlapping_reactivities.py takes a directory, or directories, containing reactivity files from a long RNA with overlapping priming sites 
and calculates the reactivity map for the entire RNA. Multiple directories should be used for experimental replicates. 

Usage:
    
overlapping_reactivities.py [options] <reactivities dir> ...
 

General options:
-h, --help                		Opens help message
-v, --version                		Displays version number
-o, --output                            Root name for file output (default: output, 'dinosaur' -> dinosaur_gammas.txt)

'''

def get_version():
    return "0.1.0"

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    def __init__(self):
        pass

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:",
                                       ["version",
                                        "help",
                                        "output"])

        except getopt.error, msg:
            raise Usage(msg)
        
        root_name = 'output'
        for option, value in opts:
            if option in ("-v", "--version"):
                print "%s v%s" % (name,get_version())
                exit(0)
            elif option in ("-h", "--help"):
                raise Usage(help_message)
            elif option in ("-o", "--output"):
                root_name = value
                
        if len(args) < 1:
            raise Usage(help_message)        
            
        return args,root_name

    def check(self):
        pass

def cumsum(list_in):
    total = 0; list_out = []
    for x in list_in:
        total += x
        list_out.append(total)
    return list_out

def import_data(directory):
    
    data = []
    #Look for any file that doesn't start with '.'
    reactivity_files = glob.glob(directory+"*?.*")
    for filename in reactivity_files:
        with open(filename,'rU') as file1:
            lines = file1.readlines()
            data.append([x.strip().split('\t') for x in lines[1:]]) #skip the headers
     
    return data

def coverages_calculate(cumsum_lists):
    
    coverages = []
    for i in range(0,len(max(cumsum_lists,key=len))): #index should be adjusted +1 for true coverages by position, keeping as is for easy math, dropped later at return to _main_
            counts = []
            for x in cumsum_lists:
                try:
                    count = x[i]
                except IndexError:
                    count = 0
                counts.append(count)
            coverages.append(sum(counts))
    return coverages

def calc_gammas_betas(data):
    
#1. First, calculate the number of reads recorded at each position and the local coverages
    positive_reads = []; negative_reads = []; pos_reads_total = []; neg_reads_total = []
    #This block will sum all of the reads by position across a
    for dataset in data:
        pos = [int(position[4]) for position in dataset]
        neg = [int(position[5]) for position in dataset]
        positive_reads.append(pos)
        negative_reads.append(neg)
        pos_reads_total = [sum(x) for x in itertools.izip_longest(pos_reads_total,pos,fillvalue=0)]
        neg_reads_total = [sum(x) for x in itertools.izip_longest(neg_reads_total,neg,fillvalue=0)]
    
    #Get the cumulative sums of each reactivity list (results in a matrix of cumulative read counts by priming site)
    pos_cumsum = [cumsum(dataset) for dataset in positive_reads]
    neg_cumsum = [cumsum(dataset) for dataset in negative_reads]
                   
    #Now get the coverage values for each position. This includes any position upstream of base n, base n 
    positive_coverages = coverages_calculate(pos_cumsum)
    negative_coverages = coverages_calculate(neg_cumsum)
    
#2. Second, calculate the gammas and betas for each position
    gammas = []; betas = []
    for position in range(0,len(pos_reads_total)-1):
        #calculate Xk- (Xkn), Xk+ (Xkp), Ck+ (Ckp), C+k-1 (Ckpm1), etc. using the summed read counts
        Xkp = pos_reads_total[position+1]; Xkn = neg_reads_total[position+1]; Xk = Xkp + Xkn
        Ckp = positive_coverages[position]; Ckn = negative_coverages[position]; Ck = Ckp + Ckn
        Ckpm1 = positive_coverages[position+1]; Cknm1 = negative_coverages[position+1]
        
        #if positive_reads[position+1]/positive_coverages[position+2] >= negative_reads[position+1]/negative_coverages[position+2]:  #Xk+/C+k-1 >= Xk-/C-k-1
        if Xkp/Ckpm1 >= Xkn/Cknm1:
            gamma = Xkn / Cknm1
            beta = ( Xkp / Ckpm1 - gamma ) / ( 1 - gamma )
            betas.append(beta)
        else:
            gamma = Xk / (Xk + Ck)  
            gamma = Xkn / Cknm1
            betas.append(0)
        gammas.append(gamma)
        
    return positive_coverages[1:], negative_coverages[1:], gammas, betas

def output_data(array,root_name='output',datatype=""):
    
    filename = root_name + '_' + datatype + '.txt'
    print_data = zip(*array)
    with open(filename,'w') as file1:
        for data in print_data:
            file1.write("\t".join([str(x) for x in data]) + "\n")
       
def main(argv=None,):

    params = Params()

    try:
        if argv is None:
            argv = sys.argv
            
            args,root_name = params.parse_options(argv)
            params.check()
        
        pos_coverages = []; neg_coverages = []; gamma_sets = []; beta_sets = []
        for arg in args:
            data = import_data(arg)
            positive_coverages, negative_coverages, gammas, betas = calc_gammas_betas(data)
            pos_coverages.append(positive_coverages)
            neg_coverages.append(negative_coverages)
            gamma_sets.append(gammas)
            beta_sets.append(betas)
        
        types_dict = {"pos_coverages":pos_coverages, "neg_coverages":neg_coverages, "gammas":gamma_sets, "betas":beta_sets}
        for key,value in types_dict.iteritems():    
            output_data(value,root_name,key)
        
    
    except Usage, err:
        print sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg) 
        print >> sys.stderr, "" 
        return 2

if __name__ == "__main__":
    sys.exit(main())

