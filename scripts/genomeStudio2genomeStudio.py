'''
genomeStudio2genomeStudio.py
=========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Manipulate genomeStudio output.

Usage
-----

Example::

   python genomeStudio2genomeStudio.py

Type::

   python genomeStudio2genomeStudio.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections
import CGAT.Experiment as E

####################################################
####################################################
####################################################
def readSampleSheet(sampleSheet, chop = False):
    '''
    return dictionary mapping labid to
    sample name - has to be a tab separated
    file of type labid\tsamplecondition (no header)
    '''
    inf = open(sampleSheet)
    mapping = {}
    for line in inf.readlines():
        data = line[:-1].split("\t")
        labid, description = data[0], data[1]

        # no spaces in description please
        description = description.replace(" ", "")
        mapping[labid] = description

    reps = collections.defaultdict(int)
    for labid, description in mapping.iteritems():
        reps[description] += 1        
        mapping[labid] = description + ".R%i" % reps[description]
        
    return mapping
        

####################################################
####################################################
####################################################
def outputMapping(probe_profile, mapping, type = "probe"):
    '''
    output renamed sample probe profile file
    '''
    inf = probe_profile
    outf = open("originalid2newid.tsv", "w")
    for line in inf.readlines():
        if line.startswith("TargetID"):
            header = line.strip().split("\t")
            break
    found = set()
    for h in header[2:]:
        k = h.split("-")[1]
        if k in found:
            continue
        else:
            found.add(k)
            outf.write("\t".join([k, mapping[k]]) + "\n")
    outf.close()

####################################################
####################################################
####################################################
def renameColumns(probe_profile, mapping, type = "probe"):
    '''
    output renamed sample probe profile file
    '''
    inf = probe_profile
    for line in inf.readlines():
        if line.startswith("TargetID"):
            header = line.strip().split("\t")
            if type == "probe":
                new_header = ['TargetID', 'ProbeID'] +  [x.split("-")[0] + "-" + mapping[x.split("-")[1]] for x in header[2:]]
            else:
                new_header = ['TargetID', 'ProbeID'] +  [mapping[x.split(".")[0]] + "." + x.split(".")[1] for x in header[2:]]
            sys.stdout.write("\t".join(new_header) + "\n")
        else:
            data = line.strip().split("\t")
            sys.stdout.write("\t".join(data) + "\n")

####################################################
####################################################
####################################################
def remove_columns(probe_profile, sample_list):
    '''
    remove columns according to sample list
    '''
    indices_to_remove = []
    inf = probe_profile
    for line in inf.readlines():
        if line.startswith("TargetID"):
            header = line.strip().split("\t")
            for sample in sample_list:
                to_remove = [i for i in range(len(header)) if re.match(".*%s$" % sample, header[i])]
                indices_to_remove.extend(to_remove)
            new_header = [header[i] for i in range(len(header)) if i not in indices_to_remove]
            sys.stdout.write("\t".join(new_header) + "\n")
        else:
            if indices_to_remove:
                data = line.strip().split("\t")
                data = [data[i] for i in range(len(data)) if i not in indices_to_remove]
                sys.stdout.write("\t".join(data) + "\n")
            else:
                data = line.strip().split("\t")
                sys.stdout.write("\t".join(data) + "\n")

####################################################
####################################################
####################################################
def reorder_columns(probe_profile, sample_list):
    '''
    place the "sample_list" at the end of the file
    used as a control in downstream analysis
    '''
    indices_for_end = []
    inf = probe_profile
    for line in inf.readlines():
        if line.startswith("TargetID"):
            header = line.strip().split("\t")
            for sample in sample_list:
                for_end = [i for i in range(len(header)) if header[i].find("-"+sample) != -1]
                indices_for_end.extend(for_end)
            new_header = [header[i] for i in range(len(header)) if i not in indices_for_end] + [header[i] for i in range(len(header)) if i in indices_for_end]
            sys.stdout.write("\t".join(new_header) + "\n")
        else:
            if indices_for_end:
                data = line.strip().split("\t")
                data = [data[i] for i in range(len(data)) if i not in indices_for_end] + [data[i] for i in range(len(data)) if i in indices_for_end]
                sys.stdout.write("\t".join(data) + "\n")
            else:
                data = line.strip().split("\t")
                sys.stdout.write("\t".join(data) + "\n")

####################################################
####################################################
####################################################
def combine_data(file1, file2):
    '''
    combine two genome studio files 
    '''
    inf1, inf2 = file1, open(file2)
    probe2data = {}
    probe2gene = {}
    
    # store info for file 1
    for line in inf1.readlines():
        if line.startswith("TargetID"):
            header = line.strip().split("\t")
        else:
            data = line.strip().split("\t")
            if len(data) == 1:
                continue
            else:
                gene, probe = data[0], data[1]
                data = data[2:]
                probe2gene[probe] = gene
                probe2data[probe] = data
            

    # itreate over file 2 and output merged rows
    for line in inf2.readlines():
        if line.startswith("TargetID"):
            # extend previous header
            header += line.strip().split("\t")[2:]
            sys.stdout.write("\t".join(header) + "\n")
        else:
            data = line.strip().split("\t")
            if len(data) == 1:
                sys.stdout.write(data[0] + "\n")
            else:
                gene, probe = data[0], data[1]
                data = data[2:]
                sys.stdout.write("\t".join(
                    [gene, probe] + probe2data[probe] + data) + "\n")

####################################################
####################################################
####################################################
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-d", "--do", dest="do", type="choice",
                      choices = ("rename-columns", "reorder-by", "remove", "combine-files"),
                      help="what do you want to do?"  )
    parser.add_option("-l", "--sample-list", dest="sample_list", type="string",
                      help="list of samples to remove"  )
    parser.add_option("-m", "--sample-sheet", dest="sample_sheet", type="string",
                      help="rename columns according to sample sheet")
    parser.add_option("-c", "--chop", dest="chop", action="store_true",
                      help="chop suffix '_A etc from samples'")
    parser.add_option("--output-mapping", dest="output_mapping", action="store_true",
                      help="output original to new identifier mapping")
    parser.add_option("-t", "--type", dest="type", type = "choice",
                      choices = ("control", "probe"), help="is it control or probe profile?")
    parser.add_option("--file2", dest="file_2", type = "string",
                      help="file to combine to stdin")


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # at the moment the mapping file is output with a default
    # name
    if options.output_mapping:
        assert options.sample_sheet, "must supply a sample sheet to output mapping"
        probe_profile, sample_sheet = options.stdin, options.sample_sheet
        samples = readSampleSheet(sample_sheet, options.chop)
        outputMapping(probe_profile, samples, type="probe")
        E.warn("cannot perform %s if --output-mapping is specified" % options.do)
        # cannot do anything else if this is specified
        return
        
    if options.do == "remove":
        assert options.sample_list, "must provide a list of sample to remove"
        samples = options.sample_list.split(",")
        remove_columns(options.stdin, samples)

    elif options.do == "rename-columns":
        assert options.sample_sheet, "must supply a sample sheet to rename columns"
        probe_profile, sample_sheet = options.stdin, options.sample_sheet
        renameColumns(probe_profile, readSampleSheet(sample_sheet, options.chop), type = options.type)
    
    elif options.do == "reorder-by":
        assert options.sample_list, "must specify a sample list to place in last columns"
        samples = options.sample_list.split(",")
        reorder_columns(options.stdin, samples)

    elif options.do == "combine-files":
        assert options.file_2, "must specify a second file to combine"
        combine_data(options.stdin, options.file_2)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


