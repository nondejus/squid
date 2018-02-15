#!/usr/bin/env python3
import re, sys
import argparse
import gzip
from collections import Counter      
from datetime import datetime

def generate_pri_mRNA(input_gtf, region_id, region_identifier):
    """ 
    iterate through gtf. at regions (i.e. transcript) extract
    out identifier (i.e. transcript_id) and generate full
    pre-mRNA records. Append primary_transcript to identifier.
    """
     
    # build dicitonary of transcripts
    exon_dict = count_exons(input_gtf, region_identifier)

    # keep transcripts with more than 1 exon
    multi_exon_dict = {k: v for k, v in exon_dict.items() if v > 1}

    input_gtf.seek(0)
    
    # add header to output gtf indicating that it's been modified
    print("##primary transcripts were added to the annotations using add_primary_transcripts.py")
    print("##date: " + datetime.now().strftime('%Y-%m-%d')) 

    for index, record in enumerate(input_gtf):

        # preserve all original records
        print(record, end = "")
        
        if record.startswith("#"):
            continue
        
        line = record.split("\t")
        
        region = line[2]
        if region != region_id: 
            continue 


        attributes = line[8]
        attributes = attributes.rstrip('\n'); 

        pattern = region_identifier + ' \"([^;]+)\";'
        
        m = re.search(pattern, attributes)

        if m is None: 
            print("unable to find {} at line {}".format(region_identifier,
                str(index + 1)), file = sys.stderr)
            continue

        else:
            attr_id = m.group(1)
                
        # only make primary_transcript record if multi-exonic
        
        if attr_id not in multi_exon_dict:
            continue
        
        
        # output a new transcript record
        new_attr_id = "pre_" + attr_id
        
        if region_id == "gene":
            new_attributes = attributes.replace(";", 
                    '; transcript_id \"' + new_attr_id + '";', 
                    1) 
        else:
            new_attributes = attributes.replace(attr_id, new_attr_id) 

        line[2] = 'transcript'
        line[8] = new_attributes
        print("{}".format("\t".join(line)))
        
        # then make an exon record
        line[2] = "exon"
        new_attributes = new_attributes.rstrip();  
        new_attributes = new_attributes + " exon_number \"1\";"
        new_attributes = new_attributes + " exon_id \"" + new_attr_id + ".exon1\"" + ";"

        line[8] = new_attributes
        print("{}".format("\t".join(line)))

def count_exons(input_gtf, region_identifier):
    """ take open file object,
    iterate through and build a dictionary of number 
    of exons per transcript id 
    """

    exon_dict = Counter()
    
    for index, record in enumerate(input_gtf):

        if record.startswith("#"):
            continue
        
        line = record.split("\t")
        
        region = line[2]
        if region != "exon": 
            continue 

        attributes = line[8]

        pattern = region_identifier + ' \"([^;]+)\";'
        
        m = re.search(pattern, attributes)

        if m is None: 
            print("unable to find {} at line {} while counting exons".format(region_identifier,
                str(index + 1)), file = sys.stderr)
            continue

        else:
            exon_dict[m.group(1)] += 1

    return exon_dict

def main():

    parser = argparse.ArgumentParser(description = """
           Add primary transcript (pre-mRNA) records 
           into GTF annotations. Used for estimating
           mature vs. pre-mRNA levels using RSEM/Salmon """)

    parser.add_argument('-i','--input_gtf', 
            help = """input gtf file""", required = True)
    
    parser.add_argument('-r','--region', 
            help = """region to use for defining 
            region coordinates, (gene or transcript) (default =
            'transcript') """, 
            default = "transcript",
            required = False)
    
    parser.add_argument('-a','--identifier', 
            help = """attribute to extract to use for defining 
            region name (default = 'transcript_id') """, 
            default = "transcript_id",
            required = False)
    
    args = parser.parse_args()
    
    # use either gzipped or plain text input
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith('.gz') else open(f) 
    
    input_gtf = gzopen(args.input_gtf)
    
    generate_pri_mRNA(input_gtf, args.region, args.identifier)

    input_gtf.close() 

if __name__ == '__main__': main()
