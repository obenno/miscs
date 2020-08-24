#! /usr/bin/env python


## This is for converting genome-based GFF to transcript-based GFF
## in case you used transcripts fasta as mapping reference, for example
## nanopore direct RNA sequencing project

import os
import sys, getopt
import re
import datetime
import gffutils

def parse_args (argv):
    File_in = ""
    File_out = ""
    outfmt = ""
    try:
        opts, args = getopt.getopt(argv,"hi:o:f:",["ifile=","ofile=","outfmt="])
        if len(opts)==0:
            print("genomeGFF_to_transcriptGFF.py -i <input> -o <output> -f <gtf/gff>")
            sys.exit(2)
        for opt, arg in opts:
            if opt == "-h":
                print ("Useage: input cufflinks gtf, output gtf or gff\n", "\n",
                       "ExtractALEandAFE.py -i <input> -o <output> -f <outfmt>")
                sys.exit()
            elif opt in ("-i"):
                File_in = arg
            elif opt in ("-o"):
                File_out = arg
            elif opt in ("-f"):
                outfmt = arg

        print ("Input file is %s" % File_in)
        print ("Output file is %s" % File_out)
        print ("Output format is %s" % outfmt)
    except getopt.GetoptError:
        print("genomeGFF_to_transcriptGFF.py -i <input> -o <output> -f <gtf/gff>")
        sys.exit(2)

    return File_in, File_out, outfmt


def make_GFFdb (File_in):
    # Define input as gtf or gff
    if re.search("\.gff.", File_in):
        # Change id_spec if your GFF doesn't follow default ID setting
        # disable_infer_transcripts was set to True
        db = gffutils.create_db(data=File_in,
                                dbfn="tmpGFF.db",
                                force=True,
                                id_spec="ID",
                                keep_order=True,
                                merge_strategy="error", # Throw an error when encountering duplicate ID
                                sort_attribute_values=True,
                                disable_infer_transcripts=True)
    elif ".gtf" in File_in:
        print("For now, only phytozome gff3 file is supported")
        sys.exit(1)

def infer_transcriptGFF (dbinput):
    ## Read input db file
    db = gffutils.FeatureDB(dbinput, keep_order=True)
    outGFF=[]
    for transcript in db.all_features(featuretype="mRNA", order_by="seqid"):
        transcript_id=transcript.attributes['ID'][0]
        transcript_len=0
        five_UTR_len=0
        CDS_len=0
        three_UTR_len=0
        outline=""
        anno_id=""
        parent_id=""
        for exon in db.children(transcript, featuretype="exon", order_by="start"):
            transcript_len += exon.end-exon.start+1
        for five_utr in db.children(transcript, featuretype="five_prime_UTR", order_by="start"):
            five_UTR_len += five_utr.end-five_utr.start+1
        for cds in db.children(transcript, featuretype="CDS", order_by="start"):
            CDS_len += cds.end-cds.start+1
        three_UTR_len = transcript_len - five_UTR_len - CDS_len
        # Gene line
        anno_id= "".join(["ID=", "Gene.", transcript_id])
        outline="\t".join([transcript_id, "Transcript", "gene", "1", str(transcript_len), ".", "+", ".", anno_id])
        outGFF.append(outline)
        # mRNA line
        anno_id= "".join(["ID=", "mRNA.", transcript_id])
        parent_id = "".join(["Parent=", "Gene.", transcript_id])
        outline="\t".join([transcript_id, "Transcript", "mRNA", "1", str(transcript_len), ".", "+", ".", ";".join([anno_id, parent_id])])
        outGFF.append(outline)
        # 5UTR line
        anno_id= "".join(["ID=", "fiveUTR.", transcript_id])
        parent_id = "".join(["Parent=", "mRNA.", transcript_id])
        outline="\t".join([transcript_id, "Transcript", "five_prime_UTR", "1", str(five_UTR_len), ".", "+", ".", ";".join([anno_id, parent_id])])
        outGFF.append(outline)
        # CDS line
        anno_id= "".join(["ID=", "cds.", transcript_id])
        parent_id = "".join(["Parent=", "mRNA.", transcript_id])
        outline="\t".join([transcript_id, "Transcript", "CDS", str(five_UTR_len+1), str(five_UTR_len+CDS_len), ".", "+", ".", ";".join([anno_id, parent_id])])
        outGFF.append(outline)
        # 3UTR line
        anno_id= "".join(["ID=", "threeUTR.", transcript_id])
        parent_id = "".join(["Parent=", "mRNA.", transcript_id])
        outline="\t".join([transcript_id, "Transcript", "three_prime_UTR", str(five_UTR_len+CDS_len+1), str(transcript_len), ".", "+", ".", ";".join([anno_id, parent_id])])
        outGFF.append(outline)
    
    return outGFF

def main ():
    File_in, File_out, outfmt = parse_args(sys.argv[1:])
    make_GFFdb(File_in)
    outGFF = infer_transcriptGFF("tmpGFF.db")
    outfile=open(File_out,"w")
    for i in outGFF:
        print(i, file = outfile)
    outfile.close()
    os.remove("tmpGFF.db")

if __name__ == "__main__":
    main()