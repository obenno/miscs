#! /usr/bin/python3

import sys
import numpy as np
import gffutils
import getopt

def parse_args (argv):
    File_in = ""
    File_out = ""
    outfmt = ""
    try:
        opts, args = getopt.getopt(argv,"hc:g:o:",["cdsfile=","gtffile=","output="])
    except getopt.GetoptError:
        print("mapCDS_to_Gtf.py -c <cds_info> -g <gtf> -o <output>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print ("Useage: input cufflinks gtf, transdecoder cds info\n", "\n",
                   "mapCDS_to_Gtf.py -c <cds_info> -g <gtf> -o <output>", "\n")
            sys.exit()
        elif opt in ("-c"):
            cds_file = arg
        elif opt in ("-g"):
            gtf_file = arg
        elif opt in ("-o"):
            out_file = arg

    print ("CDS file is %s" % cds_file)
    print ("GTF file is %s" % gtf_file)
    print ("Output is %s" % out_file)
    return cds_file, gtf_file, out_file

def make_GFFdb (gtf_file):
    # Define input as gtf or gff
    if ".gtf" in gtf_file:
        # gtf of gffcompare output already has transcript feature,
        # so disable_infer_transcripts was set to True
        db = gffutils.create_db(data=gtf_file,
                                dbfn="tmpGFF.db",
                                force=True,
                                id_spec={"gene": "gene_id", "transcript": "transcript_id"},
                                keep_order=True, merge_strategy="merge",
                                sort_attribute_values=True,
                                disable_infer_transcripts=True)
    elif ".gff" in gtf_file:
        sys.exit("Only GTF file was supported")

def make_cdsList (cds_file):
    cds = open(cds_file, "rt")
    cdsinfo = {}
    for x in cds:
        trans, cds_start, cds_end = x.split("\t")[0], x.split("\t")[3], x.split("\t")[4]
        cdsinfo[trans] = [cds_start, cds_end]
    cds.close()
    return cdsinfo

def get_exonList (GFFdb, TranscriptFeature):
    ## This function return exons list of input isoform
    exons = list(GFFdb.children(TranscriptFeature,
                                featuretype = ("exon"),
                                order_by = "start"))
        #if TranscriptFeature.strand == "-":
    #   exons.reverse()
    return exons

def make_newFeature(chrom, source, featureType, start, end, strand, transcript, gene):
    transcript = "\"" + transcript + "\""
    gene = "\"" + gene + "\""
    annotation = "transcript_id " + transcript + ";" + " " + "gene_id " + gene + ";"
    data = (chrom, source, featureType, str(start), str(end), ".", strand, ".", annotation)
    newline = "\t".join(data)
    f = gffutils.feature.feature_from_line(newline)
    return f

def calc_cds_frame(exonsList, strand):
    ## This functin used input CDS exons list,
    ## calculate frame column for each CDS exon
    ordered_cdsExon = []
    if strand == "+":
        ordered_cdsExon = exonsList
    elif strand == "-":
        ordered_cdsExon = exonsList[::-1]

    frame = 0
    cdsLen = 0
    for i in ordered_cdsExon:
        i.frame = str(frame)
        exonLen = i.stop - i.start +1
        cdsLen = cdsLen + exonLen
        if cdsLen%3 == 0:
            frame = 0
        else:
            frame = 3 - (cdsLen%3)
    return ordered_cdsExon

def Get_cds (exonsList, cds_start, cds_end, strand):
    ## This function is using exons and cds info
    ## to generate cds regions

    fiveUTR_exons = []
    left_UTR = []
    threeUTR_exons = []
    right_UTR = []
    cds_exons = []
    # Get transcript length
    transLen = 0
    for i in exonsList:
        transLen = transLen + i.stop - i.start + 1
    print("Transcript length is ", transLen)
    ## If the strand is -, change CDS start/end to
    ## positive orentation
    print("CDS start is ", cds_start)
    print("CDS end is ", cds_end)
    if strand == "-":
        start = transLen - int(cds_end) + 1
        end = transLen - int(cds_start) + 1
    elif strand == "+":
        start = int(cds_start)
        end = int(cds_end)
    else:
        sys.exit("The strand should be either + or -, please check")

    ## Condition1: start/end not located in the same exon
    ## Get utr/cds boundary exon
    scanedLen = 0
    for i in exonsList:
        exonLen = i.stop-i.start + 1
        scanedLen = scanedLen + exonLen
        if scanedLen-exonLen+1 <= start and \
           scanedLen >= start and \
           scanedLen < end:
            lastUTR_start = i.start
            firstCDS_start = i.start + start - (scanedLen - exonLen) -1
            lastUTR_end = firstCDS_start - 1
            firstCDS_end = i.stop
            lastUTR = make_newFeature(i.chrom, "TransDecoder", "five_prime_UTR",
                                      lastUTR_start, lastUTR_end, strand,
                                      i["transcript_id"][0], i["gene_id"][0])
            left_UTR.append(lastUTR)
            firstCDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                       firstCDS_start, firstCDS_end, strand,
                                       i["transcript_id"][0], i["gene_id"][0])
            cds_exons.append(firstCDS)
            idx_left = exonsList.index(i)
            break

    # Get right_UTR
    scanedLen = 0
    for i in exonsList:
        exonLen = i.stop-i.start + 1
        scanedLen = scanedLen + exonLen
        if scanedLen >= end and \
           scanedLen-exonLen+1 > start and \
           scanedLen-exonLen+1 <= end:
            lastCDS_start = i.start
            lastCDS_end = i.start + exonLen - (scanedLen - end) -1
            firstUTR_start = lastCDS_end +1
            firstUTR_end = i.stop
            firstUTR = make_newFeature(i.chrom, "TransDecoder", "three_prime_UTR",
                                       firstUTR_start, firstUTR_end, strand,
                                       i["transcript_id"][0], i["gene_id"][0])
            right_UTR.append(firstUTR)
            lastCDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                      lastCDS_start, lastCDS_end, strand,
                                      i["transcript_id"][0], i["gene_id"][0])
            cds_exons.append(lastCDS)
            idx_right = exonsList.index(i)
            break
    ## Condition2: start/end located in the same exon
    scanedLen = 0
    for i in exonsList:
        exonLen = i.stop-i.start + 1
        scanedLen = scanedLen + exonLen
        if scanedLen >= end and scanedLen -exonLen+1 <= start:
            lastUTR_start = i.start
            lastUTR_end = i.start + start - (scanedLen - exonLen) -1 -1
            CDS_start = lastUTR_end +1
            CDS_end = CDS_start + (end -start +1) -1
            firstUTR_start = CDS_end +1
            firstUTR_end = i.stop
            lastUTR = make_newFeature(i.chrom, "TransDecoder", "five_prime_UTR",
                                      lastUTR_start, lastUTR_end, strand,
                                      i["transcript_id"][0], i["gene_id"][0])
            firstUTR = make_newFeature(i.chrom, "TransDecoder", "three_prime_UTR",
                                       firstUTR_start, firstUTR_end, strand,
                                       i["transcript_id"][0], i["gene_id"][0])
            CDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                  CDS_start, CDS_end, strand,
                                  i["transcript_id"][0], i["gene_id"][0])
            cds_exons.append(CDS)
            left_UTR.append(lastUTR)
            right_UTR.append(firstUTR)
            idx_left = exonsList.index(i)
            idx_right = exonsList.index(i)
            break
    # Complement CDS exonsList and right_UTR
    for i in exonsList:
        idx = exonsList.index(i)
        if idx < idx_left:
            i.source = "TransDecoder"
            i.featuretype = "five_prime_UTR"
            left_UTR.insert(len(left_UTR)-1, i)
        if idx > idx_left and idx < idx_right:
            i.source = "TransDecoder"
            i.featuretype = "CDS"
            cds_exons.insert(len(cds_exons)-1, i)
        if idx > idx_right:
            i.source = "TransDecoder"
            i.featuretype = "three_prime_UTR"
            right_UTR.append(i)

    if strand == "+":
        fiveUTR_exons = left_UTR
        threeUTR_exons = right_UTR
    elif strand == "-":
        for i in left_UTR:
            i.featuretype = "three_prime_UTR"
        threeUTR_exons = left_UTR
        for i in right_UTR:
            i.featuretype = "five_prime_UTR"
        fiveUTR_exons = right_UTR

    # These codes is to get cds feature frame column
    cds_exons = calc_cds_frame(cds_exons, strand)

    return fiveUTR_exons, threeUTR_exons, cds_exons

def main():
    cds_file, gtf_file, out_file = parse_args(sys.argv[1:])
    cdsinfo = make_cdsList(cds_file)
    outfn = open(out_file, "w+")
    # Make GFFdb, store db in tmpGFF.db file
    make_GFFdb(gtf_file)
    db = gffutils.FeatureDB('tmpGFF.db', keep_order=True)
    for mRNA in db.features_of_type("transcript", order_by="start"):
        # print transcript_id
        print("Processing transcript_id: ", mRNA["transcript_id"][0])
        # print transcript
        print(mRNA, file = outfn)

        # print original exons
        for i in db.children(mRNA, featuretype = "exon",
                             order_by = "start"):
            print(i, file = outfn)

        transcript = mRNA.id
        exons = get_exonList(db, mRNA)
        # Some transcript doesn't have CDS
        # Need to check first
        if transcript in cdsinfo:
            cds_start = cdsinfo[transcript][0]
            cds_end = cdsinfo[transcript][1]
            fiveUTR, threeUTR, cds_exons = Get_cds(exons, cds_start, cds_end, mRNA.strand)
            
            #Disable printing fiveUTR and threeUTR
            #for i in fiveUTR:
            #    print(i, file = outfn)
            #for i in threeUTR:
            #    print(i, file = outfn)
            for i in cds_exons:
                print(i, file = outfn)
    outfn.close()
if __name__ == "__main__":
    main()
