#! /usr/bin/python3

## This script is a modified version of mapCDS_to_Gtf.py to map
## peptides to transcript genomic regions

import sys
import numpy as np
import gffutils
import getopt

def parse_args (argv):
    File_in = ""
    File_out = ""
    outfmt = ""
    try:
        opts, args = getopt.getopt(argv,"hp:g:o:",["pepfile=","gtf_file=","output="])
    except getopt.GetoptError:
        print("mapPep_to_CDS.py -p <pep_info> -g <gtf> -o <output>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print ("Useage: input cufflinks gtf, transdecoder pep info\n", "\n",
                   "mapPep_to_CDS.py -p <pep_info> -g <gtf> -o <output>", "\n")
            sys.exit()
        elif opt in ("-p"):
            pep_file = arg
        elif opt in ("-g"):
            gtf_file = arg
        elif opt in ("-o"):
            out_file = arg

    print ("Peptides file is %s" % pep_file)
    print ("GTF file is %s" % gtf_file)
    print ("Output is %s" % out_file)
    return pep_file, gtf_file, out_file

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
                                disable_infer_transcripts=False)
    elif ".gff" in gtf_file:
        sys.exit("Only GTF file was supported")

def make_pepList (pep_file):
    pepFile = open(pep_file, "rt")
    pepinfo = {}
    ## Please be noted that peptide start and end is the start/end of protein sequences
    ## Need to be converted to transcript genomic coordinates
    for x in pepFile:
        x = x.splitlines()[0]
        pep, trans, pep_start, pep_end = x.split("\t")[0], x.split("\t")[1], x.split("\t")[2], x.split("\t")[3]
        pepinfo[trans] = {}
        pepinfo[trans][pep] = [pep_start, pep_end]
    pepFile.close()
    return pepinfo

def get_exonList (GFFdb, TranscriptFeature, Featuretype):
    ## This function return exons list of input isoform
    exons = list(GFFdb.children(TranscriptFeature,
                                featuretype = Featuretype,
                                order_by = "start"))
        #if TranscriptFeature.strand == "-":
    #   exons.reverse()
    return exons

#def get_cdsList (GFFdb, TranscriptFeature):
#    ## This function return cds list of input isoform
#    cds = list(GFFdb.children(TranscriptFeature,
#                              featuretype = ("CDS"),
#                              order_by = "start"))
#    return cds

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

    ## Frame column is the 8th column of gff3/gtf
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

def Get_cds (exonsList, pep_ID, cds_start, cds_end, strand):
    ## This is cds generating function inherited from
    ## mapCDS_to_Gtf.py, we can use this function to
    ## generate peptide exons like each peptide is
    ## a small CDS
    
    ## This function is using exons and cds info
    ## to generate cds regions;
    ## cds_start and cds_end is position relative
    ## to transcript not genomic coordinates.

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
                                      pep_ID, i["gene_id"][0])
            left_UTR.append(lastUTR)
            firstCDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                       firstCDS_start, firstCDS_end, strand,
                                       pep_ID, i["gene_id"][0])
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
                                       pep_ID, i["gene_id"][0])
            right_UTR.append(firstUTR)
            lastCDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                      lastCDS_start, lastCDS_end, strand,
                                      pep_ID, i["gene_id"][0])
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
                                      pep_ID, i["gene_id"][0])
            firstUTR = make_newFeature(i.chrom, "TransDecoder", "three_prime_UTR",
                                       firstUTR_start, firstUTR_end, strand,
                                       pep_ID, i["gene_id"][0])
            CDS = make_newFeature(i.chrom, "TransDecoder", "CDS",
                                  CDS_start, CDS_end, strand,
                                  pep_ID, i["gene_id"][0])
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

def get_CDS_relativePos(exonsList, cdsList, strand):
    ## This function is used to get relative position of
    ## cds on transcript, return relative start and end

    ## exonsList and cdsList should be ordered according
    ## to transcript orentation
    cds_start = 0
    cds_end = 0
    scanedLen = 0
    cds_relative_end = 0
    cds_relative_start = 0
    if strand == "+":
        cds_start = cdsList[0].start
        cds_end = cdsList[-1].end
        cdsLen = 0
        for i in range(len(exonsList)):
            exon_start = exonsList[i].start
            exon_end = exonsList[i].stop
            exonLen = exon_end - exon_start + 1
            if cds_start>=exon_start and cds_start<=exon_end:
                cds_relative_start = scanedLen + cds_start - exon_start + 1
                cdsLen = exon_end - cds_start + 1
                for p in range(i,len(exonsList)):
                    if cds_end>=exonsList[p].start and cds_end<=exonsList[p].stop:
                        if p == i:
                            cdsLen = cds_end - cds_start + 1
                        elif p > i:
                            cdsLen = cdsLen + cds_end - exonsList[p].start + 1
                        cds_relative_end = cds_relative_start + cdsLen -1
                        break
                    else:
                        if p > i:
                            cdsLen = cdsLen + exonsList[p].stop - exonsList[p].start + 1
                break
            else:
                scanedLen = scanedLen + exonLen
        print("cdsLen", cdsLen)
    elif strand == "-":
        cds_start = cdsList[-1].end
        cds_end = cdsList[0].start
        cdsLen = 0
        exonsList = exonsList[::-1]
        for i in range(len(exonsList)):
            exon_start = exonsList[i].start
            exon_end = exonsList[i].stop
            exonLen = exon_end - exon_start + 1
            if cds_start >= exon_start and cds_start <= exon_end:
                cds_relative_start = scanedLen + exon_end - cds_start + 1
                cdsLen = cds_start - exon_start + 1
                for p in range(i,len(exonsList)):
                    if cds_end >= exonsList[p].start and cds_end <= exonsList[p].stop:
                        if p == i:
                            cdsLen = cds_start - cds_end + 1
                        elif p > i:
                            cdsLen = cdsLen + exonsList[p].stop - cds_end + 1
                        cds_relative_end = cds_relative_start + cdsLen -1
                        break
                    else:
                        if p > i:
                            cdsLen = cdsLen + exonsList[p].stop - exonsList[p].start + 1
                break
            else:
                scanedLen = scanedLen + exonLen
        print ("cdsLen", cdsLen)
    return cds_relative_start, cds_relative_end

def pep_position_to_trans_relativePOS (pep_start, pep_end, cds_relative_start, cds_relative_end):
    pep_start = int(pep_start)
    pep_end = int(pep_end)
    cds_relative_start = int(cds_relative_start)
    cds_relative_end = int(cds_relative_end)
    trans_start = cds_relative_start+(pep_start - 1)*3
    trans_end = trans_start + (pep_end - pep_start + 1)*3 - 1
    return trans_start, trans_end

def main():
    pep_file, gtf_file, out_file = parse_args(sys.argv[1:])
    pepinfo = make_pepList(pep_file)
    #print(pepinfo)
    outfn = open(out_file, "w+")
    # Make GFFdb, store db in tmpGFF.db file
    make_GFFdb(gtf_file)
    db = gffutils.FeatureDB('tmpGFF.db', keep_order=True)
    for mRNA in db.features_of_type("transcript", order_by="start"):
        # print transcript_id
        print("Processing transcript_id: ", mRNA["transcript_id"][0])
        # print transcript
        #print(mRNA, file = outfn)
        transcript = mRNA.id
        exons = get_exonList(db, mRNA, "exon")
        cds = get_exonList(db, mRNA, "CDS")
        # Some transcript doesn't have CDS
        # Need to check first
        cds_relative_start, cds_relative_end = get_CDS_relativePos(exons, cds, mRNA.strand)
        print("relative position: ", cds_relative_start, " ", cds_relative_end)
        if transcript in pepinfo:
            for pep in pepinfo[transcript]:
                pep_start = pepinfo[transcript][pep][0]
                pep_end = pepinfo[transcript][pep][1]
                print("pep_start is ", pep_start)
                print("pep_end is", pep_end)
                ## Covert pep start and end to relative position on transcript
                trans_start, trans_end = pep_position_to_trans_relativePOS(pep_start, pep_end, cds_relative_start, cds_relative_end)
                print("trans_start", trans_start)
                print("trans_end", trans_end)
                fiveUTR, threeUTR, pep_exons = Get_cds(exons, pep, trans_start, trans_end, mRNA.strand)

                #Disable printing fiveUTR and threeUTR
                #for i in fiveUTR:
                #    print(i, file = outfn)
                #for i in threeUTR:
                #    print(i, file = outfn)
                print(pep_exons)
                for i in pep_exons:
                    print(i, file = outfn)
    outfn.close()
if __name__ == "__main__":
    main()
