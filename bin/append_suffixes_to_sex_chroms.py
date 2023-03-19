#!/usr/bin/env python3

import argparse
#from tqdm import tqdm
from g2gtools import gtf
from collections import namedtuple, defaultdict

gtfInfoFields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
GTFRecord = namedtuple('GTFRecord', gtfInfoFields)

ATTRIBUTES_TO_ALTER = {'gene_id':'gene_id', 'transcript_id':'transcript_id', 'exon_id':'exon_id', 'protein_id':'protein_id', 'ccds_id':'ccds_id'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Standardize naming of X/Y gene to ensure overlapping genes are treated appropriately by downstream software."""
    )
    parser.add_argument("-i", type=str, help="haploid_gtf")
    parser.add_argument('-o', type=str, help='output gtf filename')
    parser.add_argument('-a', type=bool, help='add L/R indicator to X and Y ids in gtf')
    parser.add_argument("-X", type=str, help="X chromosome ID")
    parser.add_argument("-Y", type=str, help="Y chromosome ID")
    parser.add_argument('--progress_bar', action='store_true')

    args = parser.parse_args()
    return args

def GTFRecord_to_str(r):
    string_r = [str(x) if x is not None else '.' for x in r]
    string_r[-1] = '; '.join([f"{str(k)} \"{str(v)}\"" for k,v in r.attributes.items()])
    return '\t'.join(string_r)+'\n'


def parse_gtf(gtf_file, XY_suffix, progress_bar=False):
    print('parsing gtf...')
    gtf_db = gtf.GTF(gtf_file)
    relations = defaultdict(lambda: defaultdict(lambda: {'exon':list(), 'CDS':list()}))
    entries = {'gene':dict(), 'transcript':dict(), 'exon':dict(), 'CDS':dict()}
    entriestosort = {'gene':list(), 'transcript':list()}
    Y = list(XY_suffix.keys())[1]
    all_records = [x for x in gtf_db]
    
    for record in all_records:
        d = record._asdict()
        if record.seqid == Y:
            for key, value in record.attributes.items():
                if isinstance(value, str):
                    if value[-2:] == '_1':
                        record.attributes[key] = value[:-2]
                    if '_1-' in value:
                        record.attributes[key] = value.replace('_1-', '-')
        # exon_id fix - every exon needs an exon_id, but in most cases the exon_id is not included
        if set(('exon_number', 'transcript_id')).issubset(set(record.attributes.keys())):
            if d['type'] in ('exon', 'CDS'):
                record.attributes['exon_id'] = record.attributes['transcript_id'] + '-' + record.attributes['exon_number']
            if d['type'] == 'CDS':
                record.attributes['ccds_id'] = 'cds-' + record.attributes['exon_id']
        for field in ['gene_id', 'transcript_id', 'exon_id', 'protein_id', 'ccds_id']:
            if (field in record.attributes.keys()) and (record.seqid in XY_suffix.keys()) and (record.attributes[field][-2:] not in XY_suffix.values()):
                record.attributes[field] = record.attributes[field] + XY_suffix[record.seqid]
        if record.type == 'gene':
            id = record.attributes['gene_id']
        elif record.type == 'transcript':
            id = record.attributes['transcript_id']
        elif record.type == 'exon':
            id = record.attributes['exon_id']
        elif record.type == 'CDS':
            id = record.attributes['ccds_id']
        else:
            continue
        
        entries[record.type][id] = record
        if record.type in ('exon', 'CDS'):
            relations[record.attributes['gene_id']][record.attributes['transcript_id']][record.type].append(id)
        elif record.type == 'gene':
            entriestosort[record.type].append((id, record.start))
    return relations, entries, entriestosort


def sort_gtf(gtf, progress_bar):
    print ('sorting gtf...')
    relations, entries, entriestosort = gtf
    
    for k,v in entriestosort.items():
        entriestosort[k] = sorted(v, key=lambda x: x[1])

    records_to_print = list()
#    print ('sorting gene records...')
#    if progress_bar: 
#       genes = tqdm(entriestosort['gene'])
#    else:
    genes = entriestosort['gene']
    for gene_id, gene_start in genes:
        records_to_print.append(entries['gene'][gene_id])
        is_negative_strand = entries['gene'][gene_id].strand=='-'
        transcripts = [entries['transcript'][id] for id in relations[gene_id].keys()]
        transcripts = sorted(transcripts, key = lambda x: x.start, reverse=is_negative_strand)
        for transcript_record in transcripts:
            transcript_id = transcript_record.attributes['transcript_id']
            records_to_print.append(transcript_record)
            exons = [entries['exon'][id] for id in relations[gene_id][transcript_id]['exon']]
            exons = sorted(exons, key = lambda x: x.attributes['exon_number'], reverse=is_negative_strand)
            for exon in exons:
                records_to_print.append(exon)
            CDSs = [entries['CDS'][id] for id in relations[gene_id][transcript_id]['CDS']]
            CDSs = sorted(CDSs, key = lambda x: x.attributes['exon_number'], reverse=is_negative_strand)
            for CDS in CDSs:
                records_to_print.append(CDS)
    return records_to_print

def main():
    '''
    Going to try to handle sex chromosomes by treating Y chrom as effectively an X chromosome w/ many very large structural variants    
    If I were to do that, I'd need to modify the GTF as so below
    But I'd also need to convert the X and Y VCI into one record. That... honestly seems hard.
    Maybe I can try something closer to my original plan, and just "fix" the PAR X/Y trancripts that have the same ID but different contigs (which is confusing extract?)
    Yes. Proceed as currently doing. But annotation (GTF file) is going to need to ensure that genes in the PAR have _L and _R suffixes to allow for downstream processing in manner similar to autosomes.
    Once we hit FISHPOND and/or SEESAW, we can ensure that SEESAW only imports allelic counts from PAR by filtering for genes/transcripts with _L and _R suffixes (which SEESAW might even do already, but should be straightforward to do on my own.)
    ONLY trick is - do I need to change the contig info? I think that's the one part I don't want to change.
    '''
    args = parse_args()

    gtf_file = args.i
    add_suffixes = args.a
    if add_suffixes:
        XY_suffix = {args.X:'_L', args.Y:'_R'} # will have to convert back before use, unless SEESAW expects a similar setup? 
    else:
        XY_suffix = {args.X:'', args.Y:''}

    gtf_db = parse_gtf(gtf_file, XY_suffix, progress_bar=args.progress_bar)
    with open(args.o, 'w') as out_gtf:
        for record in sort_gtf(gtf_db, args.progress_bar):
            _=out_gtf.write(GTFRecord_to_str(record))



# X = 'NC_060947.1'; Y = 'NC_060948.1'; XY_suffix = {X:'_L', Y:'_R'}

if __name__=='__main__':
    main()
