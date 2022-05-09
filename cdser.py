#imports
import sys

import antismash
from Bio import SeqIO

def strip_duplicate_cds(record):
    locations = set()
    features = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        if str(feature.location) in locations:
            continue
        features.append(feature)
        locations.add(str(feature.location))
    record.features = features
    return record

def find_cds(record):
    cdses = []
    record = antismash.common.secmet.Record.from_biopython(record, "bacteria")
    for cds in record.get_cds_features():
        cdses.append(("%s|%s" % (record.id, cds.get_name()), cds))
    with open("cdses.fasta", "a+") as handle:
        for name, cds in cdses:
            handle.write(">%s\n%s\n" % (name, cds.location.extract(record.seq)))
    with open("cdses.faa", "a+") as handle:
        for name, cds in cdses:
            handle.write(">%s\n%s\n" % (name, cds.translation))

def main(genbank):
    record = strip_duplicate_cds(list(SeqIO.parse(genbank, "genbank"))[0])
    find_cds(record)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: genbank.gbk')
    else:
        main(sys.argv[1])
