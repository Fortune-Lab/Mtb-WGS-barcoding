from needletail import parse_fastx_file, NeedletailError, reverse_complement
import re

#(GCGGCCGCGAATTCCG[AG])([ATCG]{16,19})([GC]AATTCGATGGCCT)/)

barcodes = {}

try:
    for record in parse_fastx_file("LIB065162_GEN00299530_G2G_S1_L001_R1_001.fastq.gz"):
        pattern = r'(GCGGCCGCGAATTCCG[AG])([ATCG]{16,19})([GC]AATTCGATGGCCT)'
        match = re.search(pattern, record.seq)
        if match:
            if match.group(2) in barcodes:
                barcodes[match.group(2)] += 1
                #print(match.group(2), barcodes[match.group(2)])
            else:
                barcodes[match.group(2)] = 1
        #print(record.id)
        #print(record.seq)
        #print(record.qual)
except NeedletailError:
    print("Invalid Fastq file")

barcodes_sorted = dict(sorted(barcodes.items(), key=lambda item: item[1], reverse=True))

for key, value in barcodes_sorted.items():
    print(f"key: {key}, value: {value}")
