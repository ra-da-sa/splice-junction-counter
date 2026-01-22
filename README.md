# splice-junction-counter
This script processes mRNA sequencing alignments (.sam) to identify splice junctions. The output is a tab-delimited file (output.txt) containing gene IDs, junction start coordinates, junction end coordinates, and supporting read counts.

## input
This script takes two input files:
1. .sam file containing aligned mRNA sequencing reads. The file should have the following SAM columns: Column 3 (RNAME): chromosome of the alignment, Column 4 (POS): leftmost alignment positionColumn 6 (CIGAR): CIGAR string describing the alignment, Last column (NH:i:x): number of alignment locations for the read
3. .txt file containing gene information. Thefile should have gene IDs and their genomic coordinates.

## usage
The script is run using the following command format:
```
python splice-junction-py <SamFile.sam> <GeneLocation_Summary.txt>
```
## output
The script produces a single output file called output.txt, which contains the following columns:
1. Gene ID
2. Junction start (gene coordinate)
3. Junction end (gene coordinate)
3. Number of reads supporting the junction

## known issues
1. The script only checks file extensions and does not validate the actual content of the input files, assuming that the user has saved each file with the correct extension.
2. The script should implement argparse to handle command-line arguments.
3. The script should implement Logging to verify each step of the calculation and to record any errors.
