import sys
import re 

# =========================
# part 1: get user input at the command line 
# =========================
# initialise placeholders for the input files that will be used (.sam and .txt)
sam_file = None 
gene_locations_file = None

if len(sys.argv) != 3: # ensure the user provides exactly two input files at the command line
    print('python3 script.py <SamFile.sam> <GeneLocation_Summary.txt>') # show user the expected format of an input
    sys.exit(1)

# loop through the input arguments and assign them to the correct placeholder based on the file extension
for arg in sys.argv[1:]:  
    # check if the input ends with .sam  and assign it to sam_file
    if arg.endswith('.sam'): 
        sam_file = arg
    # check if the input ends with .txt and assign it to gene_locations_file
    elif arg.endswith('.txt'):
        gene_locations_file = arg
    else:
        print(f"Error: Unknown file type '{arg}'. Expected a .sam or .txt file.") # print an error if the input is not in .sam or .txt file type
        print("Usage: python3 script.py <SamFile.sam> <GeneLocation_Summary.txt>") # ask the user to provide a new input with valid file type
        sys.exit(1)

# =========================
# part 2A: define a function to extract intron start and end positions from a CIGAR string
# =========================
# extract_intron takes 2 arguments: 
#   1. cigar: a CIGAR string of a single alignment (e.g., "10M5N20M")
#   2. start_position: a position on the chromosome where the alignment starts
def extract_intron(cigar, start_position):
    introns = []  # create a list to store (start, end) positions of each intron

    matches = re.finditer(r'(\d+)([MDN])', cigar) # split the CIGAR string into length and operation groups: "10M" → (10)(M) (see line 42-3)
    current_position = start_position # set start_position as the current_position to track our location on the chromosome as we process the CIGAR string
    
    # loop through each CIGAR operation in the string
    for match in matches: 
        cigar_length = match.group(1) # group 1: operation length
        cigar_type = match.group(2) # group 2: operation type (M, D, or N)

        # if the operation is M (match) or D (deletion),
        # add the length of operation to the current_position to move along the chromosome
        if cigar_type == 'M' or cigar_type == 'D':
            current_position = current_position + int(cigar_length)
        
        # If the operation is N (skipped region), 
        # record the start and end positions of the intron
        if cigar_type == 'N':
            intron_start = current_position # set current_position as the start position of intron
            intron_end = current_position + int(cigar_length) # add operation length to current_position to get the end position of intron
            introns.append((intron_start, intron_end)) # record the intron start and end positions as a tuple
            current_position = intron_end # set intron_end position as the current_position to move along the chromosome
    
    return introns # return a list of (start, end) positions for all introns in a CIGAR string

# =========================
# part 2B: test extract_intron with assertions
# =========================
# test a simple case with one intron (M and D mixed)
assert extract_intron("10M5N20D", 100) == [(110, 115)]

# test a complex case with multiple introns (M and D mixed)
assert extract_intron("10D5N20M3N7D", 100) == [(110, 115), (135, 138)]

# test a complex case with no introns (M and D mixed with no N)
assert extract_intron("10M5D20M", 100) == []

# =========================
# part 3: extract junctions from the CIGAR string of all reads in the sam_file
# =========================
# note 1: only extract junctions from reads that meet two conditions:
#   condition 1: the read only aligns once (NH:i:1)
#   condition 2: the read is split (contains at least one 'N' in the CIGAR string) 
# note 2: when a read is split (contains 'N'), each split region (intron) can be referred to as a junction

junctions = [] # create a list to store (chromosome, start position, end position) for each junction

# open sam_file
try: 
    with open(sam_file) as input: 

        # loop through each read (line) in the sam_file
        for line in input: 
            if line.startswith('@'): # skip header rows starting with '@'
                continue
                
            else:
                line = line.rstrip() # removes whitespace at the end of line
                line = line.split('\t')  # split the line into columns; .sam file is tab-delimited (\t)

                if line[-1] == 'NH:i:1': # check if the read passes condition 1
                    chromosome = line[2] # column 3 contains the chromosome where the read aligned; store in variable 'chromosome'
                    alignment_start = int(line[3]) # column 4 contains the position where the alignment starts; store in variable 'alignment_start'
                    cigar = str(line[5]) # column 6 contains the CIGAR string describing the alignment; store in variable 'cigar'
                    
                    if 'N' in cigar: # check if the read pass condition 2
                        junctions_per_cigar = extract_intron(cigar, alignment_start) # extract junction(s) from each CIGAR string
                        
                        # record the chromosome, start, and end positions of each junction in the list 'junctions' (line 81)
                        for start_position, end_position in junctions_per_cigar: 
                            junctions.append((chromosome, start_position, end_position))

# if sam_file cannot be opened, print an error message and exit
except OSError:
    sys.exit(f'File {sam_file} cannot be opened. Check if the file exists and is readable. Please try again.')

# check the output
# print(junctions)  
# expected output: [('TGME49_chrVIII', 6947299, 6947987)]

# =========================
# part 4: count the number of reads supporting each unique junction
# =========================
junctions_counts_dict = {} # create a dict to store junctions as keys and their counts as values {(chromosome, start, end): count}

# loop through each junction in the list 'junctions'
for junction in junctions: 
    # if the junction is not in the dict, store junction in the dict and set its count to 1
    if junction not in junctions_counts_dict:
        junctions_counts_dict[junction] = 1

    # if the junction is already in the dict, increase its count by 1
    else:
        junctions_counts_dict[junction] += 1 

# check the output
# print(junctions_counts_dict)
# expected output: {('TGME49_chrVIII', 6947299, 6947987): 2}

# =========================
# part 5: extract the location (chromosome, start position, end position) for each gene
# =========================
genes_dict = {} # create a dict to store gene as keys and their locations as values {gene: (chromosome, start position, end position)}

# open gene_locations_file
try:
    with open(gene_locations_file) as input: 

        # loop through each gene (line) in the file
        for line in input: 
            line = line.rstrip()
            line = line.split('\t') # split the line into columns; gene_locations_file is tab-delimited (\t)
            gene_id = line[0] # column 1 contains the gene ID; store in variable 'gene_id'
            gene_location = line[2] # column 3 contains the gene location; store in variable 'gene_location'

            # parse the gene location string into 5 groups: "TGME49_chrVIII:6,631,349..6,636,865(+)" → (chromosome)(:)(start)(..)(end) 
            match = re.search(r'(\w+)([:])([\d,]+)([.]{2})([\d,]+)', gene_location) 
            if match:
                chromosome = match.group(1) # group 1: chromosome
                start_position = int(match.group(3).replace(',', '')) # group 2: start position of the gene
                end_position = int(match.group(5).replace(',', '')) # group 3: end position of the gene

                #  store each gene and its location (chromosome, start position, end position) in the dict 'genes_dict'
                genes_dict[gene_id] = (chromosome, start_position, end_position) 

# if gene_locations_file cannot be opened, print an error message and exit
except FileNotFoundError:
    sys.exit(f'File {gene_locations_file} cannot be opened. Check if the file exists and is readable. Please try again.')

# check the output
# print(genes_dict)
# expected output: {'TGME49_200480': ('TGME49_chrVIII', 6916849, 6917687)}

# =========================
# part 6: write an output table recording gene ID, junction start, junction end, and number of reads supporting each junction (4-column table)
# =========================
with open('output.txt', 'w') as output: # open output file (hardcode as student number)
    current_gene_id = None  # create a variable 'current_gene_id' to track the gene ID while looping through junctions
    
    # loop through each unique junction and its read count in the dict 'junctions_counts_dict'
    for (junction_chr, junction_start, junction_end), read_count in junctions_counts_dict.items(): 

        # loop through each gene and its location in the dict 'genes_dict' 
        for gene_id, (gene_chr, gene_start, gene_end) in genes_dict.items(): 
            # check if:
            #   1. the junction is on the same chromosome as the gene
            #   2. the junction lies within the gene boundaries
            if junction_chr == gene_chr and junction_start >= gene_start and junction_end <= gene_end:

                # add an empty line in the output before continuing to the next gene
                if current_gene_id is not None and current_gene_id != gene_id:
                    output.write('\n')

                # write gene ID, junction start, junction end, and read count as a tab-delimited line
                output.write(f'{gene_id}\t{junction_start}\t{junction_end}\t{read_count}\n')

                # update the current_gene_id tracker
                current_gene_id = gene_id
