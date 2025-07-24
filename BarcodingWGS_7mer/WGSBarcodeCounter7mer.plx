use strict;
use warnings;


# Get command line arguments: run name, read1 fastq, read2 fastq
my($run, $fastq1, $fastq2) = @ARGV;

# Define the expected L5 sequences and their labels
my %L5 = (
        "CCTAGTCTTCCAAACTAGCTACGCGGGTTCGATTCCCGTCGCCCGCTCCGCTGGTCAGA" => 'insertion_site',
        "GAGCCCTAGTCTTCCAAACTAGCGACGCGGGTTCGATTCCCGTCGCCCGCTCGGGCCATG" => 'minus_insertion',
);

# Build regex alternation for L5 motif
my $L5 = "(CCTAGTCTTCCAAACTAGCTACGCGGGTTCGATTCCCGTCGCCCGCTCCGCTGGTCAGA|GAGCCCTAGTCTTCCAAACTAGCGACGCGGGTTCGATTCCCGTCGCCCGCTCGGGCCATG)";

# Regex pattern to capture barcode motif (barcode extracted as $2 by match)
my $BARCODE_MOTIF = "(TACCCC?GA([ACTG]+)AATTCGATGGCCTA)";

# Concatenated regex for qtag sequence IDs
my $QTAGS= "(TGGTGTTCAAGCTTTCGGCTAGATGT|TGGTGTTCAAGCTTAGGAACACCAAG|TGGTGTTCAAGCTTTCGCCGAGCAGT|TGGTGTTCAAGCTTCGAGCGCGAGGA|TGGTGTTCAAGCTTTGGCGAATATGG|TGGTGTTCAAGCTTTCTTCTACAACA|TGGTGTTCAAGCTTAGCACGCCTTGT|TGGTGTTCAAGCTTGCAACTTCTTCA|TGGTGTTCAAGCTTAAGAAGTCCAAC)";
my %qtags = (
        "TGGTGTTCAAGCTTTCGGCTAGATGT" => "19",
        "TGGTGTTCAAGCTTAGGAACACCAAG" => "23",
        "TGGTGTTCAAGCTTTCGCCGAGCAGT" => "22",
        "TGGTGTTCAAGCTTCGAGCGCGAGGA" => "24",
        "TGGTGTTCAAGCTTTGGCGAATATGG" => "25",
        "TGGTGTTCAAGCTTTCTTCTACAACA" => "26",
        "TGGTGTTCAAGCTTAGCACGCCTTGT" => "27",
        "TGGTGTTCAAGCTTGCAACTTCTTCA" => "26_2",
        "TGGTGTTCAAGCTTAAGAAGTCCAAC" => "17",
);

#------------------------------------------------------------#
# Calculate Hamming distance of two equal-length DNA strings
#------------------------------------------------------------#
sub hamming{
     # String length is assumed to be equal
     my ($p,$o) = @_;
     return 0 if $p eq $o;
     my $len = length ($p);
     my $num_mismatch = 0;
     for (my $i=0; $i<$len; $i++){
        ++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
     }
     return $num_mismatch;
}

#------------------------------------------------------------#
# Unused: Converts tab-separated sequences into a hash structure
# Used nowhere in current script.
#------------------------------------------------------------#
sub get_files {
        my %hash = ();
        my $string = shift;
        my @array = split/\t/, $string;
        foreach my $seq (@array[1,2]) {
                push @{$hash{$array[0]}}, $seq;
        }
        return %hash;
}

#------------------------------------------------------------#
# Main file handle: Open FASTQ files converted on-the-fly to FASTA via seqtk.
# Store filehandles in $fh1 (read1), $fh2 (read2)
#------------------------------------------------------------#
open my ($fh1), "seqtk seq -A $fastq1 |" or die $!;
open my ($fh2), "seqtk seq -A $fastq2 |" or die $!;
my $header = "undef";

#------------------------------------------------------------#
# MAIN LOOP: Iterate through both files simultaneously, per sequence entry
#------------------------------------------------------------#
while (!eof($fh1) and !eof($fh2)) {
        my $read1 = <$fh1>;
        my $read2 = <$fh2>;
        # If line is a FASTA header, parse for ID and continue
        if ($read1 =~ /^\>/) {
                $header = get_header_info($read1);
                next;
        }
        # Identify features for each read
        my @read1 = feature_matches($read1);
        my @read2 = feature_matches($read2);

        # If both reads had 'undef' for all three features, skip this entry
        next if (count_item(@read1) == '3' and count_item(@read2) == '3');

        # Check if both reads found a junction feature, try to merge
        if (my $jnc = merge_junction($read1[2],$read2[2])) {
                # Output in special format for junctions
                print join("\t", $run, $header, $jnc), "\n";
                next;
        }
        # Combine barcodes and qtags from both reads
        my %hash = combine_reads($read1[0], $read1[1], $read2[0], $read2[1]);

        # Flag if there are two different barcodes from the pair
        my $mm = "FALSE";
        if (scalar(keys %hash) == '2') {
                $mm = "TRUE";
        }
        # Output: for each barcode->[qtag] combo
        foreach my $bc (keys %hash) {
                if (scalar(keys %{$hash{$bc}}) == '2') {
                        delete $hash{$bc}{'undef'};
                }
                foreach my $qtag (keys %{$hash{$bc}} ) {
                        print join("\t", $run, $header, $bc, $qtag, $mm), "\n";
                }
        }
}

#------------------------------------------------------------#
# Given four values (read1_barcode, read1_qtag, read2_barcode, read2_qtag)
# return a hash: { barcode => { qtag => count } }
#------------------------------------------------------------#
sub combine_reads {
        my @array = @_;
        my %hash = ();
        if ($array[0] ne 'undef') {
                $hash{$array[0]}{$array[1]}++;
        }
        if ($array[2] ne 'undef') {
                $hash{$array[2]}{$array[3]}++;
        }
        return %hash;
}

#------------------------------------------------------------#
# Merge two barcode calls -- handle cases for presence/absence or conflict.
# Only actually used for barcodes, but currently not called in the main code!
#------------------------------------------------------------#
sub merge_barcode_reads {
        my $a = shift;
        my $b = shift;
        my $c = shift;
        my $d = shift;
        if ($a eq 'undef' and $b =~ /[ATGC]+/) {
                return $b;
        }elsif ($a =~ /[ATCG]+/ and $b eq 'undef') {
                return $a;
        }elsif (($a =~ /[ATCG]+/ and $b =~ /[ATCG]+/) and ( $a eq $b) ) {
                return $a;
        }elsif (($a =~ /[ATCG]+/ and $b =~ /[ATCG]+/) and ( $a ne $b) ) {
                return $a . '_' . $b;
        }
        return "";
}

#------------------------------------------------------------#
# Merge L5 junction calls: return a non-undef value if only one is present,
# or return call if both are the *same* insertion/label type.
#------------------------------------------------------------#
sub merge_junction{
        my $a = shift;
        my $b = shift;
        if ($a eq 'undef' and $b =~ /insertion/) {
                return $b;
        }elsif ($a =~ /insertion/ and $b eq 'undef') {
                return $a;
        }elsif ($a =~ /insertion/ and $b =~ /insertion/) {
                return $a;
        }
        return "";
}

#------------------------------------------------------------#
# Count number of 'undef' items in an array (used for skipping)
#------------------------------------------------------------#
sub count_item {
        my @array = @_;
        my $c = grep {$_ eq 'undef'} @array;
        return $c;
}

#------------------------------------------------------------#
# Parse FASTA header line (strip ">", retain only ID part)
#------------------------------------------------------------#
sub get_header_info {
        my $line = shift;
        my @header  = split/ /, $line;
        $header[0] =~ s/\>//;
        chomp $header[0];
        return $header[0];
}

#------------------------------------------------------------#
# Get DNA reverse complement
#------------------------------------------------------------#
sub revcomp {
        my $string = shift;
        my $revcomp = reverse $string;
        $revcomp =~ tr/ATGCatgc/TACGtacg/;
        return $revcomp;
}

#------------------------------------------------------------#
# Extract features (barcode, qtag, L5) from a DNA string --
# tries both forward and reverse complement. Returns
#  (barcode, qtag, L5 label) -- or "undef" for each missing.
#------------------------------------------------------------#
sub feature_matches {
        my $line = shift;
        my $bc = 'undef';
        my $qtag = 'undef';
        my $l5 = 'undef';

        # Forward strand
        if ($line =~ /$BARCODE_MOTIF/) {
                $bc = $2;
        }
        if ($line =~ /$QTAGS/) {
                $qtag = $qtags{$1};
        }
        if ($line =~ /$L5/) {
                $l5 = $L5{$1};
        }
        # Reverse complement
        if (revcomp($line) =~ /$QTAGS/) {
                $qtag = $qtags{$1};
        }
        if (revcomp($line) =~ /$BARCODE_MOTIF/) {
                $bc = $2;
        }
        if (revcomp($line) =~ /$L5/) {
                $l5 = $L5{$1};
        }
        return ($bc, $qtag, $l5);
}
