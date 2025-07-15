use strict;
use warnings;
use Data::Dumper;
use File::Find;
use File::Basename;
use Cwd;

# Hash to count q-tags
my %qtag_counts = ();

# Motif (regex) for barcodes (16-19 bases of A/T/C/G)
my $BARCODE_MOTIF = "[ATCG]{16,19}";

# Motif (regex) for mcount; C + 3N + C + 3N + C + 3N + _constant_
my $MCOUNT_MOTIF = "[ATCG]{0,3}C[ACTG]{3}C[ACTG]{3}C[ACTG]{3}GCGCAACGCG";

# Hash for barcode counts by read set (file prefix)
my %BARCODE_COUNTS = ();

# Mapping known qtag sequences to sample numbers/labels
my %qtags = (
    "TCGGCTAGATGT" => "19",
    "AGGAACACCAAG" => "23",
    "TCGCCGAGCAGT" => "22",
    "CGAGCGCGAGGA" => "24",
    "TGGCGAATATGG" => "25",
    "TCTTCTACAACA" => "26",
    "AGCACGCCTTGT" => "27",
    "GCAACTTCTTCA" => "26_2",
    "AAGAAGTCCAAC" => "17",
    "AGGTGTCGTCAT" => '29',
);

# Log2 utility function
sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

# Calculate Hamming distance between two equal-length strings
sub hamming{
    my ($p,$o) = @_;
    return 0 if $p eq $o; # no mismatch
    my $len = length($p);
    my $num_mismatch = 0;
    for (my $i=0; $i<$len; $i++){
        ++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
    }
    return $num_mismatch;
}

# Get list of files, indexed by the first three underscore-separated name fields
my %FILES = get_files();

sub get_files {
    my %hash = ();
    my @files = glob("*.fastq.gz");
    foreach my $file (@files) {
        my @array = split/\_/, $file;
        my $string = join("_", @array[0,1,2]);
        push @{$hash{$string}}, $file;
    }
    return %hash;
}

# ----- Main per-file loop -----
foreach my $reads (keys %FILES) {
    print "Processing $reads\n";
    my %UNKNOWN = ();    # To count unknown/ambiguous qtags
    my %data = ();       # To store extracted tags
    my $data_file = $reads . "_" . 'reads.tsv';        # Output: all reads with parsed tags
    my $chimera_file = $reads . "_" . 'chimera_data.txt';  # Output: chimeric barcode-QT assignment info

    open my($of1), ">$chimera_file";
    open my($of6), ">$data_file";
    # Write header for chimera file
    print $of1 join("\t", 'run', 'total_reads', 'match_all_features_counts'), "\n";

    # Read FASTQ, convert to FASTA using seqtk; read in stream mode
    open my($fh1), "seqtk seq -A ${$FILES{$reads}}[0] | ";

    my $counter = 0;
    my %feature_counts = ();    # Types of feature patterns seen
    my %summary = ();           # Per-run summary
    my %chimera = ();           # Chimera detection: barcode => qtag counts
    my %umi = ();               # Tracks UMIs per barcode
    my %alt = ();               # Not clearly described - for possible alternate records?
    while(defined (my $line1 = <$fh1>) ) {
        next if $line1 =~ /^\>/;         # Skip fasta header
        chomp $line1;
        my @tags = ();                   # Extracted tags for this read
        my @features = ();               # Feature descriptors
        my $qtag = 'undef';              # To store sample qtag label if assigned
        $counter++;
        #---------------------------------------
        # Try to extract all motifs by regex
        if ( $line1 =~ /($MCOUNT_MOTIF)(TGCGGCCGCGAATTCCG[AG])($BARCODE_MOTIF)([GC]AATTCGATGGC)[ACTG]+GGTGTTCAAGCTT([ATCG]{12})/) {
            # If regex matches main sequence motif (all tags found)
            my @test = bin_qtags($5, keys %qtags);       # Bin by qtag with ≤2 mismatches
            if (scalar(@test) == 1) {
                # Matching only one qtag

                push @tags, $1;           # mcount region
                push @features, 'mc';     # feature type: mc
                # my @bc = split('', $3); # Unused: for extracting barcode bases
                # my $string = join('', @bc[3,4,5,7,8,9,10]);
                my $string = $3;          # extracted barcode
                push @tags, $string;      
                push @features, 'bc';     # feature type: barcode
                push @tags, $test[0];     # found (binned) qtag
                push @features, 'qtag';   # feature type: qtag
            } else {
                # No unambiguous qtag: flag this
                $UNKNOWN{$5}++;
                push @tags, 'MISSING';
                push @features, 'na';
                push @tags, 'MISSING';
                push @features, 'na';
                push @tags, 'MISSING';
                push @features, 'na';
                push @tags, 'MISSING';
                push @features, 'na';
            }
        } else {
            # Main motif not matched -- flag as MISSING
            push @tags, 'MISSING';
            push @features, 'na';
            push @tags, 'MISSING';
            push @features, 'na';
            push @tags, 'MISSING';
            push @features, 'na';
            push @tags, 'MISSING';
            push @features, 'na';
        }
        #---------------------------------------
        # If successfully identified final tag, get group label from %qtags map
        if (defined $qtags{$tags[$#tags]} ) {
            $qtag = $qtags{$tags[$#tags]};
        }

        # Update missing/total counters for report
        if ( grep { $_ eq 'MISSING'} @tags ) {
            $summary{'hits'}++;
        }
        $summary{'total'}++;
        my $string = join(",", @features);
        $feature_counts{$string}++;    # Tally feature pattern for report

        # Skip further downstream processing if tags are missing
        next if grep { $_ eq 'MISSING'} @tags;
        $BARCODE_COUNTS{$tags[1]}{$reads}++;   # Count how often barcode occurs in this run
        my $umi_bc = $tags[0] . '_' . $tags[1]; # For later, possibly
        $umi{$tags[0]}{$tags[1]}++;            # Count UMIs per barcode
        $alt{$umi_bc}{$qtag}++;                # Counts barcodes per qtag/umi combination (might be for cross-check)
        # Write result to read-data table
        print $of6 join("\t", $reads, $qtag, $tags[1], $tags[0], $line1), "\n";
        $data{$qtag}{$tags[1]}{$tags[0]}++;    # Full triple record: sample, barcode, umi
        $chimera{$tags[1]}{$qtags{$tags[2]}}++;# For each barcode, track all qtag assignments (for chimera detection)
    }

    #--------------------------
    # UMI histogram report -- seems experimental, not output by default
    my %COUNTS = ();
    foreach my $umi ( keys %umi) {
        foreach my $barcode ( keys %{$umi{$umi}} ) {
            $COUNTS{$umi{$umi}{$barcode}}++;
        }
    }

    #--------------------------
    # Chimera detection output: if a barcode is observed with more than one sample (qtag)
    foreach my $bc (keys %chimera) {
        my @q = keys %{$chimera{$bc}};
        next if scalar(@q) == '1';          # only care if barcode seen assigned to >1 qtag
        for my $item (@q) {
            print $of1 join("\t", $reads, $bc, $item, $chimera{$bc}{$item}), "\n";
        }
    }

    #--------------------------
    # Filter out barcodes that appear more than once with a UMI (not used by default, just shown as QA)
    foreach my $qtag ( keys %data) {
        foreach my $barcode (keys %{$data{$qtag}}) {
            foreach my $umi (keys %{$data{$qtag}{$barcode}} ) {
                if (defined $umi{$umi}{$barcode}) {
                    if ($umi{$umi}{$barcode} > 1) {
                        # Can delete; shown commented-out in original
                        # delete $data{$qtag}{$barcode}{$umi};
                    }
                }
            }
        }
    }
}

#-------------------------------
# Make summary data of UMI counts per qtag/barcode
sub get_mcounts {
    my %data = @_;
    my %clusters = ();
    foreach my $qtag (keys %data) {
        foreach my $barcode (keys %{$data{$qtag}}) {
            $clusters{$qtag}{$barcode} = scalar(keys %{$data{$qtag}{$barcode}});
        }
    }
    return %clusters;
}

# Utility: extract 8th field (field 7, 0-based) from each array-of-arrays, sorted numerically
sub return_list {
    my @array = @_;
    my @list = ();
    for my $item (@array) {
        push @list, $item->[7];
    }
    @list = sort {$a <=> $b} @list;
    return @list;
}

# Detect and remove barcodes within 1 mismatch of each other (potential chimeras/index hopping)
sub remove_barcodes {
    my $run = shift;
    my @barcodes = @_;
    my %cache = ();
    my $out = $run . '.chimera_one_off.txt';
    open my($df), ">$out";
    my @barcodes_filtered = ();
    for (my $i = 0; $i <= $#barcodes; $i++) {
        next if $cache{$i};
        for (my $j = $i +1; $j <= $#barcodes; $j++) {
            next if $cache{$j};
            my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
            if ($mm <= 1) {
                print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
                $cache{$j}++;
            }
        }
    }
    for my $x (0.. $#barcodes) {
        next if $cache{$x};
        push @barcodes_filtered, $barcodes[$x];
    }
    close($out);
    return @barcodes_filtered;
}

# Variant of above, not used in main script at present
sub remove_barcodes_edit {
    my $add = shift;
    my $run = shift;
    my @barcodes = @_;
    my %cache = ();
    my $out = $run . '.chimera_one_off.txt';
    open my($df), ">$out";
    my @barcodes_filtered = ();
    for (my $i = 0; $i <= $#barcodes; $i++) {
        next if $cache{$i};
        for (my $j = $i +1; $j <= $#barcodes; $j++) {
            next if $cache{$j};
            my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
            if ($mm <= 1) {
                print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
                $cache{$j}++;
            }
        }
    }
    for my $x (0.. $#barcodes) {
        next if $cache{$x};
        push @barcodes_filtered, $barcodes[$x];
    }
    close($out);
    return @barcodes_filtered;
}

# Attempt to match test sequence to any known qtag (≤2 mismatches)
sub bin_qtags {
    my $test = shift;
    my @qtags = @_;
    my %cache = ();
    my @array = ();
    for (my $i = 0; $i <= $#qtags; $i++) {
        my $mm = hamming($qtags[$i], $test);
        if ($mm <= 2) {
            push @array, $qtags[$i];
        }
    }
    return @array;
}

exit;
