#!/usr/bin/env perl

# Add a custom site_perl library path for Math::Derivative
use lib '/n/boslfs02/LABS/sfortune_lab/Lab/envs/barcoding/lib/perl5/site_perl/';
use Math::Derivative qw(Derivative1 Derivative2); # For numeric derivatives
use strict;
use warnings;
use Data::Dumper;
use File::Find;
use File::Basename;
use Cwd;
use Getopt::Long 'HelpMessage';

# Parse command-line options and set program parameters
GetOptions(
    'min_reads=i'      => \(my $MIN_READS = '10000'),    # minimum number of reads
    'percent_cutoff=i' => \(my $PERCENT_CUTOFF = '1'),   # percent cutoff for thresholding
    'add_chimeras=s'   => \(my $DELETE = 'FALSE'),       # whether to merge one-off barcodes
    'filter_umi=s'     => \(my $FILTER = 'FALSE'),       # whether to remove non-unique UMIs
    'run_mode=s'       => \(my $ASSAY = 'TRUE'),         # mode for thresholding (assay vs library)
    'help'             =>   sub { HelpMessage(0) },      # print usage/help
) or HelpMessage(1);

=head1 NAME

=head1 SYNOPSIS

  --percent_cutoff,-p   percent cutoff for threshold (defaults to 1).
  --min_reads,-m        minimum reads (defaults to 10000).
  --add_chimeras,-a     add chimeras and one offs to parent barcode (defaults to FALSE).
  --filter_umi,-f       filter non unique umis.
  --run_mode,-r         default assay mode. Set to FALSE and --percentcutoff to 0 for library mode.
  --help,-h             Print this help.

=head1 VERSION

0.01

=cut

# Get current working directory and use last path element as the 'run' name
my $dir = getcwd;
my $wrkdir = (split /\//, $dir)[-1];

# Utility: log base 2
sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

# Utility: Calculate Hamming distance (different bases count) between two equal-length strings
sub hamming {
    my ($p,$o) = @_;
    return 0 if $p eq $o; # perfect match
    my $len = length($p);
    my $num_mismatch = 0;
    for (my $i=0; $i<$len; $i++){
        next if(length($len) < $i); # prevent off-end errors
        ++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
    }
    return $num_mismatch;
}

# Get all *reads.tsv files and index by first three _-separated tokens
my %FILES = get_files();

sub get_files {
    my %hash = ();
    my @files = glob("*reads.tsv");
    foreach my $file (@files) {
        my @array = split /_/, $file;
        my $string = join("_", @array[0,1,2]);
        $hash{$string} = $file;
    }
    return %hash;
}

# Read chimera data (from prior analysis)
my %chimeras = get_chimeras();

# Calculate fraction of reads (percent) per barcode that are chimeric (multiple qtag/sample IDs)
my %chimera_percent = calculate_percent_chimera(%chimeras);

# Output summary CSV for thresholding
my $out_file = $wrkdir . '_' . 'threshold_data.csv';
open my($of), ">$out_file";

# Output CSV header, fields separated by commas
my @header = qw(run index qbid counts norm dydx dy2dx2 dy2dx2_cutoff percent percent_chimera);
print $of join(",", @header), "\n";

#-------------------- MAIN LOOP: process each run of reads --------------------
foreach my $run (keys %FILES) {
    print "Processing $FILES{$run}\n";
    my %data = (); # qtag -> barcode -> umi -> count
    my %umi = ();  # umi -> barcode -> count (to filter UMIs)

    # Read table of reads for this run
    open my($fh), "$FILES{$run}";
    while(my $line =<$fh> ) {
        chomp $line;
        my @array = split /\t/, $line;
        $umi{$array[3]}{$array[2]}++;
        $data{$array[1]}{$array[2]}{$array[3]}++;
    }

    # Optionally filter UMIs such that only unique UMI-barcode per qtag are kept
    sub filter_umi {
        my $hash_ref = shift;
        my $umi_ref = shift;
        my $flag    = shift;
        foreach my $qtag ( keys %{$hash_ref}) {
            foreach my $barcode (keys %{$$hash_ref{$qtag}}) {
                foreach my $umi (keys %{$$hash_ref{$qtag}{$barcode}} ) {
                    if (defined ${$umi_ref}{$umi}{$barcode}) {
                        if (${$umi_ref}{$umi}{$barcode} > 1) {
                            if ($flag eq 'TRUE') {
                                delete ${$hash_ref}{$qtag}{$barcode}{$umi};
                            }
                        }
                    }
                }
            }
        }
    }
    if ($FILTER eq 'TRUE') {
        filter_umi(\%data, \%umi, $FILTER);
    }

    # Count number of distinct UMIs per (qtag,barcode)
    my %mcounts = get_mcounts(%data);

    # Process all counts into an array of [count, barcode, qtag]
    my @values = ();
    my @x = ();  # x-axis: index (rank/order)
    my @y = ();  # y-axis: log2(norm_count)
    my @xt = ();
    my @yt = ();
    my $sum = 0;
    foreach my $qtag (keys %mcounts) {
        foreach my $barcodes ( keys %{$mcounts{$qtag}} ) {
            push @values, [$mcounts{$qtag}{$barcodes}, $barcodes, $qtag] unless $mcounts{$qtag}{$barcodes} == 0;
        }
    }
    # Sort by decreasing count
    @values = sort {$b->[0] <=> $a->[0] } @values;

    # Optionally collapse similar barcodes (deletes or merges) as a chimera fix
    if ($ASSAY eq 'TRUE') {
        @values = remove_barcodes($DELETE, $run, @values);
        @values = sort {$b->[0] <=> $a->[0] } @values;
    }

    # Sum total unique molecules for this run
    for my $counts (@values) {
        $sum += $counts->[0];
    }
    # Insert a copy of top value at front, for "reference"
    unshift @values, $values[0];

    # If not enough reads for this sample, skip output
    if ($sum <= $MIN_READS) {
        print $of $run, "\t", "FAILED THRESHOLD <= $MIN_READS" . "\n";
        next;
    }

    # Calculate % of total for each barcode/qtag combo
    my $index = 0;
    my @percent = ();
    my $flag = 0;
    for my $item (@values) {
        $index++;
        my $percent = 100 * ($item->[0] / $sum);
        my $rounded = sprintf("%.4f", $percent);
        push @percent, $rounded;
        push @x,  $index;
        push @y, log2($item->[0] / $sum);
        next if $flag;
        if ($percent >= $PERCENT_CUTOFF) {
            push @xt,  $index;
            push @yt, log2($item->[0] / $sum);
        } else {
            # Stop adding 'true' thresholded barcodes on first below-cutoff
            push @xt,  $index;
            push @yt, log2(1/ $sum);
            $flag++;
        }
    }

    #---- Calculate numeric derivatives for elbow/inflection-point threshold ----
    my @dydx       = Derivative1(\@x,\@y);
    my @d2ydx2     = Derivative2(\@x,\@y);
    my @dydxt      = Derivative1(\@xt,\@yt);
    my @d2ydx2t    = Derivative2(\@xt,\@yt);

    my @results = ();
    # For each barcode: store run, index, qbid, count, normalized log2-count, derivatives, percent, barcode string
    for my $j (1..$#dydx) {
        if (defined $d2ydx2t[$j]) {
            my $qbid = $values[$j][2] . $values[$j][1];
            my $rounded = sprintf "%.0f", $percent[$j];
            push @results, [$run, $j, $qbid, $values[$j][0], $y[$j], $dydx[$j], $d2ydx2[$j], $d2ydx2t[$j], $rounded, $values[$j][1]];
        } else {
            my $qbid = $values[$j][2] .  $values[$j][1];
            my $rounded = sprintf "%.0f", $percent[$j];
            push @results, [$run, $j, $qbid, $values[$j][0], $y[$j], $dydx[$j], $d2ydx2[$j], '0', $rounded, $values[$j][1]];
        }
    }

    # z and u: sorted lists of second derivative and thresholded cutoff for later
    my @z = return_list('6', @results);
    my @u = return_list('7', @results);

    # Save output, print summary for each barcode found, flagging thresholds
    my $flag2 = 'TRUE';
    for my $stack ( @results) {
        my $chimera = '0';
        # If we have per-barcode chimera percent, note it
        if (defined $chimera_percent{$stack->[0]}{$stack->[9]}) {
            $chimera =  $chimera_percent{$stack->[0]}{$stack->[9]};
        }

        # Output threshold barcode(s) - main, alternate, and 'other'
        if ($stack->[7] == $u[0] and $stack->[6] == $z[0]) {
            print $of join(",", @$stack[0..8], $chimera, $flag2), "\n";
            $flag2 = 'FALSE';
        } elsif ( $stack->[7] == $u[0] or $stack->[6] == $z[0]) {
            print $of join(",", @$stack[0..8], $chimera, $flag2, 'ALT'), "\n";
            $flag2 = 'FALSE';
        } else {
            print $of join(",", @$stack[0..8], $chimera, $flag2), "\n";
        }
    }
}

### -------------------------- SUBROUTINES --------------------------

# Returns hash of number of UMIs per qtag/barcode
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

# Returns a list of a specified field number from a list of arrayrefs,
# sorted numerically ascending
sub return_list {
    my $value = shift;
    my @array = @_;
    my @list = ();
    for my $item (@array) {
        push @list, $item->[$value];
    }
    @list = sort {$a <=> $b} @list;
    return @list
}

# Merge/remove one-off barcodes (distance 1 in Hamming)
sub remove_barcodes {
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
            if ($mm <= '1') {
                print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
                $cache{$j}++;
                next if $add eq 'FALSE'; # Only merge counts if option set
                $barcodes[$i][0] += $barcodes[$j][0];
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

# Not used in main, but provided for matching qtags by up to 2 mismatches
sub bin_qtags {
    my $test = shift;
    my @qtags = @_;
    my %cache = ();
    my @array = ();
    for (my $i = 0; $i <= $#qtags; $i++) {
        my $mm = hamming($qtags[$i], $test);
        if ($mm <= '2') {
            push @array, $qtags[$i];
        }
    }
    return @array;
}

# Parse all *_chimera_data.txt files to build hash:
# $hash{run}{barcode} = [list of qtag counts]
sub get_chimeras {
    my %hash = ();
    my @files = glob("*_chimera_data.txt");
    foreach my $file (@files) {
        open my($fh), "$file";
        while ( my $line = <$fh> ) {
            next if $. == 1; # skip header
            chomp $line;
            my @array = split /\t/, $line;
            push @{$hash{$array[0]}{$array[1]}}, $array[3]; # [run][barcode] = list of qtag-counts
        }
        close($fh);
    }
    return %hash;
}

# For each barcode in each run, compute percent of total reads not assigned to
# the predominant qtag (artifact rate)
sub calculate_percent_chimera {
    my %results = ();
    my %hash = @_;
    foreach my $run (keys %hash) {
        foreach my $barcode (keys %{$hash{$run}}) {
            my @list = sort {$b <=> $a} @{$hash{$run}{$barcode}};
            my $sum = 0;
            my $rest = 0;
            for my $i (0..$#list) {
                $sum += $list[$i];
                next if $i == 0; # skip top
                $rest += $list[$i];
            }
            my $frac = ($rest / $sum) * 100;
            my $rounded = sprintf("%.4f", $frac);
            $results{$run}{$barcode} = $rounded; # percent chimeric
        }
    }
    return %results;
}

exit;
