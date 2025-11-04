#!/usr/bin/env perl

#use lib '/n/boslfs02/LABS/sfortune_lab/Lab/envs/barcoding/lib/perl5/site_perl/';
use lib '/n/home11/mchase/.conda/envs/barcoding/lib/perl5/site_perl';
use Math::Derivative qw(Derivative1 Derivative2);
use strict;
use warnings;
use Data::Dumper;
use File::Find;
use File::Basename;
use Cwd;
use Getopt::Long 'HelpMessage';

GetOptions(
	'min_reads=i' => \(my $MIN_READS = '10000'),
	'percent_cutoff=i' => \(my $PERCENT_CUTOFF = '1'),
	'add_chimeras=s' => \(my $DELETE = 'FALSE'),
	'filter_umi=s' => \(my $FILTER = 'FALSE'),
	'run_mode=s' => \(my $ASSAY = 'TRUE'),
	'help'     =>   sub { HelpMessage(0) },
) or HelpMessage(1);

#die "$0 requires the license holder argument (--holder)\n" unless $holder_name;
=head1 NAME

=head1 SYNOPSIS

  --percent_cutoff,-p	percent cutoff for threshold (defaults to 1).
  --min_reads,-m	minimum reads (defaults to 10000).
  --add_chimeras,-a     add chimeras and one offs to parent barcode (defaults to FALSE).
  --filter_umi,-f       filter non unique umis.
  --run_mode,-r         default assay mode. Set to FALSE and --percentcutoff to 0 for library mode.
  --help,-h		Print this help.

=head1 VERSION

0.01

=cut

my $dir = getcwd;
my $wrkdir = (split/\//, $dir)[-1]; 
sub log2 {
        my $n = shift;
        return log($n)/log(2);
    }

#sub hamming{
#     #String length is assumed to be equal
#     my ($p,$o) = @_;
#     return 0 if $p eq $o; 
#     my $len = length ($p);
#     my $num_mismatch = 0;
#     for (my $i=0; $i<$len; $i++){
#     	++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
#     }
#
#     return $num_mismatch;
#}

sub hamming {
    my ($p, $o) = @_;
    return 0 if $p eq $o;
    return undef if length($p) != length($o);  # <<<< returns nothing if length does not match
    my $len = length($p);
    my $num_mismatch = 0;
    for (my $i = 0; $i < $len; $i++) {
        ++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
    }
    return $num_mismatch;
}



my %FILES = get_files();##############

sub get_files {
	my %hash = ();
	my @files = glob("*reads.tsv");
	foreach my $file (@files) {
		my @array = split/\_/, $file;
		#my $string = join("_", @array[0,1]);
		$hash{$array[0]} = $file;
	}
	return %hash;
}
my %chimeras = get_chimeras();
my %chimera_percent = calculate_percent_chimera(%chimeras);

my $out_file = $wrkdir . '_' . 'threshold_data.csv';

open my($of), ">$out_file";

my @header = qw( run    index   qbid    counts  norm    dydx dy2dx2 dy2dx2_cutoff percent percent_chimera);
print $of join(",", @header), "\n";

foreach my $run (keys %FILES) {
	print "Processing $FILES{$run}\n";
	my %data = ();
	my %umi = ();
	open my($fh), "$FILES{$run}";
	while(my $line =<$fh> ) {
		chomp $line;
		my @array = split/\t/, $line;
		$umi{$array[3]}{$array[2]}++;
		$data{$array[1]}{$array[2]}{$array[3]}++;
	}
sub filter_umi {
	my $hash_ref = shift;
	my $umi_ref = shift;
	my $flag = shift;
      	foreach my $qtag ( keys %{$hash_ref}) {

               	foreach my $barcode (keys %{$$hash_ref{$qtag}}) {
                       	foreach my $umi (keys %{$$hash_ref{$qtag}{$barcode}} ) {
                               	if (defined ${$umi_ref}{$umi}{$barcode}) {
                                       	if (${$umi_ref}{$umi}{$barcode} > '1') {
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
	my %mcounts = get_mcounts(%data);
	#my %mcounts = %data;
	my @values = ();
	my @x = ();
	my @y = ();
	my @xt = ();
	my @yt = ();
	my $sum = '0';
	foreach my $qtag (keys %mcounts) {
		foreach my $barcodes ( keys %{$mcounts{$qtag}} ) {
			push @values, [$mcounts{$qtag}{$barcodes}, $barcodes, $qtag] unless $mcounts{$qtag}{$barcodes} == '0';
		}
	}
########################CALCULATE DERIVITIVES AND DETERMINE THRESHOLD
	@values = sort {$b->[0] <=> $a->[0] } @values;
	if ($ASSAY eq 'TRUE') {
		@values = remove_barcodes($DELETE, $run, @values);##change to reference
		@values = sort {$b->[0] <=> $a->[0] } @values;
	}

	for my $counts (@values) {
		$sum += $counts->[0];
	}
	unshift @values, $values[0];
	if ($sum <= $MIN_READS) {
		print $of $run, "\t", "FAILED THRESHOLD <= $MIN_READS" . "\n";
		next;
	}
	
	my $index = '0';
	my @percent = ();
	my $flag = '0';
	for my $item (@values) {
        	$index++;
		my $percent = '100' * ($item->[0] / $sum);
		my $rounded = sprintf("%.4f", $percent);
		push @percent, $rounded;
                push @x,  $index;
                push @y, log2($item->[0] / $sum);
		next if $flag;
		if ($percent >= $PERCENT_CUTOFF) {
			push @xt,  $index;
                	push @yt, log2($item->[0] / $sum);
		}else {
			push @xt,  $index;
                        push @yt, log2('1'/ $sum);
			$flag++;
		}
        }
	
	
		
        my @dydx=Derivative1(\@x,\@y);
        my @d2ydx2=Derivative2(\@x,\@y,);
	my @dydxt=Derivative1(\@xt,\@yt);                                                                                                                                                                          	  my @d2ydx2t=Derivative2(\@xt,\@yt,);

	my @results = ();
	
        for my $j (1..$#dydx) {
		if (defined $d2ydx2t[$j]) {
			my $qbid = $values[$j][2] . $values[$j][1];	
			my $rounded = sprintf "%.0f", $percent[$j];
			push @results, [$run, $j, $qbid, $values[$j][0], $y[$j], $dydx[$j], $d2ydx2[$j], $d2ydx2t[$j], $rounded, $values[$j][1]];
		}else {
			my $qbid = $values[$j][2] .  $values[$j][1];
			my $rounded = sprintf "%.0f", $percent[$j];
			push @results, [$run, $j, $qbid, $values[$j][0], $y[$j], $dydx[$j], $d2ydx2[$j], '0', $rounded, $values[$j][1]];
		}
        }
	my @z = return_list('6', @results);
	my @u = return_list('7', @results);

	my $flag2 = 'TRUE';
	for my $stack ( @results) {
		###add conditional for $stack->[8] == $u[0] or $stack->[7] == $z[0]
		my $chimera = '0';
		##########BROKEN need barcode
                if (defined $chimera_percent{$stack->[0]}{$stack->[9]}) {
                        $chimera =  $chimera_percent{$stack->[0]}{$stack->[9]};
                }

		if ($stack->[7] == $u[0] and $stack->[6] == $z[0]) {
			print $of join(",", @$stack[0..8], $chimera, $flag2), "\n";
			#print $of join(",", @$stack[0..8], $chimera, $flag2, $stack->[1]), "\n";
			$flag2 = 'FALSE';
		}elsif ( $stack->[7] == $u[0] or $stack->[6] == $z[0]) {
			 print $of join(",", @$stack[0..8], $chimera, $flag2, 'ALT'), "\n";
			$flag2 = 'FALSE';
		}else {
			print $of join(",", @$stack[0..8], $chimera, $flag2), "\n";
		}
	}
}


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

sub return_list {
	my $value = shift;
	my @array = @_;
	my @list = ();
	for my $item (@array) {
		##7
		push @list, $item->[$value];
	}
	@list = sort {$a <=> $b} @list;
	return @list
}
		

#sub remove_barcodes {
#	my $add = shift;
#	my $run = shift;
#	my @barcodes = @_;
#	my %cache = ();
#	my $out = $run . '.chimera_one_off.txt';
#	open my($df), ">$out";
#	my @barcodes_filtered = ();
#	for (my $i = 0; $i <= $#barcodes; $i++) {
#		next if $cache{$i};
 #               for (my $j = $i +1; $j <= $#barcodes; $j++) {
 #               	next if $cache{$j};
  #                      my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
#                        if ($mm <= '1') {
#				print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
#				$cache{$j}++;
#				next if $add eq 'FALSE';
#				$barcodes[$i][0] += $barcodes[$j][0];
#                        }
#                }
#        }
#
#	for my $x (0.. $#barcodes) {
#        	next if $cache{$x};
#                push @barcodes_filtered, $barcodes[$x];
#       	}
#	close($out);
#	return @barcodes_filtered;
#}

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
        for (my $j = $i + 1; $j <= $#barcodes; $j++) {
            next if $cache{$j};
            my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
            next unless defined $mm; # Skip if hamming() failed (unequal length)
            # Optionally, uncomment the next line for debugging:
            # unless (defined $mm) { warn "Skipping barcodes with unequal length: $barcodes[$i][1] vs $barcodes[$j][1]\n"; next; }
            if ($mm <= 1) {
                print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
                $cache{$j}++;
                next if $add eq 'FALSE';
                $barcodes[$i][0] += $barcodes[$j][0];
            }
        }
    }

    for my $x (0..$#barcodes) {
        next if $cache{$x};
        push @barcodes_filtered, $barcodes[$x];
    }
    close($df);
    return @barcodes_filtered;
}



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


#my %hash = get_chimeras();

sub get_chimeras {
        my %hash = ();
        my @files = glob("*_chimera_data.txt");
        foreach my $file (@files) {
                open my($fh), "$file";
                while ( my $line = <$fh> ) {
                        next if $. == '1';
                        chomp $line;
                        my @array = split/\t/, $line;
                        push @{$hash{$array[0]}{$array[1]}}, $array[3];
                }
                close($fh);
        }
        return %hash;
}

#calculate_percent_chimera();

sub calculate_percent_chimera {
        my %results = ();
	my %hash = @_;
        foreach my $run (keys %hash) {
                foreach my $barcode (keys %{$hash{$run}}) {

                        my @list = sort {$b <=> $a} @{$hash{$run}{$barcode}};
                        #next unless $list[0] >= '1000';
                        my $sum = '0';
                        my $rest = '0';
                        for my $i (0..$#list) {
                                $sum += $list[$i];
                                next if $i == '0';
                                $rest += $list[$i];
                        }
                        my $frac = ($rest / $sum) * 100;
                        my $rounded = sprintf("%.4f", $frac);
                        #my $d = Statistics::Diversity::Shannon->new( data => \@list);
                        #my $H = $d->index;
                        #my $E = $d->evenness;
                        $results{$run}{$barcode} = $rounded;
                        #print join("\t", $run, $barcode, $rounded, $H, $E, @list), "\n";
                }
        }
        return %results;
}



exit;
