#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;
use File::Find;
use File::Basename;
use Cwd;

my %qtag_counts = ();

my $BARCODE_MOTIF = "GCGGCCGCGAATTCCG[AG]([ATCG]{16,19})[GC]AATTCGATGGC";
my $MCOUNT_MOTIF = "[ATCG]{0,3}C[ACTG]{3}C[ACTG]{3}C[ACTG]{3}GCGCAACGCG";
my %BARCODE_COUNTS = ();

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
	"AGGTGTCGTCAT" => "29",
);

sub log2 {
        my $n = shift;
        return log($n)/log(2);
    }

sub hamming{
     #String length is assumed to be equal
     my ($p,$o) = @_;
     return 0 if $p eq $o; 
     my $len = length ($p);
     my $num_mismatch = 0;
     for (my $i=0; $i<$len; $i++){
     	++$num_mismatch if substr($p, $i, 1) ne substr($o, $i, 1);
     }

     return $num_mismatch;
}

my %FILES = get_files();

sub get_files {
	my %hash = ();
	my @files = glob("*.fastq.gz");
	foreach my $file (@files) {
		my @array = split/\_/, $file;
		#my $string = join("_", @array[0,1]);
		push @{$hash{$array[0]}}, $file;
	}
	return %hash;
}


foreach my $reads (keys %FILES) {
	print "Processing $reads\n";
	my %UNKNOWN = ();
	my %data = ();
	my $data_file = $reads . "_" . 'reads.tsv';
	my $chimera_file = $reads . "_" . 'chimera_data.txt';

	open my($of1), ">$chimera_file";
	open my($of6), ">$data_file";
	print $of1 join("\t", 'run', 'total_reads', 'match_all_features_counts'), "\n";
 
	open my($fh1), "seqtk seq -A ${$FILES{$reads}}[0] | ";

	my $counter = '0';
	my %feature_counts = ();
	my %summary = ();
	my %chimera = ();
	my %umi = ();
	my %alt = ();
	while(defined (my $line1 = <$fh1>) ) {
		next if $line1 =~ /^\>/;
		chomp $line1;
		my @tags = ();
		my @features = ();
		my $qtag = 'undef';
		$counter++;
		#my $BARCODE_MOTIF = "(GCGGCCGCGAATTCCG[AG]([ATCG]{16,19})([GC]AATTCGATGGC))";
		if ( $line1 =~ /($MCOUNT_MOTIF)([ATCG]+)($BARCODE_MOTIF)([ATCG]+TGGTGTTCAAGCTT)([ATCG]{12})/) {
			#print join("\t", $1, $2, $3, $4, $5, $6), "\n";
			my @test = bin_qtags($6, keys %qtags);
			#print $test[0], "\n";
			if (scalar(@test) == '1') {
				
				push @tags, $1;
				push @features, 'mc';
				#my @bc = split('', $3);
                        	#my $string = join('', @bc[3,4,5,7,8,9,10]);
                        	push @tags, $4;
                        	push @features, 'bc';
				push @tags, $test[0];
                        	push @features, 'qtag';
			}else {
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
		}else {
			push @tags, 'MISSING';
			push @features, 'na';
                    	push @tags, 'MISSING';
                        push @features, 'na';
			push @tags, 'MISSING';
                        push @features, 'na';	
                        push @tags, 'MISSING';
                        push @features, 'na';
		}
	

		if (defined $qtags{$tags[$#tags]} ) {
			$qtag = $qtags{$tags[$#tags]};
		}

		if ( grep { $_ eq 'MISSING'} @tags ) {
			$summary{'hits'}++;
		}
		$summary{'total'}++;
		my $string = join(",", @features);
		$feature_counts{$string}++;
		next if grep { $_ eq 'MISSING'} @tags;
		$BARCODE_COUNTS{$tags[1]}{$reads}++;
		my $umi_bc = $tags[0] . '_' . $tags[1];
		$umi{$tags[0]}{$tags[1]}++;
		$alt{$umi_bc}{$qtag}++;## change name of variable.
		print $of6 join("\t", $reads, $qtag, $tags[1], $tags[0], $line1, $6), "\n";
		$data{$qtag}{$tags[1]}{$tags[0]}++;
		$chimera{$tags[1]}{$qtags{$tags[2]}}++;
	}

	my %COUNTS = ();
	foreach my $umi ( keys %umi) {
		foreach my $barcode ( keys %{$umi{$umi}} ) {
			$COUNTS{$umi{$umi}{$barcode}}++;
		}
	}

	foreach my $bc (keys %chimera) {
		my @q = keys %{$chimera{$bc}};
		next if scalar(@q) == '1';
		for my $item (@q) {
       			print $of1 join("\t", $reads, $bc, $item, $chimera{$bc}{$item}), "\n";
        	}
	}
	###umi barcode filter
	foreach my $qtag ( keys %data) {
		foreach my $barcode (keys %{$data{$qtag}}) {
			foreach my $umi (keys %{$data{$qtag}{$barcode}} ) {
				if (defined $umi{$umi}{$barcode}) {
					if ($umi{$umi}{$barcode} > '1') {
					#	delete $data{$qtag}{$barcode}{$umi};
					}
				}
			}
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
	my @array = @_;
	my @list = ();
	for my $item (@array) {
		push @list, $item->[7];
	}
	@list = sort {$a <=> $b} @list;
	return @list
}
		

sub remove_barcodes {
	my $run = shift;
	my @barcodes = @_;
	my %cache = ();
	my $out = $run . '.chimera_one_off.txt';
	open my($df), ">$out";
	my @barcodes_filtered = ();
	for (my $i = 0; $i <= $#barcodes; $i++) {
		next if $cache{$i};
		#print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], '0'), "\n";
                for (my $j = $i +1; $j <= $#barcodes; $j++) {
                	next if $cache{$j};
                        my $mm = hamming($barcodes[$i][1], $barcodes[$j][1]);
			if ($mm <= '1') {
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

##############################
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
                        if ($mm <= '1') {
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




exit;
