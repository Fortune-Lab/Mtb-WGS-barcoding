#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# CT C ATT C ATA C TTA GCGCAACGCGT 
#CTG C TGT C AGA C GGT GCGCAACGCGT

my $LINK = "GCGCAACGCGTGCGGCCGCGAATTCCGA";
#my $LINK = "GCGCAACGCGTG[ATCG]{16}";
my $MCOUNT_MOTIF = "C[ACTG]{3}C[ACTG]{3}C[ACTG]{3}";
my %BARCODE_COUNTS = ();

my %FILES = get_files();

#my %hash = get_barcodes();

open my($of), ">", "log_file";

sub get_barcodes {
	my %hash = ();
	my $file = "raw_hits_edited_curated.tsv";
	open my($fh), $file;
	while (my $line =<$fh> ) {
		chomp $line;
		my @array = split/\t/, $line;
		$hash{$array[2]} = $array[0];
	}
	return %hash;
}


sub get_files {
	my %hash = ();
	#_R1_001.fastq.gz
	my @files = glob("*_R1_001.fastq.gz");
	foreach my $file (@files) {
		next if $file =~ /Undetermined/;
		my @array = split/\_/, $file;
		#B6_2_14_S14_R1_001.fastq.gz
		my $string = join("_", @array[0,1,2]);
		$hash{$string} = $file;
	}
	return %hash;
}
#@LH00697:42:22WT2HLT3:1:1107:2433:1048 1:N:0:ATTACTCG+AGGCTATA
foreach my $reads (keys %FILES) {
	#print "Processing $reads\n";
	my $total = '0';
	my $counter = '0';
	open my($fh1), "seqtk seq -A $FILES{$reads} | ";
	while(defined (my $line1 = <$fh1>) ) {
		next if $line1 =~ /^\>/;
		chomp $line1;
		$total++;
		#if ( $line1 =~ /($MCOUNT_MOTIF)($LINK)([ATCG]{17})/) {
		#if ( $line1 =~ /($MCOUNT_MOTIF)($LINK)([ATCG]+$)/) {
		if ( $line1 =~ /(GCGGCCGCGAATTCCG[AG])([ATCG]{16,19})([GC]AATTCGATGGCCT)/) {
		#if ( $line1 =~ /($MCOUNT_MOTIF).+(GCGGCCGCGAATTCCG[AG])([ATCG]+)([GC]AATTCGATGGCCT)/) {
			my $link1 = $1;
			my $bc = $2;
			my $link2 = $3;
			#if (defined $hash{$bc}) {
				$counter++;
			#	print join("\t", $reads, $mc, $link1, $bc, $link2, $hash{$bc}), "\n"; 		
			#}else {
				print join("\t", $reads, $link1, $bc, $link2), "\n";	
			#}
		}
	}
	print $of join("\t", $reads, $total, $counter), "\n";
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
			#print join("\t", @{$barcodes[$i]}). "\n";	
				print $df join("\t", $run, $barcodes[$i][1], $barcodes[$i][0], $barcodes[$j][1], $barcodes[$j][0], $mm), "\n";
				#$barcodes[$i][0] += $barcodes[$j][0];
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
