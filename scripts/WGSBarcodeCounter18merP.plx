use strict;
use warnings;

my %FILES = get_files(shift @ARGV);

my %L5 = (
        "CCTAGTCTTCCAAACTAGCTACGCGGGTTCGATTCCCGTCGCCCGCTCCGCTGGTCAGA" => 'insertion_site',
        "GAGCCCTAGTCTTCCAAACTAGCGACGCGGGTTCGATTCCCGTCGCCCGCTCGGGCCATG" => 'minus_insertion',
);


my $L5 = "(CCTAGTCTTCCAAACTAGCTACGCGGGTTCGATTCCCGTCGCCCGCTCCGCTGGTCAGA|GAGCCCTAGTCTTCCAAACTAGCGACGCGGGTTCGATTCCCGTCGCCCGCTCGGGCCATG)";


#my $BARCODE_MOTIF = "(GCGGCCGCGAATTCCG[AG]([ATCG]{16,19})([GC]AATTCGATGGC))";
my $BARCODE_MOTIF = "(GCGGCCGCGAATTCCG[AG]([ATCG]+)([GC]AATTCGATGGC))";
my $QTAGS= "(TGGTGTTCAAGCTTTCGGCTAGATGT|TGGTGTTCAAGCTTAGGAACACCAAG|TGGTGTTCAAGCTTTCGCCGAGCAGT|TGGTGTTCAAGCTTCGAGCGCGAGGA|TGGTGTTCAAGCTTTGGCGAATATGG|TGGTGTTCAAGCTTTCTTCTACAACA|TGGTGTTCAAGCTTAGCACGCCTTGT|TGGTGTTCAAGCTTGCAACTTCTTCA|TGGTGTTCAAGCTTAAGAAGTCCAAC|TGGTGTTCAAGCTTAGGTGTCGTCAT)";
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
	"TGGTGTTCAAGCTTAGGTGTCGTCAT" => "29",
);

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

sub get_files {
        my %hash = ();
	my $string = shift;
	my @array = split/\t/, $string;
	foreach my $seq (@array[1,2]) {
                push @{$hash{$array[0]}}, $seq;
        }
        return %hash;
}



foreach my $run ( keys %FILES){ 

	open my ($fh1), "seqtk seq -A $FILES{$run}[0] |" or die $!;
	open my ($fh2), "seqtk seq -A $FILES{$run}[1] |" or die $!;
	my $header = "undef";
	while (!eof($fh1) and !eof($fh2)) {
		my $read1 = <$fh1>;
		my $read2 = <$fh2>;
		if ($read1 =~ /^\>/) {
			$header = get_header_info($read1);
			next;
		}
		my @read1 = feature_matches($read1);
		my @read2 = feature_matches($read2);
		next if (count_item(@read1) == '3' and count_item(@read2) == '3');
		if (my $jnc = merge_junction($read1[2],$read2[2])) {
			print join("\t", $run, $header, $jnc), "\n";
			next;
		}
                my %hash = combine_reads($read1[0], $read1[1], $read2[0], $read2[1]);
		my $mm = "FALSE";
		if (scalar(keys %hash) == '2') {
			$mm = "TRUE";
		}	
		foreach my $bc (keys %hash) {
			if (scalar(keys %{$hash{$bc}}) == '2') {
				delete $hash{$bc}{'undef'};
			}
			foreach my $qtag (keys %{$hash{$bc}} ) {
				print join("\t", $run, $header, $bc, $qtag, $mm), "\n";
			}
		}	
	}
}


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



sub count_item {
	my @array = @_;
	my $c = grep {$_ eq 'undef'} @array;
	return $c;
}

sub get_header_info {
	my $line = shift;
	my @header  = split/ /, $line;
        #my $r = (split/\:/, $header[1])[0];
        #$header = $header[0];
        $header[0] =~ s/\>//;
	chomp $header[0];
	return $header[0];

}

sub revcomp {
        my $string = shift;
        my $revcomp = reverse $string;
        $revcomp =~ tr/ATGCatgc/TACGtacg/;
        return $revcomp;
}

sub feature_matches {
	my $line = shift;
	my $bc = 'undef';
	my $qtag = 'undef';
	my $l5 = 'undef';
       	if ($line =~ /$BARCODE_MOTIF/) {
        	$bc = $2;
    	}
        if ($line =~ /$QTAGS/) {
             	$qtag = $qtags{$1};
      	}
	if ($line =~ /$L5/) {
		$l5 = $L5{$1};
       	}
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

