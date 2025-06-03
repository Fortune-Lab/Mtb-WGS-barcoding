use strict;
use warnings;
use Data::Dumper;

my %link = ();
my %mm = ();
my %bc = ();
my %l5 = ();

my $file = shift @ARGV;
open my($fh), $file;

while(my $line = <$fh>) {
	chomp $line;
	my @array = split/\t/, $line;
	@array = @array[0,1,4,5,6];
	if ($array[2] =~ /insertion/) {
		$l5{$array[0]}{$array[2]}++;
	}elsif ($array[2] =~ /[ATCG]+/ and $array[3] eq 'undef' and $array[4] eq 'FALSE') {
		$bc{$array[0]}{$array[2]}++;
	}elsif ($array[2] =~ /[ATCG]+/ and $array[4] eq 'TRUE') {
		push @{$mm{$array[0]}{$array[1]}}, [@array[2..$#array]];
	}elsif ($array[2] =~ /[ATCG]+/ and $array[3] =~ /\d+/) {
		push @{$link{$array[0]}{$array[2]}}, $array[3];
		$bc{$array[0]}{$array[2]}++;
	}

}


print Dumper(\%bc);
	
#fix mismatched barcodes
#add logic to catch qtag if present
foreach my $run ( keys %mm) {
	my $qtag = 'undef';
	my $bc = "";
	foreach my $read (keys %{$mm{$run}} ) {
		my $score1 = '0';
		my $score2 = '0';
		my $qtag = 'undef';
		my $bc = "";
		if (${$mm{$run}{$read}[0]}[1] =~ /\d+/ ) {
			$qtag = ${$mm{$run}{$read}[0]}[1];
		}
		if (${$mm{$run}{$read}[1]}[1] =~ /\d+/ ) {
			$qtag = ${$mm{$run}{$read}[1]}[1];
                }
		 
		if (defined $bc{$run}{${$mm{$run}{$read}[0]}[0]} ) {
			$score1 = $bc{$run}{${$mm{$run}{$read}[0]}[0]};
		}
		if (defined $bc{$run}{${$mm{$run}{$read}[1]}[0]} ) {
                        $score2 = $bc{$run}{${$mm{$run}{$read}[1]}[0]};
                }
		if ($score1 > $score2 ) {
			$bc{$run}{${$mm{$run}{$read}[0]}[0]} += '1';
			$bc = ${$mm{$run}{$read}[0]}[0];
		} else { 
			$bc{$run}{${$mm{$run}{$read}[1]}[0]} += '1';
			$bc = ${$mm{$run}{$read}[1]}[0];
		}
	}
	if ($qtag =~ /\d+/) {
		push @{$link{$run}{$bc}}, $qtag;
	}
}


open my($of), ">Barcode_Qtag_18_counts.tsv";

foreach my $run (keys %bc) {
	foreach my $bc (keys %{$bc{$run}}) {
		if (defined $link{$run}{$bc}) {
			print $of join("\t", $run, $bc, $bc{$run}{$bc}, reduce_array(@{$link{$run}{$bc}})), "\n";
		}else {
			print $of join("\t", $run, $bc, $bc{$run}{$bc}, 'undef'), "\n";
		}
	}
}

close($of);

open my($if), ">L5_integation_18_data.tsv";

print $if join("\t", 'sample', 'minus', 'plus', 'fraction'), "\n";

foreach my $run ( keys %l5) {
	my @counts = ();
	for my $junction ('minus_insertion', 'insertion_site') {
		if (defined $l5{$run}{$junction}) {
			push @counts, $l5{$run}{$junction};
		}else {
			 push @counts, '0';
		}
	}
	my $prct = $counts[1] / ($counts[0] + $counts[1]);
	print $if join("\t", $run, @counts, $prct), "\n";		
}

close($if);

sub reduce_array {
	my @array = @_;
	my %hash = ();
	foreach my $item (@array) {
		$hash{$item}++;
	}
	return keys %hash;
}


