use strict;
use warnings;
use Data::Dumper;


#220803Sha_D22-221150    CTACGCAT        CTAGCAT 1       0.0121951219512195      undef
#220803Sha_D22-221150    CTACGCAA        CTAGCAA 79      0.963414634146341       25
#220803Sha_D22-221150    CTCCGCAA        CTCGCAA 1       0.0121951219512195      undef
#220803Sha_D22-221150    CTAGGCAA        CTAGCAA 1       0.0121951219512195      undef
#220803Sha_D22-221281    GGACGGGC        GGAGGGC 1       0.0116279069767442      undef
#220803Sha_D22-221281    GTGCTGGC        GTGTGGC 20      0.232558139534884       27
#220803Sha_D22-221281    ATCCTTGA        ATCTTGA 1       0.0116279069767442      undef
#220803Sha_D22-221281    ATTCTTGA        ATTTTGA 64      0.744186046511628       25
#185_26  TTCCGAAA        31      19
#185_26  TTCCTAAA        1       undef
#185_26  GTCCGAAA        1       undef


my %hash = ();

my $file = 'Barcode_Qtag_counts_7mer.tsv';

open my($fh), $file;
while(my $line = <$fh>) {
	chomp $line;
	my @array = split/\t/, $line;
	push @{$hash{$array[0]}}, [$array[1], $array[2], $array[3]];
}
close($fh);

foreach my $run ( keys %hash) {
	my @data = sort {$b->[1] <=> $a->[1] } @{$hash{$run}};
	my @filtered = remove_barcodes(@data);
	@filtered = sort {$b->[1] <=> $a->[1] } @filtered;
	for my $result (@filtered) {
		next if $$result[1] eq '1';
		print join("\t", $run, @$result), "\n";
	}
}




sub remove_barcodes {
        my @barcodes = @_;
        my %cache = ();
        my @barcodes_filtered = ();
        for (my $i = 0; $i <= $#barcodes; $i++) {
                next if $cache{$i};
                for (my $j = $i +1; $j <= $#barcodes; $j++) {
                        next if $cache{$j};
			next unless length($barcodes[$i][0]) == length($barcodes[$j][0]);
                        my $mm = hamming($barcodes[$i][0], $barcodes[$j][0]);
                        if ($mm <= '1') {
				$barcodes[$i][1] += $barcodes[$j][1];#####
                                $cache{$j}++;
                          
                        }
                }
        }

        for my $x (0.. $#barcodes) {
                next if $cache{$x};
                push @barcodes_filtered, $barcodes[$x];
        }
        return @barcodes_filtered;
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

