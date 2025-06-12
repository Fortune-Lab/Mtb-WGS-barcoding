use strict;
use warnings;
use Data::Dumper;

#S       0       2227    *       .       *       *       *       P585_DMuga_A10v1_A710523_S169_espK      2227    *
#S       1       1111    *       .       *       *       *       GCF_030566105   1111    *
#H       1       2190    99.9    +       0       0       438MD673M1078D  SAMN13051719    2190    GCF_030566105   1111
#H       1       2190    99.9    +       0       0       438MD673M1078D  GCF_900520315   2190    GCF_030566105   1111

my %hash = ();
my $file = 'clusters.uc';
open my($fh), "$file";
while(my $line = <$fh> ){
	chomp $line;
	my @array = split/\t/, $line;
	next if $array[0] eq 'C';
	push @{$hash{$array[1]}{$array[0]}}, $array[8]; 
}


foreach my $cluster (sort {$a <=> $b} keys %hash) {
	if ( defined $hash{$cluster}{'H'} ) {
		print join("\t", $cluster, @{$hash{$cluster}{'S'}}, @{$hash{$cluster}{'H'}}), "\n";
	}else {
		print join("\t", $cluster, @{$hash{$cluster}{'S'}}), "\n"; 
	}
}


