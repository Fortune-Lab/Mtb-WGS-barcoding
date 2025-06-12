use strict;
use warnings;
use Bio::SeqIO;

my $seq_obj = Bio::SeqIO->new(-file => 'nr.fasta', -format => 'fasta');
while(my $seq  = $seq_obj->next_seq()) {
	my @header = split/\;/, $seq->primary_id();
	$header[1] =~ s/size\=//;
	print join("\t", @header, $seq->seq()), "\n";
}
