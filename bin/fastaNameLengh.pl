#!/usr/bin/perl -w

## Help text in perldoc format
=head1 SYNOPSIS

fastaNamesSizes.pl - returns a tab list with name and size of sequences 
                     in a fasta file on the standard output, and some stats
                     on the standard error.

=head1 USAGE

fastaNamesSizes.pl [-f format] fasta_file

=head1 AUTHOR

Lionel Guy (lionel.guy@ebc.uu.se)

=cut

## load perl packages
use strict;
use Getopt::Std;
use Bio::SeqIO;

## handle input parameters
our $opt_f = 'fasta';
getopts('f:');

usage() unless @ARGV;

## instatiate variables
my ($count, $sum) = (0, 0);
my ($min, $max);
my $file = shift @ARGV;

## use Bio::SeqIO package to read input file
my $seq_io = Bio::SeqIO->new(-file => "$file", -format => $opt_f);

## Iterate over sequences in input file to 
## generate idividual sequence and summary 
## statistics.
while (my $seq = $seq_io->next_seq){ ## selects one sequence at a time
    ## set variables for THIS sequence
    my $id = $seq->display_id;
    my $len = $seq->length;
    ## print output for THIS sequence 
    print "$id\t$len\n";
    
    ## set stats for ALL sequences
    $count++; ## seq count
    $sum += $len; ## add length of this sequence to all other sequences
    $min = $len if (!$min || $len < $min); ## find shortest sequence
    $max = $len if (!$max || $len > $max); ## find longest sequence
}
my $average = int($sum/$count+0.5); ## Compute average

## print summary statistics to standard error
print STDERR "# $count sequences, total length $sum.\n";
print STDERR "# Minimum len: $min. Max: $max. Average: $average\n";

## subroutines
## Display helptext and exit program
sub usage{
    system("perldoc $0");
    exit;
}
