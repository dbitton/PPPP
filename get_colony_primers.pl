#!/usr/local/bin/perl -w

#
# Author: cjp
# Maintainer: cjp
# Created:
# Last Modified: 2006-08-29 cjp: initial version
#

use strict;

use SangerPaths qw(core);
use SangerWeb;

use CGI qw/:cgi :standard/;
use CGI::Carp qw(fatalsToBrowser);

#use Data::Dumper;

use SDBM_File;
use Fcntl;

use Website::Utilities::IdGenerator;

use vars qw(%genehash $datadir $TMPDIR);

$datadir = '/nfs/WWWdev/SANGER_docs/data/PostGenomics/S_pombe/PPPP/';
$TMPDIR = '/nfs/WWWdev/SANGER_docs/htdocs/tmp/spge';

$ENV{PATH} = '';

&main();
1;

sub file_to_primer {

    my ($out_file, $gene, $gene_seq, $inc, $opt_len, $min_len, $max_len,
        $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $length) = @_;

    my $excl_start = 1000 - 3*$inc;
    my $excl_len   = $length + 1 + 6*$inc;

    my $max_amp_len = length $gene_seq;

    if ($opt_len < $min_len) {
        $opt_len = $min_len;
    }

    open PRIMER, ">$out_file"
        or die "Couldn't write $out_file: $!";

    print PRIMER
        "PRIMER_SEQUENCE_ID=$gene\n",
        "SEQUENCE=$gene_seq\n",
        "PRIMER_OPT_SIZE=$opt_len\n",
        "PRIMER_MIN_SIZE=$min_len\n",
        "PRIMER_MAX_SIZE=$max_len\n",
        "PRIMER_OPT_TM=$opt_tm\n",
        "PRIMER_MIN_TM=$min_tm\n",
        "PRIMER_MAX_TM=$max_tm\n",
        "PRIMER_MIN_GC=$min_gc\n",
        "PRIMER_MAX_GC=$max_gc\n",
        "PRIMER_PRODUCT_SIZE_RANGE=$excl_len-$max_amp_len\n",
        "EXCLUDED_REGION=$excl_start,$excl_len\n",
        "PRIMER_SALT_CONC=50.0\n",
        "PRIMER_DNA_CONC=50.0\n",
        "PRIMER_SELF_ANY=8\n",
        "PRIMER_SELF_END=3\n",
        "PRIMER_GC_CLAMP=0\n",
        "PRIMER_NUM_RETURN=10\n",
        "=\n";

    close PRIMER
        or die "Couldn't close for writing $out_file: $!";

}

sub get_seq {

    my ($primer3, $gene, $inc, $add_seq, $opt_len, $min_len, $max_len,
        $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc) = @_;

    tie (%genehash, 'SDBM_File', $datadir . 'pombe_data_cjp.dbr', O_RDONLY, 0666);

    my ($chr, $complement, $start, $end) = (-1, -1, -1, -1);

    foreach my $key (keys %genehash) {

        if (lc $key eq lc $gene) {

            ($chr, $complement, $start, $end) = split /,/, $genehash{$key};
        }

    }

    my $length = ($end - $start) + 1;

    my $chr_seq = $datadir . 'all_chromosomes.seq';

    open CHR_SEQ, $chr_seq
        or die "Cannot open file $chr_seq: $!";


    my $seq = '';
    my $gene_seq = '';

    while (defined (my $line = <CHR_SEQ>) ) {

        if (substr($line, 0, 5) eq ">CHR${chr}") {

            $seq = <CHR_SEQ>;

        }

        $gene_seq = substr($seq, $start - 1 - $add_seq, $length + 2 * $add_seq);

        if ($complement == 1) {

            $gene_seq =~ tr/ACGT/TGCA/;
            $gene_seq = reverse $gene_seq;

        }

    }

    #print "[$gene_seq]\n";

    close CHR_SEQ
        or die "Cannot open file $chr_seq: $!";

    my $out_file     = $TMPDIR . '/primer.tmp';
    my $results_file = $TMPDIR . '/results.tmp';

    file_to_primer($out_file, $gene, $gene_seq, $inc, $opt_len, $min_len, $max_len,
                   $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $length);

    my $primer3 = $datadir . 'primer3_core';

    system("$primer3 < $out_file > $results_file");

    open PRIMER_OUT, $results_file
        or die "Cannot open file $results_file: $!";

    while (defined (my $line = <PRIMER_OUT>) ) {

        next
            unless $line =~ /^PRIMER_(LEFT|RIGHT)(_\d+)?_?(SEQUENCE|=|TM|GC_PERCENT)/;

        $line =~ s/^PRIMER_(LEFT|RIGHT)_?(\d+)?_(SEQUENCE).*=/(ucfirst(lc($1))) . " Primer: " . ($1 eq "LEFT"?" ":"") . ($2?$2+1:1) . "\t" . (ucfirst(lc($3))) . ": \t\t"/e;

        $line =~ s/^PRIMER_(LEFT|RIGHT)_?(\d+)?_(TM|GC).*=([\d\.]+)/(ucfirst(lc($1))) . " Primer: " . ($1 eq "LEFT"?" ":"") . ($2?$2+1:1) . "\t" . (uc($3)) . ": \t\t" . (sprintf "%.1f", $4)/e;

        $line =~ s/^PRIMER_(LEFT)_?(\d+)?=(\d+),(\d+)/"Left Primer:  " . ($2?$2+1:1) . "\tDistance from ORF: " . (1000-$3) . "\t\tPrimer length: " . $4/e;
        $line =~ s/^PRIMER_(RIGHT)_?(\d+)?=(\d+),(\d+)/"Right Primer: " . ($2?$2+1:1) . "\tDistance from ORF: " . ($3-1000-$length) . "\t\tPrimer length: " . $4/e;

        print "<pre>$line</pre>";

        print "<pre>\n\n</pre>"
            if $line =~ /^Right.+GC/i;

    }

    close PRIMER_OUT
        or die "Cannot open file $results_file: $!";


}

sub main {

    my $title = qq(Get Pombe Primers For Colony Checking of Genetic Modifications);
    my $cgi   = CGI->new();
    my $sw    = SangerWeb->new({'title'   => $title,
                                'banner'  => $title,
                                'inifile' => qq($ENV{'DOCUMENT_ROOT'}/PostGenomics/S_pombe/header.ini),
                               });

    #my $primer3 = '/usr/local/pubseq/bin/primer3_core';

    my $primer3 = $sw->server_root() . "/bin-offline/primer3";

    print $sw->header();

    print
        start_form(),
        p(qq(Gene Name @{[textfield('name')]})),
        p(qq(Original primer increment @{[textfield('inc', '40', 3)]})),
        p(qq(Amount of up- or downstream sequence to use @{[textfield('add_seq', '1000', 3)]})),
        p(qq(Optimum primer length @{[textfield('opt_len', '20', 3)]})),
        p(qq(Minimum primer length @{[textfield('min_len', '20', 3)]})),
        p(qq(Maximum primer length @{[textfield('max_len', '28', 3)]})),
        p(qq(Optimum Tm @{[textfield('opt_tm', '60.0', 3)]})),
        p(qq(Minimum Tm @{[textfield('min_tm', '57.0', 3)]})),
        p(qq(Maximum Tm @{[textfield('max_tm', '63.0', 3)]})),
        p(qq(Minimum GC content @{[textfield('min_gc', '30', 3)]})),
        p(qq(Maximum GC content @{[textfield('max_gc', '70', 3)]})),
        submit(),
        reset(),
        end_form(),
        hr();

    if ($cgi->param()) {

        my $gene    = $cgi->param('name') || "";
        ($gene)     = $gene =~ /([a-zA-Z0-9_\-\.]+)/;

        my $inc     = $cgi->param('inc') || "";
        ($inc)      = $inc =~ /(\d+)/;

        my $add_seq = $cgi->param('add_seq') || "";
        ($add_seq)  = $add_seq =~ /(\d+)/;

        my $opt_len = $cgi->param('opt_len') || "";
        ($opt_len)  = $opt_len =~ /(\d+)/;

        my $min_len = $cgi->param('min_len') || "";
        ($min_len)  = $min_len =~ /(\d+)/;

        my $max_len = $cgi->param('max_len') || "";
        ($max_len)  = $max_len =~ /(\d+)/;

        my $opt_tm  = $cgi->param('opt_tm') || "";
        ($opt_tm)   = $opt_tm =~ /([\d\.]+)/;

        my $min_tm  = $cgi->param('min_tm') || "";
        ($min_tm)   = $min_tm =~ /([\d\.]+)/;

        my $max_tm  = $cgi->param('max_tm') || "";
        ($max_tm)   = $max_tm =~ /([\d\.]+)/;

        my $min_gc  = $cgi->param('min_gc') || "";
        ($min_gc)   = $min_gc =~ /([\d\.]+)/;

        my $max_gc  = $cgi->param('max_gc') || "";
        ($max_gc)   = $max_gc =~ /([\d\.]+)/;

        my $seq =
            get_seq($primer3, $gene, $inc, $add_seq, $opt_len, $min_len, $max_len,
                    $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc);

        #print "[$cmd]\n"; debugging print statement

        #system($cmd);

    }

    print $sw->footer();

}
