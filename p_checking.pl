#!/usr/local/bin/perl -wT
#####################################################################
#
# Author: cjp
# Maintainer: cjp
# Created: 2006-08-29
# Last Modified: 2006-08-29 cjp: initial version
#
#####################################################################


use lib qw(/var/www/lib/core);
use SangerWeb;

use CGI qw/:cgi :standard/;
use CGI::Carp qw(fatalsToBrowser);

#use Data::Dumper;

use SDBM_File;
use Fcntl;

use Website::Utilities::IdGenerator;

use vars qw(%genehash $datadir $bindir $TMPURI $TMPDIR $act_gene_seq);

$act_gene_seq = '';

$datadir = SangerWeb->data_root() . '/PPPP/';
$datadir = SangerWeb->data_root() . '/PPPP/';

#####################
use strict;

#use SangerPaths qw(core);
#use SangerWeb;

#use CGI qw/:cgi :standard/;
#use CGI::Carp qw(fatalsToBrowser);

#use Data::Dumper;

#use SDBM_File;
#use Fcntl;

#use Website::Utilities::IdGenerator;

#use vars qw(%genehash $datadir $bindir $TMPURI $TMPDIR $act_gene_seq);

#$act_gene_seq = '';

#$datadir = SangerWeb->data_root() . '/PostGenomics/S_pombe/PPPP/';
#$datadir = SangerWeb->data_root() . '/PostGenomics/S_pombe/PPPP/';

$TMPURI = '/tmp/spge';
$TMPDIR = $ENV{DOCUMENT_ROOT} . $TMPURI;

$ENV{PATH} = '';

&main();
1;

sub rev_comp {

  my $primer = shift;

  $primer = reverse $primer;
  $primer =~ tr/ACGT/TGCA/;

  return $primer;

}


sub file_to_primer {

    my ($out_file, $gene, $add_seq, $design_seq, $len, $inc, $opt_len, $min_len, $max_len,
        $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $gene_length, $for_sel, $rev_sel) = @_;

    #print "[$for_sel] [$rev_sel] [$inc]\n";

    my ($for_val) = $for_sel =~ /^([F\d]+)_[ACGT]+_[ACGT]*$/;
    my ($rev_val) = $rev_sel =~ /^([R\d]+)_[ACGT]+_[ACGT]*$/;

    if (not $for_val) {

        $for_val = 0;

    }

    if (not $rev_val) {

        $rev_val = 0;

    }

    my $excl_start = 0;
    my $excl_len   = 0;

    my $max_amp_len = length $design_seq;

    if ($for_val eq 'F') {

        $excl_start = $add_seq - ($rev_val*$inc + $len + 1);

    } elsif ($rev_val eq 'R') {

        $excl_start = $add_seq - ($for_val*$inc + $len + 1);

    } else {

        $excl_start = $add_seq - ($for_val*$inc + $len + 1);

    }

    if ($for_val eq 'F') {

        $excl_len = 2*$len + $rev_val*$inc + 5;

    } elsif ($rev_val eq 'R') {

        $excl_len = 2*$len + $for_val*$inc + 5;

    } else {

        $excl_len = $gene_length + 2*$len + ($for_val + $rev_val)*$inc + 2;

    }

    if ($opt_len < $min_len) {

        $opt_len = $min_len;

    }

    my $min_amp_len = $excl_len;

    if ($min_amp_len < $excl_start) {

        $min_amp_len = $excl_start + 1;

    } elsif ($min_amp_len < ($max_amp_len - ($excl_start + $excl_len) ) ) {

        $min_amp_len =$max_amp_len - ($excl_start + $excl_len) + 1;

    }

    open PRIMER, ">$out_file"
        or die "Couldn't write to $out_file: $!";

    print PRIMER
        "PRIMER_SEQUENCE_ID=$gene\n",
        "SEQUENCE=$design_seq\n",
        "PRIMER_OPT_SIZE=$opt_len\n",
        "PRIMER_MIN_SIZE=$min_len\n",
        "PRIMER_MAX_SIZE=$max_len\n",
        "PRIMER_OPT_TM=$opt_tm\n",
        "PRIMER_MIN_TM=$min_tm\n",
        "PRIMER_MAX_TM=$max_tm\n",
        "PRIMER_MIN_GC=$min_gc\n",
        "PRIMER_MAX_GC=$max_gc\n",
        "PRIMER_PRODUCT_SIZE_RANGE=$min_amp_len-$max_amp_len\n",
        "EXCLUDED_REGION=$excl_start,$excl_len\n",
        "PRIMER_SALT_CONC=50.0\n",
        "PRIMER_DNA_CONC=50.0\n",
        "PRIMER_SELF_ANY=8\n",
        "PRIMER_SELF_END=3\n",
        "PRIMER_GC_CLAMP=0\n",
        "PRIMER_NUM_RETURN=1\n",
        "=\n";

    close PRIMER
        or die "Couldn't close for writing $out_file: $!";

}

sub get_seq {

    my ($primer3, $gene, $len, $inc, $add_seq, $opt_len, $min_len, $max_len,
        $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $for_sel, $rev_sel) = @_;

    my $dbm_file = $datadir . 'pombe_data_cjp.dbr';

    tie %genehash, 'SDBM_File', $dbm_file, O_RDONLY, 0666
        or die "Cannot tie DMB file $dbm_file: $!";

    my ($chr, $complement, $start, $end) = (-1, -1, -1, -1);

    foreach my $key (keys %genehash) {

        if (lc $key eq lc $gene) {

            ($chr, $complement, $start, $end) = split /,/, $genehash{$key};

        }

    }

    untie %genehash
        or die "Cannot untie DMB file $dbm_file: $!";

    my ($for_val, $for_seq, $for_tag) = $for_sel =~ /^([F\d]+)_([ACGT])+_([ACGT]*)$/;
    my ($rev_val, $rev_seq, $rev_tag) = $rev_sel =~ /^([R\d]+)_([ACGT])+_([ACGT]*)$/;

    my $gene_length = ($end - $start) + 1;

    my $chr_file = $datadir . 'all_chromosomes.seq';

    open CHR_FILE, $chr_file
        or die "Cannot open file $chr_file: $!";


    my $chr_seq    = '';
    my $gene_seq   = '';
    my $design_seq = '';

    while (defined (my $line = <CHR_FILE>) ) {

        if (substr($line, 0, 5) eq ">CHR${chr}") {

            $chr_seq = <CHR_FILE>;

            $chr_seq = uc $chr_seq;

        }

        $gene_seq = substr $chr_seq, $start - ($add_seq + 1), $gene_length + 2 * $add_seq;

        if ( ($for_val eq 'F' and not $complement) or
             ($rev_val eq 'R' and $complement) ) {

            $design_seq =
                substr
                    $chr_seq,
                    $start + $gene_length - ($add_seq + 1),
                    2 * $add_seq;

        } elsif ( ($rev_val eq 'R' and not $complement) or
                  ($for_val eq 'F' and $complement) ) {

            $design_seq =
                substr
                    $chr_seq,
                    $start - ($add_seq + 1),
                    2 * $add_seq;

        } else {

            $design_seq = $gene_seq;

        }

        $act_gene_seq = substr $chr_seq, $start - 1, $gene_length;

        if ($complement) {

            $gene_seq = rev_comp($gene_seq);

            $act_gene_seq = rev_comp($act_gene_seq);

            $design_seq = rev_comp($design_seq);

        }

        #print "length: ", length $gene_seq, "<br />\n";

    }

    close CHR_FILE
        or die "Cannot open file $chr_file: $!";

    my ($fpath) = $TMPDIR =~ m|([a-zA-Z0-9_\-/\.]+)|;

    my $id = Website::Utilities::IdGenerator->get_unique_id();

    my $out_file     = $fpath . "/${id}_primer.tmp";
    my $results_file = $fpath . "/${id}_results.tmp";

    #print "[$out_file] [$results_file]\n";

    file_to_primer($out_file, $gene, $add_seq, $design_seq, $len, $inc, $opt_len, $min_len, $max_len,
                   $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $gene_length, $for_sel, $rev_sel);

    system("$primer3 < $out_file > $results_file");

    open PRIMER_OUT, $results_file
        or die "Cannot open file $results_file: $!";

    while (defined (my $line = <PRIMER_OUT>) ) {

        next
            unless $line =~ /^PRIMER_(LEFT|RIGHT)(_\d+)?_?(SEQUENCE|=|TM|GC_PERCENT)/;

        $line =~ s/^PRIMER_(LEFT|RIGHT)_?(\d+)?_(SEQUENCE).*=/(ucfirst(lc($1))) . " Primer: " . ($1 eq "LEFT"?" ":"") . "\t" . (ucfirst(lc($3))) . ": \t"/e;

        $line =~ s/^PRIMER_(LEFT|RIGHT)_?(\d+)?_(TM|GC).*=([\d\.]+)/(ucfirst(lc($1))) . " Primer: " . ($1 eq "LEFT"?" ":"") . "\t" . (uc($3)) . ": \t\t" . (sprintf "%.1f", $4)/e;

        # ***** TODO *****

        $line =~ s/^PRIMER_(LEFT)_?(\d+)?=(\d+),(\d+)/"Left Primer:  \tLength: \t" . $4/e;
        $line =~ s/^PRIMER_(RIGHT)_?(\d+)?=(\d+),(\d+)/"Right Primer: \tLength: \t" . $4/e;

        print "<pre>$line</pre>";

        print "<pre>\n</pre><br />"
            if $line =~ /^Right.+GC/i;

    }

    close PRIMER_OUT
        or die "Cannot open file $results_file: $!";


    $for_sel =~ s/^[F\d]+_([ACGT]+)_[ACGT]*$/$1/;
    $rev_sel =~ s/^[R\d]+_([ACGT]+)_[ACGT]*$/$1/;

    $rev_sel = rev_comp($rev_sel);

    $gene_seq =~ s/$act_gene_seq/<font color='yellow'>$act_gene_seq<\/font>/;

    $gene_seq =~ s/$for_sel/<font color='red'>$for_sel<\/font>/;
    $gene_seq =~ s/$rev_sel/<font color='blue'>$rev_sel<\/font>/;

    my $print_flag = 0;

    my $error_flag = 0;

    open PRIMER_OUT, $results_file
        or die "Cannot open file $results_file: $!";

    while (defined (my $line = <PRIMER_OUT>) ) {

        $error_flag = 1, last
            if $line =~ /^PRIMER_ERROR/;

        if ($line =~ /^PRIMER_LEFT_(\d+_)?SEQUENCE/) {

            my ($left_seq) = $line =~ /^PRIMER_LEFT_(?:\d+_)?SEQUENCE.*=([ACGT]+)/;

            $gene_seq =~ s/$left_seq/<font color='green'>$left_seq<\/font>/s;

            $print_flag = 1;

        }

        if ($line =~ /^PRIMER_RIGHT_(\d+_)?SEQUENCE/) {

            my ($right_seq) = $line =~ /^PRIMER_RIGHT_(?:\d+_)?SEQUENCE.*=([ACGT]+)/;

            $right_seq = rev_comp($right_seq);

            $gene_seq =~ s/$right_seq/<font color='green'>$right_seq<\/font>/s;

        }

    }

    $gene_seq =~ s/([ACGT]{1,80})/$1\n/g;
    $gene_seq =~ s/([ACGT]{1,10})/$1 /g;

    print
        "<em><font color='red'>NOTE:</font> Use the sequence above to order the reverse (right) checking ",
        "primer; the diagram below is just to show context.  <br />A possible choice of primer to go ",
        "with the left checking primer (in the kanamycin plasmid) is: CGGATGTGATGTGAGAACTGTATCCTAGC ",
        "(Kan reverse) and one for the right checking primer is: CGCTATACTGCTGTCGATTCG (Kan forward).",
#        "<br />Also, the above targeting primers do not contain any extra plasmid or tag sequences and ",
#        "shouldn't be used to order primers.  ",
        "</em><br /><br />Gene sequences are shown with the ORF in yellow, the forward gene deletion primer ",
        "in red, the reverse primer in blue (reverse complemented), and the checking primers in green ",
        "(with the reverse primer also being reverse complemented):<br />\n<pre>$gene_seq</pre>"
            if $print_flag;

    print
        "No primers found; expand your amount of up- and downstream sequence to use ",
        "and/or loosen your PCR parameters.<br />"
            if $error_flag;

    close PRIMER_OUT
        or die "Cannot open file $results_file: $!";

}

sub main {

    my $title = qq(PPPP Checking Primers for Correct Integration);

    my $cgi   = CGI->new();

    my $sw    = SangerWeb->new({'title'   => $title,
                                'banner'  => $title,
                                'inifile' => qq($ENV{'DOCUMENT_ROOT'}/PostGenomics/S_pombe/header.ini),
                               });

#    my $primer3 = $sw->server_root() . '/bin-offline/primer3';
my $primer3 = $sw->server_root() . '/cgi-bin/primer3-2.0.0-alpha/src/primer3_core';

    print $sw->header();

    my $for_sel = 3;
    my $rev_sel = 3;

    if ($cgi->param() ) {

        $for_sel   = $cgi->param('for_sel') || "";
        ($for_sel) = $for_sel =~ /([F\d]+_[ACGT]+_[ACGT]*)/;

        $rev_sel   = $cgi->param('rev_sel') || "";
        ($rev_sel) = $rev_sel =~ /([R\d]+_[ACGT]+_[ACGT]*)/;

    }

    print
        start_form(),
        p(qq(Recommended defaults are provided.<br />) ),
        p(qq(Gene Name @{[textfield('gene')]}) ),
        p(qq(Targeting primer length @{[textfield('length', '80', 3)]}) ),
        p(qq(Targeting primer increment @{[textfield('increment', '40', 3)]}) ),
        p(qq(Number of basepairs up- or downstream of <br />the target sequence to use for primer search
             @{[textfield('add_seq', 400, 3)]}) ),
        p(qq(Optimum primer length @{[textfield('opt_len', '22', 3)]}) ),
        p(qq(Minimum primer length @{[textfield('min_len', '20', 3)]}) ),
        p(qq(Maximum primer length @{[textfield('max_len', '28', 3)]}) ),
        p(qq(Optimum Tm @{[textfield('opt_tm', '60.0', 3)]}) ),
        p(qq(Minimum Tm @{[textfield('min_tm', '57.0', 3)]}) ),
        p(qq(Maximum Tm @{[textfield('max_tm', '63.0', 3)]}) ),
        p(qq(Minimum GC content @{[textfield('min_gc', '30', 3)]}) ),
        p(qq(Maximum GC content @{[textfield('max_gc', '70', 3)]}) ),
        hidden('for_sel', $for_sel),
        hidden('rev_sel', $rev_sel),
        submit(),
        reset(),
        end_form(),
        hr();

    if ($cgi->param() ) {

        my $gene    = $cgi->param('gene') || "";
        ($gene)     = $gene =~ /([a-zA-Z0-9_\-\.]+)/;

        my $len     = $cgi->param('length') || "";
        ($len)      = $len =~ /(\d+)/;

        my $inc     = $cgi->param('increment') || "";
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

        my ($for_seq, $for_tag) = $for_sel =~ /^[F\d]+_([ACGT]+)_([ACGT]*)$/;
        my ($rev_seq, $rev_tag) = $rev_sel =~ /^[R\d]+_([ACGT]+)_([ACGT]*)$/;

        print
            "Press 'Submit Query' button to get checking picking primers.",
            "<pre>Targeting forward primer:\n$for_seq - $for_tag\n\n",
            "Targeting reverse primer:\n$rev_seq - $rev_tag\n\n</pre>"
                if $for_sel;

        get_seq($primer3, $gene, $len, $inc, $add_seq, $opt_len, $min_len, $max_len,
                $opt_tm, $min_tm, $max_tm, $min_gc, $max_gc, $for_sel, $rev_sel)
            if defined $opt_tm;

    }

    print $sw->footer();

}
