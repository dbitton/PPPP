#!/usr/local/bin/perl -wT
#####################################################################
#
# Author: cjp
# Maintainer: cjp
# Created: 2006-10-26
# Last Modified: 2006-10-26 cjp: re-write
#
#####################################################################

use strict;
use CGI qw/:cgi :standard/;
use CGI::Carp qw(fatalsToBrowser);

use SDBM_File;
use Fcntl;

use lib qw(/var/www/lib/core);
use SangerWeb;

use vars qw($na_conc);

$na_conc = 0.2;


$ENV{'PATH'} = "";

&main();
1;

sub main {

  my $title = qq(PPPP Primers for C-terminal tagging);

  my $cgi   = CGI->new();

  my $sw    = SangerWeb->new({'title'   => $title,
                              'banner'  => $title,
                             # 'inifile' => qq($ENV{'DOCUMENT_ROOT'}/header.ini),
                             });

#  print $sw->header();
print qq(Content-type: text/html\n\n);

print qq (<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Bahler Lab</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link href="/css/other.css" rel="stylesheet" type="text/css" />
</head>
<body>
<div id="header">
<div style="text-align:left; margin-left:70px; margin-bottom:0px; margin-top:0px;"><span style="text-shadow:#FFF; font-family:'Arial Black', Gadget, sans-serif; font-size:32pt; color:#EBF1DE; ">b&auml;hler</span><span style="text-shadow:auto; font-family:'Arial Black', Gadget, sans-serif; font-size:32pt; color:#77933C; ">lab&#13;</span></div>
  <div style="text-align:left; margin-left:30px; margin-top:0px;"><span style="font-family:'Arial Black', Gadget, sans-serif; font-size:20pt; color:#4F6228; "> Genome Regulation</span></div>
  <div class="O">
  <div style="text-align:center; margin-bottom:0; margin-top:-26px; margin-left:120px;margin-right:-20px;"><span style="font-family:Arial, Helvetica, sans-serif; font-size:55pt; color:#D7E4BD; "><strong><em>primer programs&#13;</em></strong></span></div>
<div style="text-align:right; margin-bottom:0px; margin-top:-127px; margin-left:150px;"><span style="font-family:'Arial Black', Helvetica, sans-serif; font-size:36pt; color:white; ">
<img src="/images/uclogo.gif";/></span></div>
  
<div></div>
</div>


</div>
<div id="sidebar">
                <div id="menu">
                        <ul>
                                <li ><a href="/index.html" title="">Home</a></li>
                                <li><a href="/people" title="">People</a></li>
                                <li><a href="/research.htm" title="">Research</a></li>
                                <li><a href="/publications" title="">Publications</a></li>

                                <li class="active first"><a href="/resources" title="">Resources</a></li>

                                <li ><a href="/contact.htm" title="">Contact</a></li>
                        </ul>
                </div>


        </div>

         <div id="content">

                <div class="feature bg7">
        </div>

                        <div class="content" >
<table class="violet2" width="100%"><tr><td>

);

  print
    start_form(),
    p(qq(The help file for this program is available <A HREF='/PPPP/help_file_C_term.shtml' TARGET='resource window'> here</A>) ),
    p(qq(Enter the <B>name of gene</B> to delete @{[textfield('gene')]}) ),
    p(qq(Enter desired <B>length of primer target sequence</B> @{[textfield('length', 80, 3)]}) ),
    p(qq(Which plasmid will you use as a <B>PCR template</B>?) ),
    radio_group(-name    => 'plasmid',
                -default => 'pFA6a',
                -values  => ['pFA6a', 'other']),
    p(qq(The forward primer is designed to allow C-terminal tagging of full length proteins.
           Five reverse primers sequences will be suggested to you. The first primer is directly
           downstream of the stop codon. The remaining reverse primers are positioned away from the
           ORF in user-defined increments.) ),
    p(qq(How far away (base pairs from ORF) would you like each subsequent primer to be?) ),
    p(qq(<B>Primer increment</B> @{[textfield ('increment', 40, 3)]}) ),
    p(qq(Long strings of identical bases are highlighted in <FONT COLOR ='red'>colour</FONT>.) ),
    submit(),
    reset(),
    end_form(),
 p(qq(</td> <td><img src="/gfx/pppp.jpg" alight=right>
</td></tr><tr>)),
 p(qq(</tr></table>)),

    hr();

  if ($cgi->param() ) {

    my $gene      = $cgi->param('gene') || "";
    ($gene)       = $gene =~ /([a-zA-Z0-9_\-\.]+)/;

    my $plasmid   = $cgi->param('plasmid') || "";
    ($plasmid)    = $plasmid =~ /([a-zA-Z0-9_\-\.]+)/;

    my $length    = $cgi->param('length') || "";
    ($length)     = $length =~ /(\d+)/;

    my $increment = $cgi->param('increment') || "";
    ($increment)  = $increment =~ /(\d+)/;

    run_primer_c_term($gene, $plasmid, $length, $increment);

  }

  print "</pre>";
  print '</div></div>'; 
# print $sw->footer();

}


sub run_primer_c_term {

  my $gene      = shift;
  my $plasmid   = shift;
  my $length    = shift;
  my $increment = shift;

  my $datadir = SangerWeb->data_root() . '/PPPP/';

  $gene = 0
    unless defined $gene;

  $plasmid = 0
    unless defined $plasmid;

  $length = 80
    unless defined $length;

  $increment = 40
    unless defined $increment;

  my $start = 0;

  undef my %genehash;
  my $dbm_file = $datadir . 'pombe_data_cjp.dbr';

  tie %genehash, 'SDBM_File', $dbm_file, O_RDONLY, 0666
    or die "Cannot tie DMB file $dbm_file: $!";

  my ($chr, $complement, $start, $end);

  foreach my $key (keys %genehash) {

    if (lc $key eq lc $gene) {

      ($chr, $complement, $start, $end) = split /,/, $genehash{$key};

    }

  }

  untie %genehash
    or die "Cannot untie DMB file $dbm_file: $!";

  if (not $gene or $gene eq '-p' or not $start) {

    $gene = ''
      if $gene == 0;

    print "The gene $gene was not found in the database.\n";

    return 0;

  }

  my $gene_length = ($end - $start) + 1;

  print "<pre>The gene is on chromosome $chr.\n\n";

  if (not $complement) {

    print "The gene is on the forward strand.\n\n";

  } else {

    print "The gene is on the reverse strand.\n\n";

  }

  print
    "The length of the gene is $gene_length (including all introns and exons).\n\n",
    "The first pair of forward primers is adjacent to the ORF.\n",
    "Each subsequent pair of primers is $increment bp further from the ORF than the last.\n",
    "The reverse primer is internal to the ORF and ends at the start codon.\n\n";

  my $primer_dist2 = 2*$increment;
  my $primer_dist3 = 3*$increment;
  my $primer_dist4 = 4*$increment;

  my @primer_dist = (0, $increment, $primer_dist2, $primer_dist3, $primer_dist4);

  my $chr_file = $datadir . 'all_chromosomes.seq';

  open CHR_FILE, $chr_file
    or die "Cannot open sequence file $chr_file: $!";

  my $chr_seq = '';

  my $primer_forward = '';

  undef my @primer_reverse;

  my $start_codon = '';

  while (defined (my $line = <CHR_FILE>) ) {

    if (substr ($line, 0, 5) eq ">CHR${chr}") {

      $chr_seq = <CHR_FILE>;

      $chr_seq = uc $chr_seq;

      if (not $complement) {

        foreach my $i (0 .. 4) {

          $primer_reverse[$i] = substr $chr_seq, $end + $primer_dist[$i], $length;
          $primer_reverse[$i] = rev_comp($primer_reverse[$i]);

        }

        $primer_forward = substr $chr_seq, $end - $length - 3, $length;

      } else {

        foreach my $i (0 .. 4) {

          $primer_reverse[$i] = substr $chr_seq, $start - ($length + $primer_dist[$i] + 1), $length;

        }

        $primer_forward = substr($chr_seq, $start + 2, $length);
        $primer_forward = rev_comp($primer_forward);

      }

    }

  }

  close CHR_FILE
    or die "Cannot close sequence file $chr_file: $!";

  my $primer_length = $length;

  my $forward_tag = '';
  my $reverse_tag = '';

  if (lc $plasmid eq 'pfa6a') {

    print "Primers for $gene using plasmid $plasmid.\n";

    $primer_length = 20 + $length;

    $forward_tag = 'CGGATCCCCGGGTTAATTAA';
    $reverse_tag = 'GAATTCGAGCTCGTTTAAAC';

  } else {

    print "Primers for $gene with no attached plasmid sequence.\n";

  }

  print
      "\nTo design checking primers, select your forward and reverse targeting ",
      "primers \nand then press the button at the bottom.";

  undef my @gc_tm_reverse;

  my $primer_forward2 = $primer_forward . $forward_tag;

  my ($gc_forward, $primer_tm) = gc_content($primer_forward2);

  my $gc_tm_forward = [$gc_forward, $primer_tm];

  foreach my $primer (@primer_reverse) {

    my $primer_reverse = $primer . $reverse_tag;

    my ($gc_reverse, $primer_tm) = gc_content($primer_reverse);

    push @gc_tm_reverse, [$gc_reverse, $primer_tm];

  }

  #print "[$gene]\n";

  print
    start_form(-action => "/cgi-bin/PPPP/pppp_checking.pl"),
    hidden('gene', $gene),
    hidden('length', $length),
    hidden('increment', $increment);

  my $for = 'F_' . $primer_forward . '_' . $forward_tag;;

  $primer_forward =~ s/(A{4,})/<FONT COLOR='red'>$1<\/FONT>/g;
  $primer_forward =~ s/(C{4,})/<FONT COLOR='blue'>$1<\/FONT>/g;
  $primer_forward =~ s/(G{4,})/<FONT COLOR='pink'>$1<\/FONT>/g;
  $primer_forward =~ s/(T{4,})/<FONT COLOR='green'>$1<\/FONT>/g;

  my $check = "checked='checked'";

  my $gc_content_for = sprintf "%.1f", $gc_tm_forward->[0];
  my $primer_tm_for  = sprintf "%.3f", $gc_tm_forward->[1];

  print
    "\n<label><input type='radio' name='for_sel' value='$for' $check tabindex='3'/>",
    "Forward Primer = $primer_forward - $forward_tag\t",
    "% GC Content = $gc_content_for\t",
    "TM = $primer_tm_for</label>";

  print "\n";

  foreach my $i (0 .. 4) {

    my $rev_i = $i . '_' . $primer_reverse[$i] . '_' . $reverse_tag;

    $primer_reverse[$i] =~ s/(A{4,})/<FONT COLOR='red'>$1<\/FONT>/g;
    $primer_reverse[$i] =~ s/(C{4,})/<FONT COLOR='blue'>$1<\/FONT>/g;
    $primer_reverse[$i] =~ s/(G{4,})/<FONT COLOR='pink'>$1<\/FONT>/g;
    $primer_reverse[$i] =~ s/(T{4,})/<FONT COLOR='green'>$1<\/FONT>/g;

    my $check = '';

    $check = "checked='checked'"
      if $i == 0;

    my $gc_content_rev = sprintf "%.1f", $gc_tm_reverse[$i][0];
    my $primer_tm_rev  = sprintf "%.3f", $gc_tm_reverse[$i][1];

    my $i1 = $i + 1;

    print
      "\n<label><input type='radio' name='rev_sel' value='$rev_i' $check tabindex='3'/>",
      "Reverse Primer $i1 = $primer_reverse[$i] - $reverse_tag\t",
      "% GC Content = $gc_content_rev\t",
      "TM = $primer_tm_rev\t",
      "Distance from ORF = $primer_dist[$i]</label>";

  }

  print
    "\n\n",
    submit('Get checking primers'),
    reset(),
    end_form(),
    hr();

}


sub gc_content {

  my $primer = shift;

  my $primer_length = length $primer;

  my @primer = split '', $primer;

  my $count_of_A = 0;
  my $count_of_T = 0;
  my $count_of_C = 0;
  my $count_of_G = 0;

  my $errors = 0;

  foreach my $res (@primer) {

    if (uc $res eq 'A') {

      ++$count_of_A;

    } elsif (uc $res eq 'G') {

      ++$count_of_G;

    } elsif (uc $res eq 'C') {

      ++$count_of_C;

    } elsif (uc $res eq 'T') {

      ++$count_of_T;

    } else {

      ++$errors;

    }

  }

  my $total = $count_of_A + $count_of_C + $count_of_G + $count_of_T;

  my $percent_C = ($count_of_C/$total)*100;
  my $percent_G = ($count_of_G/$total)*100;

  my $G_C_content = $percent_G + $percent_C;

  my $primer_tm = tm_calculate($G_C_content, $primer_length);

  return $G_C_content, $primer_tm;

}


sub tm_calculate {

  my $G_C_content   = shift;
  my $primer_length = shift;

  my $log_Na_conc = log($na_conc)/log(10);

  my $primer_tm = (81.5 + (16.6*$log_Na_conc) + (0.41*($G_C_content) ) - (600/$primer_length) );

  return $primer_tm;

}


sub rev_comp {

  my $primer = shift;

  $primer = uc $primer;
  $primer = reverse $primer;
  $primer =~ tr/ACGT/TGCA/;

  return $primer;

}
