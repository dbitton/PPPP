#!/usr/local/bin/perl -wT
#########
# Author: zeb
# Maintainer: cjp
# Created: 
# Last Modified: 2003-09-26 rmp: tidying
#
#####################################################################
# get_pombe_deletion_primers by Zoe Birtle zoebirtle@hotmail.com
# For use with S.pombe primer design programs
# cgi script to call  "run_pombe_primers_deletion_final.pl"
# Tuesday 10th June 2003
#####################################################################

####################################################################
# the perl CGI module is used here to provide a user-friendly web-interface
# four values can be specified by the user.
# gene name, plasmid, primer length and increment size.
#####################################################################

use strict;
use CGI qw/:cgi :standard/;
use SangerPaths qw(core);
use SangerWeb;
use vars qw($EXE);

$EXE         = SangerWeb->data_root() . "/PostGenomics/S_pombe/PPPP/run_pombe_primers_deletion_final2.pl";
#$EXE         = qq(/nfs/WWWdev/SANGER_docs/data/PostGenomics/S_pombe/PPPP/run_pombe_primers_deletion_final2.pl);
$ENV{'PATH'} = "";

&main();
1;

sub main {

  my $title = qq(Get Pombe Primers For Gene Deletion);
  my $cgi   = CGI->new();
  my $sw    = SangerWeb->new({
			      'title'   => $title,
			      'banner'  => $title,
			      'inifile' => qq($ENV{'DOCUMENT_ROOT'}/PostGenomics/S_pombe/header.ini),
			     });

  print $sw->header();

  print
    start_form(),
    p(qq(The help file for this program is available <A HREF='/PostGenomics/S_pombe/PPPP/help_file_deletion.shtml' TARGET='resource window'> here</A>)),
    p(qq(Enter the <B>name of gene</B> to delete @{[textfield('name')]})),
    p(qq(Enter desired <B>length of primer target sequence</B> @{[textfield('length', '80', 3)]})),
    p(qq(i.e. the length of primer excluding plasmid specific sequence.)),
    p(qq(Which plasmid will you use as a <B>PCR template?</B>)),
    radio_group(-name    => 'plasmid',
		-default => 'pFA6a',
		-values  => ['pFA6a', 'KS-ura4', 'Other'],),
    p(qq(Five pairs of primers will be suggested to you. The primers are directly
	 upstream of the ORF (not including the ORF sequence). The remaining pairs are positioned
	 away from the ORF in user-defined increments.
	 How far away (base pairs from ORF) would you like each subsequent primer set to be?
	 Minus values acceptable, which will result in primers encroaching into the ORF.)),
    p(qq(<B> Primer increment</B>  @{[textfield('increment', '40', 3)]})),
    p(qq(Long strings of identical bases are highlighted in<FONT COLOR ='red'> colour</FONT>.)),
    submit(),
    reset(),
    end_form(),
    hr();
    
##########################################################################
# user input values, as collected from the cgi web interface, are sent to the 
# program run_pombe_primers_deletion_final.pl as options
#
##########################################################################

  if ($cgi->param()) {
    my $plasmid   = $cgi->param('plasmid') || "";
    ($plasmid)    = $plasmid =~ /([a-zA-Z0-9_\-\.]+)/;

    my $gene      = $cgi->param ('name') || "";
    ($gene)       = $gene =~ /([a-zA-Z0-9_\-\.]+)/;

    my $length    = $cgi->param ('length') || "";
    ($length)     = $length =~ /(\d+)/;

    my $increment = $cgi->param ('increment') || "";
    ($increment)  = $increment =~ /(\d+)/;

    my $cmd = "$EXE -g $gene  -p $plasmid -h $length -i $increment";
    #print "[$cmd]\n"; debugging print statement
    system($cmd);
  }

  print $sw->footer();
}
