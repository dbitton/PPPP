#!/usr/local/bin/perl -wT
#########
# Author: zeb
# Maintainer: cjp
# Created: 
# Last Modified: 2003-09-26 rmp: tidying
#
#####################################################################
# get_pombe_n_term by Zoe Birtle zoebirtle@hotmail.com
# For use with S.pombe primer design programs
# cgi script to call  "run_pombe_primers.pl"
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

$EXE         = SangerWeb->data_root().qq(/PostGenomics/S_pombe/PPPP/run_pombe_primers_ntag_final2.pl);
$ENV{'PATH'} = "";

&main();
1;

sub main {

  my $title = qq(Get Pombe Primers for Controlling Genes by Regulatable Promoter and for N-terminal Tagging);
  my $cgi   = CGI->new();
  my $sw    = SangerWeb->new({
			      'title'   => $title,
			      'banner'  => $title,
			      'inifile' => qq($ENV{'DOCUMENT_ROOT'}/PostGenomics/S_pombe/header.ini),
			     });

  print $sw->header();
  
  
  print
    start_form(),
    p(qq(The help file for this program is available  <A HREF='/PostGenomics/S_pombe/PPPP/help_file_N_term.shtml' TARGET='resource window'> here</A>\n)),
    p(qq(Enter <B>name of gene</B> to tag @{[textfield('name')]})),
    p(qq(Enter desired <B>length of primer target sequence</B> @{[textfield('lengthf', 80, 3)]})),
    p(qq(The reverse primer is designed to allow controlled expression or N-terminal tagging of
	 full length proteins. Five forward primer sequences will be suggested to you. The first
	 primer is directly upstream of the stop codon. The remaining forward primers are
	 positioned away from the ORF in user-defined increments.)),
    p(qq(How far away (base pairs from ORF) would you like each subsequent primer to be?)),
    p(qq(<B>Primer increment </B> @{[textfield('increment', 40, 3)]})),
    p(qq(Which <B>N-Terminal tag</B> would you like to use?)), 
    radio_group(
                -name    => 'tag',
		-default => 'None',
		-values  => ['None','3HA','GST','GFP'],),
    p(qq(Long Strings of identical bases are highlighted in <FONT COLOR ='red'>colours</FONT>.\n)),
    submit(),
    reset(),
    end_form(),
    hr();

##########################################################################
# user input values, as collected from the cgi web interface, are sent to the 
# program run_pombe_primers_ntag_final.pl as options
#
##########################################################################
  
  if ($cgi->param()) {
    my $tag       = $cgi->param('tag') || "";
    ($tag)        = $tag =~ /([a-zA-Z0-9_\-\.]+)/;

    my $gene      = $cgi->param ('name') || "";
    ($gene)       = $gene =~ /([a-zA-Z0-9_\-\.]+)/;

    my $length    = $cgi->param ('length') || "";
    ($length)     = $length =~ /(\d+)/;

    my $increment = $cgi->param ('increment') || "";
    ($increment)  = $increment =~ /(\d+)/;

    my $cmd = "$EXE  -g $gene  -p $tag -l $length -i $increment";
    #print "[$cmd]\n"; # debugging print statement
    system($cmd);
  }

  print $sw->footer();
}
