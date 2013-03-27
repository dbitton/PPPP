#!/usr/local/bin/perl -wT
#########
# Author: zeb
# Maintainer: cjp
# Created: 
# Last Modified: 2003-09-26 rmp: tidying
#
######################################################################
# get_pombe_c_term by Zoe Birtle zoebirtle@hotmail.com
# For use with S.pombe primer design programs
# cgi script to call  "run_pombe_primers_ctag_final.pl"
# Tuesday 10th June 2003
#####################################################################

#####################################################################
# the perl CGI module is used here to provide a user-friendly web-interface
# four values can be specified by the user.
# gene name, plasmid, primer length and increment size.
#####################################################################

use strict;
use CGI qw/:cgi :standard/;
#use lib "/nfs/WWW/SANGER_docs/perl";
use SangerPaths qw(core);
use SangerWeb;
use vars qw($EXE);

#$EXE         = "/nfs/WWWdev/SANGER_docs/data/PostGenomics/S_pombe/PPPP/run_pombe_primers_ctag_final2.pl";
$EXE         = SangerWeb->data_root().qq(/PostGenomics/S_pombe/PPPP/run_pombe_primers_ctag_final2.pl);
$ENV{'PATH'} = "";

&main();
1;

sub main {

  my $title = qq(Get Pombe Primers for C-terminal tagging);
  my $cgi   = CGI->new();
  my $sw    = SangerWeb->new({
			      'title'   => $title,
			      'banner'  => $title,
			      'inifile' => qq($ENV{'DOCUMENT_ROOT'}/PostGenomics/S_pombe/header.ini),
#                              'inifile' => qq(/nfs/WWW/SANGER_docs/htdocs/PostGenomics/S_pombe/header.ini),
			     });

  print $sw->header();

  print
    start_form(),
    p(qq(The help file for this program is available <A HREF='/PostGenomics/S_pombe/PPPP/help_file_C_term.shtml' TARGET='resource window'>here</A>)),
    p(qq(Enter <B>name of gene</B> to tag @{[textfield('name')]})),
    p(qq(Enter desired <B>length of target sequence</B> @{[textfield('length', '80', 3)]}
	   (i.e. length of primer excluding plasmid specific sequence))),
    p(qq(Which plasmid did you use as a <B>PCR template</B>?)),
    radio_group(-name    => 'plasmid',
		-default => 'pFA6a',
		-values  => ['pFA6a', 'other'],),
    p(qq(The forward primer is designed to allow C-terminal tagging of full length proteins.
	   Five reverse primers sequences will be suggested to you. The first primer is directly
	   downstream of the stop codon. The remaining reverse primers are positioned away from the
	   ORF in user-defined increments.)),
    p(qq(How far away (base pairs from ORF) would you like each subsequent primer to be?)),
    p(qq(<B>Primer increment</B> @{[textfield ('increment', 40, 3)]})),
    p(qq(Long strings of identical bases are highlighted in <FONT COLOR ='red'>colour</FONT>.)),
    submit(),
    reset(),
    end_form(),
    hr();

##########################################################################
# user input values, as collected from the cgi web interface, are sent to the 
# program run_pombe_primers_ctag_final.pl as options
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
    
    my $cmd = "$EXE  -g $gene  -p $plasmid -l $length -i $increment";
    #print "[$cmd]\n"; # debugging print statement
    system($cmd);
  }

  print $sw->footer();
}
