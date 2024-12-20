#!/usr/bin/env perl

# Basics: is there a syntax error where already bash hiccups on?

use strict;
use Test::More;
use File::stat;

my $tests = 0;
my $fileout="";
my $prg="./testssl.sh";
my $out="";

# Try to detect remainders from debugging:
my $debug_regexp='^(\s)*set (-|\+)x';
# Blacklists we use to trigger an error:
my $error_regexp1='(syntax|parse) (e|E)rror';
my $error_regexp2='testssl.sh: line';
my $error_regexp3='bash: warning';
my $error_regexp4='command not found';
my $error_regexp5='(syntax error|unexpected token)';

printf "\n%s\n", "Testing whether just calling \"./testssl.sh\" produces no error ...";
my $info    = stat($prg);
my $retMode = $info->mode;

is($retMode & 0400, 0400, "Checking \"./testssl.sh\" for read permission");
$tests++;

is($retMode & 0100, 0100, "Checking \"./testssl.sh\" for execute permission");
$tests++;

$fileout = `timeout 10 bash $prg 2>&1`;
my $retval=$?;

unlike($fileout, qr/$error_regexp1/, "regex 1");
$tests++;

unlike($fileout, qr/$error_regexp2/, "regex 2");
$tests++;

unlike($fileout, qr/$error_regexp3/, "regex 3");
$tests++;

unlike($fileout, qr/$error_regexp4/, "regex 4");
$tests++;

unlike($fileout, qr/$error_regexp5/, "regex 5");
$tests++;

is($retval, 0, "return value should be equal zero: \"$retval\"");
$tests++;

$out=`grep -E "$debug_regexp" $prg`;
unlike($out, qr/$debug_regexp/, "Debug RegEx");
$tests++;

printf "\n\n";
done_testing($tests);


#  vim:ts=5:sw=5:expandtab

