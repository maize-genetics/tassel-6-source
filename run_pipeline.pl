#!/usr/bin/perl -w

use strict;
use File::Basename;

my $top = dirname(__FILE__);
$top //= '.';

# Scan @ARGV for Java memory arguments, and put rest in @args
my $java_mem_min_default = "-Xms512m";
my $java_mem_max_default = "-Xmx1536m";
my $java_mem_min = "";
my $java_mem_max = "";
my @args;
for (my $i=0; $i<=$#ARGV; $i++){
   if ($ARGV[$i] =~ m/Xms/) {
      $java_mem_min .= "$ARGV[$i]";
      $java_mem_min=~s/–/-/g;
   }
   elsif ($ARGV[$i] =~ m/Xmx/) {
      $java_mem_max .= "$ARGV[$i]";
      $java_mem_max=~s/–/-/g;
   }
   else{
      push(@args, $ARGV[$i]);
   }
}

if ($java_mem_min eq "") { $java_mem_min = $java_mem_min_default; }
if ($java_mem_max eq "") { $java_mem_max = $java_mem_max_default; }

$ENV{'_JAVA_OPTIONS'} = "$java_mem_min $java_mem_max";

print "Memory Settings: $java_mem_min $java_mem_max\n";
print "Tassel Pipeline Arguments: " . "@args\n";

system "time ./gradlew -q run --args='-debug @args'";
