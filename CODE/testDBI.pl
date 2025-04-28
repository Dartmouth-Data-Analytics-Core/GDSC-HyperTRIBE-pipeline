#!/usr/bin/perl
use strict;
use warnings;
use DBI;

# Changes to get SQL to work (potentially)
##### To login interactively
#mysql -h dmseq-f11b-db.c.dartmouth.edu -u admin -p --ssl-mode=REQUIRED
#CREATE DATABASE dmseq;
#SET GLOBAL local_infile = 1; 
    # The above line is for this error: DBD::mysql::st execute failed: Loading local data is disabled; this must be enabled on both the client and server sides at CODE/load_matrix_data.pl line 60.


#my $host = "hypertribe-d174-db.c.dartmouth.edu";
#my $database = "dmseq";
#my $user = "gdsc"; #mysql username
#my $password = "gdsc1227"; #mysql password, if any


my $host = "dmseq-f11b-db.c.dartmouth.edu";
my $database = "dmseq";
my $user = "admin"; #mysql username
my $password = "gdscPass";

#my $dsn = "DBI:mysql:$database:$host:3306"; 
#my $dbh = DBI->connect($dsn, $user, $password, {RaiseError=>1} ) or die $DBI::errstr;

my $dsn = "DBI:mysql:database=$database;host=$host;port=3306;mysql_ssl=1";
my $dbh = DBI->connect($dsn, $user, $password, { RaiseError => 1 }) or die $DBI::errstr;

print "✅ Connected to MySQL database '$database' on host '$host'\n";


