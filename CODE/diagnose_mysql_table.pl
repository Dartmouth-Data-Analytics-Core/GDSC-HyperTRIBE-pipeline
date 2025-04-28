#!/usr/bin/env perl
use DBI; 
use strict; 
use warnings;
use Getopt::Std;

#Author Reazur Rahman 2017
#argument 1 tablename
#argument 2: list of matrix files already loaded in that MySQL table. The format of this argument is: filename1 filename2 and so on.

my $USAGE = "diagnose_mysql_table.pl -t tablename list_of_matrix_filenames_separate_by_space\n";
my %option;
getopts( 't:h', \%option );
my ($tablename);
if ( $option{t}) {
    $tablename = $option{t};
   
} else {
    die "proper parameters not passed\n$USAGE";
}

#MYSQL CONFIG VARIABLES
#my $host = "172.16.1.40"; # if mysql is hosted in a different machine
my $host = "dmseq-f11b-db.c.dartmouth.edu";
my $database = "dmseq";
my $user = "admin"; #mysql username
my $password = "gdscPass";

# connect to the mysql database
my $dsn = "DBI:mysql:$database:$host:3306;mysql_local_infile=1";
#my $dbh = DBI->connect($dsn, $user, $password, { RaiseError => 1 }) or die $DBI::errstr;

my $dbh = DBI->connect($dsn, $user, $password, { 
    RaiseError => 1, 
    mysql_read_default_file => '/CODE/mysql.cnf',
    mysql_ssl => 1,  
}) or die $DBI::errstr;

my ($sth);
#create the table
eval {
    $sth = $dbh->prepare("select COUNT(*) from $tablename");
    $sth->execute; 
};

if ($@) {
    print "Error in database creation: $@";
}
my $total=0;


print "\ttablename: $tablename: ";
my $tab_size;
while ( my @result = $sth->fetchrow_array ) { 
    print "$result[0]\n"; 
    $tab_size = $result[0];
} 
#die;
foreach my $file (@ARGV) {
    my $num_lines = `wc -l <$file`;
    chomp $num_lines;
    print "file: $file: $num_lines;\n";
    $total += $num_lines;
}


my $text = "";
if ($total == $tab_size) {
    $text = "SUCCESSFUL";
} else {
    $text = "FAILED";
}

print "\ttotal lines: $total\n\ttable load status: $text\n-------------------------\n";

#close connection to database
$sth->finish; 
$dbh->disconnect;
