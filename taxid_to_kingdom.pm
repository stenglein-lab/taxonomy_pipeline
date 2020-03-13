#!/usr/bin/perl
#
# This script converts a taxid at species level or lower
# into a kingdom 
#
# It uses the information stored in the NCBI categories.dmp file, which is part of the NCBI taxonomy database.
# 
# This information is retreived from a local sqlite database 
#
# Mark Stenglein - June 8, 2011
#
# Updates:
#
#   Feb 25, 2020: convert to using sqlite db
#

package taxid_to_kingdom;

use DBI;
use strict;

use base 'Exporter';
our @EXPORT = qw(taxid_to_kingdom);

# my $cat_filename="/home/databases/NCBI_Taxonomy/categories.dmp";
my $cat_db="/home/databases/NCBI_Taxonomy/categories_dmp.sqlite3";


sub taxid_to_kingdom
{
   my @taxids = @_;
   my @kingdoms = ();
	my %taxid_kingdom_map = ();


   # connect to sqlite database
   my $dbh = DBI->connect("DBI:SQLite:dbname=$cat_db",
                          '','',
                          {'RaiseError' => 1}) or die $DBI::errstr;

   # determine kingdoms for these taxids
   my $num_taxids = scalar @taxids;
   my $qs = ("?");
   if ($num_taxids > 1)
   {
      my @qs_array = ("?") x (scalar @taxids);
      $qs = join (", ", @qs_array);
   }
   my $sql_string = "SELECT kingdom,taxid FROM categories_dmp where taxid in ( $qs )";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute(@taxids);

   # iterate through rows of mysql output
   # need to do this because:
   #
   # (1) some taxids might not return kingdoms
   # (2) results not necessarily in order of input
   #
   while ( my ($kingdom, $taxid) = $sth->fetchrow_array())
   {
      $taxid_kingdom_map{$taxid} = $kingdom;
   }

   # create array to return
   # will return an array that is 1:1 with input array
   foreach my $taxid (@taxids)
   {
      push @kingdoms, $taxid_kingdom_map{$taxid};
   }

	$dbh->disconnect;

   return @kingdoms;
}


# Perl packages must return true value
1; 


# The taxcat dump contains a single file -
# 
#   categories.dmp
# 
# categories.dmp contains a single line for each node
# that is at or below the species level in the NCBI 
# taxonomy database.
# 
# The first column is the top-level category -
# 
#   A = Archaea
#   B = Bacteria
#   E = Eukaryota
#   V = Viruses and Viroids
#   U = Unclassified and Other
# 
# The third column is the taxid itself,
# and the second column is the corresponding
# species-level taxid.
# 
# These nodes in the taxonomy -
# 
#   242703 - Acidilobus saccharovorans
#   666510 - Acidilobus saccharovorans 345-15 
# 
# will appear in categories.dmp as -
# 
# A       242703  242703
# A       242703  666510
# 
