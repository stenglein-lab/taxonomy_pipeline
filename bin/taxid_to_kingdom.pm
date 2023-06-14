#!/usr/bin/env perl
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

# this will be a connection to a database
my $_dbh = undef;

# my $cat_db="/home/databases/NCBI_Taxonomy/categories_dmp.sqlite3";

# connect to SQLite database if not already connected
sub _open_db_connection
{
   my $db_path = shift @_;

   # check that the SQlite database file exist
   if (! -e $db_path)
   {
      die("ERROR: NCBI categories.dmp SQLite file ($db_path) does not exist or is not accessible.\n");
   }

   # only open connection if not already connected
   if (!defined $_dbh)
   {
      # warn "opening dbi connection to $db_path.\n";
      $_dbh = DBI->connect("DBI:SQLite:dbname=$db_path",
                          '','',
                          {'RaiseError' => 1}) or die $DBI::errstr;
   }
   return $_dbh;
}

# map taxids to kingdoms as defined in the NCBI taxonomy categories.dmp file
sub taxid_to_kingdom
{
   my $db_file = shift @_;
   my $dbh = _open_db_connection($db_file);

   my @taxids = @_;
   my %taxid_kingdom_map = ();

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

   # check that all taxids are accounted for
   foreach my $tid (@taxids)
   {
      if (!$taxid_kingdom_map{$tid})
      {
         # warn "WARNING: no kingdom for taxid $tid\n";
         # create a place for taxid in hash, with undefined value
         $taxid_kingdom_map{$tid} = undef;
      }
   }

   # return a ref to hash with results
   return \%taxid_kingdom_map;
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
