#!/usr/bin/env perl

# 
# The functions in this package return names
# corresponding to NCBI TAXIDs
# # Mark Stenglein, June 9, 2011
#

package taxid_to_name;

use DBI;
use strict;

use base 'Exporter';
our @EXPORT = qw(taxid_to_description taxid_to_scientific_name taxid_to_common_name);

# this will be a connection to a database
my $_dbh = undef;

# connect to SQLite database if not already connected
sub _open_db_connection
{
   my $db_path = shift @_;

   # check that the SQlite database file exist
   if (! -e $db_path)
   {
      die("ERROR: NCBI names.dmp SQLite file ($db_path) does not exist or is not accessible.\n");
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

sub taxid_to_scientific_name
{
   my $db_file = shift @_;
   my @taxids = @_;

   return taxid_to_name($db_file, "scientific name", @taxids);
}

sub taxid_to_common_name
{
   my $db_file = shift @_;
   my @taxids = @_;

   return taxid_to_name($db_file, "common name", @taxids);
}

sub taxid_to_name
{
   # first argument: path to sqlite database 
   my $db_file = shift @_;
   my $dbh = _open_db_connection($db_file);

   # second argument: type of name
   # typically: "scientific name" or "common name"
   my $name_class = shift @_;

   # rest of arguments: taxids
   my @taxids = @_;
   my %taxid_name_map = ();

   # determine names for these taxids
   my $num_taxids = scalar @taxids;
   my $qs = ("?");
   if ($num_taxids > 1)
   {
      my @qs_array = ("?") x (scalar @taxids);
      $qs = join (", ", @qs_array);
   }
   my $sql_string = "SELECT taxid,name_text FROM names_dmp where taxid in ( $qs ) AND name_class=?";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute(@taxids, $name_class);

   # iterate through rows of mysql output
   # need to do this because:
   #
   # (1) some taxids might not return names
   # (2) results not necessarily in order of input
   #
   my $counter = 1;
   while ( my ($taxid, $name) = $sth->fetchrow_array())
   {
      $taxid_name_map{$taxid} = $name;
      $counter++;
   }

   # check that all taxids are accounted for
   foreach my $tid (@taxids)
   {
      if (!$taxid_name_map{$tid})
      {
         # warn "WARNING: no $name_class for taxid $tid\n";
         # create a place for taxid in hash, with undefined value
         $taxid_name_map{$tid} = undef;
      }
   }

   # return a ref to hash with results
   return \%taxid_name_map;

}


# PERL packages must return a true value
1;
