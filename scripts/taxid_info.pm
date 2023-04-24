#!/usr/bin/env perl
#
# This module returns taxonomic information about NCBI taxids: child nodes, parent node, rank, etc
#
# It uses the information stored in the NCBI nodes.dmp file, which is part of the NCBI taxonomy database.
# 
# This information is retreived from a local sqlite database 
#
# Mark Stenglein - March 27, 2023
#
#

package taxid_info;

use DBI;
use strict;

use base 'Exporter';
our @EXPORT = qw(taxid_descendents taxid_parent taxid_rank);

# this will be a connection to a database
my $_dbh = undef;

# connect to SQLite database if not already connected
sub _open_db_connection
{
   my $db_path = shift @_;

   # check that the SQlite database file exist
   if (! -e $db_path)
   {
      die("ERROR: NCBI nodes.dmp SQLite file ($db_path) does not exist or is not accessible.\n");
   }

   # only open connection if not already connected
   if (!defined $_dbh)
   {
      # warn "opening dbi connection to $db_path.\n";
      $_dbh = DBI->connect("DBI:SQLite:dbname=$db_path",
                          '','',
                          {'RaiseError' => 1}) or die $DBI::errstr;

      if(!$_dbh)
      {
         die "Error: failed to connect to MySQL database $db_path. $DBI::errstr";
      }
   }


   return $_dbh;
}

# return all descendent taxids of one taxid
sub taxid_descendents
{
   my $db_file = shift @_;
   my $dbh = _open_db_connection($db_file);

   my $taxid = shift @_;
   my @descendents = ();

   # retreive all descendents of this taxid 
   my $sql_string = "SELECT taxid FROM nodes_dmp where parent_taxid=?";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute($taxid);

   # iterate through rows of mysql output
   while (my ($child_taxid) = $sth->fetchrow_array())
   {
      # don't output a child if it's the same as the parent 
      # (this is true for root node: 1, which has itself for a parent)
      if ($child_taxid != $taxid)
      {
         push @descendents, $child_taxid
      }
   }

   # return ref to array with results
   return \@descendents;
}

# return parent taxid of one taxid
sub taxid_parent
{
   my $db_file = shift @_;
   my $dbh = _open_db_connection($db_file);

   my $taxid = shift @_;
   my $parent_taxid = undef;

   # retreive all descendents of this taxid 
   my $sql_string = "SELECT parent_taxid FROM nodes_dmp where taxid=?";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute($taxid);

   # iterate through rows of mysql output
   while (my ($parent_tid) = $sth->fetchrow_array())
   {
      if ($parent_taxid) 
      {
         die ("ERROR: was only expecting one parent taxid for taxid $taxid, but got >1\n");
      }
      $parent_taxid = $parent_tid;
   }

   # warn "$taxid --P--> $parent_taxid\n";

   # return the parent taxid
   return $parent_taxid;
}

# return the taxonomic rank of one taxid
sub taxid_rank
{
   my $db_file = shift @_;
   my $dbh = _open_db_connection($db_file);

   my $taxid = shift @_;
   my $rank = undef;

   # retreive rank of this taxid
   my $sql_string = "SELECT rank FROM nodes_dmp where taxid=?";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute($taxid);

   # iterate through rows of mysql output
   # make sure only one returned
   while (my ($this_rank) = $sth->fetchrow_array())
   {
      if ($rank) 
      {
         die ("ERROR: was only expecting one rank for taxid $taxid, but got >1\n");
      }
      $rank = $this_rank;
   }

   # return the rank
   return $rank;
}


# Perl packages must return true value
1; 


