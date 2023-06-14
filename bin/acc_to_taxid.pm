# 
# The functions in this package return NCBI Taxids for NCBI accessions #
# Depends on sorted_file_search_by_field.pm package
# Depends also on local copy of NCBI taxonomy db (see below)
# Depends on perl DBI package (to install: cpan DBI)
#
# Mark Stenglein, June 7, 2011
# Updated: 2/25/2020 to use sqlite3 db instead of mysql
#

package acc_to_taxid;

use DBI;
use strict;
use Data::Dumper;

use base 'Exporter';
our @EXPORT = qw(acc_to_taxid);

# this will be a connection to a database
my $_dbh = undef;

# Map NCBI accessions -> NCBI Taxids
# 
# This function takes 2 arguments:
# 
# 1) A path to a SQLite3 database containing a map of NCBI accessions to taxids
# 2) One or more NCBI accessions  
# 
# It returns a hash mapping accessions -> taxids 
# 
sub acc_to_taxid
{
   my $db_file = shift @_;

   my $dbh = _open_db_connection($db_file);

   # $taxids is a hash ref
   my $taxids = _acc_to_taxid($dbh, @_);

   return $taxids;
}

# connect to SQLite database if not already connected
sub _open_db_connection
{
   my $db_path = shift @_;

   # check that the SQlite database file exist
   if (! -e $db_path)
   {
      die("ERROR: NCBI Taxonomy SQLite file ($db_path) does not exist or is not accessible.\n");
   }

   # only open connection if not already connected
   if (!defined $_dbh) 
   {
      warn "opening dbi connection to $db_path.\n";
      $_dbh = DBI->connect("DBI:SQLite:dbname=$db_path",
                          '','',
                          {'RaiseError' => 1}) or die $DBI::errstr;
   }
   return $_dbh;
}

# map accessions -> taxids
sub _acc_to_taxid
{
   my $dbh = shift @_;
   my @accs = @_;
   my %acc_taxid_map = ();

   # trivial case: no input
   if (scalar @accs == 0) { return %acc_taxid_map; }

   # check for duplicates
   my %dups = ();
   foreach my $a (@accs)
   {
      $dups{$a} += 1;
      if ($dups{$a} > 1)
      {
         warn "WARNING: input to acc_to_taxid: accession $a is duplicated\n";
      }
   }

   # this will keep track of inputs passes as accession.version
   # will only do one lookup, for accession without version
   my %acc_ver_map = ();

   # check to see if accession in gi|XXXXX| format 
   for (my $i = 0; $i < (scalar @accs); $i++)
   {
      # if this gi is in the form it is in NCBI BLAST results
      if ($accs[$i] =~ /gi\|(\d+)\|/)
      {
         die ("error converting NCBI accessions to taxids: unsupported gi| format")
      }

      # convert from acc.version format to just acc
      if ($accs[$i] =~ /(\S+)\.\d+/)
      {
         # keep track of accession.versions for which we dropped version
         $acc_ver_map{$1} = $accs[$i];

	 # shorten to just accession (drop version)
         $accs[$i] = $1;
      }
   }

   # lookup taxids: prepare SQL
   my $num_accs = scalar @accs;
   my $qs = ("?");
   if ($num_accs > 1)
   {
      # this creates a bunch of ? that DBI will fill in with values below in prepare/execute
      my @qs_array = ("?") x (scalar @accs);
      $qs = join (", ", @qs_array);
   }
   
   my $sql_string = "SELECT acc, taxid FROM acc_taxid_map where acc in ( $qs )";
   my $sth = $dbh->prepare( $sql_string );
   $sth->execute(@accs);
   
   # iterate through rows of mysql output
   # need to do this because:
   #
   # (1) some ACCs might not return TAXIDs
   # (2) results not necessarily in order of input
   #
   while ( my ($acc, $taxid) = $sth->fetchrow_array())
   {
      if ($acc_taxid_map{$acc}) 
      {
         # not expecting >1 taxid per accession (expecting 1:1 mapping)
         die ("ERROR: more than one taxid defined for accession $acc.\n");
      }
      $acc_taxid_map{$acc} = $taxid;
   }

   # check that all accessions are accounted for
   foreach my $acc (@accs)
   {
      if (!$acc_taxid_map{$acc}) 
      { 
         warn "WARNING: no taxid for accession: $acc\n";
	 # create a place for accession in hash, with undefined value
         $acc_taxid_map{$acc} = undef; 
      }
      
      # if someone passed in accession.version, return acccession.version
      if ($acc_ver_map{$acc}) 
      { 
         my $acc_ver = $acc_ver_map{$acc}; 
         $acc_taxid_map{$acc_ver} = $acc_taxid_map{$acc}; 
	 # get rid of no-longer needed accession (without version) hash element
         delete($acc_taxid_map{$acc}); 
      }
   }
   
   # return a reference to the hash mapping acc->taxid
   return \%acc_taxid_map;
}

# PERL packages must return a true value
1;
