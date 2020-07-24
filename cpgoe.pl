#! /usr/local/bin/perl
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "\nUsage: $0 -f <fasta file>  -h\n";
$usage .= "            -f   <input fasta file>\n";
$usage .= "            -h	Displays this help and exit\n\n";

my $fasta;
my $help;

GetOptions(
           'f=s'	  => \$fasta, 
           'h'        => \$help,
          );

if($help){
    print $usage;
    exit 0;
}

unless($fasta && -f $fasta){
    print STDERR "cannot access [$fasta] file !\n";
    print STDERR $usage;
    exit 1;
}



sub countbases_ {
	my $seq = shift;
	my $str=$seq->seq;
	my $nr_G=($str =~ tr/G/G/);		
	my $nr_C=($str =~ tr/C/C/);		
	my $nr_A=($str =~ tr/A/A/);		
	my $nr_T=($str =~ tr/T/T/);		
	my $nr_N=($str =~ tr/N/N/);	
 	my $G=sprintf("%.4f",$nr_G/$seq->length);
  	my $C=sprintf("%.4f",$nr_C/$seq->length);
  	my $A=sprintf("%.4f",$nr_A/$seq->length);
  	my $T=sprintf("%.4f",$nr_T/$seq->length);
  	my $N=sprintf("%.4f",$nr_N/$seq->length);	
	my @tot_cg=($str =~ /CG/g);   	#gets all occurences of "CG" in the sequence in a table
	my @tot_gc=($str =~ /GC/g);	#gets all occurences of "GC" in the sequence in a table
	my $nr1=@tot_cg;  			#scalar interpretation of the table, returns, the number of occurrences of "CG"
	my $nr2=@tot_gc;			#scalar interpretation of the table, returns, the number of occurrences of "GC"
	my $lg=length($str); 
	my $lg_no=length($str)-$nr_N;  #total number of bases A,C,G,T  (without Ns)
	my $cpg_oe;
	my $gpc_oe;
    	
	if ($nr_C!=0 && $nr_G!=0 ) {
    		$cpg_oe=sprintf("%.3f",$nr1/$nr_G/$nr_C*$lg_no);
		$gpc_oe=sprintf("%.3f",$nr2/$nr_G/$nr_C*$lg_no);
   	 } 
    	else {  $cpg_oe="na"; $gpc_oe="na" };    # cannot estimate these ratios if there are no Cs or no Gs
    
  	my $string="$A\t$C\t$G\t$T\t$N\t$cpg_oe\t$gpc_oe";
  	return($string)  
}
  




#MAIN
my $count=0; 
my $tab=$fasta.".basecomp_and_ratios.txt";
#&cluster_info_;
open COMP,">$tab";
print COMP join ("\t","snumber","id","ssize","#A","#C","#G","#T","#N","cpg_oe","gpc_oe","\n");

my $current=Bio::SeqIO->new(-file => $fasta, -format =>  'fasta'); 	
 	while (my $seq=$current->next_seq) {	
		$count++;
		print COMP join ("\t",$count,$seq->display_id,length($seq->seq),countbases_($seq),"\n");		

	}	
print "\n counted $count genes ";





