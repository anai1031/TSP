
## This script calculates the MAF for all the SNPs appearing in the simulated segments for both A.th and C.ru

## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn

use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;
use Array::Utils qw(:all);

my $out_file_a = "XXXXXXXXXXX/MAF-Ath";
my $out_file_c = "XXXXXXXXXXX/MAF-Cru";

open(OUTA,">$out_file_a") or die "Cannot open $out_file_a to write.\n";
open(OUTC,">$out_file_c") or die "Cannot open $out_file_a to write.\n";

print OUTA "File\t Num\t MAFa\n";
print OUTC "File\t Num\t MAFc\n";

for ($i = 1; $i < 1000001; $i++){
	@seqA = ();
  @seqC = ();
	my $in = "XXXXXXXXXXXX/fa/$i.fa"; 
 
  $seqio_object = Bio::SeqIO -> new(-file => $in);  
  $seq_object = $seqio_object -> next_seq;
  $len = $seq_object -> length;       # get the length;
  
  for ($a = 0; $a < 2; $a ++){
    for ($b = 0; $b < $len; $b++){
        my %alphabet;
        $freq[$a][$b] = \%alphabet;
    }
  }
  
  $seqio_object = Bio::SeqIO -> new(-file => $in);
	while($seq_object = $seqio_object -> next_seq){      
      $title = $seq_object -> display_id;
      $seq = $seq_object -> seq;
          
      $ini = substr($title,0,2);
      
      if ($ini eq "1_"){
      	$species = 0;
        push @seqA, $seq;
      }
      if ($ini eq "2_") {
      	$species = 1;
        push @seqC, $seq;
      }
      
      for ($m = 0; $m < $len; $m++){   	  
    	  $char = substr($seq,$m,1);
      	$hashRef = $freq[$species][$m];
      	&accumu($char,$hashRef);
      }
  }
  $countA = 0;
  $countC = 0;
  
    for ($k = 0; $k < $len; $k ++){
      	
          @keysAth = keys %{$freq[0][$k]};
          @keysCru = keys %{$freq[1][$k]};
        
          @valuesAth = values %{$freq[0][$k]};  
          @valuesCru = values %{$freq[1][$k]};     

          $athSize = @keysAth;   
          $cruSize = @keysCru;
        
          if($athSize == 2){   # both of size 2
        	  $countA ++;
        	  if ($valuesAth[0] < $valuesAth[1]){
        	  	$maf = $valuesAth[0] / 81;       	  	
        	  }else{
        	    $maf = $valuesAth[1] / 81;  	
        	  }
        	  print OUTA "$i\t$countA\t".sprintf("%.2f", $maf)."\n";
          }
        
          if($cruSize == 2){ 
          	$countC ++;
          	if ($valuesCru[0] < $valuesCru[1]){
        	  	$maf = $valuesCru[0] / 23;       	  	
        	  }else{
        	    $maf = $valuesCru[1] / 23;  	
        	  }
        	  print OUTC "$i\t$countC\t".sprintf("%.2f", $maf)."\n";
          }  
        	
    } # end of analyzing the specific region of the current sequence.    
	
		
}

sub accumu(){
    #my %hTable = %{$_[1]};
    my $hRef = $_[1];

    if (exists $hRef -> {$_[0]}){
       $$hRef{$_[0]} = $$hRef{$_[0]} + 1;
    }else{
       $$hRef{$_[0]} = 1;
    }
}

close OUT;
