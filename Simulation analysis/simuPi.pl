
## This script calculates Pi for every simulated segment for both A.th and C.ru

## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn

use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;
use Array::Utils qw(:all);

my $out_file = "XXXXXXXXXXXXXXXX/simuPiFst-100";
open(OUT,">$out_file") or die "Cannot open $out_file to write.\n";
print OUT "File\t Pi_A\t Pi_C\t Fst\n";

for ($i = 1; $i < 1000001; $i++){
	@seqA = ();
  @seqC = ();
	my $in = "XXXXXXXXXXX./fa/$i.fa"; 
 
  $seqio_object = Bio::SeqIO -> new(-file => $in);  
  $seq_object = $seqio_object -> next_seq;
  $len = $seq_object -> length;       

  
  $seqio_object = Bio::SeqIO -> new(-file => $in);
	while($seq_object = $seqio_object -> next_seq){      
      $title = $seq_object -> display_id;
      $seq = $seq_object -> seq;
          
      $ini = substr($title,0,2);
      
      if ($ini eq "1_"){
      	#$species = 0;
        push @seqA, $seq;
      }
      if ($ini eq "2_") {
      	#$species = 1;
        push @seqC, $seq;
      }
      
  }
  
  ## Calculation of Pi starts here....
	$countA = 0;
	$countC = 0;
	$count_between = 0;
	
  for ($a1 = 0; $a1 < 81; $a1 ++){
  	for ($a2 = 0; $a2 < $a1; $a2 ++){
  		$countA += &compare ($seqA[$a1], $seqA[$a2]);
  		
  	}    	
  }
  for ($c1 = 0; $c1 < 23; $c1 ++){
  	for ($c2 = 0; $c2 < $c1; $c2 ++){
  		$countC += &compare ($seqC[$c1], $seqC[$c2]);
  	}
  	
  }

  $Pi_A = $countA / ((80+1)*80/2) / 100;
  $Pi_C = $countC / ((22+1)*22/2) / 100;
  
  
  print OUT "$i\t". sprintf("%.5f", $Pi_A)."\t".sprintf("%.5f", $Pi_C)."\n";

	
				
}

sub compare(){
    #my %hTable = %{$_[1]};
    $seq1 = $_[0];
    $seq2 = $_[1];
    $num = 0;

    for ($n = 0; $n < 100; $n ++){
    	$char1 = substr($seq1, $n, 1);
    	$char2 = substr($seq2, $n, 1);
    	if ($char1 ne $char2) {
    		$num ++;    		
    	}
    }
    return $num;
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
