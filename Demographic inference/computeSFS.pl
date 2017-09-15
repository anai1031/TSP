
## This script computes the SFS (MAF) of the two populations to feed into the fastsimcoal software, for both the joint and Multi versions. 
## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn


use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;
#use Array::Utils qw(:all);

$countSites = 0;
@count = ();

for ($a = 0; $a <= 81; $a ++){
  for ($b = 0; $b <= 23; $b ++){
  	$count[$a][$b] = 0;
  }
}

$outFreq = "XXXXXXXXX/AthCru_jointMAFpop1_0.obs";
open(OUT,">$outFreq") or die "Cannot open $outFreq to read.\n";

$outFreq12 = "XXXXXXXXX/AthCru_MSFS.obs";
open(OUT12,">$outFreq12") or die "Cannot open $outFreq12 to read.\n";
print OUT12 "1 observations. No. of demes and sample sizes are on next line\n2\t81\t23\n";

print OUT "1 observation\n\t";
for ($a = 0; $a < 81; $a ++){
	print OUT "d0_".$a."\t";
}
print OUT "d0_81\n";


for ($i = 1; $i < 16046; $i ++) {
	
	$file = "XXXXXXXXXXXXXXX/$i.fa"; ####### input file of 4-fold degenerate sites
	$faExist = -e $file;
	if (!$faExist){
		next;
	}
		
	$seqio_object = Bio::SeqIO -> new(-file => $file);
	$seq_object = $seqio_object -> next_seq;
	$sequence = $seq_object -> seq;
	$len = length($sequence);
	$countSites = $countSites + $len;
	
	for ($j = 0; $j < $len; $j++){
       my %alphabet;
       $freq[$j] = \%alphabet;
       
       my %alphabetAth;
       $freqAth[$j] = \%alphabetAth;
       
       my %alphabetCru;
       $freqCru[$j] = \%alphabetCru;
  }
 
	
	$seqio_object = Bio::SeqIO -> new(-file => $file);		   
	while ($seq_object = $seqio_object -> next_seq){	    		  
    $sequence = $seq_object -> seq;
    $title = $seq_object -> display_id;
    $ini = substr($title,0,3);
    for ($j = 0; $j < $len; $j++) {
    	if ($ini eq "Ath" || $ini eq "COL" || $ini eq "MTE" || $ini eq "Cru"){    		
    	  $char = substr($sequence,$j,1);
        &accumu($char,$freq[$j]); 

        if ($ini eq "Ath" || $ini eq "COL"){
        	&accumu($char,$freqAth[$j]);
        }
        if ($ini eq "MTE" || $ini eq "Cru"){
        	&accumu($char,$freqCru[$j]);
        }
      }
    }
	} 
	
  for ($k = 0; $k < $len; $k ++){
	  @ks = keys %{$freq[$k]};
	  @vs = values %{$freq[$k]};
	  
        if ($vs[0] == $vs[1]){
            $char = $ks[1];
            $c1 = $freqAth[$k]{$char};
            $c2 = $freqCru[$k]{$char};
            if ($c1 eq ""){$c1 = 0};
            if ($c2 eq ""){$c2 = 0};

            $count[$c1][$c2] = $count[$c1][$c2]+0.5;
            $count[81-$c1][23-$c2] = $count[81-$c1][23-$c2] + 0.5
        }else{

    	      if ($vs[0] > $vs[1]){
    	        $char = $ks[1];		  	
    	      }else {
    	        $char = $ks[0];	
    	      }
    	  
    	      $c1 = $freqAth[$k]{$char};
    	      $c2 = $freqCru[$k]{$char};
    
    	      if ($c1 eq ""){$c1 = 0};
    	      if ($c2 eq ""){$c2 = 0};
    	 
    	      $count[$c1][$c2]++;
        }
	  	
  }
	

}





for ($b = 0; $b <= 23; $b ++){
	print OUT "d1_".$b."\t";
	for ($a = 0; $a < 81; $a ++){
		print OUT $count[$a][$b]."\t";
	}
	print OUT $count[81][$b]."\n";
}

for ($a = 0; $a <= 81; $a ++){
  for ($b = 0; $b <= 23; $b ++){
  	print OUT12 $count[$a][$b]."\t";
  }	
}


close OUT;
close OUT12;


sub accumu(){
    #my %hTable = %{$_[1]};
    my $hRef = $_[1];

    if (exists $hRef -> {$_[0]}){
       $$hRef{$_[0]} = $$hRef{$_[0]} + 1;
    }else{
       $$hRef{$_[0]} = 1;
    }
}
