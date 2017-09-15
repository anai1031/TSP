
## This script extracts the orthologous 4-fold sites from the exiting aligned fasta files including Ath, Aly, Cru.
## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn

use Bio::SeqIO;
use Bio::Seq;
#use Bio::Tools::GFF;

for ($i = 1; $i < 16047; $i ++) {
	$file = "XXXXXXX/$i.max";	 # check whether there is any coding error in the eff files.
	$maxExist = -e $file;
	if (!$maxExist){
		next;
	}
	
  $seqio_object = Bio::SeqIO -> new(-file => $file);		   
	while ($seq_object = $seqio_object -> next_seq){	    		  
	  $title = $seq_object -> display_id;
	  if ($title eq "MN47"){
	    $lengthAly = $seq_object -> length;
	  }  
	  
	  if ($title eq "MTE"){
	    $lengthCru = $seq_object -> length;
	  }
	  
	  if ($title eq "COL"){
	    $lengthAth = $seq_object -> length;
	  }	   
	}
	
	if (($lengthAth % 3 != 0) || ($lengthCru % 3 != 0) || ($lengthAly % 3 != 0)){ ## check if any of the CDS region is of length not multiple of 3
		print $i.": Coding error!!\n";
		next;
	}
	
	$fileAthCru = "XXXXX/$i.fa"; ## Ath and Cru examples
	$fileAly = "XXXXX/$i.fa"; ## Aly examples

	
	$out = "XXXXXXXXX/$i.fa"; ## the output to save the 4-fold degenerate sites
	
	open(OUT,">$out") or die "Cannot open $out to read.\n";		
	
	$seqio_object = Bio::SeqIO -> new(-file => $fileAthCru);
	while ($seq_object = $seqio_object -> next_seq){	    		  
	  $title = $seq_object -> display_id;
	  if ($title eq "MN47"){
	    $seqAly = $seq_object -> seq;
	  }  
	  
	  if ($title eq "MTE"){
	    $seqCru = $seq_object -> seq;
	  }
	  
	  if ($title eq "COL"){
	    $seqAth = $seq_object -> seq;
	  }	   
	}
	@posAth = ();
	@posAthAly = ();
	@posAthAlyCru = ();
	@codonP = ();
	
	for ($j = 0; $j < length($seqAth); $j++){
		 $charAth = substr($seqAth,$j,1);
	   if ($charAth ne "-"){
	     	push @posAth, $j;
	   } 	   	
	}
	while (@posAth){
		$p = shift @posAth;
		$charAly = substr($seqAly,$p,1);
		if ($charAly ne "-"){
		  push @posAthAly, $p;	
	  }		   	
	} 
	while (@posAthAly){
		$p = shift @posAthAly;
		$charCru = substr($seqCru,$p,1);
		if ($charCru ne "-"){
		  push @posAthAlyCru, $p;	
	  }		   	
	} 
	
	while(@posAthAlyCru){
	  $p = 	shift @posAthAlyCru;	  
	  $pAth = $p;
	  $pAly = $p;
	  $pCru = $p;
	  for ($k = 0; $k < $p; $k ++){
	  	$charAth = substr($seqAth,$k,1);
	  	$charAly = substr($seqAly,$k,1);
	  	$charCru = substr($seqCru,$k,1);
	    if ($charAth eq "-"){
	    	$pAth--;
	    }
	    if ($charAly eq "-"){
	    	$pAly--;
	    }
	    if ($charCru eq "-"){
	    	$pCru--;
	    }
	    
	  }
	  
	  if ((($pAth % 3) == ($pAly % 3)) && (($pAth % 3) == ($pCru % 3))){
	    if (($pAth % 3) == 0){
	    	if ((substr($seqAth, $p+1, 1) ne "-")&&(substr($seqAth, $p+2, 1) ne "-")&&(substr($seqAly, $p+1, 1) ne "-")&&(substr($seqAly, $p+2, 1) ne "-")&&(substr($seqCru, $p+1, 1) ne "-")&&(substr($seqCru, $p+2, 1) ne "-")){
	    		push @codonP, $p;
	    		shift @posAthAlyCru;
	    		shift @posAthAlyCru;
	      }	    	
	    }		  	
	  }		
	}

	

     for ($b = 0; $b < @codonP; $b++){
       my %alphabet;
       $freq[$b] = \%alphabet;
     }

	
	$seqio_object = Bio::SeqIO -> new(-file => $fileAthCru); ## A.th & C.ru & MN47

	while ($seq_object = $seqio_object -> next_seq){
      $seq = $seq_object -> seq;  
      $ini = substr($title,0,3);
      if ($ini eq "Ath" || $ini eq "COL" || $ini eq "MTE" || $ini eq "Cru" || $ini eq "MN4"){ 
      }else {
      	next;
      }           
          
      for ($k = 0; $k < @codonP; $k ++){
      	$start = $codonP[$k];
    	  $str = substr($seq,$start,2);
      	$hashRef = $freq[$k];
      	
      	if (($str eq "GT")||($str eq "TC")||($str eq "CT")||($str eq "CC")||($str eq "CG")||($str eq "AC")||($str eq "GC")||($str eq "GG")){
      		&accumu($str,$hashRef);   		
      	}else {
      	  &accumu("XX",$hashRef);	
      	}
      	
      }
      
  } ## record the information of the first two codons for each amnio acid in each sequence;
  
  	$seqio_object = Bio::SeqIO -> new(-file => $fileAly); ## A.ly
	  while ($seq_object = $seqio_object -> next_seq){
      $seq = $seq_object -> seq;            
          
      for ($k = 0; $k < @codonP; $k ++){
      	$start = $codonP[$k];
    	  $str = substr($seq,$start,2);
      	$hashRef = $freq[$k];
      	
      	if (($str eq "GT")||($str eq "TC")||($str eq "CT")||($str eq "CC")||($str eq "CG")||($str eq "AC")||($str eq "GC")||($str eq "GG")){
      		&accumu($str,$hashRef);   		
      	}else {
      	  &accumu("XX",$hashRef);	
      	}
      }
      
  } ## record the information of the first two codons for each amnio acid in each sequence;
  
  
  
  @arrPrint = ();
  
  for ($k = 0; $k < @codonP; $k ++){
      	
        @keysCodon = keys %{$freq[$k]};  

        $keySize = @keysCodon;   #force the array into a scalar context
        
        if ($keySize == 1){
        	if ($keysCodon[0] ne "XX"){
        	  push @arrPrint, $codonP[$k];
        	}
          
        }

  }
       	      
  $seqio_object = Bio::SeqIO -> new(-file => $fileAthCru); 
  while ($seq_object = $seqio_object -> next_seq){
    	$seq = $seq_object -> seq;
    	$title = $seq_object -> display_id;   
    	                    
      $ini = substr($title,0,3);
      if ($ini eq "Ath" || $ini eq "COL" || $ini eq "MTE" || $ini eq "Cru" || $ini eq "MN4"){
      	print OUT ">".$title."\n";      	
      }else{
      	next;
      }
      
      for ($k = 0; $k < @arrPrint; $k ++){
      	$count = $arrPrint[$k];
      	$char = substr ($seq, $count+2, 1);
      	print OUT $char;
      }
      print OUT "\n";
	}
	
	$seqio_object = Bio::SeqIO -> new(-file => $fileAly);
  while ($seq_object = $seqio_object -> next_seq){
    	$seq = $seq_object -> seq;
    	$title = $seq_object -> display_id;   
      print OUT ">".$title."\n"; 
           	      
      for ($k = 0; $k < @arrPrint; $k ++){
      	$count = $arrPrint[$k];
      	$char = substr ($seq, $count+2, 1);
      	print OUT $char;
      }
      print OUT "\n";
	}
	
	
	close OUT;
}



sub accumu(){
    my $hRef = $_[1];

    if (exists $hRef -> {$_[0]}){
       $$hRef{$_[0]} = $$hRef{$_[0]} + 1;
    }else{
       $$hRef{$_[0]} = 1;
    }
}
