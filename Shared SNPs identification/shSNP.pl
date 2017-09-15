
## This script looks for the orthologous shared bi-allelic SNPs between A.th and C.ru, requiring both SNP with MAF > 0.05

## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn

use Array::Utils qw(:all);

my ($line_next, $numSites, $remainder, $quotient, $j, $i, $k, $m, $n, $curr10, $pos, $species, $ini, $arrSize, $tf);
my (@l_next, @taxaArray, @freq, @keysAth, @keysAly, @keysCru, @keysCor, @keysCgr, @diff);
my %hashTable;

my $out_file = "XXXXXX"; ### The output file
open(OUT,">$out_file") or die "Cannot open $out_file to write.\n";
print OUT "Number\t File\t Position\t Distribution\n";


our $count = 0;

for ($i = 1; $i < 16048; $i ++){ ## process each orthologue
    my $in = "XXXXXXX/$i.phyi.phy"; ### The aligned sequence output from MUSCLE
    eval{
        open(ALN,"<$in") or die "Cannot open $in to read.\n";
    };
    if ($@){
        print "An error occurred ($@), continuing\n";
        next;
    }   
    my $filename = $in;
    open(PHY,"<$filename") or die "Cannot open $filename to read.\n";
    
    
    $line_next = <PHY>;
    @l_next = split(/\s+/,$line_next);

    $numSites = $l_next[1];
    $remainder = ($numSites+10)%60;    ## The first block of a .phyi.phy file has 50 sites; the last block has $remainder sites; all the others have 60.
    $quotient =($numSites+10-$remainder)/60;

    for ($c = 0; $c < 2; $c ++){
        for ($d = 0; $d < 60; $d ++){
            my %alphabet;
            $freq[$c][$d] = \%alphabet;
        }
    }
    for ($j = 0; $j < 104; $j ++){
        $line_next = <PHY>;
        @l_next = split(/\s+/,$line_next);
        
        $taxaArray[$j] = $l_next[0];
        #print $taxaArray[$j]."\n";
        $ini = substr($l_next[0],0,3);
    
        if ($ini eq "Ath" || $ini eq "COL"){
          $species = 0;
        } else{
        	  if ($ini eq "MTE" || $ini eq "Cru") {
              $species = 1;
            }         
        }
       
        if ($quotient > 0){
            for ($ii = 0; $ii < 5; $ii ++){
                $curr10 = $l_next[$ii+1];
                for ($a = 0; $a < 10; $a ++){
                     $char = substr($curr10,$a,1);
                     $pos = $ii*10 + $a;
                     
                     $hashRef = $freq[$species][$pos];
    
                     if ($char eq "A"){
                         &accumu("A",$hashRef);
                     }
                     if ($char eq "T"){
                         &accumu("T",$hashRef);
                     }
                     if ($char eq "C"){
                         &accumu("C",$hashRef);
                     }
                     if ($char eq "G"){
                         &accumu("G",$hashRef);
                     }
    
                }
            }
        }
        
     
    }

    for ($j = 0; $j < 50; $j ++){
        @keysAth = keys %{$freq[0][$j]};
        @keysCru = keys %{$freq[1][$j]};
        @valuesAth = values %{$freq[0][$j]};  
        @valuesCru = values %{$freq[1][$j]};     

        $arrSize = @keysAth;   
        
        $tf = 0;
        if($arrSize == 2){ ## the SNP in A.th is bi-allelic 
          if (($valuesAth[0]>4) && ($valuesAth[1]>4) && $valuesCru[0]>1 && $valuesCru[1]>1){ ## the SNP in both species is with MAF>0.05;
          	
          	@diff = &array_diff(\@keysAth, \@keysCru); ## The SNPs in both species share the same allele pair
            if ($#diff == -1){
            	  if (!$tf){
            	  	$count ++;
            	  }
            	  $tf =1;
                print OUT $count."\t";
            	  print OUT $i."\t";
                print OUT ($j+1)."\t";
                print OUT "Ath: ";
                while (($key, $value) = each %{$freq[0][$j]}){
                    print OUT "$key => $value, ";
                }
             
                print OUT "Cru: ";
                while (($key, $value) = each %{$freq[2][$j]}){
                    print OUT "$key => $value, ";
                }
                print OUT "\n";
            }
          }
        }
    }
    
    
    for ($c = 0; $c < 2; $c ++){
        for ($d = 0; $d < 50; $d ++){
            %{$freq[$c][$d]} = ();
        }
    
    }  
    
    ##################################################################
    for ($k = 0; $k < $quotient-1; $k ++){
        $line_next = <PHY>;
        for ($n = 0; $n < 104; $n ++){
            $line_next = <PHY>;
            @l_next = split(/\s+/,$line_next);
    
            $ini = substr($taxaArray[$n],0,3);
            if ($ini eq "Ath" || $ini eq "COL"){
                $species = 0
            } else{
            	if ($ini eq "MTE" || $ini eq "Cru") {
                 $species = 2;
              }
            }
        
            for ($ii = 0; $ii < 6; $ii ++){
                $curr10 = $l_next[$ii];
                for ($a = 0; $a < 10; $a ++){
                     $char = substr($curr10,$a,1);
                     $pos = $ii*10 + $a;
    
                     $hashRef = $freq[$species][$pos];
    
                     if ($char eq "A"){
                         &accumu("A",$hashRef);
                     }
                     if ($char eq "T"){
                         &accumu("T",$hashRef);
                     }
                     if ($char eq "C"){
                         &accumu("C",$hashRef);
                     }
                     if ($char eq "G"){
                         &accumu("G",$hashRef);
                     }
    
                }
    
    
            }
    
        }
    
        for ($j = 0; $j < 60; $j ++){
            @keysAth = keys %{$freq[0][$j]};
            @keysCru = keys %{$freq[1][$j]};
            @valuesAth = values %{$freq[0][$j]};
            @valuesCru = values %{$freq[1][$j]};
            
            $arrSize = @keysAth;
            
            $tf = 0;
    
            if($arrSize == 2){

              if (($valuesAth[0]>4) && ($valuesAth[1]>4) && $valuesCru[0]>1 && $valuesCru[1]>1){
              	
              	@diff = &array_diff(\@keysAth, \@keysCru);
                if ($#diff == -1){
                	  if (!$tf){
            	  	      $count ++;
            	      }
            	      $tf =1;
                    print OUT $count."\t";
            	      print OUT $i."\t";
                    print OUT (60*$k +50 +$j +1)."\t";
                    print OUT "Ath: ";
                    while (($key, $value) = each %{$freq[0][$j]}){
                        print OUT "$key => $value, ";
                    }
                    
                    print OUT "Cru: ";
                    while (($key, $value) = each %{$freq[2][$j]}){
                        print OUT "$key => $value, ";
                    }
                    print OUT "\n";                            
                }
              }
            }
    
        }
     
        for ($c = 0; $c < 2; $c ++){
            for ($d = 0; $d < 60; $d ++){
                %{$freq[$c][$d]} = ();
            }
    
        }
    }
     ##################################################################   
    if ($remainder > 0){
        $line_next = <PHY>;
       
        $tail = $remainder % 10;
        $num = ($remainder - $tail)/10;
    
        for ($n = 0; $n < 104; $n ++){
            $line_next = <PHY>;
            @l_next = split(/\s+/,$line_next);
            $ini = substr($taxaArray[$n],0,3);
            if ($ini eq "Ath" || $ini eq "COL"){
                $species = 0
            } else{
            	if ($ini eq "MTE" || $ini eq "Cru") {
                $species = 1;
              }               
            }
            for ($ii = 0; $ii < $num; $ii ++){
                $curr10 = $l_next[$ii];
                for ($a = 0; $a < 10; $a ++){
                     $char = substr($curr10,$a,1);
                     $pos = $ii*10 + $a;
    
                     $hashRef = $freq[$species][$pos];
    
                     if ($char eq "A"){
                         &accumu("A",$hashRef);
                     }
                     if ($char eq "T"){
                         &accumu("T",$hashRef);
                     }
                     if ($char eq "C"){
                         &accumu("C",$hashRef);
                     }
                     if ($char eq "G"){
                         &accumu("G",$hashRef);
                     }
                }
    
            }
            
            $curr10 = $l_next[$num]; # the rest chars, not necessarily 10;
    
            for ($ii = 0; $ii < $tail; $ii ++){
                $char = substr($curr10,$ii,1);
                $pos = $num*10 + $ii;
                $hashRef = $freq[$species][$pos];
    
                if ($char eq "A"){
                    &accumu("A",$hashRef);
                }
                if ($char eq "T"){
                    &accumu("T",$hashRef);
                }
                if ($char eq "C"){
                    &accumu("C",$hashRef);
                }
                if ($char eq "G"){
                    &accumu("G",$hashRef);
                }              
            }
        }
    
        for ($j = 0; $j < $remainder; $j ++){
            @keysAth = keys %{$freq[0][$j]};
            @keysCru = keys %{$freq[1][$j]};
            @valuesAth = values %{$freq[0][$j]};
            @valuesCru = values %{$freq[1][$j]};

            $arrSize = @keysAth;
            $tf = 0;
    
            if($arrSize == 2){
              if (($valuesAth[0]>4) && ($valuesAth[1]>4) && $valuesCru[0]>1 && $valuesCru[1]>1){
              	
              	@diff = &array_diff(\@keysAth, \@keysCru);
                if ($#diff == -1){
                    if (!$tf){
            	  	      $count ++;
            	      }
            	      $tf =1;
                    print OUT $count."\t";
            	      print OUT $i."\t";
                    print OUT (60*$quotient -10 +$j +1)."\t";
                    print OUT "Ath: ";
                    while (($key, $value) = each %{$freq[0][$j]}){
                        print OUT "$key => $value, ";
                    }
                   
                    print OUT "Cru: ";
                    while (($key, $value) = each %{$freq[2][$j]}){
                        print OUT "$key => $value, ";
                    }
                    print OUT "\n";

                }
              }
            }
    
        }  
      
    } # end of remainder;

} # end of 16046 loop

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




