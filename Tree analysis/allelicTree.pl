## automatically scan all the 100 bp windows covering each shared SNP in the genic region of the candidate genes and check whether all the samples of both Ath and Cru in the given tree can be clustered into a separate group, respectively. If not, this tree is not a species tree and can be an allelic tree, then further check will be done to see if there are 2 or more shared SNPs in such window to meet the searching criterion described in the paper. All the qualifying trees will be checked manually to confirm to be an allelic tree.
## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn


use Array::Utils qw(:all);
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;

my $i;
my ($line, $next);
my (@splits, @arr, @arrTmp);
my $count = 0;

my $out = "XXXXXXXXXXXXX/treeResult-shSNP2";
open(OUT,">$out") or die "Cannot open $out to read.\n";

@ath = (COL);
@rubella = (MTE);
@all = ();

for ($a = 1; $a < 81; $a++){
  push @ath, "Ath_$a";
}

for ($b = 1; $b < 23; $b++){
  push @rubella, "Cru_$b";
}   

push @all, @ath;
push @all, @rubella; 

for ($i = 1; $i < 976; $i ++){ ## search among the 975 shared SNPs in the genic region of the 433 candidate genes
	#print $i."\n";
	
	for ($j = 1; $j < 101; $j ++){ ## search among all the windows covering each site;
		$tsp = 0;
		$shTSP = 0;
    my $in = "XXXXXXXXXXXXXXXXXXXX/$i/$j.tree.nwk";
    $treeExist = -e $in;
    if ($treeExist){
      open(TREE,"<$in") or die "Cannot open $in to read.\n";
      $line = "";    
      while ($next = <TREE>){
          chomp ($next);
          $line = $line.$next;
      }
      
      for ($g = 0; $g < length($line); $g ++){
          $char = substr($line, $g, 1);
          if ($char eq ':'){
              substr($line, $g, 1) = "";
              while ((substr($line, $g, 1) ne ")") && (substr($line, $g, 1) ne ",")){
                  substr($line, $g, 1) = "";
              }
          }  
      }
      
      @a = ();
      @splits = ();
      $c = 0;
      for ($k = 0; $k < length($line); $k ++){
          @arrTmp = (); # Be careful!!!! This is the way to define an empty array. @array = "" defines an array which hasone element ---- "";
          $lineTmp = "";
          $char = substr($line, $k, 1);
          if ($char eq "("){
              push @arr, "(";
          }
  
          if ($char eq ")"){
              $ele = pop @arr;
              while ($ele ne "("){
                  push @arrTmp, $ele;
                  $ele = pop @arr;
              }
              $splits[$c] = [@arrTmp];             
              
              foreach $e (@arrTmp){
                  push @arr, $e;
              }
              $c++;
          
          }
  
         if ($char =~ /[a-zA-Z]/){
             while (($char ne ",") && ($char ne ")")){
                 $lineTmp = $lineTmp.$char;
                 $k ++;
                 $char = substr($line, $k, 1);
             }
             push @arr, $lineTmp;
             if ($char eq ")"){
                 $k --;
             }
         }
      }
      pop @splits;  

      if ((!(&isIn(\@ath)))&&(!(&isIn(\@rubella)))){ ## This tree is not a species tree; then check whether there are 2 or more shared SNP in this region
      	 $fa = "XXXXXXXXXXXXXXXXXXXXXXX/$i/$j.fa"; ##### The aligned sequence if fasta format
      	 $fileExist = -e $fa;
         if (!$fileExist){
         	 next;
         }
      	 for ($a = 0; $a < 2; $a ++){
           for ($b = 0; $b < 100; $b++){
             my %alphabet;
             $freq[$a][$b] = \%alphabet;
           }
         }
         $seqio_object = Bio::SeqIO -> new(-file => $fa);      
         while ($seq_object = $seqio_object -> next_seq){
         	
         	 $title = $seq_object -> display_id;
           $seq = $seq_object -> seq;  
               
           $ini = substr($title,0,3);
           
           if ($ini eq "Ath" || $ini eq "COL"){
             $species = 0;
           }
           if ($ini eq "MTE" || $ini eq "Cru") {
             $species = 1;
           }
               
           for ($m = 0; $m < 100; $m++){
         	  
         	  $char = substr($seq,$m,1);
           	$hashRef = $freq[$species][$m];
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
         
         for ($k = 0; $k < 100; $k ++){
      	
           @keysAth = keys %{$freq[0][$k]};
           @keysCru = keys %{$freq[1][$k]};
         
           @valuesAth = values %{$freq[0][$k]};  
           @valuesCru = values %{$freq[1][$k]};     
 
           $athSize = @keysAth;   #force the array into a scalar context
           $cruSize = @keysCru;
         
           if(($athSize == 2)&&($cruSize == 2)){   # both of size 2
         	  $shTSP++;
         	  @diff = &array_diff(\@keysAth, \@keysCru);
             if ($#diff == -1){
             	  $tsp ++;
             }
           }       
          
         }
         
         if ($tsp>1){ ## Yes, the region surrounding the shared (2 or more) SNPs clusters by allele, then print out the tree.
           $count ++;
           print OUT $count.": tree $i-$j\n$line\n\n";
         }
      } 
      
      close TREE;    
   }    
  }
}

close OUT;

sub accumu(){
    #my %hTable = %{$_[1]};
    my $hRef = $_[1];

    if (exists $hRef -> {$_[0]}){
       $$hRef{$_[0]} = $$hRef{$_[0]} + 1;
    }else{
       $$hRef{$_[0]} = 1;
    }
}

 
sub isIn {
    @test = @{$_[0]};
    $TorF = 0;
   
    for ($n = 0; $n <= $#splits; $n++){
        @tmp = @{$splits[$n]};
        @diff = &array_diff(\@test, \@tmp);
        if ($#diff == -1){
            $TorF = 1;
            last;
        }else {
            @minus = &array_minus(\@all, \@test);
            @diff = &array_diff(\@minus, \@tmp);
            if ($#diff == -1){
                $TorF = 1;
                last;
            }
        }

    }
    return $TorF;   
}





