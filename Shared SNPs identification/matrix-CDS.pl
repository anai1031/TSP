
## This script extracts the coding sequences of an orthologue of all the samples under study
## Author: Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science
## Email: qiongwu@ibcas.ac.cn

use Bio::SeqIO;
use Bio::Seq;
#use Bio::Tools::GFF;
use strict;

my $filename = $ARGV[0];

my $Ath_gff = "Athaliana_167_gene.gff3"; ## the GFF file for species 1
my $Cru_gff = "Crubella_183_gene.gff3"; ## the GFF file for species 2

my $Ath_ref = "Athaliana_167.fa"; ## the reference file for species 1
my $Cru_ref = "Crubella_183.fa"; ## the reference file for species 2


my $i = 0;
my $run = 1;
my @arrCDS;
my ($name, $key, $title, $line, $lineCDS, $line_next, $search, $seq);
my (@l_next, @l_CDS, @lin);
my ($seqio_object, $seq_object, $seq_ob);
our ($sequenceATH, $sequenceCRU, $sequenceR418, $seqCru);
my @s;
my $j;
my @sequenceCru;
my (@array, @current, @next, $currentLine, $nextLine, $true, $seqR418);


open(SQL,"<$filename") or die "Cannot open $filename to read.\n"; ## the file listing 16014 orthologues


while ($run){

 $line_next = <SQL>; 

  if (eof) {
    $run = 0;
        
  } else { 
    $i++; 
    my $CDSsequence = "";
    @l_next = split(/\s+/,$line_next);

     while (@l_next){
    
       $name = shift @l_next;
       $key = shift @l_next;
      
       if ($name eq "Ath"){ ##extract the coding sequences of this orthologue of all the A.thaliana samples from the downloaded SNP matrix
           
           $search = substr($key, length($key)-8, 8);
           open(ATHGFF,"<$Ath_gff") or die "Cannot open $Ath_gff to read.\n"; 
 
           while ($line = <ATHGFF>) {
             @lin = split(/\s+/,$line);
             if (($lin[8] =~/$search/) && ($lin[2] eq "CDS")) {
               push @arrCDS, $line;     
             } 
             
           }  #end of while <ATHgff>
 
           close ATHGFF;
           $sequenceATH = "";
  
           $lineCDS = $arrCDS[0];
           @l_CDS = split(/\s+/,$lineCDS);
            
           @s = &thaliana(); 
 
           if ($l_CDS[6] eq "-") {
               @arrCDS = reverse (@arrCDS);
           }
  
           while (@arrCDS){
                               
               $lineCDS = shift @arrCDS;
               @l_CDS = split(/\s+/,$lineCDS);
               print $lineCDS;
 
               $seqio_object = Bio::SeqIO -> new(-file => $Ath_ref);
               $seq_object = $seqio_object -> next_seq;
               $title = $seq_object -> display_id;              
 
               while ($l_CDS[0] ne $title) {
                 $seq_object = $seqio_object -> next_seq;
                 $title = $seq_object -> display_id;
                 print $title;
               }
               $seq = $seq_object -> subseq($l_CDS[3],$l_CDS[4]);
         
               $sequenceATH = $sequenceATH.$seq;              
               
           } #end of while @arrCDs
           
           
           if ($l_CDS[6] eq "-") {
              $seq_ob = Bio::Seq->new(-seq => $sequenceATH);
              $seq_ob = $seq_ob  -> revcom;
              $sequenceATH = $seq_ob -> seq;
           }                                     
       } #end of if Ath
 
       if ($name eq "Cru"){ ##extract the coding sequences of all the C.rubella samples from the called .VCF files
         $search = substr($key, length($key)-8, 8);
         open(CRUGFF,"<$Cru_gff") or die "Cannot open $Cru_gff to read.\n";

          while ($line = <CRUGFF>) {
            @lin = split(/\s+/,$line);

            if (($lin[8] =~/$search/) && ($lin[2] eq "CDS")) {
               push @arrCDS, $line;
            }
          }           
          close CRUGFF;
         
           $lineCDS = $arrCDS[0];

           @l_CDS = split(/\s+/,$lineCDS);
           if ($l_CDS[6] eq "-") {
             @arrCDS = reverse (@arrCDS);
           }
           
           our $sequenceCRU = "";
           
           for ($j = 0; $j < 22; $j++){
              $sequenceCru[$j] = "";
           }

           while (@arrCDS){

              $lineCDS = shift @arrCDS;
              @l_CDS = split(/\s+/,$lineCDS);
             
              print $lineCDS;
              
              $seqio_object = Bio::SeqIO -> new(-file => $Cru_ref);

              $seq_object = $seqio_object -> next_seq;

              $title = $seq_object -> display_id;

              while ($l_CDS[0] ne $title) {
                  $seq_object = $seqio_object -> next_seq;
                  $title = $seq_object -> display_id;
                  print $title;
              }
              $seq = $seq_object -> subseq($l_CDS[3],$l_CDS[4]);
              $sequenceCRU = $sequenceCRU.$seq;
               
              for ($j =1; $j < 23; $j++){
                 $seqCru = &checkVCF (".XXXXXXXXX/$j.vcf"); ### SNP information for a Cru sample saved in a .vcf file
                 $sequenceCru[$j-1] = $sequenceCru[$j-1].$seqCru;
              }                

           }
         if ($l_CDS[6] eq "-") {
            $seq_ob = Bio::Seq->new(-seq => $sequenceCRU);
            $seq_ob = $seq_ob  -> revcom;
            $sequenceCRU = $seq_ob -> seq;
            
            for ($j =1; $j < 23; $j++){
               $seq_ob = Bio::Seq->new(-seq => $sequenceCru[$j-1]);
               $seq_ob = $seq_ob  -> revcom;
               $sequenceCru[$j-1] = $seq_ob -> seq;
            }

         }

       } # end of if Cru

     }  # end of while l_next

  
  my $out_file = "XXXXXXXXXXX/".($i+$ARGV[1]).".max"; ## The output file for the coding sequences of all the samples on this orthologue
  open(OUT,">$out_file") or die "Cannot open $out_file to write.\n";  

  print OUT ">COL\n".$sequenceATH."\n";
  for ($j = 0; $j < 80; $j++){
      print OUT ">Ath_".($j+1)."\n".$s[$j]."\n";
  }
 
  print OUT ">MTE\n".$sequenceCRU."\n"; 
  for ($j = 1; $j < 23; $j++){
      print OUT ">Cru_".$j."\n".$sequenceCru[$j-1]."\n";
  }

  close OUT;
  
 } # end of else



} # end of while run

sub checkVCF (){ ## This function extracts the SNP information for each sample from the .VCF file;
    @array = @l_CDS;
    open(VCF,"<$_[0]") or die "Cannot open the file to read.\n";
    $sequenceR418 = "";
  
    $currentLine = <VCF>;
    @current = split(/\s+/,$currentLine);
   
    while ($current[0] ne $array[0]){
      $currentLine = <VCF>;
      @current = split(/\s+/,$currentLine);     
    }
    
    if (eof){
      return $seq_object -> subseq($array[3],$array[4]);
    } # no snp falls in this scaffold

    while ($current[1] < $array[3]){
      $currentLine = <VCF>;
      @current = split(/\s+/,$currentLine);

      if (eof || $current[0] ne $array[0]){
        return $seq_object -> subseq($array[3],$array[4]);
      }
    } # the snps falling in this scaffold does not lie in the range
    @current = split(/\s+/,$currentLine);
    if ($current[1] > $array[4]){
        return $seq_object -> subseq($array[3],$array[4]);
    } # out of range, too

    if ($current[1] > $array[3]){
      $sequenceR418 = $seq_object -> subseq($array[3],$current[1]-1);
    }
    
    $true = 1;
    while ($true){ # finally, appropriate SNPs!!!
      if ((length($current[4]) != 1) || (length($current[3]) != 1)){
         $sequenceR418 = $sequenceR418.$seq_object -> subseq($current[1],$current[1]);  
      } else {
         $sequenceR418 = $sequenceR418.$current[4];
      }
     
      $nextLine = <R418>;
      @next = split(/\s+/,$nextLine);
      if (eof || $next[0] ne $array[0] || $next[1]>$array[4]){
        if ($array[4] > $current[1]){
          $sequenceR418 = $sequenceR418.$seq_object -> subseq($current[1]+1,$array[4]);
        }
        return $sequenceR418;
      }# current is the last position
      if (($next[1] <= $array[4]) && ($next[1] > $current[1]+1)){
        $sequenceR418 = $sequenceR418.$seq_object -> subseq($current[1]+1,$next[1]-1);
      } 
      
      $currentLine = $nextLine;
      @current = split(/\s+/,$currentLine);

    }

   close VCF;
}
sub thaliana (){ ## This function extracts the CDS sequences from the SNP matrix of 80 A.thaliana samples
    my @arr = @arrCDS;    
    my @names;
    my $k;
    my $tmp;
    push @names, "Tair_1.txt";    #The original matrix file was separated into 5 files corresponding to the 5 chromosomes.
    push @names, "Tair_2.txt";
    push @names, "Tair_3.txt";
    push @names, "Tair_4.txt";
    push @names, "Tair_5.txt";

    if ($l_CDS[6] eq "-") {
       @arr = reverse (@arr);  
    }

    my $chrNo = substr($l_CDS[0],3,1);
    print "ChrNo: ".$chrNo."\n";
    open(TAIR,"<$names[$chrNo-1]") or die "Cannot open Tair_".$chrNo.".txt to read.\n";
    print "File: ".$names[$chrNo-1]."\n";
    
    for ($k = 0; $k < 80; $k++){
       $s[$k] = "";
    }

    while (@arr){
       my $lCDS = shift @arr;
       
       print $lCDS."\n";
   
       my @l_cds = split(/\s+/,$lCDS);
       my $start = $l_cds[3];
       my $end = $l_cds[4];

       my $l = <TAIR>;
       my @ls = split(/\s+/,$l);

       while ($ls[1] < $start - 1){
           $l = <TAIR>;
           @ls = split(/\s+/,$l);  
       }

       while ($ls[1] < $end){
           $l = <TAIR>;
           @ls = split(/\s+/,$l);
           
           for ($k = 0; $k < 80; $k++){
            $s[$k] = $s[$k].$ls[$k+3];
           }
       }
    }
    
    if ($l_CDS[6] eq "-") {
       for ($k = 0; $k < 80; $k++){
          $seq_ob = Bio::Seq->new(-seq => $s[$k],-alphabet => 'dna');
          $seq_ob = $seq_ob  -> revcom;
          $s[$k] = $seq_ob -> seq;
       }
    }
    close TAIR;
    return @s; 
}

close SQL;

