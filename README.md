# Long-term balancing selection contributes to adaptation in Arabidopsis and its relatives

All the scripts are written in Perl. To run them, you have to install Perl, which can be downloaded from https://www.perl.org/, and make sure all relevant modules are installed as well. Meantime, in each script, the specific settings need to be modified accordingly, such as the path to the relevant files etc.

The scripts are authored by Qiong Wu, State Key Laboratory of Systematic and Evolutionary Botany, Institute of Botany, Chinese Academy of Science. They are functionally organized in directories as described below.

Shared SNPs identification
---
*	**matrix-CDS.pl**: extract the coding sequences of an orthologue of all the samples under study; similar work can be done for sequences in genic region. The Ath SNPs are in a matrix and the Cru SNPs are in .vcf files. The reference genome and the GFF annotation files have to be specified for each species in the script. <br> 

**Usage**: `perl matrix-CDS.pl <sql file> <Start file number>`<br>
**Required**: <br>
`<sql file>: this is the output file from InParanoid v2.0, with each line listing one specific orthologue for each species, based on their respective GFF annotations. Example lines are listed below:`<br>
  
    Aly	871508|PACid:16052529	Ath	AT3G02260.1|PACid:19661360	Cru	Carubv10016453m|PACid:20899174<br>
    Aly	489168|PACid:16049759	Ath	AT5G23110.1|PACid:19671139	Cru	Carubv10000018m|PACid:20908653<br>

`<Start file number>: Since there are around 16,000 orthologues, to run in parallel, the sql file was separated into small files in order, containing, e.g., 1-1000,1001-2000, 2001-3000, etc. orthologues. The start number has to be applied to keep consistent with the order in the whole list.` <br>
**Dependency**: `Before applying this script, you need to run InParanoid to get the orthologue information of different species.`<br>

-
*	**shSNP.pl**: look for the orthologous shared bi-allelic SNPs between A.th and C.ru, requiring both SNPs with MAF > 0.05; This script automatically scans the aligned sequence (output from MUSCLE v3.8.31) for each orthologue. The output of this procedure is a list of all the shared SNPs and their respective allele frequency in each species.<br>

**Usage**: `perl shSNP.pl`<br>
**Dependency**: `Before applying this script, you need to have all the sequences for each orthologue aligned using PHYLIP.`<br>

Demographic inference
****
*	**extract-4fold-tsp.pl**: automatically extract the orthologous 4-fold sites from the exiting aligned fasta files. The output of this procedure is all the matrices containing the 4-fold degenerate sites of each orthologue.<br>

**Usage**: `perl extract-4fold-tsp.pl`<br>
**Dependency**: `Before applying this script, you need to have the orthologous sequences from different species aligned using PHYLIP and converted into fasta format.`<br>

---
*	**computeSFS.pl**: automatically computes the Site Frequency Spectrum (MAF) of the two populations from all of the 4-fold degenerate sites. The output of this procedure is the SFS to feed into the fastsimcoal software, for both the joint and Multi versions. <br>

**Usage**: `perl computeSFS.pl`<br>
**Dependency**: `Before applying this script, you need to have all the matrices of the 4-fold degenerate sites for each orthologue.`<br>

Tree analysis
---
*	**allelicTree.pl**: automatically scan all the 100 bp windows covering each shared SNP in the genic region of the candidate genes and check whether all the samples of both Ath and Cru in the given tree can be clustered into a separate group, respectively. If not, this tree is not a species tree and can be an allelic tree, then further check will be done to see if there are 2 or more shared SNPs in such window to meet the searching criterion described in the paper. All the qualifying trees will be checked manually to confirm to be an allelic tree. The output is a file recording the information of each qualifying tree.<br>

**Usage**: `perl allelicTree.pl`<br>
**Dependency**: `Before applying this script, you need to build phylogenetic trees for all the 100 bp windows covering each shared SNP in the genic region of the candidate genes and the trees need to be in Newick format.`<br>

Simulation analysis
---
*	**simuPi.pl**: calculates Pi for every simulated 100 bp segment for Ath and Cru separately.<br>

**Usage**: `perl simuPi.pl`<br>

*	**simuMAF.pl**: calculates the MAF for all the SNPs appearing in the simulated segments for Ath and Cru separately.<br>

**Usage**: `perl simuPi.pl`<br>

**Dependency**: `Before applying these two scripts, we generated 1,000,000 simulated 100 bp neutral sequences for all the samples under the inferred demographic model using fastsimcoal2.`<br>
