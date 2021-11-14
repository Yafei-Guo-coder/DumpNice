(base) [xuebo@lulab1 Durum]$ cat estimate.sh 
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1/analysis
OMP_NUM_THREADS=1 smc++ estimate 6.5e-9 /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1/data/chr*.smc.gz --polarization-error 0.5 \
-o /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1/analysis --cores 30 --spline pchip --timepoints 1e3 1e5 --regularization-penalty 6.0 --knots 10

(base) [xuebo@lulab1 vcf]$ cat header.sh 
#!/bin/bash
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
sed -i '15 r /data1/home/xuebo/Projects/Speciation/smcpp/sample15/group/vcf_contifheader.txt' chr$i.vcf 
bgzip -c chr$i.vcf >  chr$i.vcf.gz
tabix chr$i.vcf.gz
done
zcat chr10.vcf.gz | less -LS

(base) [xuebo@lulab1 data]$ ls
chr1.smc.gz   chr14.smc.gz  chr19.smc.gz  chr21.smc.gz  chr26.smc.gz  chr3.smc.gz   chr33.smc.gz  chr38.smc.gz  chr40.smc.gz  chr9.smc.gz
chr10.smc.gz  chr15.smc.gz  chr2.smc.gz   chr22.smc.gz  chr27.smc.gz  chr31.smc.gz  chr34.smc.gz  chr39.smc.gz  chr7.smc.gz
chr13.smc.gz  chr16.smc.gz  chr20.smc.gz  chr25.smc.gz  chr28.smc.gz  chr32.smc.gz  chr37.smc.gz  chr4.smc.gz   chr8.smc.gz

(base) [xuebo@lulab1 data]$ zcat chr10.smc.gz  | head
# SMC++ {"version": "1.15.4.dev2+g1823009.d20200510", "pids": ["Durum"], "undist": [[["B118_B122", 0], ["B118_B122", 1], ["B119_XI_S5", 0], ["B119_XI_S5", 1], ["B116_B120", 0], ["B116_B120", 1], ["B115_B124", 0], ["B115_B124", 1], ["B114_B124", 0], ["B114_B124", 1], ["B119_B125", 0], ["B119_B125", 1], ["B122_XI_S2", 0], ["B122_XI_S2", 1], ["B120_B122", 0], ["B120_B122", 1], ["B119_B120", 0], ["B119_B120", 1], ["B115_XI_S5", 0], ["B115_XI_S5", 1], ["B117_B122", 0], ["B117_B122", 1], ["B114_B117", 0], ["B114_B117", 1], ["B113_B115", 0], ["B113_B115", 1], ["B113_B125", 0], ["B113_B125", 1]]], "dist": [[["XI_S1_XI_S4", 0], ["XI_S1_XI_S4", 1]]]}
1 0 0 28
1270 -1 0 0
240 0 0 28
1 -1 0 28
77 0 0 28
177 -1 0 0
319 0 0 28
2376 -1 0 0
19 0 0 28

(base) [xuebo@lulab1 V1]$ cat group_Durum_header.txt 
XI_S1_XI_S4,B118_B122,B119_XI_S5,B116_B120,B115_B124,B114_B124,B119_B125,B122_XI_S2,B120_B122,B119_B120,B115_XI_S5,B117_B122,B114_B117,B113_B115,B113_B125


(base) [xuebo@lulab1 Durum]$ cat getsmcformat.sh 
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1
cd /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1
sed 's/\t/_/g' /data1/home/xuebo/Projects/Speciation/smcpp/sample15/group/Durum_15.txt |  tr "\n" "," | sed s'/.$//' > group_Durum_header.txt &
  mkdir data
name=$(sed -n 1p group_Durum_header.txt)
echo $name
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/vcf/chr$i.vcf.gz \
/data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/V1/data/chr$i.smc.gz $i Durum:$name  &
  done


(base) [xuebo@lulab1 Durum]$ cat getdiploidvcf.sh 
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/vcf
cd /data1/home/xuebo/Projects/Speciation/smcpp/sample15/Durum/vcf
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr$i.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/sample15/group/Durum_15.txt --out chr$i.vcf &
  done


(base) [xuebo@lulab1 mask1]$ cat getMerged.sh 
for i in 21 30
do
zcat /data1/home/yaozhou/Projects/EVO/data/merge/lineage/chr/chr$i.masked.bed.gz > chr$i.bed
cat chr$i.bed /data1/home/yaozhou/data/ref/wheat/gff/chr$i.masked.txt /data2/yaozhou/Projects/reference/chr$i.masked.txt > chr$i.masked.txt
sort -k2n -k3n chr$i.masked.txt > chr$i.sorted.txt
bedtools merge -i chr$i.sorted.txt -d 100 > chr$i.masked.bed &
  done


(base) [xuebo@lulab1 E3removeExon]$ cat getvcf.sh 
#!/bin/bash
for chr in {1..42}
do
WGS --model vcf --type unintersect --file /data1/home/xuebo/Projects/Speciation/smcpp/E3unmasked/chr${chr}.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/gff/chr${chr}.exon --out /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr${chr}.vcf &
  done


(base) [xuebo@lulab1 gff]$ cat getAgff.sh 
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
grep -P "^$i" GeneLulabNoUn1_1.gff3 >> A.gff3
done


for i in {1,2,3,4,7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28,31,32,33,34,37,38,39,40}
do
grep -P "^$i" GeneLulabNoUn1_1.gff3 >> AB.gff3
done
(base) [xuebo@lulab1 gff]$ cat getExon.sh 
for i in {1..42}
do
WGS --model gff3 --chr $i --type gff2exon --file GeneLulabNoUn1_1.gff3 --out chr$i.exon &
  done
(base) [xuebo@lulab1 gff]$ cat getGenic.sh 
for i in {1..42}
do
WGS --model gff3 --chr $i --type gff2genic --file GeneLulabNoUn1_1.gff3 --out chr$i.genic &
  done


(base) [xuebo@lulab1 gff]$ cat get4D.pl
#!/usr/bin/perl
use warnings ;
use strict ;
=begin
This script takes a gff3 (or gtf or gff) file and a fasta file as input. 
It then identifed all four-fold degenerate positions and retains those.

=cut

### file information is passed via command line
my $gff = $ARGV[0] ; 
my $fa = $ARGV[1] ; 



### read in set of 4d degenerate codons
### this file can take any format as long as all it contains are \t and 4-fold degenerate codons in the second-nth columns on each line
my %codon ; 
open CODON, "<codon_table.txt" ;
while(<CODON>) { 
  chomp $_ ;
  my @split = split ( /\t/, $_ ) ; 
  foreach ( 1..$#split ) { 
              $codon{$split[$_]} ++ ;
}
}
close CODON ; 
### gene information is extracted from standard gff3 files.
### this also shoudl work for gtf and gff files
my %gene ; 
open GFF, "<$ARGV[0]" ;
while (<GFF>) { 
  chomp $_ ; 
  if($_ =~ /^\s*#/) { next; }
     my @split = split ( /\t/, $_ ) ; 	
     ### parse information for mRNA molecule
     if ( $split[2] eq "CDS" ) { 
       my $id = $split[0]."_".$split[3]."_".$split[4]."_".$split[7] ;
       $gene{$id}{"strand"} = $split[6] ;
       $gene{$id}{"chrom"} = $split[0] ; 
       $gene{$id}{"start"} = $split[3] - 1 ;
       $gene{$id}{"stop"} = $split[4] - 1 ; 
       $gene{$id}{"frame"} = $split[7] ;
     }
}
close GFF ; 

### sequence data 
my $chrom = "" ;
my $chr = "start" ;
open FA, "<$fa" ;
while (<FA>) { 
  chomp $_ ; 
  if ( $_ =~ m/^>(.+)/ || eof() ) {
    my $new_chr = $1 ; 
    
    if ( $chr ne "start" ) {
      
      my %syn ;                               ### this hash will store the positions of all 4-d sites on the chromosome
      my %non ;                               ### this hash will store all non-4d sites
      foreach my $mrna ( keys %gene ) {
        if ( $gene{$mrna}{"chrom"} eq $chr ) {
          
          #### extract the sequences first
          my $protein = "" ;
          my @sites ; 
          if ( $gene{$mrna}{"strand"} eq '+' ) { 
            $protein = substr( $chrom, $gene{$mrna}{"start"}, $gene{$mrna}{"stop"} - $gene{$mrna}{"start"} + 1 ) ; 
            foreach ( $gene{$mrna}{"start"}..$gene{$mrna}{"stop"} ) { 
              push @sites, $_ ;
            }
          }
          else {
            ### if on - strand, we need to reverse compliment this sequence
            $protein .= reverse_complement ( substr( $chrom, $gene{$mrna}{"start"}, $gene{$mrna}{"stop"} - $gene{$mrna}{"start"} + 1 ) ) ;
            foreach ( $gene{$mrna}{"start"}..$gene{$mrna}{"stop"} ) { 
              push @sites, $_ ;
            }
            @sites = reverse( @sites ) ; ### reverse list of sites to account for - strand
          }
          
          ### identify four fold degenerate sites and non-4d sites
          for ( my $i = $gene{$mrna}{"frame"} ; $i < length( $protein ) - 2; $i += 3 ) { 
            if ( exists( $codon{ substr( $protein, $i, 3 ) } ) ) { 
              $syn{$chr}{$sites[$i+2]} ++ ;           ### only the third base can be 4d
            }
            else {
              $non{$chr}{$sites[$i+2]} ++ ;
            }
            $non{$chr}{$sites[$i]} ++ ;                 ### first and second base are always not 4d
            $non{$chr}{$sites[$i+1]} ++ ; 
          }
        }
      }
      
      #### foreach 4d site, check to make sure it's 4d in all transcripts
      #### overlapping reading frames could cause problems here. Hence, whole chromosomes at once
      foreach my $ch ( sort ( keys %syn ) ) { 
        foreach my $pos ( sort { $a <=> $b } ( keys %{$syn{$ch}} ) ) { 
          if ( !exists( $non{$ch}{$pos} ) ) { 
            print $ch, "\t", $pos + 1, "\n" ;
          }
        }
      }
      
      ### reset
      $chrom = "" ;
      $chr = $new_chr ; 
    }
    else { 
      $chr = $new_chr ; 
    }
  }
  
  ### put everything in upper case
  else { 
    $_ =~ s/a/A/g ;
    $_ =~ s/t/T/g ;
    $_ =~ s/c/C/g ;
    $_ =~ s/g/G/g ;
    $chrom .= $_ ; 
  }
}

#### reverse compliment for - strand sequences
sub reverse_complement {
  my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}



