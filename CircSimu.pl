use 5.012;
use Getopt::Long;
my ($fq1, $fq2, $out, $gtf, $coverage, $coverage2, $rand_mode, $rand_mode2, $read_length, $seq_err, $insert_length, $insert_length2, $perc_minor, $sigma, $sigma2, $ref_dir, $help, $if_chr, $exon_skipping, $psi, $circ_read);
Getopt::Long::GetOptions(
	#'1=s'	=>	\$fq1,
	#'2=s'	=>	\$fq2,
	'O=s'	=>	\$out,
	'G=s'	=>	\$gtf,
	'C=i'	=>	\$coverage,
	'LC=i'	=>	\$coverage2,
	'R=i'	=>	\$rand_mode,
	'LR=i'	=>	\$rand_mode2,
	'L=i'	=>	\$read_length,
	'E=i'	=>	\$seq_err,
	#'I=i'	=>	\$insert_length,
	'D=s'	=>	\$ref_dir,
	'CHR1=i'	=>	\$if_chr,
	'M=i'	=>	\$insert_length,
	'M2=i'	=>	\$insert_length2,
	'PM=i'	=>	\$perc_minor,
	'S=i'	=>	\$sigma,
	'S2=i'	=>	\$sigma2,
	'SE=i'	=>	\$exon_skipping,
	'PSI=i'	=>	\$psi,
	'H!'	=>	\$help
);

srand(5);
my @die_reason;
my $if_die;
#my ($fq1, $fq2, $gtf) = (">>./simulate_80_10X_80bp_200_1.fq",">>./simulate_80_10X_80bp_200_2.fq",">>./gtf_80_10X_80bp_200.out");
#my $coverage = 10;
#my $coverage2 = 100;
#my $rand_mode = 1;
#my $rand_mode2 = 2;
#my $read_length = 80;
my $cRNA_size = 0;
#my $if_PE = 2;
#my $seq_err = 1;
#my $insert_length = 200;
#my $ref_dir = "/panfs/home/zhao/gaoyuan/bwaphage/hg19/";

if(defined($help)){
	print "This is CIRI_AS_simulator, a simulation tool for circRNAs. Welcome!\n\n";
	print "Written by Yuan Gao. Any questions please mail to gaoyuan06\@mails.ucas.ac.cn.\n\n";
	print "Arguments (all required):\n";
	#print "\t-1\t\toutput simulated PE reads file 1 name\n";
	#print "\t-2\t\toutput simulated PE reads file 2 name\n";
	print "\t-O\t\tprefix of output files\n";
	print "\t-G\t\tinput gtf formatted annotation file name\n";
	print "\t-C\t\tset coverage or max coverage (when choosing -R 2) for circRNAs\n";
	print "\t-LC\t\tset coverage or max coverage (when choosing -LR 2) for linear transcripts\n";
	print "\t-R\t\tset random mode for circRNAs: 1 for constant coverage; 2 for random coverage\n";
	print "\t-LR\t\tset random mode for linear transcripts: 1 for constant coverage; 2 for random coverage\n";
	print "\t-L\t\tread length(/bp) of simulated reads (e.g. 100)\n";
	print "\t-E\t\tpercentage of sequencing error (e.g. 2)\n";
	#print "\t-I\t\tinsertion length (should be larger than read length) (e.g. 350)\n";
	print "\t-D\t\tdirectory of reference sequence(s) (please make sure all references referred in gtf file are included in the directory)\n";
	print "\t-CHR1\t\tif only choose chr1 to simulate sequencing reads: 1 for yes; 0 for no\n";
	print "\t-M\t\taverage(mu/bp) of insert length (major normal distribution) (e.g. 320)\n";
	print "\t-M2\t\taverage(mu/bp) of insert length (minor normal distribution) (e.g. 550)\n";
	print "\t-PM\t\tpercentage of minor normal distribution in total distribution (e.g. 10; 0 for no minor distribution)\n";
	print "\t-S\t\tstandard deviation(sigma/bp) of insert length (e.g. 70)\n";
	print "\t-S2\t\tstandard deviation(sigma/bp) of insert length (e.g. 70)\n";
	print "\t-SE\t\twhether simulate exon skipping: 1 for yes; 0 for no\n";
	print "\t-PSI\t\tpercentage of splice in for skipping exon(-SE should be 1)\n";
	print "\t-H\t\tshow help information\n";
	$if_die = 1;
}elsif( !defined($out) or !defined($gtf) or !defined($coverage) or !defined($coverage2) or !defined($rand_mode) or !defined($rand_mode2) or !defined($read_length) or !defined($seq_err) or !defined($insert_length) or !defined($ref_dir) or !defined($if_chr) or !defined($sigma) or !defined($exon_skipping) or !defined($psi) ){
	$if_die = 2;
	#push @die_reason, "Please input complete arguments.\n";
}#elsif($insert_length <= $read_length){
#	$if_die = 1;
#	print "Insertion length should be larger than read length.\n";
#}
#die if $if_die == 1;



if($sigma <= 0){
	push @die_reason, "standard deviation of insert_length should be larger than 0!\n";
}
if($insert_length < $read_length){
	push @die_reason, "average of insert_length should be larger than read_length!\n";
}
if($read_length < 40){
	push @die_reason, "read_length should be larger than 40!\n";
}
if($sigma >= $insert_length){
	push @die_reason, "standard deviation of insert_length should be smaller than average of insert_length!\n";
}
if($psi>=100){
	push @die_reason, "percentage of splice in for skipping exon should be smaller than 100!\n";
}elsif($psi<=0 and $exon_skipping == 1){
	push @die_reason, "percentage of splice in for skipping exon should be larger than 0!\n";
}

if($if_die == 1){
	die "\n";
}elsif($if_die == 2){
	die "Please input complete arguments.\n";
}elsif (@die_reason > 0){
	print @die_reason;
	die;
}

$fq1 = $out."_1.fq";
$fq2 = $out."_2.fq";
$out = $out.".out";
$circ_read = $out."_circ_read.txt";
open FQ1, ">", $fq1 or die "cannot write to $fq1:$!";
open FQ2, ">", $fq2 or die;
open CIRCREAD, ">", $circ_read or die;
print CIRCREAD "circRNA_id\tread_id/bsj_read_id\n";

$insert_length -= $read_length;
my $pi = 3.14159265359;
#$out = ">>./".$out;
my %chr_gene_trsc_exon;
my %chr_gene_trsc_all;
my %chr_gene_trsc_intron;
my %chr_inter;
my %ss_exon_id;
my $ss_element_id;
my $pre_gene = '';
my @gene_anno;
my @chr;
my $seqID;
my $sim_total;
$ref_dir = $ref_dir."/" unless rindex($ref_dir, "/") == length($ref_dir) - 1;
open GTF, "<", $gtf or die "cannot open gtf file: $!";
open OUT, ">", $out or die;

while(<GTF>){
	chomp;
	next if /^#/;
	my @line = split /\t/;
	if ( ($if_chr == 1 and $line[0] ne "chr1") or $line[0] eq "chrM" ){
		last;
	}
	my @atr = split '; ', $line[8];
	if($pre_gene ne $atr[0] and $pre_gene ne ''){
		&split_transcript(@gene_anno);
		#&split_transcript_all(@gene_anno);
		@gene_anno = ();
		#print "$line[0]\t$pre_gene\n";
	}
	push @gene_anno, $_;
	$pre_gene = $atr[0];
}

#for my $chr(@chr){
#	for my $gene(keys %{$chr_$gene_trsc_exon{$chr}}){
#		for my $trsc(keys %{$chr_$gene_trsc_exon{$chr}{$gene}}){
#			for $exon(@{$chr_$gene_trsc_exon{$chr}{$gene}{$trsc}}){
#				print OUT "$chr\t$gene\t$trsc\t@$exon\n";
#			}
#		}
#	}
#}
for my $chromo(@chr){
	open CHR, "<", $ref_dir."$chromo.fa" or open CHR, "<", $ref_dir."$chromo.fasta" or die "cannot open the chr fasta file $chromo: $!";
	my $uni_seq = 0;
	my $chr_seq;
	while(<CHR>){
		chomp;
		if(/^>/ and $uni_seq == 0){
			$uni_seq = 1;
		}elsif(/^>/){
			die "There are more than one sequence in $chromo file. Please check!";
		}else{
			$chr_seq .= $_;
		}
	}

	my ($count, $gene_num, $inter_num, $circ_id);
	my @intron_gene;
	my $total_gene = keys %{$chr_gene_trsc_exon{$chromo}};

	print $chromo."\t".$total_gene."\n";
	##generating exon circ

	GEN: for my $gene(keys %{$chr_gene_trsc_exon{$chromo}}){
		$count += 1;
		if ($count > ($total_gene * 0.8)){
			push (@intron_gene, $gene);
			next GEN;
		}
		my @trscs = keys (%{$chr_gene_trsc_exon{$chromo}{$gene}});
		my $trsc_rand = int(rand($#trscs+1));
		my $trsc = $trscs[$trsc_rand];
			my %ss_exon_id;
			my ($cRNA_seq, $cRNA_seq2);
			my $trsc_seq;
			my $if_cRNA = int(rand(2));	#the probability of cRNA generated from this transcrpt
			my $if_linear = int(rand(2));
			for my $i(0 .. $#{$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}){
				my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$i];
				my $exon_seq;
				if($$exon[2] eq '+'){
					$exon_seq = substr($chr_seq, $$exon[0]-1, $$exon[1]-$$exon[0]+1);
				}else{
					$exon_seq = &comp_rev(substr($chr_seq, $$exon[0]-1, $$exon[1]-$$exon[0]+1));
				}
				push @$exon, $exon_seq;
				push @$exon, substr($chr_seq, $$exon[0]-3, 2);
				push @$exon, substr($chr_seq, $$exon[1], 2);
				$trsc_seq .= $$exon[-3];
				if ( $if_cRNA == 1 and ($$exon[1]-$$exon[0]+1>20 and $$exon[1]-$$exon[0]+1<2000) and ( ($$exon[2] eq '-' and $$exon[-2] =~ /AC/i and $$exon[-1] =~ /CT/i) or ($$exon[2] eq '+' and $$exon[-2] =~ /AG/i and $$exon[-1] =~ /GT/i) ) ){
					#$rand_exon1 ++;
					#$rand_num1{$rand_exon1} = $i;
					$ss_exon_id{$i} = 1;
				}
			}
			&simulate_reads2( $rand_mode2, $trsc_seq, $coverage2 ) if ($if_linear == 1);	#length($trsc_seq) > $insert_length and
			if( $if_cRNA == 1 and scalar(keys %ss_exon_id)>=3 ){	#and $rand_exon2>=1
				my @qualified_exon_sort = sort {$a <=> $b} (keys %ss_exon_id);
				my ($bingo_num1, $bingo_num2, $bingo_num3);
				my (%p_loci1, %p_loci2, $total_p1, $total_p2, $accu_p1, $accu_p2);
				for my $qualified_exon($qualified_exon_sort[0] .. $qualified_exon_sort[-1]){
					next GEN if !exists $ss_exon_id{$qualified_exon};
				}
				my $pre_loci1 = int($#qualified_exon_sort/2)-2*log(@qualified_exon_sort)/log(10);
				my $pre_loci2 = int($#qualified_exon_sort/2)+2*log(@qualified_exon_sort)/log(10)-1; 	########################################### add -1
				for my $i(0 .. int($#qualified_exon_sort/2)-$exon_skipping){
					$p_loci1{$i} = 1/(($i-$pre_loci1)**2+.01);
					$total_p1 += $p_loci1{$i};
				}
				for my $i(int($#qualified_exon_sort/2) .. $#qualified_exon_sort){		########################################### delete +1
					$p_loci2{$i} = 1/(($i-$pre_loci2)**2+.01);
					$total_p2 += $p_loci2{$i};
				}
				my $dice_loci1 = rand(1);
				for my $i(0 .. int($#qualified_exon_sort/2)-$exon_skipping){
					my $accu_pre = $accu_p1;
					$accu_p1 += $p_loci1{$i}/$total_p1;
					if ($dice_loci1 > $accu_pre and $dice_loci1 <= $accu_p1){
						$bingo_num1 = $i;
					}
				}
				if(!defined $bingo_num1){
					$bingo_num1 = int($#qualified_exon_sort/2)-1;
					print "1!!$dice_loci1\t", scalar(@qualified_exon_sort),"\n";
				}
				my $dice_loci2 = rand(1);
				for my $i(int($#qualified_exon_sort/2) .. $#qualified_exon_sort){
					my $accu_pre = $accu_p2;
					$accu_p2 += $p_loci2{$i}/$total_p2;
					if ($dice_loci2 > $accu_pre and $dice_loci2 <= $accu_p2){
						$bingo_num2 = $i;
					}
				}
				if(!defined $bingo_num2){
					$bingo_num2 = int($#qualified_exon_sort/2)+1;
					print "2!!$dice_loci2\t", scalar(@qualified_exon_sort),"\n";
				}
				next GEN if $bingo_num1>$bingo_num2;		########################################### add this to prevent error
				if(int(rand($exon_skipping + 1)) == 1){	#should have at least three exons  and $rand_num1{$bingo_num1} <= $rand_num1{$bingo_num2}-2
					$bingo_num3 = int(rand($bingo_num2-$bingo_num1-1)+$bingo_num1+1);
					if ($bingo_num3 < $bingo_num2 and $bingo_num3 > $bingo_num1){
						for my $j( $qualified_exon_sort[$bingo_num1] .. $qualified_exon_sort[$bingo_num2] ){
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
							$cRNA_seq2 .= $$exon[-3] unless $j == $qualified_exon_sort[$bingo_num3];
						}
						next GEN if length($cRNA_seq2) <100 or length($cRNA_seq2) >4000;
						next GEN if length($cRNA_seq) <100 or length($cRNA_seq) >4000;
							if ($chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0] < $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1]){
								print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0], "|", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1], "\texon", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0]."|".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1];
							}elsif($chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1] > $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0]){
								print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0], "|", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1], "\texon", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0]."|".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1];
							}else{
								next GEN;
							}
							for my $j( $qualified_exon_sort[$bingo_num1] .. $qualified_exon_sort[$bingo_num2] ){
								my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
								print OUT "$$exon[0]:$$exon[1]!$$exon[2],";
							}
							print OUT "\nisoform2_", length($cRNA_seq2), "\t";
							for my $j( $qualified_exon_sort[$bingo_num1] .. $qualified_exon_sort[$bingo_num2] ){
								my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
								print OUT "$$exon[0]:$$exon[1]!$$exon[2]," unless $j == $qualified_exon_sort[$bingo_num3];
							}
							print OUT "\n";
							&simulate_reads( $rand_mode, $cRNA_seq, $coverage*$psi/100, $circ_id );
							&simulate_reads( $rand_mode, $cRNA_seq2, $coverage*(100-$psi)/100, $circ_id );
							$sim_total++;

					}else{
						#print "!!!$bingo_num1\t$bingo_num2\t$bingo_num3\n";
					}
				}else{	#should have at least two exons	and $rand_num1{$bingo_num1} <= $rand_num1{$bingo_num2}-1
						for my $j( $qualified_exon_sort[$bingo_num1].. $qualified_exon_sort[$bingo_num2] ){
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
							#$cRNA_seq2 .= $$exon[-3] unless $j == $rand_num1{$bingo_num3};
						}
						next GEN if length($cRNA_seq) < 100 or length($cRNA_seq) >4000;
							if ($chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0] < $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1]){
								print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0], "|", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1], "\texon", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][0]."|".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][1];
							}elsif($chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1] > $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0]){
								print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0], "|", $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1], "\texon", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num2]][0]."|".$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$qualified_exon_sort[$bingo_num1]][1];
							}else{
								next GEN;
							}
							for my $j( $qualified_exon_sort[$bingo_num1] .. $qualified_exon_sort[$bingo_num2] ){
								my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
								print OUT "$$exon[0]:$$exon[1]!$$exon[2],";
							}
							print OUT "\n";
							&simulate_reads( $rand_mode, $cRNA_seq, $coverage, $circ_id );
							#&simulate_reads( $rand_mode, $cRNA_seq2, $coverage );
							$sim_total++;
				}
			}

	}

	##generating intron circRNA
	GEN_ALL: for my $gene(keys %{$chr_gene_trsc_all{$chromo}}){
		#print $gene."\n";

	  next unless grep /^$gene$/, @intron_gene;###next unless $gene ~~ @intron_gene;
		my ($cRNA_seq, $cRNA_seq2, $trsc_seq, %ss_element_id);
	  my @trscs = keys (%{$chr_gene_trsc_all{$chromo}{$gene}});
	  my $trsc_rand = int(rand($#trscs+1));
	  my $trsc = $trscs[$trsc_rand];
		my $if_cRNA = int(rand(2));
			### single intron circRNA
		  $gene_num += 1;
			if ($gene_num < ($#intron_gene * 0.05)){
				if ( $if_cRNA == 1 ) {
					my $i = int(rand($#{$chr_gene_trsc_intron{$chromo}{$gene}{$trsc}}));
					my $element = $chr_gene_trsc_intron{$chromo}{$gene}{$trsc}[$i];
					next GEN_ALL if $$element[1]-$$element[0]>20 and $$element[1]-$$element[0]<4000;
		      if($$element[2] eq '+'){
		        $cRNA_seq = substr($chr_seq, $$element[0]-1, $$element[1]-$$element[0]+1);
		      }else{
		        $cRNA_seq = &comp_rev(substr($chr_seq, $$element[0]-1, $$element[1]-$$element[0]+1));
		      }
					next GEN_ALL if length($cRNA_seq) < 100 or length($cRNA_seq) >4000;
					print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $$element[0], "|", $$element[1], "\tintron", "\nisoform1_", length($cRNA_seq), "\t";
					$circ_id = $chromo.":".$$element[0]."|".$$element[1];
					print OUT "\n";
					&simulate_reads( $rand_mode, $cRNA_seq, $coverage*$psi/100, $circ_id );
					$sim_total++;
					next GEN_ALL;
				} else {
					next GEN_ALL;
				}
			}

			for my $i(0 .. $#{$chr_gene_trsc_all{$chromo}{$gene}{$trsc}}){
	      my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$i];
	      my $element_seq;
	      if($$element[2] eq '+'){
	        $element_seq = substr($chr_seq, $$element[0]-1, $$element[1]-$$element[0]+1);
	      }else{
	        $element_seq = &comp_rev(substr($chr_seq, $$element[0]-1, $$element[1]-$$element[0]+1));
	      }
	      push @$element, $element_seq;
	      push @$element, substr($chr_seq, $$element[0]-3, 2);
	      push @$element, substr($chr_seq, $$element[1], 2);
	      $trsc_seq .= $$element[-3];
	      if ( $if_cRNA == 1 and ($$element[1]-$$element[0]+1>20 and $$element[1]-$$element[0]+1<2000) ){
	        $ss_element_id{$i} = 1;
	      }
	    }

			if( $if_cRNA == 1 and scalar(keys %ss_element_id)>=3 ){	#and $rand_element2>=1
	      my @qualified_element_sort = sort {$a <=> $b} (keys %ss_element_id);
	      my ($bingo_all_num1, $bingo_all_num2, $bingo_all_num3);
	      my (%p_all_loci1, %p_all_loci2, $total_all_p1, $total_all_p2, $accu_all_p1, $accu_all_p2);
	      for my $qualified_element($qualified_element_sort[0] .. $qualified_element_sort[-1]){
	        next GEN_ALL if !exists $ss_element_id{$qualified_element};
	      }
	      my $pre_all_loci1 = int($#qualified_element_sort/2)-2*log(@qualified_element_sort)/log(10);
	      my $pre_all_loci2 = int($#qualified_element_sort/2)+2*log(@qualified_element_sort)/log(10)-1; 	########################################### add -1
	      for my $i(0 .. int($#qualified_element_sort/2)-$exon_skipping){
	        $p_all_loci1{$i} = 1/(($i-$pre_all_loci1)**2+.01);
	        $total_all_p1 += $p_all_loci1{$i};
	      }
	      for my $i(int($#qualified_element_sort/2) .. $#qualified_element_sort){		########################################### delete +1
	        $p_all_loci2{$i} = 1/(($i-$pre_all_loci2)**2+.01);
	        $total_all_p2 += $p_all_loci2{$i};
	      }
	      my $dice_all_loci1 = rand(1);
	      for my $i(0 .. int($#qualified_element_sort/2)-$exon_skipping){
	        my $accu_all_pre = $accu_all_p1;
	        $accu_all_p1 += $p_all_loci1{$i}/$total_all_p1;
	        if ($dice_all_loci1 > $accu_all_pre and $dice_all_loci1 <= $accu_all_p1){
	          $bingo_all_num1 = $i;
	        }
	      }
	      if(!defined $bingo_all_num1){
	        $bingo_all_num1 = int($#qualified_element_sort/2)-1;
	        print "1!!$dice_all_loci1\t", scalar(@qualified_element_sort),"\n";
	      }
	      my $dice_all_loci2 = rand(1);
	      for my $i(int($#qualified_element_sort/2) .. $#qualified_element_sort){
	        my $accu_all_pre = $accu_all_p2;
	        $accu_all_p2 += $p_all_loci2{$i}/$total_all_p2;
	        if ($dice_all_loci2 > $accu_all_pre and $dice_all_loci2 <= $accu_all_p2){
	          $bingo_all_num2 = $i;
	        }
	      }
	      if(!defined $bingo_all_num2){
	        $bingo_all_num2 = int($#qualified_element_sort/2)+1;
	        print "2!!$dice_all_loci2\t", scalar(@qualified_element_sort),"\n";
	      }
	      next GEN_ALL if $bingo_all_num1>$bingo_all_num2;		########################################### add this to prevent error
	      if($exon_skipping == 1){	#should have at least three elements  and $rand_num1{$bingo_all_num1} <= $rand_num1{$bingo_all_num2}-2
	        $bingo_all_num3 = int(rand($bingo_all_num2-$bingo_all_num1-1)+$bingo_all_num1+1);
	        if ($bingo_all_num3 < $bingo_all_num2 and $bingo_all_num3 > $bingo_all_num1){
	          for my $j( $qualified_element_sort[$bingo_all_num1] .. $qualified_element_sort[$bingo_all_num2] ){
	            my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$j];
	            $cRNA_seq .= $$element[-3];
	            $cRNA_seq2 .= $$element[-3] unless $j == $qualified_element_sort[$bingo_all_num3];
	          }
	          next GEN_ALL if length($cRNA_seq2) <100 or length($cRNA_seq2) >4000;
						next GEN_ALL if length($cRNA_seq) <100 or length($cRNA_seq) >4000;	# or length($cRNA_seq) > 850;
	            if ($chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0] < $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1]){
	              print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0], "|", $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1], "\tintron", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0]."|".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1];
	            }elsif($chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1] > $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0]){
	              print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0], "|", $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1], "\tintron", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0]."|".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1];
	            }else{
	              next GEN_ALL;
	            }
	            for my $j( $qualified_element_sort[$bingo_all_num1] .. $qualified_element_sort[$bingo_all_num2] ){
	              my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$j];
	              print OUT "$$element[0]:$$element[1]!$$element[2],";
	            }
	            print OUT "\nisoform2_", length($cRNA_seq2), "\t";
	            for my $j( $qualified_element_sort[$bingo_all_num1] .. $qualified_element_sort[$bingo_all_num2] ){
	              my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$j];
	              print OUT "$$element[0]:$$element[1]!$$element[2]," unless $j == $qualified_element_sort[$bingo_all_num3];
	            }
	            print OUT "\n";
	            &simulate_reads( $rand_mode, $cRNA_seq, $coverage*$psi/100, $circ_id );
	            &simulate_reads( $rand_mode, $cRNA_seq2, $coverage*(100-$psi)/100, $circ_id );
	            $sim_total++;

	        }else{
	          print "!!!$bingo_all_num1\t$bingo_all_num2\t$bingo_all_num3\n";
	        }
	      }elsif($exon_skipping == 0){	#should have at least two elements	and $rand_num1{$bingo_all_num1} <= $rand_num1{$bingo_all_num2}-1
	          for my $j( $qualified_element_sort[$bingo_all_num1].. $qualified_element_sort[$bingo_all_num2] ){
	            my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$j];
	            $cRNA_seq .= $$element[-3];
	            #$cRNA_seq2 .= $$element[-3] unless $j == $rand_num1{$bingo_all_num3};
	          }
	          next GEN_ALL if length($cRNA_seq) < 100 or length($cRNA_seq) >4000;
	            if ($chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0] < $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1]){
	              print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0], "|", $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1], "\tintron", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][0]."|".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][1];
	            }elsif($chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1] > $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0]){
	              print OUT "$chromo\t$gene\t$trsc\t", $chromo, ":" , $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0], "|", $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1], "\tintron", "\nisoform1_", length($cRNA_seq), "\t";
								$circ_id = $chromo.":".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num2]][0]."|".$chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$qualified_element_sort[$bingo_all_num1]][1];
	            }else{
	              next GEN_ALL;
	            }
	            for my $j( $qualified_element_sort[$bingo_all_num1] .. $qualified_element_sort[$bingo_all_num2] ){
	              my $element = $chr_gene_trsc_all{$chromo}{$gene}{$trsc}[$j];
	              print OUT "$$element[0]:$$element[1]!$$element[2],";
	            }
	            print OUT "\n";
	            &simulate_reads( $rand_mode, $cRNA_seq, $coverage, $circ_id );
	            #&simulate_reads( $rand_mode, $cRNA_seq2, $coverage );
	            $sim_total++;
	      }
	    }

	}

  ###generating intergenic circRNA
	INTER: for my $i(0 .. $#{$chr_inter{$chromo}}){
		my $region = $chr_inter{$chromo}[$i];
	  my $inter_total_num = $#intron_gene * 0.1;
		$inter_num += 1;
		next INTER if $inter_num > $inter_total_num;
		my ($if_cRNA, $pre_length, $inter_circ_start, $inter_circ_end, $inter_circ_strand, $cRNA_seq);
		  $if_cRNA = int(rand(2));
			next INTER if $if_cRNA == 0;
		  $pre_length = length($$region[1] - $$region[0] + 1);
			$inter_circ_start = int(rand($pre_length * 0.5)) + $$region[0];
			$inter_circ_end = $inter_circ_start + 100 + int(rand(3800));
			if (int(rand(2)) == 1){
				$inter_circ_strand = "+";
			} else {
				$inter_circ_strand = "-";
			}
			if($inter_circ_strand eq '+'){
				$cRNA_seq = substr($chr_seq, $inter_circ_start-1, $inter_circ_end-$inter_circ_start+1);
			}else{
				$cRNA_seq = &comp_rev(substr($chr_seq, $inter_circ_start-1, $inter_circ_end-$inter_circ_start+1));
			}

			print OUT "$chromo\t", $chromo, ":" , $inter_circ_start, "|", $inter_circ_end, "\tintergenic", "\nisoform1_", length($cRNA_seq), "\t", $inter_circ_start, ":", $inter_circ_end, "!", $inter_circ_strand, ",";
			$circ_id = $chromo.":".$inter_circ_start."|".$inter_circ_end;
			print OUT "\n";
			&simulate_reads( $rand_mode, $cRNA_seq, $coverage*$psi/100, $circ_id );
			$sim_total++;
	}

}
print OUT "!!total: $sim_total\n";
print "circRNA totally simulated: $sim_total\n";

sub simulate_reads2{
	my $mode = shift @_;
	my $trsc_coverage;
	my $seq_length = length($_[0]);
	if ($mode == 1){
		$trsc_coverage = $_[1];
	}else{
		$trsc_coverage = rand($_[1]+1);
	}
	my ($read_num, undef) = sort{$b <=> $a}(int( $seq_length * $trsc_coverage / $read_length / 2 ),1);
	my $err_num = int( $seq_length * $trsc_coverage * $seq_err / 100 );
	my %err_read;
	for (1 .. $err_num){
		my $err_loci = int( rand( (($read_num)*2)+1 ) );
		$err_read{$err_loci} ++;
	}
	for my $x( 1 .. $read_num ){
		my $ins_len_rand = &insert_length_calculation_normal(rand(1), rand(1), rand(1));		#insert length can be simulated later
		#my $start_loci = int( rand($seq_length - $ins_len_rand - $read_length) );
		next if $ins_len_rand < 0 or $ins_len_rand >= $seq_length-$read_length;
		my ($start_loci, $start_loci2);
		$start_loci = int( rand($seq_length-$ins_len_rand-$read_length+1) );
		$start_loci2 = $start_loci + $ins_len_rand;	# - $read_length;
		next if $start_loci2 > $seq_length-$read_length;
		$seqID ++;
		my $if_1st = int(rand(2));
		my ($seq1, $seq2);
		if ($if_1st == 1){
			$seq1 = substr( $_[0], $start_loci, $read_length );
			$seq2 = &comp_rev( substr( $_[0], $start_loci2, $read_length ) );
		}else{
			$seq2 = substr( $_[0], $start_loci, $read_length );
			$seq1 = &comp_rev( substr( $_[0], $start_loci2, $read_length ) );
		}
		die "$seq1" if length($seq1) != $read_length;
		die "$seq2" if length($seq2) != $read_length;
		my @errs1;
		for (1 .. $err_read{$x * 2 - 1}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs1;
			push @errs1, $err_loci;
			my $ori_base = substr( $seq1, $err_loci, 1 );
			substr( $seq1, $err_loci, 1 ) = &simulate_seq_error($ori_base);
		}
		my @errs2;
		for (1 .. $err_read{$x * 2}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs2;
			push @errs2, $err_loci;
			my $ori_base = substr( $seq2, $err_loci, 1 );
			substr( $seq2, $err_loci, 1 ) = &simulate_seq_error($ori_base);
		}
		print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
		print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
		print FQ1 "$seq1\n";
		print FQ2 "$seq2\n";
		print FQ1 "+\n";
		print FQ2 "+\n";
		print FQ1 ("!" x $read_length) . "\n";
		print FQ2 ("!" x $read_length) . "\n";
	}
}

sub insert_length_calculation_linear{
	if($_[0]<=0.5){
		my $x = int(($_[0]*2)**.5*$insert_length);
	}else{
		my $x = int((2-(2-2*$_[0])**.5)*$insert_length);
	}
}

sub insert_length_calculation_normal{
	my $y;
	if($_[2] < $perc_minor/100){
		$y = sqrt(-2*log($_[0]))*cos(2*$pi*$_[1])*$sigma2+$insert_length2;
	}else{
		$y = sqrt(-2*log($_[0]))*cos(2*$pi*$_[1])*$sigma+$insert_length;
	}
}

sub simulate_reads{
	my $mode = shift @_;
	my $seq_length = length($_[0]);
	my $cRNA_coverage;
	my $junc_read = ();
	my $circ_id = $_[-1];
	#my $seq4substr;
	if ($mode == 1){
		$cRNA_coverage = $_[1];
	}else{
		$cRNA_coverage = rand($_[1]+1);
	}
	my $seq4substr = $_[0] x 12;
	my $read_num = int( $seq_length * $cRNA_coverage / $read_length / 2 );
	my $err_num = int( $seq_length * $cRNA_coverage * $seq_err / 100 );
	my %err_read;
	for (1 .. $err_num){
		my $err_loci = int( rand( ($read_num)*2 ) +1);
		$err_read{$err_loci} ++;
	}
	for my $x( 1 .. $read_num ){
		my $ins_len_rand = &insert_length_calculation_normal(rand(1), rand(1), rand(1));		#insert length can be simulated later
		next if $ins_len_rand < 0; #( or $ins_len_rand + $read_length > $seq_length );
		my ($start_loci, $start_loci2);
		$start_loci = int( rand($seq_length) );
		$start_loci2 = $start_loci + $ins_len_rand;	# - $read_length;
		$seqID ++;
		my $if_1st = int( rand(2) );
		my ($seq1, $seq2);

		$junc_read .= $seqID."," if ($start_loci > ($seq_length - $read_length)) and ($start_loci < $seq_length);
		for my $times (1 .. 9) {
			$junc_read .= $seqID."," if $start_loci2 > (($seq_length * $times) - $read_length) and $start_loci2 < ($seq_length * $times);
		}

		if ($if_1st == 1){
			$seq1 = substr( $seq4substr, $start_loci, $read_length );
			$seq2 = &comp_rev( substr( $seq4substr, $start_loci2, $read_length ) );
		}else{
			$seq2 = substr( $seq4substr, $start_loci, $read_length );
			$seq1 = &comp_rev( substr( $seq4substr, $start_loci2, $read_length ) );
		}
		die "$seq2" if length($seq2) != $read_length;
		die "$seq1" if length($seq1) != $read_length;
		my @errs1;
		for (1 .. $err_read{$x * 2 - 1}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs1;
			push @errs1, $err_loci;
			my $ori_base = substr( $seq1, $err_loci, 1 );
			substr( $seq1, $err_loci, 1) = &simulate_seq_error($ori_base);
		}
		my @errs2;
		for (1 .. $err_read{$x * 2}){
			my $err_loci = int( rand($read_length) );
			redo if $err_loci ~~ @errs2;
			push @errs2, $err_loci;
			my $ori_base = substr( $seq2, $err_loci, 1 );
			substr( $seq2, $err_loci, 1) = &simulate_seq_error($ori_base);
		}
		print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
		print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
		print FQ1 "$seq1\n";
		print FQ2 "$seq2\n";
		print FQ1 "+\n";
		print FQ2 "+\n";
		print FQ1 ("!" x $read_length) . "\n";
		print FQ2 ("!" x $read_length) . "\n";
		print OUT ">\t$x\t$seqID\n";
		print CIRCREAD "$circ_id\tsimulate:"."$seqID\n";
		if($if_1st == 1){
			print OUT "**\t1\n" if ($start_loci >= $seq_length-$read_length+20 and $start_loci <= $seq_length-20);
			print OUT "**\t2\n" if ($start_loci2%$seq_length >= $seq_length-$read_length+20 and $start_loci2%$seq_length <= $seq_length-20);
		}else{
			print OUT "**\t1\n" if ($start_loci2%$seq_length >= $seq_length-$read_length+20 and $start_loci2%$seq_length <= $seq_length-20);
			print OUT "**\t2\n" if ($start_loci >= $seq_length-$read_length+20 and $start_loci <= $seq_length-20);
		}
	}
	my (@junc, $tmp, %exist_junc, $junc_read_id);
	if ($junc_read =~ /\d+/){
		@junc = split /,/,$junc_read;
		for $tmp (@junc){
			$exist_junc{$tmp} = 1;
		}
		for (sort keys %exist_junc){
			$junc_read_id .= "simulate:".$_.",";
		}
		print CIRCREAD "$circ_id\tbsj:"."$junc_read_id\n";
	}
	
}
sub simulate_seq_error{
	my $ori_base = $_[0];
	my @base = ('A', 'T', 'C', 'G');
	my $err_base_index;
	for my $i( 0 .. $#base ){
		if ($base[$i] =~ /$ori_base/i){
			while(1){
				$err_base_index = int(rand(4));
				last unless $err_base_index  == $i;
			}
			last;
		}
	}
	$base[$err_base_index];
}
sub comp_rev{
	my $seq = reverse($_[0]);
	$seq =~ s/[Aa]/X/g;
	$seq =~ s/[Tt]/A/g;
	$seq =~ s/X/T/g;
	$seq =~ s/[Cc]/Y/g;
	$seq =~ s/[Gg]/C/g;
	$seq =~ s/Y/G/g;
	$seq;
}

sub split_transcript{
	my ($count, $intron_end, $intron_start, $intron_strand, $inter_start, $inter_end, $transcript, %trans_exist, %chr_inter_exist, %intron);
	for (@_){
		my @line = split /\t/;
		next if $line[0] =~ /\_/;
		if($line[2] eq 'exon'){
			my @atr = split ('; ', $line[8], 3);
			if ($atr[1] =~ /transcript_id \"(\S+?)\"/){
				##exon
				push @{$chr_gene_trsc_exon{$line[0]}{$atr[0]}{$1}}, [ $line[3], $line[4], $line[6] ];
				##exon and intron
				unless ($intron{$1}) {
					$intron{$1} = 1;
					$intron_start = $line[4];
					push @{$chr_gene_trsc_all{$line[0]}{$atr[0]}{$1}}, [ $line[3], $line[4], $line[6] ];
				} else {
					$intron_end = $line[3];
					$intron_strand = $line[6];
					push @{$chr_gene_trsc_all{$line[0]}{$atr[0]}{$1}}, [ $intron_start, $intron_end, $intron_strand ];
					push @{$chr_gene_trsc_intron{$line[0]}{$atr[0]}{$1}}, [ $intron_start, $intron_end, $intron_strand ];
					$intron_start = $line[4];
					push @{$chr_gene_trsc_all{$line[0]}{$atr[0]}{$1}}, [ $line[3], $line[4], $line[6] ];
				}
				#$trans_exist{$1} = 1 and $intron_start = () unless $chr_gene_trsc_all{$line[0]}{$atr[0]}{$1};
			}else{
				print "error: no transcript_id found for $atr[0]!\n";
			}
		}
    ##intergenic

		if($line[2] eq 'transcript'){
			$chr_inter_exist{$line[0]} = 1 and $inter_start = () unless $chr_inter_exist{$line[0]};
			unless ($inter_start){
				$inter_start = $line[4];
		  } else {
				$inter_end = $line[3];
				push @{$chr_inter{$line[0]}}, [ $inter_start, $inter_end ];
				$inter_start = $line[4];
			}
		}


	}
	#print $_[0]."\n";
	my @line2 = split (/\t/, $_[0], 2);
	push @chr, $line2[0] unless $line2[0] ~~ @chr;
}
