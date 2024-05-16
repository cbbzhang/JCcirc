#!/usr/bin/perl

=begin SUMMARY
Main steps:
1. Check and record parameters designated by user
2. Prepare files to process
3. Obtaining junction reads and junction contigs
4. Mapping junction read sequences and junction contig sequences to the concatenate circRNA genomic sequences
5. Organizing all the fragments and filtering out low-confidence fragments
6. Generating alternative splicing circRNA
7. Output
=cut

use strict;
use Thread;
use threads::shared;
use Getopt::Long;
my $version = '1.0.0';
open STDERR, ">>./JCcirc.log";
### Parameters can be designated by user
my ( $circ, $genome, $output, $help, $gtf, $thread, $read1, $read2, $contig, $difference );
Getopt::Long::GetOptions (
        'circ|C=s'                                =>      \$circ,
        'genome|G=s'                              =>      \$genome,
        'annotation|F=s'                          =>      \$gtf,
        'output|O=s'                              =>      \$output,
        'thread|P=i'                              =>      \$thread,
        'read1=s'                                 =>      \$read1,
        'read2=s'                                 =>      \$read2,
        'contig=s'                                =>      \$contig,
        'difference|D=i'                          =>      \$difference,
        'help|H!'                                 =>      \$help,
);

if ( !defined($circ) and !defined($genome) and  !defined($output) and !defined($help) and !defined($gtf) and !defined($read1)and !defined($read2)and !defined($contig)) {
  print "Please use the --help or -H option to get usage information.\n";
### Show help information to user if requested
} elsif (defined($help)) {
  print "
Program:  CIRCSeq (circRNA sequence)
Version:  $version

Usage:    perl cFLSeq.pl -C circ -G genome -F annotation -O circ_contig -P 8 --read1 read_1.fq --read2 read_2.fq --contig contig.fa -D 0

Arguments:

    -C, --circ
          input circRNA file, which including chromosome, start site, end site, host gene, and junction reads ID (required).
    -O, --output
          directory of output (required).
    -G, --genome
          FASTA file of all reference sequences. Please make sure this file is
          the same one provided to prediction tool (required).
    -F, --annotation
          gene annotation file in gtf format. Please make sure this file is
          the same one provided to prediction tool.
    -P, --thread
          set number of threads for parallel running, default is 4.
    --read1
          RNA-Seq data, read_1 (paired end, fastq format).
    --read2
          RNA-Seq data, read_1 (paired end, fastq format).
    --contig
          contig sequences (required).
    --difference
          the difference in support numbers between adjacent fragments when generating circRNA isoforms, default is 0 (recommend setting to 0, 1, or 2, the larger number means stricter).
    -H, --help
          show this help information.

";
} elsif ( !defined($circ) or !defined($genome) or !defined($output) or !defined($read1) or !defined($read2) or !defined($contig)) {
  print "Please check input file.\n";
} elsif ($difference > 2) {
  print "Error: Parameter difference|D is greater than 2, it can be set to 0, 1, or 2\n";
} else {
  $difference = 0 if (!defined($difference));
  $thread = 4 if (!defined($thread));
  mkdir( $output ) or die "can not make a new output directory" unless (-e $output);
  our( $all, $frag_res_directory, $fragment_result_all, $contig_frag_res, $frag_sort_contig, $fragment_result_sort, $map_res_directory, $fragment_anno, $fragment_final, $circ_full_seq );
  $frag_res_directory = "$output/fragment";
  mkdir( $frag_res_directory ) or die "can not make a new output directory" unless (-e $frag_res_directory);
  $map_res_directory = "$output/map";
  mkdir( $map_res_directory ) or die "can not make a new output directory" unless (-e $map_res_directory);
  $fragment_result_all = "$frag_res_directory/fragment_all.txt";
  $contig_frag_res = "$frag_res_directory/contig_frag.txt";
  $frag_sort_contig = "$frag_res_directory/fragment_sort_contig.txt";
  $fragment_result_sort = "$frag_res_directory/fragment_sort.txt";
  $fragment_anno = "$frag_res_directory/fragment_anno.txt" if $gtf;
  $fragment_final = "$output/fragment_final.txt";
  $circ_full_seq = "$output/circ_full_seq.fa";
  open FRAGCONTIG, ">", $contig_frag_res or die "cannot write to $contig_frag_res.";
  open FRAGALL, ">", $fragment_result_all or die "cannot write to $fragment_result_all.";
  open FRAGSORTCONTIG, ">", $frag_sort_contig or die "cannot write to $contig_frag_res.";
  open FRAGSORT, ">", $fragment_result_sort or die "cannot write to $fragment_result_sort.";
  open FRAGANNO, ">", $fragment_anno or die "cannot write to $fragment_anno." if $gtf;
  open FRAGFINAL, ">", $fragment_final or die "cannot write to $fragment_final.";
  open CIRCSEQ, ">", $circ_full_seq or die "cannot write to $circ_full_seq.";

  ###### Step 1: Scanning circRNA information, genome sequences, and annotation information #####
  print 		'[', scalar(localtime), "] Step 1: Scanning circRNA information, genome sequences, and annotation information.\n";
  print 		'[', scalar(localtime), "]   Scanning genome sequences.\n";
  my ( $chrom );
  our ( %genome, %gene_exon );
  open GENOME, "$genome" or (print "Can not open $genome!\n" and die);###genome.fa
  while ( <GENOME> ){
    chomp;
    $chrom = $1 and next if $_ =~ />(\S+)/;
    $genome{$chrom} .= $_;
  }

  if ($gtf){
    print 		'[', scalar(localtime), "]   Scanning annotation information.\n";
    my( @gene_info, $exon, $gene, $intron_start, $intron_end );

    open GTF, "$gtf" or (print "Can not open $gtf!\n" and die);
    while(<GTF>){
      chomp;
      next unless $_ =~ /\S+/;
      next if $_ =~ /\#/;
      @gene_info = split /\t/,$_;
      next unless $gene_info[2] =~ /exon/;
      $exon = $gene_info[3]."-".$gene_info[4];
      $gene = $1 if $_ =~ /gene_id "(\S+?)";/;
      $gene_exon{$gene} .= $exon.",";
    }
  }

  print 		'[', scalar(localtime), "]   Scaning circRNA information.\n";
  our( %circ_to_gene, %junc_exist, %junc_contig_exist );
  my( @circ_info, @junction_reads );
  my( $in3, $circ_id, $circ_gene, $circ_line );#######$in and $out is file handle
  
  open $in3,"$circ" or (print "Can not open $circ!\n" and die);
  while(<$in3>){#print $_;
    $circ_line += 1;
    @circ_info = split /\t/,$_;
    $circ_id = $circ_info[0].":".$circ_info[1]."|".$circ_info[2];
    $circ_to_gene{$circ_id} = $circ_info[3];
    @junction_reads = split /\,/,$circ_info[5];
    for my $junction_read (@junction_reads) {
      $junc_exist{$junction_read} = 1;
    }
  }

  print 		'[', scalar(localtime), "]   Total circRNA: $circ_line.\n";
  print 		'[', scalar(localtime), "] Step 1 Finished!\n";

  ###### Second step: Extracting junction reads from sequencing data#####
  print 		'[', scalar(localtime), "] Step 2: Obtaining junction reads and junction contigs.\n";
  print 		'[', scalar(localtime), "]   Obtaining junction reads.\n";
  my ($in4, $in5, $in6, $in7, $in8, $in9, $read_id, $count, $contig_id, @map_info, $map_length, $junc_id);
  our (%read_seq, %contig_seq, %read_contig, $read_length);

  open $in4, "$read1" or (print "Can not open $read1!\n" and die);
  while(<$in4>){
    chomp;
    $read_id = $1 and $junc_id = $2 and $count = 1 and next if $_ =~ /\@((\S+)\/(\S+))/;
    $count = () and $read_id = () and next unless $junc_exist{$junc_id};
    if ($count == 1){
      $read_length = length($_) unless $read_length =~ /\d+/;
      $read_seq{$read_id} = $_;
      $count = 0;
    }
    $read_id = ();
  }

  open $in5, "$read2" or (print "Can not open $read2!\n" and die);
  while(<$in5>){
    chomp;
    $read_id = $1 and $junc_id = $2 and $count = 1 and next if $_ =~ /\@((\S+)\/(\S+))/;
    $count = () and $read_id = () and next unless $junc_exist{$junc_id};
    if ($count == 1){
      $read_seq{$read_id} = $_;
      $count = 0;
    }
    $read_id = ();
  }

  print 		'[', scalar(localtime), "]   Mapping RNA-Seq data to contig sequences.\n";

  if (-s($contig) and -s($read1) and -s($read2)) {
    system "bwa index $contig";
    system "bwa mem -t $thread $contig $read1 $read2 > $output/read_to_contig.sam";
  } else {
    print "Please check reads file and contig file!\n" and die;
  }

  print 		'[', scalar(localtime), "]   Obtaining junction contigs.\n";
  open $in6, "$output/read_to_contig.sam" or (print "Can not open $output/read_to_contig.sam" and die);
  while(<$in6>){
    next if $_ =~ /\@SQ/ or $_ =~ /\@PG/;
    @map_info = split /\t/,$_;
    next if $map_info[2] =~ /\*/;
    next unless $junc_exist{$map_info[0]};
    $junc_contig_exist{$map_info[2]} = 1;
    if ($map_info[5] =~ /(\S+)M/){
      $map_length = $1;
      next if ($map_length * 2) < $read_length;
      $read_contig{$map_info[0]} .= $map_info[2].",";
    } elsif ($map_info[5] =~ /(\S+)M(\S+)D(\S+)M/) {
      $map_length = $1 + $3;
      next if ($map_length * 2) < $read_length;
      $read_contig{$map_info[0]} .= $map_info[2].",";
    } elsif ($map_info[5] =~ /(\S+)M(\S+)I(\S+)M/) {
      $map_length = $1 + $3;
      next if ($map_length * 2) < $read_length;
      $read_contig{$map_info[0]} .= $map_info[2].",";
    }
  }
  system "rm $contig\.*";
  

  open $in7, "$contig" or (print "Can not open $contig!\n" and die);##build relationship between contig id and sequences, and delete contig that shorter than read length
  while(<$in7>){
    chomp;
    $contig_id = $1 and next if $_ =~ />(\S+)/;
    $contig_id = () and next unless $junc_contig_exist{$contig_id};
    $contig_seq{$contig_id} .= $_;
  }

  print 		'[', scalar(localtime), "] Step 2 Finished!\n";

 ###Build connnection between circRNAs and contig sequences
  print 		'[', scalar(localtime), "] Step 3: Mapping junction read sequences and junction contig sequences to the concatenate circRNA genomic sequences.\n";
  my( @circ_file, @ths, $num, $t );
  our ($circ_row, $circ_number);
  our (%contig_frag);
  @circ_file = &split_circ_file($circ, $thread);
  
  @ths=();
  for my $circ_part (@circ_file){
    push @ths, threads -> create( \&circ_contig_read, $circ_part, $output );
  }
  foreach $t (@ths) {
    $t -> join();
  }

  print 		'[', scalar(localtime), "] Step 3 Finished!\n";
  system "rm -rf $output/map";

  print 		'[', scalar(localtime), "] Step 4: Organizing all the fragments and filtering out low-confidence fragments\n";
  ####sort and count fragments
  print 		'[', scalar(localtime), "]   Organizing all the fragments and related contigs.\n";
  my( %frag_count, %frag_set, %frag_exist, %contig_frag );
  my( @each_circ_frag, @each_info, @fragment_set, @fragment_set_sort, @fragment_new, @fragment_new_sort );
  my( $each_circ, $circ_id, $count, $new_count, $frag, $tmp, $tmp_start, $tmp_end, $frag_start, $frag_end, $fragment_set, $frag_start_new, $new_frag );
  if (-s($contig_frag_res)){
    open $in8, "$contig_frag_res" or die;
    while(<$in8>){
      chomp;
      $contig_frag{$1}{$2} = $3 if $_ =~ /(\S+)\s+(\S+)\s+(\S+)/;
    }
  } else {
    print "$contig_frag_res is empty!" and die;
  }
  #print "Total number: ".$all."\n";
  if (-s($fragment_result_all)) {
    open $in9, "$fragment_result_all" or die;
    while (<$in9>){
      chomp;
      $each_circ = $_;
      @each_info = split /\t/, $each_circ;
      $circ_id = $each_info[0];
      #next unless $filter1{$circ_id};#delete no cover junction site
      @fragment_set = split /\,/, $each_info[1];
      @fragment_set_sort = &sort_frag(@fragment_set);
      %frag_count=(); $count=(); $fragment_set=(); $tmp=();$frag=();@fragment_new=();
      for $frag ( @fragment_set_sort ){
        next unless $frag =~ /\d+/;
        if ($frag_count{$circ_id}{$frag}){
          $frag_count{$circ_id}{$frag} += 1;
        } else {
          push(@fragment_new, $frag);
          $frag_count{$circ_id}{$frag} += 1;
        }
      }

      @fragment_new_sort = &sort_frag(@fragment_new);
      for $frag (@fragment_new_sort) {
        next unless $frag_count{$circ_id}{$frag} =~ /\d+/;
        #print $circ_id."\t".$frag."/".$frag_count{$circ_id}{$frag}."\n";
        $fragment_set .= $frag."/".$frag_count{$circ_id}{$frag}.",";
        print FRAGSORTCONTIG "$circ_id\t$frag\t$contig_frag{$circ_id}{$frag}\n";
      }
      print FRAGSORT "$circ_id\t$fragment_set\n";
    }
  } else {
    print "$fragment_result_all is empty!" and die;
  }

  my (%frag_count, %contig_frag, %frag_gene, @fragment, @circ_frag, @fragment_anno, @circ_genes, @exons, @fragment_exon, @info );
  my ( $line, $circ_id, $circ_gene, $circ_chr, $circ_start, $circ_end, $fragment, $frag, $frag_start, $frag_end, $gene, $exon, $exon_cluster, $exon_start, $exon_end, $count);
  my (%frag_exon, $last_frag, $fragment_anno_all, $frag_anno_set);

  if (-s($frag_sort_contig)){
    open $in8, "$frag_sort_contig" or die;
    while(<$in8>){
      chomp;
      @info = split /\t/,$_;
      $contig_frag{$info[0]}{$info[1]} = $info[2];
    }
  } else {
    print "$frag_sort_contig is empty!" and die;
  }
  
  my $sort_line = `system "wc -l $fragment_result_sort"`;

  if ( -s($fragment_result_sort) ) {
    open SORT, "$fragment_result_sort" or die;
    while(<SORT>){
      chomp;
      $line = $_;
      next unless $line =~ /(\S+)/;
      %frag_count=();@fragment_anno=();@circ_genes=();@exons=();@fragment=();#@fragment_set_new=();
      %frag_gene = ();
      @circ_frag = split /\t/,$line;
      $circ_id = $circ_frag[0];
      $circ_chr = $1 and $circ_start = $2 and $circ_end = $3 if $circ_id =~ /(\S+)\:(\S+)\|(\S+)/;
      $circ_gene = $circ_to_gene{$circ_id};

      ####Use annotation file to complete circexon

      if ($gtf) {
        @circ_genes = split /,/,$circ_gene;
        @fragment = split /,/, $circ_frag[1];
        $last_frag=();@fragment_exon=();
        for (@fragment){
          next unless $_ =~ /\d+/;
          $frag = $1 and $count = $2 if $_ =~ /(\S+)\/(\S+)/;
          $frag_start = $circ_start + $1 and $frag_end = $circ_start + $2 if $frag =~ /(\S+)\-(\S+)/;
          for $gene (@circ_genes) {#print $gene."\n";
            next unless $gene =~ /\S+/;
            $exon_cluster = $gene_exon{$gene};
            @exons = split /,/, $exon_cluster;
            for $exon (@exons) {#print $exon."\n";
              $exon_start = $1 and $exon_end = $2 if $exon =~ /(\S+)\-(\S+)/;
              if ((($frag_start + 20) > $exon_start) and (($frag_start - 20) < $exon_start)) {
                $frag_start = $exon_start;
              }
              if ((($frag_end + 20) > $exon_end) and (($frag_end - 20) < $exon_end)) {
                $frag_end = $exon_end;
              }
            }
          }

          $fragment = $frag_start."-".$frag_end;

          unless ($fragment eq $frag) {
            $contig_frag{$circ_id}{$fragment} .= $contig_frag{$circ_id}{$frag};
            #print $frag."\t".$contig_frag{$circ_id}{$frag}."\t".$fragment."\t".$contig_frag{$circ_id}{$fragment}."\n";
            delete($contig_frag{$circ_id}{$frag});

          }
          if ($frag_count{$circ_id}{$fragment}){
            $frag_count{$circ_id}{$fragment} += $count;
          } else {
            $frag_count{$circ_id}{$fragment} = $count;
            push(@fragment_exon, $fragment);
          }

          #print $circ_id."\t".$fragment."\t".$frag_count{$circ_id}{$fragment}."\t".$contig_frag{$circ_id}{$fragment}."\n";
          #print $fragment."\t".$frag_exon{$fragment}."\n";
          $frag_start = $1 and $frag_end = $2 if $fragment =~ /(\S+)\-(\S+)/;
          ##Determine which exon fragment is located

          for $gene (@circ_genes) {
            next unless $gene =~ /\S+/;
            $exon_cluster = $gene_exon{$gene};
            @exons = split /,/, $exon_cluster;
            for $exon (@exons) {
              $exon_start = $1 and $exon_end = $2 if $exon =~ /(\S+)\-(\S+)/;
              if (($frag_start > ($exon_start - 3)) and ($frag_end < ($exon_end + 3))) {
                $frag_exon{$fragment} = $exon;
                last;
              }
            }

          }
        }

        for $fragment (@fragment_exon){
          #print $circ_id."\t".$fragment."\t".$frag_count{$circ_id}{$fragment}."\t".$frag_exon{$fragment}."\n";
          $frag_start = $1 and $frag_end = $2 if $fragment =~ /(\S+)\-(\S+)/;
          ##Determine which exon fragment is located
          next if $fragment ~~ @fragment_anno;
          if (($last_frag =~ /\S+/) and ($frag_exon{$fragment} =~ /\S+/) and ($frag_exon{$last_frag} =~ /\S+/) and ($frag_exon{$fragment} eq $frag_exon{$last_frag})) {

            my $now_frag_start = $1 and my $now_frag_end = $2 if $fragment =~ /(\S+)\-(\S+)/;
            my $last_frag_start = $1 and my $last_frag_end = $2 if $last_frag =~ /(\S+)\-(\S+)/;
            if ((($last_frag_start < $now_frag_start) or ($last_frag_start == $now_frag_start)) and ($now_frag_start < $last_frag_end) and (($now_frag_end > $last_frag_end) or ($now_frag_end == $last_frag_end))) {
              my $new_frag = $last_frag_start."-".$now_frag_end;
              if (($new_frag ne $last_frag) and ($new_frag ne $fragment)) {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$fragment} + $frag_count{$circ_id}{$last_frag};
                #print $frag_count{$circ_id}{$new_frag}."+=".$frag_count{$circ_id}{$fragment}."+". $frag_count{$circ_id}{$last_frag}."\n";
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$fragment});
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment}.$contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$fragment});
                delete($contig_frag{$circ_id}{$last_frag});
              } elsif ($new_frag ne $last_frag) {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$last_frag};
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$last_frag});
              } else {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$fragment};
                delete($frag_count{$circ_id}{$fragment});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment};
                delete($contig_frag{$circ_id}{$fragment});
              }
              $last_frag = $new_frag;
              $frag_exon{$last_frag} = $frag_exon{$fragment};
            } elsif (($last_frag_end < $now_frag_start) or ($last_frag_end == $now_frag_start)) {
              my $new_frag = $last_frag_start."-".$now_frag_end;
              if (($frag_count{$circ_id}{$fragment} > $frag_count{$circ_id}{$last_frag}) or ($frag_count{$circ_id}{$fragment} == $frag_count{$circ_id}{$last_frag})) {
                $frag_count{$circ_id}{$new_frag} = $frag_count{$circ_id}{$fragment};
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$fragment});
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment}.$contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$fragment});
                delete($contig_frag{$circ_id}{$last_frag});
              } else {
                $frag_count{$circ_id}{$new_frag} = $frag_count{$circ_id}{$last_frag};
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$fragment});
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment}.$contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$fragment});
                delete($contig_frag{$circ_id}{$last_frag});
              }
              $last_frag = $new_frag;
              $frag_exon{$last_frag} = $frag_exon{$fragment};
            } elsif (($now_frag_start < $last_frag_start) and ($now_frag_end > $last_frag_end)) {
              my $new_frag = $fragment;
              $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$last_frag};
              pop(@fragment_anno);
              push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
              delete($frag_count{$circ_id}{$last_frag});
              $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$last_frag};
              delete($contig_frag{$circ_id}{$last_frag});
              $last_frag = $new_frag;
              $frag_exon{$last_frag} = $frag_exon{$fragment};
            }
          } elsif (($last_frag =~ /\S+/) and ($frag_exon{$last_frag} !~ /\S+/) and ($frag_exon{$fragment} !~ /\S+/)) {
            my $now_frag_start = $1 and my $now_frag_end = $2 if $fragment =~ /(\S+)\-(\S+)/;
            my $last_frag_start = $1 and my $last_frag_end = $2 if $last_frag =~ /(\S+)\-(\S+)/;
            if ((abs($now_frag_start - $last_frag_start) < 10) or (abs($now_frag_end - $last_frag_end) < 10)){
              if (($last_frag_start < $now_frag_start) or ($last_frag_start == $now_frag_start)) {
                $new_frag = $last_frag_start."-".$now_frag_end;
              } else {
                $new_frag = $now_frag_start."-".$now_frag_end;
              }
              if (($new_frag ne $last_frag) and ($new_frag ne $fragment)) {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$fragment} + $frag_count{$circ_id}{$last_frag};
                #print $new_frag."\t".$frag_count{$circ_id}{$new_frag}."\n";
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$fragment});
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment}.$contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$fragment});
                delete($contig_frag{$circ_id}{$last_frag});
              } elsif ($new_frag ne $last_frag) {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$last_frag};
                pop(@fragment_anno);
                push(@fragment_anno, $new_frag) unless $new_frag ~~ @fragment_anno;
                delete($frag_count{$circ_id}{$last_frag});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$last_frag};
                delete($contig_frag{$circ_id}{$last_frag});
              } else {
                $frag_count{$circ_id}{$new_frag} += $frag_count{$circ_id}{$fragment};
                delete($frag_count{$circ_id}{$fragment});
                $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$fragment};
                delete($contig_frag{$circ_id}{$fragment});
              }
              $last_frag = $new_frag;
              $frag_exon{$last_frag} = $frag_exon{$fragment};
            } else {
              push(@fragment_anno, $fragment) unless $fragment ~~ @fragment_anno;
              $last_frag = $fragment;
            }
          } else {
            push(@fragment_anno, $fragment) unless $fragment ~~ @fragment_anno;
            $last_frag = $fragment;
          }
        }

        my @frag_anno_contig = ();
        $frag_anno_set = ();
        for my $frag_anno (@fragment_anno) {
          @frag_anno_contig = split /\,/, $contig_frag{$circ_id}{$frag_anno};
          delete ($contig_frag{$circ_id}{$frag_anno});
          my @frag_contig = ();
          for my $contig_id (@frag_anno_contig) {
            next unless $contig_id =~ /\S+/;
            push(@frag_contig, $contig_id) unless $contig_id ~~ @frag_contig;
          }
          for my $contig_id (@frag_contig){
            $contig_frag{$circ_id}{$frag_anno} .= $contig_id.",";
          }
          #print $circ_id."\t".$frag_anno."\t".$frag_count{$circ_id}{$frag_anno}."\t".$frag_exon{$last_frag}."\t".$contig_frag{$circ_id}{$frag_anno}."\n";
          $frag_anno_set .= $frag_anno."/".$frag_count{$circ_id}{$frag_anno}.",";
          #print $frag_anno_set."\n";
          print FRAGANNO "$circ_id\t$frag_anno\/$frag_count{$circ_id}{$frag_anno}\t$contig_frag{$circ_id}{$frag_anno}\n";
        }
        $fragment_anno_all .= $circ_id."\t".$frag_anno_set."\n";
        #print $circ_id."\t".$frag_anno_set."\n";
      } else {
        $frag_anno_set = ();
        @fragment = split /,/, $circ_frag[1];
        $last_frag=();@fragment_exon=();
        for (@fragment){
          next unless $_ =~ /\d+/;
          $frag = $1 and $count = $2 if $_ =~ /(\S+)\/(\S+)/;
          $frag_start = $circ_start + $1 and $frag_end = $circ_start + $2 if $frag =~ /(\S+)\-(\S+)/;

          $fragment = $frag_start."-".$frag_end;
          unless ($fragment eq $frag) {
            $contig_frag{$circ_id}{$fragment} .= $contig_frag{$circ_id}{$frag};
            #print $frag."\t".$contig_frag{$circ_id}{$frag}."\t".$fragment."\t".$contig_frag{$circ_id}{$fragment}."\n";
            delete($contig_frag{$circ_id}{$frag});
          }
          if ($frag_count{$circ_id}{$fragment}){
            $frag_count{$circ_id}{$fragment} += $count;
          } else {
            $frag_count{$circ_id}{$fragment} = $count;
            push(@fragment_exon, $fragment);
          }
          $frag_anno_set .= $fragment."/".$frag_count{$circ_id}{$fragment}.",";
        }
        $fragment_anno_all .= $circ_id."\t".$frag_anno_set."\n";
      }
    }
  } else {
    print "$fragment_result_sort is empty!" and die;
  }

  print 		'[', scalar(localtime), "] Step 4 Finished!\n";
  print 		'[', scalar(localtime), "] Step 5: Generating alternative splicing circRNA.\n";
  print 		'[', scalar(localtime), "]   Generating longest circRNA full length sequences.\n";
  print 		'[', scalar(localtime), "]   Generating alternative splicing circRNA.\n";
  ####delete fragment according frag_count
  my (@fragment_merge, $fragment_result_final);
  if ($fragment_anno_all =~ /\S+/) {
    my @fragment_anno_all = split /\n/, $fragment_anno_all;
    for my $circ_frag (@fragment_anno_all){
      @fragment_merge = ();
      $count=();
      my @circ_frag = split /\t/,$circ_frag;
      $circ_id = $circ_frag[0];
      my @frag_anno = split /\,/, $circ_frag[1];
      for my $frag (@frag_anno){
        $count += 1;
        if ($count == 1){
          $tmp = $frag;
          $tmp = $1 and $tmp_start = $2 and $tmp_end = $3 and $frag_count{$circ_id}{$tmp} = $4 if $tmp =~ /((\S+)-(\S+))\/(\S+)/;
          push (@fragment_merge, $tmp);
        } else {
          $frag = $1 and $frag_start = $2 and $frag_end = $3 and $frag_count{$circ_id}{$frag} = $4 if $frag =~ /((\S+)-(\S+))\/(\S+)/;
          my( $new_start, $new_end );
          if ((($tmp_start + 10) > $frag_start) and (($tmp_end < $frag_end) or ($tmp_end == $frag_end))){
            #print $circ_id."\t".$tmp."\t".$frag."\n";
            ###completely cover
            if (($frag_count{$circ_id}{$tmp} < $frag_count{$circ_id}{$frag}) or ($frag_count{$circ_id}{$tmp} == $frag_count{$circ_id}{$frag})){
              $frag_count{$circ_id}{$frag} += $frag_count{$circ_id}{$tmp};
              delete($frag_count{$circ_id}{$tmp}) unless $tmp eq $frag;
              $contig_frag{$circ_id}{$frag} .= $contig_frag{$circ_id}{$tmp};
              delete($contig_frag{$circ_id}{$tmp}) unless $tmp eq $frag;
              $tmp = $frag;
              pop(@fragment_merge);
              push (@fragment_merge, $tmp);
              $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
            } else {
              $frag_count{$circ_id}{$tmp} += $frag_count{$circ_id}{$frag};
              delete($frag_count{$circ_id}{$frag}) unless $tmp eq $frag;
              $contig_frag{$circ_id}{$tmp} .= $contig_frag{$circ_id}{$frag};
              delete($contig_frag{$circ_id}{$frag}) unless $tmp eq $frag;
              $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
            }
          } elsif (($frag_start > $tmp_start) and ($frag_end < ($tmp_end + 10))) {
            #print $circ_id."\t".$tmp."\t".$frag."\n";
            ###completely cover
            if ($frag_count{$circ_id}{$frag} > $frag_count{$circ_id}{$tmp}){
              $frag_count{$circ_id}{$frag} += $frag_count{$circ_id}{$tmp};
              delete($frag_count{$circ_id}{$tmp}) unless $tmp eq $frag;
              $contig_frag{$circ_id}{$frag} .= $contig_frag{$circ_id}{$tmp};
              delete($contig_frag{$circ_id}{$tmp}) unless $tmp eq $frag;
              $tmp = $frag;
              pop(@fragment_merge);
              push (@fragment_merge, $tmp);
              $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
            } else {
              $frag_count{$circ_id}{$tmp} += $frag_count{$circ_id}{$frag};
              delete($frag_count{$circ_id}{$frag}) unless $tmp eq $frag;
              $contig_frag{$circ_id}{$tmp} .= $contig_frag{$circ_id}{$frag};
              delete($contig_frag{$circ_id}{$frag}) unless $tmp eq $frag;
              $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
            }
          } elsif ((($frag_start > ($tmp_start + 10))) and ($frag_start < $tmp_end) and (($frag_end > ($tmp_end + 10)))){
            #print $circ_id."\t".$tmp."\t".$frag."\n";
            ##part overlap
            #merge two fragments of two fragments both large than 3
            $new_start = $tmp_start;
            $new_end = $frag_end;
            $new_frag = $new_start."-".$new_end;
            $frag_count{$circ_id}{$new_frag} = $frag_count{$circ_id}{$tmp} + $frag_count{$circ_id}{$frag};
            delete($frag_count{$circ_id}{$tmp}) unless $tmp eq $new_frag;
            delete($frag_count{$circ_id}{$frag}) unless $frag eq $new_frag;;
            $contig_frag{$circ_id}{$new_frag} .= $contig_frag{$circ_id}{$tmp}.$contig_frag{$circ_id}{$frag};
            delete($contig_frag{$circ_id}{$tmp}) unless $tmp eq $new_frag;;
            delete($contig_frag{$circ_id}{$frag}) unless $frag eq $new_frag;;
            $tmp = $new_frag;
            pop(@fragment_merge);
            push (@fragment_merge, $tmp);
            $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
          } else {
            $tmp = $frag;
            push (@fragment_merge, $tmp);
            $tmp_start = $1 and $tmp_end = $2 if $tmp =~ /(\S+)-(\S+)/;
          }
        }
      }

      ###delete fragment that was support count less than average - 1, delete if (frag_count < average - 1)
      my ($frag_num, $median, $frag_average, $start_num, $end_num, $frag_new, $aaa, @frag_del_new);
      $frag_num = scalar(@fragment_merge);
      next if $frag_num == 0;
      my $all_num = 0;
      my %contig_exist;
      for $frag (@fragment_merge) {
        my @frag_contig = split /\,/, $contig_frag{$circ_id}{$frag};
        my %contig_exist = ();
        for my $contig (@frag_contig) {
          $contig_exist{$contig} = 1;
        }
        $contig_frag{$circ_id}{$frag} = ();
        for (keys %contig_exist) {
          $contig_frag{$circ_id}{$frag} .= $_.",";
        }
        $all_num += $frag_count{$circ_id}{$frag};
      }
      $median=();
      $frag_average=(); $start_num =(); $end_num=();
      $frag_average = int($all_num / $frag_num);

      for my $frag_new (@fragment_merge){#print $frag_new."\n";
        push( @frag_del_new, $frag_new ) and next unless $frag_count{$circ_id}{$frag_new} < 3;###count >= 3,save
        push( @frag_del_new, $frag_new ) and next if (($frag_count{$circ_id}{$frag_new} > ($frag_average - 2)) or ($frag_count{$circ_id}{$frag_new} == ($frag_average - 2)));
      }

      my( $frag_num, %frag_contig, $isoform1, $isoform2, $i, $j, $one_count, $count, $tmp_frag, $fragment_result );
      $frag_num =  scalar(@frag_del_new);


      $fragment_result = ();
      $circ_start = $2 and $circ_end = $3 if $circ_id =~ /(\S+)\:(\S+)\|(\S+)/;
      for my $frag (@frag_del_new) {
        $count += 1;
        $one_count += 1 if $frag_count{$circ_id}{$frag} == 1;
        if ($count == 1) {
          if ($frag =~ /(\S+)-(\S+)/) {
            $tmp_frag = $circ_start."-".$2;
            $frag_count{$circ_id}{$tmp_frag} = $frag_count{$circ_id}{$frag};
            delete($frag_count{$circ_id}{$frag}) unless $tmp_frag eq $frag;
            $frag_del_new[0] = $tmp_frag unless $tmp_frag eq $frag;
            $frag = $tmp_frag;
          }
        }
        if ($count == $frag_num) {
          if ($frag =~ /(\S+)-(\S+)/) {
            $tmp_frag = $1."-".$circ_end;
            $frag_count{$circ_id}{$tmp_frag} = $frag_count{$circ_id}{$frag};
            delete($frag_count{$circ_id}{$frag}) unless $tmp_frag eq $frag;
            $frag_del_new[-1] = $tmp_frag unless $tmp_frag eq $frag;
            $frag = $tmp_frag;
          }
        }
        $isoform1 .= $frag."/".$frag_count{$circ_id}{$frag}.",";
        $fragment_result .= $frag."\/".$frag_count{$circ_id}{$frag}.",";
      }
      #next if (($one_count * 2) > $frag_num) or (($one_count * 2) == $frag_num);####delete circRNA if more than half fragments count is 1.
      $fragment_result_final .= $circ_id."\t".$fragment_result."\n";

      
      ##if fragment counts bigger than 2, considering alternative splicing
      if ($frag_num > 2) {
        ##generating all possibility isoforms
        for my $frag (@frag_del_new) {
          my @contigs = ();
          @contigs = split /\,/, $contig_frag{$circ_id}{$frag};
          for (@contigs) {
            $frag_contig{$frag} .= $_.",";
          }
        }
        my @isoforms = ();

        for my $frag (@frag_del_new) {
          my @contigs =  split /\,/, $frag_contig{$frag};
          for my $contig (@contigs) {
            my $isoform .= $frag.",";
            for my $inter_frag (@frag_del_new) {
              next if $inter_frag eq $frag;
              my @inter_contigs = split /\t/,$frag_contig{$inter_frag};
              my $item;
              $isoform .= $inter_frag."," if $contig ~~ @inter_contigs;
            }
            push (@isoforms, $isoform);
          }
        }

        my $first_frag = @frag_del_new[0];
        my $last_frag = @frag_del_new[-1];
        my (%isoform_uniq, @fragments_sort, $isoform_sort);
        for my $isoform (@isoforms) {
          my %frag_uniq = ();
          my $isoform_new = ();
          my @fragments = split /,/, $isoform;
          push (@fragments, $first_frag) unless $first_frag ~~ @fragments;
          push (@fragments, $last_frag) unless $last_frag ~~ @fragments;
          @fragments_sort = &sort_frag(@fragments);
          for (@fragments_sort){
            $isoform_new .= $_."\/".$frag_count{$circ_id}{$_}.",";
          }
          $isoform_uniq{$isoform_new} += 1;
        }

        ##deleting the most likely to be deleted fragment
        if ($frag_num % 2 == 0) {
          #fragments count is even
          for $i (1 .. ($frag_num/2 - 1)) {
            $j = ($frag_num - 1) - $i;
            ##fragment_i matched fragment_j
            #fragment_i_coverage - fragment_j_coverage =< 3, fragment_i_coverage < fragment_i-1_coverage, and fragment_i_coverage < fragment_i+1_coverage, delete fragment i
            if (($frag_count{$circ_id}{$frag_del_new[$i]} - $frag_count{$circ_id}{$frag_del_new[$j]}) > 2) {
              if (($frag_count{$circ_id}{$frag_del_new[$j]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$j-1]}) {
                if (($frag_count{$circ_id}{$frag_del_new[$j]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$j+1]}) {
                  #print $circ_id."\t".$frag_del_new[$i]."\n";
                  delete $frag_del_new[$j];
                }
              }
            } elsif (($frag_count{$circ_id}{$frag_del_new[$j]} - $frag_count{$circ_id}{$frag_del_new[$i]}) > 2) {
              if (($frag_count{$circ_id}{$frag_del_new[$i]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$i-1]}) {
                if (($frag_count{$circ_id}{$frag_del_new[$i]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$i+1]}) {
                  #print $circ_id."\t".$frag_del_new[$i]."\n";
                  delete $frag_del_new[$i];
                }
              }
            } else {
              next;
            }
          }
        } else {
          for $i (1 .. (int($frag_num/2) - 1)) {
            $j = ($frag_num - 1) - $i;
            if (($frag_count{$circ_id}{$frag_del_new[$i]} - $frag_count{$circ_id}{$frag_del_new[$j]}) > 2) {
              if (($frag_count{$circ_id}{$frag_del_new[$j]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$j-1]}) {
                if (($frag_count{$circ_id}{$frag_del_new[$j]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$j+1]}) {
                  delete $frag_del_new[$j];
                }
              }
            } elsif (($frag_count{$circ_id}{$frag_del_new[$j]} - $frag_count{$circ_id}{$frag_del_new[$i]}) > 2) {
              if (($frag_count{$circ_id}{$frag_del_new[$i]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$i-1]}) {
                if (($frag_count{$circ_id}{$frag_del_new[$i]} + $difference) < $frag_count{$circ_id}{$frag_del_new[$i+1]}) {
                  delete $frag_del_new[$i];
                }
              }
            } else {
              next;
            }
          }
          ######odd number
          ###coverage_median < (coverage_median-1)/2 and coverage_median < (coverage_median+1)/2, we exclude fragment median in circRNA isoforms
          $median = int($frag_num/2) + 1;
          if (($frag_count{$circ_id}{$frag_del_new[$median]} * 2) <  $frag_count{$circ_id}{$frag_del_new[$median + 1]}) {
            if (($frag_count{$circ_id}{$frag_del_new[$median]} * 2) <  $frag_count{$circ_id}{$frag_del_new[$median - 1]}) {
              delete $frag_del_new[$median];
            }
          }
        }
        $isoform2 = ();
        for (@frag_del_new){
          next unless $_ =~ /\S+/;
          $isoform2 .= $_."/".$frag_count{$circ_id}{$_}.",";
        }
        #print $circ_id."\t".$isoform1."\t".$isoform2."\n";
        unless ($isoform1 eq $isoform2) {
          for (keys %isoform_uniq) {
            $fragment_result_final .= $circ_id."\t".$_."\n" if $_ eq $isoform2;
          }
        }
      } else {
        next;
      }
    }
  } else {
    print "Not find any result!" and die;
  }

  my ($circ_seq, $final_line);
  my @all_res = split /\n/, $fragment_result_final;
  for (@all_res) {
    $final_line += 1;
    my $line = &merge_fragment($_);
    my $total_length = &total_length($line);
    next if $total_length > 5000;
    print FRAGFINAL "$line";
    $circ_seq = &get_seq($line);
    print CIRCSEQ "$circ_seq";
  }
  
  print 		'[', scalar(localtime), "] Step 5 Finished!\n";

  print 		'[', scalar(localtime), "]   Total circRNA isoforms: $final_line.\n";

  print 		'[', scalar(localtime), "] Program Finished ^v^!\n";

  sub circ_contig_read {
    my ( $circ_split, $output ) = @_;
    my ( %contig_uniq, @circ_info, @junc_read_set, @contig_cluster );
    my ( $out1, $out2, $contig_info, $junc_read, $junc_read_id, $contig_id_set, $junc_read_contig_seq, $map_res_directory, $sam, $num );
    $map_res_directory = "$output/map";
    mkdir( $map_res_directory ) or die "can not make a new output directory" unless (-e $map_res_directory);
    print 		'[', scalar(localtime), "]   Start processing $circ_split!\n";
    open FILE2, "$circ_split" or die;
    while( <FILE2> ){
      next if $_ =~ /start/;
      %contig_uniq=(); $contig_info=(); $contig_id_set=();$junc_read_contig_seq=();
      @circ_info = split /\t/,$_;
      my $circ_chr = $circ_info[0];
      my $circ_start = $circ_info[1];
      my $circ_end = $circ_info[2];
      my $circ_strand = $circ_info[3];
      my $circ_id = $circ_chr.":".$circ_start."|".$circ_end;
      my $circ_name = $circ_chr."_".$circ_start."_".$circ_end;
      my $circ_seq = substr( $genome{$circ_chr}, $circ_start - 1, $circ_end - $circ_start + 1 );
      my $circ_seq2 = $circ_seq.$circ_seq;
      $circ_seq2 = &comp_rev($circ_seq2) if $circ_strand =~ /-/;
      my $circ_info_seq = ">" . $circ_id . "\n" . $circ_seq2 . "\n";
      open $out1, ">", "$map_res_directory/$circ_name\.fa" or die;
      print $out1( $circ_info_seq );
      #######junction reads###########
      $junc_read = $circ_info[5];#print $junc_read."\n";
      @junc_read_set = split /,/,$circ_info[5];
      for $junc_read_id (@junc_read_set){
        next unless $junc_read_id =~ /\d+/;
        my $junc_read_id_1 = $junc_read_id."/1";
        my $junc_read_id_2 = $junc_read_id."/2";
        $junc_read_contig_seq .= ">".$junc_read_id_1."\n".$read_seq{$junc_read_id_1}."\n";
        $junc_read_contig_seq .= ">".$junc_read_id_2."\n".$read_seq{$junc_read_id_2}."\n";
        $contig_info .= $read_contig{$junc_read_id}.",";
      }
      ################################
      ##############contigs############
      #print $contig_info."\n";
      @contig_cluster = split /\,/,$contig_info;
      for my $contig_id (@contig_cluster){#print $contig_id."\n";
        $contig_uniq{$contig_id} = 1;
      }
      for my $contig_id (keys %contig_uniq){
        next unless $contig_id =~ /\S+/;
        next if length($contig_seq{$contig_id}) < $read_length;
        $junc_read_contig_seq .= ">".$contig_id."\n".$contig_seq{$contig_id}."\n";
      }

      #################################
      open $out2, ">", "$map_res_directory/$circ_name\_read_contig.fa" or die;##/map/
      print $out2( $junc_read_contig_seq );

      system "bwa index $map_res_directory/$circ_name\.fa";
      system "bwa mem $map_res_directory/$circ_name\.fa $map_res_directory/$circ_name\_read_contig.fa -o $map_res_directory/$circ_name\.sam";
      $sam = "$map_res_directory/$circ_name\.sam";
      &read_sam($sam) and $num += 1 if -s($sam);
      system "rm $map_res_directory/$circ_name*";
    }
    print 		'[', scalar(localtime), "]   Completing the alignment of $num circRNAs in $circ_split !\n";
  }

  sub split_circ_file {
    my ( $circ, $n ) = @_;
    my ( @circ_row, @split_circ, $circ_raw, $input_dir, $circ_row );
    open(FILE1,"$circ") or die;
    @circ_row = <FILE1>;
    $circ_row = @circ_row;
    if (rindex($circ, "/") >= 0) {
      $circ_raw = substr($circ, rindex($circ, "/")+1);
    } else {
      $circ_raw = $circ;
    }
    if ($n == 1){
      if ($circ_row < 1){
        print "Fail to get file row number for $circ.\nFatal error. Aborted.\n";
        die "Fail to get file row number for $circ: $!";
      }
      push(@split_circ, $circ);
    }
    if ($n >= 2) {
      print 		'[', scalar(localtime), "]   Requesting system to split circRNA into $n pieces.\n";
      $input_dir = $output;
      if ($circ_row < 1){
        print "Fail to get file row number for $circ.\nFatal error. Aborted.\n";
        die "Fail to get file row number for $circ: $!";
      }

      my $division;
      if ($circ_row % $n != 0) {
        $division = int($circ_row/$n) + 1;
      } else {
        $division = int($circ_row/$n);
      }

      if(!defined($division) or $division < 0.01*$circ_row){
          print "Fail to calculate divided file row number, total file row number $circ_row, thread $n.\nFatal error. Aborted.\n";
          die "Fail to calculate divided file row number, total file row number $circ_row, thread $n: $!";
      }
      mkdir "$input_dir/split_circ";
      system "split -l $division $circ $input_dir/split_circ/$circ_raw";
      @split_circ = <$input_dir/split_circ/$circ_raw*>;
      my $n2 = 0;
      for my $part (@split_circ) {
        if ($part =~ /$circ_raw[a-z]+$/){
          $n2 ++;
        }
      }
      if ($n2 == $n){
        print 		'[', scalar(localtime), "]   circRNA was divided successfully.\n";
      } else {
        print "Cannot split $circ into $n ($n2) pieces with row of $division and named them as $input_dir/split_circ/$circ_raw.\nFatal error. Aborted.\n";
        die "Cannot split $circ into $n ($n2) pieces with row of $division and named them as $input_dir/split_circ/$circ_raw: $!";
      }
    }
    return @split_circ;
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

  sub sort_frag{
    my @fragment = @_;
    my (@new_fragment, $new_frag, @fragment_sort, @fragment_set);
    for my $frag (@fragment) {
      $new_frag = $2.".".$1 if $frag =~ /(\S+)-(\S+)/;
      push ( @new_fragment, $new_frag);
    }
    @fragment_sort = sort {$a <=> $b} @new_fragment;
    for my $new_frag (@fragment_sort) {
      my $frag = $2."-".$1 if $new_frag =~ /(\S+)\.(\S+)/;
      push (@fragment_set, $frag);
    }
    return @fragment_set;
  }

  sub read_sam{
    my $sam = shift @_;
    my ($fragment_set, $frag_start, $frag_end, $frag_start_another, $frag_end_another, $part_end, $circ_length);
    my (@align_info);

    open SAM, "$sam" or (print "Can not open $sam" and die);
    while(<SAM>){
      chomp;
      $frag_start=(); $frag_end=(); $frag_start_another=(); $frag_end_another=();
      next if $_ =~ /PG/;
      $circ_id = $1 and $circ_length = $2 / 2 and next if $_ =~ /\@SQ\s+SN:(\S+)\s+LN:(\S+)/;
      @align_info = split /\t/,$_;
      next if $align_info[2] =~ /\*/;
      $frag_start = $align_info[3];
      next unless $frag_start =~ /\d+/;
      if ($align_info[5]=~ /(\d+)M(\d+)D(\d+)M/){##alignment result: 60M50D90M
        $frag_end = $frag_start + $1 - 1;
        next if ($frag_end - $frag_start) < 20;##Fragment was ignored if it shorter than 20bp
        $frag_start_another = $frag_end + $2;
        $frag_end_another = $frag_start_another + $3 - 1;
        next if ($frag_end_another - $frag_start_another) < 20;
        ####first fragment location change into circRNA boundary
        if (($frag_end < $circ_length) or ($frag_end == $circ_length)){##fragment located in circRNA boundary
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } elsif (($frag_start < $circ_length) and ($frag_end > $circ_length)) {##fragment covered circRNA junction site
          if (($circ_length - $frag_start) > 10){###fragment which located in first circRNA seqences
            $fragment = $frag_start."-".$circ_length;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
          if (($frag_end - $circ_length + 1) > 10){###fragment which located in second circRNA seqences
            $part_end = $frag_end - $circ_length + 1;
            $fragment = "1-".$part_end;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
        } elsif (($frag_start > $circ_length) or ($frag_start == $circ_length)) {##fragment located in second circRNA sequences
          $frag_start = $frag_start - $circ_length + 1;
          $frag_end = $frag_end - $circ_length + 1;
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } else {
          print "Fragment $frag_start-$frag_end not located in circRNA sequences.\n" and die;
        }
        ####second fragment location change into circRNA boundary
        if (($frag_end_another < $circ_length) or ($frag_end_another == $circ_length)){
          $fragment = $frag_start_another."-".$frag_end_another;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } elsif (($frag_start_another < $circ_length) and ($frag_end_another > $circ_length)) {
          if (($circ_length - $frag_start_another) > 10){
            $fragment = $frag_start_another."-".$circ_length;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;
          }
          if (($frag_end_another - $circ_length + 1) > 10){
            $part_end = $frag_end_another - $circ_length + 1;
            $fragment = "1-".$part_end;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;
          }
        } elsif (($frag_start_another > $circ_length) or ($frag_start_another == $circ_length)) {
          $frag_start_another = $frag_start_another - $circ_length + 1;
          $frag_end_another = $frag_end_another - $circ_length + 1;
          $fragment = $frag_start_another."-".$frag_end_another;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } else {
          print "Fragment $frag_start_another-$frag_end_another not located in circRNA sequences.\n" and die;
        }
      } elsif ($align_info[5]=~ /(\d+)M(\d+)I(\d+)M/) {##alignment result: 60M10I80M
        $frag_end = $frag_start + $1 - 1;
        next if ($frag_end - $frag_start) < 20;##Fragment was ignored if it shorter than 20bp
        $frag_start_another = $frag_end + 1;
        $frag_end_another = $frag_start_another + $3 - 1;
        next if ($frag_end_another - $frag_start_another) < 20;
        ####first fragment location change into circRNA boundary
        if (($frag_end < $circ_length) or ($frag_end == $circ_length)){##fragment located in circRNA boundary
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } elsif (($frag_start < $circ_length) and ($frag_end > $circ_length)) {##fragment covered circRNA junction site
          if (($circ_length - $frag_start) > 10){###fragment which located in first circRNA seqences
            $fragment = $frag_start."-".$circ_length;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
          if (($frag_end - $circ_length + 1) > 10){###fragment which located in second circRNA seqences
            $part_end = $frag_end - $circ_length + 1;
            $fragment = "1-".$part_end;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
        } elsif (($frag_start > $circ_length) or ($frag_start == $circ_length)) {##fragment located in second circRNA sequences
          $frag_start = $frag_start - $circ_length + 1;
          $frag_end = $frag_end - $circ_length + 1;
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } else {
          print "Fragment $frag_start-$frag_end not located in circRNA sequences.\n" and die;
        }
        ####second fragment location change into circRNA boundary
        if (($frag_end_another < $circ_length) or ($frag_end_another == $circ_length)){
          $fragment = $frag_start_another."-".$frag_end_another;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } elsif (($frag_start_another < $circ_length) and ($frag_end_another > $circ_length)) {
          if (($circ_length - $frag_start_another) > 10){
            $fragment = $frag_start_another."-".$circ_length;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;
          }
          if (($frag_end_another - $circ_length + 1) > 10){
            $part_end = $frag_end_another - $circ_length + 1;
            $fragment = "1-".$part_end;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;
          }
        } elsif (($frag_start_another > $circ_length) or ($frag_start_another == $circ_length)) {
          $frag_start_another = $frag_start_another - $circ_length + 1;
          $frag_end_another = $frag_end_another - $circ_length + 1;
          $fragment = $frag_start_another."-".$frag_end_another;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } else {
          print "Fragment $frag_start_another-$frag_end_another not located in circRNA sequences.\n" and die;
        }
      } elsif ($align_info[5]=~ /(\d+)M/){##alignment result: 60M
        $frag_end = $frag_start + $1 - 1;
        next unless $frag_start =~ /\d+/;
        next if ($frag_end - $frag_start) < 20;##Fragment was ignored if it shorter than 20bp
        if (($frag_end < $circ_length) or ($frag_end == $circ_length)){##fragment located in circRNA boundary
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } elsif (($frag_start < $circ_length) and ($frag_end > $circ_length)) {##fragment covered circRNA junction site
          if (($circ_length - $frag_start) > 10){###fragment which located in first circRNA seqences
            $fragment = $frag_start."-".$circ_length;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
          if (($frag_end - $circ_length + 1) > 10){###fragment which located in second circRNA seqences
            $part_end = $frag_end - $circ_length + 1;
            $fragment = "1-".$part_end;
            $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
            $fragment_set .= $fragment.",";
          } else {
            next;##delete fragment shorter than 10bp
          }
        } elsif (($frag_start > $circ_length) or ($frag_start == $circ_length)) {##fragment located in second circRNA sequences
          $frag_start = $frag_start - $circ_length + 1;
          $frag_end = $frag_end - $circ_length + 1;
          $fragment = $frag_start."-".$frag_end;
          $contig_frag{$circ_id}{$fragment} .= $align_info[0].",";
          $fragment_set .= $fragment.",";
        } else {
          print "Fragment $frag_start-$frag_end not located in circRNA sequences.\n" and die;
        }
      } else {
        print "Can not deal with this map information.\n$_\n";
      }
    }
    my @frag = split /\,/, $fragment_set;
    for (@frag){
      print FRAGCONTIG "$circ_id\t$_\t$contig_frag{$circ_id}{$_}\n";
    }
    print FRAGALL "$circ_id\t$fragment_set\n";

  }

  sub get_seq{
    my $line = $_;
    my (@circ_info, @frag_cluster, $circ_chr, $circ_intinal, $circ_id, $frag_start, $frag_end, $circ_sequence, $circ_length, $circ_seq_res);
    @circ_info = split /\t/,$line;
    $circ_id = $circ_info[0];
    $circ_chr = $1 and $circ_intinal = $2 if $circ_info[0] =~ /(\S+)\:(\S+)\|(\S+)/;
    @frag_cluster= split /,/,$circ_info[1];
    for my $frag (@frag_cluster){
      next unless $frag =~ /\-/;
      $frag_start = $1 and $frag_end = $2 if $frag =~ /(\S+)\-(\S+)\/(\S+)/;
      $circ_sequence .= substr($genome{$circ_chr}, $frag_start, $frag_end - $frag_start + 1);
    }
    $circ_length = length($circ_sequence);
    $circ_seq_res = ">".$circ_id."\t".$circ_length."\n".$circ_sequence."\n";
  }

  sub merge_fragment{
    my ($count, $tmp, $frag, $tmp_id, $frag_id, $tmp_start, $frag_start, $tmp_end, $frag_end, @res_frag, %frag_count, $frag_set, $result);
    my @res_line = split /\t/,$_;
    my $circ_id = $res_line[0];
    $frag_set = ();
    my @fragments = split /\,/,$res_line[1];
    for my $frag (@fragments) {
      $count += 1;
      if ($count == 1){
        $tmp = $frag;
        $tmp_id = $1 and $tmp_start = $2 and $tmp_end = $3 and $frag_count{$1} = $4 if $tmp =~ /((\S+)-(\S+))\/(\S+)/;
        push(@res_frag, $tmp_id);
      } else {
        $tmp_id = $1 and $tmp_start = $2 and $tmp_end = $3 and $frag_count{$1} = $4 if $tmp =~ /((\S+)-(\S+))\/(\S+)/;
        $frag_id = $1 and $frag_start = $2 and $frag_end = $3 and $frag_count{$1} = $4 if $frag =~ /((\S+)-(\S+))\/(\S+)/;
        if ($frag_start < $tmp_end) {
          $new_frag = $tmp_start."-".$frag_end;
          $frag_count{$new_frag} = $frag_count{$tmp_id} + $frag_count{$frag_id};
          pop(@res_frag);
          push(@res_frag,$new_frag);
          $tmp = $new_frag."/".$frag_count{$new_frag};
        } elsif (abs($frag_start - $tmp_end) < 10){
          $new_frag = $tmp_start."-".$frag_end;
          if ($frag_count{$tmp_id} > $frag_count{$frag_id}){
            $frag_count{$new_frag} = $frag_count{$tmp_id};
          } else {
            $frag_count{$new_frag} = $frag_count{$frag_id};
          }
          pop(@res_frag);
          push(@res_frag,$new_frag);
          $tmp = $new_frag."/".$frag_count{$new_frag};
        } else {
          $tmp = $frag;
          $tmp_id = $1 and $tmp_start = $2 and $tmp_end = $3 and $frag_count{$1} = $4 if $tmp =~ /((\S+)-(\S+))\/(\S+)/;
          push(@res_frag, $tmp_id);
        }
      }
    }
    for (@res_frag){
      next unless $_ =~ /\S+/;
      $frag_set .= $_."/".$frag_count{$_}.",";
    }
    $result .= $circ_id."\t".$frag_set."\n";
  }

  sub total_length{
    my @res_line = split /\t/,$_;
    my $circ_id = $res_line[0];
    my @fragments = split /\,/,$res_line[1];
    my $total_length = ();
    for my $frag (@fragments) {
      $frag_start = $1 and $frag_end = $2 if $frag =~ /(\S+)-(\S+)\/(\S+)/;
      $total_length += $frag_end - $frag_start + 1;
    }
    return $total_length;
  }
}
