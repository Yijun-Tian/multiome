#!/usr/bin/perl -w
# substr ATAC-seq BAM alignment to obtain read start parts, which contained Tn5 insertion sequencing info
# module load SAMtools/1.16.1-GCC-11.3.0

use strict;
no warnings 'uninitialized';

my $samtools="/app/eb/software/SAMtools/1.16.1-GCC-11.3.0/bin/samtools";
my $deletion_string = "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD";

sub main{
  my $bamfile = $ARGV[0];
  my $outname = $ARGV[1];
  my $id = $ARGV[2];
  my $start = time;
  open(SAM, ">$outname.sam") || die("error writing to file $outname.sam\n");
  open(HEAD, "$samtools view -H $bamfile |") || die("can't read header\n");
  while(my $line = <HEAD>){
  chomp($line);
  print SAM $line,"\n";
  }
  open(INFILE, "$samtools view -D CB:$id $bamfile | " ) || die("Error running samtools to convert BAM file\n");
  while(my $line = <INFILE>){
    chomp($line);
    my %alignment = getAlignmentInfo($line);
    my $uuid = $alignment{"name"};
    my $SEQ = $alignment{"seq"};
    my $QUAL = $alignment{"qual"};
    next if ($alignment{"flag"} & 4 || $alignment{"cigar"} =~ m/\*/);
    getReference(\%alignment);
    my $strand = $alignment{"flag"} & 16 ? "C":"W";
    my $pos;
    my @template;
    my @rev;
    my ($ciglength, $tl);
    if($strand eq "C"){
    next if ($alignment{"cigar"} !~ m/[0-9]+M$/);
    $alignment{"ins"} = substr $SEQ, -29;
    $alignment{"qs"} = substr $QUAL, -29;
    @template = $alignment{"cigar_e"};
    @rev = scalar reverse @template;
    ($ciglength, $tl) = getCIGARlength(@rev);
    $alignment{"cigarn"} = getCIGAR(substr $alignment{"cigar_e"}, -$ciglength);
    $pos = $alignment{"pos"} + $alignment{"tlen"} + $alignment{"ref_nm_i"} - $tl;
    $alignment{"ref_new"} = $pos >=1 ? $pos : "1" ;
    }
    else{
    next if ($alignment{"cigar"} !~ m/^[0-9]+M/);
    $alignment{"ins"} = substr $SEQ, 0, 29;
    $alignment{"qs"}  = substr $QUAL, 0, 29;
    ($ciglength, $tl) = getCIGARlength($alignment{"cigar_e"}); 
    $alignment{"cigarn"} = getCIGAR(substr $alignment{"cigar_e"}, 0, $ciglength);
    $pos = $alignment{"pos"} ;
    $alignment{"ref_new"} = $pos >=1 ? $pos : "1" ;
    }
    #print $alignment{"cigar"},"\t",$tl,"\t",$alignment{"tlen"},"\n";
#
    print SAM $alignment{"name"},"\t",
              $alignment{"flag"},"\t",
              $alignment{"chr"},"\t",
              $alignment{"ref_new"},"\t",
              $alignment{"mapq"},"\t",
              $alignment{"cigarn"},"\t",
	      $alignment{"rnext"},"\t",
	      $alignment{"pnext"},"\t",
	      $alignment{"tlength"},"\t",
              $alignment{"ins"},"\t",
              $alignment{"qs"},"\t",
              $alignment{"tags"},"\n";
    }
close(SAM);
}

sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = $tmp[0];
	$alignment{"flag"} = $tmp[1];
	$alignment{"chr"} = $tmp[2];
	$alignment{"pos"} = $tmp[3];
	$alignment{"mapq"} = $tmp[4];
	$alignment{"cigar"} = $tmp[5];
	$alignment{"rnext"} = $tmp[6];
	$alignment{"pnext"} = $tmp[7];
	$alignment{"tlength"} = $tmp[8];
	$alignment{"seq"} = $tmp[9];
	$alignment{"qual"} = $tmp[10];
        $alignment{"tags"} = splice(@tmp, 11);        
	return %alignment;
}

sub getCIGARlength{
	my $template = shift;
	my $i = 0;
	my $count = 0;
	my $tl = 0;

	while($i < length($template)){     
		my $current = substr $template, $i, 1;
		$count++ if($current !~ m/D/);
		$tl++ if($current =~ m/[MD]/);
		$i++;
		last if($count == 29);
		}
	return $i, $tl;
}

sub getCIGAR{
	my $string = shift;
	my $i = 0;
	my $nCIGAR;
	my $count = 0;
	while($i < length($string)){
	my $exp;
	my $current = substr($string, $i, 1);
	my $next = substr($string, $i+1, 1);
	if($current eq $next){
	$count++;		
	}
	else{
	$exp = $count + 1 . $current;
	$count = 0; 
	}
	$nCIGAR = $nCIGAR . $exp; 
	$i++;
	}
	return $nCIGAR;
}

sub getReference{
	my $alignment = shift;
	my $cigar = $alignment->{"cigar"};
	my $query = $alignment->{"seq"};
	my $qual = $alignment->{"qual"};
	my $ref_seq;
	my $ref_qual;
	my $cigar_e;
	# Adapted from Ben Langmead (BtlBio::Alignment:Util.pm)
	# CIGAR fields are always in pairs
	my $i = 0;
	my $j = 0;
	my $nm_i = 0;
	my $nm_d = 0;
	while($i < length($cigar)){
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die("Bad cigar string : '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || die("count not parse operation at pos $i: '$cigar'");
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=DIS]/);
		my $exp = $op x $runlen;
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s =~ s/[ATGCatgc]/I/g if($op eq "I"); 
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}else{
			#deletion from reference
			$nm_d += $runlen;
			length($deletion_string) > $runlen || die("deletion is too long at $runlen: '$cigar'");
			$clip_s = substr($deletion_string, 0, $runlen);
			$clip_q = substr($deletion_string, 0, $runlen);
		}
		$i++;
		$ref_seq = $ref_seq . $clip_s if($op =~ m/[MD]/);
		$ref_qual = $ref_qual . $clip_q if($op =~ m/[MD]/);
		$cigar_e = $cigar_e . $exp; 
	}
	
	$alignment->{"ref_match_seq"} = $ref_seq;
	$alignment->{"ref_match_qual"} = $ref_qual;
	$alignment->{"ref_nm_i"} = $nm_i;
	$alignment->{"ref_nm_d"} = $nm_d;
	$alignment->{"tlen"} = length($ref_seq) - $nm_i if($ref_seq);
	$alignment->{"cigar_e"} = $cigar_e;
}

main();

