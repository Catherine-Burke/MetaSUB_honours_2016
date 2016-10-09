package Phylosift::Command::demux;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use JSON;
use Carp;
use Phylosift::Utilities qw(debug ps_open miss);
use IO::Zlib;

use strict;
use warnings;

our $VERSION = "v1.0.1";
my $USEARCH = "/home/koadman/Downloads/usearch7.0.1001_i86linux32";
my $CLUSTER_IDENTITY = "0.89";
my $FLASH = "flash";
my $RANDOM_LEN;

sub description {
	return
	  "phylosift demux - process a barcoded Illumina run into a single interleaved file. Input is a directory of FastQ files generated by Illumina's pipeline."
	  ;
}

sub abstract {
	return "process a barcoded Illumina run into a single interleaved file. Input is a directory of FastQ files generated by Illumina's pipeline.";
}

sub usage_desc { "demux %o <illumina directory>" }

sub options {
	return (
			 [ "sample-map=s",    "A file with the list of samples names, their barcode sequences, and any sample metadata", { required => 1 } ],
			 [ "barcode-table=s", "A table containing barcode and adapter information",                                      { required => 1 } ],
			 [ "output=s",        "Name of the output file",                                                                 { required => 1 } ],
			 [ "barcode-pos=i",  "Position in read(s) where barcode exists (zero-based)", { required => 1, default => 0 } ],
			 [ "stats-output=s", "Name of an output file for sample summary statistics" ],
			 [ "flash",          "Assemble overlapping read pairs into a single longer read with FLASH" ],
			 [ "samplefiles",    "Output each sample separately" ],
			 [ "no-interleaving", "Save reads 1 and 2 in separate files instead of interleaving them" ],
			 [ "rev-barcode=i",  "The specified barcode is read in the reverse complement direction" ],
			 [ "swap-barcodes",  "Swap the first and second barcodes in the sample file" ],
			 [ "cluster",         "Perform clustering of randomized barcode tags" ],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
}

sub load_opt {
	my %args = @_;
	my $opt  = $args{opt};
	$Phylosift::Settings::configuration        = $opt->{config};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::my_debug             = $opt->{debug};

	$Phylosift::Utilities::debuglevel = $Phylosift::Settings::my_debug || 0;
}

my $include_numeric_read_id = 0;

sub execute {
	my ( $self, $opt, $args ) = @_;
	load_opt( opt => $opt );
	Phylosift::Command::sanity_check();
	Phylosift::Utilities::program_checks();
	my $ps;
	Phylosift::Utilities::data_checks( self => $ps );

	croak("Please provide the path to a directory of Illumina fastq files") unless defined @$args[0];

	my %sample_map = read_sample_file( sample_file => $opt->{sample_map}, rev_barcode => $opt->{rev_barcode}, swap_barcodes => $opt->{swap_barcodes} );
	my %barcode_table = read_barcodes( barcode_file => $opt->{barcode_table}, rev_barcode => $opt->{rev_barcode} );
	my $barcode_length;
	foreach my $read(keys(%barcode_table)){
		foreach my $bc(keys(%{$barcode_table{$read}})){
			$barcode_length = length($bc) unless defined($barcode_length);
			die "Error barcodes have unequal lengths" if $barcode_length != length($bc);
		}
	}
	print "Barcode length is $barcode_length\n";
	my $illumina_dir  = @$args[0];
	my @files         = <$illumina_dir/*_R1_001.fastq.gz>;
	croak(
		"Unable to locate any Illumina FastQ files in $illumina_dir\nPlease check that the directory contains FastQ files named according to Illumina's convention, e.g. names ending with R1_001.fastq.gz, R2_001.fastq.gz, etc"
	) unless @files > 0;

	my %outstream;
	if( $opt->{samplefiles} ) {
		foreach my $barcode(keys(%sample_map)){
			my $sample = $sample_map{$barcode};
			if ( $opt->{flash} ) {
				$outstream{$barcode}{both} = ps_open("| $FLASH -I -z -o $sample --max-overlap=250 --min-overlap=15 -");
			}elsif( $opt->{no_interleaving} ){
				$outstream{$barcode}{r1} = new IO::Zlib;
				$outstream{$barcode}{r2} = new IO::Zlib;
				$outstream{$barcode}{r1}->open( "$sample.r1.fastq.gz", "wb9" );
				$outstream{$barcode}{r2}->open( "$sample.r2.fastq.gz", "wb9" );
			}else{
				$outstream{$barcode}{both} = new IO::Zlib;
				$outstream{$barcode}{both}->open( "$sample.fastq.gz", "wb9" );
			}
		}
	} else {
		my $gzout_name = $opt->{output};
		if ( $opt->{flash} ) {
			$outstream{both} = ps_open("| $FLASH -I -z -o $gzout_name --max-overlap=250 --min-overlap=15 -");
		}elsif( $opt->{no_interleaving} ){
			$outstream{r1} = new IO::Zlib;
			$outstream{r2} = new IO::Zlib;
			$outstream{r1}->open( "$gzout_name.r1.fastq.gz", "wb9" );
			$outstream{r2}->open( "$gzout_name.r2.fastq.gz", "wb9" );
		}else{
			$outstream{both} = new IO::Zlib;
			$gzout_name .= ".fastq.gz" unless $gzout_name =~ /\.gz$/;
			$outstream{both}->open( $gzout_name, "wb9" );
		}
	}

	foreach my $file (@files) {
		print STDERR "Processing $file\n";
		$file =~ m/^(\S+)_\S\S_(\d+).fastq.gz/;
		process_barcodes(
						  barcode     => \%barcode_table,
						  sample_map  => \%sample_map,
						  gzout       => \%outstream,
						  in_file     => $file,
						  barcode_pos => $opt->{barcode_pos},
						  barcode_len => $barcode_length,
						  cluster => $opt->{cluster},
						  output_base => $opt->{output}
		);
	}
}

my $r1_defined = 0;
my $r2_defined = 0;

sub process_barcodes {
	my %args        = @_;
	my $barcode     = $args{barcode} || miss("barcode");
	my $sample_map  = $args{sample_map} || miss("sample_map");
	my $gzout       = $args{gzout} || miss("gzout");
	my $in_file     = $args{in_file} || miss("in_file");
	my $barcode_pos = $args{barcode_pos};
	my $barcode_len = $args{barcode_len} || miss("barcode_len");
	my $output_base = $args{output_base};

	$in_file =~ m/^(\S+)_\S\S_(\d+).fastq.gz/;
	my $core  = $1;
	my $index = $2;

	# logic: try four-read with names R1,I1,I2,R2, then R1,R2,R3,R4, then R1,R2,R3, finally R1,R2 w/inline barcodes
	my $i1_file = "$core"."_I1_$index.fastq.gz";
	$i1_file = "$core"."_R2_$index.fastq.gz" unless -e $i1_file;
	$i1_file = undef unless -e $i1_file;
	my $i2_file = "$core"."_I2_$index.fastq.gz";
	$i2_file = "$core"."_R3_$index.fastq.gz" unless -e $i2_file;
	$i2_file = undef unless -e $i2_file;
	my $r2_file = "$core"."_R4_$index.fastq.gz";
	$r2_file = "$core"."_R3_$index.fastq.gz" unless -e $r2_file;
	$r2_file = "$core"."_R2_$index.fastq.gz" unless -e $r2_file;
	$i2_file = undef if $i1_file eq $r2_file;
	$i1_file = undef if $i1_file eq $r2_file;
	print STDERR "i1_file $i1_file\ni2_file $i2_file\n" if defined($i1_file);
	print STDERR "r2_file $r2_file\n\n";

	my %STREAMS;
	$STREAMS{1} = ps_open("gzip -cd $in_file |");
	$STREAMS{2} = ps_open("gzip -cd $i1_file |") if defined $i1_file;
	$STREAMS{3} = ps_open("gzip -cd $i2_file |") if defined $i2_file;
	$STREAMS{4} = ps_open("gzip -cd $r2_file |");

	my $process_limit = 99999999999999;
#	my $process_limit = 1000000;
	my %bc_counts;
	my %rbc_counts;
	my %sample_counts;
	my $counter  = 0;
	my $defcount = 0;
	my $scount   = 0;
	my $READ_CLUSTER_LEN = 32 - $RANDOM_LEN - $barcode_len;

	# if we've done nextera-style dual barcodes then the barcode coupled to read 1 comes in on the 3rd read
	my ($rbc1id, $rbc2id) = get_random_bc_read_id(i1_file=>$i1_file, i2_file => $i2_file);
	my $RSCLUSTFA = ps_open(">$output_base.random_seq.fa") if defined($args{cluster});
	my $RSR1FA = ps_open(">$output_base.random_read1.fa") if defined($args{cluster});
	my $RSR2FA = ps_open(">$output_base.random_read2.fa") if defined($args{cluster});
	while (1) {
		$counter++;
		last if $counter > $process_limit;
		my %bcdata =
		  get_barcoded_read(
							 STREAMS       => \%STREAMS,
							 barcode       => $barcode,
							 barcode_pos   => $barcode_pos,
							 barcode_len   => $barcode_len,
							 sample_map    => $sample_map,
							 sample_counts => \%sample_counts
		  );
		last if !defined($bcdata{reads});
		my $readname = $bcdata{reads}->{1}[0];
		$readname =~ s/^@/>/g;
		$readname =~ s/ .+//g;

		if(defined($args{cluster})){
			if($bcdata{r1d} && $bcdata{rbc}->[$rbc1id] ne ""){
				print $RSR1FA $readname;
				print $RSR1FA $bcdata{bc}->[0].$bcdata{rbc}->[$rbc1id].substr($bcdata{reads}->{1}[1], 0, $READ_CLUSTER_LEN)."\n";
				$r1_defined++;
			}
			if($bcdata{r2d} && $bcdata{rbc}->[$rbc2id] ne ""){
				print $RSR2FA $readname;
				print $RSR2FA $bcdata{bc}->[1].$bcdata{rbc}->[$rbc2id].substr($bcdata{reads}->{4}[1], 0, $READ_CLUSTER_LEN)."\n";
				$r2_defined++;
			}
		}
		if ($bcdata{defined}) {
			if(defined($args{cluster}) && $bcdata{sample}){
				print $RSCLUSTFA $readname;
				print $RSCLUSTFA $bcdata{bc}->[0].$bcdata{rbc}->[$rbc1id].substr($bcdata{reads}->{1}[1], 0, $READ_CLUSTER_LEN);
				print $RSCLUSTFA $bcdata{bc}->[1].$bcdata{rbc}->[$rbc2id].substr($bcdata{reads}->{4}[1], 0, $READ_CLUSTER_LEN)."\n";
			}
			my $jbc  = join( ":", @{$bcdata{bc}} );
			my $jrbc = join( ":", @{$bcdata{rbc}} );
			$bc_counts{$jbc} = 0 unless defined( $bc_counts{$jbc} );
			$bc_counts{$jbc}++;
			$rbc_counts{$jbc}{$jrbc} = 0 unless defined( $rbc_counts{$jbc} ) && defined( $rbc_counts{$jbc}{$jrbc} );
			$rbc_counts{$jbc}{$jrbc}++;
			$defcount++;

			my $r1_out = defined($gzout->{both}) ? $gzout->{both} : $gzout->{r1};
			my $r2_out = defined($gzout->{both}) ? $gzout->{both} : $gzout->{r2};
			$r1_out = defined($r1_out) ? $r1_out : defined($gzout->{$jbc}{both}) ? $gzout->{$jbc}{both} : $gzout->{$jbc}{r1};
			$r2_out = defined($r2_out) ? $r2_out : defined($gzout->{$jbc}{both}) ? $gzout->{$jbc}{both} : $gzout->{$jbc}{r2};
			print $r1_out @{ $bcdata{reads}->{1} } if defined $r1_out;
			print $r2_out @{ $bcdata{reads}->{4} } if defined $r2_out;
		}
		$scount += $bcdata{sample};
		print "$scount / $defcount / $counter reads: matching samples / matching barcodes / total\n" if ( $counter % 1000 == 0 );
	}
	if(defined($args{cluster})){
		close $RSR1FA;
		close $RSR2FA;
		close $RSCLUSTFA;
	}
	foreach my $STRM(values(%STREAMS)){
		close $STRM;
	}

	print "r1 defined $r1_defined\tr2 defined $r2_defined\n";

	if(defined($args{cluster})){
		cluster_molecular_tags(output_base => $output_base, in_file => $in_file, i1_file => $i1_file, i2_file => $i2_file,
r2_file => $r2_file, barcode_pos => $barcode_pos, barcode_len => $barcode_len, process_limit => $process_limit, barcode => $barcode, 
sample_map => $sample_map, random_barcode_len => $RANDOM_LEN);
	}

	my $maxbc;
	open(RSOUT, ">$output_base.random_stats.txt");
	foreach my $jbc ( sort { $bc_counts{$a} <=> $bc_counts{$b} } ( keys %bc_counts ) ) {
		print "$sample_map->{$jbc}\t$bc_counts{$jbc}\n" if defined($sample_map->{$jbc});
		print "$jbc\t$bc_counts{$jbc}\n" unless defined($sample_map->{$jbc});
		$maxbc = $jbc;

		foreach my $jrbc(sort {$rbc_counts{$maxbc}{$a} <=> $rbc_counts{$maxbc}{$b}} (keys %{$rbc_counts{$maxbc}})){
			print RSOUT "$maxbc\t$jrbc\t$rbc_counts{$maxbc}{$jrbc}\n";
		}
	}
}

sub get_random_bc_read_id {
	my %args        = @_;
	my $i1_file     = $args{i1_file};
	my $i2_file     = $args{i2_file};
	my $rbc1id = defined $i1_file && defined $i2_file ? 3 : 1;
	my $rbc2id = defined $i1_file && defined $i2_file ? 2 : 4;
	return ($rbc1id, $rbc2id);
}

sub cluster_molecular_tags {
	my %args        = @_;
	my $output_base = $args{output_base};
	my $in_file     = $args{in_file};
	my $i1_file     = $args{i1_file};
	my $i2_file     = $args{i2_file};
	my $r2_file     = $args{r2_file};
	my $barcode_pos = $args{barcode_pos};
	my $barcode_len = $args{barcode_len};
	my $process_limit = $args{process_limit};
	my $barcode = $args{barcode};
	my $sample_map = $args{sample_map};
	my $random_barcode_len = $args{random_barcode_len} || miss("random_barcode_len");

	# now cluster the random tags
	# cluster separately for r1, r2, and both together
	my $usearch = "$USEARCH -threads 4  -cluster_fast $output_base.random_seq.fa -id $CLUSTER_IDENTITY -sizeout -consout $output_base.random_clust.fa";
	system($usearch);
	my $usearch_r1 = "$USEARCH -threads 4  -cluster_fast $output_base.random_read1.fa -id $CLUSTER_IDENTITY -sizeout -consout $output_base.random_clust1.fa";
	system($usearch_r1);
	my $usearch_r2 = "$USEARCH -threads 4  -cluster_fast $output_base.random_read2.fa -id $CLUSTER_IDENTITY -sizeout -consout $output_base.random_clust2.fa";
	system($usearch_r2);

	# read the full length template random tags back in and
	# split them in half into the random tags found on each end of the template
	my $RCLUSTBOTH = ps_open("$output_base.random_clust.fa");
	my %rtag_combos1;
	my %rtag_combos2;
	my $count;
	my $name;
	my %groupings;
	my %combo_counts;
	my $seq = "";
	while(my $line = <$RCLUSTBOTH>){
		if($line =~ /^>(.+);size=(\d+)/){
			if(length($seq) > 0){
				my $half_len = length($seq)/2;
				my $r1_tag = substr($seq, 0, $half_len);
				my $r2_tag = substr($seq, $half_len, $half_len);
				$rtag_combos1{$r1_tag}{$count} = $name;
				$rtag_combos2{$r2_tag}{$count} = $name;
	
				my $rbc1 = substr($seq,$barcode_len, $half_len - $barcode_len);
				my $rbc2 = substr($seq,$half_len+$barcode_len);
				my $cname = $name;
				$cname = substr($cname, 1); # trim off leading >
				chomp $cname;
				$combo_counts{1}{$r1_tag}{$cname} = $count; 
				$combo_counts{2}{$r2_tag}{$cname} = $count;
				$groupings{$cname}{rbc1} = $rbc1;
				$groupings{$cname}{rbc2} = $rbc2;
				$groupings{$cname}{r1_tag} = $r1_tag;
				$groupings{$cname}{r2_tag} = $r2_tag;
				$groupings{$cname}{count} = $count;
				$groupings{$cname}{name} = $name;
				$seq = "";
			}

			$name = $line;
			$count = $2;
		}else{
			chomp $line;
			$seq .= $line;
		}
	}
	close $RCLUSTBOTH;
	
	# identify the true templates as those with the most frequently observed combination of each individual barcode 
	# those templates containing a barcode that also occurs in a more abundant template are likely recombinants
	# then write out the tags found in valid full length templates to a representatives file
	# create one of these files for each tag
	my $recombinant = 0;
	open(R1CR, ">$output_base.random_clust_rep1.fa");
	open(R2CR, ">$output_base.random_clust_rep2.fa");
	open(CLUSTERSIZES, ">$output_base.cluster_sizes.txt");
	print CLUSTERSIZES "Cluster_name\tCount\tMax with r1\tMax with_r2\tRecombinant1\tRecombinant2\n";
	foreach my $group(keys(%groupings)){
		my $r1_max = $combo_counts{1}{$groupings{$group}{r1_tag}}{$group} == max(values($combo_counts{1}{$groupings{$group}{r1_tag}}));
		my $r2_max = $combo_counts{2}{$groupings{$group}{r2_tag}}{$group} == max(values($combo_counts{2}{$groupings{$group}{r2_tag}}));
		my $is_max = $r1_max && $r2_max; 
		$groupings{$group}{valid} = $is_max;
		$recombinant++ unless $is_max;
		print CLUSTERSIZES "$group\t".$combo_counts{1}{$groupings{$group}{r1_tag}}{$group}."\t".max(values($combo_counts{1}{$groupings{$group}{r1_tag}}))."\t".max(values($combo_counts{2}{$groupings{$group}{r2_tag}}))."\t".($r1_max?0:1)."\t".($r2_max?0:1)."\n";
		next unless $groupings{$group}{valid};
		print R1CR "$groupings{$group}{name}$groupings{$group}{r1_tag}\n";
		print R2CR "$groupings{$group}{name}$groupings{$group}{r2_tag}\n";
	}
	close CLUSTERSIZES;
	close R1CR;
	close R2CR;
	print STDERR "Found $recombinant / ".scalar(keys(%groupings))." potentially recombinant full length templates\n";
	
	# find tags that exist only in R1 and R2 and add them to the representatives
	# then match all reads against those tags to define read groupings
	# quack is a hash mapping individual read names to group names
	my $quack = add_and_match_representatives(rid=>1, output_base => $output_base);
	my $quack2 = add_and_match_representatives(rid=>2, output_base => $output_base);
	
	my $min_read_pairs = 3; # don't write a cluster unless it has at least this many sequences in the chunk
	my $clusters = 0;
	`mkdir -p $output_base.clusters`;

	my $counter = 1;
	my %STREAMS;
	my %sample_counts;
	my $chunk_size = 1000000; # process reads in batches of this size to avoid running out of memory
	my ($rbc1id, $rbc2id) = get_random_bc_read_id(i1_file=>$i1_file, i2_file => $i2_file);
	$STREAMS{1} = ps_open("gzip -cd $in_file |");
	$STREAMS{2} = ps_open("gzip -cd $i1_file |") if defined $i1_file;
	$STREAMS{3} = ps_open("gzip -cd $i2_file |") if defined $i2_file;
	$STREAMS{4} = ps_open("gzip -cd $r2_file |");
	print STDERR "rbc1id $rbc1id rbc2id $rbc2id\n";
	while (1) {
		my %quick; # stores the reads, grouped together by template molecule
		my %bcdata;
		while($counter % $chunk_size > 0){
			$counter++;
			last if $counter > $process_limit;
			%bcdata = get_barcoded_read(
								 STREAMS       => \%STREAMS,
								 barcode       => $barcode,
								 barcode_pos   => $barcode_pos,
								 barcode_len   => $barcode_len,
								 sample_map    => $sample_map,
								 sample_counts => \%sample_counts
			  );
			last if !defined($bcdata{reads});
			my $readname = $bcdata{reads}->{1}[0];
			$readname =~ s/^@//g;
			$readname =~ s/ .+//g;
			chomp $readname;
			if(defined($quack->{$readname}) && defined($groupings{$quack->{$readname}})){
				$bcdata{reads}->{1}[0] =~ s/:/:r1pair:/; # tag the read with the side that it belongs on
				$bcdata{reads}->{4}[0] =~ s/:/:r1pair:/;
				$quick{$quack->{$readname}}{$readname}{reads} = $bcdata{reads};
				$quick{$quack->{$readname}}{$readname}{side} = 1;
			}
			if(defined($quack2->{$readname}) && defined($groupings{$quack2->{$readname}})){
				$bcdata{reads}->{1}[0] =~ s/:/:r2pair:/;
				$bcdata{reads}->{4}[0] =~ s/:/:r2pair:/;
				$bcdata{reads}->{1}[0] =~ s/:r2pair:r1pair:/:/; # untag it if it's a LEnd+REnd read!
				$bcdata{reads}->{4}[0] =~ s/:r2pair:r1pair:/:/;
				$quick{$quack2->{$readname}}{$readname}{reads} = $bcdata{reads};
				$quick{$quack2->{$readname}}{$readname}{side} = 2;
			}
		}
		print STDERR "Processing ".scalar(keys(%quick))." read groups\n";
		$counter++;

		foreach my $readgroup(keys(%quick)){
			my $cname = $readgroup;		
			$cname =~ s/\;size.+//g;
			$cname =~ s/\:/_/g;
			$cname =~ s/centroid\=//g;
			$cname =~ s/\;seqs=\d+//g;
			next unless scalar(keys(%{$quick{$readgroup}})) > $min_read_pairs;
	
			$cname .= ".".$groupings{$readgroup}{rbc1}.".".$groupings{$readgroup}{rbc2};
#			print STDERR "Processing $cname\n";
			`mkdir -p $output_base.clusters/$cname`;
			open(CLUSTER, ">>$output_base.clusters/$cname/$cname.fq");
			my $side1 = 0; my $side2 = 0;
			foreach my $read(keys(%{$quick{$readgroup}})){
				next if $read eq "rbc1" || $read eq "rbc2";
				my $reads = $quick{$readgroup}{$read}{reads};
				print CLUSTER @{ $reads->{1} } if defined($reads->{1});
				print CLUSTER @{ $reads->{4} } if defined($reads->{4});
				$side1++ if $quick{$readgroup}{$read}{side} == 1;
				$side2++ if $quick{$readgroup}{$read}{side} == 2;
			}
			close CLUSTER;
			unless(-x "$output_base.clusters/$cname/$cname.lib"){
				open(CLIB, ">$output_base.clusters/$cname/$cname.lib");
				print CLIB "[LIB]\nshuf=$cname.fq\nins=600\n";
				close CLIB;
				$clusters++;
			}
		}
		last if !defined($bcdata{reads});
	}
}


sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

sub write_cluster_reps {
	my %args          = @_;
	my $rtag_combos = $args{rtag_combos};
	my $filename = $args{filename};
	my $groupings = $args{groupings};
	my $RCLUST = ps_open(">$filename");	
	my $recomb = 0;
	foreach my $tag(keys(%{$rtag_combos})){
		my $iter = 0;
		foreach my $name ( sort { $b <=> $a } (keys %{$rtag_combos->{$tag}}) ) {
			print $RCLUST "$rtag_combos->{$tag}{$name}$tag\n";
			$recomb++ if $iter > 0;
			$iter++;
		}
	}
	close $RCLUST;
	return $recomb;
}

sub add_and_match_representatives {
	my %args          = @_;
	my $rid = $args{rid};
	my $output_base = $args{output_base};
	
	# search clustered tags from just this read against the tags found in full length templates
	my $usearch = "$USEARCH -threads 4  -usearch_global $output_base.random_clust$rid.fa -db $output_base.random_clust_rep$rid.fa -strand plus -id $CLUSTER_IDENTITY -uc $output_base.readmap_rep_r$rid.uc";
	system($usearch);
	my $RUC = ps_open("$output_base.readmap_rep_r$rid.uc");
	# collect a list of tags that were not found in full length templates -- these will be added
	my %keepers;
	while(my $line = <$RUC>){
		chomp $line;
		my @dat = split(/\t/, $line);
		next unless $dat[9] eq "*";
		$keepers{$dat[8]} = 1;
	}
	close $RUC;
	
	# add any new tags that were not found in the full length templates
	my $RCLUSTREP = ps_open(">>$output_base.random_clust_rep$rid.fa");
	my $RCLUST = ps_open("$output_base.random_clust$rid.fa");
	my $name;
	while(my $line = <$RCLUST>){
		if($line =~ /^>(.+)\n/){
			$name = $1;
		}elsif(defined($keepers{$name})){
			print $RCLUSTREP ">$name\n$line";
		}
	}
	close $RCLUSTREP;
	close $RCLUST;

	# match all the reads against the complete set of random tags on this side
	$usearch = "$USEARCH -threads 4  -usearch_global $output_base.random_read$rid.fa -db $output_base.random_clust_rep$rid.fa -strand plus -id $CLUSTER_IDENTITY -uc $output_base.readmap_r$rid.uc";
	system($usearch);

	# now parse reads by tag & store in %quack for later assembly
	my %quack;
	$RUC = ps_open("$output_base.readmap_r$rid.uc");
	while(my $line = <$RUC>){
		chomp $line;
		my @dat = split(/\t/, $line);
		next if $dat[9] eq "*";
		$quack{$dat[8]} = $dat[9];
	}
	close $RUC;
	return \%quack;
}

sub get_barcoded_read {
	my %args          = @_;
	my $STREAMS       = $args{STREAMS} || miss("STREAMS");
	my $barcode       = $args{barcode} || miss("barcode");
	my $barcode_pos   = $args{barcode_pos};
	my $barcode_len   = $args{barcode_len} || miss("barcode_len");
	my $sample_map    = $args{sample_map} || miss("sample_map");
	my $sample_counts = $args{sample_counts} || miss("sample_counts");

	my %reads;
	for ( my $i = 1; $i < 5; $i++ ) {
		$reads{$i} = ();
		for ( my $j = 0; $j < 4; $j++ ) {
			my $si = $STREAMS->{$i};
			$reads{$i}[$j] = <$si> if defined $si;
			return if ( $i == 1 && !defined( $reads{$i}[$j] ) );
		}
	}
	my @bc = ("","");
	my @rbc = ("","","","","");
	my $defined_bc = 1;    # 1 if all barcodes were defined in the table
	my $sample_bc  = 0;    # 1 if the combination of barcodes matches a sample
	my $linker     = 1;    # 1 if the linker sequence was found in both reads
	my $r1d        = 0;
	my $r2d        = 0;
	foreach my $rid ( sort {$a <=> $b} keys(%$barcode) ) {
		my $b = substr( $reads{$rid}[1], $barcode_pos, $barcode_len );
		
#		my @qual = split( //, substr( $reads{$rid}[3], $barcode_pos, $barcode_len ) );
#		my $qavg = 0;
#		foreach my $qq(@qual){ $qavg += ord($qq); }
#		$qavg /= @qual;
#		print STDERR "qavg too low: $qavg\n" if $qavg <= 60;
#		next unless $qavg > 60;

		# try to error-correct the barcode
		my $ecb = $b;
		$ecb = $barcode->{$rid}{$b}{parent} if defined $barcode->{$rid}{$b};
#		print "barcode $b parent .".$barcode->{$rid}{$b}{parent}."\n"; # if $rid == 4;
		$defined_bc = 0 unless defined $barcode->{$rid}{$ecb};
		$r1d = $rid if defined $barcode->{$rid}{$ecb} && $rid == 1 || $rid == 3;
		$r2d = $rid if defined $barcode->{$rid}{$ecb} && $rid == 4 || $rid == 2;
		$bc[0] = $ecb if $r1d == $rid;
		$bc[1] = $ecb if $r2d == $rid;
#		print "read $rid $ecb undefined\n" unless $defined_bc; # || $rid != 4;

		# now capture any randomized barcodes
		my $arbc = "";
		if ( defined( $barcode->{$rid}{$ecb} ) && defined( $barcode->{$rid}{$ecb}{random_pos} ) ) {
			$arbc = substr( $reads{$rid}[1], $barcode->{$rid}{$ecb}{random_pos}, $barcode->{$rid}{$ecb}{random_len} );
		}

		# look for linker sequence, if known
		my $any_matched = 0;
		if( defined( $barcode->{$rid}{$ecb} ) && defined( $barcode->{$rid}{$ecb}{linker} ) ) {
			eval {
			require String::Approx;
			my $linker_start = $barcode->{$rid}{$ecb}{trim_len} - length($barcode->{$rid}{$ecb}{linker}) - 5;
			$any_matched = String::Approx::amatch($barcode->{$rid}{$ecb}{linker}, [ "I2","D2","S25%" ], substr( $reads{$rid}[1], $linker_start, length($barcode->{$rid}{$ecb}{linker}) + 10 ));
#			print STDERR "Found linker ".$barcode->{$rid}{$ecb}{linker}." in ".substr( $reads{$rid}[1], $linker_start, length($barcode->{$rid}{$ecb}{linker}) + 10 )."\n" if $any_matched;
			$linker = $any_matched && $linker;
			$defined_bc = 0 unless $any_matched;
			$arbc = "" unless $any_matched;
			$r1d = 0 if !$any_matched && $rid == 1;
			$r2d = 0 if !$any_matched && $rid == 4;
		}
		}
		$rbc[$rid] = $arbc;

		# trim the read if the linker was found
		if ( $any_matched && defined( $barcode->{$rid}{$ecb} ) && defined( $barcode->{$rid}{$ecb}{trim_len} ) ) {
			$reads{$rid}[1] = substr( $reads{$rid}[1], $barcode->{$rid}{$ecb}{trim_len} );
			$reads{$rid}[3] = substr( $reads{$rid}[3], $barcode->{$rid}{$ecb}{trim_len} );
		}
	}

	# have the barcodes. tag the sequence with the barcodes and sample names
	my $joined_bc  = join( ":", @bc );
#	print "Joined bc $joined_bc\n";
	my $joined_rbc = join( ":", @rbc );
	chomp $reads{1}[0];
	$reads{1}[0] .= " barcode=$joined_bc";
	$reads{1}[0] .= " random=$joined_rbc" if $defined_bc && @rbc > 0;
	chomp $reads{4}[0];
	$reads{4}[0] .= " barcode=$joined_bc";
	$reads{4}[0] .= " random=$joined_rbc" if $defined_bc && @rbc > 0;
	if ( $defined_bc && defined( $sample_map->{$joined_bc} && $linker) ) {
		$sample_counts->{ $sample_map->{$joined_bc} } = 0 unless defined $sample_counts->{ $sample_map->{$joined_bc} };
		$sample_counts->{ $sample_map->{$joined_bc} }++;
		my $readid = $sample_map->{$joined_bc}."_".$sample_counts->{ $sample_map->{$joined_bc} };
		my $readid1 = $readid.($include_numeric_read_id ? "/1 " : " ");
		my $readid2 = $readid.($include_numeric_read_id ? "/2 " : " ");
		$reads{1}[0] =~ s/^\@/\@$readid1/;
		$reads{4}[0] =~ s/^\@/\@$readid2/;
		$sample_bc = 1;
	}
	$reads{1}[0] .= "\n";
	$reads{4}[0] .= "\n";

	return ( reads => \%reads, bc => \@bc, rbc => \@rbc, defined => $defined_bc, sample => $sample_bc, r1d => $r1d, r2d => $r2d );
}

# file format
# barcode	bc_name	read_id trim_len random_pos random_len
#
# barcode: nucleotide sequence of the barcode
# bc_name: a name for the barcode
# read_id: the read containing the barcode -- For Illumina this will usually be 1, 2, 3, or 4
# trim_len: if not empty, then the amount to trim from the read containing the barcode
# random_pos: if not empty, position in the barcode read of a randomly synthesized read identifier
# random_len: if not empty, length in the read of a randomly synthesized read identifier
sub read_barcodes {
	my %args = @_;
	my $barcode_file = $args{barcode_file} || miss("barcode_file");
	my %barcode;
	my $INBC = ps_open( $args{barcode_file} );
	while ( my $line = <$INBC> ) {
		chomp($line);
		next if $line =~ /^#/;    # skip a header line
		my @line = split( /\t/, $line );
		next if @line < 2; # probably an empty line, skip it.
		my $rid = $line[2];
		if(defined($args{rev_barcode}) && $rid == $args{rev_barcode}){
			# reverse complement the DNA sequence
			$line[0] = reverse($line[0]);
		    $line[0] =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;			
		}
		warn("Read ID of $rid may be out of range for barcode $line[0]") unless $rid > 0 && $rid <= 4;
		$barcode{$rid}{ $line[0] }{name}       = $line[1];
		$barcode{$rid}{ $line[0] }{parent}     = $line[0];
		$barcode{$rid}{ $line[0] }{trim_len}   = $line[3] if defined $line[3];
		$barcode{$rid}{ $line[0] }{random_pos} = $line[4] if defined $line[4];
		$barcode{$rid}{ $line[0] }{random_len} = $line[5] if defined $line[5];
		$barcode{$rid}{ $line[0] }{linker}     = $line[6] if defined $line[6];
		
		die "Error, can not handle multiple random barcode lengths\n" if defined($RANDOM_LEN) && defined($line[5]) && $RANDOM_LEN != $line[5];
		$RANDOM_LEN = $barcode{$rid}{ $line[0] }{random_len} if defined($line[5]) && !defined($RANDOM_LEN);

		# insert all single-error barcodes
		for ( my $i = 0; $i < length( $line[0] ); $i++ ) {
			my @chars = ( "A", "C", "G", "T", "N" );
			my $s = $line[0];
			foreach my $ck (@chars) {
				substr( $s, $_, 1 ) =~ s/[ACGT]/$ck/ for $i;
				print STDERR "Barcode collision! $s => $line[0] was already mapped to $barcode{$rid}{$s}{parent} for read $line[2]!!\n" if defined $barcode{$rid}{$s} && $s ne $line[0];
				$barcode{$rid}{$s}{parent} = $line[0];
			}
		}
	}
	die "Error, no barcodes found!\nPlease check the barcode file for proper tab-delimited formatting\n" if(keys(%barcode)==0);
	close($INBC);
	return %barcode;
}

# read barcodes to names mapping in QIIME mapping format
# column 1: sample name
# column 2: barcode, in XXXX:YYYY format if dual barcode
# columns 3+: arbitrary metadata
sub read_sample_file {
	my %args = @_;
	my $sample_file = $args{sample_file} || miss("sample_file");

	my %barcode_map;
	debug "Reading name mapping\n";
	my $INMAP = ps_open($sample_file);
	while ( my $line = <$INMAP> ) {
		chomp($line);
		next if $line =~ /^#/;
		my @line = split( /\t/, $line );
		next if @line < 2;
		# handle tab-delimited
		$line[1] .= ":$line[2]" if($line[1] !~ /\:/ && defined($line[2]));
		if(defined($args{swap_barcodes})){
			my @bc = split(/:/, $line[1]);
			$line[1] = $bc[1].":".$bc[0];
		}
		if(defined($args{rev_barcode})){
			my @bc = split(/:/, $line[1]);
			if($args{rev_barcode}==3){
				$bc[0] = reverse($bc[0]);
			        $bc[0] =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;			
			}
			if($args{rev_barcode}==2){
				$bc[1] = reverse($bc[1]);
			        $bc[1] =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;			
			}
			$line[1] = $bc[0].":".$bc[1];
		}

		print STDERR "Mapping barcode $line[1] to $line[0]\n";
		$barcode_map{ $line[1] } = $line[0];
	}
	close($INMAP);
	return %barcode_map;
}

# clean up a fastq header line so it's compatible with MG-RAST
sub clean_line {
	my %args = @_;
	$args{line} =~ s/ /:/g;
	chomp $args{line};
	$args{line} .= "/".$args{num}."\n";
	return $args{line};
}

1;