#!/usr/bin/env perl
use Modern::Perl;
use autodie;


##############################################
# Import external libraries
use Getopt::Long;
use File::Path qw(make_path);
use File::Temp;


##############################################
# Read command options
my ($help, $chr_size_file, $ifile, $ofile_prefix);
my $memory = 500000000;
GetOptions(
        'h'               => \$help,
        'chr_size_file=s' => \$chr_size_file,
        'i=s'             => \$ifile,
        'o_prefix=s'      => \$ofile_prefix,
        'm=i'             => \$memory,
) or usage();
usage() if $help;

# Sanity check
map {defined $_ or usage()} ($chr_size_file, $ifile, $ofile_prefix);


############################
# Initialize
my $time = time;


############################
# Extra sanity checks
die "Input and output can not be the same file\n" if ($ifile eq $ofile_prefix);
die "Input file does not exists. $ifile\n" if (!-e $ifile);


############################
# Create output directory if needed
my (undef, $output_directory) = File::Spec->splitpath($ofile_prefix); make_path($output_directory);


############################
# Create a temporary directory for intermediate files
my $temp_dir = File::Temp->newdir();
warn "Temporary dir created: ".$temp_dir->dirname."\n";
my $temp_bam = File::Spec->catfile($temp_dir->dirname, 'file.bam');
my $temp_sorted_bam_no_suffix = File::Spec->catfile($temp_dir->dirname, 'file.sorted');
my $temp_bedgraph_plus = File::Spec->catfile($temp_dir->dirname, 'file.sorted.plus.bedgraph');
my $temp_bedgraph_minus = File::Spec->catfile($temp_dir->dirname, 'file.sorted.minus.bedgraph');


############################
# Prepare command
my $cmd = "samtools view -Sb -o $temp_bam $ifile;".
          "samtools sort -m $memory $temp_bam $temp_sorted_bam_no_suffix;".
          "genomeCoverageBed -split -bg -strand \"+\" -ibam $temp_sorted_bam_no_suffix.bam -g $chr_size_file > $temp_bedgraph_plus;".
          "genomeCoverageBed -split -bg -strand \"-\" -ibam $temp_sorted_bam_no_suffix.bam -g $chr_size_file > $temp_bedgraph_minus;".
          "bedGraphToBigWig $temp_bedgraph_plus $chr_size_file $ofile_prefix.plus.bw;".
          "bedGraphToBigWig $temp_bedgraph_minus $chr_size_file $ofile_prefix.minus.bw";


############################
# Run job
warn "Running command: ".$cmd."\n";
system "$cmd";

warn "Time elapsed:\t".(time-$time)."sec\n";

###########################################
# Subroutines used
###########################################
sub usage {
	print "\nUsage:   $0 [options]\n\n".
	      "Convert a SAM file into BigWig. Two BigWig files, one for each strand, are created.\n".
	      "Options:\n".
	      "        -chr_size_file <Str>   file that has the chromosome sizes\n".
	      "        -i <Str>               input SAM file\n".
	      "        -o_prefix <Str>        the output SAM filepaths will be created by appending the sufixes .plus.bw and .minus.bw to this value\n".
	      "        -m <Int>               the maximum number of memory in bytes to be used during the sorting step (Default: 500000000)\n".
	      "        -h                     print this help\n\n";
	exit;
}
