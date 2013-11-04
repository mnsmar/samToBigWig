#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Path qw(make_path);
use File::Temp;

use lib '/home/mns/lib/perl/class/v7.0';

use MyBio::NGS::Sample;

############################
# Read options and arguments
my $help;
my $memory = 500000000;
GetOptions(
        'h'            => \$help,
        'm=i'          => \$memory,
) or usage();
usage() if $help;

my $sample_folder = shift;
my $ifile_relative_path = shift;
my $ofile_relative_path = shift;

# Sanity check
if (!defined $sample_folder or !defined $ifile_relative_path or !defined $ofile_relative_path) {
	usage();
}


############################
# Initialize
my $time = time;
my $params_file = File::Spec->catfile($sample_folder, 'params.xml');
my $sam_input_file = File::Spec->catfile($sample_folder, $ifile_relative_path);
my $bigwig_output_file = File::Spec->catfile($sample_folder, $ofile_relative_path);


############################
# Extra sanity checks
if ($sam_input_file eq $bigwig_output_file) {
	die "Input file and output file should not be the same.\nInput: $sam_input_file\nOutput: $bigwig_output_file\n";
}
if (!-e $sam_input_file) {
	die "Input file does not exists. $sam_input_file\n";
}
if (!-e $params_file) {
	die "File with params does not exists. $params_file\n";
}


############################
# Create output directory if needed
my (undef, $output_directory) = File::Spec->splitpath($bigwig_output_file); make_path($output_directory);


############################
# Create a temporary directory for intermediate files
my $temp_dir = File::Temp->newdir();
warn "Temporary dir created: ".$temp_dir->dirname."\n";
my $temp_bam = File::Spec->catfile($temp_dir->dirname, 'file.bam');
my $temp_sorted_bam_no_suffix = File::Spec->catfile($temp_dir->dirname, 'file.sorted');
my $temp_bedgraph_plus = File::Spec->catfile($temp_dir->dirname, 'file.sorted.plus.bedgraph');
my $temp_bedgraph_minus = File::Spec->catfile($temp_dir->dirname, 'file.sorted.minus.bedgraph');


############################
# Read sample info
my $sample = MyBio::NGS::Sample->new({XML => $params_file});
my $chr_size_file = '/store/data/UCSC/'.$sample->align->{'bwa'}->{'assembly'}.'/chrom.sizes';


############################
# Prepare command
my $cmd = "samtools view -Sb -o $temp_bam $sam_input_file;".
          "samtools sort -m $memory $temp_bam $temp_sorted_bam_no_suffix;".
          "genomeCoverageBed -split -bg -strand \"+\" -ibam $temp_sorted_bam_no_suffix.bam -g $chr_size_file > $temp_bedgraph_plus;".
          "genomeCoverageBed -split -bg -strand \"-\" -ibam $temp_sorted_bam_no_suffix.bam -g $chr_size_file > $temp_bedgraph_minus;".
          "bedGraphToBigWig $temp_bedgraph_plus $chr_size_file $bigwig_output_file.plus.bw;".
          "bedGraphToBigWig $temp_bedgraph_minus $chr_size_file $bigwig_output_file.minus.bw";


############################
# Run job
warn "Running command: ".$cmd."\n";
system "$cmd";

warn "Time elapsed:\t".(time-$time)."sec\n";

###########################################
# Subroutines used
###########################################
sub usage {
	print "\nUsage:   $0 [options] sample_folder input_file_relative_path output_file_relative_path\n\n".
	      "Options:\n".
	      "        -t             number of threads to be used\n".
	      "        -devel         enables development mode\n".
	      "        -h             print this help\n\n";
	exit;
}
