#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 [options] genome.fasta > genome.RD.fasta

    --length <int>    default: 100000
    将长度低于此设定值的基因组序列和全基因组序列进行BLAST比对。

    --coverage <float>    default: 0.90
    --identity <float>    default: 0.95
    --evalue <float>    default: 1e-10
    若序列和长度大于该序列的基因组其它序列进行比对时，其结果的evalue、identity和coverage值优于以上阈值，则认为该序列属于基因组中的重复序列，对其予以剔除。

    --CPU <int>    default: 8
    设置blastn命令运行时使用CPU线程数。

    --tmp <string>    default: tmp
    设置临时文件夹路径。

USAGE
if (@ARGV==0){die $usage}

my ($length, $coverage, $identity, $evalue, $CPU, $tmp);
GetOptions(
    "length:s" => \$length,
    "coverage:f" => \$coverage,
    "identity:f" => \$identity,
    "evalue:f" => \$evalue,
    "CPU:i" => \$CPU,
    "tmp:s" => \$tmp,
);
$length ||= 100000;
$coverage ||= 0.90;
$identity ||= 0.95;
$evalue ||= 1e-10;
$CPU ||= 8;
$tmp ||= "tmp";

mkdir $tmp unless -e $tmp;
my $pwd = `pwd`;
chomp($pwd);
my $input = $ARGV[0];
unless ($input =~ m/^\//) {
    $input = "$pwd/$input";
}
unless ($tmp =~ m/^\//) {
    $tmp = "$pwd/$tmp";
}
chdir $tmp;

my $cmdString = "makeblastdb -in $input -dbtype nucl -title genome -parse_seqids -out genome -logfile makeblastdb.log";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

open IN,$input or die "Can not open file $input, $!\n";
my (%seq, $seq_id, %length, @seq);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) { $seq_id = $1; push @seq, $1; }
    else { $seq{$seq_id} .= $_; $length{$seq_id} += length($_); }
}
close IN;

open OUT, ">", "short_sequences.fasta" or die "Can not create file $tmp/short_sequences.fasta, $!\n";
my @short_seq_id;
foreach (sort {$length{$b} <=> $length{$a}} keys %seq) {
    if ($length{$_} < $length) {
        print OUT ">$_\n$seq{$_}\n";
        push @short_seq_id, $_;
    }
}
close OUT;

my $perc_identity = $identity * 100;
$cmdString = "blastn -query $tmp/short_sequences.fasta -db $tmp/genome -out $tmp/blast.out -evalue 1e-10 -perc_identity $perc_identity -outfmt 7 -num_threads $CPU";
print STDERR "CMD: $cmdString\n";
(system $cmdString) == 0 or die "Failed to execute: $cmdString\n";

open IN, "$tmp/blast.out" or die "Can not open file $tmp/blast.out, $!\n";
my (%out, %align);
while (<IN>) {
    next if m/^#/;
    @_ = split /\t/;
    next if $_[0] eq $_[1];
    next if $length{$_[0]} > $length{$_[1]};
    my $align_length = $_[7] - $_[6] + 1;
    my $cov = $align_length / $length{$_[0]};
    if ($cov >= 0.01 or $align_length >= 5000) {
        $out{$_[0]} .= "$_[0]\t$_[1]\t$_[6]\t$_[7]\t$_[8]\t$_[9]\t$cov\n";
        $align{$_[0]}{"$_[6]\t$_[7]}"} = 1;
    }
}
close IN;

my %cov;
foreach my $id (keys %align) {
    my @cov = sort {$a <=> $b} keys %{$align{$id}};
    my $match_length = &match_length(@cov);
    $cov{$id} = $match_length / $length{$id};
}

print STDERR "Coverage of short sequences:\n";
foreach (sort {$length{$b} <=> $length{$a}} @short_seq_id) {
    my $cov = 0;
    $cov = $cov{$_} if $cov{$_};
    print STDERR "$_\t$cov\n";
}

open OUT, ">", "alignment.tab" or die "Can not open file $tmp/alignment.tab, $!\n";
print OUT "query\tsubject\tquery_start\tquery_end\tsubject_start\tsubject_end\tcoverage\n";
print STDERR "Duplicated sequences:\n";
foreach (sort {$length{$a} <=> $length{$b}} @short_seq_id) {
    if ($cov{$_} > $coverage) {
        print OUT $out{$_};
        print STDERR "$_\n";
    }
}
close OUT;
print STDERR "The detail alignment infomation for these duplicated sequences: $tmp/alignment.tab\n";

foreach (@seq) {
    print ">$_\n$seq{$_}\n" unless $cov{$_} > $coverage;
}

sub match_length {
    my @inter_sorted_site;
    foreach (@_) {
        my @aaa = $_ =~ m/(\d+)/g;
        @aaa = sort { $a <=> $b } @aaa;
        push @inter_sorted_site, "$aaa[0]\t$aaa[1]";
    }
    @inter_sorted_site = sort { $a <=> $b } @inter_sorted_site;

    my $out_site_number;
    my $former_region = shift @inter_sorted_site;
    my @aaa = $former_region =~ m/(\d+)/g;
    $out_site_number += ($aaa[1] - $aaa[0] + 1);
    foreach (@inter_sorted_site) {
        my @former_region = $former_region =~ m/(\d+)/g;
        my @present_region = $_ =~ m/(\d+)/g;
        
        if ($present_region[0] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $present_region[0] + 1);
            $former_region = $_;
        }
        elsif ($present_region[1] > $former_region[1]) {
            $out_site_number += ($present_region[1] - $former_region[1]);
            $former_region = $_;
        }
        else {
            next
        }
    }
    return $out_site_number;
}
