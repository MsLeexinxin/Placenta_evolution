#!usr/bin/perl;
my %seq;
my %species;
my ($aa_fa,$pep2sym)=@ARGV;
open IN,"$aa_fa";
while (<IN>){
chomp;
if (/^>/){
my @list=split/\s+/,$_;
$id=$list[0];
}else{
$seq{$id}.=$_;
}
}
open IN,"$pep2sym";
while (<IN>){
chomp;
my @list=split/\t/,$_;
my $key=">$list[0]";
$species{$key}=$list[1];
}
close IN;

my %hash;
foreach my $key (keys %seq){
if ((exists $seq{$key}) && (not exists $hash{$species{$key}})){
my @list=split/>/,$key;
print ">$list[1] gene:$list[1] gene_symbol:$species{$key}\n$seq{$key}\n";
$hash{$species{$key}}=$key;
}
}


