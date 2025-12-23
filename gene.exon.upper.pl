#!usr/bin/perl -w
open List1,"$ARGV[0]" || die"$!";
open List2,"$ARGV[1]" || die"$!";
open List3,"$ARGV[2]" || die"$!";
my %hash;

while (<List1>) {
        chomp;
        my @tt=split;
        $chrom=$tt[0];
		$star=$tt[1];
		$end=$tt[2];
		for($i=$star;$i<=$end;$i++){
		$tag=$chrom.":".$i;
		$hash{$tag}=$i;
		}
}

while (<List2>) {
        chomp;
		if(/^>/){
		s/>//;
        @aa=split;
		$name=$aa[0];}
		
	   else{$seq{$name}.=$_;}

}

while (<List3>) {
        chomp;
        @pp=split;
		$string=$seq{$pp[0]};
		print ">$pp[0]\n";
    for($i=1;$i<=$pp[1];$i++){
	$flag=$pp[0].":".$i;
	$num=$i-1;
	$s=substr($string,$num,1);
	if(exists $hash{$flag}){
	$t=uc($s);
	print "$t";}
	else{print "$s";}
	}
	   print "\n";
}
