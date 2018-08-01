#!/usr/bin/perl

if(!$cutoff){
$cutoff = 30;
}

while(<>){
chomp($_);
	if(/>.*/){
	print $_,"\n";
	}
	else{
		while($_ =~ /N+/g){
			if(length($&) < $cutoff){
			$_ = $`.$';
			}
		}
	
	$_ =~ s/^N+//;
	$_ =~ s/N+$//;	
	print $_,"\n";
	}
}
