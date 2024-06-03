#!/usr/bin/perl

while ($line=<>){
    chop($line);
    $line =~ s/^\s+//;
    $line =~ s/D-/E-/g;
#    print STDERR $line,"\n";
    @fields=split /\s+/,$line;
#    print STDERR "Fields ";
#    for ($i=0;$i<=$#fields; $i++) {print STDERR "|",$fields[$i];}
#    print STDERR "\n";
#    if ($line=~/methd iprec istopc   :/){$prec=$fields[4];$ovr=$fields[5];}
    if ($line=~/Generating Matrix/){
	$size=$fields[2];
	$size =~ s/\(size=//;
	$size =~ s/\)//;
    } 
    if ($line=~/Grid dimensions/){
	$grid=$fields[3];
	$grid=~s/x//gi;
    } 
    if ($line=~/Unit cube discretized with dim/){$idim=$fields[5];}
    if ($line=~/Overall matrix creation time :/){$matbld=$fields[5];}
    if ($line=~/Number of levels   :/){$nlev=$fields[4];}
    if ($line=~/variant:/){$pvariant=$fields[1];}
    if ($line=~/Degree:/){$pdegree=$fields[1];}
    if ($line=~/Smoother:/){
	$line=<>;
	chop($line);
	$line =~ s/^\s+//;
	$line =~ s/D-/E-/g;
	$smoother=$line;
	if ($smoother =~/Point Jacobi/) {
	    $pvariant ="L1-Jacobi";
	}
    }
    if ($line=~/Operator complexity:/){$opcmpl=$fields[2];}
    if ($line=~/Average coarsening :/){$avgcrs=$fields[3];}
    if ($line=~/Coarse Matrix: Global size:/){$crssz=$fields[4];}
    if ($line=~/Linear system size                 :/){$size=$fields[4];} 
    if ($line=~/Time to build hierarchy            :/){$thier=$fields[5];}
    if ($line=~/Time to build smoothers            :/){$tsmth=$fields[5];}
    if ($line=~/Total nonzeros          for A:/){$anz=$fields[4];}
    if ($line=~/Total nonzeros          for PREC:/){$pnz=$fields[4];}
    if ($line=~/Total memory occupation for A/){$asizeof=$fields[6];}
    if ($line=~/Total memory occupation            :/){$szpp=$fields[5];}
    if ($line=~/Time to solve system               :/){$tslv=$fields[5];}
    if ($line=~/Time per iteration                 :/){$titer=$fields[4];}
    if ($line=~/Number of RHS  /){$nrhs=$fields[3];}
    if ($line=~/Iterations to convergence          :/){$nit=$fields[4];}
    if ($line=~/Storage format for A               :/){ $afmt=$fields[5];}
    if ($line=~/Storage type for AGPU:/){ $agmt=$fields[4];}
    if ($line=~/MFLOPS                       \(GPU.\)/){
	$gflps=$fields[3]/1000.0;
    }
    
    if ($line=~/Total memory occupation for DESC_A/){
	printf("%-24s %-18s  %5d   %10.7f \n",
	       $smoother,$pvariant,$pdegree,$titer);
	
    }

    
}

