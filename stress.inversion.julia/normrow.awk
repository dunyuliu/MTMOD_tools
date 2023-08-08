{
    sum=0.0;
    if(($1 != "") && (substr($1,1,1) != "#")){
	for(i=1;i <= NF;i++)
	    sum += ($i)*($i);
	print(sqrt(sum));
    }
}
