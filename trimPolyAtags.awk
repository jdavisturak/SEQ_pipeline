{
if(NR%4==1){a=$0} if(NR%4==2) {b=$0} if(NR%4==3) {c=$0}   

if(NR%4==0){
	d=$0
  if (b ~ /^TTTTT/){ gsub(/^T*/,"",b);}
  if (b ~ /AAAAA$/){ gsub(/A*$/,"",b)}  
  d=substr(d,length(d)-length(b)+1,length(b));
  
  printf("%s\n%s\n%s\n%s\n",a,b,c,d)
} 
}

