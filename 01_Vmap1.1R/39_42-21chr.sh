cat $1 | awk 'BEGIN{FS="\t";OFS="\t"}{
         if($0~/^Chr/){
			print $0;
			}
		else{
			if($1==1){
				$1="chr1A";
				print $0;
				}
			else if($1==2){
				$1="chr1A";
				$2=$2+471304005;
				print $0;
				}
			else if($1==3 ){$1="chr1B";print $0;}
			else if($1==4 ){$1="chr1B";$2=$2+438720154;print $0}
			else if($1==5 ){$1="chr1D";print $0}
			else if($1==6 ){$1="chr1D";$2=$2+452179604;print $0}
			else if($1==7 ){$1="chr2A";print $0}
			else if($1==8 ){$1="chr2A";$2=$2+462376173;print $0}
			else if($1==9 ){$1="chr2B";print $0}
			else if($1==10 ){$1="chr2B";$2=$2+453218924;print $0}
			else if($1==11 ){$1="chr2D";print $0}
			else if($1==12 ){$1="chr2D";$2=$2+462216879;print $0}
			else if($1==13 ){$1="chr3A";print $0}
			else if($1==14 ){$1="chr3A";$2=$2+454103970;print $0}
			else if($1==15 ){$1="chr3B";print $0}
			else if($1==16 ){$1="chr3B";$2=$2+448155269;print $0}
			else if($1==17 ){$1="chr3D";print $0}
			else if($1==18 ){$1="chr3D";$2=$2+476235359;print $0}
			else if($1==19 ){$1="chr4A";print $0}
			else if($1==20 ){$1="chr4A";$2=$2+452555092;print $0}
			else if($1==21 ){$1="chr4B";print $0}
			else if($1==22 ){$1="chr4B";$2=$2+451014251;print $0}
			else if($1==23 ){$1="chr4D";print $0}
			else if($1==24 ){$1="chr4D";$2=$2+451004620;print $0}
			else if($1==25 ){$1="chr5A";print $0}
			else if($1==26 ){$1="chr5A";$2=$2+453230519;print $0}
			else if($1==27 ){$1="chr5B";print $0}
			else if($1==28 ){$1="chr5B";$2=$2+451372872;print $0}
			else if($1==29 ){$1="chr5D";print $0}
			else if($1==30 ){$1="chr5D";$2=$2+451901030;print $0}
			else if($1==31 ){$1="chr6A";print $0}
			else if($1==32 ){$1="chr6A";$2=$2+452440856;print $0}
			else if($1==33 ){$1="chr6B";print $0}
			else if($1==34 ){$1="chr6B";$2=$2+452077197;print $0}
			else if($1==35 ){$1="chr6D";print $0}
			else if($1==36 ){$1="chr6D";$2=$2+450509124;print $0}
			else if($1==37 ){$1="chr7A";print $0}
			else if($1==38 ){$1="chr7A";$2=$2+450046986;print $0}
			else if($1==39 ){$1="chr7B";print $0}
			else if($1==40 ){$1="chr7B";$2=$2+453822637;print $0}
			else if($1==41 ){$1="chr7D";print $0}
			else if($1==42 ){$1="chr7D";$2=$2+453812268;print $0}
			else{}
			}
		}'