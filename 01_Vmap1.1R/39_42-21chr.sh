#42-21.sh
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
		
		
		
#21-42.sh
cat $1 | awk 'BEGIN{FS="\t";OFS="\t"}{
         if($0~/^chromo/){
			print $0;
			}
		else{
			if($1~/chr1A/ && $4 <= 471304005){
				$1="1";
				print $0;
				}
			else if($1~/chr1A/ && $4 > 471304005){
				$1="2";
				$4=$4-471304005;
				print $0;
				}
			else if($1~/chr1B/ && $4 <= 438720154){$1="3";print $0;}
			else if($1~/chr1B/ && $4 > 438720154){$1="4";$4=$4-438720154;print $0}
			else if($1~/chr1D/ && $4 <= 452179604){$1="5";print $0}
			else if($1~/chr1D/ && $4 > 452179604){$1="6";$4=$4-452179604;print $0}
			else if($1~/chr2A/ && $4 <= 462376173){$1="7";print $0}
			else if($1~/chr2A/ && $4 > 462376173){$1="8";$4=$4-462376173;print $0}
			else if($1~/chr2B/ && $4 <= 453218924){$1="9";print $0}
			else if($1~/chr2B/ && $4 > 453218924){$1="10";$4=$4-453218924;print $0}
			else if($1~/chr2D/ && $4 <= 462216879){$1="11";print $0}
			else if($1~/chr2D/ && $4 > 462216879){$1="12";$4=$4-462216879;print $0}
			else if($1~/chr3A/ && $4 <= 454103970){$1="13";print $0}
			else if($1~/chr3A/ && $4 > 454103970){$1="14";$4=$4-454103970;print $0}
			else if($1~/chr3B/ && $4 <= 448155269){$1="15";print $0}
			else if($1~/chr3B/ && $4 > 448155269){$1="16";$4=$4-448155269;print $0}
			else if($1~/chr3D/ && $4 <= 476235359){$1="17";print $0}
			else if($1~/chr3D/ && $4 > 476235359){$1="18";$4=$4-476235359;print $0}
			else if($1~/chr4A/ && $4 <= 452555092){$1="19";print $0}
			else if($1~/chr4A/ && $4 > 452555092){$1="20";$4=$4-452555092;print $0}
			else if($1~/chr4B/ && $4 <= 451014251){$1="21";print $0}
			else if($1~/chr4B/ && $4 > 451014251){$1="22";$4=$4-451014251;print $0}
			else if($1~/chr4D/ && $4 <= 451004620){$1="23";print $0}
			else if($1~/chr4D/ && $4 > 451004620){$1="24";$4=$4-451004620;print $0}
			else if($1~/chr5A/ && $4 <= 453230519){$1="25";print $0}
			else if($1~/chr5A/ && $4 > 453230519){$1="26";$4=$4-453230519;print $0}
			else if($1~/chr5B/ && $4 <= 451372872){$1="27";print $0}
			else if($1~/chr5B/ && $4 > 451372872){$1="28";$4=$4-451372872;print $0}
			else if($1~/chr5D/ && $4 <= 451901030){$1="29";print $0}
			else if($1~/chr5D/ && $4 > 451901030){$1="30";$4=$4-451901030;print $0}
			else if($1~/chr6A/ && $4 <= 452440856){$1="31";print $0}
			else if($1~/chr6A/ && $4 > 452440856){$1="32";$4=$4-452440856;print $0}
			else if($1~/chr6B/ && $4 <= 452077197){$1="33";print $0}
			else if($1~/chr6B/ && $4 > 452077197){$1="34";$4=$4-452077197;print $0}
			else if($1~/chr6D/ && $4 <= 450509124){$1="35";print $0}
			else if($1~/chr6D/ && $4 > 450509124){$1="36";$4=$4-450509124;print $0}
			else if($1~/chr7A/ && $4 <= 450046986){$1="37";print $0}
			else if($1~/chr7A/ && $4 > 450046986){$1="38";$4=$4-450046986;print $0}
			else if($1~/chr7B/ && $4 <= 453822637){$1="39";print $0}
			else if($1~/chr7B/ && $4 > 453822637){$1="40";$4=$4-453822637;print $0}
			else if($1~/chr7D/ && $4 <= 453812268){$1="41";print $0}
			else if($1~/chr7D/ && $4 > 453812268){$1="42";$4=$4-453812268;print $0}
			else{}
			}
		}'
