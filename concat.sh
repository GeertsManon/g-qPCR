#!/bin/bash


usage() {
	echo -e "\n\tusage\t\t./concat -t <filetype> -o <filename>"
        echo -e "\t-t\t\tfile type"
        echo -e "\tfiletypes\tlog"
        echo -e "\t\t\tkDNA/minicircles/*.mapping.stats.txt"
        echo -e "\t\t\tkDNA/minicircles/*.depth.stats.txt"
        echo -e "\t\t\tsummary"
        echo -e "\t-o\t\toutput file name\n"
}

while getopts "t:o:" opt; do
    case $opt in
        t) filetype=${OPTARG};;
	o) filename=${OPTARG};;
    esac
done

if [[ $# -ne 4 ]] ; then
	usage
	exit
else


####	A LOG FILES

	if [[ $filetype == "log" ]]; then
		ls *log | while read file; do
			echo $(echo $file | sed 's/.log//g') > tmp.$file
			cat $file | awk -F ':  ' '{print $2}' >> tmp.$file
		done
		paste <(echo 'sample'; cat 108AT.log | awk -F ':' '{print $1}') tmp.* > $filename
		rm tmp.*


####	B SUMMARY VARIANT FILTERS

        elif [[ $filetype == "summary" ]]; then
                ls *.summary | while read file; do
                        echo $(echo $file | sed 's/.summary//g') > tmp.$file
                        cat $file | awk -F ': ' '{print $2}' >> tmp.$file
                done
                paste <(echo 'variant filter'; cat QUAL.summary | awk -F ':' '{print $1}') tmp.* > $filename
                #rm *.summary



#####	C MAPPING STATS KOMICS

	elif [[ $filetype == "kDNA/minicircles/*.mapping.stats.txt" ]]; then
		ls kDNA/minicircles/*.mapping.stats.txt |  while read file; do 
			sample=$(echo $file | sed 's/.mapping.stats.txt//g' | sed 's/kDNA\/minicircles\///g')
			echo $sample > tmp.$sample
			awk -F ': ' '{print $2}' $file >> tmp.$sample
		done
		paste <(echo 'sample'; awk -F ':' '{print $1}' kDNA/minicircles/108AT.mapping.stats.txt) tmp.* > $filename
		rm tmp.* 


#####	C DEPTH STATS KOMICS

#	elif [[ $filetype == "kDNA/minicircles/*.depth.stats.txt" ]]; then
#	sample=$(echo $file | sed 's/.depth.stats.txt//g' | sed 's/kDNA\/minicircles\///g')
#	echo $sample $(awk '{ total += $3; count++ } END { print total/count }' $file)
#done | column -t

	else
		echo -e "\n\tERROR: file type doesn't exist"

	fi
fi
