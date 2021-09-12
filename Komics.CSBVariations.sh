#!/bin/bash


####################################################
#
#
# CHECK FOR CSB3 VARIATIONS
# =========================
#
# Description
#
# 	KOMICS retrieves contigs based on 
#		GGGGTTG[G/A]TGTA
# 	Check whether other variations exists
#
#	a) grep CSB3 variation in CSB1 containing contigs
#	b) check whether CSB3 variation is close to CSB1
#	c) check whether contig is of expected length
#		(700-1300 bp)
#	d) count occurence of CSB3 variation
#
# Usage
#
#
#
####################################################


#####	A MAKE LIST OF CSB3 VARIATIONS

CSB1='GGGCGT[T|G]C'
CSB1compl='G[C|A]ACGCCC'
CSB3=GGGGTTGGTGTA
COUNTER=0

CSB3var=$(for i in {1..12}; do COUNTER=$(( COUNTER + 1 )); echo $CSB3 | sed "s/././$COUNTER"; done | paste -s -d '|')

CSB3='GGGGTTG[G|A]TGTA'
CSB3compl='TACA[C|T]CAACCCC'


#####	B FUNCTIONS

function fasta_oneline {

  tr ' ' '_' | awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' | tr ' ' '\n'

}


##### 	C CHECK FOR CSB3 VARIATIONS

echo -e '\n\n\n---------------------- CSB3 variations ----------------------\n\n\n'

for file in kDNA/minicircles/*contigs.fasta; do

	# check if CSB3 variations other than GGGGTTG[G/A]TGTA
	# exists in contigs having CSB1 
        # var=($(cat $file | fasta_oneline | grep -E -B1 "$CSB1|$CSB1compl" | grep -Eo -B1 "$CSB3var|$(./complement.sh $CSB3var)" | sort | uniq | grep -vE "$CSB3|$CSB3compl"))
	var=($(cat $file | fasta_oneline | grep -E "$CSB1" | grep -Eo "$CSB3var" | sort | uniq | grep -vE "$CSB3"))

	# if CSB3 variations exists
        if [ ${#var[@]} -gt 0 ]; then

		# print file name
                echo -e "In file" $file "found following variations on CSB3: \n"

        fi

	# for every CSB3 variations
	for i in "${var[@]}"; do

		# print whether max 100 nucleotides from CSB1
       	        closetocsb1=$(cat $file | fasta_oneline | grep -o '.\{0,0\}GGGCGT.\{0,100\}' | grep $i)
       	        if [ -z "$closetocsb1" ]; then
			echo -e "\t" $i "\t --> not close to CSB1" "\n"
		else
			echo -e "\t" $i "\t --> close to CSB1"
			length=$(cat $file | fasta_oneline | grep $i | wc -c)
			contig=$(cat $file | fasta_oneline | grep -B1 $i | head -n1)

			# print whether contig is of expected length
			if [ $length -gt 1300 ]; then
				echo -e "\t\t\t" "-->" $contig "is to big:" $length "\n"
			elif [ $length -lt 700 ]; then
				echo -e "\t\t\t" "-->" $contig "is to small:" $length "\n"
			else
				# if yes, print occurenca
				echo -e "\t\t\t" "-->" $contig "is of expected length:" $length
				occurence=$(cat $file | grep -c $i)
				echo -e "\t\t\t" "-->" "occurence of CSB3 variation:" $occurence "\n"
			fi

		fi

		# closetocsb1compl=$(cat $file | fasta_oneline | grep -o '.\{0,100\}ACGCCC.\{0,0\}' | grep $i)
		# if [ ! -z $closetocsb1compl ]; then
		#	echo -e "\t" $i "\t --> check compl"
		#	echo $closetocsb1compl
		# fi

	done

done


#####   D MAKE LIST OF CSB1 VARIATIONS

CSB1='GGGCGTGC'
COUNTER=0

CSB1var=$(for i in {1..8}; do COUNTER=$(( COUNTER + 1 )); echo $CSB1 | sed "s/././$COUNTER"; done | paste -s -d '|')

CSB1='GGGCGT[T|G]C'
CSB1compl='G[C|A]ACGCCC'


#####   D CHECK FOR CSB3 VARIATIONS

echo -e '\n\n---------------------- CSB1 variations ----------------------\n\n\n'

for file in kDNA/minicircles/*contigs.fasta; do

        # check if CSB1 variations other than GGGCGT[T|G]C
        # exists in contigs having CSB3
        var=($(cat $file | fasta_oneline | grep -E "$CSB3" | grep -Eo "$CSB1var" | sort | uniq | grep -vE "$CSB1"))

        # if CSB1 variations exists
        if [ ${#var[@]} -gt 0 ]; then

                # print file name
                echo -e "In file" $file "found following variations on CSB1: \n"

        fi

        # for every CSB1 variations
        for i in "${var[@]}"; do

                # print whether max 100 nucleotides from CSB3
                closetocsb3=$(cat $file | fasta_oneline | grep -o '.\{0,100\}GGGGTTG.\{0,0\}' | grep $i)
                if [ -z "$closetocsb3" ]; then
                        echo -e "\t" $i "\t --> not close to CSB3" "\n"
                else
                    	echo -e "\t" $i "\t --> close to CSB3"
                        length=$(cat $file | fasta_oneline | grep $i | wc -c)
                        contig=$(cat $file | fasta_oneline | grep -B1 $i | head -n1)

                        # print whether contig is of expected length
                        if [ $length -gt 1300 ]; then
                                echo -e "\t\t\t" "-->" $contig "is to big:" $length "\n"
                        elif [ $length -lt 700 ]; then
                                echo -e "\t\t\t" "-->" $contig "is to small:" $length "\n"
                        else
                                # if yes, print occurenca
                                echo -e "\t\t\t" "-->" $contig "is of expected length:" $length
                                occurence=$(cat $file | grep -c $i)
                                echo -e "\t\t\t" "-->" "occurence of CSB3 variation:" $occurence "\n"
                        fi

                fi

        done

done


#####	D EXTRA

# grep '>' kDNA/minicircles/tmp.108AT.othercontigs.fasta | awk -F '_' '$2<1200 && $2>800 {print $0}' | while read line; do
# 		awk -v line="$line" '/^>/ { p = ($0 ~ line)} p' kDNA/minicircles/tmp.108AT.csb3contigs.fasta | \
# done




