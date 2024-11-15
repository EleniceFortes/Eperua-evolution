#################
###11 May 2023###
#################

# create directory Raw_data2
	cd Eperua
	mkdir Raw_data2

# copy raw data of plate 2 from Jacob`s server to my server:

	cp /data/jbl256/Eperua/raw_data_plate2/*fq.gz /data/eaf236/Eperua/Raw_data2/
	
# unzip reads in Raw_data2 (128 files in the directory)
	cd Raw_data2
	gunzip *.gz
	
# for specimens repeated in both plates, copy raw reads from Raw_data to Raw_data2 to be able to concatenate everything
	#each of this specimens have 8 reads, 4 foward and 4 reverse (OBTU207 6 F and 6 R), the other have just 2 foward and 2 reverse
	#repeated specimens: CERRA01 CARDE2632 MANA194b FROE1458 PRAE869 OBTU207 
	
	cp /data/eaf236/Eperua/Raw_data/CARD2632* /data/eaf236/Eperua/Raw_data2/
	cp /data/eaf236/Eperua/Raw_data/CERR01* /data/eaf236/Eperua/Raw_data2/
	cp /data/eaf236/Eperua/Raw_data/MANA194B* /data/eaf236/Eperua/Raw_data2/
	cp /data/eaf236/Eperua/Raw_data/FROE1458* /data/eaf236/Eperua/Raw_data2/
	cp /data/eaf236/Eperua/Raw_data/PRAE869* /data/eaf236/Eperua/Raw_data2/
	cp /data/eaf236/Eperua/Raw_data/OBTU207* /data/eaf236/Eperua/Raw_data2/
	
	#DONE 156 files (128 + (5*4) + 8) = 156

# Concatenate fowards and reverse reads

	cd Eperua 
	mkdir concatenated_raw_reads2
	
	cd Raw_data2
	
	#cat REPLACE*1.fq > ../concatenated_raw_reads2/REPLACE_combined.R1.fq
	#cat REPLACE*2.fq > ../concatenated_raw_reads2/REPLACE_combined.R2.fq
	
	#Test run
	
	cat CERR*1.fq > ../concatenated_raw_reads2/CERR01_combined.R1.fq
	cat CERR*2.fq > ../concatenated_raw_reads2/CERR01_combined.R2.fq
	
	cat PR*869*1.fq > ../concatenated_raw_reads2/PRAE869_combined.R1.fq
	cat PR*869*2.fq > ../concatenated_raw_reads2/PRAE869_combined.R2.fq
	
		
	#put all others together in one run
	
	cat AUGO1496*1.fq > ../concatenated_raw_reads2/AUGO1496_combined.R1.fq
	cat AUGO1496*2.fq > ../concatenated_raw_reads2/AUGO1496_combined.R2.fq
	
	cat AUGO2974*1.fq > ../concatenated_raw_reads2/AUGO2974_combined.R1.fq
	cat AUGO2974*2.fq > ../concatenated_raw_reads2/AUGO2974_combined.R2.fq
	
	cat AUGOA4*1.fq > ../concatenated_raw_reads2/AUGOA4_combined.R1.fq
	cat AUGOA4*2.fq > ../concatenated_raw_reads2/AUGOA4_combined.R2.fq
	
	cat CARD*1.fq > ../concatenated_raw_reads2/CARD2632_combined.R1.fq
	cat CARD*2.fq > ../concatenated_raw_reads2/CARD2632_combined.R2.fq
	
	cat COPA1661*1.fq > ../concatenated_raw_reads2/COPA1661_combined.R1.fq
	cat COPA1661*2.fq > ../concatenated_raw_reads2/COPA1661_combined.R2.fq
	
	cat DUCK190b*1.fq > ../concatenated_raw_reads2/DUCK190b_combined.R1.fq
	cat DUCK190b*2.fq > ../concatenated_raw_reads2/DUCK190b_combined.R2.fq
	
	cat FALC165*1.fq > ../concatenated_raw_reads2/FALC165_combined.R1.fq
	cat FALC165*2.fq > ../concatenated_raw_reads2/FALC165_combined.R2.fq
	
	cat FRO32418*1.fq > ../concatenated_raw_reads2/FRO32418_combined.R1.fq
	cat FRO32418*2.fq > ../concatenated_raw_reads2/FRO32418_combined.R2.fq
	
	cat FROE1458*1.fq > ../concatenated_raw_reads2/FROE1458_combined.R1.fq
	cat FROE1458*2.fq > ../concatenated_raw_reads2/FROE1458_combined.R2.fq
	
	cat GLBI192b*1.fq > ../concatenated_raw_reads2/GLBI192b_combined.R1.fq
	cat GLBI192b*2.fq > ../concatenated_raw_reads2/GLBI192b_combined.R2.fq
	
	cat MANA194*1.fq > ../concatenated_raw_reads2/MANA194b_combined.R1.fq
	cat MANA194*2.fq > ../concatenated_raw_reads2/MANA194b_combined.R2.fq
	
	cat OBTU207*1.fq > ../concatenated_raw_reads2/OBTU207_combined.R1.fq
	cat OBTU207*2.fq > ../concatenated_raw_reads2/OBTU207_combined.R2.fq
	
	cat OLE17008*1.fq > ../concatenated_raw_reads2/OLE17008_combined.R1.fq
	cat OLE17008*2.fq > ../concatenated_raw_reads2/OLE17008_combined.R2.fq
	
	cat OLEI1901*1.fq > ../concatenated_raw_reads2/OLEI1901_combined.R1.fq
	cat OLEI1901*2.fq > ../concatenated_raw_reads2/OLEI1901_combined.R2.fq
	
	cat PRAE4758*1.fq > ../concatenated_raw_reads2/PRAE4758_combined.R1.fq
	cat PRAE4758*2.fq > ../concatenated_raw_reads2/PRAE4758_combined.R2.fq
	
	cat RUBI168*1.fq > ../concatenated_raw_reads2/RUBI168_combined.R1.fq
	cat RUBI168*2.fq > ../concatenated_raw_reads2/RUBI168_combined.R2.fq
	
	cat RUBI446*1.fq > ../concatenated_raw_reads2/RUBI446_combined.R1.fq
	cat RUBI446*2.fq > ../concatenated_raw_reads2/RUBI446_combined.R2.fq
	
	cat RUBI57579*1.fq > ../concatenated_raw_reads2/RUBI57579_combined.R1.fq
	cat RUBI57579*2.fq > ../concatenated_raw_reads2/RUBI57579_combined.R2.fq
	
	cat RUBI6012*1.fq > ../concatenated_raw_reads2/RUBI6012_combined.R1.fq
	cat RUBI6012*2.fq > ../concatenated_raw_reads2/RUBI6012_combined.R2.fq
	
	cat RUBI8047*1.fq > ../concatenated_raw_reads2/RUBI8047_combined.R1.fq
	cat RUBI8047*2.fq > ../concatenated_raw_reads2/RUBI8047_combined.R2.fq
	
	cat RUBI8184*1.fq > ../concatenated_raw_reads2/RUBI8184_combined.R1.fq
	cat RUBI8184*2.fq > ../concatenated_raw_reads2/RUBI8184_combined.R2.fq
	
	cat SCHO2588*1.fq > ../concatenated_raw_reads2/SCHO2588_combined.R1.fq
	cat SCHO2588*2.fq > ../concatenated_raw_reads2/SCHO2588_combined.R2.fq
	
	cat SCHO3258*1.fq > ../concatenated_raw_reads2/SCHO3258_combined.R1.fq
	cat SCHO3258*2.fq > ../concatenated_raw_reads2/SCHO3258_combined.R2.fq
	
	cat SCHO4231*1.fq > ../concatenated_raw_reads2/SCHO4231_combined.R1.fq
	cat SCHO4231*2.fq > ../concatenated_raw_reads2/SCHO4231_combined.R2.fq
	
	cat SCHO7708*1.fq > ../concatenated_raw_reads2/SCHO7708_combined.R1.fq
	cat SCHO7708*2.fq > ../concatenated_raw_reads2/SCHO7708_combined.R2.fq
	
	cat STEM1203*1.fq > ../concatenated_raw_reads2/STEM1203_combined.R1.fq
	cat STEM1203*2.fq > ../concatenated_raw_reads2/STEM1203_combined.R2.fq
	
	cat STEMONO8*1.fq > ../concatenated_raw_reads2/STEMONO8_combined.R1.fq
	cat STEMONO8*2.fq > ../concatenated_raw_reads2/STEMONO8_combined.R2.fq
	
	cat TESS696*1.fq > ../concatenated_raw_reads2/TESS696_combined.R1.fq
	cat TESS696*2.fq > ../concatenated_raw_reads2/TESS696_combined.R2.fq
	
	cat UNI2122*1.fq > ../concatenated_raw_reads2/UNI2122_combined.R1.fq
	cat UNI2122*2.fq > ../concatenated_raw_reads2/UNI2122_combined.R2.fq
	
	cat VEN53514*1.fq > ../concatenated_raw_reads2/VEN53514_combined.R1.fq
	cat VEN53514*2.fq > ../concatenated_raw_reads2/VEN53514_combined.R2.fq
	
# zip raw reads
	
	cd Raw_data2
	nohup gzip *.fq &
	
	cd Raw_data
	nohup gzip *.fq &
	
# zip concatenated reads

	cd concatenated_raw_reads2
	nohup gzip *.fq &

	
# Clean with Fastp

	cd Eperua
	mkdir cleaned_reads2
	
	# test run

	for file in /data/eaf236/Eperua/concatenated_raw_reads2/*R1.fq.gz
	do
		name=`basename $file .R1.fq.gz`
		echo "Running fastp on $name"
		forward=$name".R1.fq.gz"
		reverse=$name".R2.fq.gz"
		fastp -i concatenated_raw_reads2/$forward -o cleaned_reads2/$forward -I concatenated_raw_reads2/$reverse -O cleaned_reads2/$reverse -z 4 -q 20 --trim_poly_g --length_required 75 --thread 4
	done
	
	# exclude what was made
	
	# do nohup, for this, first create an .sh file with the comand above and put it in Eperua directory
	
	cd Eperua
	nohup sh fastp3.sh > fastp.log &
	
	# doing for PRAE869 (something went wrong in the concatenation)
	fastp -i PRAE869_combined.R1.fq.gz -o ../cleaned_reads2/PRAE869_combined.R1.fq.gz -I PRAE869_combined.R2.fq.gz -O ../cleaned_reads2/PRAE869_combined.R2.fq.gz -z 4 -q 20 --trim_poly_g --length_required 75 --thread 4

	cd cleaned_reads2
	gunzip PRAE869*.gz
	
	
# unzip cleaned reads cleaned_reads2
	cd cleaned_reads2
	gunzip *.gz

	
#################
###12 May 2023###
#################

# Hybpiper

	# first exclude CERRA01 CARDE2632 MANA194b FROE1458 PRAE869 OBTU207 from directory Hybpiper
		rm -r CARD2632
        rm -r CERR01
        rm -r MANA194B
        rm -r FROE1458
        rm -r PRAE869
        rm -r OBTU207
	# put the Detarioideae_baits.fasta file in Hybpiper directory
	# create the name list and put it in the Hybpiper directory
		ls > namelist2.txt # (it was cleaned manually)
		
	# running hybpiper
			
		cd Hybpiper
		conda activate hybpiper
					
		while read name; 
	do hybpiper assemble -t_dna Detarioideae_baits.fasta  -r /data/eaf236/Eperua/cleaned_reads2/$name*.fq --prefix $name --bwa --cpu 12;
	done < namelist2.txt
	
	#VEN53514 didin`t work, I don`t know why, try again individually
	hybpiper assemble -t_dna Detarioideae_baits.fasta -r /data/eaf236/Eperua/cleaned_reads2/VEN53514_combined.R*.fq --prefix VEN53514 --bwa --cpu 12
	
	#PRAE869 didin`t work, error in concatenation and fastp
	hybpiper assemble -t_dna Detarioideae_baits.fasta -r /data/eaf236/Eperua/cleaned_reads2/PRAE869_combined.R*.fq --prefix PRAE869 --bwa --cpu 12
	
	
	# hybpiper stats DO for all samples
	
	hybpiper stats -t_dna Detarioideae_baits.fasta gene namelist_all.txt
	
	# recovery heatmap for all samples
	
	hybpiper recovery_heatmap seq_lengths.tsv 
	
	# retrieve sequences all samples
	
	mkdir retrieved_sequences_unaligned_all
	
	hybpiper retrieve_sequences dna -t_dna Detarioideae_baits.fasta --sample_names namelist_all.txt --fasta_dir /data/eaf236/Eperua/retrieved_sequences_unaligned_all
	
	
# align with Maft

		# testing the comand 
		# remove all files in directory Maft
		
		cd Eperua
		conda activate Eperua_tree
		
		for file in /data/eaf236/Eperua/retrieved_sequences_unaligned_all/DN*.FNA 
		do
			name=`basename $file .FNA`
			mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > /data/eaf236/Eperua/Maft/$name.mafft.fasta
		done
	
		# nohup create a .sh. file with the comand above and put in directory Eperua
		
		nohup sh maft_all.sh > maft_all.log & 
	
#run trimal to cleanup the alignments
		
		#remove all files in directory Trimal 
		#remove exons with 0 coverage (exclude files with 0 size) in directory Maft: 48 files

	for file in /data/eaf236/Eperua/Maft/*.mafft.fasta
	do
		name=`basename $file .mafft.fasta`
		trimal -in $file -out /data/eaf236/Eperua/Trimal/$name.trimAL.mafft.fasta -automated1
	done
	
		#everything went ok, but the message Error: the symbol 'n' accesing the matrix is not defined in this object
		
		for file in /data/eaf236/Eperua/Maft/*.mafft.fasta
	do
		name=`basename $file .mafft.fasta`
		trimal -in $file -out /data/eaf236/Eperua/Trimal/$name.trimAL.mafft.fas -automated1
	done
	
	
			for file in /data/eaf236/Eperua/Maft/*.mafft.fasta
	do
		name=`basename $file .mafft.fasta`
		trimal -in $file -out /data/eaf236/Eperua/Trimal/$name.trimAL.mafft.nexus -automated1
	done
	
# Raxml to reconstruct trees from all exons individually

	for file in /data/eaf236/Eperua/Trimal/*.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxml-ng --all --msa /data/eaf236/Eperua/Trimal/$name --msa-format FASTA --model GTR+G --prefix $name --threads 6 --bs-trees 100;
        	done
        
     #nohup: create a .sh file with the command above, put in Trees_raxml
     
     nohup sh raxml_individual_genes_nohup.sh > raxml_individual_genes_nohup.log &
 
     
	
#################
###14 May 2023###
#################

# catfasta2phyml for all exons concatenated

	catfasta2phyml -c /data/eaf236/Eperua/Trimal/DN*.fasta > concatenated_exons.phy 2> partitions.txt
	
	
#partition finder
	#preparing CFG file:
		#exclude "/data/eaf236/Eperua/Trimal/" before exons name in file 'partitions.txt'
		#add a `;` after each partition using excel, exclude .maft.trimal
		#copy and paste the partitions block in the CFG file
	#paste this CFG file in the same folder as your alignment
	
	#Test without nohup
	
	conda activate PartitionFinder2
	python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/ -p 20 -- raxml
	
	# Nohup
	  
	nohup python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/ -p 20 -- raxml & 

		
#################
###15 May 2023###
#################

# Raxml to reconstruct trees from all concatenated exons without defining partitions
	mkdir Trees_raxml_concat
	
	raxml-ng --all --msa /data/eaf236/Eperua/Trees_raxml_concat/concatenated_exons.phy --msa-format PHYLIP --model GTR+G --prefix concatenated_exons --threads 12 --bs-trees 1000	#didn`t work
	
	#reduce threads because it is in the 80%
	raxml-ng --all --msa /data/eaf236/Eperua/Trees_raxml_concat/concatenated_exons.phy --model GTR+G --prefix concatenated_exons --threads 6 --bs-trees 1000	#worked, but need check the algiment
	
	#nohup: create a .sh file with the command above, put in Trees_raxml_concat
	nohup  sh raxml_concat_nohup.sh > raxml_concat_nohup.log &
	
# Species tree with Astral

	conda install -c bioconda astral-tree
	
	mkdir Trees_astral
	
	cd Trees_raxml
	cat *.bestTree > /data/eaf236/Eperua/Trees_astral/all_individual_BestTrees.tre
	ls *.bootstraps > /data/eaf236/Eperua/Trees_astral/bootstrap_trees.txt
	
	#put the .bootstraps in the Trees_astral directory
	 cp /data/eaf236/Eperua/Trees_raxml/*.bootstraps /data/eaf236/Eperua/Trees_astral/
	
	cd Trees_astral
	astral -i all_individual_BestTrees.tre -b bootstrap_trees.txt -o Astral_bootstrap_analysis.tre #100 bootstraps by defaulty
	
	#nohup
	 #nohup: create a .sh file with the command above, put in Trees_astral
	nohup  sh astral_nohup.sh > astral_nohup.log &
	
	#astral no bp
	astral -i all_individual_BestTrees.tre -o Astral_all_exon_all_taxa_analysis.tre
	
# GOTREE
	mkdir Trees_support-raxml
	cp /data/eaf236/Eperua/Trees_raxml/*.support /data/eaf236/Eperua/Trees_support-raxml/
	
 	mkdir summary_stats
 	mkdir collapsed_trees
 	mkdir rooted_trees
 	mkdir well_supported
 	mkdir larger_than70_nodes

	#reroot trees
	
	cd Trees_support-raxml
	
	for file in /data/eaf236/Eperua/Trees_support-raxml/*.support
	do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "GUIB16798"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
	done

	#trees that don`t have  GUIB16798
	# get the files with very low size in rooted_trees directory, copy to my computer, make a list with the exons names, in excel write the comand usinf the function CONCAT, after copy to to here
	# Order to choose an outgroup availalble:  HYME10482, COPA1661, SIND837, Augoardia, Stemonocoleus, Eurypetalum, Eperua purpurea (checked in seqlenghts [hybpiper])
		
		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN13548_c0_g1_i1_m45854_exon1.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN17269_c0_g2_i1_m19008_exon2.trimAL.mafft.fasta.raxml.support 
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok! (very weird exon)
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN18086_c0_g2_i2_m30440_exon2.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN19054_c0_g1_i1_m2815_exon11.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		# Ok! Few tips!!!!!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN19084_c0_g1_i1_m2866_exon6.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "PURP1973"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done                         

		#Ok!	
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN19502_c0_g1_i1_m17231_exon19.trimAL.mafft.fasta.raxml.support                  
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN19647_c0_g2_i1_m40560_exon6.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGOA4"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN19886_c0_g1_i1_m46583_exon3.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN20944_c0_g1_i1_m40434_exon4.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN21520_c0_g1_i3_m44313_exon2.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN21520_c0_g1_i3_m44313_exon4.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN21902_c0_g1_i1_m30132_exon6.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN21902_c0_g1_i1_m30132_exon7.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN22546_c0_g1_i1_m60856_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN22546_c0_g1_i1_m60856_exon3.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN23800_c0_g1_i1_m17566_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN24663_c0_g1_i1_m33180_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN24846_c0_g1_i1_m32341_exon3.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN25570_c0_g1_i2_m58578_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN26351_c0_g1_i1_m23568_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN26351_c0_g1_i1_m23568_exon3.trimAL.mafft.fasta.raxml.support                   	
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN26461_c0_g1_i1_m50737_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN27486_c0_g1_i1_m22242_exon4.trimAL.mafft.fasta.raxml.support                   	
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "STEM695"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN28031_c0_g1_i2_m59063_exon1.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok! Few taxa! Just Eperua!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN28528_c0_g1_i1_m11259_exon4.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "OLEI3209"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN28725_c0_g2_i1_m46401_exon8.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN28783_c0_g2_i2_m46468_exon2.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29041_c0_g1_i1_m6521_exon1.trimAL.mafft.fasta.raxml.support                    
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29365_c1_g2_i3_m40788_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29365_c1_g2_i3_m40788_exon2.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29386_c0_g1_i1_m40845_exon5.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29501_c0_g3_i1_m39522_exon7.trimAL.mafft.fasta.raxml.support                
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29585_c2_g1_i2_m39452_exon7.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29746_c0_g1_i1_m43048_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN29865_c0_g3_i4_m54356_exon10.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "STEM695"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN30117_c0_g3_i2_m57312_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "SIND837"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN30459_c0_g1_i1_m16208_exon10.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN30845_c0_g1_i2_m7227_exon1.trimAL.mafft.fasta.raxml.support                    
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN30871_c0_g2_i2_m7374_exon12.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN30871_c0_g2_i2_m7374_exon19.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31690_c1_g10_i1_m18498_exon1.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31690_c1_g10_i1_m18498_exon5.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31703_c1_g1_i5_m31328_exon4.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31870_c0_g3_i1_m26604_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31922_c0_g1_i1_m1232_exon13.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN31922_c0_g1_i1_m1232_exon18.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN32007_c3_g1_i1_m19432_exon2.trimAL.mafft.fasta.raxml.support                 
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "STEM695"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN32117_c1_g2_i1_m56954_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO2974"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok! 
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN328_c0_g1_i1_m12305_exon1.trimAL.mafft.fasta.raxml.support                     
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN34771_c0_g1_i1_m28451_exon10.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok! 
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN3533_c0_g1_i1_m33975_exon3.trimAL.mafft.fasta.raxml.support                    
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
		#Ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN38841_c0_g1_i1_m50501_exon8.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done

		#ok!
		for file in /data/eaf236/Eperua/Trees_support-raxml/DN47367_c0_g1_i1_m32833_exon1.trimAL.mafft.fasta.raxml.support
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "UNIJ1892"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees/$name.rooted.treefile.tre 
		done
		
#collapse weak nodes
		cd Eperua
		
		for file in /data/eaf236/Eperua/rooted_trees/*.treefile.tre
		do
			name=`basename $file .treefile.tre`
			echo "Collapsing weakly supported nodes for $name"
			gotree collapse support -i $file -s 70 -o /data/eaf236/Eperua/collapsed_trees/$name.collapsed.treefile.tre 
		done

		#Gotree: basic command for one tree using a 70% cutoff
		for file in /data/eaf236/Eperua/collapsed_trees/*.tre
		do
			name=`basename $file .tre`
			echo "Collapsing weakly supported nodes for $name"
			gotree stats edges -i $file > /data/eaf236/Eperua/summary_stats/$name.summary_file.txt
		done

		#Gotree: select only internal nodes (not polytomy collapsed nodes)
		for file in /data/eaf236/Eperua/summary_stats/*.summary_file.txt
		do
			name=`basename $file .min4.fasta.one_accession.summary_file.txt`
			awk '$4 != "N/A" && $5 == "false"' $file > /data/eaf236/Eperua/well_supported/$name.wellsupported.txt
		done

		#Gotree: double check that the correct nodes that are well-supported were saved (in this case anything with a bootstrap value of 70 or higher)
		for file in /data/eaf236/Eperua/well_supported/*.wellsupported.txt
		do
			name=`basename $file .wellsupported.txt`
			awk '$4 > 69' $file > /data/eaf236/Eperua/larger_than70_nodes/$name.larger_than70nodes.txt
		done

		#Gotree: Count number of nodes that are well-supported by counting the lines with text, each node will be its own line
		for file in /data/eaf236/Eperua/larger_than70_nodes/*.larger_than70nodes.txt
		do
			wc -l $file >> larger_than_70_nodes.txt
		done

# Look for recombination with PhiPack

		mkdir good_exons_alignments
		
		#do a list of the good exons (more than 30 well-supported nodes) using list from gotree (to write the comand I used CONCAT function in excel)
		# get the good exons from Trimal directory and put in good_exons_alignments
		
		cp /data/eaf236/Eperua/Trimal/DN13752_c0_g1_i1_m46903_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN17269_c0_g2_i1_m19008_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN17760_c0_g1_i1_m21669_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN17760_c0_g1_i1_m21669_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN18377_c0_g1_i1_m51274_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN18731_c0_g1_i1_m25776.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN19192_c0_g1_i1_m38288_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN19342_c0_g1_i1_m23149.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN19502_c0_g1_i1_m17231_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN20864_c0_g1_i1_m41593_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN21145_c0_g1_i2_m47760_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN21885_c0_g1_i1_m39862_exon4.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN22484_c0_g2_i1_m29446_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN23392_c0_g1_i1_m41339_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN23392_c0_g1_i1_m41339_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN23406_c0_g1_i2_m34414_exon8.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN23425_c0_g1_i1_m34571.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN23800_c0_g1_i1_m17566_exon9.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN24512_c1_g1_i1_m28144.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN24820_c0_g1_i1_m32355_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN24958_c0_g1_i1_m15171_exon4.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN24986_c0_g1_i1_m15066_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25026_c0_g1_i1_m45398_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25037_c0_g2_i2_m45353_exon12.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25038_c0_g1_i2_m45244_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25216_c0_g1_i1_m16555_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25586_c0_g1_i1_m58685.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN25707_c0_g1_i1_m49445_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26153_c0_g1_i1_m6090_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26377_c0_g1_i1_m23549_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26464_c0_g1_i1_m50719_exon12.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26464_c0_g1_i1_m50719_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26562_c0_g1_i1_m43840_exon6.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26693_c0_g2_i1_m35133_exon4.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN26821_c0_g1_i1_m50297_exon14.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN27315_c1_g1_i1_m52862.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN27437_c0_g1_i1_m22231.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN27578_c0_g1_i1_m4004.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN27817_c0_g1_i3_m29148.*.fasta /data/eaf236/Eperua/good_exons_alignments
		
		cp /data/eaf236/Eperua/Trimal/DN27948_c0_g1_i1_m13099_exon6.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN27992_c0_g1_i2_m12824_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN28417_c0_g1_i1_m25084_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN28518_c0_g1_i1_m11012_exon6.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN28715_c0_g1_i1_m46539.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN28749_c0_g1_i2_m46219_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN28989_c0_g2_i1_m1747_exon8.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29041_c0_g1_i1_m6521_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29099_c1_g2_i3_m6395_exon6.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29203_c0_g1_i1_m20275_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29295_c0_g2_i1_m20177_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29365_c1_g2_i3_m40788_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29442_c1_g1_i1_m37877_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29553_c0_g1_i1_m39616_exon9.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29620_c0_g1_i1_m58040_exon11.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29746_c0_g1_i1_m43048_exon5.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29763_c0_g1_i1_m43193.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN29972_c3_g2_i1_m9303_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30085_c1_g1_i1_m4862_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30085_c1_g1_i1_m4862_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30209_c0_g1_i3_m7946_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30315_c0_g1_i1_m59988.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30316_c1_g1_i1_m60090.*.fasta /data/eaf236/Eperua/good_exons_alignments
		
		cp /data/eaf236/Eperua/Trimal/DN30640_c0_g1_i1_m13938_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN30794_c2_g6_i1_m38726_exon5.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31125_c0_g1_i1_m41749_exon4.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31150_c2_g1_i1_m42164_exon10.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31203_c2_g2_i2_m3107_exon15.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31371_c0_g2_i1_m44822_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31739_c1_g3_i1_m31072_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31739_c1_g3_i1_m31072_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31870_c0_g3_i1_m26604_exon6.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN31870_c0_g3_i1_m26604_exon8.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN32114_c2_g1_i1_m56266_exon9.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN32152_c0_g1_i1_m56670_exon5.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN32205_c0_g2_i1_m36481_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN32229_c1_g5_i1_m37111.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN32328_c0_g1_i1_m5338.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN39025_c0_g1_i1_m8629_exon3.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN3971_c0_g2_i1_m30057_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN42923_c0_g1_i1_m42674_exon1.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN42923_c0_g1_i1_m42674_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN47285_c0_g1_i1_m29628.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN47306_c0_g1_i1_m32839_exon7.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN5812_c0_g1_i1_m13648_exon2.*.fasta /data/eaf236/Eperua/good_exons_alignments
		cp /data/eaf236/Eperua/Trimal/DN9107_c0_g2_i1_m9012_exon4.*.fasta /data/eaf236/Eperua/good_exons_alignments

		#run phipack
		mkdir phipack
		cd phipack
		
		for file in /data/eaf236/Eperua/good_exons_alignments/*.fasta
		do
			name=`basename $file fasta`
			echo "running PhiPack on $name"
			Phi -f $file -o -p
			mv  Phi.log $name.Phi.log
		done
		
		#significance values 
		awk 'FNR==27 {print FILENAME,$0}' /data/eaf236/Eperua/phipack/*.log > PhiPack_signficance_thresholds.txt #FNR 27 means the line 27
		
		                                                                                                                         
#################
###16 May 2023###
#################

# Analising the good exons (>=30 nodes well-supported)

	cd Trees_raxml_support
	
	#used gotree result to get the list of exons with >=30 nodes well-supported
	cat DN13752_c0_g1_i1_m46903_exon2.trimAL.mafft.fasta.raxml.support DN17269_c0_g2_i1_m19008_exon3.trimAL.mafft.fasta.raxml.support DN17760_c0_g1_i1_m21669_exon1.trimAL.mafft.fasta.raxml.support DN17760_c0_g1_i1_m21669_exon2.trimAL.mafft.fasta.raxml.support DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta.raxml.support DN18731_c0_g1_i1_m25776.trimAL.mafft.fasta.raxml.support DN19192_c0_g1_i1_m38288_exon1.trimAL.mafft.fasta.raxml.support DN19342_c0_g1_i1_m23149.trimAL.mafft.fasta.raxml.support DN19502_c0_g1_i1_m17231_exon1.trimAL.mafft.fasta.raxml.support DN20864_c0_g1_i1_m41593_exon1.trimAL.mafft.fasta.raxml.support DN21145_c0_g1_i2_m47760_exon1.trimAL.mafft.fasta.raxml.support DN21885_c0_g1_i1_m39862_exon4.trimAL.mafft.fasta.raxml.support DN22484_c0_g2_i1_m29446_exon1.trimAL.mafft.fasta.raxml.support DN23392_c0_g1_i1_m41339_exon1.trimAL.mafft.fasta.raxml.support DN23392_c0_g1_i1_m41339_exon2.trimAL.mafft.fasta.raxml.support DN23406_c0_g1_i2_m34414_exon8.trimAL.mafft.fasta.raxml.support DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta.raxml.support DN23800_c0_g1_i1_m17566_exon9.trimAL.mafft.fasta.raxml.support DN24512_c1_g1_i1_m28144.trimAL.mafft.fasta.raxml.support DN24820_c0_g1_i1_m32355_exon2.trimAL.mafft.fasta.raxml.support DN24958_c0_g1_i1_m15171_exon4.trimAL.mafft.fasta.raxml.support DN24986_c0_g1_i1_m15066_exon1.trimAL.mafft.fasta.raxml.support DN25026_c0_g1_i1_m45398_exon1.trimAL.mafft.fasta.raxml.support DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta.raxml.support DN25038_c0_g1_i2_m45244_exon3.trimAL.mafft.fasta.raxml.support DN25216_c0_g1_i1_m16555_exon2.trimAL.mafft.fasta.raxml.support DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta.raxml.support DN25707_c0_g1_i1_m49445_exon3.trimAL.mafft.fasta.raxml.support DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta.raxml.support DN26377_c0_g1_i1_m23549_exon3.trimAL.mafft.fasta.raxml.support DN26464_c0_g1_i1_m50719_exon12.trimAL.mafft.fasta.raxml.support DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta.raxml.support DN26562_c0_g1_i1_m43840_exon6.trimAL.mafft.fasta.raxml.support DN26693_c0_g2_i1_m35133_exon4.trimAL.mafft.fasta.raxml.support DN26821_c0_g1_i1_m50297_exon14.trimAL.mafft.fasta.raxml.support DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta.raxml.support DN27437_c0_g1_i1_m22231.trimAL.mafft.fasta.raxml.support DN27578_c0_g1_i1_m4004.trimAL.mafft.fasta.raxml.support DN27817_c0_g1_i3_m29148.trimAL.mafft.fasta.raxml.support DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta.raxml.support DN27992_c0_g1_i2_m12824_exon7.trimAL.mafft.fasta.raxml.support DN28417_c0_g1_i1_m25084_exon1.trimAL.mafft.fasta.raxml.support DN28518_c0_g1_i1_m11012_exon6.trimAL.mafft.fasta.raxml.support DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta.raxml.support DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta.raxml.support DN28989_c0_g2_i1_m1747_exon8.trimAL.mafft.fasta.raxml.support DN29041_c0_g1_i1_m6521_exon7.trimAL.mafft.fasta.raxml.support DN29099_c1_g2_i3_m6395_exon6.trimAL.mafft.fasta.raxml.support DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta.raxml.support DN29295_c0_g2_i1_m20177_exon1.trimAL.mafft.fasta.raxml.support DN29365_c1_g2_i3_m40788_exon3.trimAL.mafft.fasta.raxml.support DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta.raxml.support DN29553_c0_g1_i1_m39616_exon9.trimAL.mafft.fasta.raxml.support DN29620_c0_g1_i1_m58040_exon11.trimAL.mafft.fasta.raxml.support DN29746_c0_g1_i1_m43048_exon5.trimAL.mafft.fasta.raxml.support DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta.raxml.support DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta.raxml.support DN30085_c1_g1_i1_m4862_exon1.trimAL.mafft.fasta.raxml.support DN30085_c1_g1_i1_m4862_exon2.trimAL.mafft.fasta.raxml.support DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta.raxml.support DN30315_c0_g1_i1_m59988.trimAL.mafft.fasta.raxml.support DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta.raxml.support DN30640_c0_g1_i1_m13938_exon1.trimAL.mafft.fasta.raxml.support DN30794_c2_g6_i1_m38726_exon5.trimAL.mafft.fasta.raxml.support DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta.raxml.support DN31150_c2_g1_i1_m42164_exon10.trimAL.mafft.fasta.raxml.support DN31203_c2_g2_i2_m3107_exon15.trimAL.mafft.fasta.raxml.support DN31371_c0_g2_i1_m44822_exon2.trimAL.mafft.fasta.raxml.support DN31739_c1_g3_i1_m31072_exon1.trimAL.mafft.fasta.raxml.support DN31739_c1_g3_i1_m31072_exon7.trimAL.mafft.fasta.raxml.support DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta.raxml.support DN31870_c0_g3_i1_m26604_exon8.trimAL.mafft.fasta.raxml.support DN32114_c2_g1_i1_m56266_exon9.trimAL.mafft.fasta.raxml.support DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta.raxml.support DN32205_c0_g2_i1_m36481_exon3.trimAL.mafft.fasta.raxml.support DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta.raxml.support DN32328_c0_g1_i1_m5338.trimAL.mafft.fasta.raxml.support DN39025_c0_g1_i1_m8629_exon3.trimAL.mafft.fasta.raxml.support DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta.raxml.support DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta.raxml.support DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta.raxml.support DN47285_c0_g1_i1_m29628.trimAL.mafft.fasta.raxml.support DN47306_c0_g1_i1_m32839_exon7.trimAL.mafft.fasta.raxml.support DN5812_c0_g1_i1_m13648_exon2.trimAL.mafft.fasta.raxml.support DN9107_c0_g2_i1_m9012_exon4.trimAL.mafft.fasta.raxml.support > /data/eaf236/Eperua/good_exons.raxml.support.tre
	

#################
###19 May 2023###
#################
	
#Next steps
	#Exons with > 30 nodes well-supported: Tree (concatenated, Astral) from the good exons (signal for no recombination) and from the "bad" exons (signal for recombination): 
		#after run Phyparts to quantify the signal.
	# Concatenated and species tree from all exons excluding those taxa with low exon coverage and exons with low coverage:
			#GLBA499
			#PRAE1000
			#RUBI23804
			#TESS3457
			
	#Testing exons with >= 30 no nodes well-supported

# Preparing directory for next analysis

	#Directory for Astral
		mkdir good_exons_trees #with at least 30 nodes well-supported
		cd Trees_raxml
		
		cp DN13752_c0_g1_i1_m46903_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN17269_c0_g2_i1_m19008_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN17760_c0_g1_i1_m21669_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN17760_c0_g1_i1_m21669_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN18731_c0_g1_i1_m25776.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN19192_c0_g1_i1_m38288_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN19342_c0_g1_i1_m23149.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN19502_c0_g1_i1_m17231_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN20864_c0_g1_i1_m41593_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN21145_c0_g1_i2_m47760_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN21885_c0_g1_i1_m39862_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN22484_c0_g2_i1_m29446_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN23392_c0_g1_i1_m41339_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN23392_c0_g1_i1_m41339_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN23406_c0_g1_i2_m34414_exon8.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN23800_c0_g1_i1_m17566_exon9.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN24512_c1_g1_i1_m28144.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN24820_c0_g1_i1_m32355_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN24958_c0_g1_i1_m15171_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN24986_c0_g1_i1_m15066_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25026_c0_g1_i1_m45398_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25038_c0_g1_i2_m45244_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25216_c0_g1_i1_m16555_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN25707_c0_g1_i1_m49445_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26377_c0_g1_i1_m23549_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26464_c0_g1_i1_m50719_exon12.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		
		cp DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26562_c0_g1_i1_m43840_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26693_c0_g2_i1_m35133_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN26821_c0_g1_i1_m50297_exon14.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27437_c0_g1_i1_m22231.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27578_c0_g1_i1_m4004.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27817_c0_g1_i3_m29148.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN27992_c0_g1_i2_m12824_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN28417_c0_g1_i1_m25084_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN28518_c0_g1_i1_m11012_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN28989_c0_g2_i1_m1747_exon8.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29041_c0_g1_i1_m6521_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29099_c1_g2_i3_m6395_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29295_c0_g2_i1_m20177_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29365_c1_g2_i3_m40788_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29553_c0_g1_i1_m39616_exon9.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29620_c0_g1_i1_m58040_exon11.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29746_c0_g1_i1_m43048_exon5.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30085_c1_g1_i1_m4862_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30085_c1_g1_i1_m4862_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30315_c0_g1_i1_m59988.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30640_c0_g1_i1_m13938_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN30794_c2_g6_i1_m38726_exon5.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		
		cp DN31150_c2_g1_i1_m42164_exon10.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31203_c2_g2_i2_m3107_exon15.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31371_c0_g2_i1_m44822_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31739_c1_g3_i1_m31072_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31739_c1_g3_i1_m31072_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN31870_c0_g3_i1_m26604_exon8.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN32114_c2_g1_i1_m56266_exon9.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN32205_c0_g2_i1_m36481_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN32328_c0_g1_i1_m5338.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN39025_c0_g1_i1_m8629_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN47285_c0_g1_i1_m29628.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN47306_c0_g1_i1_m32839_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN5812_c0_g1_i1_m13648_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees
		cp DN9107_c0_g2_i1_m9012_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees

		cd good_exons_trees
		mkdir exons_trees_signal_for_recombination #exons with p-value less than 0.00059 moving to a subdirectory
		
				mv DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta.raxml.bestTree /data/eaf236/Eperua/good_exons_trees/exons_trees_signal_for_recombination

		
		# Directory for concatenated analysis
	
		# In the directory 'good_exons_alignments' I have all the exons with at least 30 nodes well-supported.
		# I`ll create a subdirectory called 'exons_signal_for_recombination' and move the ones with signal for recombination to this sub-directory.
		
		cd good_exons_alignments
		mkdir exons_signal_for_recombination
				
		mv DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination

# Species tree with Astral: good exons (>=30 nodes well-supported) with recombination AND good exons (>=30 nodes well-supported) no recombination

		mkdir Trees_astral_good_exons
	
	#copy the .bootstraps of the good exons (all with >=30 nodes well-supported) to the Trees_astral_good_exons directory
	
		cd Trees_raxml
		
		cp DN13752_c0_g1_i1_m46903_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN17269_c0_g2_i1_m19008_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN17760_c0_g1_i1_m21669_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN17760_c0_g1_i1_m21669_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN18731_c0_g1_i1_m25776.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN19192_c0_g1_i1_m38288_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN19342_c0_g1_i1_m23149.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN19502_c0_g1_i1_m17231_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN20864_c0_g1_i1_m41593_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN21145_c0_g1_i2_m47760_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN21885_c0_g1_i1_m39862_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN22484_c0_g2_i1_m29446_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN23392_c0_g1_i1_m41339_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN23392_c0_g1_i1_m41339_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN23406_c0_g1_i2_m34414_exon8.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN23800_c0_g1_i1_m17566_exon9.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN24512_c1_g1_i1_m28144.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN24820_c0_g1_i1_m32355_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN24958_c0_g1_i1_m15171_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN24986_c0_g1_i1_m15066_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN25026_c0_g1_i1_m45398_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN25038_c0_g1_i2_m45244_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN25216_c0_g1_i1_m16555_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		
		cp DN25707_c0_g1_i1_m49445_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26377_c0_g1_i1_m23549_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26464_c0_g1_i1_m50719_exon12.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26562_c0_g1_i1_m43840_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26693_c0_g2_i1_m35133_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN26821_c0_g1_i1_m50297_exon14.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27437_c0_g1_i1_m22231.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27578_c0_g1_i1_m4004.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27817_c0_g1_i3_m29148.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN27992_c0_g1_i2_m12824_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN28417_c0_g1_i1_m25084_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN28518_c0_g1_i1_m11012_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN28989_c0_g2_i1_m1747_exon8.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29041_c0_g1_i1_m6521_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29099_c1_g2_i3_m6395_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29295_c0_g2_i1_m20177_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29365_c1_g2_i3_m40788_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29553_c0_g1_i1_m39616_exon9.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29620_c0_g1_i1_m58040_exon11.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29746_c0_g1_i1_m43048_exon5.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		
		cp DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30085_c1_g1_i1_m4862_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30085_c1_g1_i1_m4862_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30315_c0_g1_i1_m59988.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30640_c0_g1_i1_m13938_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN30794_c2_g6_i1_m38726_exon5.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31150_c2_g1_i1_m42164_exon10.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31203_c2_g2_i2_m3107_exon15.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31371_c0_g2_i1_m44822_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31739_c1_g3_i1_m31072_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31739_c1_g3_i1_m31072_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN31870_c0_g3_i1_m26604_exon8.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN32114_c2_g1_i1_m56266_exon9.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN32205_c0_g2_i1_m36481_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN32328_c0_g1_i1_m5338.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN39025_c0_g1_i1_m8629_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN47285_c0_g1_i1_m29628.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN47306_c0_g1_i1_m32839_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN5812_c0_g1_i1_m13648_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
		cp DN9107_c0_g2_i1_m9012_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons
	 
	#move the bootstraps from the exons with recombination to a new subdirectory
	
		cd Trees_astral_good_exons
		mkdir Trees_astral_good_exons_signal_for_recombination
		
		mv DN18377_c0_g1_i1_m51274_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN23425_c0_g1_i1_m34571.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN25037_c0_g2_i2_m45353_exon12.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN25586_c0_g1_i1_m58685.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN26153_c0_g1_i1_m6090_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN26464_c0_g1_i1_m50719_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN27315_c1_g1_i1_m52862.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN27948_c0_g1_i1_m13099_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN28715_c0_g1_i1_m46539.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN28749_c0_g1_i2_m46219_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN29203_c0_g1_i1_m20275_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN29442_c1_g1_i1_m37877_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN29763_c0_g1_i1_m43193.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN29972_c3_g2_i1_m9303_exon3.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN30209_c0_g1_i3_m7946_exon7.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN30316_c1_g1_i1_m60090.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN31125_c0_g1_i1_m41749_exon4.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN31870_c0_g3_i1_m26604_exon6.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN32152_c0_g1_i1_m56670_exon5.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN32229_c1_g5_i1_m37111.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN3971_c0_g2_i1_m30057_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon1.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination
		mv DN42923_c0_g1_i1_m42674_exon2.trimAL.mafft.fasta.raxml.bootstraps /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination

		#preparing bootstrap.txt files
		#no recombination
		cd Trees_astral_good_exons
		ls *.bootstraps > /data/eaf236/Eperua/Trees_astral_good_exons/good_exons_no_recombination_bootstrap_trees.txt
		#with recombination
		cd Trees_astral_good_exons_signal_for_recombination
		ls *.bootstraps > /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination/good_exons_with_recombination_bootstrap_trees.txt
		
		#catting trees no recombination
		cd good_exons_trees 
		cat *.bestTree > /data/eaf236/Eperua/Trees_astral_good_exons/good_exons_no_recombination_BestTrees.tre
		#catting trees with recombination
		cd  exons_trees_signal_for_recombination #(in directory good_exons_trees) 
		cat *.bestTree > /data/eaf236/Eperua/Trees_astral_good_exons/Trees_astral_good_exons_signal_for_recombination/good_exons_WITH_recombination_BestTrees.tre
	
#################
###22 May 2023###
#################

		#Test running Astral
		#running Astral for good exons/no recombination
		cd Trees_astral_good_exons
		astral -i good_exons_no_recombination_BestTrees.tre -b good_exons_no_recombination_bootstrap_trees.txt -o good_exons_no_recombination_Astral_bootstrap_analysis.tre #100 bootstraps by defaulty
		#no bootstrap
		astral -i good_exons_no_recombination_BestTrees.tre -o good_exons_no_recombination_Astral_analysis.tre #100 bootstraps by defaulty
		#running Astral for good exons/no recombination
		cd Trees_astral_good_exons_signal_for_recombination # (in directory Trees_astral_good_exons)
		astral -i good_exons_WITH_recombination_BestTrees.tre -b good_exons_with_recombination_bootstrap_trees.txt -o good_exons_WITH_recombination_Astral_bootstrap_analysis.tre #100 bootstraps by defaulty
		#no bootstrap
		astral -i good_exons_WITH_recombination_BestTrees.tre -o good_exons_WITH_recombination_Astral_analysis.tre 
		#nohup Astral
		#good exons no recombination 
		 #nohup: create a .sh file with the command above, put in Trees_astral_good_exons
		nohup  sh astral_good_exons_NO_recombination_nohup.sh > astral_good_exons_NO_recombination_nohup.log &
		#good exons WITH recombination 
		 #nohup: create a .sh file with the command above, put in Trees_astral_good_exons_signal_for_recombination 
		nohup  sh astral_good_exons_WITH_recombination_nohup.sh > astral_good_exons_WITH_recombination_nohup.log &

#Concatenated tree	
	
	# catfasta2phyml for good exons with no recombination
	cd good_exons_alignments
	catfasta2phyml -c /data/eaf236/Eperua/good_exons_alignments/DN*.fasta > concatenated__good_exons_NO_recombination.phy 2> partitions_good_exons_NO_recombination.txt
	
	# catfasta2phyml for good exons with no recombination
	cd exons_signal_for_recombination
	catfasta2phyml -c /data/eaf236/Eperua/good_exons_alignments/exons_signal_for_recombination/DN*.fasta > concatenated__good_exons_WITH_recombination.phy 2> partitions_good_exons_WITH_recombination.txt
	
	
#partition finder
	#preparing CFG file:
		#exclude "/data/eaf236/Eperua/Trimal/" before exons name in file 'partitions.txt'
		#add a `;` after each partition using excel, exclude .maft.trimal
		#copy and paste the partitions block in the CFG file
	#paste this CFG file in the same folder as your alignment
	#create sub folders in the directory 'partitionfinder'
	
	#CFG file for good_exons_aligments
	#exclude FROE32452 and 	GLBA499 (no data)
	
	cd partitionfinder
	mkdir good_exons_no_recombination
	#copy the aligment 
	cd good_exons_aligments #DONE
	cp concatenated__good_exons_NO_recombination.phy > /data/eaf236/Eperua/partition_finder/good_exons_no_recombination #did manually because I opened it in my computer
	
	mkdir good_exons_with_recombination #DONE
	#copy the aligment 
	cp concatenated__good_exons_WITH_recombination.phy > /data/eaf236/Eperua/partition_finder/good_exons_with_recombination #did manually because I opened it in my computer
	
	
	#Test without nohup
	#good exons no recombination 
	conda activate PartitionFinder2
	python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/good_exons_no_recombination/ -p 20 -- raxml
	#good exons with recombination 
	python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/good_exons_with_recombination/ -p 20 -- raxml
	
	# Nohup
	#good exons no recombination  
	nohup python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/good_exons_no_recombination/ -p 20 -- raxml & 
	nohup python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/good_exons_with_recombination/ -p 20 -- raxml &
	

#Concatenated tree good exons raxml-ng

	mkdir Trees_raxml_concat_good_exons #DONE
	# copy the concatenated aligment to this directory #did manually because I opened it in my computer (partition finder)
		cd good_exons_aligments
		cp concatenated__good_exons_NO_recombination.phy > /data/eaf236/Eperua/Trees_raxml_concat_good_exons/
		
		cd exons_signal_for_recombination #(inside good_exons_aligments)
		cp concatenated__good_exons_WITH_recombination.phy > /data/eaf236/Eperua/Trees_raxml_concat_good_exons/
	
	#good exons no recombination
	#add info of partitions
	#put .part files in Trees_raxml_concat_good_exons #worked with --msa-format PHYLIP
	raxml-ng --all --msa /data/eaf236/Eperua/Trees_raxml_concat_good_exons/concatenated__good_exons_NO_recombination.phy --msa-format PHYLIP --model concatenated__good_exons_NO_recombination.part --prefix concatenated__good_exons_NO_recombination --threads 12 --bs-trees 1000	
	 #model  GTR+FC+G4m+B
	#good exons with recombination
	#put .part files in Trees_raxml_concat_good_exons #didn`t work with --msa-format PHYLIP
	raxml-ng --all --msa /data/eaf236/Eperua/Trees_raxml_concat_good_exons/concatenated__good_exons_WITH_recombination.phy --model concatenated__good_exons_WITH_recombination.part --prefix concatenated__good_exons_WITH_recombination --threads 12 --bs-trees 1000
		#model  GTR+FC+G4m+B
	#nohup: create a .sh file with the command above, put in Trees_raxml_concat
	#No recombination 
	nohup  sh raxml_concat_good_exons_NO_recombination_nohup.sh > raxml_concat_good_exons_NO_recombination_nohup.log &
	#with recombination
	nohup  sh raxml_concat_good_exons_WITH_recombination_nohup.sh > raxml_concat_good_exons_WITH_recombination_nohup.log &
	

# Trees cleaned (excluding bad samples and bad exons)

	# retrieve sequences all samples minus  #FROE32452
											 #TESS3457 
											 #GLBA499
											 #PRAE1000
											 #RUBI23804
						
	mkdir retrieved_sequences_unaligned_minus_five
	#make new taxa list excluding the five taxa, put in 'Hybpiper' directory
	cd Hybpiper
	conda activate hybpiper
	hybpiper retrieve_sequences dna -t_dna Detarioideae_baits.fasta --sample_names namelist_minus_five.txt --fasta_dir /data/eaf236/Eperua/retrieved_sequences_unaligned_minus_five
	#exclude exons with low coverage (<50% specimens = less than 57 specimens (<=57)) #total number of taxa 114 (excluding the 5 above)
		#138 exons	excluded (used stats Hybpiper to do the list) (48 zero coverage, 90 less than 50% samples more than 1)
	cd retrieved_sequences_unaligned_minus_five
	rm DN20832_c1_g1_i1_m41674_exon10.FNA DN31284_c0_g1_i1_m3006_exon5.FNA DN15361_c0_g1_i1_m23054_exon3.FNA DN29093_c0_g3_i1_m6303_exon1.FNA DN29985_c0_g2_i1_m9227_exon7.FNA DN25657_c0_g1_i2_m24593_exon3.FNA DN25216_c0_g1_i1_m16555_exon1.FNA DN32117_c1_g2_i1_m56954_exon1.FNA DN3533_c0_g1_i1_m33975_exon3.FNA DN31460_c0_g1_i1_m55429_exon7.FNA DN26259_c0_g1_i1_m9831_exon4.FNA DN21520_c0_g1_i3_m44313_exon4.FNA DN47117_c0_g1_i1_m17272_exon1.FNA DN30871_c0_g2_i2_m7374_exon12.FNA DN21902_c0_g1_i1_m30132_exon7.FNA DN29041_c0_g1_i1_m6521_exon1.FNA DN26426_c0_g1_i3_m50699_exon5.FNA DN28430_c0_g2_i1_m25299_exon6.FNA DN27439_c0_g1_i1_m22151_exon1.FNA DN28200_c0_g1_i1_m45732_exon7.FNA DN29491_c0_g1_i1_m38000_exon2.FNA DN25037_c0_g2_i2_m45353_exon2.FNA DN26461_c0_g1_i1_m50737_exon1.FNA DN31284_c0_g1_i1_m3006_exon3.FNA DN29491_c0_g1_i1_m38000_exon1.FNA DN19175_c0_g1_i1_m38298_exon3.FNA DN30000_c1_g3_i2_m4922_exon7.FNA DN19192_c0_g1_i1_m38288_exon5.FNA DN26426_c0_g1_i3_m50699_exon4.FNA DN27372_c1_g2_i2_m52979_exon3.FNA DN26598_c0_g2_i1_m43989_exon2.FNA DN32007_c3_g1_i1_m19432_exon8.FNA DN28683_c0_g1_i1_m12370_exon1.FNA DN31690_c1_g10_i1_m18498_exon14.FNA DN31147_c0_g1_i1_m42169_exon6.FNA DN25216_c0_g1_i1_m16555_exon3.FNA DN19502_c0_g1_i1_m17231_exon19.FNA DN28725_c0_g2_i1_m46401_exon8.FNA DN32007_c3_g1_i1_m19432_exon2.FNA DN29501_c0_g3_i1_m39522_exon8.FNA DN21797_c1_g1_i1_m5766_exon6.FNA DN26627_c0_g1_i1_m35035_exon8.FNA DN51253_c0_g1_i1_m13283_exon3.FNA DN28518_c0_g1_i1_m11012_exon1.FNA DN30640_c0_g1_i1_m13938_exon3.FNA DN28133_c0_g1_i2_m29758_exon1.FNA DN7229_c0_g1_i1_m41138_exon2.FNA DN34771_c0_g1_i1_m28451_exon7.FNA DN30628_c0_g1_i1_m13947_exon7.FNA DN27486_c0_g1_i1_m22242_exon4.FNA DN31150_c2_g1_i1_m42164_exon4.FNA DN20944_c0_g1_i1_m40434_exon2.FNA DN25570_c0_g1_i2_m58578_exon19.FNA DN28518_c0_g1_i1_m11012_exon7.FNA DN21902_c0_g1_i1_m30132_exon5.FNA DN31872_c0_g3_i1_m26493_exon14.FNA DN29295_c0_g2_i1_m20177_exon6.FNA DN28031_c0_g1_i2_m59063_exon12.FNA DN9107_c0_g2_i1_m9012_exon2.FNA DN26821_c0_g1_i1_m50297_exon9.FNA DN13779_c0_g1_i1_m46957_exon3.FNA DN28528_c0_g1_i1_m11259_exon4.FNA DN21520_c0_g1_i3_m44313_exon2.FNA DN29553_c0_g1_i1_m39617_exon4.FNA DN20049_c0_g1_i1_m28600_exon3.FNA DN47367_c0_g1_i1_m32833_exon1.FNA DN27362_c0_g1_i1_m52842_exon1.FNA DN26426_c0_g1_i3_m50699_exon1.FNA DN26426_c0_g1_i3_m50699_exon6.FNA DN21902_c0_g1_i1_m30132_exon6.FNA DN5812_c0_g1_i1_m13648_exon1.FNA DN31922_c0_g1_i1_m1232_exon13.FNA DN31150_c2_g1_i1_m42164_exon9.FNA DN23406_c0_g1_i2_m34414_exon4.FNA DN31690_c1_g10_i1_m18498_exon5.FNA DN25570_c0_g1_i2_m58578_exon5.FNA DN13779_c0_g1_i1_m46957_exon5.FNA DN30209_c0_g1_i3_m7946_exon4.FNA DN26464_c0_g1_i1_m50719_exon6.FNA DN38841_c0_g1_i1_m50501_exon8.FNA DN32203_c2_g1_i1_m36628_exon1.FNA DN31703_c1_g1_i5_m31328_exon4.FNA DN24808_c0_g2_i1_m32440_exon3.FNA DN17269_c0_g2_i1_m19008_exon2.FNA DN24450_c0_g1_i1_m48340_exon5.FNA DN13548_c0_g1_i1_m45854_exon6.FNA DN26570_c0_g1_i1_m43984_exon2.FNA DN29985_c0_g2_i1_m9227_exon3.FNA DN29501_c0_g3_i1_m39522_exon1.FNA DN23800_c0_g1_i1_m17566_exon1.FNA DN24846_c0_g1_i1_m32341_exon3.FNA DN29198_c0_g1_i3_m16666_exon4.FNA DN26351_c0_g1_i1_m23568_exon1.FNA DN22289_c0_g1_i1_m40222_exon8.FNA DN19054_c0_g1_i1_m2815_exon4.FNA DN26377_c0_g1_i1_m23549_exon1.FNA DN7229_c0_g1_i1_m41138_exon9.FNA DN31159_c0_g1_i1_m42190_exon1.FNA DN30459_c0_g1_i1_m16208_exon7.FNA DN26627_c0_g1_i1_m35035_exon4.FNA DN24808_c0_g2_i1_m32440_exon1.FNA DN32007_c3_g1_i1_m19432_exon1.FNA DN21457_c0_g1_i1_m6708_exon1.FNA DN22546_c0_g1_i1_m60856_exon1.FNA DN30117_c0_g3_i2_m57312_exon1.FNA DN29172_c1_g2_i1_m16971_exon7.FNA DN30845_c0_g1_i2_m7227_exon1.FNA DN30000_c1_g3_i2_m4922_exon3.FNA DN19175_c0_g1_i1_m38298_exon4.FNA DN29865_c0_g3_i4_m54356_exon10.FNA DN29198_c0_g1_i3_m16666_exon11.FNA DN28725_c0_g2_i1_m46401_exon12.FNA DN28749_c0_g1_i2_m46219_exon4.FNA DN26259_c0_g1_i1_m9831_exon1.FNA DN18086_c0_g2_i2_m30440_exon2.FNA DN19736_c0_g1_i1_m8351_exon2.FNA DN20832_c1_g1_i1_m41674_exon9.FNA DN31371_c0_g2_i1_m44822_exon10.FNA DN28783_c0_g2_i2_m46468_exon2.FNA DN19084_c0_g1_i1_m2866_exon6.FNA DN27345_c0_g1_i1_m52944_exon1.FNA DN20832_c1_g1_i1_m41674_exon5.FNA DN27948_c0_g1_i1_m13099_exon5.FNA DN21520_c0_g1_i3_m44313_exon9.FNA DN26355_c0_g1_i1_m23542_exon3.FNA DN19054_c0_g1_i1_m2815_exon5.FNA DN29501_c0_g3_i1_m39522_exon7.FNA DN26242_c0_g1_i1_m9810_exon5.FNA DN24450_c0_g1_i1_m48340_exon12.FNA DN27849_c0_g1_i1_m29151_exon3.FNA DN24663_c0_g1_i1_m33180_exon1.FNA DN29377_c2_g1_i1_m41038_exon1.FNA DN31640_c3_g2_i1_m18201_exon3.FNA DN26080_c0_g1_i1_m27870_exon8.FNA DN38841_c0_g1_i1_m50501_exon3.FNA DN27372_c1_g2_i2_m52979_exon4.FNA DN51374_c0_g1_i1_m21506_exon1.FNA DN26430_c0_g1_i1_m50610_exon6.FNA 
 
	#Maft
		cd Eperua
		mkdir Maft_minus_five_bad_exons
		conda activate Eperua_tree
				for file in /data/eaf236/Eperua/retrieved_sequences_unaligned_minus_five/DN*.FNA 
		do
			name=`basename $file .FNA`
			mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > /data/eaf236/Eperua/Maft_minus_five_bad_exons/$name.mafft.fasta
		done
			
	#run trimal to cleanup the alignments
		cd Eperua
		mkdir Trimal_minus_five_bad_exons
		for file in /data/eaf236/Eperua/Maft_minus_five_bad_exons/*.mafft.fasta
		do
			name=`basename $file .mafft.fasta`
			trimal -in $file -out /data/eaf236/Eperua/Trimal_minus_five_bad_exons/$name.trimAL.mafft.fasta -automated1
		done	
			#Error: the symbol 'n' accesing the matrix is not defined in this object
		
		#Concatenate tree raxml
		#concatenate the exons
		# catfasta2phyml for all exons/samples minus minus_five_bad_exons
		cd Trimal_minus_five_bad_exons
		catfasta2phyml -c /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN*.fasta > concatenated_all_minus_five_bad_exons.phy 2> partitions_concatenated_all_minus_five_bad_exons.txt
	   
		#Partitionfinder
		cd partitionfinder
		mkdir all_minus_five_samples_bad_exons
		#copy the aligment to minus_five_bad_exons
		#prepare the CFG file and put in minus_five_bad_exons
		#Test without nohup
		conda activate PartitionFinder2
		python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/all_minus_five_samples_bad_exons/ -p 12 -- raxml
		# Nohup
		nohup python /data/eaf236/installed_programs/partitionfinder-2.1.1/PartitionFinder.py /data/eaf236/Eperua/partition_finder/all_minus_five_samples_bad_exons/ -p 12 -- raxml & 
			
	#raxml-ng to resconstruct individual exons trees
		mkdir Trees_raxml_concat_minus
		cd Trees_raxml_concat_minus # after 25 May 2023, I changed the name of the directory to Trees_raxml_minus
		     	
		for file in /data/eaf236/Eperua/Trimal_minus_five_bad_exons/*.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxml-ng --all --msa /data/eaf236/Eperua/Trimal_minus_five_bad_exons/$name --msa-format FASTA --model GTR+G --prefix $name --threads 6 --bs-trees 100;
        done
      #nohup: create a .sh file with the command above, put in Trees_raxml
      nohup sh raxml_individual_genes_some_excluded_nohup.sh > raxml_individual_genes_some_excluded_nohup.log &
     
      #879 resulting trees (before 883 aligments!) 4 aligments not run well! I saw it in the log file and doing ls *.TMP (aligments with error)
      #running trees individually for the aliments with error
      
      #Ok!   
      raxml-ng --all --msa /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN20058_c0_g1_i1_m28642_exon1.trimAL.mafft.fasta --msa-format FASTA --model GTR+G --prefix DN20058_c0_g1_i1_m28642_exon1.trimAL.mafft.fasta --threads 6 --bs-trees 100;
			
	  #ok!
	  raxml-ng --all --msa /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN24830_c0_g2_i1_m32410_exon5.trimAL.mafft.fasta --msa-format FASTA --model GTR+G --prefix DN24830_c0_g2_i1_m32410_exon5.trimAL.mafft.fasta --threads 6 --bs-trees 100;
	  
	  #ok
	  raxml-ng --all --msa /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN30871_c0_g2_i2_m7374_exon2.trimAL.mafft.fasta --msa-format FASTA --model GTR+G --prefix DN30871_c0_g2_i2_m7374_exon2.trimAL.mafft.fasta --threads 6 --bs-trees 100;
		
 	  #ok!
 	  raxml-ng --all --msa /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN31203_c2_g2_i2_m3107_exon2.trimAL.mafft.fasta --msa-format FASTA --model GTR+G --prefix DN31203_c2_g2_i2_m3107_exon2.trimAL.mafft.fasta --threads 6 --bs-trees 100;
		
 	  # all changed model to Model: GTR+FO+G4m

			
#################
###25 May 2023###
#################

		#Species tree Astral
			mkdir Trees_astral_all_minus_five_bad_exons #(in directory Eperua)
			#preparing bootstrap.txt files
			cd Trees_raxml_minus  
			ls *.bootstraps > /data/eaf236/Eperua/Trees_astral_all_minus_five_bad_exons/all_exons_minus_five_bad_exons_bootstrap_trees.txt
			#copy the bootstrap files 
			cd Trees_raxml_minus
			cp DN*.bootstraps /data/eaf236/Eperua/Trees_astral_all_minus_five_bad_exons/
			#preparing the bestTree file
			cd Trees_raxml_minus
			cat *.bestTree > /data/eaf236/Eperua/Trees_astral_all_minus_five_bad_exons/all_exons_minus_five_bad_exons_bestTree.tre
			#running Astral testing comand	
			cd Trees_astral_all_minus_five_bad_exons
			astral -i all_exons_minus_five_bad_exons_bestTree.tre -b all_exons_minus_five_bad_exons_bootstrap_trees.txt -o all_exons_minus_five_bad_exons_Astral_bootstrap_analysis.tre #100 bootstraps by defaulty
			#nohup Astral
			#nohup: create a .sh file with the command above, put in Trees_astral_all_minus_five_bad_exons
			nohup  sh astral_all_exons_minus_five_bad_exons_nohup.sh > astral_all_exons_minus_five_bad_exons_nohup.log &
			
#################
###28 May 2023###
#################

#Astral no bp
			astral -i all_exons_minus_five_bad_exons_bestTree.tre  -o all_exons_minus_five_bad_exons_Astral_analysis.tre
			#nohup
			nohup  sh astral_all_exons_minus_five_bad_exons_no_bp_nohup.sh > astral_all_exons_minus_five_bad_exons_no_bp_nohup.log &
			
			
			# root trees
			mkdir rooted_trees_all_minus_five_bad_exons

			cd Trees_raxml_minus
	
		for file in /data/eaf236/Eperua/Trees_raxml_minus/*.support
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "GUIB16798"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#no GUIB16798
		
		#ok! DN13548_c0_g1_i1_m45854_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN13548_c0_g1_i1_m45854_exon1.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN19054_c0_g1_i1_m2815_exon11.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19054_c0_g1_i1_m2815_exon11.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN19647_c0_g2_i1_m40560_exon6.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19647_c0_g2_i1_m40560_exon6.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "AUGOA4"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done

		# Ok! DN19886_c0_g1_i1_m46583_exon3.
			for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19886_c0_g1_i1_m46583_exon3.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN20944_c0_g1_i1_m40434_exon4.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN20944_c0_g1_i1_m40434_exon4.trimAL.mafft.fasta.raxml.support                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.support`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN22546_c0_g1_i1_m60856_exon3.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN22546_c0_g1_i1_m60856_exon3.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN25570_c0_g1_i2_m58578_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN25570_c0_g1_i2_m58578_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN26351_c0_g1_i1_m23568_exon3.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN26351_c0_g1_i1_m23568_exon3.trimAL.mafft.fasta.raxml.support                   	
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN28031_c0_g1_i2_m59063_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN28031_c0_g1_i2_m59063_exon1.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		#Ok! DN29365_c1_g2_i3_m40788_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29365_c1_g2_i3_m40788_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN29365_c1_g2_i3_m40788_exon2.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29365_c1_g2_i3_m40788_exon2.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN29386_c0_g1_i1_m40845_exon5.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29386_c0_g1_i1_m40845_exon5.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN29585_c2_g1_i2_m39452_exon7.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29585_c2_g1_i2_m39452_exon7.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN29746_c0_g1_i1_m43048_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29746_c0_g1_i1_m43048_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN30459_c0_g1_i1_m16208_exon10.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN30459_c0_g1_i1_m16208_exon10.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN30871_c0_g2_i2_m7374_exon19.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN30871_c0_g2_i2_m7374_exon19.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN31690_c1_g10_i1_m18498_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31690_c1_g10_i1_m18498_exon1.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN31870_c0_g3_i1_m26604_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31870_c0_g3_i1_m26604_exon1.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN31922_c0_g1_i1_m1232_exon18.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31922_c0_g1_i1_m1232_exon18.trimAL.mafft.fasta.raxml.support                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done
		
		# Ok! DN328_c0_g1_i1_m12305_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN328_c0_g1_i1_m12305_exon1.trimAL.mafft.fasta.raxml.support                     
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done

		
		# Ok! DN34771_c0_g1_i1_m28451_exon10.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN34771_c0_g1_i1_m28451_exon10.trimAL.mafft.fasta.raxml.support                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.support`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons/$name.rooted.treefile.tre 
		done

		#Phyparts
		#copyed the 883 rooted trees tres
		cat *.tre > Rooted_gene_trees.tre # put in directoy  rooted_trees_all_minus_five_bad_exons_ inside phyparts directory
		# root species tree (didn`t work!) root manually
		 #for file in /data/eaf236/Eperua/Trees_astral_all_minus_five_bad_exons/all_exons_minus_five_bad_exons_Astral_analysis.tre                  
		 #do
		 #name=`basename $file .trimAL.mafft.fasta.raxml.support`
		 #echo "rooting tree for $name"
		 #echo "GUIB16798"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/Trees_astral_all_minus_five_bad_exons/$name.rooted.treefile.tree 
			#done
		#copy the rooted trees an cat.tree to my computer in phyparts directory
		# put the cat file with the 883 trees in the -d / rooted tree in figtree exporting as newick format (but it keep with .tre format
		 java -jar ../Phyparts-blackrim-37067bc3bc14/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees_all_minus_five_bad_exons/ -m all_exons_minus_five_bad_exons_Astral_analysis_rooted.tre -o out_nuclear
		#worked with 883 trees in 883 files (not cat)
		java -jar ../Phyparts-blackrim-37067bc3bc14/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees_all_minus_five_bad_exons_/ -m all_exons_minus_five_bad_exons_Astral_analysis_rooted.tre -o out_exons
		#put everything in the server nohup
		nohup  sh phyparts_nohup.sh > phyparts_nohup.log &
		
		
		#Phyparts plot
		#copied the file phyparts plot from Lab Jacob directory to the phypart directory
		#changed
		#Species Tree
			#sptree_file = "all_exons_minus_five_bad_exons_Astral_analysis_rooted.tre"
			#total_genes = 883
		conda install -y matplotlib
		python phyparts_plot.py
		
	
#################
###01 June 2023###
#################	

			#Raxml-ng for the concatenated tree
			mkdir Trees_raxml_concat_all_minus_five_bad_exons
			#copy the aligment to Trees_raxml_concat_all_minus_five_bad_exons
			#prepare the .part file with name best_partitions_concatenated_all_minus_five_bad_exons.part
			conda activate Eperua_tree
			raxml-ng --all --msa /data/eaf236/Eperua/Trees_raxml_concat_all_minus_five_bad_exons/concatenated_all_minus_five_bad_exons.phy --model best_partitions_concatenated_all_minus_five_bad_exons.part --prefix concatenated_all_minus_five_bad_exons --threads 24 --bs-trees 1000
			#nohup
			nohup  sh concat_minus_nohup.sh > concat_minus_nohup.log &
		

#################
###09 June 2023###
#################		
		
	#Concatenate tree all exon minus five to FASTA file
		#concatenate the exons
		# catfasta2phyml for all exons/samples minus minus_five_bad_exons
		cd Trimal_minus_five_bad_exons
		catfasta2phyml -c /data/eaf236/Eperua/Trimal_minus_five_bad_exons/DN*.fasta > concatenated_all_minus_five_bad_exons.fasta 2> partitions_concatenated_all_minus_five_bad_exons_fasta.txt
	   	#worked; partition scheme is the same generated for the concatenated PHY file used in Partition Finder.
		
	
	
#################
###13 June 2023###
#################		
	
#Phyparts
		#Root bestTre and not the support tree
		
		# root trees
			mkdir rooted_trees_all_minus_five_bad_exons_bestTree

			cd Trees_raxml_minus
	
		for file in /data/eaf236/Eperua/Trees_raxml_minus/*.bestTree
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "GUIB16798"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		
		#no GUIB16798
		#ok! DN13548_c0_g1_i1_m45854_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN13548_c0_g1_i1_m45854_exon1.trimAL.mafft.fasta.raxml.bestTree                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN19054_c0_g1_i1_m2815_exon11.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19054_c0_g1_i1_m2815_exon11.trimAL.mafft.fasta.raxml.bestTree                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN19647_c0_g2_i1_m40560_exon6.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19647_c0_g2_i1_m40560_exon6.trimAL.mafft.fasta.raxml.bestTree                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "AUGOA4"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN19886_c0_g1_i1_m46583_exon3.
			for file in /data/eaf236/Eperua/Trees_raxml_minus/DN19886_c0_g1_i1_m46583_exon3.trimAL.mafft.fasta.raxml.bestTree                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN20944_c0_g1_i1_m40434_exon4.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN20944_c0_g1_i1_m40434_exon4.trimAL.mafft.fasta.raxml.bestTree                   
		do
			name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
			echo "rooting tree for $name"
			echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN22546_c0_g1_i1_m60856_exon3.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN22546_c0_g1_i1_m60856_exon3.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN25570_c0_g1_i2_m58578_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN25570_c0_g1_i2_m58578_exon1.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN26351_c0_g1_i1_m23568_exon3.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN26351_c0_g1_i1_m23568_exon3.trimAL.mafft.fasta.raxml.bestTree                   	
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN28031_c0_g1_i2_m59063_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN28031_c0_g1_i2_m59063_exon1.trimAL.mafft.fasta.raxml.bestTree                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		#Ok! DN29365_c1_g2_i3_m40788_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29365_c1_g2_i3_m40788_exon1.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN29365_c1_g2_i3_m40788_exon2.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29365_c1_g2_i3_m40788_exon2.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN29386_c0_g1_i1_m40845_exon5.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29386_c0_g1_i1_m40845_exon5.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN29585_c2_g1_i2_m39452_exon7.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29585_c2_g1_i2_m39452_exon7.trimAL.mafft.fasta.raxml.bestTree                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "AUGO1496"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN29746_c0_g1_i1_m43048_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN29746_c0_g1_i1_m43048_exon1.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN30459_c0_g1_i1_m16208_exon10.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN30459_c0_g1_i1_m16208_exon10.trimAL.mafft.fasta.raxml.bestTree                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN30871_c0_g2_i2_m7374_exon19.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN30871_c0_g2_i2_m7374_exon19.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "COPA1661"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN31690_c1_g10_i1_m18498_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31690_c1_g10_i1_m18498_exon1.trimAL.mafft.fasta.raxml.bestTree                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN31870_c0_g3_i1_m26604_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31870_c0_g3_i1_m26604_exon1.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN31922_c0_g1_i1_m1232_exon18.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN31922_c0_g1_i1_m1232_exon18.trimAL.mafft.fasta.raxml.bestTree                   
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN328_c0_g1_i1_m12305_exon1.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN328_c0_g1_i1_m12305_exon1.trimAL.mafft.fasta.raxml.bestTree                     
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		# Ok! DN34771_c0_g1_i1_m28451_exon10.
		for file in /data/eaf236/Eperua/Trees_raxml_minus/DN34771_c0_g1_i1_m28451_exon10.trimAL.mafft.fasta.raxml.bestTree                  
		do
		name=`basename $file .trimAL.mafft.fasta.raxml.bestTree`
		echo "rooting tree for $name"
		echo "HYME10482"| gotree reroot outgroup --input $file -l - > /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/$name.rooted.treefile.tre 
		done
		
		#didn`t work!
		cd installed_programs
		cd Phyparts-blackrim-37067bc3bc14
		#exclude old things in Phyparts directory rooted trees directory
		mkdir rooted_trees
		#copy rooted BetTree to directory rooted_trees
		cp /data/eaf236/Eperua/rooted_trees_all_minus_five_bad_exons_bestTree/*.tre /data/eaf236/installed_programs/Phyparts-blackrim-37067bc3bc14/rooted_trees/
		# rooted species tree in the directory Phyparts-blackrim-37067bc3bc14/
		#run phyparts
		cd Phyparts-blackrim-37067bc3bc14
		java -jar ../Phyparts-blackrim-37067bc3bc14/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees/ -m all_exons_minus_five_bad_exons_Astral_analysis_rooted.tre -o out_exons
		
#################
###14 June 2023###
#################	
		
#again with newick format
		cd Phyparts-blackrim-37067bc3bc14
		java -jar ../Phyparts-blackrim-37067bc3bc14/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees/ -m astral_rooted.nwk -o out_exons
		#nohup
		nohup  sh phyparts_nohup.sh > phyparts_nohup.log &
			#started at 12:23 06/12/23
			##error
			#Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
        #at java.base/java.util.HashMap$KeySet.iterator(HashMap.java:913)
        #at java.base/java.util.AbstractCollection.toArray(AbstractCollection.java:184)
        #at org.opentree.bitarray.LongBitSet$1.<init>(LongBitSet.java:219)
        #at org.opentree.bitarray.LongBitSet.iterator(LongBitSet.java:217)
        #at org.opentree.bitarray.CompactLongSet.iterator(CompactLongSet.java:368)
        #at org.opentree.bitarray.CompactLongSet.containsAll(CompactLongSet.java:270)
        #at phyparts.CladeMapper.mapConcordanceConflict(CladeMapper.java:156)
        #at phyparts.CmdLineInt.decon(CmdLineInt.java:291)
        #at phyparts.CmdLineInt.parse(CmdLineInt.java:79)
        #at phyparts.Main.main(Main.java:7)
        
        #Phyparts 
			#C:\Users\forte\OneDrive\Área de Trabalho\Phyparts-blackrim-37067bc3bc14\target
			java -jar phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees -m astral_rooted.nwk -o out_exons			
		
		#trying again server
		java -Xmx60g -jar ../Phyparts-blackrim-37067bc3bc14/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rooted_trees/ -m astral_rooted.nwk -o out_exons	
		#nohup
		nohup  sh phyparts_nohup.sh > phyparts_nohup.log &
		#started at 12:23 06/20/23
		
		
#################
###26 July 2023###
#################	
	
	#Species tree Astral
	#plot the quartet support
			mkdir Trees_astral_all_minus_five_bad_exons #(in directory Eperua)
			#copy all_exons_minus_five_bad_exons_bestTree.tre to Astral directory in my computer
			#running Astral 
			cd desktop/Atral 5.7.1
			java -jar astral.5.7.1.jar -i all_exons_minus_five_bad_exons_bestTree.tre -t2 -o all_exons_minus_five_bad_exons_Astral_quartet.tre 


			
#################
###16 August 2023###
#################		

#Seqkit

	
	seqkit grep -f selected_accessions.txt DN328_c0_g1_i1_m12305_exon1.FNA -o filtered_DN328_c0_g1_i1_m12305_exon1.FNA
	
	seqkit grep -f selected_accessions.txt DN997_c0_g1_i1_m32554_exon1.FNA -o filtered_DN997_c0_g1_i1_m32554_exon1.FNA
	
	seqkit grep -f selected_accessions.txt DN997_c0_g1_i1_m32554_exon2.FNA -o filtered_DN997_c0_g1_i1_m32554_exon2.FNA
	
	seqkit grep -f selected_accessions.txt DN997_c0_g1_i1_m32554_exon4.FNA -o filtered_DN997_c0_g1_i1_m32554_exon4.FNA
	
	seqkit grep -f selected_accessions.txt DN3533_c0_g1_i1_m33975_exon1.FNA -o filtered_DN3533_c0_g1_i1_m33975_exon1.FNA
	
	
#Maft

	#did manually additing the follow parameters:
			mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > $name.mafft.fasta
	
#Raxml 

	for file in \Users\Elenice\Desktop\standard-RAxML-master\WindowsExecutables_v8.2.4\*.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxmlHPC --all --msa /data/eaf236/Eperua/Trimal/$name --msa-format FASTA --model GTR+G --prefix $name --threads 6 --bs-trees 100;
        	done
        	
        	
        raxmlHPC --all --msa filtered_DN328_c0_g1_i1_m12305_exon1_mafft.fas --msa-format FASTA --model GTR+G --prefix filtered_DN328_c0_g1_i1_m12305_exon1_mafft --threads 6 --bs-trees 100;
        
#################
###16 August 2023###
#################
	
#Maft 

	for file in /home/elenicefortes/Documents/Test_five_exons/seqkit/filtered_*.FNA 
		do
			name=`basename $file .FNA`
			mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > /home/elenicefortes/Documents/Test_five_exons/maft/$name.mafft.fasta
		done
#Trimal 

		for file in /home/elenicefortes/Documents/Test_five_exons/maft/*.mafft.fasta
			do
				name=`basename $file .mafft.fasta`
				trimal -in $file -out /home/elenicefortes/Documents/Test_five_exons/trimal/$name.trimAL.mafft.fasta -automated1
			done
#Raxml

	for file in /home/elenicefortes/Documents/Test_five_exons/trimal/*.trimAL.mafft.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxml-ng --all --msa /home/elenicefortes/Documents/Test_five_exons/trimal/$name --msa-format FASTA --model GTR+G --prefix $name --bs-trees 100;
        	done
#Cat besttree
	cd trees_raxml
	cat *.bestTree > /home/elenicefortes/Documents/Test_five_exons/cat_bestTree/five_exons_besttrees.tre
	
###############################################################################################################################################################################	
	
####SegKit for all exons
	for file in /home/elenicefortes/Documents/Seqkit/retrieved_sequences_unaligned_minus_five/*.FNA
				do
					name=`basename $file .FNA`
					seqkit grep -f selected_accessions.txt /home/elenicefortes/Documents/Seqkit/retrieved_sequences_unaligned_minus_five/$name.FNA -o /home/elenicefortes/Documents/Seqkit/filtered_$name.FNA
	done

#Maft U

	for file in /home/elenicefortes/Documents/Seqkit/filtered_*.FNA 
		do
			name=`basename $file .FNA`
			mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > /home/elenicefortes/Documents/Maft/$name.mafft.fasta
		done
#Trimal 

		for file in /home/elenicefortes/Documents/Maft/*.mafft.fasta
			do
				name=`basename $file .mafft.fasta`
				trimal -in $file -out /home/elenicefortes/Documents/Trimal/$name.trimAL.mafft.fasta -automated1
			done
#Raxml ng 

	for file in /home/elenicefortes/Documents/Trimal/*.trimAL.mafft.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxml-ng --all --msa /home/elenicefortes/Documents/Trimal/$name --msa-format FASTA --model GTR+G --prefix $name --bs-trees 100;
        	done
#Cat besttree ####used to run SNaQ
	cd trees_raxml
	cat *.bestTree > /home/elenicefortes/Documents/Cat_BestTrees/filtered_883exons_besttrees.tre
	
#astral

	java -jar /home/elenicefortes/anaconda3/envs/Eperua_tree/share/astral-tree-5.7.8-0/astral.5.7.8.jar -i filtered_883exons_besttrees.tre -t2 -o filtered_883exons_besttrees_Astral.tre
	###run again with no quartets


	
		
#######################
###10 November 2023###
#######################	

# Redoing trees with 62 exons with no signal for recombination

	#First exclude manually the four taxa from the gene aligments
			#>GLBA499
			#>PRAE1000
			#>RUBI23804
			#>FROE32452
			#>TESS3457
			
	#Astral tree:
		#reconstruct gene trees #nobootstrap
		
		for file in /home/elenicefortes/Documents/TREE62exons/DN*.fasta
	do
        	name=`basename $file`
        	echo "Pruning fasta file to just screening set for $name"
        	raxml-ng --all --msa /home/elenicefortes/Documents/TREE62exons/$name --msa-format FASTA --model GTR+G --prefix $name --bs-trees 1;
        	done
        	
		#reconstruct species tree
		
		cat *.bestTree > /home/elenicefortes/Documents/TREE62exons/astral_tree/62exons_besttrees.tre
		
		#astral

	java -jar /home/elenicefortes/anaconda3/envs/Eperua_tree/share/astral-tree-5.7.8-0/astral.5.7.8.jar -i 62exons_besttrees.tre -t2 -o 62exons_Astral.tre
	
	java -jar /home/elenicefortes/anaconda3/envs/Eperua_tree/share/astral-tree-5.7.8-0/astral.5.7.8.jar -i 62exons_besttrees.tre -o 62exons_Astral2.tre
	
	#Concatenated tree
		#concat sequences
		conda install -c bioconda catfasta2phyml
		catfasta2phyml -c /home/elenicefortes/Documents/TREE62exons/DN*.fasta > concatenated_62exons.phy 2> partitions.txt
		
		#reconstruct concatenated tree ###not needed, similar result when I included praes3187
		conda activate Eperua_tree
		raxml-ng --all --msa /home/elenicefortes/Documents/TREE62exons/raxml_concat_tree/concatenated_62exons.phy --model GTR+G --prefix concatenated_exons --bs-trees 1000		
			
		#use this tree because it has 114 tips