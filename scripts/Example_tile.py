import sys
import os
from Bio import SeqIO
import glob
from collections import defaultdict
import numpy as np

# usage: python Example1.py <prefix of files R1 and R3, e.g.I7_i72_S5_L001> <output Folder (needs to be created before running the script)>
# inputs=sys.argv
# inputFile=inputs[1]
# outFolder=inputs[2]


for file in glob.glob("*R1_001.fastq"): #global function to search for every fastq file, python cannot do this by default 

    if str.split(file, "_")[0] != "Undetermined":  # if begins with undetermined by splitting the name according to _ and using the first element
        inputFile = str.split(file, "_R1")[0]
        outFolder = inputFile + ".Results"
        print inputFile
        os.makedirs(os.getcwd() + "/" + outFolder)

        good_tags = ["AAGTGCCG", "AATAATGT", "AGAGGACC", "AGTTAGAG", "CACGATGA", "CCGCGGGA", "GATGGCTC", "GCCCTCCG","GTACCGAG", "GTGGCTGC"]
        dict_tags = defaultdict(list)
        # each time we find a tag, we create a fastq file and append this read to the file, put the tag in used tags list, so if we face the same tag we will just find the file
        for record in SeqIO.parse((inputFile + "_R1_001.fastq"),"fastq"):  # each read is a record, which has an attribute
            Header = record.id  # id is the name
            Seq = record.seq  # seq is the sequence
            Tag = str(Seq[:8])
            coordinate = Header.split(":")[3]+Header.split(":")[4]+":"+Header.split(":")[5]+":"+Header.split(":")[6]  # these are the coordinates and tile number
            [coordinate].append(Tag)
            if Tag in good_tags:
                # print coordinate
                filename = outFolder + "/" + inputFile + "." + str(Tag) + "_R1.fastq"  # naming output file
                if os.path.exists(os.getcwd() + "/" + filename):  # if file already exists
                    output_handle = open(filename, "a")  # a means append the read to the output file
                else:
                    output_handle = open(filename, "w")  # otherwise create a new file
                output_handle.write(record.format("fastq"))  # write the read to the new file
                output_handle.close()
            else:
                filename_unas = outFolder + "/Unassigned" + "." + inputFile + ".R1.fastq"
                output_handle = open(filename_unas, "a")
                output_handle.write(record.format("fastq"))
                output_handle.close()
        print "R1 complete"

        for record2 in SeqIO.parse((inputFile +"_R3_001.fastq"),"fastq"):  # each read is a record, which has an attribute
            Header2 = record2.id  # id is the name
            coordinate2 = Header2.split(":")[3]+Header2.split(":")[4]+":" +Header2.split(":")[5] + ":" + Header2.split(":")[6]  # these are the coordinates
            Tag2 = dict_tags[coordinate2][0]  # Tag is assigned to the coordinates in the dictionary
            if Tag2 in good_tags:
                filename2 = outFolder + "/" + inputFile + "." + str(Tag2) + "_R2.fastq"
                if os.path.exists(os.getcwd() + "/" + filename2):  # if it exists
                    output_handle = open(filename2, "a")  # then append
                else:  # otherwise
                    output_handle = open(filename2, "w")  # make a new file
                output_handle.write(record2.format("fastq"))
                output_handle.close()
            else:
                filename_unas2 = outFolder + "/Unassigned" + "." + inputFile + ".R2.fastq"
                output_handle = open(filename_unas2, "a")
                output_handle.write(record2.format("fastq"))
                output_handle.close()
        print "R2 complete"




#for k in dict_tags.keys():
#	if len(dict_tags[k])>1:
#		print k