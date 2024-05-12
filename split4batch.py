import sys
import os

# Because cd-batch accepts 1000 sequences max per file
# this function will spilt the input file into small files of 999 sequences max
# by default the output folder is named splits/ !! 

# python split4batch.py dp-proteome.faa output_folder/

def split4batch(input_file, output_folder = "./splits" ):

	if output_folder[-1] == "/":
		output_folder = output_folder[:-1]

	if not os.path.exists(output_folder):
	    os.makedirs(output_folder)

	with open(input_file) as f:

		file_num = 0
		seq_num = 0

		output_prefix = output_folder + "/part-"

		output = open(output_prefix + str(file_num) + ".faa", "w+")
		
		for line in f:

			if line[0] == ">":
				seq_num += 1
				if seq_num > 999:
					seq_num = 1
					file_num += 1
					output.close()
					output = open(output_prefix + str(file_num) + ".faa", "w+")
			
			print(line, end="", file = output)

		output.close()



if __name__ == "__main__":

	if len(sys.argv) < 2:
		print("input file missing")
		print("ex: python split4batch.py input.faa output_folder")
	elif len(sys.argv) == 2:
		split4batch(sys.argv[1])
	else: 
		split4batch(sys.argv[1], sys.argv[2])
