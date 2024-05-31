import sys
import time
import requests


CD_BATCH_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";


def split4batch(input_file):
	# splits large file to batches of 999 sequences max

	batches = list()
	current_batch = ""

	with open(input_file) as f:

		seq_num = 0
		
		for line in f:

			if line[0] == ">":
				seq_num += 1
				if seq_num > 999:
					seq_num = 1
					batches.append(current_batch)
					current_batch = ""
			
			current_batch += line
		
		batches.append(current_batch)

	return batches


def submit_to_cdbatch(sequences):


	PARAMS = {
	    'useid1': "true",
	    'maxhit': 250,
	    'filter': "true",
	    'db': 'pfam',
	    'evalue': 0.01,
	    'cddefl': "false",
	    'qdefl': "false",
	    'dmode': "rep",
	    'clonly': "false",
	    'tdata': "hits",
	    'queries': "".join(sequences)
	}

	req = requests.post(CD_BATCH_URL, data=PARAMS)

	if req.status_code == 200:
		res = req.text.split("\n")
		for row in res:
			if "#status" in row:
				if "3" in row:
					return res[1].split("\t")[1] # return the job id
				else:
					break

	raise ValueError("Something went wrong on job "+ job_id)


def check_job_status(job_id):

	req = requests.post(CD_BATCH_URL, data={'tdata': "hits", 'cdsid': job_id})

	if req.status_code == 200:
		for row in req.text.split("\n"):
			if "#status" in row:
				return "success" in row or "0" in row

	raise ValueError("Something went wrong on job "+ job_id)


def display_job_results(job_id):

	req = requests.post(CD_BATCH_URL, data={
		'tdata': "hits",
		'cddefl': "false",
		'qdefl': "false",
		'dmode': "rep",
		'clonly': "false",
		'cdsid': job_id
	})

	if req.status_code == 200:		
		print(req.text)
		return
	
	raise ValueError("Something went wrong on job "+ job_id)
	



if __name__ == "__main__":

	if len(sys.argv) < 2:
		print("usage: python", sys.argv[0], "input.faa")
		exit()

	batches = split4batch(sys.argv[1])

	for batch in batches:

		job_id = submit_to_cdbatch(batch)
			
		time.sleep(5) # To not overload the server :]

		while not check_job_status(job_id):
			time.sleep(5)

		display_job_results(job_id)
		