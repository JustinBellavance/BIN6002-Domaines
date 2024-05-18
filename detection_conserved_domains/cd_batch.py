# Inspired from the perl script from the ncbi cd-batch help page
# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBWebAPI_samplePERL

import requests
import sys
import time


CD_BATCH_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";


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

	print("sopmething went wrong")
	exit()


def check_job_status(job_id):

	req = requests.post(CD_BATCH_URL, data={'tdata': "hits", 'cdsid': job_id})

	if req.status_code == 200:
		for row in req.text.split("\n"):
			if "#status" in row:
				return "success" in row or "0" in row

	print("something went wrong", job_id)
	exit()


def display_results(job_id):

	req = requests.post(CD_BATCH_URL, data={
		'tdata': "hits",
		'cddefl': "false",
		'qdefl': "false",
		'dmode': "rep",
		'clonly': "false",
		'cdsid': job_id
	})

	if req.status_code != 200:		
		print("something went wrong", job_id)
		exit()

	print(req.text)


if __name__ == '__main__':
	
	if len(sys.argv) < 2:
		print("usage: python cd-batch.py input.faa")
		print("note: the input file should be less than 1000 sequences.")
		exit()

	sequences = open(sys.argv[1]).readlines()

	job_id = submit_to_cdbatch(sequences)

	time.sleep(10)
	
	while not check_job_status(job_id):
		time.sleep(10)

	display_results(job_id)