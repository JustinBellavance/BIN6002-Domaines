import os
import sys
import subprocess




def get_tm_score(dippa_pdb_path, ref_pdb_path):
	tmscore_process = subprocess.run(["./TMscore {} {}".format(dippa_pdb_path, ref_pdb_path)], shell=True, capture_output=True, text=True)

	if len(tmscore_process.stderr):
		#print("Error")
		#print(tmscore_process.stderr)
		
		return -1

	for line in tmscore_process.stdout.split("\n"):
		if "TM-score" == line[:8]:
			return (line.split()[2])


if __name__ == "__main__":

	if len(sys.argv) < 3:

		print("usage: python3 calc_tm_score.py DIPPA_XXXX.pdb ref_pdbs/")
		sys.exit()

	dippa_pdb = sys.argv[1]
	pdbs_folder = sys.argv[2]

	for pdb in os.listdir(pdbs_folder):
		tmscore = get_tm_score(dippa_pdb, pdbs_folder+pdb)
		print(dippa_pdb[:-4],pdb[:-4],tmscore)
