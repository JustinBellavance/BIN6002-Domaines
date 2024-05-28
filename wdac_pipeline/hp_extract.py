import sys


# hypothetical proteins (hp) extraction
# python hp-extract dp-proteome.faa > hp-sub.faa
def extract_hypthetical_proteins(input_file, all_hp = False):

	with open(input_file) as f:

		print_line = False

		for line in f:
			if line[0] == ">":
				print_line = False
				if all_hp and "product=hypothetical protein" in line:
					print_line = True
				elif "product=hypothetical protein ;;" in line and ("motif:PFAM" in line or "domain" in line):
					print_line = True

			if print_line:
				print(line, end="")



if __name__ == "__main__":

	if len(sys.argv) < 2:
		print("input file missing")
		print("usage: python", sys.argv[0], "dp-proteome.faa > hp-sub.faa")
		print("usage: python", sys.argv[0], "dp-proteome.faa --all-hp > hp-all.faa")
	else: 
		extract_hypthetical_proteins(sys.argv[1], "--all-hp" in sys.argv)

