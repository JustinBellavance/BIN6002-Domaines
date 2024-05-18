import sys



# hypothetical proteins (hp) extraction
# python hp-extract dp-proteome.faa > hp-sub.faa
def extract_hypothetical_proteins_with_domains(input_file):

	with open(input_file) as f:

		print_line = False

		for line in f:
			if line[0] == ">":
				print_line = False
				#if "product=hypothetical protein" in line:
				#if "product=hypothetical protein" in line and ";;" not in line:
				if "product=hypothetical protein ;;" in line and ("motif:PFAM" in line or "domain" in line):
					print_line = True

			if print_line:
				print(line, end="")

def extract_hypothetical_proteins(input_file):

	with open(input_file) as f:

		print_line = False

		for line in f:
			if line[0] == ">":
				print_line = False
				if "product=hypothetical protein" in line:
					print_line = True

			if print_line:
				print(line, end="")

def extract_hypothetical_proteins_no_domains(input_file):

	with open(input_file) as f:

		print_line = False

		for line in f:
			if line[0] == ">":
				print_line = False
				if "product=hypothetical protein" in line.lower() and (";;" not in line and "domain" not in line):
					print_line = True

			if print_line:
				print(line, end="")


if __name__ == "__main__":

	if len(sys.argv) < 2:
		print("input file missing")
		print("ex: python hp-extract dp-proteome.faa > hp-sub.faa")
	elif len(sys.argv) == 3 and sys.argv[2] == "--all":
		extract_hypothetical_proteins(sys.argv[1])
	elif len(sys.argv) == 3 and sys.argv[2] == "--empty":
		extract_hypothetical_proteins_no_domains(sys.argv[1])
	else: 
		extract_hypothetical_proteins_with_domains(sys.argv[1])

