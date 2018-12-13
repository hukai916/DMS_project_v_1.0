AAList		= ['Gly', 'Ala', 'Ser', 'Thr', 'Cys', 'Pro',
	  			'Asn', 'Gln', 'Asp', 'Glu', 'His', 'Arg',
	  			'Val', 'Ile', 'Leu', 'Met', 'Phe', 'Tyr',
	  			'Trp', 'Ter', '???']

AA_GROUPS 	= [('Polar',0,7), ('Charged',8,12), ('Non-polar',13,19)]

Property2AA = dict([("Polar", ('Gly', 'Ala', 'Ser', 'Thr', 'Cys', 'Pro', 'Asn', 'Gln')),
		    		("Charged", ('Asp', 'Glu', 'His', 'Arg')),
		    		("Non-polar", ('Val', 'Ile', 'Leu', 'Met', 'Phe', 'Tyr', 'Trp'))])

CodonList 	= ["TTT","TTC","TTA","TTG","CTT","CTC","CTA",
				"CTG","ATT","ATC","ATA","ATG","GTT","GTC",
				"GTA","GTG","TCT","TCC","TCA","TCG","AGT",
				"AGC","CCT","CCC","CCA","CCG","ACT","ACC",
				"ACA","ACG","GCT","GCC","GCA","GCG","TAT",
				"TAC","TAA","TAG","TGA","CAT","CAC","CAA",
				"CAG","AAT","AAC","AAA","AAG","GAT","GAC",
				"GAA","GAG","TGT","TGC","TGG","CGT","CGC",
				"CGA","CGG","AGA","AGG","GGT","GGC","GGA",
				"GGG"]

CODON_GROUPS = [("Phe",0,1),("Leu",2,7),("Ile",8,10),
				("Met",11,11),("Val",12,15),("Ser",16,21),
				("Pro",22,25),("Thr",26,29),("Ala",30,33),
				("Tyr",34,35),("Stop",36,38),("His",39,40),
                ("Gln",41,42),("Asn",43,44),("Lys",45,46),
				("Asp",47,48),("Glu",49,50),("Cys",51,52),
				("Trp",53,53),("Arg",54,59),("Gly",60,63)]

AA2CODON 	= dict([("Phe", ("TTT","TTC")),
					("Leu", ("TTA","TTG","CTT","CTC","CTA","CTG")),
                	("Ile", ("ATT","ATC","ATA")),
					("Met", ("ATG")),
                	("Val", ("GTT","GTC","GTA","GTG")),
					("Ser", ("TCT","TCC","TCA","TCG","AGT","AGC")),
                	("Pro", ("CCT","CCC","CCA","CCG")),
					("Thr", ("ACT","ACC","ACA","ACG")),
                	("Ala", ("GCT","GCC","GCA","GCG")),
					("Tyr", ("TAT","TAC")),
                	("Stop", ("TAA","TAG","TGA")),
					("His", ("CAT","CAC")),
                	("Gln", ("CAA","CAG")),
					("Asn", ("AAT","AAC")),
                	("Lys", ("AAA","AAG")),
					("Asp", ("GAT","GAC")),
                	("Glu", ("GAA","GAG")),
					("Cys", ("TGT","TGC")),
                	("Trp", ("TGG")),
					("Arg", ("CGT","CGC","CGA","CGG","AGA","AGG")),
                	("Gly", ("GGT","GGC","GGA","GGG")),
					("Ter", ("TAA","TAG","TGA"))])

CODON2AA 	= dict([("TTT",("Phe","F")),("TTC",("Phe","F")),("TTA",("Leu","L")),
                 	("TTG",("Leu","L")),("CTT",("Leu","L")),("CTA",("Leu","L")),
                 	("CTG",("Leu","L")),("CTC",("Leu","L")),("ATT",("Ile","I")),
                 	("ATC",("Ile","I")),("ATA",("Ile","I")),("ATG",("Met","M")),
                 	("GTT",("Val","V")),("GTC",("Val","V")),("GTA",("Val","V")),
	                ("GTG",("Val","V")),("TCT",("Ser","S")),("TCC",("Ser","S")),
	                ("TCA",("Ser","S")),("TCG",("Ser","S")),("CCT",("Pro","P")),
	                ("CCC",("Pro","P")),("CCA",("Pro","P")),("CCG",("Pro","P")),
	                ("ACT",("Thr","T")),("ACC",("Thr","T")),("ACA",("Thr","T")),
	                ("ACG",("Thr","T")),("GCT",("Ala","A")),("GCC",("Ala","A")),
	                ("GCA",("Ala","A")),("GCG",("Ala","A")),("TAT",("Tyr","Y")),
	                ("TAC",("Tyr","Y")),("TAA",("Stop","*")),("TAG",("Stop","*")),
	                ("TGA",("Stop","*")),("CAT",("His","H")),("CAC",("His","H")),
	                ("CAA",("Gln","Q")),("CAG",("Gln","Q")),("AAT",("Asn","N")),
	                ("AAC",("Asn","N")),("AAA",("Lys","K")),("AAG",("Lys","K")),
	                ("GAT",("Asp","D")),("GAC",("Asp","D")),("GAA",("Glu","E")),
	                ("GAG",("GlT","E")),("TGT",("Cys","C")),("TGC",("Cys","C")),
	                ("TGG",("Trp","W")),("CGT",("Arg","R")),("CGC",("Arg","R")),
	                ("CGA",("Arg","R")),("CGG",("Arg","R")),("AGT",("Ser","S")),
	                ("AGC",("Ser","S")),("AGA",("Arg","R")),("AGG",("Arg","R")),
	                ("GGT",("Gly","G")),("GGC",("Gly","G")),("GGA",("Gly","G")),
	                ("GGG",("Gly","G"))])

CODON321 	= dict([("Ala", "A"), ("Gln", "Q"), ("Lys", "K"), ("Trp", "W"),
                 	("Arg", "R"), ("Gly", "G"), ("Met", "M"), ("Tyr", "Y"),
                 	("Asn", "N"), ("His", "H"), ("Phe", "F"), ("Val", "V"),
                 	("Cys", "C"), ("Ile", "I"), ("Ser", "S"), ("Pro", "P"),
                 	("Glu", "E"), ("Leu", "L"), ("Thr", "T"), ("Asp", "D"),
                 	("Ter", "*"), ("Stop", "*")])
