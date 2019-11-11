import os
import glob

energies = [30,50,70,90,150,250]

path = "/home/software/Calo/results/jet_highstat/jetscan30k/"

for e in energies:
	string = "jetscan_"+str(e)
	f = "jetscan_"+str(e)+".root"
	cmnd1 = "time ./dohist"+" "+str(string)+" "+str(path+string)
	os.system(cmnd1)	
