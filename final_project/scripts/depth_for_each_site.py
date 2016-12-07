### get full site

finaler = open("scripts/temp/patient1.final.bed","r");
full = open("scripts/temp/patient1.final.full.txt","w");

site = [];
dsite = [];
chrom = "chr1";
esite = 0;
header = 1;

for line in finaler:
	if header == 1:
		site.append(int(line.split("\t")[1]));
		esite = int(line.split("\t")[2]);
		dsite.append(int(line.split("\t")[3]));
		chrom = str(line.split("\t")[0]);
		header += 1;
		continue;
	if chrom == str(line.split("\t")[0]) and int(line.split("\t")[1]) == esite :
		site.append(int(line.split("\t")[1]));
		esite = int(line.split("\t")[2]);
		dsite.append(int(line.split("\t")[3]));
	else:
		#print
		k = 0;
		for j in range(site[0],site[-1]):
			if(int(j) < site[k+1]):
				full.write(chrom + "\t" + str(j) + "\t" + str(dsite[k]) + "\n");
			else:
				k += 1;
				full.write(chrom + "\t" + str(j) + "\t" + str(dsite[k]) + "\n");
		#reinitial
		chrom = str(line.split("\t")[0]);
		site = [int(line.split("\t")[1])];
		dsite = [int(line.split("\t")[3])];
		esite = int(line.split("\t")[2]);

#print the last one.
k = 0;
for j in range(site[0],site[-1]):
	if(int(j) < site[k+1]):
		full.write(chrom + "\t" + str(j) + "\t" + str(dsite[k]) + "\n");
	else:
		k += 1;
		full.write(chrom + "\t" + str(j) + "\t" + str(dsite[k]) + "\n");

	

finaler.close();
full.close();

exit();
