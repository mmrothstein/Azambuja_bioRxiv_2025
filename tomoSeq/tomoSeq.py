import os, subprocess
from collections import defaultdict

root = "/data/Megan/tomoAnalysis/tomoSeq"

filekey = {}

if not os.path.exists('%s/trimmed_fastq' % (root)):
	os.mkdir('%s/trimmed_fastq' % (root))

if not os.path.exists('%s/BAM' % (root)):
	os.mkdir('%s/BAM' % (root))

if not os.path.exists('%s/counts' % (root)):
	os.mkdir('%s/counts' % (root))

for fastq in os.listdir('%s/fastq/' % (root)):
	if fastq[:4] == '1073':
		sample = '_'.join(fastq.split('_')[4:8])
		filekey[sample]=fastq.rstrip()
	if fastq[:4] == '9898':
		sample = '_'.join(fastq.split('_')[4:6])
		filekey[sample]=fastq.rstrip()

for f in filekey.keys():
	I = '%s/fastq/%s' % (root, filekey[f])
	O = '%s/trimmed_fastq/trimmed_%s.fastq' % (root, f)
	cmd = "/home/mmr234/Programs/bbmap/bbduk.sh in=%s out=%s ref=/home/mmr234/Programs/bbmap/resources/polyA.fa.gz,/home/mmr234/Programs/bbmap/resources/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20" % (I, O)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out,err) = proc.communicate()
	print(out)
	print(err)

for f in filekey.keys():
	I = '%s/trimmed_fastq/trimmed_%s.fastq' % (root, f)
	O = '%s/BAM/%s_' % (root, f)
	cmd = "STAR --runThreadN 8 --genomeDir /data/Megan/tomoAnalysis/tomoSeq/genome/ --readFilesIn %s --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s" % (I, O)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out,err) = proc.communicate()
	print(out)
	print(err)

for f in filekey.keys():
	I = '%s/BAM/%s_Aligned.sortedByCoord.out.bam' % (root, f)
	cmd =  "samtools index %s" % (I)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out,err) = proc.communicate()
	print(out)
	print(err)

for f in filekey.keys():
	I = '%s/BAM/%s_Aligned.sortedByCoord.out.bam' % (root, f)
	cmd = "htseq-count -m intersection-nonempty -s yes -f bam -r pos %s /data/Megan/tomoAnalysis/tomoSeq/genome/Gallus_gallus.GRCg6a.99_SIRV.gtf > %s/counts/%s_counts.txt" % (I, root, f)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out,err) = proc.communicate()
	print(out)
	print(err)



counts = defaultdict(int)
genes = set()

for f in filekey.keys():
	for line in open('%s/counts/%s_counts.txt' % (root, f)):
		curline = line.rstrip().split('\t')
		counts[(f, curline[0])] = int(curline[1])
		genes.add(curline[0])


header = 'Gene'
for f in sorted(filekey.keys()):
	header = header + '\t' + f


with open("total_counts.txt","w") as out:
	out.write('%s\n' % (header))
	for gene in sorted(genes):
		toprint = gene
		for f in sorted(filekey.keys()):
			toprint = toprint + '\t' + str(counts[(f,gene)])
		out.write('%s\n' % (toprint))

















