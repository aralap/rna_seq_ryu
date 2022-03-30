import os, random
from Bio import SeqIO

# DNA codon table 
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_dict = dict(zip(codons, amino_acids))

#
def check_mkdir(directory):
    if not os.path.exists(directory):
        mkdir = "mkdir "+directory
        os.system(mkdir)
    return 0

# SRA download
def SRA(input_n):
    sra_fold = '/Users/aaptekmann/sratoolkit.2.10.8-mac64/bin/'
    os.chdir(sra_fold)
    if not os.path.isfile(cwd+'/FASTQS/'+input_n+'.fastq' ):
        print('Trying to download fastq from NCBI')
        os.system(sra_fold+'fastq-dump '+input_n+ ' -O '+ cwd+'/FASTQS/')
    else:
        print('Previous FASTQ found, skipping download')
    return 0

def trimo(input_n):
    # Trimomatic
    #This will perform the following:
    #Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
    #Remove leading low quality or N bases (below quality 3) (LEADING:3)
    #Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
    #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
    #Drop reads below the 36 bases long (MINLEN:36)

    trimo  = '/Applications/Trimmomatic-0.39/trimmomatic-0.39.jar'
    adapter = '/Applications/Trimmomatic-0.39/adapters/TruSeq3-SE.fa' 
    input_f = cwd+'/FASTQS/'+input_n+'.fastq'

    if not os.path.isfile(cwd+'/TRIMMED_FASTQS/'+input_n+'.fq'):
        print('Trimmomatic')
        os.system( 'java -jar '+ trimo + ' SE -phred33 '+input_f+' '+ cwd+'/TRIMMED_FASTQS/'+input_n+'.fq'+' ILLUMINACLIP:'+adapter+':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')
    else:
        print('Previous TRIMMED_FASTQ found, skipping trimming')
    return 0

def SPAdes(input_n):   
    # metaSPAdes {Nurk, 2017 #74} with variable kmer sizes (-k 21,33,55,77,99,127).
    spades_folder = '/Users/aaptekmann/SPAdes/SPAdes-3.14.1-Darwin/bin/'
    if not os.path.isfile(cwd+'/ASSEMBLED/'+input_n+'/scaffolds.fasta'):
        os.system(spades_folder+'spades.py --only-assembler -s '+cwd+'/TRIMMED_FASTQS/'+input_n+'.fq -o '+cwd+'/ASSEMBLED/'+input_n  )
    else:
        print('Previos ASSEMBLY found, skipping spades')
        #print('MANUAL RUN COMMAND:\n'+spades_folder+'spades.py -s '+cwd+'/TRIMMED_FASTQS/'+input_n+'.fq -o '+cwd+'/ASSEMBLED/'+input_n )
    return 0

def prokka(input_n):
# Prokka {Seemann, 2014 #75} 
    if not os.path.isdir(cwd+'/ANNOTATED/'+input_n):
        os.system('prokka --centre C --locustag L --outdir '+cwd+'/ANNOTATED/'+input_n+' '+cwd+'/ASSEMBLED/'+input_n+'/scaffolds.fasta')
        os.system('cp '+cwd+'/ANNOTATED/'+input_n+'/*.faa '+ cwd+'/PROTEOMES/'+input_n+'.faa')
    else:
        print('Previos ANNOTATION found, skipping prokka')
    return 0


def translate_six_frame(input_n):
    #Input fastq file
    translated_records = [] 
    out_file = cwd+'/TRANS_READS/'+input_n+'.fasta'
    if os.path.isfile(out_file) and os.stat(out_file).st_size > 730:
        print('Previos TRANS_READS found, skipping TRANS_READS')
    else:  
        print('Running TRANS_READS.')  
        reads_file = cwd+'/TRIMMED_FASTQS/'+input_n+'.fq'
        records = [r for r in SeqIO.parse(reads_file, "fastq")]
        if len(records) > 1000000:
            records = random.sample(records, 2000000)
        for record in records:
            rid = record.id
            
            fwd0 = record.translate(to_stop=True)
            if len(fwd0) > 11 : 
                fwd0.id = rid+'_f0'
                translated_records.append(fwd0)

            fwd1 = record[1:].translate(to_stop=True)
            if len(fwd1) > 11 : 
                fwd1.id = rid+'_f1'
                translated_records.append(fwd1)

            fwd2 = record[2:].translate(to_stop=True)
            if len(fwd2) > 11 : 
                fwd2.id = rid+'_f2'
                translated_records.append(fwd2)
            
            rev =  record.reverse_complement()
            
            rev0 = rev.translate(to_stop=True)
            if len(rev0) > 11 : 
                rev0.id = rid+'_r0'
                translated_records.append(rev0)
            
            rev1 = rev[1:].translate(to_stop=True)
            if len(rev1) > 11 : 
                rev1.id = rid+'_r1'
                translated_records.append(rev1)

            rev2 = rev[2:].translate(to_stop=True)
            if len(rev2) > 11 : 
                rev2.id = rid+'_r2'
                translated_records.append(rev2) 
        if len(translated_records) > 2000000:
            translated_records = random.sample(translated_records, 2000000)           
        SeqIO.write(translated_records, cwd+'/TRANS_READS/'+input_n+'.fasta', "fasta")
    return 0

if __name__ == "__main__":
    # just sample list, get() method TBD
    input_list = ['SRR12347146','SRR12347145','SRR12347144','SRR12347143','SRR12347142','SRR12347141','SRR12347140']
    cwd = os.getcwd()


    i = 0 
    for input_n in input_list:
        i+=1 
        print('\nProcesing:',input_n, '\tSample', i, 'of', len(input_list))
        SRA(input_n)
        trimo(input_n)
        SPAdes(input_n)
        prokka(input_n)
        translate_six_frame(input_n)
