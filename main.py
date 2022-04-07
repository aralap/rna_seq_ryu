import os, random
from snippet import get_config_parameter
from Bio import SeqIO
import pandas as pd

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

#General QC
#
def FASTQC(input_path):    
    fastqc_path = get_config_parameter('fastqc')
    output_dir = './QC_Plots' # future support for dynamic assign
    check_mkdir(output_dir)
    print(f'Running QC on {input_path}')
    os.system(f'{fastqc_path} {input_path} -o {output_dir}')


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

def TRIMO(input_path):
    # Trimomatic
    #This will perform the following:
    #Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
    #Remove leading low quality or N bases (below quality 3) (LEADING:3)
    #Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
    #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
    #Drop reads below the 36 bases long (MINLEN:36)
    basename = os.path.basename(input_path).split('.')[0]
    trimo  = '/Applications/Trimmomatic-0.39/trimmomatic-0.39.jar'
    adapter = '/Applications/Trimmomatic-0.39/adapters/TruSeq3-SE.fa' 
    output_dir = './TRIMMED_FASTQS' # future support for dynamic assign
    check_mkdir(output_dir)
    if not os.path.isfile(cwd+'/TRIMMED_FASTQS/'+basename+'.fq'):
        print('Trimmomatic')
        os.system( 'java -jar '+ trimo + ' SE -phred33 '+input_path+' '+ cwd+'/TRIMMED_FASTQS/'+basename+'.fq'+' ILLUMINACLIP:'+adapter+':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')
    else:
        print('Previous TRIMMED_FASTQ found, skipping trimming')
    return 0

def HISAT(input_path,reference):
        basename = os.path.basename(input_path).split('.')[0]
        print('Previous TRIMMED_FASTQ found, running hisat2')
        hisat_path = get_config_parameter('hisat')
        output_dir = './ALIGNED_READS' # future support for dynamic assign
        check_mkdir(output_dir)
        print(f'Running Hisat on {input_path}')
        print(f'{hisat_path} -p 8 --dta -x {reference} -U {input_path} -S {output_dir}/{basename}.sam')
        os.system(f'{hisat_path} -p 8 --dta -x {reference} -U {input_path} -S {output_dir}/{basename}.sam')
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

def STAR_INDEX():
    star = get_config_parameter('star')
    # Folder to save index
    d = '/Users/aaptekmann/Desktop/CDI/Reference_Genomes/SC'
    # Single file genome fasta (i.e. from ensemble refs)
    f = '/Users/aaptekmann/Desktop/CDI/Reference_Genomes/Glabrata/Candida_glabrata.GCA000002545v2.dna.toplevel.fa'
    # Gene annotations (could do without, but much much slower and worst quality result)
    gtf = '/Users/aaptekmann/Desktop/CDI/Reference_Genomes/Glabrata/Candida_glabrata.GCA000002545v2.52.gtf'
    l = f'{star} --runThreadN 16 --runMode genomeGenerate --genomeDir {d} --genomeFastaFiles {f} --sjdbGTFfile {gtf} --sjdbOverhang 149'
    os.system(l)

def STAR_ALIGN(gz_file, ref_genome, output_folder):
    # Index folder (previously generated)
    d = ref_genome#'/Users/aaptekmann/Desktop/CDI/Reference_Genomes/SC'
    gz = gz_file#'/Users/aaptekmann/Desktop/CDI/Erika_Schor/NC2-1_R1_001.fastq.gz'
    o =output_folder# '/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts/'
    # for some reason in mac, uncompress works funky, so, uncompress before run
    basename = os.path.basename(gz).split('.')[0]
    if os.path.exists(f'{o}/{basename}/ReadsPerGene.out.tab'):
        os.system(f'rm {o}/{basename}/Aligned.out.bam')
        os.system(f'rm {o}/{basename}/Aligned.sortedByCoord.out.bam')
        return 0
    check_mkdir(f'{o}/{basename}')
    os.system(f'gunzip -k {gz}')
    fq = gz[:-3]
    #l = f'STAR --runThreadN 16 --genomeDir {d} --readFilesIn {fq} --outFileNamePrefix {o}/{basename}/ '
    l = f'STAR --runThreadN 16 --outSAMtype BAM Unsorted SortedByCoordinate --quantMode GeneCounts --genomeDir {d} --readFilesIn {fq} --outFileNamePrefix {o}/{basename}/ '

    print(l)
    os.system(l)
    os.system(f'rm {o}/{basename}/Aligned.out.bam')
    os.system(f'rm {o}/{basename}/Aligned.sortedByCoord.out.bam')
    os.system(f'rm {fq}')

def SAMTOOLS(bname,o):
    l = f'samtools sort -o {o}/{bname}/sorted.bam {o}/{bname}/Aligned.out.sam'
    os.system(l)

def stringtie(bname,gtf,o):
    stringtie = get_config_parameter('stringtie')
    l = f'{stringtie} -o {o}/{bname}/counted.gtf -G {gtf} {o}/{bname}/sorted.bam'
    os.system(l)

def create_count_matrix(counts_dir): 
    import code
    directory = counts_dir #'/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts/'
    mtrx=0
    for filename in os.listdir(directory):
        pat = os.path.join(directory, filename)
        if os.path.exists(f'{pat}/ReadsPerGene.out.tab'):
            df = pd.read_csv(f'{pat}/ReadsPerGene.out.tab',names=["Gene", "CountU", "Count+", "Count-"], sep='\t', header=4)
            df[filename] = df.max(axis=1)
            df = df.drop(columns=['CountU', 'Count+', 'Count-'])
            if type(mtrx) == type(0):
                mtrx = df
            else:
                mtrx = mtrx.set_index('Gene').join(df.set_index('Gene'))    
            print(f'Yes {filename}')
    mtrx.to_csv(f'{directory}/count_matrix.tsv', sep='\t')        
    return mtrx

if __name__ == "__main__":
    # just sample list, get() method TBD
    cwd = os.getcwd()
    directory = '/Users/aaptekmann/Desktop/CDI/Erika_Schor/RNA_SEQ_DATA'

    i = 0 
    for filename in os.listdir(directory):
        print('\nProcesing:',filename)
        if filename.endswith('gz'):
            i+=1 
            #FASTQC(input_path)
            #TRIMO(input_path)
            reference = "/Users/aaptekmann/Desktop/CDI/Reference_Genomes/hg38_tran/genome_tran.1.ht2"
            gtf = '/Users/aaptekmann/Desktop/CDI/Reference_Genomes/Glabrata/Candida_glabrata.GCA000002545v2.52.gtf'

            o = '/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts'
            bname = os.path.basename(filename).split('.')[0]

            # HISAT(input_path, reference)
            # Get ref from https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz
            #STAR_INDEX()
            ref_genome = '/Users/aaptekmann/Desktop/CDI/Reference_Genomes/SC'
            gz_file = f'{directory}/{filename}'
            output_folder = '/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts/'
            STAR_ALIGN(gz_file,ref_genome,output_folder)
        #SAMTOOLS(bname,o)
        #stringtie(bname,gtf,o)
        # SRA(input_n)
        # SPAdes(input_n)
        # prokka(input_n)
        # translate_six_frame(input_n)
    counts_dir = '/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts/'
    create_count_matrix(counts_dir)
    
