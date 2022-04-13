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
    d = get_config_parameter('ref_genome_folder')
    # Single file genome fasta (i.e. from ensemble refs)
    f = get_config_parameter('genome_fasta')
    # Gene annotations (could do without, but much much slower and worst quality result)
    gtf = get_config_parameter('genome_gtf')
    l = f'{star} --runThreadN 16 --runMode genomeGenerate --genomeDir {d} --genomeFastaFiles {f} --sjdbGTFfile {gtf} --sjdbOverhang 149 '
    os.system(l)

def STAR_ALIGN(gz_file, ref_genome, output_folder):
    # Index folder (previously generated)
    d = ref_genome#'/Users/aaptekmann/Desktop/CDI/Reference_Genomes/SC'
    gz = gz_file#'/Users/aaptekmann/Desktop/CDI/Erika_Schor/NC2-1_R1_001.fastq.gz'
    o =output_folder# '/Users/aaptekmann/Desktop/CDI/Erika_Schor/Counts/'
    # for some reason in mac, uncompress works funky, so, uncompress before run
    basename = os.path.basename(gz).split('.')[0]
    if os.path.exists(f'{o}/{basename}/ReadsPerGene.out.tab'):
        if os.path.exists(f'{o}/{basename}/Aligned.out.bam'):
            os.system(f'rm {o}/{basename}/Aligned.out.bam')
        if os.path.exists(f'{o}/{basename}/Aligned.sortedByCoord.out.bam'):    
            os.system(f'rm {o}/{basename}/Aligned.sortedByCoord.out.bam')
        return 0
    check_mkdir(f'{o}/{basename}')
    os.system(f'gunzip -k {gz}')
    fq = gz[:-3]
    #l = f'STAR --runThreadN 16 --genomeDir {d} --readFilesIn {fq} --outFileNamePrefix {o}/{basename}/ '
    l = f'STAR --runThreadN 16 --outSAMtype BAM Unsorted SortedByCoordinate --quantMode GeneCounts --genomeDir {d} --readFilesIn {fq} --outFileNamePrefix {o}/{basename}/ '

    print(l)
    os.system(l)
    # Add checks for existance...
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
            print(f'Reading gene counts for {filename}')
            df = pd.read_csv(f'{pat}/ReadsPerGene.out.tab',names=["Gene", "CountU", "Count+", "Count-"], sep='\t', header=4)
            df[filename] = df.max(axis=1)
            df = df.drop(columns=['CountU', 'Count+', 'Count-'])
            if type(mtrx) == type(0):
                mtrx = df
                print('Initialized matrix')
            else:
                mtrx[filename] = df[filename]   
    mtrx.to_csv(f'{directory}/count_matrix.tsv', sep='\t',  index=False)        
    return mtrx

def create_exp_design(counts_dir):
    for line in open(f'{counts_dir}/count_matrix.tsv','r'):
        names = line
        break
    out = open(f'{counts_dir}exp_design.tsv','w')
    for name in names.split():
        if 'Gene' in name:
            out.write(f'run\tcondition\ttype\ttime\n')
            continue
        elif 'H' in name:
            NA1='treatment'
            NA2="oxidative"
        elif 'M' in name: 
            NA1='treatment'
            NA2="dna_damage"
        elif 'NN' in name: 
            NA1='treatment'
            NA2="N_starvation" 
        elif 'P3' in name: 
            NA1='treatment'
            NA2="acidic"   
        elif 'Na' in name: 
            NA1='treatment'
            NA2="osmotic" 
        elif 'CR' in name: 
            NA1='treatment'
            NA2="congo_red"
        elif 'NC' in name: 
            NA1='treatment'
            NA2="C_starvation" 
        elif 'CT' in name: 
            NA1='control'
            NA2="time_series"
        elif 'NS' in name: 
            NA1='control'
            NA2="no_stress" 
        else:
            'ERROR'       
        if '-4' in name:
            NA3 = 4
        elif '-2' in name:
            NA3 = 2 
        elif '-1' in name:
            NA3 = 1    
        else:
            'ERROR'                                        
        out.write(f'{name}\t{NA1}\t{NA2}\t{NA3}\n')
    out.close()

if __name__ == "__main__":
    STAR_INDEX()
    # # just sample list, get() method TBD
    # cwd = os.getcwd()
    # directory = get_config_parameter('data_dir')
    # reference = get_config_parameter('genome_fasta')
    # ref_genome_folder = get_config_parameter('ref_genome_folder')

    # i = 0 
    # for filename in os.listdir(directory):
    #     if filename.endswith('gz'):
    #         print('\nProcesing:',filename)
    #         i+=1 
    #         #FASTQC(input_path)
    #         #TRIMO(input_path)

    #         o = get_config_parameter('out_dir')
    #         bname = os.path.basename(filename).split('.')[0]

    #         # HISAT(input_path, reference)
    #         gz_file = f'{directory}/{filename}'
    #         STAR_ALIGN(gz_file,ref_genome_folder,o)
    #     #SAMTOOLS(bname,o)
    #     #stringtie(bname,gtf,o)
    #     # SRA(input_n)
    #     # SPAdes(input_n)
    #     # prokka(input_n)
    #     # translate_six_frame(input_n)
    # create_count_matrix(o)
    # create_exp_design(o)
