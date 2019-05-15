# -*- coding: future_fstrings -*-

import cProfile
import decimal
def num_lines_in_file(filename):
    '''
    Finds the number of lines in the file provided

    filename = filename, example: filename.txt

    how to call: num_of_lines = num_lines_in_file('filename.txt')
    '''
    num_lines = sum(1 for line in open(filename))
    return num_lines

def generate_sequence_file(filename,filewritename):
    seq_file = open(filewritename,'w')
    with open(filename, 'r') as f:
        for count, line in enumerate(f, start = 1):
            if (count+2)%4 == 0:  #DNA sequence is in every 4 lines (2nd one in the list)
                seq_file.write(line)

def generate_quality_file(filename,filewritename):
    seq_file = open(filewritename,'w')
    with open(filename, 'r') as f:
        for count, line in enumerate(f, start = 1):
            if count%4 == 0:  #Quality sequence is in every 4 lines (4th one in the list)
                seq_file.write(line)

def quality_quantification(filename, filewritename,num_values_quantized):
    """
    num_values_quantized = number of quality values you want to take in account because quality gets worse as the sequence progresses, look at sequence file to coordinate this
    """

    quality_file = list(open(filename,'r')) #converted to list for efficiency in future computations

    # dictionary denoted by curly braces
    fastq_quality = {}
    qualities = '\n!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~' #from fastq documentation arranged from lowest quality to highest quality

    #
    for index in range(len(qualities)): #iterates through each char in qualities
        fastq_quality[qualities[index]] = index #assigns values for each char in qualities identical to their index

    quality_quant = open(filewritename,'w') #creates the file that will receive the written quality number totals per sequence

    for i in range(len(quality_file)): #iterates through the each element of the quality_file list
        qual_sum = 0
        for j in range(1,num_values_quantized+2): #ignore first character of each line
            qual = fastq_quality[quality_file[i][j]] #dictionary value of char j in the quality line i of quality_file
            qual_sum += qual #sums up qualities for each line
        quality_quant.write(str(qual_sum)+'\n')

def remove_bad_quality_seq(sequence_only_file, quality_quantified_file,filewritename):
    import statistics
    seq_file_modified = open(filewritename,'w')
    seq_file = list(open(sequence_only_file,'r'))
    quality_quant_file = list(open(quality_quantified_file,'r'))
    quality_quant_file = [int(i) for i in quality_quant_file] #quality_quant_file = list(map(int,list(open(quality_quantified_file,'r')))) also works but a little slower
    median_quality = statistics.median(quality_quant_file)
    std_quality = statistics.stdev(quality_quant_file)

    for index in range(len(quality_quant_file)):
        if quality_quant_file[index] < (median_quality - (0.25*std_quality)): #Formula created from finding the least background and keeping that as best result
            pass
        else:
            seq_file_modified.write(seq_file[index])



def percent_indels(filename):
    #barcode should always be the 2-6 bases in the sequence
    #potential inputs
    #decimal.setcontext(decimal.Context(prec=15))

    quality_ensuring_sequence = 'TCTCTCATCGTGGGGG' #should be present with the barcode if the sequencing for the specific sequence was successful
    start_basenum_of_qual_ens_seq = 40
    end_basenum_of_qual_ens_seq = 56
    CRISPR_cut_site_6_bases = "ACCTAC"
    location_wild_type_non_frame_shift = 70 #CATCTT is located at 70
    wild_type_non_frame_shift = 'CATCTT'

    f = open(filename, 'r')
    sequence_list = f.readlines()

    reg = [0]*15
    mut = [0]*15
    frame_shift = [0]*15

    #dictionary
    mouse_barcodes = {"TCACG":0,"TAGGC":1,"CAGTG":2,"AGATC":3,"ATCAG":4,"GCTAC":5,"GTCAA":6,"TGTCA":7,"TAGAG":8,"TGAAA":9,"TTTCG":10,"AGTGG":11,"CTGAT":12,"TTCCT":13,"AACTA":14}
    for sequence in sequence_list:
        for barcode in mouse_barcodes:
            i = mouse_barcodes.get(barcode,0)
            if (sequence[1:6] == barcode) and (sequence[start_basenum_of_qual_ens_seq:end_basenum_of_qual_ens_seq] == quality_ensuring_sequence):
                if CRISPR_cut_site_6_bases in sequence:
                    reg[i]+=1
                else:
                    mut[i]+=1
                    for start in range(location_wild_type_non_frame_shift - (location_wild_type_non_frame_shift//3)*3,location_wild_type_non_frame_shift+34,3):
                        if (sequence[start:start+6] == wild_type_non_frame_shift): #location of CATCTT is 70 and location of ACCTAC is 59 so to determine whether it's frameshift, I have to scan 3 base multiples of where CATCTT is
                            pass
                        else:
                            frame_shift[i]+=1
                            break

    mouse_names = ['RP1','RP10','RP11','RPNM','RQ1','RQ8','RQ10','RQ11','RQNM','SR2','SR10','SR11','SY1','SY2','SY10']
    indel_percents = {} #dict
    frame_shift_percents = {}

    for index in range(len(mouse_names)):
        percent_indels = (float(mut[index])/(mut[index]+reg[index]))*100
        percent_frame_shifts_if_mutants = (float(frame_shift[index])/mut[index])*100
        indel_percents[mouse_names[index]] = percent_indels
        frame_shift_percents[mouse_names[index]] = percent_frame_shifts_if_mutants



    for mouse in indel_percents:
        print(f"Percent Indels in {mouse}: {indel_percents[mouse]}%")
        # print(f"Percent Frameshift Mutations if Mutant in {mouse}: {frame_shift_percents[mouse]}%")


def main():

    print('Number of lines in file:', num_lines_in_file('RPRQ_1.txt'))
    generate_sequence_file('RPRQ_1.txt', 'RPRQ_1_seq_only.txt')
    generate_quality_file('RPRQ_1.txt', 'RPRQ_1_quality_only.txt')
    quality_quantification('RPRQ_1_quality_only.txt','RPRQ_1_quality_quantification_first_86.txt',86)
    remove_bad_quality_seq('RPRQ_1_seq_only.txt','RPRQ_1_quality_quantification_first_86.txt','RPRQ_1_seq_only_modified_for_quality_first_86.txt')
    percent_indels('RPRQ_1_seq_only_modified_for_quality_first_86.txt')

if __name__ == "__main__":
    cProfile.run('main()')
