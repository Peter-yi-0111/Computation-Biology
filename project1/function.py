import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import sys
from Levenshtein import distance as lev
from Levenshtein import editops as edi
from termcolor import colored

def genome_dictionary(salmonella_enterica_file , k):

    dictionary = {}
    count = 0
    #read gene sequences in file
    for genome_sequence in Bio.SeqIO.parse(salmonella_enterica_file, "fasta"):

        #count = count + 1
        genome_sequence = str(genome_sequence.seq)

        for i in range(len(genome_sequence) - k + 1 ):
            dictionary_key = genome_sequence[i:i+k]

            if dictionary_key not in dictionary:
                dictionary[dictionary_key] = 1
            else:
                dictionary[dictionary_key] = dictionary[dictionary_key] + 1

        if count == 100:
            break

    return dictionary

def sort_dictionary_frequency(dictionary):
    dictionary_value = list(dictionary.values())
    dictionary_frequency = dict(Counter(dictionary_value))
    #sort
    sort_dict = {}
    for i in range(1,300):
        if i in dictionary_frequency.keys():
            sort_dict[i] = dictionary_frequency[i]
    return sort_dict

def plot_dictionary_frequency(sort_dict, k ):
    x = list(sort_dict.keys())
    y = list(sort_dict.values())
    plt.figure()
    plt.plot(x, y,  label = f'k = {k}')
    plt.yscale('log')
    plt.legend(loc = 'upper right')
    plt.show()

def multi_plot(sort_dict, k ):
    x = list(sort_dict.keys())
    y = list(sort_dict.values())
    plt.plot(x, y,  label = f'k = {k}')
    plt.yscale('log')
    plt.legend(loc = 'upper right')

def delete_error(dictionary, error):
    for key,value in dict(dictionary).items():
        if value < error:
            del dictionary[key]
    return dictionary

def concate_kmers(dictionary):
    orginal_list = list(dictionary)
    concate_list = []
    for i in range(len(orginal_list)):
        list_control = 0
        string = orginal_list[i]
        string_tail = string[1:]
        string_head = string[:-1]
        count = 1
        for j in range(i+1,len(orginal_list)):
            temp = orginal_list[j]
            temp_tail = temp[1:]
            temp_head = temp[:-1]
            if string_tail == temp_head:
                string = string + temp[-1]
                count = count + 1
                string_tail = string[count:]
                string_head = string[:-count]

            elif string_head == temp_tail:
                string = temp[0] + string
                count = count + 1
                string_tail = string[count:]
                string_head = string[:-count]

        for k in range(len(concate_list)):
            if string in concate_list[k]:
                list_control = list_control + 1
        if list_control == 0:
            concate_list.append(string)

    return concate_list

def min_distances(orginal_seqs, variant_seqs):
    #Levenshtein
    min_distances = {}
    for orginal in orginal_seqs:
        closest_variant = ''
        min_distance = sys.maxsize
        for variant in variant_seqs:

            distance = lev(orginal, variant)
            if distance <= min_distance:
                min_distance = distance
                closest_variant = variant

        if closest_variant != '':
            if min_distance <= 10:
                min_distances[(orginal,closest_variant)] = min_distance

    return min_distances

def fixed_length(min_distances):

    fix_length_min_dinstance = {}

    for key, val in min_distances.items():

        if len(key[0]) != len(key[1]):
            if len(key[0])>len(key[1]):
                step = edi(key[1],key[0])
                newstr = key[1]

            elif len(key[0])<len(key[1]):
                step = edi(key[0],key[1])
                newstr = key[0]

            for i in range(len(step)):
                if step[i][0] == 'insert':
                    newstr = newstr[:step[i][2]] + " " + newstr[step[i][2]:]
                elif step[i][0] == 'delete':
                    newstr = newstr[:step[i][2]]+newstr[step[i][2]+1:]

            if len(key[0])>len(key[1]):
                fix_length_min_dinstance[(key[0],newstr)] = val

            elif len(key[0])<len(key[1]):
                fix_length_min_dinstance[(newstr,key[1])] = val

        if len(key[0]) == len(key[1]):
            fix_length_min_dinstance[key] = val

    return fix_length_min_dinstance

def print_colored_SNP(min_distances):
    for key, val in min_distances.items():

        if len(key[0]) == len(key[1]):
            diff = [i for i in range(len(key[0])) if key[0][i] != key[1][i]]
            real = ''
            variant = ''

            for i in range(len(key[0])):
                if i in diff:
                    real += colored(key[0][i], 'red')
                    variant += colored(key[1][i], 'red')
                else:
                    real += key[0][i]
                    variant += key[1][i]

            print('SNP')
            print(f'\toriginal: {real}')
            print(f'\tvariant: {variant}')
            print("\n---------------------------------------------------------------------------------------------------\n")


