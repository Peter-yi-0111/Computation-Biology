import function
import matplotlib.pyplot as plt
import pdb

#---------------------------parameter change---------------------------

original_file = "salmonella-enterica/salmonella-enterica.reads.fna"
variant_file = "salmonella-enterica/salmonella-enterica-variant.reads.fna"
k = 20
error = 5

#---------------------------main--------------------------

##salmonella-enterica-dictionary
print("Create orginial file dictionary")
original_dictionary = function.genome_dictionary(original_file, k)

##plot salmonella-enterica k-mers graph
#print("Plot orginial file k-mers,k = {}".format(k))
#sort_dict = function.sort_dictionary_frequency(original_dictionary)
#function.plot_dictionary_frequency(sort_dict,k)

##variant salmonella-enterica-dictionary
print("Create variant file dictionary")
variant_dictionary = function.genome_dictionary(variant_file, k)

##plot salmonella-enterica k-mers graph
#print("Plot variant file k-mers,k = {}".format(k))
#sort_dict = function.sort_dictionary_frequency(variant_dictionary)
#function.plot_dictionary_frequency(sort_dict,k)

#delete sequence error
print("Delete sequence error")
original_dictionary = function.delete_error(original_dictionary, error)
variant_dictionary = function.delete_error(variant_dictionary, error)

#concate k_mers sequence
print("Concate k_mers sequence")
original_concat = function.concate_kmers(original_dictionary)
variant__concat = function.concate_kmers(variant_dictionary)

#calculate min distance between two sequence
print("Calculate min distance between two sequence")
mindistance = function.min_distances(original_concat, variant__concat)

#fix length
print("Fix length")
fix_length = function.fixed_length(mindistance)

#print SNP
print("\nResult:")
print("SNP")
function.print_colored_SNP(fix_length)

#---------------------------multi plot--------------------------
'''
#multi plot
multi_plot = True #Do you need graph merging for multiple k values? will take a lot of time！
k_init = 17
k_end = 77
k_interval = 5

while multi_plot == True :
    multi_plot = False
    plt.figure()
    for k in range(k_init,k_end,k_interval):
        print(k)
        orginal_dictionary = function.genome_dictionary(orginal_file, k)
        sort_dict = function.sort_dictionary_frequency(orginal_dictionary)

        function.multi_plot(sort_dict,k)
    plt.show()
'''







