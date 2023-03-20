Biology Project 1

Please download the package first

Download the package:
$ pip install -r requirements.txt

Execute program:
$ python computational_biology.py

All result images are stored in "result"！

##Notice
original_file : Put the fna file of salmonella-enterica
variant_file : Put the fna file of salmonella-enterica-variant

If you want to display a graph of Kmers, uncomment the following programs:

    print("Plot orginial file k-mers,k = {}".format(k))
    sort_dict = function.sort_dictionary_frequency(original_dictionary)
    function.plot_dictionary_frequency(sort_dict,k)

    print("Plot variant file k-mers,k = {}".format(k))
    sort_dict = function.sort_dictionary_frequency(variant_dictionary)
    function.plot_dictionary_frequency(sort_dict,k)

If you want to display Kmers plots for multiple K values, uncomment the following program:

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
