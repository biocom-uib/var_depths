# var_depths
Variance of leaves' depths related algorithms. Scripts associated to the paper "On the variance of the leaves' depths
in a rooted tree I: Extremal values".

## Contents

### R folder

It contains the following *R Mardown* files:

* `Figure11.Rmd`: Produces Figure 11 in the main text
* `Figure12.Rmd`: Produces Figure 11 in the main text
* `MinVar.Rmd`: Illustrates the function `MinVar` to compute the minimum value of $V$ for trees with $n$ leaves and the multiset(s) of depths of the tree(s) achieving it.
* `Lists_of_MinVar.Rmd`: Computes the tables in folder *Tables*
* `UniqueMinima.Rmd`: Checks the uniqueness of the type of trees achieving the minimum $V$ value for $n$ up to 2^20^

It also contains the HTML files obtained knitting them.

### Python folder
Python script defining the functions
* `var_depths`: that computes the variance of the leaves' depths of a given tree
* `max_var_depths`: that returns the tuple formed by the maximum variance of a given number of leaves and the only tree that attains it
* `min_var_depths`: that returns the tuple formed by the minimum variance of a given number of leaves and ONE tree that attains it
* `min_var_depths_vector`: that returns the tuple formed by the minimal variance of a given number of leaves and the vector of l_i's  

### Tables folder

It contains the types of trees with minimum V value for n from 2^8^ to 2^20^. Each file *Resultsk_k+1.txt* contains the results from 2^k^to 2^k+1^. The last file, from 2^19^ to 2^20^, is split into two files.

### Table 3 folder

It contains the R Markdown file and the files used in the computations yielding Table 3 in the main text
