##A.4.1 Computation of $E_Y(\Phi_n)$ {#eyphi}


We shall compute the values of $E_Y(V_n)$, for $n=3,\ldots,8$, from the $V$ values of all trees in the corresponding $\mathcal{BT}_n$. Then, on the one hand, we have computed the probability of each tree in each $\mathcal{BT}_n$ under the Yule model with the following function: 
```{r}
yule.prob = function(tree){
  if (class(tree)=="phylo") 
    tree=graph.edgelist(tree$edge, directed=TRUE)  
  sp = shortest.paths(tree,mode = "out")
  deg = degree(tree,mode="out")
  leaves = which(deg==0)
  n = length(leaves) 
  k.node = function(node){
    subtree=which(sp[node,]<Inf)
    return(length(intersect(leaves,subtree)))
  } 
  kappas = sapply(which(deg>0), k.node) 
  value = (2^(n-1)/as.numeric(big.factorial(n)))*
              prod(1/(kappas-1))
  return(value)
} 
```

On the other hand, we have computed the $V$ value of each
tree with the function `V.index` and, finally, we have computed the desired expected value for each
n as the sum over all phylogenetic trees in $\mathcal{BT}_n$ of the product of their $V$-value and their probability:

```{r,eval=FALSE}
exp.V.yule = c()
for(n in 3:8){ 
  trees=read.tree(file=paste("./bintrees-n",n,".txt",sep=""))
  indices = sapply(trees, cophen.index)
  probs=sapply(trees, yule.prob)
  exp.yule[n]=sum(indices*probs)
} 
exp.V.yule
```
```{r,echo=FALSE}
read.table("./C2-EYPhi.txt")[,1]
```

So, the results agree with the figures given by our formula.



##A.4.2 Computation of $E_U(\Phi_n)$ {#euphi}
The formula in Theorem 2.28 can be computed with the following function (it
uses the function `big.double.factorial` explained in Section [$\color{blue}{\text{A.3}}$](#gefun)): 
```{r}
EUPhi = function(n){
  return(as.numeric((n*(n-1)/4)*
            (big.double.factorial(2*n-2)/
              big.double.factorial(2*n-3)-2)))
}
```

For $n=3,...,20$ the results are:
```{r}
sapply(3:20,EUPhi)
```

To double-check the formula, we have computed the values of $E_U(\Phi)$, for $n=3,\ldots,8$, as the arithmetic mean of the total cophenetic indices of all
phylogenetic trees in $\mathcal{BT}_n$ computed with our function `cophen.index`.
```{r,eval=FALSE}
exp.uni = c()
for(n in 3:8){ 
  trees=read.tree(file=paste("./bintrees-n",n,".txt",sep=""))
  indices = sapply(trees, cophen.index) 
  exp.uni[n]=mean(indices)
} 
exp.uni
```
```{r,echo=FALSE}
read.table("./C2-EUPhi.txt")[,1] 
```

Again, the results agree with the figures given by our formula.
