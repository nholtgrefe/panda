# `phypanda`
`phypanda` is a Python package that contains PaNDA (Phylogenetic Network Diversity Algorithms): a software framework for exploring, visualizing and maximizing phylogenetic diversity in phylogenetic networks.


## List of important features
- Maximizing all-paths diversity in a rooted phylogenetic network
- Computing the all-paths diversity for a given set of taxa in a rooted phylogenetic network
- Greedily select a set of taxa with large all-paths diversity

## Installation
If you have an up-to-date version of [Python](https://www.python.org/downloads/) installed on your device, the standard package manager `pip` should come pre-installed. Then, you can install `phypanda` from [PyPI](https://pypi.org/project/phypanda/) by simply using the following command in a terminal:

```
python -m pip install phypanda
```

  
## Example usage

### Importing the package
To get started with `phypanda`, open a Python shell and import the package with:

```
import phypanda as ppa
```

### Maximizing all-paths diversity
To maximize the all-paths diversity of a rooted network (given in `eNewick` format) for a set of `k` taxa, use the function `solve_MAPPD`. For example, when finding a set of 10 taxa with maximum all-paths diversity in the [Xiphophorus network](https://github.com/nholtgrefe/panda/blob/main/data/exp2_xiphophorus_network.txt) from the PaNDA paper, run:

```
enewick = '((((((((((Xgordoni:1.3295084631587457,Xmeyeri:1.3295084631587457):0.0,Xcouchianus:1.329508093234352):6.999730834529853,Xvariatus:8.329238927764205):2.1769451514229345,Xevelynae:10.50618407918714):1.118605313770228,(Xxiphidium:7.2210504457107145,#H24:0.0):4.403738947246653):0.0,Xmilleri:11.624787067955268):4.296868586395352,Xandersi:15.92165565435062):0.9486610416497712,Xmaculatus:16.87031669600039):0.5723386247384958,((((Xmontezumae:7.221055986870681,(Xcortezi:5.485599585171238,((Xmalinche:5.485605240002155,Xbirchmanni:5.485605240002155):0.0)#H26:0.0):1.7354564016994427):0.0,((Xnigrensis:2.4303498026154564,Xmultilineatus:2.4303498026154564):0.19174715477323678,(Xpygmaeus:1.347820846400494,Xcontinens:1.347820846400494):1.2742761109881993):4.598960284156991):0.0,#H26:1.7354540549075645):0.0)#H24:10.2216024589192):2.1886232296055894,((Xclemenciae:11.254572014210282,Xmonticolus:11.254572014210282):6.4012991117391564,(#H25:1.6332001759073602,(Xsignum:10.266850863153604,((Xhellerii:8.633649976742058)#H25:1.6332013685506936,(Xalvarezi:8.362082334652573,Xmayae:8.362082334652573):1.9047690106401785):0.0):0.0):7.3890209733000205):1.975407424395037);'
k = 10
pd, taxa = ppa.solve_MAPPD(enewick, k)
```

To print the resulting maximum all-paths diversity and the selected taxa, run:

```
print(f"Maximum all-paths diversity for k = {k} is {pd}")
print("Selected taxa:", taxa)
```

        
For a complete overview of different methods and extra parameter options, please check the method descriptions in the [source code](https://github.com/nholtgrefe/panda/tree/main/phypanda/src/phypanda) of `phypanda`.


## Citation
If you use `phypanda`, please cite the corresponding paper:

> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks**.
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*
> bioRxiv, 2025. doi: [10.1101/2025.11.14.688467](https://www.biorxiv.org/content/10.1101/2025.11.14.688467)