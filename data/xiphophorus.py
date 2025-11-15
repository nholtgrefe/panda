from panda import DirectedNetwork, all_paths_pd

xipho_newick = "((((((((((Xgordoni:1.3295084631587457,Xmeyeri:1.3295084631587457):0.0,Xcouchianus:1.329508093234352):6.999730834529853,Xvariatus:8.329238927764205):2.1769451514229345,Xevelynae:10.50618407918714):1.118605313770228,(Xxiphidium:7.2210504457107145,#H24:0.0):4.403738947246653):0.0,Xmilleri:11.624787067955268):4.296868586395352,Xandersi:15.92165565435062):0.9486610416497712,Xmaculatus:16.87031669600039):0.5723386247384958,((((Xmontezumae:7.221055986870681,(Xcortezi:5.485599585171238,((Xmalinche:5.485605240002155,Xbirchmanni:5.485605240002155):0.0)#H26:0.0):1.7354564016994427):0.0,((Xnigrensis:2.4303498026154564,Xmultilineatus:2.4303498026154564):0.19174715477323678,(Xpygmaeus:1.347820846400494,Xcontinens:1.347820846400494):1.2742761109881993):4.598960284156991):0.0,#H26:1.7354540549075645):0.0)#H24:10.2216024589192):2.1886232296055894,((Xclemenciae:11.254572014210282,Xmonticolus:11.254572014210282):6.4012991117391564,(#H25:1.6332001759073602,(Xsignum:10.266850863153604,((Xhellerii:8.633649976742058)#H25:1.6332013685506936,(Xalvarezi:8.362082334652573,Xmayae:8.362082334652573):1.9047690106401785):0.0):0.0):7.3890209733000205):1.975407424395037);"
N = DirectedNetwork()
N.load_from_enewick(xipho_newick)

#%%
sols = []
for k in range(1, len(N.leaves)+1):
    val, leaves = all_paths_pd(N, k, tree_extension="optimal_XP")
    sols.append(leaves)

# Print the marginal species for every increasing k and print No if a species is removed
prev_sol = set()
for i, sol in enumerate(sols):
    sol_set = set(sol)
    #print(sol_set)
    diff = sol_set - prev_sol
    if len(diff) == 0:
        print(f"k={i+1}: Yes")
    else:
        print(f"k={i+1}: Added {diff}")
    prev_sol = sol_set


#%%

indices = dict()
for leaf in N.leaves:
    index = 0
    for k in range(1, len(N.leaves)):
        s1, _ = all_paths_pd(N, k, tree_extension="optimal_XP")
        N_minus_leaf = N.subnetwork(N.leaves - {leaf})
        s2, _ = all_paths_pd(N_minus_leaf, k, tree_extension="optimal_XP")
        index += (s1 - s2)
    index = index / (len(N.leaves)-1)
    indices[leaf] = index

#%%

sorted_indices = sorted(indices.items(), key=lambda x:x[1])
for key, value in sorted_indices:
    print(key, value)

