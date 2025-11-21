
if __name__ == "__main__":
    from phypanda import *
    import pandas as pd
    import time, csv, os

    start = time.perf_counter()

    # Path to input file
    folder1 = "/path/to/folder/"
    file_path1 = folder1 + "exp1_simulated_networks.csv"
    
    # Path to output file
    folder2 = "/path/to/folder/"
    file_path2 = folder2 + "exp1_results.csv"

    # Load as dataframe (tab-delimited)
    df = pd.read_csv(file_path1, sep=" ", header=0)

    # open a results file once, before the loop
    if not os.path.exists(file_path2):
        with open(file_path2, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")

            # Header row
            header = ['row_id', 
                    'leaf_nr',
                    'scanwidth',
                    'level',
                    'ret_number',
                    'time_load_net',
                    'time_sw_load',
                    'time_sw_compute',
                    'time_sw_canonical',
                    'time_panda1',
                    'time_panda2',
                    'time_panda3',
                    'k1',
                    'k2',
                    'k3',
                    'diversity1',
                    'diversity2',
                    'diversity3']
            
            writer.writerow(header)

            for _, row in df.iterrows():

                # Load information
                row_id = row["id"]
                newick = row["newick"]
                leaf_nr = row["ntips"]
                level = row["level"]

                # Load network
                t1 = time.perf_counter()
                N = DirectedNetwork()
                N.load_from_enewick(newick)
                t2 = time.perf_counter()

                k1 = int(leaf_nr / 10)
                k2 = int(leaf_nr / 2)
                k3 = leaf_nr

                ret_number = len(N.reticulation_nodes())

                # Compute scanwidth
                t3 = time.perf_counter()
                G = DAG(nx.DiGraph(N))
                t4 = time.perf_counter()
                scanwidth, extension = G.optimal_scanwidth()
                t5 = time.perf_counter()
                tree_extension = extension.canonical_tree_extension()

                # Compute PD 1
                t6 = time.perf_counter()
                diversity1, solution = solve_MAPPD(N, k1, tree_extension=tree_extension)
                t7 = time.perf_counter()

                # Compute PD 2
                diversity2, solution = solve_MAPPD(N, k2, tree_extension=tree_extension)
                t8 = time.perf_counter()

                # Compute PD 3
                diversity3, solution = solve_MAPPD(N, k3, tree_extension=tree_extension)
                t9 = time.perf_counter()

                # Double check
                assert N.is_binary() == True
                assert len(N.leaves) == leaf_nr
                assert N.level() == level

                # Create statistics
                time_load_net = t2 - t1
                time_sw_load = t4 - t3
                time_sw_compute = t5 - t4
                time_sw_canonical = t6 - t5
                time_panda1 = t7 - t6
                time_panda2 = t8 - t7
                time_panda3 = t9 - t8

                info = [row_id, 
                        leaf_nr,
                        scanwidth,
                        level,
                        ret_number,
                        time_load_net,
                        time_sw_load,
                        time_sw_compute,
                        time_sw_canonical,
                        time_panda1,
                        time_panda2,
                        time_panda3,
                        k1,
                        k2,
                        k3,
                        diversity1,
                        diversity2,
                        diversity3]

                # Write one row to file
                writer.writerow(info)

    end = time.perf_counter()
    print(f"Finished in {end-start} seconds.")
