from NetworkInferenceWrapper import NetworkInferenceWrapper
import pandas as pd
import numpy as np
import subprocess
import os
import sys
import preprocessing as prp
import csv
test_suite = os.path.join(os.path.dirname(__file__))
sys.path.append(test_suite)

class ARACNEWrapper(NetworkInferenceWrapper):

    def _infer_network(self, expression_data):
        """Method to infer a network from expression data using the GENIE3 algorithm.

        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression data stored in a data frame with sample identifiers as indices and 
            gene symbols as column names.

        Returns
        -------
        inferred_network : pd.DataFrame
            A data frame whose entries in the column 'score' correspond to edge scores in inferred network. For undirected networks,
            columns 'node_lower' and 'node_upper' contain the gene symbols of the nodes that are connected by the edge.
            Fr directed networks, these columns are named 'source' and 'target'.
        """
        main = os.path.join(test_suite, '..')
        prefix = 'aracne'

        # -o: set output folder
        out_dir = os.path.join(main, 'temp', prefix)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_path = os.path.join(out_dir, 'network.txt')

        # remove columns with zero standard deviation and normalize columns to unit variance
        expression_data = expression_data.loc[:, (expression_data.std() != 0)]
        expression_data = prp.normalizeToUnitVariance(expression_data)

        # -e: save expression_data to csv
        expression_data = expression_data.T # ARACNe expects gene x sample data set
        gene_dict = dict(zip(expression_data.index, range(len(expression_data.index)))) # ARACNe expects numeric gene identifiers in first column
        expression_data.insert(loc=0, column='gene', value=[gene_dict[i] for i in expression_data.index])
        data_path = os.path.join(main, 'temp', 'aracne', f'{prefix}_expression_data.txt')
        expression_data.to_csv(data_path, sep='\t', index=False)

        # get regulators and remove such genes that are not present in expression_data
        ktf_path = os.path.join(main, 'data', 'regulators.csv')
        regulators = np.loadtxt(ktf_path, delimiter='\t', dtype=str)
        regulators = regulators[np.isin(regulators, expression_data.index)]
        regulators = [gene_dict[i] for i in regulators]
        regulator_path = os.path.join(main, 'temp', 'aracne', f'{prefix}_regulators.txt')
        pd.DataFrame(regulators).to_csv(regulator_path, sep='\t', index=False, header=False)

        # set parameters seed and p-value
        p = '1E-4'
        seed = '1'

        # run ARACNe:
        cur = os.getcwd()
        os.chdir(os.path.join(main, 'algorithms', 'ARACNe-AP'))
        exe = os.path.join('dist','aracne.jar')
        thresholdCommand = f'java -Xmx5G -jar {exe} -e {data_path}  -o {out_dir} --tfs {regulator_path} --pvalue {p} --seed {seed} --calculateThreshold'
        subprocess.run(thresholdCommand, shell=True)
        for i in range(3): # TODO: multiple or not?
            command = f'java -Xmx5G -jar {exe} -e {data_path}  -o {out_dir} --tfs {regulator_path} --pvalue {p} --seed {seed}'
            subprocess.run(command, shell=True)
        command = f'java -Xmx5G -jar {exe} -o {out_dir} --consolidate --nobonferroni' # TODO: why is the network empty if we do not add --nobonferroni ?
        subprocess.run(command, shell=True)
        os.chdir(cur)

        # get results
        network = pd.read_csv(out_path, sep='\t', index_col=False)
        inv_gene_dict = {v: k for k, v in gene_dict.items()}
        network['Regulator'] = [inv_gene_dict[i] for i in network['Regulator']]
        network['Target'] = [inv_gene_dict[i] for i in network['Target']]
        network = network.sort_values(by=['MI'], axis=0, ascending=False)

        network = network.rename({'Regulator': 'source', 'Target': 'target', 'MI': 'score'}, axis='columns')
        network['type'] = 'directed'
        
        # remove temporary directory
        subprocess.call('rm -r '+str(out_dir), shell=True)
        return network

    def _get_top_k_edges(self, i, k):
            """Method to return the top k edges for the inferred network for block i. 
            
            Parameters
            ----------
            i : int
                ID of partion block. Must fall in range(0, len(self.partition)).
            k : Maximal number of edges to be returned.
            
            Returns
            -------
            top_k_edges : set
                Set of tuples encoding edges. For edges without a sense, use tuples of form (<gene_1>, <gene_2>),
                where <gene_1> and <gene_2> are gene symbols. For edges with a sense (e.g., positive or negative
                correlation), use tuples of form (<gene_1>, <gene_2>, <sense>), where <sense> is either -1 or 1.
                For undirected edges, ensure that <gene_1> <= <gene_2> for all tuples contained in edge set.
            """
            block = self._inferred_networks[i]
            top_k_edges = []
            for j in range(k):
                source = block.iloc[j, 0]
                target = block.iloc[j, 1]
                top_k_edges.append((source, target)) # TODO p-vals?
            return set(top_k_edges)
