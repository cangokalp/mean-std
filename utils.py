import os
import sys
import logging
import subprocess
import time
import pdb
import pickle
import scipy
import numpy as np
import networkx as nx
import mosek
import mosek.fusion as mf
from mosek.fusion import Expr, Set, Domain, Var, ObjectiveSense, Matrix

def gen_instance(**kwargs):
    """takes network parameters as input and outputs a file that can be read by the generator"""
    pass


def read_metadata(lines, generator='netgen'):
    """
    Read metadata tags and values from a TNTP file, returning a dictionary whose
    keys are the tags (strings between the <> characters) and corresponding values.
    The last metadata line (reading <END OF METADATA>) is stored with a value giving
    the line number this tag was found in.  You can use this to proceed with reading
    the rest of the file after the metadata.
    """
    if generator == 'goto':
        metadata = dict()
        lineNumber = 0
        for line in lines:
            lineNumber += 1
            line.strip()
            if line.find('Grid example:') >= 0:
                nodeTagPos = line.find('N = ')
                arcTagPos = line.find('M = ')
                maxcapTagPos = line.find('MAXCAP = ')
                maxcostTagPos = line.find('MAXCOST = ')
                seedTagPos = line.find('SEED = ')

                metadata['Number of nodes'] = line[
                    nodeTagPos + 4: arcTagPos].strip()
                metadata['Number of arcs'] = line[
                    arcTagPos + 4: maxcapTagPos].strip()
                metadata['Maximum arc cost'] = line[
                    maxcapTagPos + len('MAXCAP = '): maxcostTagPos].strip()
                metadata['Maximum arc capacity'] = line[
                    maxcostTagPos + len('MAXCOST = '): seedTagPos].strip()

            if line.find('p') == 0:
                metadata['END OF METADATA'] = lineNumber
                return metadata

    if generator == 'netgen':
        metadata = dict()
        lineNumber = 0
        for line in lines:
            lineNumber += 1
            line.strip()

            commentPos1 = line.find("NETGEN")
            commentPos2 = line.find("Problem")
            commentPos3 = line.find("-")

            if commentPos1 >= 0 or commentPos2 >= 0 or commentPos3 >= 0:  # strip comments
                # line = line[:commentPos]
                continue
            if len(line) == 1:
                continue

            startTagPos = line.find("c")
            endTagPos = line.find(":")

            if line.find('*** Minimum cost flow ***') >= 0:
                metadata['END OF METADATA'] = lineNumber
                return metadata

            if startTagPos < 0 or endTagPos < 0 or startTagPos >= endTagPos:
                print("Error reading this metadata line, ignoring: '%s'" % line)
                continue
            metadataTag = line[startTagPos + 1: endTagPos].strip()
            metadataValue = line[endTagPos + 1:]

            metadata[metadataTag] = metadataValue.strip()

        print("Warning: END OF METADATA not found in file")

    return metadata


def load_data(networkFileName, G, generator='netgen'):
    
    try:
        with open(networkFileName, "r") as networkFile:
            fileLines = networkFile.read().splitlines()

            # Set default parameters for metadata, then read

            metadata = read_metadata(fileLines, generator)

            n = int(metadata['Number of nodes'])
            m = int(metadata['Number of arcs'])
            max_arc_cost = int(metadata['Maximum arc cost'])
            max_arc_cap = int(metadata['Maximum arc capacity'])

            var = np.zeros(m)
            mu = np.zeros(m)
            sigma = np.zeros(m)
            cap = np.zeros(m)
            b = np.zeros(n)

            rows = []
            cols = []
            values = []
            i = 0
            arc_dict = {}
        
            if generator == 'netgen':
                min_arc_cost = int(metadata['Minimum arc cost'])
                total_supply = int(metadata['Total supply'])
                skltn_max_arc_perc = int(metadata['With max cost'][:-1])
                skltn_cap_arc_perc = int(metadata['Capacitated'][:-1])
                random_seed = int(metadata['Random seed'])

            firstlinedone = False
            extradone = False

            for line in fileLines[metadata['END OF METADATA']:]:
                # Ignore n and p and blank lines

                if generator == 'goto':

                    if line.find("n") == 0 and not extradone:
                        data = line.split()
                        if not firstlinedone:
                            source = int(data[1])
                            total_supply = int(data[2])
                            firstlinedone = True
                        else:
                            sink = int(data[1])
                            extradone = True

                line = line.strip()
                pPos = line.find("p")
                nPos = line.find("n")

                if nPos >= 0 or pPos >= 0:  # strip comments
                    # line = line[:commentPos]
                    continue
                if len(line) == 1:
                    continue

                data = line.split()

                u = int(data[1])
                v = int(data[2])
                mu[i] = float(data[5])
                if mu[i] == 0:
                    mu[i] = np.random.uniform(0, max_arc_cost)
                cov_coef = np.random.uniform(0.15, 0.3)
                sigma[i] = mu[i] * cov_coef
                var[i] = sigma[i]**2
                cap[i] = float(data[4])
                arc_dict[i] = (u, v)
                
                rows.append(u - 1)
                cols.append(i)
                values.append(1)
                rows.append(v - 1)
                cols.append(i)
                values.append(-1)
                
                G.nxg.add_edge(u, v, capacity=cap[i], mu=mu[i], var=var[i], std=sigma[i])
                i += 1

        nx.set_node_attributes(G.nxg, 0, 'demand')
        G.nxg.node[source]['demand'] = total_supply
        G.nxg.node[sink]['demand'] = -total_supply
        G.d_top = total_supply
        G.mu = mu
        G.sigma = sigma
        G.var = var
        G.cap = cap
        G.A_msk = Matrix.sparse(n, m, rows, cols, values)
        G.A = scipy.sparse.csc_matrix((values, (rows, cols)))

        b[source-1] = total_supply
        b[sink-1] = -total_supply
        G.b = b
        G.n = n
        G.m = m
        G.arc_dict = arc_dict

        ##add feasible flow for nmcc algo
        # nx.set_edge_attributes(G.nxg, 0, 'flow')
        # G.nxg.add_edge(source, sink, capacity=total_supply, flow=total_supply, mu=mu*1e3 ,var=(std**2)*1e3)
        return G

    except IOError:
        print("\nError reading network file %s" % networkFile)
        traceback.print_exc(file=sys.stdout)


class Logger(object):
    """
    Creates a class that will both print and log any
    output text. See https://stackoverflow.com/a/5916874
    for original source code. Modified to add date and
    time to end of file name.
    """

    def __init__(self, **kwargs):
        num_nodes = kwargs['num_nodes']
        # pname = 'saved_runs/' + str(num_nodes) + '/' + experiment_name
        pname = 'test'
        # lamBar = kwargs['lambar'] * kwargs['scale']
        lamBar = kwargs['lambar']

        lamstr = str(lamBar)
        lamstr = lamstr.replace(".", "-")
        narcs = kwargs['num_arcs']
        mu_top = kwargs['mu_top']
        var_top = kwargs['var_top']
        seed = kwargs['seed']
        save_extension = lamstr + '_' + \
            str(mu_top) + str(var_top) + '_' + str(narcs) + '_' + str(seed)
        # self.filename = pname + save_extension + '.txt'
        self.filename = pname + '.txt'

        if not os.path.exists(pname):
            os.makedirs(pname)

        self.log = open(self.filename, "w")

    def write(self, message):
        # self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass


def read_nmcc(filename):
    f = open(filename, 'r+')
    nmcc = []

    i = 0
    for line in f.readlines():
        i += 1
        if i == 1:
            nmcc_cost = line.strip()
        else:
            if line.find('@time_spent') >= 0:
                nmcc_time_ms = float(
                    line[line.find(':') + 1:]) * 1e-6
            else:
                nmcc.append(int(line) + 1)
                # nmcc.append(int(line))

    try:
        nmcc.append(nmcc[0])
    except:
        pass

    f.close()

    return nmcc, nmcc_time_ms, nmcc_cost


def solve_nmcc():
    
    args = ("../dev/msmcf/msmcf/hello_lemon /Users/cgokalp/repos/dev/msmcf/residual_graph.lgf", "-c")
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
    popen.wait()
    output = popen.stdout.read()


def output_graph(filename, g, G, ptype='reg'):
    g = g.nxg
    nodelist = [n for n in g.nodes]
    f = open(filename, 'w+')
    f.write('@nodes\n')
    f.write('label\n')
    for node in nodelist:
        f.write(str(node) + '\n')
    f.write('@arcs\n')
    f.write('                weight       label\n')
    jj = 0

    for u, v, e in g.edges(data=True):
        jj += 1
        try:
            flow = G.nxg[u][v]['flow']
            if ptype == 'sa':
                w = G.mu_scalar * G.nxg[u][v]['mu'] + \
                    2 * G.lam * flow * G.nxg[u][v]['var']
            elif ptype == 'xi':
                w = 2 * curvar * G.nxg[u][v]['flow'] + \
                    2 * G.lam * G.nxg[u][v]['xi'] * curvar
            else:
                w = G.nxg[u][v]['mu']
            w = w / 1e6
            f.write(str(u) + '       ' + str(v) + '       ' +
                    str((1.0) * w) + '       ' + str(u) + '\n')
        except:
            flow = G.nxg[v][u]['flow']
            if ptype == 'sa':
                w = G.mu_scalar * G.nxg[v][u]['mu'] + \
                    2 * G.lam * flow * G.nxg[v][u]['var']
            elif ptype == 'xi':
                w = 2 * curvar * G.nxg[v][u]['flow'] + \
                    2 * G.lam * G.nxg[v][u]['xi'] * curvar
            else:
                w = G.nxg[v][u]['mu']
            w = w / 1e6
            f.write(str(u) + '       ' + str(v) + '       ' + str((-1.0)
                                                                  * w) + '       ' + str(u) + '\n')
    f.write('@attributes\n')
    f.write('source ' + str(nodelist[0]) + '\n')
    f.write('target ' + str(nodelist[-1]) + '\n')
    f.close()


def save(list_to_save, item, num_nodes, experiment_name=''):

    pname = 'saved_runs/' + str(num_nodes) + '/' + experiment_name
    fname = pname + item + '.pickle'

    if not os.path.exists(pname):
        os.makedirs(pname)

    with open(fname, 'wb') as f:
        pickle.dump(list_to_save, f)


def load(item, num_nodes, experiment_name=''):
    pname = 'saved_runs/' + str(num_nodes) + '/' + experiment_name
    fname = pname + item + '.pickle'

    with open(fname, 'rb') as f:
        item = pickle.load(f)
    return item

