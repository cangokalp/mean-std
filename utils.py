import os
import sys
import logging
import subprocess
import time
import pdb
import pickle
import scipy
import numpy as np
import gzip
import math

# import networkx as nx
# import mosek
# import mosek.fusion as mf
# from mosek.fusion import Expr, Set, Domain, Var, ObjectiveSense, Matrix

# Where to save the figures
GRAPH_LOC = '/Users/cgokalp/repos/dev/msmcf/residual_graph.lgf'
NMCC_LOC = '/Users/cgokalp/repos/dev/msmcf/msmcf/nmcc.txt'

PROJECT_ROOT_DIR = "."

PLOTS_PATH = os.path.join(PROJECT_ROOT_DIR, "plots")
os.makedirs(PLOTS_PATH, exist_ok=True)

EXPERIMENT_PATH = os.path.join(PROJECT_ROOT_DIR, "saved_runs")
os.makedirs(EXPERIMENT_PATH, exist_ok=True)

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

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
    np.random.seed(seed=9)
    try:
        # with open(networkFileName, "r") as networkFile:
        #         pdb.set_trace()
        with gzip.open(networkFileName + '.gz', 'rt') as networkFile:

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
            arc_dict = {}

            rows = []
            cols = []
            values = []
            i = 0

            if generator == 'netgen':
                min_arc_cost = int(metadata['Minimum arc cost'])
                total_supply = int(metadata['Total supply'])
                skltn_max_arc_perc = int(metadata['With max cost'][:-1])
                skltn_cap_arc_perc = int(metadata['Capacitated'][:-1])
                random_seed = int(metadata['Random seed'])

            for line in fileLines[metadata['END OF METADATA']:]:
                # Ignore n and p and blank lines

                if line.find("n") == 0:
                    data = line.split()
                    b[int(data[1]) - 1] = int(data[2]) / 10.0
                    continue

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

                mu[i] = float(data[5]) / 100.0
                if mu[i] == 0:
                    mu[i] = np.random.uniform(
                        max_arc_cost / 4.0, max_arc_cost / 1.5)
                cov_coef = np.random.uniform(0.15, 0.3)
                sigma = mu[i] * cov_coef
                var[i] = sigma**2
                cap[i] = float(data[4]) / 10.0
                arc_dict[i] = (u, v)

                rows.append(u - 1)
                cols.append(i)
                values.append(1.0)
                rows.append(v - 1)
                cols.append(i)
                values.append(-1.0)

                G.nxg.add_edge(u, v, capacity=cap[i], mu=mu[i], var=var[i])
                i += 1

        # nx.set_node_attributes(G.nxg, 0, 'demand')
        G.mu = mu
        G.sigma = sigma
        G.var = var
        G.cap = cap
        # G.A_msk = Matrix.sparse(n, m, rows, cols, values)
        G.A = scipy.sparse.csc_matrix((values, (rows, cols)))
        G.b = b
        G.n = n
        G.m = m
        G.rows = rows
        G.cols = cols
        G.values = values
        G.arc_dict = arc_dict

    except IOError:
        print("\nError reading network file %s" % networkFile)
        traceback.print_exc(file=sys.stdout)

    # for solving the mcf with networkx library
    # for line in fileLines[metadata['END OF METADATA']:]:
    #     # Ignore n and p and blank lines

    #     if line.find("n") == 0:
    #         data = line.split()
    #         G.nxg.node[int(data[1])]['demand'] = -int(data[2])

    return G


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


def output_graph(filename, g, G, ptype='sa'):
    g = g.nxg
    nodelist = g.nodes()
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
                w = G.nxg[u][v]['mu'] + \
                    2 * G.lam * flow * G.nxg[u][v]['var']
            elif ptype == 'xi':
                w = 2 * curvar * G.nxg[u][v]['flow'] + \
                    2 * G.lam * G.nxg[u][v]['xi'] * curvar
            else:
                w = G.nxg[u][v]['mu']
            w = w
            f.write(str(u) + '       ' + str(v) + '       ' +
                    str((1.0) * w) + '       ' + str(u) + '\n')
        except:
            flow = G.nxg[v][u]['flow']
            if ptype == 'sa':
                w = G.nxg[v][u]['mu'] + \
                    2 * G.lam * flow * G.nxg[v][u]['var']
            elif ptype == 'xi':
                w = 2 * curvar * G.nxg[v][u]['flow'] + \
                    2 * G.lam * G.nxg[v][u]['xi'] * curvar
            else:
                w = G.nxg[v][u]['mu']
            w = w 
            f.write(str(u) + '       ' + str(v) + '       ' + str((-1.0)
                                                                  * w) + '       ' + str(u) + '\n')
    f.write('@attributes\n')
    f.write('source ' + str(nodelist[0]) + '\n')
    f.write('target ' + str(nodelist[-1]) + '\n')
    f.close()


def save_run(fname, data, extension='pickle'):
    path = os.path.join(EXPERIMENT_PATH, fname + "." + extension)

    with open(path, 'wb') as f:
        pickle.dump(data, f)


def load_run(fname, extension='pickle'):
    path = os.path.join(EXPERIMENT_PATH, fname + "." + extension)

    with open(path, 'rb') as f:
        item = pickle.load(f)

    return item
