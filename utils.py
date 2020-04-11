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

PROJECT_ROOT_DIR = "."

PLOTS_PATH = os.path.join(PROJECT_ROOT_DIR, "plots")
os.makedirs(PLOTS_PATH, exist_ok=True)

EXPERIMENT_PATH = os.path.join(PROJECT_ROOT_DIR, "saved_runs")
os.makedirs(EXPERIMENT_PATH, exist_ok=True)


def read_metadata(lines, generator='netgen'):
    """
    Read metadata tags and values from a netgen files
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

    elif generator == 'netgen':
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
    np.random.seed(seed=99)
    try:
        
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
                    b[int(data[1]) - 1] = int(data[2]) 
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

                mu[i] = float(data[5])

                cov_coef = np.random.uniform(0.15, 0.3)

                sigma = mu[i] * cov_coef
                var[i] = sigma**2
                cap[i] = float(data[4]) 
                arc_dict[i] = (u, v)

                rows.append(u - 1)
                cols.append(i)
                values.append(1.0)
                rows.append(v - 1)
                cols.append(i)
                values.append(-1.0)

                G.nxg.add_edge(u, v, capacity=cap[i], mu=mu[i], var=var[i])
                i += 1

        G.mu = mu
        G.sigma = sigma
        G.var = var
        G.cap = cap
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

    return G


def save_run(fname, data, extension='pickle', prefix=None):

    if prefix is not None:
        path = os.path.join(EXPERIMENT_PATH, prefix)
        path = os.path.join(path, fname + "." + extension)
    else:
        path = os.path.join(EXPERIMENT_PATH, fname + "." + extension)

    with open(path, 'wb') as f:
        pickle.dump(data, f)


def load_run(fname, extension='pickle'):
    path = os.path.join(EXPERIMENT_PATH, fname + "." + extension)

    with open(path, 'rb') as f:
        item = pickle.load(f)

    return item
