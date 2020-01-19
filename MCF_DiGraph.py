import numpy as np
import networkx as nx
import pdb
import queue
from collections import deque
import math

class MCF_DiGraph:

    """
    Extend networkx DiGraph class for our stochastic arc cost purposes
    """
    def __init__(self, lambar=0.7, residual=False):
        self.nxg = nx.DiGraph()
        self.residual = residual
        self.lam = 0
        self.initial_flow = []
        self.nmcc = []
        self.v_star = None
        self.pi = None
        self.K = None
        self.M = None
        self.nmcc_exists = None
        self.varcost = None
        self.xicost = None
        self.mu_scalar = 1
        self.lambar = lambar

    def reset(self):
        self.set_flow_attributes(self.initial_flow)


    def set_muscalar(self, muscalar):
        self.mu_scalar = muscalar

    def set_nmcc_lemon(self, nmcc):
        self.nmcc = nmcc
        # nmcc_edges = zip(nmcc[:-1], nmcc[1:])
        # self.subedges = list(nmcc_edges)
        # self.alledges = list(self.nxg.edges)
        #
        # diff = list(set(self.subedges).difference(set(self.alledges)))
        # diff = [(b, a) for a, b in diff]
        #
        # self.subedges.extend(diff)
        #
        # self.subg = self.nxg.edge_subgraph(self.subedges).copy()
        # diffn = self.nxg.edge_subgraph(diff).copy()
        #
        # mudict = nx.get_edge_attributes(diffn, 'mu')
        # for key in mudict:
        #     mudict[key] = -mudict[key]
        #
        # vardict = nx.get_edge_attributes(diffn, 'var')
        # for key in vardict:
        #     vardict[key] = -vardict[key]
        #
        # nx.set_edge_attributes(self.subg, mudict, 'mu')
        # nx.set_edge_attributes(self.subg, vardict, 'var')

    def set_lambda(self, lam=0, xi=False):
        self.lam = lam
        # if not xi:
        #     self.set_weights()
        # if xi:
        #     self.set_weights_xi()

    def get_lambda(self):
        return self.lam

    # def set_weights(self):
    #     var_cost = 0
    #     for u, v, e in self.nxg.edges(data=True):
    #         var_cost +=  np.square(e.get('flow', 0)) * e.get('var', 0)
    #         self.nxg[u][v]['weight'] = self.nxg[u][v]['mu'] + 2*self.lam*self.nxg[u][v]['flow']*self.nxg[u][v]['var']
    #     self.varcost = var_cost

    # def set_weights_xi(self):
    #     try:
    #         cost = 0
    #         for u, v, e in self.nxg.edges(data=True):
    #             self.nxg[u][v]['weight_xi'] = 2*self.nxg[u][v]['var']*self.nxg[u][v]['flow'] + 2*self.lam*self.nxg[u][v]['xi']*self.nxg[u][v]['var']
    #             cost += e.get('xi', 0) * e.get('var', 0) * e.get('flow', 0)
    #         self.xicost = cost
    #     except:
    #         pass

    def find_feasible_flow(self, d_top, mu_top, var_top):
        """Establish a feasible flow
        """
        # solve max flow problem

        G_MaxFlow = self.nxg.copy()
        G_MaxFlow.add_node('s')
        G_MaxFlow.add_node('t')
        for i, data in G_MaxFlow.nodes(data=True):
            if (i != 's') and (i != 't'):
                try:
                    b_i = G_MaxFlow.node[i]['demand']
                except:
                    pdb.set_trace()
                if b_i > 0:
                    G_MaxFlow.add_edge('s', i, capacity=b_i)
                elif b_i < 0:
                    G_MaxFlow.add_edge(i, 't', capacity=-b_i)

        cost, flow = nx.maximum_flow(G_MaxFlow, 's', 't')
        # pdb.set_trace()
        # cost, flow = nx.capacity_scaling(self.nxg)
        # print("The minimum cost is:", cost)

        # The data structure of flow is not consistent with dictionary data structure
        # needed for printing the flow vector
        feasible_flow = {}
        cap = {}
        var= {}
        mu ={}
        for i in self.nxg.nodes():
            for j in flow[i].keys():
                if j != 't':
                    if abs(round(flow[i][j],4) - flow[i][j]) <= 1e-10:
                        flow[i][j] = round(flow[i][j],4)
                    feasible_flow[i,j] = flow[i][j]

                    # if flow[i][j] == 0:
                    #     mu[i,j] = 0
                    #     var[i,j]= 0
                    #     cap[i,j] = 0
                    # else:
                    #     mu[i,j] = mu_top*10000
                    #     var[i,j]= var_top*10000
                    #     cap[i,j] = d_top

                    if feasible_flow[i,j] <=1e-12 and feasible_flow[i,j] >0:
                        pdb.set_trace()
        # print(feasible_flow)
        self.initial_flow = feasible_flow
        self.set_flow_attributes(self.initial_flow)
        # self.set_parameters(cap, mu, var)


        # n = len(self.nxg.nodes())
        # m = len(self.nxg.edges())

        # self.var = {}
        # self.cap = {}
        # self.mu = {}

        # self.A = np.zeros((n,m))
        # i = 0
        # for u, v, e in self.nxg.edges(data=True):
        #     self.var[(u,v)] = e.get('var', 0)
        #     self.cap[(u,v)] = e.get('capacity', 0)
        #     self.mu[(u,v)] = e.get('mu', 0)
        #     self.A[u-1, i] = 1
        #     self.A[v-1, i] = -1
        #     i += 1

        # self.b = {}
        # i = 0
        # for node, dat in self.nxg.nodes(data=True):
        #     self.b[node] = dat['demand']
        #     i += 1

        # self.set_weights()

    def find_feasible_flow_xi(self):
        """Establish a feasible flow
        """
        # solve max flow problem

        G_MaxFlow = self.nxg.copy()
        G_MaxFlow.add_node('s')
        G_MaxFlow.add_node('t')
        M = 5
        for i, data in G_MaxFlow.nodes(data=True):
            if (i != 's') and (i != 't'):
                ## fix
                b_i = 0
                # plus = 0
                # for v in G_MaxFlow.neighbors(i):
                #     plus += G_MaxFlow[i][v]['flow']
                    # if self.nxg[i][v]['flow'] < 1e-6:
                    #     plus += 0
                    # else:
                    #     plus += self.nxg[i][v]['flow']



                # minus = 0
                # for e in G_MaxFlow.in_edges(i):
                #     minus += G_MaxFlow[e[0]][e[1]]['flow']

                    # if self.nxg[e[0]][e[1]]['flow'] < 1e-6:
                    #     minus += 0
                    # else:
                    #     minus += self.nxg[e[0]][e[1]]['flow']


                # b_i = b_i + plus - minus

                if b_i > 0:
                    G_MaxFlow.add_edge('s', i, capacity=b_i)
                elif b_i < 0:
                    G_MaxFlow.add_edge(i, 't', capacity=-b_i)

        cost, flow = nx.maximum_flow(G_MaxFlow, 's', 't')

        # cost, flow = nx.capacity_scaling(self.nxg)
        # print("The minimum cost is:", cost)

        # The data structure of flow is not consistent with dictionary data structure
        # needed for printing the flow vector
        feasible_flow = {}
        for i in self.nxg.nodes():
            for j in flow[i].keys():
                if j != 't':
                    if abs(round(flow[i][j],4) - flow[i][j]) < 1e-10:
                        flow[i][j] = round(flow[i][j],4)
                    feasible_flow[i,j] = flow[i][j]


        self.initial_flow = feasible_flow
        self.set_flow_attributes_xi(self.initial_flow)
        # self.set_weights_xi()

    def set_parameters(self, cap, mu, var):
        """Set flow variables for the graph
        """
        nx.set_edge_attributes(self.nxg, cap, 'capacity')
        nx.set_edge_attributes(self.nxg, mu, 'mu')
        nx.set_edge_attributes(self.nxg, var, 'var')


    def set_flow_attributes(self, feasible_flow):
        """Set flow variables for the graph
        """
        nx.set_edge_attributes(self.nxg, feasible_flow, 'flow')

    def set_flow_attributes_xi(self, feasible_flow):
        """Set flow variables for the graph
        """
        nx.set_edge_attributes(self.nxg, feasible_flow, 'xi')


    def draw_graph(self):
        """Visualize the graph.
        """

        pos = nx.spring_layout(self.nxg)

        # draw nodes, edges and labels
        nx.draw_networkx_nodes(self.nxg, pos, node_size=500, node_color='green', alpha=0.5)
        nx.draw_networkx_edges(self.nxg, pos, width=3, alpha=0.5, edge_color='skyblue')
        nx.draw_networkx_labels(self.nxg, pos, font_size=12, font_family='sans-serif')

        if not self.residual:
            # show the current flow on the arcs
            edge_labels = nx.get_edge_attributes(self.nxg, 'capacity')
        else:
            # show r_ij if it's a residual graph
            edge_labels = nx.get_edge_attributes(self.nxg, 'r')

        nx.draw_networkx_edge_labels(self.nxg, pos, edge_labels=edge_labels)

    def get_soln(self):
        x = []
        for u,v,e in self.nxg.edges(data=True):
            x.append(e.get('flow', 0))
        return x

    def tot_cost(self):
        return self.mu_cost() + self.lambar * np.sqrt(self.var_cost())

    def mu_cost(self):

        # flow = np.array(list(nx.get_edge_attributes(self.nxg,'flow').values()))
        # mu = np.array(list(nx.get_edge_attributes(self.nxg,'mu').values()))
        # return np.dot(flow, mu)
        cost = 0
        for u,v,e in self.nxg.edges(data=True):
            cost += e.get('flow', 0) *e.get('mu', 0)
        return cost


    def var_cost(self):

        # flow = np.square(np.array(list(nx.get_edge_attributes(self.nxg,'flow').values())))
        # var = np.array(list(nx.get_edge_attributes(self.nxg,'var').values()))
        # return np.dot(flow, var)

        cost = 0
        for u,v,e in self.nxg.edges(data=True):
            cost += np.square(e.get('flow', 0)) * e.get('var', 0)
        return cost

    def xi_cost(self):

        # flow = np.square(np.array(list(nx.get_edge_attributes(self.nxg,'flow').values())))
        # var = np.array(list(nx.get_edge_attributes(self.nxg,'var').values()))
        # xi = np.array(list(nx.get_edge_attributes(self.nxg,'xi').values()))
        # return np.sum(np.multiply(np.multiply(xi,var), flow))
        M = 50000
        cost = 0
        for u,v,e in self.nxg.edges(data=True):

                # e.get('xi',0)-e.get('flow',0)
            # print()
            #     xi = e.get('xi',0)
            #     cost += xi * e.get('var', 0) * e.get('flow', 0)
            

            # else:
            #     xi = e.get('xi',0)
            #     cost += (xi-M) * e.get('var', 0) * e.get('flow', 0)



            cost += e.get('xi',0) * e.get('var', 0) * e.get('flow', 0)

        return cost

    # def tot_cost(self):

    #     cost = 0
    #     for u,v,e in self.nxg.edges(data=True):
    #         cost += e.get('flow', 0) *e.get('mu', 0) + self.lam*np.square(e.get('flow', 0)) * e.get('var', 0)

    #     return cost
    def set_weight_uc(self, lam):
        for u,v,e in self.nxg.edges(data=True):
            e['weight'] = e['flow'] + 2*lam*e['var']*e['flow']
            
    def FindSource(self, g):

        F = np.zeros(len(self.nxg.nodes()))
        L = np.zeros(len(self.nxg.nodes()))
        for i in range(10):
            for v in g.nodes():
                L[int(v-1)] = g.degree(v)
                for u in g.neighbors(v):
                    L[int(v-1)] = L[int(v-1)] + F[int(u-1)]
            F = L
        minval = np.min(L[np.nonzero(L)])
        L = L.tolist()
        s = L.index(minval) + 1

        return s

    def build_res_network(self, capacity='capacity', flow='flow', jw=False):
        """Build a residual network and initialize a zero flow.
        """
        R = MCF_DiGraph(residual=True)

        for u, v, e in self.nxg.edges(data=True):
            r_vu = e.get(flow, 0)
            r_uv = e.get(capacity,0) - e.get(flow,0)

            if r_vu <= 1e-14 and r_uv==e.get(capacity,0):
                r_vu = 0
            if abs(round(r_uv,4) - r_uv) <= 1e-10:
                r_uv = round(r_uv,4)

            # r_uv = fix_frac(e.get(capacity,0) - e.get(flow,0))

            # mu_val = e.get(mu, 0)
            # var_val = e.get(var, 0)
            # flow_val = e.get(flow, 0)
            # cap_val = e.get(capacity, 0)
            # weight_val = e.get(weight,0)

            if not jw:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, r=r_uv)

                if r_vu != 0:
                    R.nxg.add_edge(v, u, r=r_vu)
            else:

                if r_uv != 0:
                    R.nxg.add_edge(u, v, r=r_uv)
                if r_vu != 0:
                    R.nxg.add_edge(v, u, r=r_vu)
        return R

    def build_res_network_xi(self, capacity='capacity', oldflow='flow', flow='xi', jw=False):
        """Build a residual network and initialize a zero flow.
        """
        R = MCF_DiGraph(residual=True)
        ub = 50000
        lb = -50000

        for u, v, e in self.nxg.edges(data=True):
 
            if abs(e.get(oldflow, 0) - e.get(capacity,0))<= 1e-10:
                # r_vu = e.get(flow, 0)
                r_uv = 0 - e.get(flow,0)
                r_vu = e.get(flow,0) - lb
                # r_uv = e.get(oldflow,0) - e.get(flow,0)
                # r_uv = M - e.get(flow,0)

            elif e.get(oldflow,0) < 1e-6:
                r_uv = 0 - e.get(flow,0)
                r_vu = e.get(flow,0)
                # r_uv = ub - e.get(flow,0)
                # r_vu = e.get(flow,0)
                # print(e.get(flow,0))
                # r_uv = e.get(capacity,0) - e.get(flow,0)
                # r_uv = e.get(flow,0)

            else:
                r_uv = ub - e.get(flow,0)
                r_vu = e.get(flow,0) - lb
                # r_vu = e.get(flow, 0)
                # r_uv = e.get(capacity,0) + e.get(oldflow,0) - e.get(flow,0)
                # r_uv = 2*M - e.get(flow,0)


            if abs(r_vu) <= 1e-9 and abs(r_uv-e.get(capacity,0))<=1e-9:
                r_vu = 0
            if abs(round(r_uv,4) - r_uv) <=1e-10:
                r_uv = round(r_uv,4)

            # mu_val = e.get(mu, 0)
            # var_val = e.get(var, 0)
            # flow_val = e.get(flow, 0)
            # cap_val = e.get(capacity, 0) + e.get(oldflow,0)
            # weight_val = e.get(weight,0)

            if not jw:
                if r_uv != 0:
                    # R.nxg.add_edge(u, v, capacity=cap_val, r=r_uv)
                    R.nxg.add_edge(u, v, r=r_uv)
                if r_vu != 0:
                    # R.nxg.add_edge(v, u, capacity=cap_val, r=r_vu)
                    R.nxg.add_edge(v, u, r=r_vu)
            else:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, r=r_uv)
                if r_vu != 0:
                    R.nxg.add_edge(v, u, r=r_vu)
        return R

    def get_nmcc(self):

        return self.nmcc

    def find_delta(self):
        try:
            min_r = np.inf
            for i in range(len(self.nmcc)-1):
                u = self.nmcc[i]
                v = self.nmcc[i+1]
                rij = self.nxg[u][v]['r']
                min_r = min(min_r, rij)
        except:
            pdb.set_trace()

        # min_r = min(np.array(list(nx.get_edge_attributes(self.subg,'r').values())))
        return min_r

    def find_xsi_star_xi(self, G):

        nominator = 0
        denominator = 0

        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]

            if G.nxg.has_edge(u,v):
                flow = G.nxg[u][v]['flow']
                # weight = G.nxg[u][v]['weight_xi']
                var = G.nxg[u][v]['var']
                nominator -= 2*var*G.nxg[u][v]['flow'] + 2*G.lam*G.nxg[u][v]['xi']*var

            else:
                flow = G.nxg[v][u]['flow']
                # weight = -G.nxg[v][u]['weight_xi']
                var = G.nxg[v][u]['var']
                nominator += 2*var*G.nxg[v][u]['flow'] + 2*G.lam*G.nxg[v][u]['xi']*var

            # nominator -= weight
            # denominator += 2*G.lam*abs(self.nxg[u][v]['var'])
            denominator += 2*G.lam*abs(var)
        try:
            cc = nominator/denominator
            # cc = math.floor(cc * 1000)/1000.0
        except:
            pdb.set_trace()


        return cc

    def find_xsi_star(self, G):

        # denominator = 2*G.lam*(sum(abs(np.array(list(nx.get_edge_attributes(G.subg,'var').values())))))
        # nominator = sum(np.array(list(nx.get_edge_attributes(G.subg,'mu').values()))) + 2*G.lam*(sum(np.array(list(nx.get_edge_attributes(G.subg,'flow').values()))* np.array(list(nx.get_edge_attributes(G.subg,'var').values()))))
        # cc= -nominator/denominator

        # pdb.set_trace()

        nominator = 0
        denominator = 0
        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]

            if G.nxg.has_edge(u,v):
                flow = G.nxg[u][v]['flow']
                # mu = self.mu[(u,v)] 
                # var = self.var[(u,v)] 

                # weight = G.nxg[u][v]['weight']
                var = G.nxg[u][v]['var']
                nominator -= self.mu_scalar*G.nxg[u][v]['mu'] + 2*G.lam*flow*var
                # nominator -= self.mu_scalar*mu + 2*G.lam*flow*var


            else:
                # weight = -G.nxg[v][u]['weight']
                flow = G.nxg[v][u]['flow']
                # mu = self.mu[(v,u)] 
                # var = self.var[(v,u)] 
                var = G.nxg[v][u]['var']
                nominator += self.mu_scalar*G.nxg[v][u]['mu'] + 2*G.lam*flow*var

            # nominator -= weight
            # denominator += 2*G.lam*abs(self.nxg[u][v]['var'])
            denominator += 2*G.lam*abs(var)
        try:
            cc = nominator/denominator
        except:
            pdb.set_trace()

        return cc

    def remove_no_r(self, u, v):
        if self.nxg[u][v]['r'] <= 1e-9:
            self.nxg.remove_edge(u,v)

    def augment_flow(self, min_r):
        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]

            cur_r = self.nxg[u][v]['r']
            if abs(round(cur_r - min_r, 4) - (cur_r- min_r)) <=1e-10:
                self.nxg[u][v]['r'] = round(cur_r - min_r, 4)
            else:
                self.nxg[u][v]['r'] = cur_r - min_r

            self.remove_no_r(u, v)

            if self.nxg.has_edge(v, u):
                cur_r = self.nxg[v][u]['r']
                if abs(round(cur_r + min_r, 4) - (cur_r + min_r)) <=1e-10:
                    self.nxg[v][u]['r'] = round(cur_r + min_r, 4)
                else:
                    self.nxg[v][u]['r'] = cur_r + min_r

                self.remove_no_r(v, u)

            else:
                self.nxg.add_edge(v, u, r=min_r)
                self.remove_no_r(v, u)


    def adjust_flow(self, nmcc, delta):

        for i in range(len(nmcc)-1):
            u = nmcc[i]
            v = nmcc[i+1]

            if self.nxg.has_edge(u, v):
                flow = self.nxg[u][v]['flow']
                if abs(round(flow + delta,4) - (flow + delta)) <=1e-10:
                    self.nxg[u][v]['flow'] = round(flow + delta, 4)
                else:
                    self.nxg[u][v]['flow'] = flow + delta

            else:
                flow = self.nxg[v][u]['flow']
                if abs(round(flow - delta,4) - (flow - delta)) <=1e-10:
                    self.nxg[v][u]['flow'] = round(flow - delta, 4)
                else:
                    self.nxg[v][u]['flow'] = flow - delta

    def adjust_flow_xi(self, nmcc, delta):

        for i in range(len(nmcc)-1):
            u = nmcc[i]
            v = nmcc[i+1]

            # if self.nxg.has_edge(u, v):
            #     self.nxg[u][v]['xi'] = self.nxg[u][v]['xi'] + delta
            # else:
            #     self.nxg[v][u]['xi'] = self.nxg[v][u]['xi'] - delta

            if self.nxg.has_edge(u, v):
                xi = self.nxg[u][v]['xi']
                if abs(round(xi + delta,4) - (xi + delta)) <=1e-10:
                    self.nxg[u][v]['xi'] = round(xi + delta, 4)
                else:
                    self.nxg[u][v]['xi'] = xi + delta

            else:
                xi = self.nxg[v][u]['xi']
                if abs(round(xi - delta,4) - (xi - delta)) <=1e-10:
                    self.nxg[v][u]['xi'] = round(xi - delta, 4)
                else:
                    self.nxg[v][u]['xi'] = xi - delta


    def check_nmcc(self, G, lam, g=None ,cycle_len=0):
        """Temporary solution
        """

        tol = 1e-10
        valid = []
        # for c in nx.simple_cycles(self.nxg):
        for c in nx.simple_cycles(g):
            # if len(c) == cycle_len:
            if len(c) != 2:
                # pass
            # else:
                c.append(c[0])
                valid.append(c)

        if len(valid) == 0:
            return False
        else:

            path_cost = []
            for c in valid:
                cost = 0
                for i in range(len(c)-1):
                    u = c[i]
                    v = c[i+1]

                    if G.has_edge(u,v):
                        cost += G[u][v]['weight']
                    else:
                        cost += -G[v][u]['weight']

                    # try:


                        # cost += self.nxg[u][v]['mu'] + 2*lam*G[u][v]['flow']*self.nxg[u][v]['var']
                        # cost += g[u][v]['mu'] + 2*lam*G[u][v]['flow']*g[u][v]['var']

                    # except:
                        # cost += self.nxg[u][v]['mu'] + 2*lam*G[v][u]['flow']*self.nxg[u][v]['var']
                        # cost += g[u][v]['mu'] + 2*lam*G[v][u]['flow']*g[u][v]['var']

                path_cost.append(cost/(len(c)-1))

        W = valid[np.argmin(path_cost)]
        print('=======ENUMERATION==========')
        print('============================')
        print('| Cycle: ', W, ' | Cost: ', np.min(path_cost), ' |')
        print('============================')

        # cycle_cost = []
        possible = []
        for i in range(len(valid)):
        #     cycle_cost[i] = path_cost[i]/len(valid[i])

        # idx = argmin(cycle_cost)
        # print('cycle: ', valid[idx])
        # print('mean cost: ', path_cost[idx])
        # print('=====')

            if path_cost[i]/len(valid[i]) < 0 and len(valid[i])==(cycle_len+1):
                possible.append((valid[i], path_cost[i]))
                # mincost = min(mincost, cycle_cost)

        print('possible cycles with correct length')
        for cycle, cost in possible:
            print(cycle, cost)
        #         print(valid[i])
        #         print(path_cost[i])
        #         print('=====')

        # print('=======now pos======')
        # for i in range(len(valid)):

        #     if path_cost[i] >= 0:
        #         print(valid[i])
        #         print(path_cost[i]/(len(valid[i])-1))
        #         print('=====')


            # if  np.min(path_cost) < 0:
            #     W = valid[np.argmin(path_cost)]
            #     self.nmcc_cost = np.min(path_cost)
            #     self.nmcc = W

            #     pdb.set_trace()
            #     print(W)

            #     delta = self.find_delta()
            #     if delta <= 0:
            #         return False

            #     if abs(self.nmcc_cost) < tol:
            #         self.nmcc_cost = 0.0
            #         return False
            #     else:
            #         return True
            # else:
            #     return False

    def DG_3(self, G):
        g = self.nxg
        n = len(g.nodes())
        m = len(g.edges())

        if not nx.is_strongly_connected(g):
            scc = [c for c in sorted(nx.strongly_connected_components(g),key=len, reverse=True)]
            scc_count = len(scc)
            min_lam = []
            M_l = []
            K_l = []
            v_l =[]
            pi_l = []
            glen_l =[]
            curg_l = []
            dkv_l = []
            for i in range(scc_count):
                nodelist = scc[i]
                if len(nodelist) <= 1:
                    continue

                cur_g = self.nxg.copy()

                allnodes = np.arange(len(self.nxg.nodes)) + 1
                removals = list(set(allnodes)-set(nodelist))
                cur_g.remove_nodes_from(removals)
                s = self.FindSource(cur_g)
                nodelist = list(nodelist)

                len_cur_g = len(nodelist)
                d_kv = np.full((n+1,n), np.inf)
                pi = np.full((n+1, n), -1)
                M = np.zeros(n)
                K = np.zeros(n)

                valid = np.zeros((n+1, n))
                lastlevel = np.full(n, -1)
                d_kv[0, s-1] = 0
                valid[0, s-1] = -1
                lastlevel[s-1] = 0
                q = deque()
                q.append((0,s))
                v_star = None

                ## Body
                k, v = q.popleft()
                while k < len_cur_g:
                    for u in cur_g.neighbors(v):
                        if lastlevel[u-1] < k+1:
                            q.append((k+1,u))
                            valid[k+1, u-1] = lastlevel[u-1]
                            lastlevel[u-1] = k+1

                        if G.nxg.has_edge(v,u):
                            w_vu = G.nxg[v][u]['weight']
                        else:
                            w_vu = -G.nxg[u][v]['weight']
                        pdb.set_trace()
                        if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                            d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                            pi[k+1, u-1] = v
                    k,v = q.popleft()
                ## Tail

                lam = np.inf
                for v in cur_g.nodes():
                    if lastlevel[v-1] == n:
                        M[v-1] = -np.inf
                        k = valid[len_cur_g, v-1]
                        while k > -1:
                            if M[v-1] < (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k):
                                M[v-1] = (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k)
                                K[v-1] = k
                            k = valid[k, v-1]
                        if lam > M[v-1]:
                            lam = M[v-1]
                            v_star = v

                min_lam.append(lam)
                M_l.append(M)
                K_l.append(K)
                v_l.append(v_star)
                pi_l.append(pi)
                glen_l.append(len_cur_g)
                curg_l.append(cur_g)
                dkv_l.append(d_kv)

            min_ind = np.argmin(min_lam)
            lam = min_lam[min_ind]
            M = M_l[min_ind]
            K = K_l[min_ind]
            v_star = v_l[min_ind]
            pi = pi_l[min_ind]
            glen = glen_l[min_ind]
            g = curg_l[min_ind]
            self.d_kv = dkv_l[min_ind]
            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            pdb.set_trace()
            print(lam)
            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

        else:

            d_kv = np.full((n+1,n), np.inf)
            pi = np.full((n+1, n), -1)
            valid = np.zeros((n+1, n))
            M = np.zeros(n)
            K = np.zeros(n)

            lastlevel = np.full(n, -1)
            s = self.FindSource(g)
            d_kv[0, s-1] = 0
            valid[0, s-1] = -1
            lastlevel[s-1] = 0

            q = deque()
            q.append((0,s))
            ## Body
            k,v = q.popleft()
            while k < n:
                for u in g.neighbors(v):
                    if lastlevel[u-1] < k+1:
                        q.append((k+1,u))
                        valid[k+1, u-1] = lastlevel[u-1]
                        lastlevel[u-1] = k+1

                    if G.nxg.has_edge(v,u):
                        w_vu = G.nxg[v][u]['weight']
                    else:
                        w_vu = -G.nxg[u][v]['weight']

                    if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                        d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                        pi[k+1, u-1] = v
                k,v = q.popleft()

            ## Tail
            lam = np.inf
            for v in g.nodes():
                if lastlevel[v-1] == n:
                    M[v-1] = -np.inf
                    k = int(valid[n, v-1])
                    while k > -1:
                        if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                            M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                            K[v-1] = k
                        k = int(valid[k, v-1])
                    if lam > M[v-1]:
                        lam = M[v-1]
                        self.v_star = v

            self.M = M
            self.K = K
            self.pi = pi
            self.lam_star = lam
            self.g = g
            self.d_kv = d_kv
            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

    def DG_2(self, G):
        g = self.nxg

        if not nx.is_strongly_connected(g):
            scc = [c for c in sorted(nx.strongly_connected_components(g),key=len, reverse=True)]
            scc_count = len(scc)
            min_lam = []
            M_l = []
            K_l = []
            v_l =[]
            pi_l = []
            glen_l =[]
            curg_l = []
            dkv_l = []
            for i in range(scc_count):
                nodelist = scc[i]
                if len(nodelist) <= 1:
                    continue

                cur_g = self.nxg.copy()

                allnodes = np.arange(len(self.nxg.nodes)) + 1
                removals = list(set(allnodes)-set(nodelist))
                cur_g.remove_nodes_from(removals)

                s = self.FindSource(cur_g)
                nodelist = list(nodelist)
                n = len(g.nodes())
                m = len(g.edges())
                len_cur_g = len(nodelist)
                d_kv = np.full((n+1,n), np.inf)
                pi = np.full((n+1, n),-1)
                Visit = np.zeros((n, 2),dtype=bool)
                M = np.zeros(n)
                K = np.zeros(n)
                v_star = None

                firstv = nodelist[0]
                lastv = nodelist[-1]


                d_kv[0, s-1] = 0
                # pi[0, s-1] = None
                turn = 0
                Visit[s-1, turn] = True
                ## Body

                for k in range(len_cur_g):
                # for k in range(firstv-1,lastv-firstv+1):
                # for k in range(lastv-firstv+1):
                    for v in nodelist:
                        if Visit[v-1, turn] == True:
                            Visit[v-1, turn] = False
                            for u in cur_g.neighbors(v):

                                if G.nxg.has_edge(v,u):
                                    w_vu = G.nxg[v][u]['weight']
                                else:
                                    w_vu = -G.nxg[u][v]['weight']

                                if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                    d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                    pi[k+1, u-1] = v
                                    Visit[u-1, 1-turn] = True
                    turn = 1 - turn
                ## Tail
                lam = np.inf
                for v in nodelist:
                    if Visit[v-1, turn] == True:
                        M[v-1] = -np.inf
                        for k in range(len_cur_g):
                        # for k in range(firstv-1,lastv-firstv+1):
                        # for k in range(lastv-firstv+1):

                            if M[v-1] < (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k):
                                M[v-1] = (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k)
                                K[v-1] = k

                        if lam > M[v-1]:
                            lam = M[v-1]
                            v_star = v
                min_lam.append(lam)
                M_l.append(M)
                K_l.append(K)
                v_l.append(v_star)
                pi_l.append(pi)
                glen_l.append(len_cur_g)
                curg_l.append(cur_g)
                dkv_l.append(d_kv)

            min_ind = np.argmin(min_lam)
            lam = min_lam[min_ind]
            M = M_l[min_ind]
            K = K_l[min_ind]
            v_star = v_l[min_ind]
            pi = pi_l[min_ind]
            glen = glen_l[min_ind]
            g = curg_l[min_ind]
            self.d_kv = dkv_l[min_ind]
            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

        else:

            n = len(g.nodes())
            m = len(g.edges())
            d_kv = np.full((n+1,n), np.inf)
            pi = np.full((n+1, n), -1)
            Visit = np.zeros((n, 2),dtype=bool)
            M = np.zeros(n)
            K = np.zeros(n)

            s = self.FindSource(g)
            d_kv[0, s-1] = 0

            # pi[0, s-1] = None
            turn = 0
            Visit[s-1, turn] = True

            ## Body
            for k in range(n):
                for v in g.nodes():
                    if Visit[v-1, turn] == True:
                        Visit[v-1, turn] = False
                        for u in g.neighbors(v):

                            if G.nxg.has_edge(v,u):
                                w_vu = G.nxg[v][u]['weight']
                            else:
                                w_vu = -G.nxg[u][v]['weight']

                            if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                pi[k+1, u-1] = v
                                Visit[u-1, 1-turn] = True
                turn = 1 - turn

            ## Tail
            lam = np.inf
            for v in g.nodes():
                if Visit[v-1, turn] == True:
                    M[v-1] = -np.inf
                    for k in range(n):
                        if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                            M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                            K[v-1] = k

                    if lam > M[v-1]:
                        lam = M[v-1]
                        v_star = v

            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            self.d_kv = d_kv
            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

    def DG_2_dict(self, G):
        g = self.nxg
        n = len(g.nodes())
        m = len(g.edges())
        orig_dict = nx.to_dict_of_dicts(G.nxg)
        if not nx.is_strongly_connected(g):
            scc = [c for c in sorted(nx.strongly_connected_components(g),key=len, reverse=True)]
            scc_count = len(scc)
            min_lam = []
            M_l = []
            K_l = []
            v_l =[]
            pi_l = []
            glen_l =[]
            curg_l = []
            dkv_l = []
            for i in range(scc_count):
                nodelist = scc[i]
                if len(nodelist) <= 1:
                    continue

                cur_g = self.nxg.copy()

                allnodes = np.arange(len(self.nxg.nodes)) + 1
                removals = list(set(allnodes)-set(nodelist))
                cur_g.remove_nodes_from(removals)

                s = self.FindSource(cur_g)
                nodelist = list(nodelist)

                cur_r_dict = nx.to_dict_of_dicts(g)

                len_cur_g = len(nodelist)
                d_kv = np.full((n+1,n), np.inf)
                pi = np.full((n+1, n),-1)
                Visit = np.zeros((n, 2),dtype=bool)
                M = np.zeros(n)
                K = np.zeros(n)
                v_star = None

                # firstv = nodelist[0]
                # lastv = nodelist[-1]


                d_kv[0, s-1] = 0
                # pi[0, s-1] = None
                turn = 0
                Visit[s-1, turn] = True
                ## Body

                for k in range(len_cur_g):
                # for k in range(firstv-1,lastv-firstv+1):
                # for k in range(lastv-firstv+1):
                    for v in nodelist:
                        if Visit[v-1, turn] == True:
                            Visit[v-1, turn] = False
                            for u in cur_r_dict[v].keys():

                                try:
                                    w_vu = orig_dict[v][u]['weight']
                                except:
                                    w_vu = -orig_dict[u][v]['weight']

                                if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                    d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                    pi[k+1, u-1] = v
                                    Visit[u-1, 1-turn] = True
                    turn = 1 - turn
                ## Tail
                lam = np.inf
                for v in nodelist:
                    if Visit[v-1, turn] == True:
                        M[v-1] = -np.inf
                        for k in range(len_cur_g):
                        # for k in range(firstv-1,lastv-firstv+1):
                        # for k in range(lastv-firstv+1):

                            if M[v-1] < (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k):
                                M[v-1] = (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k)
                                K[v-1] = k

                        if lam > M[v-1]:
                            lam = M[v-1]
                            v_star = v
                min_lam.append(lam)
                M_l.append(M)
                K_l.append(K)
                v_l.append(v_star)
                pi_l.append(pi)
                glen_l.append(len_cur_g)
                curg_l.append(cur_g)
                dkv_l.append(d_kv)

            min_ind = np.argmin(min_lam)
            lam = min_lam[min_ind]
            M = M_l[min_ind]
            K = K_l[min_ind]
            v_star = v_l[min_ind]
            pi = pi_l[min_ind]
            glen = glen_l[min_ind]
            g = curg_l[min_ind]
            self.d_kv = dkv_l[min_ind]
            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            if lam < 0:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

        else:

            d_kv = np.full((n+1,n), np.inf)
            pi = np.full((n+1, n), -1)
            Visit = np.zeros((n, 2),dtype=bool)
            M = np.zeros(n)
            K = np.zeros(n)

            s = self.FindSource(g)
            d_kv[0, s-1] = 0

            # pi[0, s-1] = None
            turn = 0
            Visit[s-1, turn] = True

            cur_r_dict = nx.to_dict_of_dicts(g)
            nodelist = list(g.nodes())

            ## Body
            for k in range(n):
                for v in nodelist:
                    if Visit[v-1, turn] == True:
                        Visit[v-1, turn] = False
                        for u in cur_r_dict[v].keys():

                            try:
                                w_vu = orig_dict[v][u]['weight']
                            except:
                                w_vu = -orig_dict[u][v]['weight']

                            if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                pi[k+1, u-1] = v
                                Visit[u-1, 1-turn] = True
                turn = 1 - turn

            ## Tail
            lam = np.inf
            for v in nodelist:
                if Visit[v-1, turn] == True:
                    M[v-1] = -np.inf
                    for k in range(n):
                        if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                            M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                            K[v-1] = k

                    if lam > M[v-1]:
                        lam = M[v-1]
                        v_star = v

            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            self.d_kv = d_kv
            if lam < 0:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

    def set_nmcc(self, G):

        glen = len(self.g.nodes())
        cycle_len = int(glen - self.K[self.v_star - 1])
        path_to_v = [self.v_star]
        next_v = self.v_star

        loopv = self.v_star
        for i in range(glen, 0, -1):

            next_v = int(self.pi[i,next_v-1])
            path_to_v.append(next_v)
            # if i!=glen and next_v == loopv:
            #     break
            if len(path_to_v) != len(set(path_to_v)):
                st = path_to_v[-1]
                ind = path_to_v.index(st)
                path_to_v = path_to_v[ind:]
                break

        if len(path_to_v) == len(set(path_to_v)):
            path_to_v = [self.v_star] + path_to_v

        self.nmcc = path_to_v[::-1]
        #############
        # print('----Cycle Found-----')
        # print(self.nmcc)


        ####
        # if (len(self.nmcc)-1) != cycle_len and self.lam_star<-1e-6:
        #     print('hooop length is different')
        #     pdb.set_trace()

        # nmcc_cost = 0
        # for i in range(len(self.nmcc)-1):
        #     u = self.nmcc[i]
        #     v = self.nmcc[i+1]
        #     if G.nxg.has_edge(u,v):
        #         nmcc_cost += G.nxg[u][v]['weight']
        #     else:
        #         nmcc_cost += -G.nxg[v][u]['weight']

        # if round(self.lam_star,2) != round(nmcc_cost/(len(self.nmcc)-1),2):
        #     print('hooop cost is different')
        #     pdb.set_trace()

        # print('=======GUPTA_MMC===========')
        # print('============================')
        # print('| Cycle: ', self.nmcc, ' | Cost: ', nmcc_cost/(len(self.nmcc)-1), ' | Mmc_true: ', self.lam_star, ' |')
        # print('| True cycle lenght: ', cycle_len)
        # print(' ')
        ######


        # self.check_nmcc(G.nxg,G.lam,self.g,cycle_len)


        ##########################
        # cycle_len = int(glen - self.K[self.v_star-1])
        # print('cycle length: ', cycle_len)
        # k = glen
        # P = np.zeros(glen+1)
        # P[k] = int(self.v_star)
        # while k>0:
        #     P[k-1] = int(self.pi[k,int(P[k])-1])
        #     k = k-1

        # W = np.zeros(glen+1)

        # for i in range(1,glen+1):
        #     mymin = np.inf
        #     for j in range(glen+1):
        #         mymin = min(mymin, self.d_kv[j,i-1] -j*self.lam_star)
        #     W[i-1] = mymin

        # k = 0
        # W[k] = 0
        # while k < glen:
        #     k = k+1
        #     if G.nxg.has_edge(P[k-1],P[k]):
        #         w = G.nxg[P[k-1]][P[k]]['weight']
        #     else:
        #         w = -G.nxg[P[k]][P[k-1]]['weight']
        #     W[k] = W[k-1] + w

        # k = 0
        # cycle = []
        # while k <= self.K[self.v_star-1]:
        #     if (P[k] == P[k+cycle_len]) and (round((W[k+cycle_len] - W[k])/cycle_len,3) == round(self.lam_star,3)):
        #         # if len(P[k:k+cycle_len+1]) -1 == len(set(P[k:k+cycle_len+1])):
        #         cycle = P[k:k+cycle_len+1]
        #         cycle = cycle.astype(int)
        #     k = k+1
        # cycle = [n for n in cycle]
        # try:
        #     nx = cycle[1:].index(cycle[0])
        # except:
        #     pdb.set_trace()

        # cycle = cycle[:nx+2]
        # while len(cycle) == 0:
        #     mycopy = self.M.copy()
        #     minind = np.argmin(mycopy)
        #     mycopy[minind] = np.inf
        #     v_star = np.argmin(mycopy)
        #     cycle_len = int(glen - self.K[v_star])
        #     k=0
        #     cycle = []
        #     while k < self.K[v_star]:
        #         if (P[k] == P[k+cycle_len]) and (round((W[k+cycle_len] - W[k])/cycle_len,3) == round(self.lam_star,3)):
        #             if len(P[k:k+cycle_len+1]) -1 == len(set(P[k:k+cycle_len+1])):
        #                 cycle = P[k:k+cycle_len+1]
        #                 cycle = cycle.astype(int)
        #         k = k+1

        # self.nmcc = cycle

    def DG_3_dict(self, G):
        g = self.nxg
        n = len(g.nodes())
        m = len(g.edges())
        orig_dict = nx.to_dict_of_dicts(G.nxg)

        if not nx.is_strongly_connected(g):
            scc = [c for c in sorted(nx.strongly_connected_components(g),key=len, reverse=True)]
            scc_count = len(scc)
            min_lam = []
            M_l = []
            K_l = []
            v_l =[]
            pi_l = []
            glen_l =[]
            curg_l = []
            dkv_l = []

            for i in range(scc_count):
                nodelist = scc[i]
                if len(nodelist) <= 1:
                    continue

                cur_g = self.nxg.copy()
                allnodes = np.arange(len(self.nxg.nodes)) + 1
                removals = list(set(allnodes)-set(nodelist))
                cur_g.remove_nodes_from(removals)
                s = self.FindSource(cur_g)
                nodelist = list(nodelist)

                len_cur_g = len(nodelist)
                d_kv = np.full((n+1,n), np.inf)
                pi = np.full((n+1, n), -1)
                M = np.zeros(n)
                K = np.zeros(n)

                valid = np.zeros((n+1, n))
                lastlevel = np.full(n, -1)
                d_kv[0, s-1] = 0
                valid[0, s-1] = -1
                lastlevel[s-1] = 0
                q = deque()
                q.append((0,s))
                v_star = None

                cur_r_dict = nx.to_dict_of_dicts(cur_g)
                ## Body
                k, v = q.popleft()
                while k < len_cur_g:
                    for u in cur_r_dict[v].keys():
                        if lastlevel[u-1] < k+1:
                            q.append((k+1,u))
                            valid[k+1, u-1] = lastlevel[u-1]
                            lastlevel[u-1] = k+1

                        try:
                            w_vu = orig_dict[v][u]['weight']
                        except:
                            w_vu = -orig_dict[u][v]['weight']

                        if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                            d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                            pi[k+1, u-1] = v
                    k,v = q.popleft()
                ## Tail
                lam = np.inf
                for v in cur_r_dict.keys():
                    if lastlevel[v-1] == n:
                        M[v-1] = -np.inf
                        k = valid[len_cur_g, v-1]
                        while k > -1:
                            if M[v-1] < (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k):
                                M[v-1] = (d_kv[len_cur_g, v-1] - d_kv[k, v-1])/float(len_cur_g-k)
                                K[v-1] = k
                            k = valid[k, v-1]
                        if lam > M[v-1]:
                            lam = M[v-1]
                            v_star = v

                min_lam.append(lam)
                M_l.append(M)
                K_l.append(K)
                v_l.append(v_star)
                pi_l.append(pi)
                glen_l.append(len_cur_g)
                curg_l.append(cur_g)
                dkv_l.append(d_kv)

            min_ind = np.argmin(min_lam)
            lam = min_lam[min_ind]
            M = M_l[min_ind]
            K = K_l[min_ind]
            v_star = v_l[min_ind]
            pi = pi_l[min_ind]
            glen = glen_l[min_ind]
            g = curg_l[min_ind]
            self.d_kv = dkv_l[min_ind]
            self.M = M
            self.K = K
            self.v_star = v_star
            self.pi = pi
            self.lam_star = lam
            self.g = g
            print('lam_star: ', self.lam_star)
            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []

        else:

            d_kv = np.full((n+1,n), np.inf)
            pi = np.full((n+1, n), -1)
            valid = np.zeros((n+1, n))
            M = np.zeros(n)
            K = np.zeros(n)

            lastlevel = np.full(n, -1)
            s = self.FindSource(g)
            d_kv[0, s-1] = 0
            valid[0, s-1] = -1
            lastlevel[s-1] = 0
            q = deque()
            q.append((0,s))

            cur_r_dict = nx.to_dict_of_dicts(g)
            ## Body
            k,v = q.popleft()
            while k < n:
                for u in cur_r_dict[v].keys():
                    if lastlevel[u-1] < k+1:
                        q.append((k+1,u))
                        valid[k+1, u-1] = lastlevel[u-1]
                        lastlevel[u-1] = k+1

                    try:
                        w_vu = orig_dict[v][u]['weight']
                    except:
                        w_vu = -orig_dict[u][v]['weight']

                    if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                        d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                        pi[k+1, u-1] = v
                k,v = q.popleft()

            ## Tail
            lam = np.inf
            for v in cur_r_dict.keys():
                if lastlevel[v-1] == n:
                    M[v-1] = -np.inf
                    k = int(valid[n, v-1])
                    while k > -1:
                        if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                            M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                            K[v-1] = k
                        k = int(valid[k, v-1])
                    if lam > M[v-1]:
                        lam = M[v-1]
                        self.v_star = v

            self.M = M
            self.K = K
            self.pi = pi
            self.lam_star = lam
            self.g = g
            self.d_kv = d_kv
            print('lam_star: ', self.lam_star)

            if lam < -1e-6:
                self.nmcc_exists = True
            else:
                self.nmcc_exists = False
                self.nmcc = []
