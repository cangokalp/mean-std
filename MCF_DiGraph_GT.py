import numpy as np
import networkx as nx
import pdb
from itertools import islice
import cProfile
import pstats

class MCF_DiGraph_GT:
        
    """
    Extend networkx DiGraph class for our stochastic arc cost purposes
    """
    def __init__(self, residual=False):
        self.nxg = nx.DiGraph()
        self.residual = residual
        self.lam = 0
        self.initial_flow = []
        self.nmcc = []
        self.v_star = None
        self.pi = None
        self.K = None
        self.M = None
        
    def set_lambda(self, lam=0):
        self.lam = lam
        # self.set_weights()
    
    def get_lambda(self):
        
        return self.lam
        
    def set_weights(self):
        for u, v, e in self.nxg.edges(data=True):
            self.nxg[u][v]['weight'] = self.nxg[u][v]['mu'] + 2*self.lam*self.nxg[u][v]['var']
            
    def find_feasible_flow(self):
        """Establish a feasible flow
        """
        # solve max flow problem
        G_MaxFlow = self.nxg.copy()
        G_MaxFlow.add_node('s')
        G_MaxFlow.add_node('t')
        for i, data in G_MaxFlow.nodes(data=True):
            if (i != 's') and (i != 't'):
                b_i = G_MaxFlow.node[i]['demand']
                if b_i > 0:
                    G_MaxFlow.add_edge('s', i, capacity=b_i)
                elif b_i < 0:
                    G_MaxFlow.add_edge(i, 't', capacity=-b_i)
        
        cost, flow = nx.maximum_flow(G_MaxFlow, 's', 't')
        
        # cost, flow = nx.capacity_scaling(self.nxg)
        print("The minimum cost is:", cost)

        # The data structure of flow is not consistent with dictionary data structure 
        # needed for printing the flow vector
        feasible_flow = {}
        for i in self.nxg.nodes(): 
            for j in flow[i].keys():
                if j != 't':
                    feasible_flow[i,j] = flow[i][j]
        
        
        self.initial_flow = feasible_flow
        return feasible_flow
    
    def set_flow_attributes(self, feasible_flow):
        """Set flow variables for the graph
        """
        nx.set_edge_attributes(self.nxg, feasible_flow, 'flow')
            
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

    def mu_cost(self):

        flow = np.array(list(nx.get_edge_attributes(self.nxg,'flow').values()))
        mu = np.array(list(nx.get_edge_attributes(self.nxg,'mu').values()))

        # cost = 0
        # for u,v,e in self.nxg.edges(data=True):
        #     cost += e.get('flow', 0) * e.get('mu', 0)
            
        return np.dot(flow, mu)

    def var_cost(self):

        flow = np.square(np.array(list(nx.get_edge_attributes(self.nxg,'flow').values())))
        var = np.array(list(nx.get_edge_attributes(self.nxg,'var').values()))

        # cost = 0
        # for u,v,e in self.nxg.edges(data=True):
        #     cost += np.square(e.get('flow', 0)) * e.get('var', 0)
            
        return np.dot(flow, var)
    
    def xi_cost(self):

        flow = np.square(np.array(list(nx.get_edge_attributes(self.nxg,'flow').values())))
        var = np.array(list(nx.get_edge_attributes(self.nxg,'var').values()))
        xi = np.array(list(nx.get_edge_attributes(self.nxg,'xi').values()))

        # cost = 0
        # for u,v,e in self.nxg.edges(data=True):
        #     cost += e.get('xi', 0) * e.get('var', 0) * e.get('flow', 0)
       
        return np.sum(np.multiply(np.multiply(xi,var), flow))

    def tot_cost(self):

        cost = 0
        for u,v,e in self.nxg.edges(data=True):
            cost += e.get('flow', 0) * e.get('mu', 0) + self.lam*np.square(e.get('flow', 0)) * e.get('var', 0)

        return cost

    def FindSource(self, g):
        for i in range(len(g.nodes)):
            if g.has_node(i+1):
                break
        return i+1
    
    def build_res_network(self, demand='demand', capacity='capacity', weight='weight', mu='mu', var='var', flow='flow', jw=False):
        """Build a residual network and initialize a zero flow.
        """
        if sum(self.nxg.nodes[u].get(demand, 0) for u in self.nxg) != 0:
            raise nx.NetworkXUnfeasible("Sum of the demands should be 0.")

        R = MCF_DiGraph(residual=True)
        R.nxg.add_nodes_from((u, {'excess': -self.nxg.nodes[u].get(demand, 0),
                              'potential': 0}) for u in self.nxg)

        inf = float('inf')
        # Detect selfloops with infinite capacities and negative weights.
        for u, v, e in nx.selfloop_edges(self.nxg, data=True):
            if e.get(weight, 0) < 0 and e.get(capacity, inf) == inf:
                raise nx.NetworkXUnbounded(
                    'Negative cost cycle of infinite capacity found. '
                    'Min cost flow may be unbounded below.')

        # Extract edges with positive capacities. Self loops excluded.
        if self.nxg.is_multigraph():
            edge_list = [(u, v, k, e)
                         for u, v, k, e in self.nxg.edges(data=True, keys=True)
                         if u != v and e.get(capacity, inf) > 0]
        else:

            edge_list = [(u, v, e) for u, v, e in self.nxg.edges(data=True)
                         if u != v and e.get(capacity, inf) > 0]

        inf = max(sum(abs(R.nxg.nodes[u]['excess']) for u in R.nxg), 2 * sum(e[capacity] for u, v, e in edge_list if capacity in e and e[capacity] != inf)) or 1

        for u, v, e in edge_list:
            # pdb.set_trace()
            r_vu = e.get(flow, 0)
            r_uv = round(e.get(capacity,0) - e.get(flow,0),8)

            mu_val = e.get(mu, 0)
            var_val = e.get(var, 0)
            flow_val = e.get(flow, 0)
            cap_val = e.get(capacity, 0)
            try:
                weight_val = e.get(weight,0)
            except:
                weight_val = mu_val + self.lam * var_val


            if not jw:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, weight=weight_val, capacity=cap_val, mu=mu_val, var=var_val, r=r_uv)

                if r_vu != 0:
                    R.nxg.add_edge(v, u, weight=-weight_val, capacity=cap_val, mu=-mu_val, var=-var_val, r=r_vu)
            else:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, r=r_uv, weight=weight_val)
                if r_vu != 0:
                    R.nxg.add_edge(v, u, r=r_vu, weight=-weight_val)

        R.nxg.graph['inf'] = inf
        return R

    def find_delta(self):
        
        min_r = np.inf
        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]
            rij = self.nxg[u][v]['r']
            min_r = min(min_r, rij)

        return min_r

    def find_xsi_star(self, G):
        nominator = 0
        denominator = 0
        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]

            if G.nxg.has_edge(u,v):
                flow = G.nxg[u][v]['flow']
            else:
                flow = G.nxg[v][u]['flow']
            
            nominator -= self.nxg[u][v]['mu'] + 2*G.lam*flow*self.nxg[u][v]['var']
            
            denominator += 2*G.lam*abs(self.nxg[u][v]['var'])

        return nominator/denominator

    def remove_no_r(self,u, v):
        if self.nxg.has_edge(u,v):
            if self.nxg[u][v]['r'] == 0:
                self.nxg.remove_edge(u,v)

    def augment_flow(self, min_r):
        for i in range(len(self.nmcc)-1):
            u = self.nmcc[i]
            v = self.nmcc[i+1]
            
            # weight = self.nxg[u][v]['weight']
            try:
                mu = self.nxg[u][v]['mu']
                var = self.nxg[u][v]['var']
            except:
                pass

            if self.nxg.has_edge(u, v):
                self.nxg[u][v]['r'] -= min_r
                
            self.remove_no_r(u, v)

            if self.nxg.has_edge(v, u):
                self.nxg[v][u]['r'] += min_r
            else:
                # self.nxg.add_edge(v, u, r=min_r, weight=-weight, mu=-mu, var=-var)
                self.nxg.add_edge(v, u, r=min_r, mu=-mu, var=-var)

            self.remove_no_r(v, u)

    def adjust_flow(self, g, nmcc, delta):

        for i in range(len(nmcc)-1):
            u = nmcc[i]
            v = nmcc[i+1]

            try:         
                g.ep['flow'][u,v] += delta
            except:
                g.ep['flow'][v,u] -= delta

    def nmcc_exists(self, G):
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
                n = g.num_vertices(ignore_filter=True)
                m = g.num_edges(ignore_filter=True)
                len_cur_g = cur_g.num_vertices(ignore_filter=True)
                vertices = cur_g.vertices()
                d_kv = np.ones((n+1, n))
                d_kv = np.inf * d_kv
                pi = -1*np.ones((n+1, n))
                Visit = np.ones((n, n))
                M = np.zeros(n)
                K = np.zeros(n)
                v_star = None

                for v in vertices:
                    Visit[v-1, 0] = False
                    Visit[v-1, 1] = False

                d_kv[0, s-1] = 0

                pi[0, s-1] = None
                turn = 0
                Visit[s-1, turn] = True
                ## Body
                for k in range(len_cur_g):
                    for v in vertices:
                        if Visit[v-1, turn] == True:
                            Visit[v-1, turn] = False
                            for u in cur_g.get_out_neighbors(v):


                                try:         
                                    flow = g.ep['flow'][v,u]
                                except:
                                    flow = g.ep['flow'][u,v]

                                w_vu= r.ep['mu'][v,u] + 2*g.lam*flow*r.ep['var'][v,u]
 
                                if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                    d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                    pi[k+1, u-1] = v
                                    Visit[u-1, 1-turn] = True
                    turn = 1 - turn
                ## Tail
                lam = np.inf
                for v in vertices:
                    if Visit[v-1, turn] == True:
                        M[v-1] = -np.inf
                        for k in range(len_cur_g):

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
            self.set_nmcc(glen, G)             
            return True if lam < -1e-4 else False

        else:       

            n = g.num_vertices(ignore_filter=True)
            m = g.num_edges(ignore_filter=True)
            d_kv = np.ones((n+1, n))
            d_kv = np.inf * d_kv
            pi = -1*np.ones((n+1, n))
            Visit = np.ones((n, n))
            M = np.zeros(n)
            K = np.zeros(n)
            v_star = None
            vertices = g.vertices()

            for v in vertices:
                Visit[v-1, 0] = False
                Visit[v-1, 1] = False

            s = self.FindSource(g)
            d_kv[0, s-1] = 0

            pi[0, s-1] = None
            turn = 0
            Visit[s-1, turn] = True

            ## Body
            for k in range(n):
                for v in vertices:
                    if Visit[v-1, turn] == True:
                        Visit[v-1, turn] = False
                        for u in g.get_out_neighbors(v):
                            if G.nxg.has_edge(v,u):
                                flow = G.nxg[v][u]['flow']
                            else:
                                flow = G.nxg[u][v]['flow']
 
                            w_vu= self.nxg[v][u]['mu'] + 2*G.lam*flow*self.nxg[v][u]['var']
 
                            if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                pi[k+1, u-1] = v
                                Visit[u-1, 1-turn] = True
                turn = 1 - turn

            ## Tail
            lam = np.inf
            for v in vertices:
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
            self.set_nmcc(len(g.nodes()), G)  
            return True if lam < -1e-4 else False

    def set_nmcc(self, glen, G):
        d_kv = self.d_kv
        cycle_len = int(glen - self.K[self.v_star - 1])
        path_to_v = []
        next_v = self.v_star

        if self.lam_star < -1e-4:
            j = 0
            for i in range(glen, 0, -1):
                next_v = int(self.pi[glen-j,next_v-1])
                if i == glen:
                    loopv = next_v
                path_to_v.append(next_v)
                if i!=glen and next_v == loopv:
                    break
                if len(path_to_v) != len(set(path_to_v)):
                    st = path_to_v[-1]
                    ind = path_to_v.index(st)
                    path_to_v = path_to_v[ind:]
                    break
                j += 1
        else:
            self.nmcc = []

        self.nmcc = path_to_v[::-1]
            
    def nmcc_exists_test(self, G):
            ## Head
            g = self.nxg
            
            #test#
            if not nx.is_strongly_connected(g):
                print('not strongly connected')
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
                    for node, data in self.nxg.edges():
                        if node not in nodelist:
                            cur_g.remove_node(node)

                    n = len(cur_g.nodes())
                    m = len(cur_g.edges())
                    d_kv = np.ones((n+1, n))
                    d_kv = np.inf * d_kv
                    pi = -1*np.ones((n+1, n))
                    Visit = np.ones((n, n))
                    M = np.zeros(n)
                    K = np.zeros(n)
                    v_star = None

                    for v in cur_g.nodes():
                        Visit[v-1, 0] = False
                        Visit[v-1, 1] = False

                    s = self.FindSource()
                    d_kv[0, s-1] = 0

                    pi[0, s-1] = None
                    turn = 0
                    Visit[s-1, turn] = True

                    ## Body
                    for k in range(n):
                        for v in cur_g.nodes():
                            if Visit[v-1, turn] == True:
                                Visit[v-1, turn] = False
                                for u in cur_g.neighbors(v):

                                    w_vu= self.nxg[v][u]['weight']
                                    
                                    if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                        d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                        pi[k+1, u-1] = v
                                        Visit[u-1, 1-turn] = True
                        turn = 1 - turn

                    ## Tail
                    lam = np.inf
                    for v in cur_g.nodes():
                        if Visit[v-1, turn] == True:
                            M[v-1] = -np.inf
                            for k in range(n):
                                if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k): 
                                    M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)  
                                    K[v-1] = k

                            if lam > M[v-1]: 
                                lam = M[v-1] 
                                v_star = v

                    min_lam.append(lam)
                    M_l.append(M)
                    K_l.append(K)
                    v_l.append(v_star)
                    pi_l.append(pi)
                    glen_l.append(len(cur_g.nodes()))
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
                self.set_nmcc_test(glen, G)            
                return lam
            else:       

                n = len(g.nodes())
                m = len(g.edges())
                d_kv = np.ones((n+1, n))
                d_kv = np.inf * d_kv
                pi = -1*np.ones((n+1, n))
                Visit = np.ones((n, n))
                M = np.zeros(n)
                K = np.zeros(n)
                v_star = None

                for v in g.nodes():
                    Visit[v-1, 0] = False
                    Visit[v-1, 1] = False

                s = self.FindSource()
                d_kv[0, s-1] = 0

                pi[0, s-1] = None
                turn = 0
                Visit[s-1, turn] = True

                ## Body
                # pdb.set_trace()
                for k in range(n):
                    for v in g.nodes():
                        if Visit[v-1, turn] == True:
                            Visit[v-1, turn] = False
                            for u in g.neighbors(v):
                                w_vu = self.nxg[v][u]['weight']
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
                self.set_nmcc_test(len(self.nxg.nodes()),G)   
                return lam

    def set_nmcc_test(self, glen, G):
        d_kv = self.d_kv
        cycle_len = int(glen - self.K[self.v_star - 1])
        path_to_v = []
        next_v = self.v_star
        print('cyle length: ', cycle_len)
        print('v_star: ', self.v_star)

        # pdb.set_trace()
        if self.lam_star < 0:
            j = 0
            for i in range(glen, 0, -1):
                next_v = int(self.pi[glen-j,next_v-1])
                if i == glen:
                    loopv = next_v
                path_to_v.append(next_v)
                if next_v == loopv:
                    break
                # if len(path_to_v) != len(set(path_to_v)):
                #     st = path_to_v[-1]
                #     pdb.set_trace()
                #     ind = path_to_v.index(st)
                #     path_to_v = path_to_v[ind:]
                #     break
                j += 1

            self.nmcc = path_to_v[::-1]
        else:
            self.nmcc = []


        print(self.nmcc)
        print('====')
