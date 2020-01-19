def bisect_solve_exp(**kwargs):
    print('-----BISECTION------')
    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)
    start = time.time()

    if kwargs['R'] != None:
        R = copy.deepcopy(kwargs['R'])
    else:
        R = G.build_res_network()

    lam_bar = kwargs['lam_bar']
    stop_tol = kwargs['stop_tol']
    subproblem_tol = kwargs['subproblem_tol_bs']
    precision_start_tol = kwargs['precision_start_tol_bs']
    muscalar = kwargs['mu_scalar']
    found = False
    iters = 0

    disc_tot = 0
    nmcc_time_tot = 0
    scc_time_tot = 0
    lams = []

    varcost = kwargs['start_varcost_bs']
    nmcc_time_start = kwargs['add_time_bs']
    lam_high = kwargs['lam_high_bs']
    lam = kwargs['startlam_bs']
    low = kwargs['lam_low']

    # sigmacost = np.sqrt(varcost)
    high = lam_high
    f = 100

    objs = []
    subproblem_times = []
    times = []
    G.set_muscalar(kwargs['mu_scalar'])

    # discount = time.time()
    # algo_obj = (lam_bar * sigmacost + G.mu_cost() *
    #             kwargs['mu_scalar']) * kwargs['obj_scalar']

    # objs.append(algo_obj)
    # disc_tot += time.time() - discount
    subproblem_times.append(nmcc_time_start)
    times.append(time.time() - start - disc_tot + nmcc_time_start)
    nfultol = 1e-2

    while not found:
        iters += 1
        mid = (high + low) / 2.0
        sub_start = time.time()
        G.set_lambda(lam=mid)

        if abs(f) <= kwargs['precision_start_tol_bs']:
            nfultol = 1e-4

        discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
            G, R, varcost, fully=True, nfullytol=nfultol, var_cost_ful=kwargs['varcostful'], vartop=kwargs['var_top'], mutop=kwargs['mu_top'], difftol=1e-2)
        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time

        subproblem_times.append(
            time.time() - sub_start - discount + nmcc_time - scc_time)

        sigmacost = np.sqrt(varcost)

        disc_me = time.time()
        algo_obj= (lam_bar * sigmacost + G.mu_cost())

        # algo_obj = (mid * varcost + G.mu_cost() *
        #             kwargs['mu_scalar']) * kwargs['obj_scalar']
        times.append(time.time() - start - disc_tot +
                     nmcc_time_tot - scc_time_tot + nmcc_time_start)

        objs.append(algo_obj)
        lams.append(mid)
# 
        # print('algo_obj: {}, f: {}, lam: {}'.format(algo_obj, f, mid))

        disc_tot += time.time() - disc_me
        f = mid - float(lam_bar) / float((2 * sigmacost))

        if abs(f) <= kwargs['stop_tol']:
            found = True
            break

        if np.sign(f) == np.sign(1):
            high = mid
        else:
            low = mid

        mid_prev = mid

    end = time.time()

    # cur_r_dict = nx.to_dict_of_dicts(G.nxg)
    # cc = [v for v in cur_r_dict.keys()]
    # list_ofed = []
    # endnode = len(G.nxg.nodes())
    # for c in cc:
    #     if endnode in cur_r_dict[c].keys():
    #         list_ofed.append(c)

    # flows = []
    # for e in list_ofed:
    #     flows.append(G.nxg[e][endnode]['flow'])

    # sflows = sum(np.array(flows))
    # print('totflow: ', sflows)
    # x = nx.get_edge_attributes(G.nxg, 'flow')

    return algo_obj, iters, np.array(subproblem_times), times, objs, lams, 5, G.mu_cost(), np.sqrt(G.var_cost())
# def read_log(logfilename):
#     f = open(logfilename + '_processed.txt', 'w')
#     ff = open(logfilename + '.txt', 'r+')
#     w_to_processed = False

#     for line in ff.readlines():
#         if line.find('ITE PFEAS') >= 0:
#             w_to_processed = True
#         if line.find('Optimizer terminated') >= 0:
#             w_to_processed = False
#         if w_to_processed:
#             f.write(line)
#     f.close()

#     log = pandas.read_table(logfilename + '_processed.txt',
#                             delim_whitespace=True, sep='\t', lineterminator='\n')
#     return log
# def update_weight(R, G, lam, delta=0):
#     for u, v, e in R.edges(data=True):
#         try:
#             w = (e['flow'] + delta) + 2 * lam * \
#                  G[u][v]['var'] * (e['flow'] + delta)
#         except:
#             w = (e['flow'] + delta) + 2 * lam * \
#                  G[v][u]['var'] * (e['flow'] + delta)

#         e['weight'] = w
# # __all__ = ['capacity_scaling']
# # from itertools import chain
# # from math import log
# # import networkx as nx
# # from networkx.utils import *
# def _detect_unboundedness(R):
#     """Detect infinite-capacity negative cycles.
#     """
#     s = generate_unique_node()
#     G = nx.DiGraph()
#     G.add_nodes_from(R)

#     # Value simulating infinity.
#     inf = R.graph['inf']
#     # True infinity.
#     f_inf = float('inf')
#     for u in R:
#         for v, e in R[u].items():
#             # Compute the minimum weight of infinite-capacity (u, v) edges.
#             w = f_inf
#             for k, e in e.items():
#                 if e['capacity'] == inf:
#                     w = min(w, e['weight'])
#             if w != f_inf:
#                 G.add_edge(u, v, weight=w)

#     if nx.negative_edge_cycle(G):
#         raise nx.NetworkXUnbounded(
#             'Negative cost cycle of infinite capacity found. '
#             'Min cost flow may be unbounded below.')
# def _build_residual_network(G, demand, capacity, weight):
#     """Build a residual network and initialize a zero flow.
#     """
#     if sum(G.node[u].get(demand, 0) for u in G) != 0:
#         raise nx.NetworkXUnfeasible("Sum of the demands should be 0.")

#     R = nx.MultiDiGraph()
#     R.add_nodes_from((u, {'excess': G.node[u].get(demand, 0),
#                           'potential': 0}) for u in G)

#     inf = float('inf')
#     # Detect selfloops with infinite capacities and negative weights.
#     for u, v, e in G.selfloop_edges(data=True):
#         if e.get(weight, 0) < 0 and e.get(capacity, inf) == inf:
#             raise nx.NetworkXUnbounded(
#                 'Negative cost cycle of infinite capacity found. '
#                 'Min cost flow may be unbounded below.')

#     # Extract edges with positive capacities. Self loops excluded.
#     if G.is_multigraph():
#         edge_list = [(u, v, k, e)
#                      for u, v, k, e in G.edges(data=True, keys=True)
#                      if u != v and e.get(capacity, inf) > 0]
#     else:
#         edge_list = [(u, v, 0, e) for u, v, e in G.edges(data=True)
#                      if u != v and e.get(capacity, inf) > 0]
#     # Simulate infinity with the larger of the sum of absolute node imbalances
#     # the sum of finite edge capacities or any positive value if both sums are
#     # zero. This allows the infinite-capacity edges to be distinguished for
#     # unboundedness detection and directly participate in residual capacity
#     # calculation.
#     inf = max(sum(abs(R.node[u]['excess']) for u in R),
#               2 * sum(e[capacity] for u, v, k, e in edge_list
#                       if capacity in e and e[capacity] != inf)) or 1
#     for u, v, k, e in edge_list:
#         r = min(e.get(capacity, inf), inf)
#         w = e.get(weight, 0)
#         # Add both (u, v) and (v, u) into the residual network marked with the
#         # original key. (key[1] == True) indicates the (u, v) is in the
#         # original network.
#         R.add_edge(u, v, key=(k, True), capacity=r, weight=w, flow=0)
#         R.add_edge(v, u, key=(k, False), capacity=0, weight=-w, flow=0)

#     # Record the value simulating infinity.
#     R.graph['inf'] = inf

#     _detect_unboundedness(R)

#     return R
# def _build_flow_dict(G, R, capacity, weight):
#     """Build a flow dictionary from a residual network.
#     """
#     inf = float('inf')
#     flow_dict = {}
#     if G.is_multigraph():
#         for u in G:
#             flow_dict[u] = {}
#             for v, es in G[u].items():
#                 flow_dict[u][v] = dict(
#                     # Always saturate negative selfloops.
#                     (k, (0 if (u != v or e.get(capacity, inf) <= 0 or
#                                e.get(weight, 0) >= 0) else e[capacity]))
#                     for k, e in es.items())
#             for v, es in R[u].items():
#                 if v in flow_dict[u]:
#                     flow_dict[u][v].update((k[0], e['flow'])
#                                            for k, e in es.items()
#                                            if e['flow'] > 0)
#     else:
#         for u in G:
#             flow_dict[u] = dict(
#                 # Always saturate negative selfloops.
#                 (v, (0 if (u != v or e.get(capacity, inf) <= 0 or
#                            e.get(weight, 0) >= 0) else e[capacity]))
#                 for v, e in G[u].items())
#             flow_dict[u].update((v, e['flow']) for v, es in R[u].items()
#                                 for e in es.values() if e['flow'] > 0)
#     return flow_dict
# def capacity_scaling(G, demand='demand', capacity='capacity', weight='weight',heap=BinaryHeap, lam=0.7):

#     R = _build_residual_network(G, demand, capacity, weight)

#     inf = float('inf')
#     # Account cost of negative selfloops.
#     flow_cost = sum(
#         0 if e.get(capacity, inf) <= 0 or e.get(weight, 0) >= 0
#         else e[capacity] * e[weight]
#         for u, v, e in G.selfloop_edges(data=True))

#     # Determine the maxmimum edge capacity.
#     wmax = max(chain([-inf],
#                      (e['capacity'] for u, v, e in R.edges(data=True))))
#     if wmax == -inf:
#         # Residual network has no edges.
#         return flow_cost, _build_flow_dict(G, R, capacity, weight)

#     R_node = R.node
#     R_succ = R.succ

#     delta = 2 ** int(log(wmax, 2))
#     update_weight(R, G, lam)
#     while delta >= 0.0025:
#         # Determine the Δ-active nodes.
#         S = set()
#         T = set()
#         S_add = S.add
#         S_remove = S.remove
#         T_add = T.add
#         T_remove = T.remove
#         for u in R:
#             excess = R_node[u]['excess']
#             if excess >= delta:
#                 S_add(u)
#             elif excess <= -delta:
#                 T_add(u)
#         # Repeatedly augment flow from S to T along shortest paths until
#         # Δ-feasibility is achieved.
#         print(delta)
#         while S and T:
#             s = next(iter(S))
#             t = None
#             # Search for a shortest path in terms of reduce costs from s to
#             # any t in T in the Δ-residual network.
#             d = {}
#             pred = {s: None}
#             h = heap()
#             h_insert = h.insert
#             h_get = h.get
#             h_insert(s, 0)
#             while h:
#                 u, d_u = h.pop()
#                 d[u] = d_u
#                 if u in T:
#                     # Path found.
#                     t = u
#                     break
#                 p_u = R_node[u]['potential']
#                 for v, es in R_succ[u].items():
#                     if v in d:
#                         continue
#                     wmin = inf
#                     # Find the minimum-weighted (u, v) Δ-residual edge.
#                     for k, e in es.items():
#                         if e['capacity'] - e['flow'] >= delta:
#                             w = e['weight']
#                             if w < wmin:
#                                 wmin = w
#                                 kmin = k
#                                 emin = e
#                     if wmin == inf:
#                         continue
#                     # Update the distance label of v.
#                     d_v = d_u + wmin - p_u + R_node[v]['potential']
#                     if h_insert(v, d_v):
#                         pred[v] = (u, kmin, emin)
#             # pdb.set_trace()
#             if t is not None:
#                 # Augment Δ units of flow from s to t.
#                 while u != s:
#                     v = u
#                     u, k, e = pred[v]
#                     e['flow'] += delta
#                     R_succ[v][u][(k[0], not k[1])]['flow'] -= delta
#                 # Account node excess and deficit.
#                 R_node[s]['excess'] -= delta
#                 R_node[t]['excess'] += delta
#                 if R_node[s]['excess'] < delta:
#                     S_remove(s)
#                 if R_node[t]['excess'] > -delta:
#                     T_remove(t)
#                 # Update node potentials.
#                 d_t = d[t]
#                 for u, d_u in d.items():
#                     R_node[u]['potential'] -= d_u - d_t
#             else:
#                 # Path not found.
#                 S_remove(s)
#             update_weight(R, G, lam)

#         # Saturate Δ-residual edges with negative reduced costs to achieve
#         # Δ-optimality.
#         for u in R:
#             p_u = R_node[u]['potential']
#             for v, es in R_succ[u].items():
#                 for k, e in es.items():
#                     flow = e['capacity'] - e['flow']
#                     if e['weight'] - p_u + R_node[v]['potential'] < 0:
#                         flow = e['capacity'] - e['flow']
#                         if flow >= delta:
#                             e['flow'] += flow
#                             R_succ[v][u][(k[0], not k[1])]['flow'] -= flow
#                             R_node[u]['excess'] -= flow
#                             R_node[v]['excess'] += flow
#         update_weight(R, G, lam)

#         delta = delta / 2.0

#     # pdb.set_trace()
#     if any(R.node[u]['excess'] > 2 * delta for u in R):
#         raise nx.NetworkXUnfeasible('No flow satisfying all demands.')

#     # Calculate the flow cost.
#     for u in R:
#         for v, es in R_succ[u].items():
#             for e in es.values():
#                 flow = e['flow']
#                 if flow > 0:
#                     flow_cost += flow * e['weight']

#     return flow_cost, _build_flow_dict(G, R, capacity, weight)
# def push_flow(G, P, delta, g):
#     for i in range(len(P)-1):
#         u = P[i]
#         v = P[i+1]

#         if G.nxg.has_edge(u, v):
#             flow =  G.nxg[u][v]['flow']
#             G.nxg[u][v]['flow'] = flow + delta
#         else:
#             flow =  G.nxg[v][u]['flow']
#             G.nxg[v][u]['flow'] = flow - delta
#         g[u-1] -= delta
#         g[v-1] += delta

#     return g
# def ASSP(G, lam):
#     nx.set_node_attributes(G.nxg, 0, 'potential')
#     nx.set_edge_attributes(G.nxg, 0, 'flow')

#     theta = 0.5
#     rho = 1.0 / 3.0

#     m = len(G.nxg.edges())

#     mu = np.zeros(m)
#     var = np.zeros(m)
#     sigma = np.zeros(m)
#     cap = np.zeros(m)
#     d = np.zeros(len(G.nxg.nodes()))
#     rows = []
#     cols = []
#     values = []
#     i = 0
#     arc_dict = {}
#     for u, v, e in G.nxg.edges(data=True):
#         pdb.set_trace()
#         mu[i] = e.get('mu', 0)
#         sigma[i] = np.sqrt(e.get('var', 0))
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         arc_dict[i] = (u, v)
#         rows.append(u - 1)
#         cols.append(i)
#         values.append(1)
#         rows.append(v - 1)
#         cols.append(i)
#         values.append(-1)
#         i += 1

#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         d[i] = dat['demand']
#         i += 1


#     max_d = max(d)
#     max_d_vec = np.ones(len(cap)) * max_d
#     x_carrot = np.minimum(cap, max_d_vec)

#     eps_max = np.inf
#     eps = min(eps_max, max(1.0 / (len(G.nxg.nodes()) + 1),
#               rho * np.max(mu + var * x_carrot)))

#     # x = np.zeros(len(cap))
#     # for i,j,e in G.nxg.edges(data=True):
#     #     G.nxg[i][j]['flow'] = max(0, min(G.nxg[i][j]['capacity'], - (G.nxg[i][j]['mu'] - (G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'] )/(2*lam*G.nxg[u][v]['var']))))

#     g = np.zeros(len(G.nxg.nodes()))
#     for i in G.nxg.nodes():
#         outgoing = 0
#         for j, es in G.nxg._succ[i].items():
#             outgoing += G.nxg[i][j]['flow']
#         incoming = 0
#         for j, es in G.nxg._pred[i].items():
#             incoming += G.nxg[j][i]['flow']

#         g[i-1] = outgoing - incoming + G.nxg._node[i]['demand']

#     surp = g[g > 0]
#     pos_surp = [ list(g).index(i) for i in list(surp)]
#     s = pos_surp[0] + 1
#     P = [s]
#     while len(g[g > 0]) > 0:

#         print('-------')
#         print(eps)
#         print(P)
#         print(g)
#         '''step 1'''
#         i = P[-1]
#         step2time = True

#         for j, es in G.nxg._succ[i].items():
#             if G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'] > G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + theta * eps:
#                 step2time = False
#                 status = 'active'
#                 addme = j
#                 break

#         if step2time:
#             for j, es in G.nxg._pred[i].items():
#                 if G.nxg._node[j]['potential'] - G.nxg._node[i]['potential'] < G.nxg[j][i]['mu'] + 2*lam*G.nxg[j][i]['var']*G.nxg[j][i]['flow'] - theta * eps:
#                     step2time = False
#                     status = 'inactive'
#                     addme = j
#                     break

#         if not step2time:
#             P.append(addme)
#             if len(P) == len(list(set(P))):
#                 '''step 3a: push flow'''
#                 if g[addme-1] >= 0:
#                     continue
#                 else:

#                     fmarginP = np.inf
#                     for k in range(len(P)):
#                         if k == len(P)-1:
#                             break
#                         u = P[k]
#                         v = P[k+1]

#                         try:
#                             var = 2*lam*G.nxg[u][v]['var']
#                             deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] - (G.nxg[u][v]['mu'] + var*G.nxg[u][v]['flow']))/(var)
#                             # deltaij = min(deltaij, G.nxg[u][v]['capacity'] - G.nxg[u][v]['flow'])
                            
#                         except:
#                             var = 2*lam*G.nxg[v][u]['var']
#                             deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] + (G.nxg[v][u]['mu'] + var*G.nxg[v][u]['flow']))/(var)
#                             # deltaij = min(deltaij, G.nxg[v][u]['flow'])
                      
#                         if deltaij <= fmarginP:
#                             fmarginP = deltaij

#                     delta = min(fmarginP, g[s-1], -g[addme-1])
#                     g = push_flow(G, P, delta, g)
#                     surp = g[g > 0]
#                     pos_surp = [ list(g).index(o) for o in list(surp)]
#                     s = pos_surp[0] + 1
#                     P = [s]

#             else:
#                 '''step 3b: augment the flow on the cycle'''
#                 cs = set([x for x in P if P.count(x) > 1])
#                 C = P[P.index(cs.pop()):]
#                 fmarginC = np.inf

#                 for k in range(len(C)):
#                     if k == len(C)-1:
#                         break
#                     u = C[k]
#                     v = C[k+1]

#                     try:
#                         var = 2*lam*G.nxg[u][v]['var']
#                         deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] - (G.nxg[u][v]['mu'] + var*G.nxg[u][v]['flow']))/(var)
#                         deltaij = min(deltaij, G.nxg[u][v]['capacity'] - G.nxg[u][v]['flow'])
                        
#                     except:
#                         var = 2*lam*G.nxg[v][u]['var']
#                         deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] + (G.nxg[v][u]['mu'] + var*G.nxg[v][u]['flow']))/(var)
#                         deltaij = min(deltaij, G.nxg[v][u]['flow'])

#                     if deltaij <= fmarginC:
#                         fmarginC = deltaij

#                 g = push_flow(G, P, fmarginC, g)
#                 surp = g[g > 0]
#                 pos_surp = [ list(g).index(o) for o in list(surp)]
#                 s = pos_surp[0] + 1
#                 P = [s]

#         else:
#             '''step 2: price rise on i'''
#             minpi = np.inf
#             for j, es in G.nxg._succ[i].items():
#                 curpi = G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + eps + G.nxg._node[j]['potential']
#                 if curpi <= minpi:
#                     minpi = curpi
#             for j, es in G.nxg._pred[i].items():
#                 curpi = -G.nxg[j][i]['mu'] - 2*lam*G.nxg[j][i]['var']*G.nxg[j][i]['flow'] + eps + G.nxg._node[j]['potential']
#                 if curpi <= minpi:
#                     minpi = curpi

#             # delta_1 = np.inf
#             # for j, es in G.nxg._succ[i].items():
#             #     curpi = G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + eps - (G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'])
#             #     if curpi <= delta_1:
#             #         delta_1 = curpi
#             # delta_2 = np.inf
#             # for j, es in G.nxg._succ[i].items():
#             #     curpi = -(G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow']) + eps + (G.nxg._node[j]['potential'] - G.nxg._node[i]['potential'])
#             #     if curpi <= delta_2:
#             #         delta_2 = curpi

#             # delta_3 = np.inf
#             # minpi = min(delta_1, delta_2, delta_3)

#             G.nxg._node[i]['potential'] = minpi
#             if i != s:
#                 P.remove(i)

#     pdb.set_trace()
#     print('wait here')
# def dcinv(arc, t, lam=0.7, rev=False):
#     if not rev:
#         return (t - arc['mu'])/(2.0*lam*arc['var'])
# def tseng_ASSG(G):
#     nx.set_node_attributes(G.nxg, 0, 'potential')
#     nx.set_edge_attributes(G.nxg, 0, 'flow')

#     theta = 0.5
#     rho = 1.0 / 3.0

#     startn = {}
#     endn = {}
#     fin = {}
#     fout = {}
#     nxtin = {}
#     nxtou = {}
#     finalin = {}
#     finalou = {}

#     m = len(G.nxg.edges())
#     n = len(G.nxg.nodes())

#     for node in range(1,n+1):
#         fin[node] = 0
#         fout[node] = 0
#         finalin[node] = 0
#         finalou[node] = 0

#     for s,e,arc in G.nxg.edges(data=True):
#         if fout[s] != 0:
#             nxtou[finalou[s]] = arc
#         else:
#             fout[s] = (arc, e)
#         if fin[e] != 0:
#             nxtin[finalin[e]] = arc
#         else:
#             fin[e]=arc
#         finalou[s]=arc
#         finalin[e]=arc
#         nxtin[arc]=0
#         nxtou[arc]=0

#     for u,v,e in G.nxg.edges(data=True):
#         qfpushf[e] = (0, u, v, e)
#         qfpushb[e] = (0, u, v, e)
#         p[e] = 0

#     mu = np.zeros(m)
#     var = np.zeros(m)
#     sigma = np.zeros(m)
#     cap = np.zeros(m)
#     d = np.zeros(len(G.nxg.nodes()))
#     rows = []
#     cols = []
#     values = []
#     i = 0
#     arc_dict = {}
#     for u, v, e in G.nxg.edges(data=True):
#         mu[i] = e.get('mu', 0)
#         sigma[i] = np.sqrt(e.get('var', 0))
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         arc_dict[i] = (u, v)
#         rows.append(u - 1)
#         cols.append(i)
#         values.append(1)
#         rows.append(v - 1)
#         cols.append(i)
#         values.append(-1)
#         i += 1

#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         d[i] = dat['demand']
#         i += 1

#     max_d = max(d)
#     max_d_vec = np.ones(len(cap)) * max_d
#     x_carrot = np.minimum(cap, max_d_vec)

#     eps_max = np.inf
#     eps = min(eps_max, max(1.0 / (len(G.nxg.nodes()) + 1),
#               rho * np.max(mu + var * x_carrot)))

#     DFCT1 = {}
#     PRDCSR = {}
#     g = np.zeros(len(G.nxg.nodes()))
#     for i in G.nxg.nodes():
#         outgoing = 0
#         for j, es in G.nxg._succ[i].items():
#             outgoing += G.nxg[i][j]['flow']
#         incoming = 0
#         for j, es in G.nxg._pred[i].items():
#             incoming += G.nxg[j][i]['flow']

#         g[i-1] = outgoing - incoming + G.nxg._node[i]['demand']
#         DFCT1[i] = -g[i-1]
#         PRDCSR[i] = 0
#         SURPLUS[u] = g[i-1]

#     LARGE=2000000000
#     REF_THR = e-10
#     MAXPRICE = LARGE
#     TOTALITER = 0
#     niter = 0


#     #reduce arc capacities
#     for node in G.nxg.nodes():
#         CAPOUT = 0
#         arc = fout[node]
#         if arc > 0:
#             CAPOUT = min(LARGE - CAPOUT - arc['capacity'])
#             arc = nxtou[arc]
#         CAPOUT = min(LARGE - CAPOUT - SURPLUS[i])
#         if CAPOUT < 0:
#             print('infeasible')
#         CAPIN = 0
#         arc = fin[node]
#         if arc > 0:
#             if arc['capacity'] > CAPOUT:
#                 arc['capacity'] = CAPOUT
#             CAPIN = min(LARGE, CAPIN+arc['capacity'])
#             arc = ntxtin[arc]
#         CAPIN = min(LARGE, CAPIN + SURPLUS[node])
#         if CAPIN < 0:
#             print('infeasible')
#         arc = fout[node]
#         if arc > 0:
#             if arc['capacity'] > CAPIN:
#                 arc['capacity'] = CAPIN
#             arc = nxtou[arc]

#     #set arc flows to satisfy e-cs
#     eps2 = eps/2.0
#     eps3 = eps2/2.0

#     # for u,v,arc in G.nxg.edges(data=True):
#     #     start = u
#     #     end = v
#     #     PSTART = P[start]
#     #     PEND = P[end]
#     ###########

#     for i in G.nxg.nodes():
#         PRDCSR[i] = 0

#     NXTQUEUE = {}
#     for i in range(1, n):
#         NXTQUEUE[i] = i+1

#     NXTQUEUE[n]=1
#     NODE = 1
#     PREVNODE = NODE
#     LASTQUEUE = n

#     TOTALITER += TOTALITER
#     PRDCSR[NODE] = m + 1
#     i = NODE

#     qarcf = qfpushf[i]
#     qarcb = qfpushb[i]
#     price = p[i]

#     ###Forward arcs
#     if qarcf > 0:
#         j = qarcf[2]
#         if price - p[j] - (qarcf[3]['mu'] + 2*lam*qarcf[3]['var']*qarcf[3]['flow']) > eps2:
#             resid = min(qarcf['capacity'], dcinv(quarcf, price-p[j], lam=lam)) - qarcf[3]['flow']
#             if resid > REF_THR:
#                 qfpushf[i] = qarcf
#                 if -SURPLUS[j] > REF_THR:
#                     path_push1(j, qarcf[0])
#                     # go to 620
#                 else:
#                     if PRDCSR[j] == 0:
#                         PRDCSR[j] = qarcf
#                         i = j
#                         #go to 110
#                     else:
#                         cycle_push1(i,j,qarcf[0])
#                         #go to 620
#             else:
#                 qarcf = qnxtpushf[qarcf]
#         else:
#             qarcf = qnxtpushf[qarcf]
#         #go to 122
#     qfpushf[i] = qarcf

#     ###backward arcs
#     if qarcb > 0:
#         j = qarcb[2]
#         if price - p[j] - (qarcb[3]['mu'] + 2*lam*qarcb[3]['var']*qarcb[3]['flow']) > eps2:
#             resid = qarcb[3]['flow']  - max(0, dcinv(qarcb, p[j] - price, lam=lam)) 
#             if resid > REF_THR:
#                 qfpushb[i] = qarcb
#                 if -SURPLUS[j] > REF_THR:
#                     path_push1(j, -qarcb[0])
#                     # go to 620
#                 else:
#                     if PRDCSR[j] == 0:
#                         PRDCSR[j] = -qarcb
#                         i = j
#                         #go to 110
#                     else:
#                         cycle_push1(i,j,-qarcb[0])
#                         #go to 620
#             else:
#                 qarcb = qnxtpushf[qarcb]
#         else:
#             qarcb = qnxtpushf[qarcb]
#         #go to 123
#     qfpushb[i] = qarcb

#     if SURPLUS[i] > 0:
#         print(i, qarcf, qarcb, price)

#     if qarcb[0] + qarcf[0] == 0:
#         refprice = price
#         price = LARGE
#         niter = niter + 1
#         arc = fout[i][0]
#         if arc > 0:
#             up = arc['capacity']
#             if (arc['flow'] - up) < -REF_THR:
#                 xp = p[fout[i][1]] + (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -eps:
#                     price = xp + eps

#                     qarcf = arc
#                     qnxtpushf[arc] = 0
                    
#                 else:
#                     if (xp - price) < -eps2:
#                         if qarcf == 0:
#                             qarcf = arc
#                             qnxtpush[arc] = 0
#                         else:
#                             qnxtpushf[arc] = qarcf
#                             qarcf = arc
#             arc = nxtou[arc]
#             #goto 111

#         arc = fin[i][0]
#         if arc > 0:
#             if arc['flow'] > REF_THR:
#                 xp = p[fin[i][1]] - (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -eps:
#                     price = xp + eps
#                     qarcf = 0
#                     qarcb = arc
#                     qnxtpushb[arc] = 0
#                 else:
#                     if (xp-price) < -eps2:
#                         if qarcb == 0:
#                             qarcb = arc
#                             qnxtpushb[qarcb] = 0
#                         else:
#                             qnxtpushb[arc] = qarcb
#                             qarcb = arc

#             arc = nxtin[arc]
#             # goto 112


#         if price > LARGE:
#             qarcf = 0
#             qarcb = 0
#             arc = fout[i][0]
#             if arc > 0:
#                 xp = p[fout[i][1]] + (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -REF_THR:
#                     price = xp
#                 else:
#                     if xp > LARGE:
#                         print('price exceeds the range')
#                 arc = nxtou[arc]
#                 #goto 211
#             arc = fin[i]

#             if arc > 0:
#                 xp = p[fout[i][1]] - (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -REF_THR:
#                     price = xp
#                 else:
#                     if xp > LARGE:
#                         print('price exceeded the range')
#                 arc = nxtin[arc]
#                 #go to 212

#             if price > int(LARGE/10)
#                 price = int(LARGE/10)

#         p[i] = price
#         qfpushf[i] = qarcf
#         qfpushb[i] = qarcb


#         arc = PRDCSR[i]
#         if arc != m+1:
#             PRDCSR[i] = 0
#             if arc > 0:
#                 i = 0#startnode of arc
#             else:
#                 i = 0#endnode of -arc

#         # goto 110


#     if -SURPLUS[NODE] > THRESHSURPLUS:
#         TOTALITER +=1
#         PRDCSR[NODE] = m+1
#         i = NODE
# def subgrad(G, lam):
#     n = len(G.nxg.nodes())
#     m = len(G.nxg.edges())


#     var = np.zeros(m)
#     cap = np.zeros(m)
#     mu = np.zeros(m)

#     arc_dict = {}
#     A = np.zeros((n,m))
#     i = 0
#     for u, v, e in G.nxg.edges(data=True):
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         mu[i] = e.get('mu', 0)
#         A[u-1, i] = 1
#         A[v-1, i] = -1
#         arc_dict[i] = (u,v)
#         i += 1

#     b = np.zeros(n)
#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         b[i] = dat['demand']
#         i += 1

#     sumvar = 2*lam*var*np.eye(m)
#     sumvarinv = np.linalg.inv(sumvar)
#     H = np.dot(np.dot(A, sumvarinv),A.T)
#     Hinv = np.linalg.inv(H)

#     print('starting')
#     st = time.time()
#     p = np.zeros(n)
#     alpha = 1
#     g_p = np.zeros(n)
#     g_p[0] = 100
#     g_p[-1] = -100
#     close = False
#     iters = 1
#     while not close:
#         iters += 1

#         for i in G.nxg.nodes():
#             ij_term = 0
#             ji_term = 0
#             for j in np.where(A[i-1,:]==-1)[0]:
#                 idx = j
#                 u = arc_dict[j][0]
#             # for j, es in G.nxg._pred[i].items():
#                 # x_ji = max(0, min(G.nxg[j][i]['capacity'],(p[j-1] - p[i-1] - G.nxg[j][i]['mu'])/(2*lam*G.nxg[j][i]['var'])))
#                 x_ji = max(0, min(cap[idx],(p[u-1] - p[i-1] - mu[idx])/(2*lam*var[idx])))
#                 ji_term += x_ji
            
#             # for j, es in G.nxg._succ[i].items():
#             for j in np.where(A[i-1,:]==1)[0]:
#                 idx = j
#                 v = arc_dict[j][1]
#                 x_ij = max(0, min(cap[idx],(p[i-1] - p[v-1] - mu[idx])/(2*lam*var[idx])))
#                 ij_term += x_ij
            
#             g_p[i-1] = ji_term - ij_term + b[i-1]
#         p = p -alpha*np.dot(-Hinv,g_p)
#         if np.linalg.norm(g_p) <= 0.5:
#             close = True

#     print('time: ', time.time()-st)
#     for i,j,e in G.nxg.edges(data=True):
#             G.nxg[i][j]['flow'] = max(0, min(G.nxg[i][j]['capacity'],(p[i-1] - p[j-1] - G.nxg[i][j]['mu'])/(2*lam*G.nxg[i][j]['var'])))
#     tot_cost = G.mu_cost() + lam*G.var_cost()
#     print('obj: ', tot_cost)
#     pdb.set_trace()
