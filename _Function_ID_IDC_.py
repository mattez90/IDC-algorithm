# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 20:28:39 2021

@author: matte
"""
import networkx as nx
import numpy as np
import sys

#/---------------------------------------------------------------------------\
#|                          Support functions                                |
#\---------------------------------------------------------------------------/

def condset(var,orderlist):
    '''
    The function return the conditional set of the input variable based on 
    the topological order passed with variable orderlist. For example if 
    order list is ["A","B","E","G","Z"] the conditional set of variable "E"
    will be ["A","B"]. 
    To do: based on d-separation criteria the conditional set could be 
           simplify considering the Markov-Blanket set of the variable
           
    Parameters
    ----------
    var       : array that contains the name of the variable for which 
                we want to compute the conditional set
          
    orderlist : array that contains the topological ordering 

    Returns
    -------
    condset : array with conditional set for the variable of interest
    
    Author: Matteo Zoro 2020-02

    '''

    pos = int(np.where(orderlist==var)[0])

    if pos==0:
        condset = np.array([])
    else:
        # plus one in order to get the node itself
        condset = orderlist[0:pos]
        

    return condset   


def toposort(varsarray, order):
    return np.array([orderi for orderi in order if orderi in varsarray])


def removeedge(x,G,hiddvar):
    '''
    This function remove the nodes (and all the incoming outgoing edges from 
    the graph.  Since could be possible thathe nodes that we want to remove 
    have incoming edges from hidden variable the latter have to be remove from 
    the graph  as well (F.i if the graph is X<- U -> Y
    and X-> Y with U hidden node removing node X leads to the exclusion also 
    of the node U from the graph)

    Parameters
    ----------
    x       : array of nodes that we want to remove
    G       : networkx object that represent the graph
    hiddvar : array of hidden variable present in the graph

    Returns
    -------
    G_del : networkx object that represent the graph after the nodes are deleted
    
    Author: Matteo Zoro 2020-02

    '''
 
 
    # We have to remove all the bi-derected edge connected to the edge that
    # is removed. Note that using networkx in order to find a neighbors I
    # have combine successors and precedessors functions
    G_del = G.copy()
 
    for edge in x:

        neighbors = list(G.predecessors(edge)) + list(G.successors(edge))
        
         #Get the hidden edge
        delete = [u for u in hiddvar.tolist() if u in neighbors] + [edge]
       
        #Remove node and hidden edge connected to the removed edge
        G_del.remove_nodes_from(delete)
 
    return G_del

def incomingedgedel(nodes,G):
    '''
    This function remove all the outcoming edge of a list of nodes in the graph
    this opetarions mimic the intervention in the system that leads tu a 
    mutilated netwrok. This is implemented retriving the adjacency matrix from
    the input graph and then force to zero all values i-th column that 
    corresponds to the i-th node in the matrix. Note that in the adjacency 
    matrix [i,j]=0 , [j,i]=1  means j->i....putting [j,i]=0 then i j will be 
    separated

    Parameters
    ----------
    nodes : array containing nodes for which we want to remove the incoming edges
    G     : networkx object that rapresent the graph

    Returns
    -------
    G_ : networkx object that rapresent the graph after the incoming edges
         of interest nodes are removed

    Author: Matteo Zoro 2020-02
    '''
   
  
    #Retreive the adjency matrix
    adjmtx = nx.to_numpy_matrix(G)
  
    for nodei in nodes.tolist():
        adjmtx[:,list(nx.nodes(G)).index(nodei)]=0
       
    #Convert again to a graph object
    G_bar_x = nx.from_numpy_matrix(adjmtx,create_using=nx.DiGraph)
   
    #Map edge according to the list in the previous graph
    G_ = nx.relabel_nodes(G_bar_x, {i:j for i,j in enumerate(list(nx.nodes(G)))})
      
    return G_

def outcomeedgedel(nodes,G):
    '''
    This function remove all the outcoming edge of a list of nodes in the graph
    this opetarions mimic the intervention in the system that leads tu a 
    mutilated netwrok. This is implemented retriving the adjacency matrix from
    the input graph and then force to zero all values i-th row that 
    corresponds to the i-th node in the matrix. Note that in the adjacency 
    matrix [i,j]=0 , [j,i]=1  means j->i....putting [j,i]=0 then i j will be 
    separated

    Parameters
    ----------
    nodes : array containing nodes for which we want to remove the incoming edges
    G     : networkx object that rapresent the graph

    Returns
    -------
    G_ : networkx object that rapresent the graph after the incoming edges
         of interest nodes are removed

    Author: Matteo Zoro 2020-02
    '''
   
  
    #Retreive the adjency matrix
    adjmtx = nx.to_numpy_matrix(G)
  
    for nodei in nodes.tolist():
        adjmtx[list(nx.nodes(G)).index(nodei),:]=0
       
    #Convert again to a graph object
    G_bar_x = nx.from_numpy_matrix(adjmtx,create_using=nx.DiGraph)
   
    #Map edge according to the list in the previous graph
    G_ = nx.relabel_nodes(G_bar_x, {i:j for i,j in enumerate(list(nx.nodes(G)))})
      
    return G_


def ispresent(x,graph):
    return np.array([node for node in x.tolist() if node  in list(nx.nodes(graph))])


def C_components(graph,hiddvar,orderlist):
    '''
    This function compute the C-component of the graph as defined in Tian 
    2002 "A General Identification Condition for Causal Effects"
    (https://ftp.cs.ucla.edu/pub/stat_ser/R290-A.pdf). Consider the case of 
    presence of hidden node: X <- U -> Y; Bi-directed paths is 
    define as follow X <--> Y. In order to compute the C-componet I create 
    an adjacency matrix where I trace only the bi-directed edge leavingn the 
    node not connected to a b-directed edge "alone". If the graph 
    constructed by only bi-directed edge form a spanning tree then graph is
    a C-component. Note that the hidden node can't have incoming edge

    Parameters
    ----------
    graph : netwrokx object

    Returns
    -------
    list with connected component of the graph

    '''
 
    adj_mtx  = np.zeros((len(list(nx.nodes(graph))),len(list(nx.nodes(graph)))))
    listnode = list(nx.nodes(graph))
  
 
    #Trace all the bi-derected path starting from the empty matrix
    for u in hiddvar:
        if u in listnode:
            neighbors = list(graph.neighbors(u))
            adj_mtx[listnode.index(neighbors[0]),listnode.index(neighbors[1])] = adj_mtx[listnode.index(neighbors[0]),listnode.index(neighbors[1])] + 1
            adj_mtx[listnode.index(neighbors[1]),listnode.index(neighbors[0])] = adj_mtx[listnode.index(neighbors[1]),listnode.index(neighbors[0])] + 1
 
    # Delete nn-observable node from list of nodes and row and column from the
    # adjacency matrix
    adj_mtx2 = np.delete(np.delete(adj_mtx,[listnode.index(u) for u in hiddvar if u in listnode],axis=1),[listnode.index(u) for u in hiddvar if u in listnode],axis=0)
 
    # Delete non observable component
    for u in hiddvar:
         if u in listnode:
 
            listnode.remove(u)
      
    Gsub = nx.relabel_nodes(nx.from_numpy_matrix(adj_mtx2), {i:j for i, j in enumerate(listnode)})
    
    # I have to guarantee the topological order inside each C-component and 
    # among different C-component
    
    #Get C-component
    C_comp = [np.array(list(component_i)) for component_i in nx.connected_components(Gsub)]
    
    #Order the c-component
    C_topo = []
    posmax = []
    for ii in range(0,len(C_comp)):
        posmax = [orderlist.tolist().index(Sub_C_i) for Sub_C_i in C_comp[ii]]
        C_topo.append(max(posmax))
    
    C_topo_sort = C_topo.copy()
    C_topo_sort.sort()
    
    newind = [C_topo_sort.index(ii) for ii in C_topo]
    C_comp = [C_comp[idx] for idx in newind]
    
    #Order element inside C-component
    C_comp = [toposort(C_comp[ii], orderlist) for ii in range(0,len(C_comp))]
    
    return C_comp

def probabilitystr():
      return{"01_var" : np.array([]), "02_cond" : np.array([]), "03_sumvar" : np.array([]),"05_fraction" : False, "07_sum" : False, "08_prod" : False, "09_branch" : [] ,"04_numerator" : [], "06_denominator" : []}

def an(x,G,orderlist):
    '''
    Function return the ancestors of a interest node. Note that our definition
    of acenstor of a node include the node itself

    Parameters
    ----------
    x         : array containing nodes for which compute the ancestors
    G         : networkx object
    orderlist : array with topological order

    Returns
    -------
    anc : list of ancestor ordered according to the topological order 

     Author: Matteo Zoro 2020-02
    '''

    anc = []
    for xx in x.tolist():
        anc.extend(list(nx.ancestors(G,xx)))
        anc.extend(list(x))

    anc = np.array([ordi for ordi in orderlist if ordi in anc])
    
   
    return np.array(anc)



#/---------------------------------------------------------------------------\
#|                          Main  function Level - 0                         |
#\---------------------------------------------------------------------------/

def ID(y, x, P, G, orderlist, hiddvar,**kwargs):
    '''
    Function implement the ID algorihtm as presented in Pearl 2006 
    "Identification of Joint Interventional Distributionsin Recursive 
     Semi-Markovian Causal Models" 
    (https://www.aaai.org/Papers/AAAI/2006/AAAI06-191.pdf)
    
    
    Parameters
    ----------
    y         : array of viarables on which measures the effect of intervention
    x         : array of variables on which the intervention (do) is done
    P         : probability structure defined as the output of function 
                probabilitystr
    G         : netwrokx object
    orderlist : topological order
    hiddvar   : set of hidden (confonuders) nodes  
    **kwargs  : dictionary with following jey
                debug : if True message are print

    Raises
    ------
    ValueError: if the causal effect is not identifiable then the function 
                raise an error

    Returns
    -------
    probability structure that define the identification of causal effect
    '''

    P_out  = P.copy()
    if ("debug" in kwargs):
        debug  = kwargs["debug"]
    else :
        debug = False
    if ("debug" in kwargs):
        debug2  = kwargs["debug2"]
    else :
        debug2 = False


    # Line 0 (preliminary): i) Get observed (current) graph
    #                       ii) Get the list of (currente) node
    if len(ispresent(hiddvar,G)) > 0:
        G_obs = removeedge(ispresent(hiddvar,G),G,ispresent(hiddvar,G))
        v     = np.array(list(nx.nodes(G_obs)))
        v     = np.array([ordi for ordi in orderlist if ordi in v])
    else:
        G_obs = G.copy()
        v     = np.array(list(nx.nodes(G_obs)))
        v     = np.array([ordi for ordi in orderlist if ordi in v])
        
    if debug2:
        print(" ")
        print("--------- Graph information ------------------")
        print("> The edge of the current graph are:")
        print(">>>" + str(list(nx.edges(G))) + ";")
        print(" ")
        print("> The edge of the current observed graph are:")
        print(">>>" + str(list(nx.edges(G_obs))) + ";")
        print("---------------------------------------------")
        print(" ")
    
    #Line 1 of the algorihtm
    if len(x)==0:
        if debug:
            print(">> Enter in Line 1 of the algorithm")
        
        #Correct cases in Line 6 and 7
        if P_out["08_prod"]==True or P_out["05_fraction"]==True:
            P_out["03_sumvar"] = toposort(np.union1d(np.setdiff1d(v,y),P_out["03_sumvar"]), orderlist)  #list((set(v)-set(y)) | set(P_out["03_sumvar"]))
            P_out["07_sum"]    = True
        else:
            if len(v)==1: #==0 Means just P(var)
                 P_out["07_sum"]    = True
                 P_out["03_sumvar"] = toposort(np.union1d(np.setdiff1d(v,y),P_out["03_sumvar"]), orderlist)   #list((set(v)-set(y)) | set(P_out["03_sumvar"]))
                 
            P_out["01_var"]    = y

            
        return P_out

    #>Line 2 of the algorithm
    # Get the ancestors of y (!! note that the process is done in the observed
    # graph)

    anc = an(y,G_obs,orderlist)
    if len(np.setdiff1d(v,anc)) != 0: #len(set(v) - set(anc)) != 0:
        if debug:
            print(">> Enter in Line 2 of the algorithm")

        #Remove non ancestors(Y) node form the graph
        #G_an    = G.copy()
        edgedel = np.array([u for u in list(nx.nodes(G_obs)) if u not in anc])
        G_an = removeedge(edgedel,G,ispresent(hiddvar,G))
        
        #G_an.remove_nodes_from(edgedel)
        
        if P_out["08_prod"]==True or P_out["05_fraction"]==True:
            P_out["03_sumvar"] = toposort(np.union1d(np.setdiff1d(v,anc),P_out["03_sumvar"]), orderlist)  #list((set(v)-set(anc)) | set(P_out["03_sumvar"]))
            P_out["07_sum"]    = True
        else:
            P_out["01_var"]    = anc    
        
        
        P_iter = ID(y,np.intersect1d(x,anc),P_out,G_an,orderlist,hiddvar,**kwargs)

        return P_iter


    #->Line 3
    #--> Get ancestors of Y in G_x_bar: first create G_bar_X, then compute the
    #    ancestor
    G_bar_x  = incomingedgedel(x,G_obs)
    an_bar_x = an(y,G_bar_x,orderlist)
    W        = np.setdiff1d(np.setdiff1d(v,x),an_bar_x) # list((set(v)-set(x)) - set(an_bar_x))
    if len(W) > 0:
        if debug:
            print(">> Enter in Line 3 of the algorithm")
        #x      = x + W
        P_iter = ID(y,np.union1d(x,W),P_out,G,orderlist,hiddvar,**kwargs)
        
        return P_iter


    G_less_x   = removeedge(x,G,hiddvar)
    C_G_less_x = C_components(G_less_x,hiddvar,orderlist)

    #-> Line 4
    if len(C_G_less_x) > 1:
        
        if debug:
            print(">> Enter in Line 4 of the algorithm")
            print(">>> Condition ' if C(G\X)={S1, ..., Sk} is valid';")


        P_iter = []
        for p_iter in range(0,len(C_G_less_x)):
            P_iter_i = ID(C_G_less_x[p_iter],np.setdiff1d(v,C_G_less_x[p_iter]),P_out,G,orderlist,hiddvar,**kwargs)
            P_iter.append(P_iter_i)  
                        
        if len(np.setdiff1d(v,np.union1d(y,x))) != 0:
            P_out["07_sum"]    = True
        P_out["08_prod"]   = True
                
        P_out["03_sumvar"] = toposort(np.setdiff1d(v,np.union1d(y,x)), orderlist)
        P_out["09_branch"] = P_iter
        
        return P_out

    if len(C_G_less_x) == 1:
        if debug:
            print(">> Enter in Line 4 of the algorithm")
            print(">>> Condition ' if C(G\X)={S} is valid';")
            
        C_G = C_components(G,hiddvar,orderlist)
        #   
        # Line 5 
        if len(C_G)==1: 

            if len(np.setdiff1d(v,C_G))==0:  

                raise ValueError('Graph form an Hedge: causal effect is not identifiable')
        
                
        #-> Line 6
        #Add the condition of lenght since I'm dealing with array and If they don't have the same dimension code crash;

        #if any([all(np.sort(C_G_less_x) == np.sort(C_Gi)) for C_Gi in C_G if len(C_G_less_x)==len(C_Gi)]):
        if any([C_G_less_x[0].tolist() == C_Gi.tolist() for C_Gi in C_G]):
            if debug:
                print(">> Enter in Line 6 of the algorithm")
                        
                
            P_prod_b = []
            
            for ii in range(0,len(C_G_less_x[0])):
                
                if P_out["08_prod"]==True:
                    
                    #Compute condition set 
                    cond_list = condset(C_G_less_x[0][ii],v)
                    P_prod_i  = complex_prob(P_out,C_G_less_x[0][ii],cond_list,v,orderlist)
                   
                                        
                else:
                    P_prod_i         = P_out.copy()
                    P_prod_i["01_var"]  = np.array(C_G_less_x[0][ii])
                    
                    cond_list        = condset(C_G_less_x[0][ii],v)
                    
                    P_prod_i["02_cond"] = cond_list
                P_prod_b.append(P_prod_i)                    
                    
            #If the C-component have more than one variable
            if len(C_G_less_x[0]) > 1:
                P_6 = probabilitystr()
                P_6["03_sumvar"]  = toposort(np.setdiff1d(v,np.union1d(y,x)), orderlist)
                P_6["09_branch"]  = P_prod_b.copy()   
                P_6["08_prod"]    = True
            else:
                P_6 = P_prod_b[0].copy()
                if P_6["08_prod"] or P_6["05_fraction"]:
                    P_6["03_sumvar"]  = np.union1d(P_6["03_sumvar"],np.setdiff1d(C_G_less_x[0],y))
                else:
                    P_6["01_var"] =  np.setdiff1d(P_6["01_var"], np.union1d(P_6["03_sumvar"],np.setdiff1d(C_G_less_x[0],y)))
            
            return P_6        
                    
        #-> Line 7
        #For testing if C_G_less_x is a subset of C_Gi I use the difference 
        #function with array: if length is 0 the is a subset
        if any([len(np.setdiff1d(C_G_less_x,C_Gi))==0 for C_Gi in C_G]): 
            if debug:
                print(">> Enter in Line 7 of the algorithm")

            S            = [C_Gi for C_Gi in C_G if len(np.setdiff1d(C_G_less_x,C_Gi))==0][0]
            P_           = probabilitystr()
            P_prod_b     = []

            for iii in range(0,len(S)):
                
                if P_out["08_prod"]==True:
                    
                    #Compute condition set 
                    cond_list = condset(S[iii],v)
                    P_prod_i  = complex_prob(P_out,S[iii],cond_list,v,orderlist)
                    
                else:
                    P_prod_i         = P_out.copy()
                    P_prod_i["01_var"]  = np.array(S[iii])
                    
                    cond_list        = condset(S[iii],v)
                    
                    P_prod_i["02_cond"] = cond_list
                    
                P_prod_b.append(P_prod_i) 
                     
            if len(S) > 1:
                P_["08_prod"]   = True                        
                P_["09_branch"]  = P_prod_b.copy()
            else:
                P_ = P_prod_b[0].copy()
                
            edgedel = [u for u in list(nx.nodes(G_obs)) if u not in S]
            G_s     = removeedge(edgedel,G,hiddvar)
            P_iter  = ID(y,np.intersect1d(x,S),P_,G_s,orderlist,hiddvar,**kwargs)
    
            return P_iter  
        
def complex_prob(P, var, condset, v, orderlist):
    
    if len(condset) == 0:
        P_complex = P.copy()    
        if P_complex["07_sum"]==True:
            P_complex["03_sumvar"] = toposort(np.union1d(P_complex["03_sumvar"],np.setdiff1d(v,np.union1d(var,condset))), orderlist)
        else:
            P_complex["07_sum"] = True
            P_complex["03_sumvar"] = toposort(np.setdiff1d(v,np.union1d(var,condset)), orderlist)
            
   #If we have the condition set express the conditional probability in terms of 
   # conjunction probability require a denominator
   
    else:
        P_complex = probabilitystr()    
        P_complex["05_fraction"]    = True
        P_complex["04_numerator"]   = P.copy()
        P_complex["06_denominator"] = P.copy()
        if P["07_sum"]==True:
            P_complex["04_numerator"]["03_sumvar"]   = toposort(np.union1d(P["03_sumvar"],np.setdiff1d(v,np.union1d(var,condset))), orderlist)
            P_complex["06_denominator"]["03_sumvar"] = toposort(np.union1d(P["03_sumvar"],np.setdiff1d(v,condset)), orderlist)
        else:
            P_complex["04_numerator"]["03_sumvar"]   = toposort(np.setdiff1d(v,np.union1d(var,condset)), orderlist)
            P_complex["06_denominator"]["03_sumvar"] = toposort(np.setdiff1d(v,condset), orderlist)
           
       
    return P_complex
      
    
def probabilitystr():
      return{  "01_var" : np.array([])
             , "02_cond" : np.array([])
             , "03_sumvar" : np.array([])
             , "04_numerator" : []
             , "05_fraction" : False
             , "06_denominator" : []
             , "07_sum" : False
             , "08_prod" : False
             , "09_branch" : [] }
    
    
    
def IDC(y, x, z, P, G, orderlist, hiddvar,**kwargs):
    '''
    Function implement the ID algorihtm as presented in Pearl 2006b 
    "Identification of Conditional Interventional Distributions" 
    (https://ftp.cs.ucla.edu/pub/stat_ser/r329-uai.pdf)
    
    
    Parameters
    ----------
    y         : array of viarables on which measures the effect of intervention
    x         : array of variables on which the intervention (do) is done
    z         :  array containing condition set
    P         : probability structure defined as the output of function 
                probabilitystr
    G         : netwrokx object
    orderlist : topological order
    hiddvar   : set of hidden (confonuders) nodes  
    **kwargs  : dictionary with following jey
                debug : if True message are print

    Raises
    ------
    ValueError: if the causal effect is not identifiable then the function 
                raise an error

    Returns
    -------
    probability structure that define the identification of causal effect

    '''
    
    debug  = kwargs["debug"]
    
    #Transform array y, z and z in array
    y_set = set(list(y))
    
    if all(x) != "":
        x_set = set(list(x))
    

        #Transform z in set
    z_set = set(list(z))
    #Get graph without incoming edge on x and outgoing edge to z
    G_under_x_bar_z = incomingedgedel(x,outcomeedgedel(z,G))
    for zi in z:
        c_set = x_set.copy()
        c_set | (z_set - {zi})

        if nx.d_separated(G_under_x_bar_z,y_set,{zi},c_set):
            if debug:
                print("> Iterate over IDC: Z variable is " +str(zi))
            
            P_numerator = IDC(y, np.union1d(x,np.array([zi])),np.setdiff1d(z,np.array([zi])), P, G, orderlist, hiddvar,**kwargs)
            
            return P_numerator
                   
    if debug:
        print("> Iterate over ID")
    P_numerator = ID(np.union1d(y,z), x, P, G, orderlist, hiddvar,**kwargs)

    return P_numerator        

#/---------------------------------------------------------------------------\
#|                          Main  function Level - 1                         |
#\---------------------------------------------------------------------------/
        
def ce_identification(y,x,z,P,G,orderlist,hiddvar,**kwargs):
   '''
    This function complete the usage of ID and in particular the IDC
    algorithm. When z is not empty set the final output is the fraction of
    identified distribution P' and sum(over y) wrt P'. So the output of IDC
    algorithm is manipulated in order to include the fraction in the 
    probability structure

    Parameters
    ----------
    y         : array of viarables on which measures the effect of intervention
    x         : array of variables on which the intervention (do) is done
    z         :  array containing condition set
    P         : probability structure defined as the output of function 
                probabilitystr
    G         : netwrokx object
    orderlist : topological order
    hiddvar   : set of hidden (confonuders) nodes  
    **kwargs  : dictionary with following jey
                debug : if True message are print

    Raises
    ------
    ValueError: if the causal effect is not identifiable then the function 
                raise an error

    Returns
    -------
    probability structure that define the identification of causal effect

    '''
 
   debug  = kwargs["debug"]

   if z.tolist():
       
       if debug:
           print("@@ Condition set z is not empty: IDC version is used")
           
       P_numerator = IDC(y, x, z, P, G, orderlist, hiddvar,**kwargs)
       
       P_prime = probabilitystr()
       
       #Construct the denominator
       P_prime_denominator = P_numerator.copy()
       
       if P_prime_denominator["07_sum"]:
           P_prime_denominator["03_sumvar"] = np.union1d(y,P_prime_denominator["03_sumvar"])
       else:
           P_prime_denominator["07_sum"]    = True
           P_prime_denominator["03_sumvar"] = y
           
       #Add the denominator to whole distribution
       P_prime["05_fraction"]    = True
       P_prime["06_denominator"] = P_prime_denominator
       P_prime["04_numerator"]   = P_numerator

   else:
       if debug:
           print("@@Condition set z is empty: ID version is used")
           
       P_prime = ID(y,x,P,G,orderlist,hiddvar,**kwargs)
       
       
   return P_prime
      