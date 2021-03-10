
#//=========================================================================\\
#||                        Call function: example                           ||
#\\=========================================================================//

# In this code I test the algorithm with some examples taken from the 
# literature. An important note is the follow: for a correct implementation
# when you specify the biderected node you have to make explicit the hidden
# node when you create the DAG (I identify hiddent node in my example with
# Ui). 



#/---------------------------------------------------------------------------\
#|                         Call function: example                            |
#\---------------------------------------------------------------------------/

#------->          Causality Pearl page 92  figure 3.8(f)      <--------------
#-> Construct the graph
DG_1 = nx.DiGraph()
DG_1.add_edges_from([ ("X","Y")
                   ,("Z1","Z2")
                   ,("Z1","Y")
                   ,("Z2","Y")
                   ,("Z1","Z2")
                   ,("Z2","X")
                   ,("U1","Y")
                   ,("U1","Z1")
                   ,("U2","X")
                   ,("U2","Z2")
                    ])
nx.draw(DG_1,with_labels=True)
  
#Preparing input parameters
hiddvar          = np.array(["U1","U2"])                       # declare an array containign the hidden variables
Gobs1            = removeedge(hiddvar,DG_1,hiddvar)            # gets the observed graph (so without hidden variable)
order            = np.array(list(nx.topological_sort(Gobs1)))  # retreive the topological ordering of the graph
kwargs           = dict()                                      # specify others input
kwargs["debug"]  = True                                        # diplay the enter interrogated at each step
kwargs["debug2"] = False                                       # if true, edges of graphs at each steps are diplayed
P                = probabilitystr()                            # initalize probability structure

#Compute distribution
Pi1 = ce_identification(np.array(["Y"]),np.array(["X"]),np.array([]), P,DG_1,order,hiddvar,**kwargs)           
 
#------->          Causality Pearl page 92  figure 3.8(f)      <--------------
#-> Construct the graph
DG_2 = nx.DiGraph()
DG_2.add_edges_from([ ("X","Z1")
                     ,("Z1","Y")
                     ,("Z3","Y")
                     ,("Z2","Z3")
                     ,("Z2","Z1")
                     ,("Z2","X")
                     ,("U1","X")
                     ,("U1","Z2")
                     ,("U2","X")
                     ,("U2","Z3")
                     ,("U3","X")
                     ,("U3","Y")
                     ,("U4","Y")
                     ,("U4","Z2")
                    ])

nx.draw(DG_2,with_labels=True)
  

hiddvar2          = np.array(["U1","U2","U3","U4"])             # declare an array containign the hidden variables
Gobs2             = removeedge(hiddvar,DG_2,hiddvar)            # gets the observed graph (so without hidden variable)
order2            = np.array(list(nx.topological_sort(Gobs2)))  # retreive the topological ordering of the graph
kwargs2           = dict()                                      # specify others input
kwargs2["debug"]  = True                                        # diplay the enter interrogated at each step
kwargs2["debug2"] = False                                       # if true, edges of graphs at each steps are diplayed
P2                = probabilitystr()                            # initalize probability structure

Pi2 = ce_identification(np.array(["Y"]),np.array(["X"]),np.array([]), P2, DG_2, order2, hiddvar2, **kwargs2)     

#------->          Paper Pearl 2006a Figure 1(b)      <--------------
# Paper in which the ID algorith is presented 
# https://www.aaai.org/Papers/AAAI/2006/AAAI06-191.pdf
#-> Construct the graph
DG_3 = nx.DiGraph()
DG_3.add_edges_from([ ("W1","X")
                     ,("X","Y1")
                     ,("W2","Y2")
                     ,("U1","W1")
                     ,("U1","Y1")
                     ,("U2","W1")
                     ,("U2","Y2")
                     ,("U3","W1")
                     ,("U3","W2")
                     ,("U4","X")
                     ,("U4","W2")
                     ,("W1","W2")
                    ])
# !! This graph is not indentifiable. Remove last edge in order to get an
# identifiable graph

nx.draw(DG_3,with_labels=True)
  

hiddvar3          = np.array(["U1","U2","U3","U4"])             # declare an array containign the hidden variables
Gobs3             = removeedge(hiddvar,DG_3,hiddvar)            # gets the observed graph (so without hidden variable)
order3            = np.array(list(nx.topological_sort(Gobs3)))  # retreive the topological ordering of the graph
kwargs3           = dict()                                      # specify others input
kwargs3["debug"]  = True                                        # diplay the enter interrogated at each step
kwargs3["debug2"] = False                                       # if true edges of graphs at each steps are diplayed
P3                = probabilitystr()                            # initalize probability structure

Pi3 = ce_identification(np.array(["Y1","Y2"]),np.array(["X"]),np.array([]), P3, DG_3, order3, hiddvar3, **kwargs3)     



