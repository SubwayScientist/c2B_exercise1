This directory contains a few test sample figures corresponding to different cluster sizes and parameters.



For each cluster :
    fig_a = M_per fermi arc at n=0.95      (except for 2x2, where n =1.05 for some reason)
    fig_b = M_per self for n=0.95          (except for 2x2, where n =1.05       ,,       )
    fig_c = M_per self for n=0.97          (except for 2x2, where n =1.03       ,,       )
    fig_d = M_per self for n=0.99          (except for 2x2, where n =1.01       ,,       )
    


Hamiltonian parameters that are the same for every figures :
    t = 1,   tp = -0.2,   tpp = 0,   delta = 1.1,   eta = 0.5
    


    | 2x2 Cluster     |    Nc=8 Cluster    |    4x3 Cluster   |   4x4 Cluster   |         
 a) |   MU =  -0.85   |     MU =  -0.90    |     MU =  -0.90  |    MU =  -1.05  |
 b) |   MU =  -0.95   |     MU =    -      |     MU =  -0.90  |    MU =  -1.05  |
 c) |   MU =  -1.00   |     MU =    -      |     MU =  -0.80  |    MU =  -0.90  |     
 d) |   MU =  -1.06   |     MU =    -      |     MU =  -0.74  |    MU =  -0.70  |


The Nc = 8 cluster needs to be updated as the M_per self-energy was completely asymmetrical.
The 4x3 cluster was also a bit asymmetrical, hence it would need to be looked at.



These figures were made on commit : 98d0ae5
