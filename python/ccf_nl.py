'''
Created on April, 2017

@author: Jose_Cobena
'''
import numpy as np
import itertools
import copy
from itertools import count, zip_longest
from collections import OrderedDict

def main():
    n_atoms = 756
    cutoff = 3.5
    T = 240
    initial_step = 0


    file_data_time_t0 = 'uni-file-{:08}.txt'.format(initial_step)
    output_file_ccf = 'ccf-results-{}.txt'.format(T)
    

    #-------------------------------------------------------
    with open(file_data_time_t0, 'r') as file_t0, open(output_file_ccf, 'w') as output_file:
        
        oxygen_list = my_ox_list(file_t0)
        #print(oxygen_list)
        neigh_list_t0 = my_nl_list(oxygen_list, cutoff)
        #print(neigh_list_t0)
        oxy_for_dp_t0 = bool_list_evolution(neigh_list_t0, neigh_list_t0)
        #print(oxy_for_dp_t0)
        
 
        #dot product at t =0
        ccf_t0 =0
        for vector_t0, vector_t0, in zip_longest(oxy_for_dp_t0, oxy_for_dp_t0):
                ccf_t0 += np.dot(vector_t0[1], vector_t0[1])
        ccf_t00 = ccf_t0/len(oxy_for_dp_t0)
        
        #print(ccf_t00)
#---------------------------------------------------------------
        ccf_t2 = ccf_t(0, 10, 2,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t10 = ccf_t(10, 100, 10, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t100 = ccf_t(100, 1000, 100,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t1000 = ccf_t(1000, 10000, 1000,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t10000 = ccf_t(10000, 100000, 10000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t100000 = ccf_t(100000, 1000000, 100000,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t1000000 = ccf_t(1000000, 10000001, 1000000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        ccf_t10000000 = ccf_t(10000000, 50000001, 1000000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
         
        for x in ccf_t2:
            print(*x, file = output_file)
             
        for x in ccf_t10:
            print(*x, file = output_file)
             
        for x in ccf_t100:
            print(*x, file = output_file)
             
        for x in ccf_t1000:
            print(*x, file = output_file)
             
        for x in ccf_t10000:
            print(*x, file = output_file)
             
        for x in ccf_t100000:
            print(*x, file = output_file)
             
        for x in ccf_t1000000:
            print(*x, file = output_file)
#             
        for x in ccf_t10000000:
            print(*x, file = output_file)
#             
            
#-------------------------------------------------------------------
def ccf_t(ini_step, final_step, increments, neigh_list_t0  ,oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff):
    ccf_list = []
    
    for step in range(ini_step, final_step, increments):
           #fc = 'ccf-in%.8d.txt' % i    #old format
            fc = 'uni-file-{:08d}.txt'.format(step)
            file_tt = open(fc, 'r')
            
            oxygen_list = my_ox_list(file_tt)
            #print(oxygen_list)
            neigh_list = my_nl_list(oxygen_list, cutoff)
            #print(neigh_list)
            oxy_for_dp_tt = bool_list_evolution(neigh_list_t0, neigh_list)
            #print(oxy_for_dp_tt)
            
    
            #print(oxy_for_dp_tt)
            ccf_it = 0.0
            for vector_t0, vector_tt, in zip_longest(oxy_for_dp_t0, oxy_for_dp_tt):
                #print(vector_t0, vector_tt)
                ccf_it += np.dot(vector_t0[1], vector_tt[1])
            ccf_tave = ccf_it/n_atoms
            ccf_tt = ccf_tave/ccf_t00
            to_return = [step, ccf_tt]
            ccf_list.append(to_return)
            #print(ccf_tt)
    return ccf_list

def bool_list_evolution(oxy_ids_lists_t0, oxy_ids_lists_tt):
    neigh_atoms_time0 = copy.deepcopy(oxy_ids_lists_t0)
    neigh_atoms_timet = copy.deepcopy(oxy_ids_lists_tt)
    dic_oxy = {}
    #boolean_neighs_big_tt = []
    #loop simultaneously over the list at t=0 to list at t=t:
    for neighs_t0, neighs_tt in zip_longest(neigh_atoms_time0, neigh_atoms_timet):#, fillvalue = 0):
        keys = int(neighs_t0[0])
        
        if (neighs_t0[0] == neighs_tt[0]):    
           boolean_neighs_tt = []
           if neighs_t0[1]:
               
               for neighs_t00 in neighs_t0[1]:
                   if neighs_tt[1]:
                       for neighs_ttt in neighs_tt[1]:

                           if neighs_t00 in neighs_tt[1]:
                               boolean_neighs_tt.append(int(1))
                               break
                           else:
                               boolean_neighs_tt.append(int(0)) 
                               #print('yes')
                               break
                   else:
                         boolean_neighs_tt.append(int(0))        
                         #print('yes1')
           else:
               boolean_neighs_tt.append(int(0))
               #print('yes2')               
        else:
            #print('yes3')
            boolean_neighs_tt.append(int(0)) 
        val = boolean_neighs_tt
        dic_oxy[keys] = val
    Ordic = sorted(dic_oxy.items())
#     Ordic = OrderedDict(sorted(dic_oxy.items()))
#     for k,v in Ordic.items():
#         print(k,v)
    return (Ordic)
            
def my_ox_list(step0_file):
#     list_b = copy.deepcopy(oxy_working)
    #     oxy0 = copy.deepcopy(oxy_base_lists)#
    dict_atoms = {}
    coords = []
    for skip in range(9):
        next(step0_file)
    for rows in step0_file:
        values = rows.split()
        if values[1] == '4':
            keys = int(values[0])
            coords = [float(values[3]), float(values[4]), float(values[5])]
            dict_atoms[keys] = coords   
        
    positions = sorted(dict_atoms.items())    
    return positions
            
            
def my_nl_list(positions, cutoff):
    super_list = []
    dict_atoms = {}
    list_b = copy.deepcopy(positions)
    #oxy0 = copy.deepcopy(oxy_base_lists)#
    
    for i in positions:
        keys = int(i[0])
        inter_list = []
        for j in list_b:
            if i[0] != j[0]:
                r_test = (j[1][0] - i[1][0])*(j[1][0] - i[1][0]) + (j[1][1] - i[1][1])*(j[1][1] - i[1][1]) + \
                                 (j[1][2] - i[1][2])*(j[1][2] - i[1][2])
                             
                if np.sqrt(r_test) < cutoff:
                    inter_list.append(j[0])
                
        dict_atoms[keys] = inter_list
                
    positions_f = sorted(dict_atoms.items())    
    return positions_f



if __name__ == "__main__": main()
