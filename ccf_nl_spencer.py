'''
Created on September, 2018

@author: Jose_Cobena, Spencer Ortega
'''
import numpy as np
import itertools
import copy
from itertools import count, zip_longest
from collections import OrderedDict

def main():
	# Defined variables
    n_atoms = 756
    cutoff = 3.5
    T = 240
    initial_step = 0

    # file names
    file_data_time_t0 = 'uni-file-{:08}.txt'.format(initial_step)
    output_file_ccf = 'ccf-results-{}.txt'.format(T)
    

    #-------------------------------------------------------
    with open(file_data_time_t0, 'r') as file_t0, open(output_file_ccf, 'w') as output_file:
        
        # get list of oxygen atoms
        oxygen_list = my_ox_list(file_t0)
        #print(oxygen_list)

        # get neighbor list at t0
        neigh_list_t0 = my_nl_list(oxygen_list, cutoff)
        #print(neigh_list_t0)

        # compute binary neighbor list
        oxy_for_dp_t0 = bool_list_evolution(neigh_list_t0, neigh_list_t0)
        #print(oxy_for_dp_t0)
        
 
        # dot product at t = 0
        ccf_t0 = 0
        for vector_t0, vector_t0, in zip_longest(oxy_for_dp_t0, oxy_for_dp_t0):
                ccf_t0 += np.dot(vector_t0[1], vector_t0[1])

        # = 1
        ccf_t00 = ccf_t0/len(oxy_for_dp_t0)
        
        #print(ccf_t00)
	#---------------------------------------------------------------
        ccf_t2 = ccf_t(0, 10, 2,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        #ccf_t10 = ccf_t(10, 100, 10, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t100 = ccf_t(100, 1000, 100,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t1000 = ccf_t(1000, 10000, 1000,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t10000 = ccf_t(10000, 100000, 10000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t100000 = ccf_t(100000, 1000000, 100000,  neigh_list_t0, oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t1000000 = ccf_t(1000000, 10000001, 1000000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
        # ccf_t10000000 = ccf_t(10000000, 50000001, 1000000, neigh_list_t0,  oxy_for_dp_t0, ccf_t00,  n_atoms, cutoff)
         
        for x in ccf_t2:
            print(*x, file = output_file)
             
        # for x in ccf_t10:
        #     print(*x, file = output_file)
             
#         for x in ccf_t100:
#             print(*x, file = output_file)
             
#         for x in ccf_t1000:
#             print(*x, file = output_file)
             
#         for x in ccf_t10000:
#             print(*x, file = output_file)
             
#         for x in ccf_t100000:
#             print(*x, file = output_file)
             
#         for x in ccf_t1000000:
#             print(*x, file = output_file)
# #             
#         for x in ccf_t10000000:
#             print(*x, file = output_file)
# #             
            
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

            ccf_it = 0.0

            # for every atom in time step t
            for vector_t0, vector_tt, in zip_longest(oxy_for_dp_t0, oxy_for_dp_tt):
            	# dot product of binary neighbor list at 0 and t
                ccf_it += np.dot(vector_t0[1], vector_tt[1])

			# average of ccf for every atom in time step t              
            ccf_tave = ccf_it/n_atoms

            # ccf at t / ccf at 0
            ccf_tt = ccf_tave/ccf_t00

            # package ccf calculation with step 
            to_return = [step, ccf_tt]

            ccf_list.append(to_return)
            #print(ccf_tt)
    return ccf_list

def bool_list_evolution(oxy_ids_lists_t0, oxy_ids_lists_tt):
	neigh_atoms_time0 = copy.deepcopy(oxy_ids_lists_t0)
	neigh_atoms_timet = copy.deepcopy(oxy_ids_lists_tt)
	dic_oxy = {}

	#loop simultaneously over the list at t=0 to list at t=t:
	count = 0
	for neighs_t0, neighs_tt in zip_longest(neigh_atoms_time0, neigh_atoms_timet):#, fillvalue = 0):
		
		keys = int(neighs_t0[0])
		
		# if neighs_t0[0] != neighs_tt[0]:
		# 	print("neighs_t0 atom=: ",neighs_t0[0])
		# 	print("neighs_tt atom=: ",neighs_tt[0])


		# neighs[0] : atom id
		# neighs[1] : list of neighbor atom ids


		# if atom id of t0 and tt are the same
		# AND
		# if neighbor list of t0 is not empty
		if neighs_t0[0] == neighs_tt[0] and neighs_t0[1]:  

        	## possible move out of if statement
			boolean_neighs_tt = []

			# if neighbor list of tt is empty
			# boolean_neighs_tt = [0,...,0] of size length(neighs_t0[1])
			if not neighs_tt[1]:
				boolean_neighs_tt.append( [0]*len(neighs_t0[1]))
				break

			# for every neighbor in t0
			# check if present in neighbor list of tt
			for neighs_t00 in neighs_t0[1]:
				# if neighbor neighs_t00 is in neighbor list of tt
				if neighs_t00 in neighs_tt[1]:
					# add 1 to binary list
					boolean_neighs_tt.append(int(1))
					#print(neighs_t00, " = ", 1)

				# if neighs_t00 is no longer in in neighbor list at tt
				else:
					# add 0 to binary list
					boolean_neighs_tt.append(int(0)) 
					#print(neighs_t00, " = ", 0)		

		# if atom id of t0 and tt are NOT the same
		# ****ask jose about this case*
		# 	why are we storing bool_neigh_tt using atom id of t0 as key
		# OR
		# if neighbor list of t0 is empty
		else:
			# add 0 to binary list
			boolean_neighs_tt.append(int(0)) 

		# store key/value pair
		# key = atom id
		# val = binary neighbor list
		val = boolean_neighs_tt
		dic_oxy[keys] = val

		#print(val)

	# sort dictionary by atom id
	Ordic = sorted(dic_oxy.items())

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
            
  
# find neighbor list lamps function            
def my_nl_list(positions, cutoff):
	super_list = []
	dict_atoms = {}
	list_b = copy.deepcopy(positions)

	for i in positions:
		keys = int(i[0])
		inter_list = []

		for j in list_b:
			# atom ids are not the same
			if i[0] != j[0]:
				# Compute distance
				# sqrt(  (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2  )
				r_test = (j[1][0] - i[1][0])*(j[1][0] - i[1][0]) + (j[1][1] - i[1][1])*(j[1][1] - i[1][1]) + (j[1][2] - i[1][2])*(j[1][2] - i[1][2])
                
				# if distace is within cutoff

                if np.sqrt(r_test) < cutoff
                    inter_list.append(j[0])
                
        dict_atoms[keys] = inter_list
                
    positions_f = sorted(dict_atoms.items())    
    return positions_f



if __name__ == "__main__": main()
