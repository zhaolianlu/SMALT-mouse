######################################################################################
## Python scripts to infer parameters by simulating 3D tumor growth and             ##
## SMALT lineage tracing system. Deme subdivision is assumed to                     ##
## model cell mixing and spatial contraint.                                         ##
##                                                                                  ##
## *Spatial model: pripheral growth                                                 ##
## ABC inference: Duo Xie                                                           ##
## Original Model Author: Zheng Hu                                                  ##
## Original Date: 11/15/2022                                                        ##
## Update Date: 09/13/2022                                                          ##
######################################################################################
###This script is used to perfom analysis of Approximate Bayesian Computation (ABC) of methods

from pyabc import *
import tempfile
import pandas as pd
from pandas import Series,DataFrame
import matplotlib.pyplot as plt
import math
import sys,os,math,random
import numpy as np
from collections import Counter
import scipy.spatial.distance as scipydis
from Bio import SeqIO, AlignIO
from decimal import Decimal,ROUND_HALF_UP
from scipy.spatial.distance import cityblock
######functions
class deme():
    def __init__(self):
        self.present= 0         ## whether the deme is empty or occupied: 0-empty;1-occupied
        self.neutral = []    	## the neutral founder lineage after tumor tranformation
        self.advant = []        ## the advantageous cells having one driver mutation
        self.empty = 27         ## number of empty sites in the neighbourhood

####round in the right way
def right_round(num,keep_n):
	if isinstance(num,float):
		num = str(num)
	return Decimal(num).quantize((Decimal('0.' + '0'*keep_n)),rounding=ROUND_HALF_UP)

def seq2stat(data):
	MB_S = data.sum(axis=1)
	mean_S = np.mean(MB_S)
	std_S = np.std(MB_S)
#####calculate SFS
	mutation_S = data.sum(axis=0)
	SFS_S = mutation_S/ncells
	quantile_25_S = np.quantile(MB_S, 0.25)
	quantile_75_S = np.quantile(MB_S, 0.75)
	SFS_5to10percent = len(SFS_S[(SFS_S>=0.05) & (SFS_S<0.1)])
	SFS_10to25percent = len(SFS_S[(SFS_S>=0.1) & (SFS_S<0.25)])
	SFS_25to50percent = len(SFS_S[(SFS_S>=0.25) & (SFS_S<0.5)])
	SFS_50to75percent = len(SFS_S[(SFS_S>=0.5) & (SFS_S<0.75)])
	PD_S = np.array([])
	for i in range(0, 1000):
		idx = np.random.randint(ncells, size=2)
		sampled_S = data[idx,:]
		D_S = cityblock(sampled_S[0], sampled_S[1])
		PD_S = np.hstack([PD_S ,D_S])
	mean_PD = np.mean(PD_S)
	std_PD = np.std(PD_S)
	PD_quantile_25 = np.quantile(PD_S, 0.25)
	PD_quantile_75 = np.quantile(PD_S, 0.75)
	sum_stat = np.hstack([mean_S, std_S, quantile_25_S,quantile_75_S, mean_PD, std_PD, PD_quantile_25,PD_quantile_75,SFS_5to10percent,SFS_10to25percent,SFS_25to50percent,SFS_50to75percent ])
	return sum_stat

def obs2stat(data):
	mutation_S = data.sum(axis=0)
	ncells = len(data)
	SFS_S = mutation_S/ncells

	# Modify seq based on SFS_S
	for col in range(data.shape[1]):
		if SFS_S[col] > 0.75:
			data[:, col] = 0
	MB_S = data.sum(axis=1)
	mean_S = np.mean(MB_S)
	std_S = np.std(MB_S)
#####calculate SFS
	mutation_S = data.sum(axis=0)
	SFS_S = mutation_S/ncells
	quantile_25_S = np.quantile(MB_S, 0.25)
	quantile_75_S = np.quantile(MB_S, 0.75)
	SFS_5to10percent = len(SFS_S[(SFS_S>=0.05) & (SFS_S<0.1)])
	SFS_10to25percent = len(SFS_S[(SFS_S>=0.1) & (SFS_S<0.25)])
	SFS_25to50percent = len(SFS_S[(SFS_S>=0.25) & (SFS_S<0.5)])
	SFS_50to75percent = len(SFS_S[(SFS_S>=0.5) & (SFS_S<0.75)])
	PD_S = np.array([])
	for i in range(0, 1000):
		idx = np.random.randint(ncells, size=2)
		sampled_S = data[idx,:]
		D_S = cityblock(sampled_S[0], sampled_S[1])
		PD_S = np.hstack([PD_S ,D_S])
	mean_PD = np.mean(PD_S)
	std_PD = np.std(PD_S)
	PD_quantile_25 = np.quantile(PD_S, 0.25)
	PD_quantile_75 = np.quantile(PD_S, 0.75)
	sum_stat = np.hstack([mean_S, std_S, quantile_25_S,quantile_75_S, mean_PD, std_PD, PD_quantile_25,PD_quantile_75,SFS_5to10percent,SFS_10to25percent,SFS_25to50percent,SFS_50to75percent ])
	return sum_stat

def readdata(phy):
	seqs = []
	align = AlignIO.read(phy, "phylip-relaxed")
	for index, record in enumerate(align._records):
		if record.id == 'ref':
			next
		else:
			seqs.append(record.seq)
	seqs = np.array(seqs)
	seqs = seqs.astype(int)
	return seqs
	
def distance(simulation, data):
####calculate mean and SD
	dis = scipydis.euclidean(simulation["X_2"], data["X_2"])
	return dis

def createLattice(d):
    """
    Create a 3D cubic lattice with side length of 2d+1 where each site contains a empty deme.
    """
    lattice = {}
    for x in range(0,2*d+1):
        for y in range(0,2*d+1):
            for z in range(0,2*d+1):
                lattice[(x,y,z)] = deme()
    return lattice


def neighbor26(xxx_todo_changeme):
    """
    Moore neighbourhood: 26 neighbour sites of (a,b,c)
    """
    (a,b,c) = xxx_todo_changeme
    neighbor = [(a+i,b+j,c+k)
                for i in [-1,0,1]
                for j in [-1,0,1]
                for k in [-1,0,1]
                if not (i==0 and j==0 and k==0)]

    return neighbor

def localNeighbor(xxx_todo_changeme1,r):
    """
    A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice.
    """
    (a,b,c) = xxx_todo_changeme1
    neighbor = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbor += [(a+x,b+y,c+z)]
    return neighbor


def traceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell.
    For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell

    mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',')  # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0]      # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace


def lowerORupper(value):
    """
    A function to choose the upper or lower integral value given a non-integral number
    """
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirstDeme_t1(maxsize,lineage,current_id,current_driver_id,sfit,barcode_muts, total_mut_rate, quies_rate, birth_rate, adv_rate, hotspot_chance, hotspot_number, non_hotspot_number): 
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    t1 - one-tier driver model
    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    current_id - the starting mutation ID
    sfit - selection fitness of advantageous mutations
    """
#    neu_list = [str(current_id)]
######ensuring there is 10 mutations at the beginning
    neu_list = ["0" for i in range(10)]
    adv_list = []
   # current_deme_size = 1
    current_deme_size = len(neu_list)+len(adv_list)
    #total_mut_rate2 = 0
    total_mut_rate2 = total_mut_rate
    while current_deme_size < maxsize:
        n1,n2 = len(neu_list),len(adv_list)
        n1_double = 0
        n2_double = 0
        neu_div_list_double = []
        adv_div_list_double = []
        neu_qui_list = []
        adv_qui_list = []
        if n1 > 0:
            neu_qui_number = lowerORupper(n1*quies_rate)
            neu_div_number = lowerORupper(n1*birth_rate)
#            neu_qui_number = int(n1*quies_rate)
#            neu_div_number = int(n1*birth_rate+1)                            #number of dividing cells in this generation
            neu_pass_number = neu_qui_number+neu_div_number
            random.shuffle(neu_list)
            neu_qui_list = neu_list[0:neu_qui_number]
            neu_div_list = neu_list[neu_qui_number:neu_pass_number]
            neu_div_list_double = neu_div_list*2
            n1_double = len(neu_div_list_double)

        if n2 > 0:
            adv_qui_number = lowerORupper(n2*(quies_rate-sfit))        #number of dividing cells in this generation
            #adv_qui_number = lowerORupper(n2*(0.9-0.1-(1+sfit)*(birth_rate-0.1)))
            #adv_div_number = lowerORupper(n2*(0.1+(1+sfit)*(birth_rate-0.1)))
            adv_div_number = lowerORupper(n2*(birth_rate+sfit))
            adv_pass_number = adv_qui_number+adv_div_number
            if adv_pass_number > n2:
                adv_qui_number = int(n2*(quies_rate-sfit))
                adv_pass_number = adv_qui_number+adv_div_number
            random.shuffle(adv_list)
            adv_qui_list = adv_list[0:adv_qui_number]
            adv_div_list = adv_list[adv_qui_number:adv_pass_number]
            adv_div_list_double = adv_div_list*2
            n2_double = len(adv_div_list_double)

        if n1_double > 0:
            new_mut1 = np.random.poisson(total_mut_rate2*n1_double)
			#assign new_mut1 mutations to n1_double cells
			#####Counter({4: 2, 0: 1, 7: 1, 5: 1, 8: 1}) #th cell: mutation number
            mut_assig1 = Counter(np.random.choice(n1_double,new_mut1))
            for x1 in list(mut_assig1.keys()):
                nmut = mut_assig1[x1]
                new_mut1 = list(range(current_id+1,current_id+1+nmut))
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_div_list_double[x1]]
                    if random.random() < hotspot_chance:
                        barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
                        barcode_muts += [barcode_mut_pos]
                    else:
                        barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
                        barcode_muts +=  [barcode_mut_pos]
                neu_div_list_double[x1] = mut_str

        if n2_double > 0:
            new_mut2 = np.random.poisson(total_mut_rate2*n2_double)
            mut_assig2 = Counter(np.random.choice(n2_double,new_mut2))
            for x2 in list(mut_assig2.keys()):
                nmut = mut_assig2[x2]
                new_mut2 = list(range(current_id+1,current_id+1+nmut))
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_div_list_double[x2]]
                    if random.random() < hotspot_chance:
                        barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
                        barcode_muts += [barcode_mut_pos]
                    else:
                        barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
                        barcode_muts += [barcode_mut_pos]
                    #cut_type = random.choice(range(0,cut_type_number))+1
                adv_div_list_double[x2] = mut_str

        if random.random() < adv_rate*n1_double:
            current_driver_id += 1
            #current_n1 = len(neu_list)
            #lineage += [str(neu_list[current_n1-1])]
            adv_div_list_double += [neu_div_list_double[n1_double-1]]
            neu_div_list_double = neu_div_list_double[0:n1_double-1]

        neu_list=neu_qui_list+neu_div_list_double
        adv_list=adv_qui_list+adv_div_list_double
        current_deme_size=len(neu_list)+len(adv_list)

    return neu_list,adv_list,current_id,current_driver_id,lineage,barcode_muts


def demeGrowthFission_t1(neu_list,adv_list,lineage,current_id,current_driver_id,current_deme_number,sfit,barcode_muts, deme_size, quies_rate, birth_rate, total_mut_rate, hotspot_chance, hotspot_number, non_hotspot_number, adv_rate):
    """
    A function to simulate deme growth and fission and keep track of the mutational lineages
    t1 - one-tier driver model
	The function takes as input the lists of neutral and advantageous mutations in the deme before growth and division, the mutation lineage, the mutation ID, the driver mutation ID, the current deme number, the coefficient for the selection coefficient of advantageous mutations, and the barcode mutation site.

	First, the function generates a random number of neutral and advantageous mutations that occur during the growth phase. The number of neutral mutations is determined by sampling from a Poisson distribution with a mean equal to the total mutation rate times the growth time. The number of advantageous mutations is determined by sampling from a binomial distribution with the number of trials equal to the number of neutral mutations and a success probability equal to the probability of an advantageous mutation. These new mutations are added to the list of neutral and advantageous mutations for the deme.

	Next, the function simulates the division of the deme into two daughter demes. The lists of neutral and advantageous mutations are divided between the two daughter demes randomly, with the possibility of some mutations being lost in the process. The function also simulates the occurrence of new mutations in the daughter demes. For each new mutation, the function generates a random number that represents the position of the mutation on the SMALT barcode sequence, and assigns a unique mutation ID to the mutation. If the mutation is advantageous, the function also assigns a unique driver mutation ID to the mutation. The function updates the mutation lineage by appending the new mutation IDs to the lineage of the parent cell. The function also updates the barcode mutation site by adding the position of the new mutation
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:
        n1,n2 = len(neu_list),len(adv_list)
        n1_double = 0
        n2_double = 0
        neu_div_list_double = []
        adv_div_list_double = []
        neu_qui_list = []
        adv_qui_list = []
        if n1 > 0:
            neu_qui_number = lowerORupper(n1*quies_rate)
            neu_div_number = lowerORupper(n1*birth_rate)                            #number of dividing cells in this generation
            neu_pass_number = neu_qui_number+neu_div_number
            if neu_pass_number > n1:
                neu_qui_number = int(n1*quies_rate)
                neu_pass_number = neu_qui_number+neu_div_number
            random.shuffle(neu_list)
            neu_qui_list = neu_list[0:neu_qui_number]
            neu_div_list = neu_list[neu_qui_number:neu_pass_number]
            neu_div_list_double = neu_div_list*2
            n1_double = len(neu_div_list_double)

        if n2 > 0:
            adv_qui_number = lowerORupper(n2*(quies_rate-sfit))        #number of dividing cells in this generation
            adv_div_number = lowerORupper(n2*(birth_rate+sfit))
#            adv_qui_number = lowerORupper(n2*(0.9-0.1-(1+sfit)*(birth_rate-0.1)))
#            adv_div_number = lowerORupper(n2*(0.1+(1+sfit)*(birth_rate-0.1)))
            adv_pass_number = adv_qui_number+adv_div_number
            if adv_pass_number > n2:
                adv_qui_number = int(n2*(quies_rate-sfit))
                adv_pass_number = adv_qui_number+adv_div_number
            random.shuffle(adv_list)
            adv_qui_list = adv_list[0:adv_qui_number]
            adv_div_list = adv_list[adv_qui_number:adv_pass_number]
            adv_div_list_double = adv_div_list*2
            n2_double = len(adv_div_list_double)

        if n1_double > 0:
            new_mut1 = np.random.poisson(total_mut_rate*n1_double)
            mut_assig1 = Counter(np.random.choice(n1_double,new_mut1))
            for x1 in list(mut_assig1.keys()):
                nmut = mut_assig1[x1]
                new_mut1 = list(range(current_id+1,current_id+1+nmut))
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_div_list_double[x1]]
                    if random.random() < hotspot_chance:
                        barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
                        barcode_muts += [barcode_mut_pos]
                    else:
                        barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
                        barcode_muts +=  [barcode_mut_pos]
                neu_div_list_double[x1] = mut_str

        if n2_double > 0:
            new_mut2 = np.random.poisson(total_mut_rate*n2_double)
            mut_assig2 = Counter(np.random.choice(n2_double,new_mut2))
            for x2 in list(mut_assig2.keys()):
                nmut = mut_assig2[x2]
                new_mut2 = list(range(current_id+1,current_id+1+nmut))
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_div_list_double[x2]]
                    if random.random() < hotspot_chance:
                        barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
                        barcode_muts += [barcode_mut_pos]
                    else:
                        barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
                        barcode_muts +=  [barcode_mut_pos]
                adv_div_list_double[x2] = mut_str

        if random.random() < adv_rate*n1_double:
            current_driver_id += 1
            #current_n1 = len(neu_list)
            #lineage += [str(neu_list[current_n1-1])]
            adv_div_list_double += [neu_div_list_double[n1_double-1]]
            neu_div_list_double = neu_div_list_double[0:n1_double-1]

        neu_list=neu_qui_list+neu_div_list_double
        adv_list=adv_qui_list+adv_div_list_double
        current_deme_size=len(neu_list)+len(adv_list)

    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]

    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]

    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,current_driver_id,lineage,barcode_muts





######model adv low adv_rate = pow(10,-4)
def SMALT_models_low(parameter):
	'''
	The SMALT_models_low function is a simulation of a tumor growth process, in which each tumor cell is represented by a deme. This function simulates tumor growth from a single transformed cell to a tumor of a certain size, specified by the final_tumor_size variable.

	The function takes one argument, parameter, which is a dictionary containing several parameters that are used in the simulation. The first two parameters, mut_rate and birth_rate, are used to compute the time until the next cell division. The mut_rate parameter is set to the value of parameter["mu"] after being converted to a floating point number and raised to the power of 10. The birth_rate parameter is set to the value of parameter["birth_rate"] converted to a floating point number.

	Next, the function sets several other variables, including rd, which is the radius of the pre-created 3D space, and final_tumor_size, which is the number of cells in the final tumor. It also creates an empty lattice using the createLattice function.

	After setting up the simulation, the function then runs a loop until the number of demes in the tumor reaches the final number of demes, final_deme_number. In each iteration of the loop, the function chooses the deme with the shortest time until the next cell division, and updates the time remaining until the next cell division for all demes.

	Next, the function determines which of the deme's 26 neighboring sites are empty, and if there are any empty neighbor sites, it chooses one of them at random and creates a new deme there. It then updates the empty site counts for the neighboring demes and removes any demes that have no empty sites remaining from the list of demes on the surface.

	The function continues to simulate tumor growth until the number of demes in the tumor reaches the final number of demes. It then returns the final number of demes and the list of mutation lineages for each deme.
	'''
	mut_rate = - parameter["mu"]
	mut_rate = pow(10, mut_rate)
	#print(mut_rate)
	birth_rate = parameter["birth_rate"]
	s_coef = parameter["s_coef"]
	rd = 25                             # the radius of the pre-created 3D space
	final_tumor_size = pow(10,7)        # the number of cells in the final tumor
	final_deme_number = int(final_tumor_size/deme_size)    # the final number of demes in the tumor
	#death_rate = parameter["death_rate"]
	death_rate = 0.1
	quies_rate = 1 - death_rate - birth_rate
	non_hotspot_number = 1000           # the number of AID editing sites in each cell
	#hotspot_number = 20                 # the number of AID hotspot editing sites in each cell
	hotspot_number = 0                 # the number of AID hotspot editing sites in each cell
	hotspot_increase = 20
	total_mut_rate = mut_rate*(non_hotspot_number+hotspot_number*hotspot_increase)
	hotspot_chance = hotspot_number*hotspot_increase*1.0/(hotspot_number*hotspot_increase+non_hotspot_number)
	adv_rate = pow(10,-4)
	mut_id = 0                          # the ids of each mutation event, 1,2,3,4,5,...
	driver_mut_id = 0                   # the ids of advantageous driver mutations (modeling clonal selection)
	mutlineage = ['0']                  # the lineage tracer
	barcode_mut_site = [0]              # the mutation positions on general SMALT barcode sequence
######################################################################################
	first_neu,first_adv,mut_id,driver_mut_id,mutlineage,barcode_mut_site = initiateFirstDeme_t1(deme_size, mutlineage, mut_id, driver_mut_id, s_coef, barcode_mut_site, total_mut_rate, quies_rate, birth_rate, adv_rate, hotspot_chance, hotspot_number, non_hotspot_number)
	space = createLattice(rd)
	space[(rd,rd,rd)].present = 1
	space[(rd,rd,rd)].empty = 26
	space[(rd,rd,rd)].neutral = list(first_neu)
	space[(rd,rd,rd)].advant = list(first_adv)
	current_keys = [(rd,rd,rd)]
	current_deme_number = 1                                 #current deme number
	surface_keys = [(rd,rd,rd)]
	surface_times = [1.0]
	surface_deme_number = 1
	current_time = 0
	while current_deme_number < final_deme_number:
		deme_index = surface_times.index(min(surface_times))
		ckey = surface_keys[deme_index]

		ctime = surface_times[deme_index]
		current_time += ctime
		surface_times = [xt-ctime for xt in surface_times]
		#if current_deme_number%100 == 0:
		#	print(current_time,current_deme_number,surface_deme_number,mut_id,len(mutlineage)-1,len(barcode_mut_site)-1,driver_mut_id)

		nei_sites = neighbor26(ckey)  # neighbor sites of (rx,ry,rz)
		empty_sites = [key for key in nei_sites if space[key].present == 0]                    # the empty neighbor sites

		if len(empty_sites) > 0:
			num_neu = len(space[ckey].neutral)
			num_advant = len(space[ckey].advant)
			faction_advant = num_advant*1.0/(num_neu+num_advant)
			next_time = np.random.exponential(1.0/(1+faction_advant))
			#next_time = np.random.exponential(1.0)
			surface_times[deme_index] = next_time

			#rand_prob = random.random()
			#if rand_prob < 1-math.exp(-len(empty_sites)*0.25): # the probability for a deme to grow and divide is proportional to the # of empty neighbor sites
			pre_neu = list(space[ckey].neutral)
			pre_adv = list(space[ckey].advant)
			post_neu_l1,post_neu_l2,post_adv_l1,post_adv_l2,mut_id,driver_mut_id,mutlineage,barcode_mut_site = demeGrowthFission_t1(pre_neu,pre_adv,mutlineage,mut_id,driver_mut_id,current_deme_number,s_coef,barcode_mut_site, deme_size, quies_rate, birth_rate, total_mut_rate, hotspot_chance, hotspot_number, non_hotspot_number, adv_rate)
			space[ckey].neutral = list(post_neu_l1)
			space[ckey].advant = list(post_adv_l1)

			nkey = random.choice(empty_sites)
			space[nkey].neutral = list(post_neu_l2)
			space[nkey].advant = list(post_adv_l2)
			space[nkey].present = 1
			current_keys += [nkey]
			current_deme_number += 1

			next_nei_sites = neighbor26(nkey)
			next_empty_sites = []
			for key in next_nei_sites:
				if space[key].present == 1:
					space[key].empty -= 1
					if space[key].empty == 0:
						remove_index = surface_keys.index(key)
						del surface_keys[remove_index]
						del surface_times[remove_index]
						surface_deme_number -= 1
				else:
					next_empty_sites += [key]
			if len(next_empty_sites) > 0:
				num_neu = len(space[nkey].neutral)
				num_advant = len(space[nkey].advant)
				faction_advant = num_advant*1.0/(num_neu+num_advant)
				next_time = np.random.exponential(1.0/(1+faction_advant))
#				next_time = np.random.exponential(1.0)

				surface_keys += [nkey]
				surface_deme_number += 1
				surface_times += [next_time]

			space[nkey].empty = len(next_empty_sites)

			#if current_deme_number == met_timing:
				#    met_ancestor_list = list(post_neu_l2+post_adv_l2+post_adv2_l2)
				#    met_founder = random.choice(met_ancestor_list)

		else:
			break


	####visulization of spatial clonal structure in the central slice###
	sample_demes_n = int(right_round(ncells*10/100, 0))
	sample_demes = random.sample(current_keys, sample_demes_n)
	sample_deme_cells = []
	for key in sample_demes:
		cells0 = space[key].neutral+space[key].advant
		sample_deme_cells += cells0
	
	sample_cells = random.sample(sample_deme_cells, ncells)
	mut_array_pos = []
	barcodeseq = np.zeros([ncells, barcodelen])
	count = 0
	for k in sample_cells:
		mut_lin = traceLineage(mutlineage,k)
		mut_position = [barcode_mut_site[x] for x in mut_lin]
		mut_position = np.array(mut_position)
		mut_position = np.unique(mut_position)
		mut_position = mut_position.astype(int)
		barcodeseq[count][mut_position] = 1
		count += 1
		#print(len(mut_lin),len(mut_position),mut_lin,mut_position)
	barcodeseq = barcodeseq.astype(int)
	simultaion_sum_stat = seq2stat(barcodeseq)
	return {"X_2": simultaion_sum_stat}
#		np.append(mut_array_pos, mut_position, axis=0)

####define the two models

#######set the priors of models
parameter_prior = Distribution(
	mu = RV("uniform", 3, 1), birth_rate = RV("uniform", 0.15, 0.4), s_coef = RV("uniform", 0, 0.1)
	)

import argparse
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--nsim',type=int, default=5, help='the population size of each ABC-SMC iteration')
parser.add_argument('--npopulation',type=int, default=2, help='The maximum number of iterations as passed to ABCSMC')
parser.add_argument('--minMB',type=int, default=24, help='The minmum number of mutation burden of barcode sequence')
parser.add_argument('--demeSize',type=int, default=2000, help='The demes recapitulate the glandular structure often found in epithelial tumors such as colorectal cancer in which the gland size is approximated at 2,000- 10,000 cells')
parser.add_argument('--nprocs',type=int, default=80, help='the number of cores used in ABC-SMC')
parser.add_argument('--barcodelen',type=int, default=3004, help='the length of barcode sequence')
parser.add_argument('--alphavaule',type=float, default=0.6, help='The alpha-quantile to be used in ABC-AMC')
parser.add_argument('--inputRdata',type=str, help='The input data for ABC-SMC')
parser.add_argument('--dbpath',type=str, help='The path to directory storing the database')
args = parser.parse_args()
nsim = args.nsim
npopulation = args.npopulation
nprocs = args.nprocs
minMB = args.minMB
deme_size = args.demeSize
alphavaule = args.alphavaule
barcodelen = args.barcodelen
inputRdata = args.inputRdata
dbpath = args.dbpath
dataname = os.path.basename(inputRdata)

#####pre-process input data
obs_seq = readdata(inputRdata)
ncells = len(obs_seq)
sum_stat_obs = obs2stat(obs_seq)
processID = os.getpid()
####The ABC-SMC algorithm
abc = ABCSMC(SMALT_models_low, parameter_prior, distance, population_size = nsim,sampler=sampler.MulticoreEvalParallelSampler(nprocs), eps=QuantileEpsilon(alpha=alphavaule))
db_path = create_sqlite_db_id(file_= f"/data/xieduo/lineage_tracing/{dataname}_single_model_nsim={nsim}_npopulations={npopulation}_{processID}_nohotspot_no75.db")
abc.new(db_path, {"X_2": sum_stat_obs})
history = abc.run(max_nr_populations=npopulation)

