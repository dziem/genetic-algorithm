import numpy as np
import random
import math

x1 = [-100,100] #x1 value range, -100 t0 100
x2 = [-100,100] #x2 value range, -100 to 100
pops_num = 100 #amount of population per generation
generatione = 10000 // pops_num #the maximum amount of generation, 10k is max total pops
selected_parent_num = 20 #and also child amount, same amount
cross_probs = 0.9 #crossover probability
mutate_probs = 0.9 #mutation probability
alfha = 0.5 #for crossover calculation
t_size = 40 #tournament size, for parent selection
best_probs = 0.33 #probability to choose tournament winner = p * (1 - p)^n , n = 0 ... t_size where 0 is best and so on
creep_step = 0.001 #for creeping, -x or +x
stop_after = 25 #stop running after x generation without improvement
run_times = 30 #x times running the GA, for best solution, worst solution, and average
#changed params = original params, best
#selected_parent_num = 20
'''
cos note
The value passed in this function should be in radians.
https://www.geeksforgeeks.org/python-math-cos-function/
and also there's no further instruction in the assignment instruction document
sooo all the cos use radian, for example cos(x) -> math.cos(math.radians(x))
'''
'''
change formula/function note
just comment one functione(x,y) and leave the other one
'''

#function 1
def functione(x, y):
	a1 = 0
	a2 = 0
	for i in range(1,6):
		a1 += i * math.cos(math.radians((i + 1) * x + 1))
		a2 += i * math.cos(math.radians((i + 1) * y + 1))
	a1 = a1 * -1
	return a1 * a2
'''
#function 2
def functione(x, y):
	a = -1 * math.cos(math.radians(x))
	b = math.cos(math.radians(y))
	c = math.exp((-1 * ((x - math.pi) ** 2)) - ((y - math.pi) ** 2))
	return a * b * c
'''
def generateParent(pops_num):
	#function to generate parent, params the amount of generated individual/parent, ex. 100
	#data structure -> pops[x] = [x1 representaion value, x2 representaion value], where x is parent number x
	pops = [] #population container
	for i in range(pops_num):
		chromosome = [] #chromosome container
		chromosome.append(random.uniform(0, 1)) #random real 0 - 1 for x1
		chromosome.append(random.uniform(0, 1)) #random real 0 - 1 for x2
		pops.append(chromosome)
	return pops

def decodeChromosome(value, x): 
	#function to decode chromosome, params chromosome value, value range, ex (0.5, [-100,100])
	#just one bit, x1 or x2, not both
	return (value * (abs(x[0]) + abs(x[1]))) - abs(x[0])

def fitness(value):
	#function to find the fitness value, params function value, ex 26
	return 2 ** (-1 * value)

def chromosomeFitness(pops, x1, x2):
	#function to find fitness value of each individual, params population, value range x1 & x2, ex. (pops,[-100,100],[-100,100])
	#data structure -> fitnesse[x] = [index in pops, fitness value]
	fitnesse = [] #fitness container
	for i in range(len(pops)):
		chromosome_fitness = [] #fitness item container
		chromosome_fitness.append(i) #individual index in pops
		chromosome_fitness.append(fitness(functione(decodeChromosome(pops[i][0], x1), decodeChromosome(pops[i][1],x2)))) #individual fitness
		fitnesse.append(chromosome_fitness)
	return fitnesse

def selectingParent(fits, num_parent, pops, best_probs, t_size):
	#function to selecte parent, params all_fitness_value, amount of parent, population, best winning probability, tournament size
	#params ex. (fitnesse, 20, pops, 0.5, 32) 
	selected_parent = [] #selected parents container
	for i in range(num_parent):
		tournament_participant = [] #tournament participante container, reset every tournament
		picked = [] #picked participant container, so that an individual can't have 2 participant in the same container, reset every tournament
		selected = False #to check if a tournament winner is selected
		for j in range(t_size):
			while True: #to avoid a double pick on an individual
				p = random.randrange(len(pops)) #pick an individual
				if p not in picked: #if individual not in tournament
					picked.append(p) #save picked individual
					break
			tournament_participant.append(fits[p]) #add tournament participant
		rand_num = random.uniform(0, 1) #pick random number to pick the winner
		tournament_participant.sort(key=lambda x: x[1], reverse=True) #sort the participant based on fitness value, sort descending
		prev = 0
		for j in range(t_size): #pick winner based on probability
			if rand_num >= prev and rand_num <= (1 - (best_probs * ((1 - best_probs) ** j))): #if and else if for probability winning p * (1 - p)^n, in loop so that it's dynamic
				selected_parent.append(pops[tournament_participant[j][0]]) #selected parent/winner
				selected = True
				break
			prev = 1 - (best_probs * ((1 - best_probs) ** j))
		if selected == False: #the else from probability p*(1 - p)^n
			selected_parent.append(pops[tournament_participant[t_size - 1][0]]) #selected parent/winner
	return selected_parent

def crossOver(parents, cross_probs, alfha, mutate_probs, creep_step): #also calling mutation
	#crossover and "mutation" function, params selected parents, crossover probability, alpha for crossover calculation, mutation probability, creep step (for mutation)
	#params ex. (parents, 0.9, 0.5, 0.1, 0.001)
	#crossover using whole arithmetic crossover
	child = [] #child container, same data structure with parent/individual
	#the pairing/marriage is based on position, ex. index 0 pair with index 0 + 1, index 2 and 2 + 1, and so on
	for i in range(0,len(parents),2): #loop 0 to parent num step 2
		z = random.uniform(0, 1) #random number
		if z >= (1 - cross_probs): #if true crossover
			child1 = [] #alpha * parent1 + (1 - alpha) * parent2
			child1.append((alfha * parents[i][0]) + ((1 - alfha) * parents[i + 1][0])) #bit x1
			child1.append((alfha * parents[i][1]) + ((1 - alfha) * parents[i + 1][1])) #bit x2
			child2 = [] #alpha * parent2 + (1 - alpha) * parent1
			child2.append((alfha * parents[i + 1][0]) + ((1 - alfha) * parents[i][0])) #bit x1
			child2.append((alfha * parents[i + 1][1]) + ((1 - alfha) * parents[i][1])) #bit x2
			child.append(mutatione(child1, mutate_probs, creep_step)) #mutate child
			child.append(mutatione(child2, mutate_probs, creep_step)) #mutate child
		else: #else no crossover
			mutatione(parents[i], mutate_probs, creep_step) #but mutate parent
			mutatione(parents[i + 1], mutate_probs, creep_step) #but mutate parent
			child.append(parents[i])	
			child.append(parents[i + 1])	
	return child
	
def mutatione(child, mutate_probs, creep_step):
	#mutation function, params individual, mutation probability, creep step, ex (child, 0.1, 0.001)
	#mutation using creeping, creeping on one bit, either x1 or x2, creeping by + or -
	z = random.uniform(0, 1) #random number
	if z <= mutate_probs: #if true mutate, if false the no mutation
		xory = random.choice([0, 1]) #pick creep x1 or x2
		step = random.choice([(-1 * creep_step), creep_step]) #pick + or -
		if(step == (-1 * creep_step)):
			if(child[xory] <= creep_step):
				child[xory] = 0 #so that the creep result ain't out of 0 - 1
			else:
				child[xory] += step
		elif(step == creep_step):
			if(child[xory] >= (1 - creep_step)):
				child[xory] = 1 #so that the creep result ain't out of 0 - 1
			else:
				child[xory] += step
	return child

def replacePops(pops, child, x1, x2):
	#function to replace population with new child, params population, children, value range x1 and x2
	#params ex. (pops, children, [-100,100], [-100,100])
	pops_fitness = chromosomeFitness(pops, x1, x2) #find chromosome fitness for current population
	pops_fitness.sort(key=lambda x: x[1]) #sort fitness value descending
	pops_fitness = pops_fitness[:len(pops_fitness)-(len(pops)-len(child))] #get worst fitness value, same amount with child amount
	pops_fitness.sort(key=lambda x: x[0], reverse=True) #sort the fitness index descending, to avoid index error/crash
	for i in pops_fitness:
		pops.pop(i[0]) #remove population based on worst index
	for i in child:
		pops.append(i) #add child to population
	return pops

def bestIndividue(pops, x1, x2):
	#function to find best individual, params population, value range x1 and x2, ex (pops, [-100,100], [-100,100])
	pops_fitness = chromosomeFitness(pops, x1, x2) #find chromosome fitness for current population
	pops_fitness.sort(key=lambda x: x[1], reverse=True) #sort fitness value descending
	return pops_fitness[0][0] #return the population index of the first (or best) fitness value

sum = 0 #for average
for z in range(run_times): #run GA run_times time
	pops = generateParent(pops_num) #generate parent
	solution = pops[bestIndividue(pops, x1, x2)] #initiate solution
	stop_counter = 0 #counter to stop the GA
	for i in range(generatione): #loop n-generation times
		fits = chromosomeFitness(pops, x1, x2) #find whole population fitness
		selected_parente = selectingParent(fits, selected_parent_num, pops, best_probs, t_size) #select parent
		child = crossOver(selected_parente, cross_probs, alfha, mutate_probs, creep_step) #crossover & mutation
		replacePops(pops, child, x1, x2) #replace population
		best = pops[bestIndividue(pops, x1, x2)] #get best individual in population
		if (functione(decodeChromosome(best[0], x1), decodeChromosome(best[1],x2)) < functione(decodeChromosome(solution[0], x1), decodeChromosome(solution[1],x2))):
			#if best individual better than best so far
			stop_counter = 0 #reset counter
			solution = best #save individual
		else: #if not better
			stop_counter += 1 #count to stop
		if stop_counter >= stop_after: #if time to stop
			break
	if z == 0: #initiate best and worst individual in GA run
		best_solution = functione(decodeChromosome(solution[0], x1), decodeChromosome(solution[1],x2))
		worst_solution = functione(decodeChromosome(solution[0], x1), decodeChromosome(solution[1],x2))
		best_chromosome = solution
		worst_chromosome = solution
		sum += best_solution
	else: #replace best and worst individual in GA run with best and worst from another GA run
		the_solution = functione(decodeChromosome(solution[0], x1), decodeChromosome(solution[1],x2))
		sum += the_solution
		if the_solution < best_solution:
			best_solution = the_solution
			best_chromosome = solution
		if the_solution > worst_solution:
			worst_solution = the_solution
			worst_chromosome = solution
print('best solution')
print('x1 =',decodeChromosome(best_chromosome[0], x1) , '& x2 =' , decodeChromosome(best_chromosome[1], x2))
print('f(x1,x2) =' , best_solution)
print('worst solution')
print('x1 =',decodeChromosome(worst_chromosome[0], x1) , '& x2 =' , decodeChromosome(worst_chromosome[1], x2))
print('f(x1,x2) =' , worst_solution)
print('average f(x1, x2) =' , (sum / run_times))