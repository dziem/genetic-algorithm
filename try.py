import numpy as np
import random
import math

x1 = [-100,100]
x2 = [-100,100]
pops_num = 100
generatione = 10000 // pops_num #10k is max total pops
selected_parent_num = 20 #and also child amount, same amount
cross_probs = 0.9
mutate_probs = 0.1
alfha = 0.5 #for crossover
cee = 2 #for scaling
'''
def functione(x, y):
	a1 = 0
	a2 = 0
	for i in range(1,6):
		a1 += i * math.cos(math.radians((i + 1) * x + 1))
		a2 += i * math.cos(math.radians((i + 1) * y + 1))
	a1 = a1 * -1
	return a1 * a2
'''
def functione(x, y):
	a = -1 * math.cos(math.radians(x))
	b = math.cos(math.radians(y))
	c = math.exp((-1 * ((x - math.pi) ** 2)) - ((y - math.pi) ** 2))
	return a * b * c

def generateParent(pops_num):
	pops = [] #whole populatione
	for i in range(pops_num):
		chromosome = [] #chromosome
		chromosome.append(random.uniform(0, 1)) #x1
		chromosome.append(random.uniform(0, 1)) #x2
		pops.append(chromosome)
	return pops

def decodeChromosome(value, x): #value, range min max
	return (value * (abs(x[0]) + abs(x[1]))) - abs(x[0])

def fitness(value):
	return 2 ** (-1 * value)

def chromosomeFitness(pops, x1, x2):
	fitnesse = []
	for i in range(len(pops)):
		chromosome_fitness = []
		chromosome_fitness.append(i) #chromosome index
		chromosome_fitness.append(fitness(functione(decodeChromosome(pops[i][0], x1), decodeChromosome(pops[i][1],x2)))) #chromosome fitness
		fitnesse.append(chromosome_fitness)
	return fitnesse

def selectingParent(fits, num_parent, pops, cee):
	selected_parent = []
	#sigma scaling (pizdec)
	'''
	meanArr = np.mean(fits, axis = 0)
	stdevArr = np.std(fits, axis = 0)
	mean = meanArr[1]
	stdev = stdevArr[1]
	for i in fits:
		i[1] = i[1] + (mean - cee * stdev)
	'''
	#window scaling
	mins = min(x[1] for x in fits)
	for i in fits:
		i[1] = (i[1] - mins) + (mins * 0.0000000000001) #plus so that the min list ain't 0
	#fits.sort(key=lambda x: x[1], reverse=True)
	#fits = fits[:len(fits)-(len(pops) - num_parent)]
	random.shuffle(fits)
	sum = 0
	for i in fits:
		sum += i[1]
	selector = np.linspace(0, sum, num_parent)
	for i in selector:
		current = 0
		for j in fits:
			current += j[1]
			if current >= i:
				selected_parent.append(pops[j[0]])
				break
	return selected_parent

def crossOver(parents, cross_probs, alfha, mutate_probs):
	child = []
	for i in range(0,len(parents),2):
		z = random.uniform(0, 1)
		if z >= (1 - cross_probs): #whole arithmetic crossover
			child_chromosome = []
			newx = (alfha * parents[i][0]) + ((1 - alfha) * parents[i + 1][0])
			newy = (alfha * parents[i][1]) + ((1 - alfha) * parents[i + 1][1])
			child_chromosome.append(newx)
			child_chromosome.append(newy)
			child.append(mutatione(child_chromosome, mutate_probs))
			child.append(mutatione(child_chromosome, mutate_probs))
		else:
			mutatione(parents[i], mutate_probs)
			mutatione(parents[i + 1], mutate_probs)
			child.append(parents[i])	
			child.append(parents[i + 1])	
	return child
	
def mutatione(child, mutate_probs):
	z = random.uniform(0, 1)
	if z <= cross_probs: #creep
		xory = random.choice([0, 1])
		step = random.choice([-0.01, 0.01])
		if(step == -0.01):
			if(child[xory] <= 0.01):
				child[xory] = 0
			else:
				child[xory] += step
		elif(step == 0.01):
			if(child[xory] >= 0.99):
				child[xory] = 1
			else:
				child[xory] += step
	return child

def replacePops(pops, child, x1, x2):
	pops_fitness = chromosomeFitness(pops, x1, x2)
	pops_fitness.sort(key=lambda x: x[1], reverse=True)
	pops_fitness = pops_fitness[:len(pops_fitness)-(len(pops)-len(child))]
	pops_fitness.sort(key=lambda x: x[0], reverse=True)
	for i in pops_fitness:
		pops.pop(i[0])
	for i in child:
		pops.append(i)
	return pops

def bestIndividue(pops, x1, x2):
	pops_fitness = chromosomeFitness(pops, x1, x2)
	pops_fitness.sort(key=lambda x: x[1], reverse=True)
	return pops_fitness[0][0]

pops = generateParent(pops_num)
for i in range(generatione):
	fits = chromosomeFitness(pops, x1, x2)
	selected_parente = selectingParent(fits, selected_parent_num, pops, cee)
	child = crossOver(selected_parente, cross_probs, alfha, mutate_probs)
	replacePops(pops, child, x1, x2)
best = bestIndividue(pops, x1, x2)
print(decodeChromosome(pops[best][0], x1))
print(decodeChromosome(pops[best][1], x2))
print(functione(decodeChromosome(pops[best][0], x1), decodeChromosome(pops[best][1],x2)))