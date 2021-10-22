import random
from tqdm import tqdm
import matplotlib.pyplot as plt
#plt.style.use('fivethirtyeight')

class World:
    cell_dict = {}
    n_cities = 20
    cities = [(random.randint(1,200),random.randint(1,200)) for _ in range(n_cities)]

class Cell:
    """
    Organism, that has a DNA.
    Rewritten for the Travelling Salesman Problem
    """
    number = 0
    def __init__(self, dna=None, n_chromosomes = 6, n_genes = 1):
        self.number = Cell.number
        Cell.number += 1

        self.mutation_chance = 0.1
        self.n_mutations = 1 # In this problem, a mutation represents a switch of cities

        if dna is not None:
            self.dna = dna
        else:
            self.dna = self._random_dna(n_chromosomes, n_genes)

        World.cell_dict[self.number] = self

    def _random_dna(self, n_chromosomes=6, n_genes=1):
        """
        Generates a random DNA for Cell.
        Rewritten for the Travelling Salesman problem.
        :param n_chromosomes: Number of chromosomes.
        :param n_genes: Number of genes inside 1 Chromosome.
        :return: random DNA; example: [[2],[1],[3],[0],[4]]
        """
        dna = [[] for _ in range(n_chromosomes)]
        cities = list(range(n_chromosomes))
        random.shuffle(cities)
        for _ in range(n_genes):
            for n in range(len(dna)):
                dna[n].append(cities[n])
        return dna

    def mate(self, partner):
        """
        Produces a new Cell individual with partner.
        :param partner: type Cell, partner with which to produce child.
        :return: DNA of the produced child. (type list)
        """
        parent1 = self.dna
        parent2 = partner.dna

        # --- Crossover
        geneA = int(random.random() * len(parent1))
        geneB = int(random.random() * len(parent1))

        startGene = min(geneA, geneB)
        endGene = max(geneA, geneB)

        childP1 = []
        for i in range(startGene, endGene):
            childP1.append(parent1[i])

        childP2 = [item for item in parent2 if item not in childP1]
        child_dna = childP2[:startGene] + childP1 + childP2[startGene:]

        # --- Mutation
        if random.random() < self.mutation_chance:
            mutation_idx = []
            for _ in range(self.n_mutations):
                mutation_idx.append((random.randint(0,len(child_dna)-1),random.randint(0,len(child_dna)-1)))
            mutated_dna = []
            for y, chrom in enumerate(child_dna):
                mutated_dna.append([])
                for x, gene in enumerate(chrom):
                    mutated_dna[y].append(gene)
            for midx1 in mutation_idx:
                mutated_dna[midx1[0]],mutated_dna[midx1[1]] = mutated_dna[midx1[1]],mutated_dna[midx1[0]]
            child_dna = mutated_dna
        child_cell = Cell(dna=child_dna)
        return child_cell

    def _get_fitness(self):
        """
        Calculates the fitness of this cell.
        Rewritten for the Travelling Salesman problem.
        :return: Fitness
        """
        total_distance = self._get_distance()
        return 1 / total_distance

    def _get_distance(self):
        """
        Calculates the fitness of this cell.
        Rewritten for the Travelling Salesman problem.
        :return: Fitness
        """
        # ((((x2 - x1 )**2) + ((y2-y1)**2) )**0.5)
        citylist = []
        for gene in self.dna:
            citylist.append(World.cities[gene[0]])
        total_distance = 0
        for num, (x1,y1) in enumerate(citylist):
            try:
                x2 = citylist[num+1][0]
                y2 = citylist[num+1][1]
            except:
                x2 = citylist[0][0]
                y2 = citylist[0][1]
            distance = (((x2 - x1 )**2) + ((y2-y1)**2))**0.5
            total_distance += distance
        return total_distance

    def _mean_normalize(self, data):
        """
        Applies mean normalization to a single dimension list of values.
        :param data: type list, filled with values
        :return: returns list with all values normalized (0-1)
        """
        newdata = []
        for d in data:
            norm = 1 / (max(data) - min(data)) * (d - max(data)) + 1 #(max'-min') / (max - min) * (value - max) + max'
            newdata.append(norm)
        return newdata

    def get_similarity(self, other_cell):
        """
        Function to calculate similarity of DNAs.
        :param other_cell: type Cell, to compare DNAs with.
        :return: Similarity of the two cells in percent.
        """
        dna1 = self.dna
        dna2 = other_cell.dna
        max_similarity = len([j for sub in dna1 for j in sub])
        n_similar = 0
        for num, chromosome1 in enumerate([j for sub in dna1 for j in sub]):
            chromosome2 = [j for sub in dna2 for j in sub][num]
            if chromosome1 == chromosome2:
                n_similar += 1
        if n_similar == 0:
            return 0
        similarity = n_similar / max_similarity
        return similarity

    def __repr__(self):
        return str(self.dna)

class Population:
    """
    Main class for controlling the genetic algorithm.
    """
    def __init__(self, n_population=50, n_generations=100):
        self.n_population = n_population
        self.n_generations = n_generations

        self.population = []
        self.fitness_plot = []
        self.distance_plot = []

    def selection(self):
        """
        Uses a selection process to select the next n_selected parents for further breeding.
        :return: The selected Cells. type list.
        """
        selection = []

        fitness_list = [c._get_fitness() for c in self.population]

        best = self.population[fitness_list.index(max(fitness_list))]
        selection.append(best)
        fitness_list.remove(max(fitness_list))
        best = self.population[fitness_list.index(max(fitness_list))]
        selection.append(best)

        return selection

    def next_generation(self):
        """
        Breeds next generation.
        :return: This generation's population, The fitness plot up to this point
        """
        selection = self.selection()

        self.population = []
        for pop_idx in range(self.n_population):
            self.population.append(selection[0].mate(selection[1]))

        # --- Statistics ---
        fitness_list = [c._get_fitness() for c in self.population]
        distance_list = [c._get_distance() for c in self.population]
        self.fitness_plot.append(sum(fitness_list) / len(fitness_list))
        self.distance_plot.append(sum(distance_list) / len(distance_list))

        return self.population, self.fitness_plot, self.distance_plot

    def plot_stats(self):
        """
        Plots the fitness and distance over all generations using matplotlib.
        """
        # --- Plotting the distance ---
        plt.plot(self.distance_plot)
        plt.title('Average Distance')
        plt.ylabel('Distance')
        plt.xlabel('Generations')
        plt.legend(['Distance'], loc='lower right')
        plt.show()

        # --- Plotting the fitness ---
        plt.plot(self.fitness_plot)
        plt.title('Average Fitness')
        plt.ylabel('Fitness')
        plt.xlabel('Generations')
        plt.legend(['Fitness'], loc='lower right')
        plt.show()

    def run(self):
        """
        Automatically runs the algorithm.
        :return: New population, Fitness data
        """
        for n in tqdm(range(self.n_generations)):
            #print(f"{n}/{self.n_generations}")
            self.next_generation()
        return self.population, self.fitness_plot, self.distance_plot

if __name__ == "__main__":
    # --- Adam and Eve
    adam = Cell(n_chromosomes=len(World.cities), n_genes=1)
    eve = Cell(n_chromosomes=len(World.cities), n_genes=1)

    pop = Population(n_population=100, n_generations=10000)
    pop.population.append(adam)
    pop.population.append(eve)
    pop.run()

    pop.plot_stats()
