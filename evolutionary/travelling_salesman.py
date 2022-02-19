from evolutionary_algorithm import EvolutionaryAlgorithm, Cell, Population
import random

class World:
    cells = {}
    n_cities = 20
    cities = [(random.randint(1,50),random.randint(1,50)) for _ in range(n_cities)]

class TravellingSalesman(EvolutionaryAlgorithm):

    def get_distance(self, dna):
        citylist = []
        for gene in dna:
            citylist.append(World.cities[gene[0]])
        total_distance = 0
        for num, (x1, y1) in enumerate(citylist):
            try:
                x2 = citylist[num + 1][0]
                y2 = citylist[num + 1][1]
            except:
                x2 = citylist[0][0]
                y2 = citylist[0][1]
            distance = (((x2 - x1) ** 2) + ((y2 - y1) ** 2)) ** 0.5
            total_distance += distance
        return total_distance

    def fitness_func(self, dna):
        total_distance = self.get_distance(dna)
        return 1 / total_distance

    def gene_func(self, n_genes, n_chroms):
        dna = [[] for _ in range(n_chroms)]
        cities = list(range(n_chroms))
        random.shuffle(cities)
        for _ in range(n_genes):
            for n in range(len(dna)):
                dna[n].append(cities[n])
        return dna

    def mate_func(self, cell1, cell2):
        parent1 = cell1.dna
        parent2 = cell2.dna

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
                mutation_idx.append((random.randint(0, len(child_dna) - 1), random.randint(0, len(child_dna) - 1)))
            mutated_dna = []
            for y, chrom in enumerate(child_dna):
                mutated_dna.append([])
                for x, gene in enumerate(chrom):
                    mutated_dna[y].append(gene)
            for midx1 in mutation_idx:
                mutated_dna[midx1[0]], mutated_dna[midx1[1]] = mutated_dna[midx1[1]], mutated_dna[midx1[0]]
            child_dna = mutated_dna
        child_cell = Cell(dna=child_dna, gene_func=self.gene_func)
        return child_cell

if __name__ == "__main__":
    n_cells_per_pop = 50
    n_pops = 4
    n_genes = 1
    n_chroms = len(World.cities)
    n_epochs = 200
    verbose = True
    n_mutations = 1
    n_mutation_chance = 0.1

    evolutionary_algorithm = TravellingSalesman(n_cells_per_pop, n_pops, n_genes, n_chroms, n_epochs, verbose, n_mutations, n_mutation_chance)
    evolutionary_algorithm.run()
    evolutionary_algorithm.plot_stats()

    print(evolutionary_algorithm.get_distance(evolutionary_algorithm.selection(n_selected=1)[0].dna))
