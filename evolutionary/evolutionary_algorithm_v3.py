import random
import matplotlib.pyplot as plt
import gc

class World:
    cells = {}
    populations = {}

class Cell:
    number = 0
    def __init__(self, gene_func, dna=None, n_genes=4, n_chroms=4):
        self.number = Cell.number
        Cell.number += 1

        self.random_dna = gene_func
        self.meta = {}
        self.cached_fitness = None

        self.dna = dna
        if self.dna is None:
            self.dna = self.random_dna(n_genes, n_chroms)

        World.cells[self.number] = self

class Population:
    number = 0
    def __init__(self, gene_func, pop=None, n_cells=50, n_genes=4, n_chroms=4):
        self.number = Population.number
        Population.number += 1

        self.population = pop
        if self.population is None:
            self.population = [Cell(n_genes=n_genes, n_chroms=n_chroms, gene_func=gene_func) for cell in range(n_cells)]

        World.populations[self.number] = self

class EvolutionaryAlgorithm:
    def __init__(self, n_cells_per_pop, n_pops, n_genes, n_chroms, n_epochs=200, verbose=True, n_mutations=1, mutation_chance=0.1):
        self.populations = [Population(n_cells=n_cells_per_pop, n_genes=n_genes, n_chroms=n_chroms, gene_func=self.gene_generator) for pop in range(n_pops)]
        self.n_epochs = n_epochs

        self.verbose = verbose
        self.fitness_history = []
        self.n_mutations = n_mutations
        self.mutation_chance = mutation_chance

        n_crossovers = n_pops-1
        self.crossover_points = [n_epochs//(n_crossovers+1)*(n+1) for n in range(n_crossovers)]

    def gene_generator(self, n_genes, n_chroms):
        dna = [[random.random() for n in range(n_chroms)] for _ in range(n_genes)]
        return dna

    def selection(self, n_selected):
        """
        Uses a selection process to select the next n_selected parents for further breeding.
        :return: The selected Cells. type list. selection(2) returns [Cell, Cell]
        """
        selection = []

        fitness_list = [[(c.number, self.fitness(c)) for c in p.population] for p in self.populations]
        fitness_list = [j for sub in fitness_list for j in sub]
        sorted_values = sorted(fitness_list, key=lambda x: x[1], reverse=True)
        for n in range(n_selected):
            selection.append(World.cells[sorted_values[n][0]])

        return selection

    def selection_from_pop(self, n_selected, pop):
        """
        Uses a selection process to select the next n_selected parents for further breeding.
        :return: The selected Cells. type list. selection(2) returns [Cell, Cell]
        """
        selection = []

        fitness_list = [(c.number, self.fitness(c)) for c in pop]
        sorted_values = sorted(fitness_list, key=lambda x: x[1], reverse=True)
        for n in range(n_selected):
            selection.append(World.cells[sorted_values[n][0]])

        return selection

    def advanced_selection(self, n_selected, pop):
        """
        Uses a advanced selection process to select the parents with the best outcome for the next generation.
        :param n_selected: Placeholder if someone used old function with same parameters.
        :return: The selected Cells. type list.
        """
        selection = []
        n_preselected = random.randint(2, 10)
        fitness_list = [(c.number, self.fitness(c)) for c in pop]
        sorted_values = sorted(fitness_list, key=lambda x: x[1], reverse=True)
        preselected = [y[0] for y in sorted_values[:n_preselected]]
        combinations = {}
        for p1 in preselected:
            p1_cell = World.cells[p1]
            for p2 in preselected:
                p2_cell = World.cells[p2]
                child = self.mate(p1_cell,p2_cell)
                combinations[f"{p1}-{p2}"] = self.fitness(child)
                del World.cells[child.number]
                del child
        gc.collect()
        sorted_values = sorted(combinations.items(), key=lambda x: x[1], reverse=True)
        best = sorted_values[0][0].split("-")
        selection.append(World.cells[int(best[0])])
        selection.append(World.cells[int(best[1])])
        return selection

    def fitness(self, cell):
        if cell.cached_fitness is not None:
            return cell.cached_fitness
        single_dim = [j for sub in cell.dna for j in sub]
        fit = sum(single_dim)
        cell.cached_fitness = fit
        return fit

    def decimator(self, pop):
        return random.choice(pop)

    def mate(self, cell1, cell2):
        dna1 = cell1.dna
        dna2 = cell2.dna

        # --- Crossover
        child_dna = self.crossover(dna1, dna2)

        # --- Mutation
        if random.random() < self.mutation_chance:
            child_dna = self.mutation(child_dna)
        child_cell = Cell(dna=child_dna, gene_func=self.gene_generator)
        self.on_child(child_cell)
        return child_cell

    def crossover(self, dna1, dna2):
        crossover_point = random.randint(0, len(dna1) - 1)
        child_dna = dna1[:crossover_point] + dna2[crossover_point:]
        return child_dna

    def mutation(self, child_dna):

        mutation_idx = []
        for _ in range(self.n_mutations):
            mutation_idx.append(random.randint(0, len(child_dna) - 1))

        mutated_dna = [g if num not in mutation_idx else [random.random() for _ in range(len(child_dna[0]))] for num, g in enumerate(child_dna)]

        return mutated_dna

    def mix_populations(self, pop1, pop2):
        half = len(pop1.population)//2
        random.shuffle(pop1.population)
        random.shuffle(pop2.population)
        newpop = pop1.population[:half]+pop2.population[half:]
        self.populations.append(Population(pop=newpop, gene_func=self.gene_generator))

    def log(self, txt, end="\n"):
        if self.verbose:
            print(txt, end=end)

    def plot_stats(self):
        plt.plot(self.fitness_history)

        plt.plot([n for n in self.crossover_points], [self.fitness_history[n] for n in self.crossover_points], "^")

        plt.title('Average Fitness')
        plt.ylabel('Fitness')
        plt.xlabel('Generations')
        plt.legend(['Fitness'], loc='lower right')
        plt.show()

    def run(self):
        for epoch in range(self.n_epochs):
            self.log(f"Epoch {epoch}/{self.n_epochs-1} | Average Fitness: ", end="")

            self.on_epoch(epoch)

            for num, population in enumerate(self.populations):
                selection = self.advanced_selection(n_selected=2, pop=population.population)
                n_children = len(population.population)//2
                for c in range(n_children):
                    to_die = self.decimator(population.population)
                    del self.populations[num].population[self.populations[num].population.index(to_die)]
                    del World.cells[to_die.number]
                    del to_die
                    self.populations[num].population.append(self.mate(selection[0],selection[1]))
            gc.collect() # collect garbage (deleted stuff)

            # --- Statistics ---
            fitness_list = [[self.fitness(c) for c in p.population] for p in self.populations]
            fitness_list = [j for sub in fitness_list for j in sub]
            avg_fitness = sum(fitness_list) / len(fitness_list)
            self.log(avg_fitness)
            self.fitness_history.append(avg_fitness)

            if epoch in self.crossover_points:
                self.log("Crossing Populations")
                n_crosses = 2
                to_delete = []
                pop_crosses = random.sample(self.populations, n_crosses)
                while True:
                    if len(pop_crosses) > 1:
                        to_delete.append(self.populations[self.populations.index(pop_crosses[0])])
                        to_delete.append(self.populations[self.populations.index(pop_crosses[-1])])
                        self.mix_populations(self.populations[self.populations.index(pop_crosses[0])], self.populations[self.populations.index(pop_crosses[-1])])
                        del pop_crosses[0]
                        del pop_crosses[-1]
                    else:
                        break
                for d in to_delete:
                    del self.populations[self.populations.index(d)]

    def on_epoch(self, epoch):
        pass

    def on_child(self, child):
        pass

if __name__ == "__main__":
    n_cells_per_pop = 50
    n_pops = 64
    n_genes = 5
    n_chroms = 5
    n_epochs = 200
    verbose = True
    n_mutations = 2
    n_mutation_chance = 0.2

    evolutionary_algorithm = EvolutionaryAlgorithm(n_cells_per_pop, n_pops, n_genes, n_chroms, n_epochs, verbose, n_mutations, n_mutation_chance)
    evolutionary_algorithm.run()
    evolutionary_algorithm.plot_stats()
