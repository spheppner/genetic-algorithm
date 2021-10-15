import random

class World:
    cell_dict = {}

class DNA:
    """
    DNA(n_chromosomes=4, n_genes=1):
    dna = [[0.1],[0.6],[0.2],[0.3]]

    DNA(n_chromosomes=5, n_genes=2):
    dna = [[0.1,0.2],[0.4,0.9],[0.8,0.5],[0.2,0.1],[0.8,0.3]]
    """
    def __init__(self, dna=None, n_chromosomes=4, n_genes=1):
        if dna is None:
            self.dna = self._random_dna(n_chromosomes, n_genes)
        else:
            self.dna = dna
        self.n_chromosomes = len(self.dna)
        self.n_genes = len(self.dna[0])

    def _random_dna(self, n_chromosomes, n_genes):
        dna = [[] for _ in range(n_chromosomes)]
        for _ in range(n_genes):
            for n in range(len(dna)):
                dna[n].append(random.random())
        return dna

    def __repr__(self):
        return str(self.dna)

class Cell:
    """
    Organism, that has a DNA.
    For understanding how n_chromosomes and n_genes work, see DNA() class above.
    """
    number = 0
    def __init__(self, dna=None, n_chromosomes = 3, n_genes = 3):
        self.number = Cell.number
        Cell.number += 1

        self.mutation_chance = 0.2
        self.n_mutations = 2

        if type(dna) == DNA:
            self.dna = dna.dna
        else:
            self.dna = DNA(dna, n_chromosomes, n_genes).dna

        World.cell_dict[self.number] = self

    def mate(self, partner):
        dna1 = self.dna
        dna2 = partner.dna

        # --- Crossover
        child_dna = []
        for n in range(len(dna1)):
            parent_chromosome = random.randint(0,1)
            if parent_chromosome == 0:
                child_dna.append(dna1[n])
            elif parent_chromosome == 1:
                child_dna.append(dna2[n])
        child_dna = DNA(dna=child_dna)

        # --- Mutation
        if random.random() < self.mutation_chance:
            mutation_idx1 = []
            mutation_idx2 = []
            for _ in range(self.n_mutations):
                mutation_idx1.append(random.randint(0,child_dna.n_chromosomes-1))
                mutation_idx2.append(random.randint(0,child_dna.n_genes-1))
            mutated_dna = []
            for y, chrom in enumerate(child_dna.dna):
                mutated_dna.append([])
                for x, gene in enumerate(chrom):
                    for id, idx in enumerate(mutation_idx1):
                        ok = True
                        if y == idx and x == mutation_idx2[id]:
                            mutated_dna[y].append(random.random())
                            ok = False
                            break
                    if ok:
                        mutated_dna[y].append(gene)
            child_dna = mutated_dna
        child_cell = Cell(dna=child_dna)
        return child_cell

    def get_similarity(self, other_cell):
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


if __name__ == "__main__":
    # --- Adam and Eve
    adam = Cell(n_chromosomes=5, n_genes=2)
    eve = Cell(n_chromosomes=5, n_genes=2)

    child = adam.mate(eve)
    print(child)