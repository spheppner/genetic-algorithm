from evolutionary_algorithm_v2 import EvolutionaryAlgorithm, Cell
import random
from pickle import dump
from PIL import Image, ImageDraw, ImageChops
from numba import cuda

class World:
    cells = {}
    original = Image.open("bird_gray.png").convert("L")

class ImageGenerator(EvolutionaryAlgorithm):

    def on_epoch(self, epoch):
        if epoch == 0:
            for population in self.populations:
                for cell in population.population:
                    self.generate_image[32, 32](cell)

        if epoch % 20 == 0:
            for p in self.populations:
                pop = p.population
                best = self.selection_from_pop(1, pop)[0]
                best.meta["image"].save(f'generated_images/epoch_{epoch}_pop_{p.number}.jpg', quality=95)
                with open(f"generators/epoch_{epoch}.pkl", "wb") as f:
                    dump(best, f)

    def on_child(self, child):
        self.generate_image[32, 32](child)

    @cuda.jit
    def generate_image(self, cell):
        dna = cell.dna
        im_size = (64, 64)
        im = Image.new('L', (im_size[0],im_size[1]), 255)
        draw = ImageDraw.Draw(im)
        for gene in dna:
            x = round(gene[0] * im_size[0])
            y = round(gene[1] * im_size[1])
            r = gene[2] * 3.5
            gray = round(gene[3] * 255)
            draw.ellipse((x-r, y-r, x+r, y+r), fill=(gray), outline=(gray))
        cell.meta["image"] = im

    @cuda.jit
    def fitness(self, cell):
        if cell.cached_fitness is not None:
            return cell.cached_fitness
        diff = ImageChops.difference(cell.meta["image"], World.original)
        diff = [d/255 for d in list(diff.getdata())]
        diff = sum(diff)
        fit = 1 / diff
        cell.cached_fitness = fit
        return fit

    def crossover(self, dna1, dna2):

        # --- Crossover
        geneA = int(random.random() * len(dna1))
        geneB = int(random.random() * len(dna1))

        startGene = min(geneA, geneB)
        endGene = max(geneA, geneB)

        childP1 = []
        for i in range(startGene, endGene):
            childP1.append(dna1[i])

        childP2 = [item for item in dna2 if item not in childP1]
        child_dna = childP2[:startGene] + childP1 + childP2[startGene:]
        return child_dna

if __name__ == "__main__":
    n_cells_per_pop = 35
    n_pops = 50
    n_genes = 350
    n_chroms = 4 # [posX, posY, radius, graytone]
    n_epochs = 1000
    verbose = True
    n_mutations = 88
    mutation_chance = 0.35

    evolutionary_algorithm = ImageGenerator(n_cells_per_pop, n_pops, n_genes, n_chroms, n_epochs, verbose, n_mutations, mutation_chance)
    evolutionary_algorithm.run()
    evolutionary_algorithm.plot_stats()
