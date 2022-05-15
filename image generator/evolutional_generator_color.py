from evolutionary_algorithm_v2 import EvolutionaryAlgorithm, Cell
import random
from pickle import dump
from PIL import Image, ImageDraw, ImageChops


class World:
    cells = {}
    original = Image.open("original_color.png").convert("RGB")

class ImageGenerator(EvolutionaryAlgorithm):

    def on_epoch(self, epoch):
        if epoch == 0:
            for population in self.populations:
                for cell in population.population:
                   self.generate_image(cell)

        if epoch % 20 == 0:
            best = self.selection(n_selected=1)[0]
            #with open(f"generators/epoch_{epoch}.pkl", "wb") as f:
            #    dump(best, f)
            best.meta["image"].save(f'generated_images/color/epoch_{epoch}.jpg', quality=95)

    def on_child(self, child):
        self.generate_image(child)

    def generate_image(self, cell):
        dna = cell.dna
        im_size = (80, 80)
        im = Image.new('RGB', (im_size[0], im_size[1]), (255, 255, 255))
        draw = ImageDraw.Draw(im)
        for gene in dna:
            x = round(gene[0] * im_size[0])
            y = round(gene[1] * im_size[1])
            radius = round(gene[2] * 5)
            r = round(gene[3] * 255)
            g = round(gene[4] * 255)
            b = round(gene[5] * 255)
            draw.ellipse((x-radius, y-radius, x+radius, y+radius), fill=(r, g, b), outline=(r, g, b))
        cell.meta["image"] = im

    def fitness(self, cell):
        if cell.cached_fitness is not None:
            return cell.cached_fitness
        self.generate_image(cell)
        diff = ImageChops.difference(cell.meta["image"], World.original)
        diff = [sum([e / 255 for e in d]) for d in list(diff.getdata())]
        diff = sum(diff)
        fit = 1 / diff
        cell.cached_fitness = fit
        return fit

if __name__ == "__main__":
    n_cells_per_pop = 50
    n_pops = 40
    n_genes = 300
    n_chroms = 6 # [posX, posY, radius, r, g, b]
    n_epochs = 1000
    verbose = True
    n_mutations = 75
    mutation_chance = 0.2

    evolutionary_algorithm = ImageGenerator(n_cells_per_pop, n_pops, n_genes, n_chroms, n_epochs, verbose, n_mutations, mutation_chance)
    evolutionary_algorithm.run()
    evolutionary_algorithm.plot_stats()
