// A gene is a short length of DNA found on a chromosome that codes for a particular characteristic or protein.
// Alleles are different forms of the same gene. For example, eye colour is the gene but blue, green, brown etc. are alleles.
// (BBC Bitesize)

package main

import (
	"cmp"
	"fmt"
	"math"
	"math/rand/v2"
	"slices"
)

type ProblemStatement struct {
	Genes            map[string]int
	ElitePopulation  int
	NormalPopulation int
	MutantPopulation int
	Generations      int
	Seed1            uint64
	Seed2            uint64
	MutationRate     float64
	Objective        func(*Chromosome) float64
}

type Definition struct {
	Genes        []string
	Max          []int
	Entropy      []float64 // measured in bits
	Alleles      []*Alleles
	TotalEntropy float64
	MutationRate float64
	Generator    *rand.Rand
	Objective    func(*Chromosome) float64
}

type Chromosome struct {
	Definition *Definition
	Values     []int
	Fitness    float64
	Origin     string
}

type Alleles struct {
	counter   int
	max       int
	alleles   []int
	generator *rand.Rand
}

func (a *Alleles) SampleWithoutReplacement() (allele int) {
	if a.counter == 0 {
		a.alleles = a.generator.Perm(a.max)
	}
	allele = a.alleles[a.counter]
	a.counter = (a.counter + 1) % a.max
	return allele // int
}

func Build_Definition(p ProblemStatement) (d Definition) {
	d.Genes = make([]string, len(p.Genes))
	d.Max = make([]int, len(p.Genes))
	d.Entropy = make([]float64, len(p.Genes))
	d.Alleles = make([]*Alleles, len(p.Genes))
	d.Generator = rand.New(rand.NewPCG(p.Seed1, p.Seed2))
	d.Objective = p.Objective
	var i int = 0
	for name, max := range p.Genes {
		d.Genes[i] = name
		d.Max[i] = max
		permutations := float64(max + 1)
		fully_committed_bits := math.Floor(math.Log2(permutations))
		permutations_in_fully_committed_bits := math.Pow(2, fully_committed_bits)
		permutations_in_remaining_bit := permutations - permutations_in_fully_committed_bits
		d.Entropy[i] = fully_committed_bits + (permutations_in_remaining_bit / permutations_in_fully_committed_bits) // entropy is measured in bits, e.g. integer of max size 256 has 8 bits of entropy
		d.Alleles[i] = &Alleles{counter: 0, max: max, generator: d.Generator}
		i++
	}
	return d // Defintion
}

func New_Chromosome(d Definition) (c *Chromosome) {
	c = &Chromosome{}
	c.Definition = &d
	c.Values = make([]int, len(d.Genes))
	for i, _ := range c.Values {
		c.Values[i] = d.Alleles[i].SampleWithoutReplacement()
	}
	c.Origin = "Randomised"
	return c // *Chromosome
}

func Crossover(a *Chromosome, b *Chromosome) (child *Chromosome) {
	child = &Chromosome{}
	child.Definition = a.Definition
	child.Values = make([]int, len(a.Values))
	for i, _ := range child.Values {
		if child.Definition.Generator.IntN(2) == 1 {
			child.Values[i] = a.Values[i]
		} else {
			child.Values[i] = b.Values[i]
		}
	}
	child.Origin = "Crossover"
	return child // *Chromosome
}

func Mutate(c *Chromosome) {
	var i int = 0
	var j float64 = 0
	stop := c.Definition.TotalEntropy * c.Definition.Generator.Float64()
	for j < stop {
		j += c.Definition.Entropy[i] // more complex genes recieve more frequent mutations
		i++
	}
	c.Values[i] = c.Definition.Alleles[i].SampleWithoutReplacement()
	c.Origin = "Mutation"
}

func Evaluate(c *Chromosome) {
	c.Fitness = c.Definition.Objective(c)
}

func Rank(population []*Chromosome) {
	comparison := func(a, b *Chromosome) int {
		return -cmp.Compare(a.Fitness, b.Fitness) // sort descending
	}
	slices.SortStableFunc(population, comparison)
}

func SelectParents(population []*Chromosome) (a, b *Chromosome) {
	// assume population has already been ranked
	pool := (len(population) / 2) + len(population)%2 // adding modulo 2 rounds up the integer division
	position_a := population[0].Definition.Generator.IntN(pool)
	position_b := population[0].Definition.Generator.IntN(pool - 1)
	if position_a == position_b {
		position_b += 1 // b is sampled from a population which does not include a
	}
	a = population[position_a]
	b = population[position_b]
	return a, b // *Chromosome
}

func SprintGenerations(generations [][]*Chromosome) (output string) {
	for i, _ := range generations { // generation number
		output += fmt.Sprintf("Generation: %d\n", i)
		if generations[i][0] == nil {
			output += "(not calculated)"
		} else {
			for j, _ := range generations[i] { // chromosome number
				output = fmt.Sprint(output, generations[i][j].Fitness, ` `, generations[i][j].Origin, ` `, generations[i][j].Values, `, `)
			}
		}
		output += "\n"
	}
	return output // string
}

func (p ProblemStatement) Solve() (solution *Chromosome, chromosomes [][]*Chromosome) {
	chromosome_definition := Build_Definition(p)
	population := p.ElitePopulation + p.NormalPopulation + p.MutantPopulation
	chromosomes = make([][]*Chromosome, p.Generations)
	all_generations := make([]*Chromosome, p.Generations*population)
	for i, _ := range chromosomes {
		chromosomes[i] = all_generations[i*population : (i+1)*population : (i+1)*population]
	}
	var parent_a *Chromosome // for crossover
	var parent_b *Chromosome // for crossover
	var i int                // generation number
	var j int                // chromosome number
	for j = 0; j < population; j++ {
		chromosomes[0][j] = New_Chromosome(chromosome_definition)
		Evaluate(chromosomes[0][j])
	}

	Rank(chromosomes[0])
	for i = 1; i < p.Generations; i++ {
		j = 0
		for j < p.ElitePopulation {
			copy := *chromosomes[i-1][j] // copy the underlying chromosome
			chromosomes[i][j] = &copy
			chromosomes[i][j].Origin = "Elite"
			j++
		}
		parents := population / 2
		remix := chromosome_definition.Generator.Perm(parents)
		for n := 0; n < p.NormalPopulation; n++ {
			parent_a = chromosomes[i-1][n%parents]
			parent_b = chromosomes[i-1][remix[n%parents]]
			chromosomes[i][j] = Crossover(parent_a, parent_b)
			Evaluate(chromosomes[i][j])
			j++
		}
		remix = chromosome_definition.Generator.Perm(parents)
		for n := 0; n < p.MutantPopulation; n++ {
			parent_a = chromosomes[i-1][n%parents]
			parent_b = chromosomes[i-1][remix[n%parents]]
			chromosomes[i][j] = Crossover(parent_a, parent_b)
			Mutate(chromosomes[i][j])
			Evaluate(chromosomes[i][j])
			j++
		} // j = population, new generation is fully spawned
		Rank(chromosomes[i])
		if i > 2 && chromosomes[i][0].Fitness <= chromosomes[i-1][0].Fitness && chromosomes[i-1][0].Fitness <= chromosomes[i-2][0].Fitness {
			break // halt early if no improvement for two generations
		}
	}
	solution = chromosomes[i-1][0]
	return solution, chromosomes // *Chromosome, [][]*Chromosome
}

func total(c *Chromosome) float64 {
	return float64(c.Values[0] + c.Values[1])
}

func main() {
	genes := map[string]int{`a`: 100, `b`: 27}
	problem := ProblemStatement{
		Genes:            genes,
		ElitePopulation:  1,
		NormalPopulation: 10,
		MutantPopulation: 10,
		Generations:      10,
		Seed1:            0,
		Seed2:            0,
		MutationRate:     0,
		Objective:        total,
	}
	solution, workings := problem.Solve()
	fmt.Println(`Solution`)
	fmt.Println(solution)
	fmt.Println(SprintGenerations(workings))
}
