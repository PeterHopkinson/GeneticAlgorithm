# GeneticAlgorithm

Optimisation problem solver written in Golang.

## How To Use - Main

Compile the script and run the binary.<br>
The script has a main method which will create and run a demonstration problem.

## How To Use - Library

If the main method is removed the script can be used as a library.<br>
If using it as a library you will need to import the ProblemStatement and Chromosome classes and do the following:
* write an objective function with the pattern: `func objective (c Chromosome) float64 {}`
* create a ProblemStatement{} literal
* run the solver: `solution, workings = ProblemStatement.Solve()`
    * solution is a pointer to Chromosome object representing the best solution the algorithm found
    * workings is a pointer to an array of chromosomes representing all the generations which were tested
    * `SprintGenerations(workings)` is a convenient way to format the array as a string

## To Do

This is the first program I've written in Go.<br>
There are several things that need tidying up:
* reduce the number of Public functions and variables
* move the main method to another file, to better delineate the library files from the demo problem
* work out a better way of inputting genes (the current approach using a map puts the genes into a random order, which isn't great!)
* improve the interface so fewer objects need to be imported
    * aim is to reduce the imports to ProblemStatement plus a set of summary and print functions
* replace instances of snake_case with camelCase
* add binary for non-Go users