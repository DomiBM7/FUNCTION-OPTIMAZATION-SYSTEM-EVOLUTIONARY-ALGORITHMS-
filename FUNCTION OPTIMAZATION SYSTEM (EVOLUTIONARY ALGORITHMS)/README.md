# Differential Evolution (DE) Optimization Project

## Overview
This project implements various Differential Evolution (DE) algorithms for optimization problems. It includes 15 different DE variants and provides comprehensive performance analysis tools.

## Features
- 15 Different DE Variants including:
  - DE/rand/1/bin
  - DE/rand/1/exp
  - DE/current/1/bin
  - DE/best/2/bin
  - DE/proportional/2/exp
  - DE/tournament/1/exp
  - DE/rand-to-best/4/bin
  - DE/current-to-best/2/bin
  - DE/best/1/writeLinearOperator
  - DE/best/1/writeHeuristicOperator
  - DE/best/1/MichalewiczArithmetic
  - DE/best/1/MichalewiczGeometric
  - DE/best/1/headlessChicken
  - DE/best/1/headlessChicken/AdaptiveMutation
  - DE/rank/1/MichalewiczArithmetic

## Performance Metrics
The implementation tracks various performance metrics including:
- Best fitness values
- Average population fitness
- Worst fitness values
- Population range
- Mutation rates
- Selective pressure
- Success rates
- Convergence metrics
- Population deviation

## Visualization
The project includes plotting capabilities through the `MultiArrayPlotter` class to visualize:
- Fitness evolution
- Population statistics
- Algorithm convergence

## Usage
1. Run the main program
2. Select one of the 15 DE variants from the menu
3. The algorithm will run for the specified number of generations
4. Results will be displayed including:
   - Best solution found
   - Worst solution found
   - Average fitness
   - Various performance metrics
   - Visualization plots

## Parameters
- Population Size: Defined by POPULATION_SIZE constant
- Dimensions: Defined by DIMENSIONS constant
- Maximum Generations: Defined by MAX_GENERATIONS constant
- Mutation Factor (M): Self-adaptive in some variants
- Crossover Rate (CR): Configurable parameter

## Requirements
- Java 8 or higher
- JFreeChart library for plotting
- See requirements.txt for detailed dependencies

## Performance Analysis
The implementation provides comprehensive performance analysis including:
- Convergence analysis
- Statistical measures
- Fitness tracking
- Population diversity metrics