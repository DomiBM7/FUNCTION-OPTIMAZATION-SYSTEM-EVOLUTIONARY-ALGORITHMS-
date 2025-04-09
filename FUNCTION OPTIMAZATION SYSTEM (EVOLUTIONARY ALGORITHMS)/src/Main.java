import java.util.Random;
import java.util.Scanner;

public class h {
// ================================================================================================== PROBLEM-SPECIFIC PARAMETERS ==================================================================================================
    static final int DIMENSIONS = 1; // Number of dimensions for the optimization problem
    static final int POPULATION_SIZE = 20; // Population size
    static final int MAX_GENERATIONS = 2000; // Maximum number of generations
    static  double M = 0.3; // Mutation factor
    static  double CR = 1; // Crossover rate
    static final double B = 10;
    static final double maxDouble = Double.MAX_VALUE;
    
    private static final double INITIAL_TEMPERATURE = 100.0;
    private static final double COOLING_RATE = 0.995;
    private static final double MAX_DOUBLE = Double.MAX_VALUE;
    private static final Random rand = new Random();
    private static int generation = 0;
    static double ss0 = 500;
    static boolean blnContinue = true;
    //int index = 0;
    static Scanner scanner = new Scanner(System.in);
    
        
    public static double temperature() {
        return INITIAL_TEMPERATURE * Math.pow(COOLING_RATE, generation++);
    }

    public static void main(String[] args) {
        double[][] population = initializePopulation(); // Initialize population
        double[][] offsprings = new double [POPULATION_SIZE][DIMENSIONS];
        
        

//================================================================================================== KEEPING TRACK OF THE INITIAL POPULATION FOR COMPARISON ==================================================================================================
        double [][] initialPopulation = new double[POPULATION_SIZE][DIMENSIONS];
        for (int i = 0; i < POPULATION_SIZE; i++) {
                for (int j = 0; j < DIMENSIONS; j++) {
                	initialPopulation[i][j] = population[i][j];
                }
         }
//================================================================================================== HANDLE USER INPUT ==================================================================================================
        double[] t = new double[DIMENSIONS];
        double sr = 0.0;
        double[] ta = new double[DIMENSIONS];
        String[] strChoices = { "1. DE/rand/1/bin",
                "2. DE/rand/1/exp",
                "3. DE/current/1/bin",
                "4. DE/best/2/bin",
                "5. DE/propotional/2/exp",
                "6. DE/tournament/1/exp",
                "7. DE/rand-to-best/4/bin",
                "8. DE/current-to-best/2/bin",
                "9. DE/best/1/writeLinearOperator",
                "10. DE/best/1/writeHeuristicOperator",
                "11. DE/best/1/MichalewiczArithmetic",
                "12. DE/best/1/MichalewiczGeometric",
                "13. DE/best/1/headlessChicken",
                "14. DE/best/1/headlessChicken/AdaptiveMutation",
                "15. DE/rank/1/MichalewiczArithmetic"
                };
        System.out.println("Chose the DE pipeline you want to implement");
        
        for(int h =0 ; h < strChoices.length; h++)
        {
        	System.out.println(strChoices[h]);
        }
        System.out.println("");
        
        int choice = scanner.nextInt();
        double[] bestIndividuals = new double[MAX_GENERATIONS];
        double[] averageIndividuals = new double[MAX_GENERATIONS];
        double[] worstIndividuals = new double[MAX_GENERATIONS];
        double[] rangeIndividuals = new double[MAX_GENERATIONS];
        double[] mutationIndividuals = new double[MAX_GENERATIONS];
        double[] selectivePressureIndividuals = new double[MAX_GENERATIONS];
        double[] successRateIndividuals = new double[MAX_GENERATIONS];
        double[] convergenceIndividuals = new double[MAX_GENERATIONS];
        double[] bestFitnessIndividuals = new double[MAX_GENERATIONS];
        double[] maxFitness = new double[MAX_GENERATIONS];
        double[] maxValues = new double[MAX_GENERATIONS];
        double[] deviationIndividuals = new double[MAX_GENERATIONS];
        
//================================================================================================== PERFORM MUTATION AND CROSSOVER ==================================================================================================
        if(stopCondition_AS(population) && stopCondition_OFD(population))
        {
        for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        	//System.out.println("Gen:" + generation);
            for (int i = 0; i < POPULATION_SIZE; i++) {
            	
            	double[] targetVector = new double[DIMENSIONS];
            	double[] finalDifferenceVector = new double[DIMENSIONS];
            	double[] finalTrial = new double[DIMENSIONS];
            	double[] finalOffspring = new double[DIMENSIONS];
            	//M = selfAdaptionLogNormal(generation);
            	switch(choice)
            	{
	            	case 1:
	            		// Mutation
	                    int a, b, c;
	                    do { a = rand.nextInt(POPULATION_SIZE); } while (a == i);
	                    do { b = rand.nextInt(POPULATION_SIZE); } while (b == i || b == a);
	                    do { c = rand.nextInt(POPULATION_SIZE); } while (c == i || c == a || c == b);
	                    
	                    targetVector = population[a];
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + M * (population[b][j] - population[c][j]);
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = binaryCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    
	                    t = finalOffspring;
	                    ta = targetVector;
	            		break;
	            	case 2:
	            		// Mutation
	                    int d, e, f;
	                    do { d = rand.nextInt(POPULATION_SIZE); } while (d == i);
	                    do { e = rand.nextInt(POPULATION_SIZE); } while (e == i || e == d);
	                    do { f = rand.nextInt(POPULATION_SIZE); } while (f == i || f == d || f == d);
	                    
	                    targetVector = population[d];
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + M * (population[e][j] - population[f][j]);
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = exponentialCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    
	                    t = finalOffspring;
	                    ta = targetVector;
	            		break;
	            	case 3:
	            		
	            		double[] d1 = randomSelection(population);
	            		double[] d2 = randomSelection(population);
	                    finalDifferenceVector[0] = M*(d1[0]-d2[0]);
	                    targetVector = population[i];
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = binaryCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    
	                    t = finalOffspring;
	                    ta = targetVector;
	            		break;
	            	case 4:
	            		double[] d1_ = randomSelection(population);
	            		double[] d2_ = randomSelection(population);
	            		double[] d3_ = randomSelection(population);
	            		double[] d4_ = randomSelection(population);
	            		
	                    finalDifferenceVector[0] = M*((d1_[0]-d2_[0])+(d3_[0]-d4_[0]));
	                    targetVector = elitismSelection(population);
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = binaryCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    
	                    t = finalOffspring;
	                    ta = targetVector;
	            		break;
	            	case 5:
	            		double[] d1_5 = randomSelection(population);
	            		double[] d2_5 = randomSelection(population);
	            		double[] d3_5 = randomSelection(population);
	            		double[] d4_5 = randomSelection(population);
	            		
	                    finalDifferenceVector[0] = M*((d1_5[0]-d2_5[0])+(d3_5[0]-d4_5[0]));
	                    targetVector = propotionalSelection(population);
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = exponentialCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    
	                    t = finalOffspring;
	                    ta = targetVector;
	            		break;
	            	
	            	case 6:
	            		double[] d1_6 = tournamentSelection(population);
	            		double[] d2_6 = tournamentSelection(population);
	            		
	            		
	                    finalDifferenceVector[0] = M*((d1_6[0]-d2_6[0]));
	                    targetVector = tournamentSelection(population);
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = exponentialCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            	case 7:
	            		double[] d1_7 = tournamentSelection(population);
	            		double[] d2_7 = tournamentSelection(population);
	            		double[] d3_7 = tournamentSelection(population);
	            		double[] d4_7 = tournamentSelection(population);
	            		double[] d5_7 = tournamentSelection(population);
	            		double[] d6_7 = tournamentSelection(population);
	            		double[] d7_7 = tournamentSelection(population);
	            		double[] d8_7 = tournamentSelection(population);
	            		
	            		
	                    finalDifferenceVector[0] = M*((d1_7[0]-d2_7[0])+ (d3_7[0]-d4_7[0])+ (d5_7[0]-d6_7[0])+ (d7_7[0]-d8_7[0]));
	                    double[] tv1 = randomSelection(population);
	                    double[] tv2 = elitismSelection(population);
	                    targetVector[0] = 0.5*tv1[0] + tv2[0];
	                    
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = binaryCrossover(population[i], finalTrial, CR);
	                  
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            	case 8:
	            		double[] d1_8 = tournamentSelection(population);
	            		double[] d2_8 = tournamentSelection(population);
	            		double[] d3_8 = tournamentSelection(population);
	            		double[] d4_8 = tournamentSelection(population);
	            		
	            		
	                    finalDifferenceVector[0] = M*((d1_8[0]-d2_8[0])+(d3_8[0]-d4_8[0]));
	                    double[] tv1_8 = population[i];
	                    double[] tv2_8 = elitismSelection(population);
	                    targetVector[0] = tv1_8[0] + tv2_8[0];
	                    
	                    
	                    for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = binaryCrossover(population[i], finalTrial, CR);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            		
	            	case 9:
	            		// Mutation
	                    
	                    targetVector = elitismSelection(population);
	                    
	                    double[] d1_9 = tournamentSelection(population);
	            		double[] d2_9 = tournamentSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_9[0]-d2_9[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);
	                    finalOffspring = writeLinearOperator(population);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            		
	            	case 10:
	            		// Mutation
	                    
	                    targetVector = elitismSelection(population);
	                    
	                    double[] d1_10 = tournamentSelection(population);
	            		double[] d2_10 = tournamentSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_10[0]-d2_10[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = writeHeuristicOperator(population);
	                  
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;

	            		break;
	            	case 11:
	            		
	            		targetVector = elitismSelection(population);
	                    
	                    double[] d1_11 = tournamentSelection(population);
	            		double[] d2_11 = tournamentSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_11[0]-d2_11[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = MichalewiczArithmetic(population);
	                   
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            		
	            	case 12:
	            		targetVector = elitismSelection(population);
	                    
	                    double[] d1_12 = elitismSelection(population);
	            		double[] d2_12 = elitismSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_12[0]+d2_12[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = MichalewiczArithmetic(population);
	                    
	                    //replaceWorst(population, trial, i);
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            		
	            	case 13:
	            		targetVector = elitismSelection(population);
	                    
	                    double[] d1_13 = tournamentSelection(population);
	            		double[] d2_13 = tournamentSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_13[0]-d2_13[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = headlessChicken(population, finalOffspring);
	                   
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            	case 14:
	            		targetVector = elitismSelection(population);
	                    double m = selfAdaptionLogNormal(generation);
	                    double[] d1_14 = tournamentSelection(population);
	            		double[] d2_14 = tournamentSelection(population);
	            		M = m;
	            		finalDifferenceVector[0] = M*((d1_14[0]-d2_14[0]));
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = headlessChicken(population, finalOffspring);
	                    
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            		
	            	case 15:
	            		targetVector = rankBasedtSelection(population);
	                    
	                    double[] d1_15 = randomSelection(population);
	            		double[] d2_15 = elitismSelection(population);
	            		
	            		finalDifferenceVector[0] = M*((d1_15[0]+d2_15[0]));
	            		
	            		if (generation % 100 == 0) {
	                            
	                                if (rand.nextDouble() < 0.15) { // 10% chance to reinitialize
	                                    for (int j = 0; j < DIMENSIONS; j++) {
	                                        population[i][j] = rand.nextDouble() * (2) - 1;
	                                    }
	                                }
	                            
	                           
	                        
	                	}
	                    
	                    System.out.println(selfAdaptionLogNormal(generation));
	                    
	                    
	            		
	            		for (int j = 0; j < DIMENSIONS; j++) {
	                    	finalTrial[j] = targetVector[j] + finalDifferenceVector[j];
	                    }
	            		
	                   
	                    finalTrial = repair(finalTrial);

	                    // Crossover
	                    //double[] trial = new double[DIMENSIONS];
	                    
	                    finalOffspring = MichalewiczGeometric(population);
	                    
	                    //replaceWorst(population, trial, i);
	                    offsprings[i] = finalOffspring;
	                    parentChild(population, finalOffspring, i);
	                    ta = targetVector;
	                    t = finalOffspring;
	            		break;
	            	default:
	            		System.out.println("Please chose a valid option");
	            		break;
            	
            	}
            	
            }
            
            for(int p = 0; p < POPULATION_SIZE; p++)
            {
            	averageIndividuals[generation] += population[p][0]; 
            }
            
            averageIndividuals[generation] = averageIndividuals[generation]/POPULATION_SIZE; 
            
            

            double[] best = findBestSolution(population);
            //System.out.println("Best Solution Found at: "+ best[0]);
            double[] worst = findWorstSolution(population);
            //System.out.println("Worst Solution Found at: "+ worst[0]);
            bestIndividuals[generation] = (best[0]);
            worstIndividuals[generation] = worst[0];
            
            
            
            
            
            for(int p = 0; p < POPULATION_SIZE; p++)
            {
            	if(fitness(t) > fitness(population[p]))
            	{
            		successRateIndividuals[generation] += 1.0;
            		sr += successRateIndividuals[generation];
            	}
            	else if(fitness(t) < fitness(population[p]))
            	{
            		successRateIndividuals[generation] += -1.0;
            		sr += successRateIndividuals[generation];
            	}
            	
            	successRateIndividuals[generation] = 0.0;
            }
            
            for(int p = 0; p < MAX_GENERATIONS; p++)
            {
            	convergenceIndividuals[p] = fitness(best) - fitness(worst); 
            }
            
            for(int p = 0; p < MAX_GENERATIONS; p++)
            {
            	bestFitnessIndividuals[p] = fitness(best); 
            }
            
            for(int p = 0; p < MAX_GENERATIONS; p++)
            {
            	maxFitness[p] = 3.962; 
            }
            for(int p = 0; p < MAX_GENERATIONS; p++)
            {
            	maxValues[p] = -1.088598; 
            }
            for(int p = 0; p < MAX_GENERATIONS; p++)
            {
            	double sum_diff = 0;
            	double sum = 0;
            	double sigma = 0;
            	for (int i = 0; i < population.length; i++) {
            		sum += population[i][0];
            	}
            	
            	double mean = sum/population.length;
            	
            	for (int i = 0; i < population.length; i++) {
            		sum_diff += Math.pow(population[i][0]-mean, 2);
            	}
            	
            	sigma = Math.sqrt(sum_diff/population.length);
            	deviationIndividuals[p] = sigma; 
            }
            
          
       
        }
        
        for(int p = 0; p < MAX_GENERATIONS; p++)
        {
        	rangeIndividuals[p] =  ( Math.abs(bestIndividuals[p]) -Math.abs(worstIndividuals[p])); 
        }
        
        for(int p = 0; p < MAX_GENERATIONS; p++)
        {
        	mutationIndividuals[p] =M; 
        }
        
        for(int p = 0; p < MAX_GENERATIONS; p++)
        {
        	selectivePressureIndividuals[p] = fitness(ta); 
        }
        
        double[] in = new double[POPULATION_SIZE];
        for(int ab = 0; ab < POPULATION_SIZE; ab++)
        {
        	in[ab] = population[ab][0];
        }
        
        printPopulation(initialPopulation, population);
       /* MultiArrayPlotter fitnessplot = new MultiArrayPlotter(strChoices[choice-1], bestIndividuals, 
        		rangeIndividuals,
        		mutationIndividuals);
        fitnessplot.plotArrays(strChoices[choice-1],worstIndividuals, 
        		averageIndividuals, 
        		bestIndividuals, 
        		rangeIndividuals,
        		mutationIndividuals,
        			selectivePressureIndividuals,
        			successRateIndividuals,
        			convergenceIndividuals,
        			bestFitnessIndividuals,
        			maxFitness,
        			maxValues,
        		deviationIndividuals);*/
       
        double[] best = findBestSolution(population);
        System.out.println("Best Solution Found: "+ best[0]);
        double[] worst = findWorstSolution(population);
        System.out.println("Worst Solution Found: "+ worst[0]);
        System.out.println("Average Individual: "+averageIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Range of Individuals: "+rangeIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Selective pressure: "+selectivePressureIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Success Rate Individuals: "+sr/MAX_GENERATIONS);
        System.out.println("Deviation of Individuals: "+deviationIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Best Fitness Individuals: "+bestFitnessIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Convergence Individuals: "+convergenceIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Mutation Individuals: "+mutationIndividuals[MAX_GENERATIONS -1]);
        System.out.println("Absolute Difference Individuals: "+ Math.abs(Math.abs(averageIndividuals[MAX_GENERATIONS -1]) - Math.abs(-1.088598)));
        double[] f = {averageIndividuals[MAX_GENERATIONS -1]};
        System.out.println("Absolute Difference Fitness: "+ Math.abs(Math.abs(fitness(f)) - Math.abs(3.962)));
        
    }
        ;
}
            

//================================================================================================== PRINT THE SOLUTION ==================================================================================================
    static void printBestSolution(double[][] population)
    {
    	// Output the best solution
        double[] best = findBestSolution(population);
        int k = findBestSolutionIndex(population);
        System.out.println("Best solution found: " + java.util.Arrays.toString(best)+ "at "+ k);
    }
    
    
//================================================================================================== PRINT THE POPULATION ==================================================================================================
    static void printPopulation(double[][] population, double[][] population2)
    {
    	System.out.println("===================================================================== INITIAL POPULATION");
        
        for(int l = 0; l < POPULATION_SIZE; l++)
        {
        	for (int j = 0; j < DIMENSIONS; j++) {
        		System.out.print( population[l][j]+" ");
            }
        	System.out.println(" ");
        }
        
        
        System.out.println("  ");
        System.out.println("===================================================================== FINAL POPULATION");
        
        for(int l = 0; l < POPULATION_SIZE; l++)
        {
        	for (int j = 0; j < DIMENSIONS; j++) {
        		System.out.print( population2[l][j]+" ");
            }
        	System.out.println(" ");
        }
    }

//================================================================================================== PRINT THE POPULATION FITNESS ==================================================================================================
    static void printPopulationFitness(double[][] population, double[][] population2)
    {
    	System.out.println("===================================================================== INITIAL POPULATION FITNESS");
        
        for(int l = 0; l < POPULATION_SIZE; l++)
        {
        	for (int j = 0; j < DIMENSIONS; j++) {
        		System.out.print( fitness(population[l])+" ");
            }
        	System.out.println(" ");
        }
        
        
        System.out.println("  ");
        System.out.println("===================================================================== FINAL POPULATION FITNESS");
        
        for(int l = 0; l < POPULATION_SIZE; l++)
        {
        	for (int j = 0; j < DIMENSIONS; j++) {
        		System.out.print( fitness(population2[l])+" ");
            }
        	System.out.println(" ");
        }
    }

//================================================================================================== INITIALISE THE POPULATION WITH RANDOM VALUES ==================================================================================================
    static double[][] initializePopulation() 
{
        double[][] population = new double[POPULATION_SIZE][DIMENSIONS];
        Random rand = new Random();
        double minRange = -50;
        double maxRange = 50;
        for (int i = 0; i < POPULATION_SIZE; i++) {
            for (int j = 0; j < DIMENSIONS; j++) {
                population[i][j] = rand.nextDouble(-10000,10000);
            }
        }
        return population;
    }
//================================================================================================== CALCULATING FITNESS ==================================================================================================
    static double fitness(double[] solution) 
    {
    	double sum = 0;
        for (int i = 0; i < solution.length; i++) {
        	
            sum = -solution[i]*solution[i]*solution[i]*solution[i]*solution[i]*solution[i] + 4*solution[i]*solution[i] + Math.sin(-solution[i]) ;//+ Math.log(solution[i]); // Fitness function
        }
        //System.out.println("Fitness: "+sum);
        return sum;
    }
    
    //METHOD TO enforce range
    static double[] repair(double[] solution) {
        for (int i = 0; i < solution.length; i++) {
          //# if (solution[i] < -10) solution[i] = 0;
          //  if (solution[i] > 10) solution[i] = 0;
        }
        return solution;
    }

    // method to print details of DE reproduction
    static void logIndividualDetails(double[] parent, double[] trial, double[] child)
    {
    	System.out.print(" Parent: [ ");
        for (double v : parent) System.out.print(v + ", ");
        System.out.println(" ] ");

        System.out.print(" Trial: [ ");
        for (double v : trial) System.out.print(v + ", ");
        System.out.println(" ] ");

        System.out.print(" Child: [ ");
        for (double v : child) System.out.print(v + ", ");
        System.out.println(" ] "); 
    }

    // Find the best solution in the population
    static double[] findBestSolution(double[][] population) {
        double[] best = population[0];
        double bestFitness = fitness(best);
        for (int i = 1; i < population.length; i++) {
            double currentFitness = fitness(population[i]);
            if (currentFitness > bestFitness) { // Maximizing
                best = population[i];
                bestFitness = currentFitness;
                
            }
        }
        return best;
    }
    
 // Find the index of the best solution in the population
    static int findBestSolutionIndex(double[][] population) {
        double[] best = population[0];
        double bestFitness = fitness(best);

        int p = 0; 
        for (int i = 1; i < population.length; i++) {
        	System.out.println(fitness(population[i]));
            double currentFitness = fitness(population[i]);
            if (currentFitness > bestFitness) { // Maximizing
                best = population[i];
                bestFitness = currentFitness;
                p = i; 
            }
        }
        
        
        return p; // Return the index of the best solution after the loop
    }

    //find the worst solution in the population
    static double[] findWorstSolution(double[][] population) {
        double[] best = population[0];
        double bestFitness = fitness(best);
        for (int i = 1; i < population.length; i++) {
            double currentFitness = fitness(population[i]);
            if (currentFitness < bestFitness) { // Minimizing
                best = population[i];
                bestFitness = currentFitness;
            }
        }
        return best;
    }
    
 // Find the index of the worst solution in the population
    static int findWorstSolutionIndex(double[][] population) {
        double[] best = population[0];
        double bestFitness = fitness(best);

        int p = 0; 
        for (int i = 1; i < population.length; i++) {
        	System.out.println(fitness(population[i]));
            double currentFitness = fitness(population[i]);
            if (currentFitness < bestFitness) { // Maximizing
                best = population[i];
                bestFitness = currentFitness;
                p = i; 
            }
        }
        return p; // Return the index of the worst solution after the loop
    }
    
//================================================================================================== SELECTION ==================================================================================================
//-Random
    static double[] randomSelection(double[][] population) {
    	Random rand = new Random();
    	int j = rand.nextInt(POPULATION_SIZE); 
        double[] randIndividual = population[j];
        return randIndividual;
    }
    
//-Proportional Selection
    static double[] propotionalSelection(double[][] population) {
    	double totalFitness = 0;
    	
    	for (double[] v : population)
    	{
    		double scaledFitness = fitness(v) + maxDouble;
    		totalFitness += scaledFitness;
    	}
//-Roulette Wheel Selection
    	int i = 0;
    	double sum = (fitness(population[i])+maxDouble)/totalFitness;
    	double r = rand.nextDouble();
    	int index = 0;
    	
    	while(sum < r && i < population.length -1)
    	{
    		i = i +1;
    		double probSelection = (fitness(population[i])+maxDouble)/totalFitness;
    		sum += probSelection;
    	}
    	
    	return population[i];
    }
    
//-Tournament
    static double[] tournamentSelection(double[][] population) {
    	double tournamentSize = 5;
    	double[][] tournamentPopulation = new double[POPULATION_SIZE][DIMENSIONS];
    	Random rand = new Random();
    	for(int i = 0; i<tournamentSize; i++)
    	{
    		int j = rand.nextInt(POPULATION_SIZE); 
            double[] randIndividual = population[j];
            tournamentPopulation[i] = randIndividual;
    	}
    	
        return findBestSolution(tournamentPopulation);
    }
    
//-Rank-Based
    static double[] rankBasedtSelection(double[][] population) {
    	
    	double[][] rankPopulation = population;
    	Random rand = new Random();
    	for(int i = 0; i<POPULATION_SIZE - 1; i++)
    	{
    		if (fitness(population[i]) > fitness(population[i + 1]))
    		{
    			double[] temp = rankPopulation[i];
    			rankPopulation[i] = rankPopulation[i + 1];
    			rankPopulation[i + 1] = temp;
    		}
    	}
    	
    	int k = rand.nextInt(0, POPULATION_SIZE); // this works fine as long as POPULATION_SIZE > 0
    	int l;
    	if (k > 0) {
    	    l = rand.nextInt(0, k); // when k > 0
    	} else {
    	    l = 0; // Handle case when k = 0
    	}
    	
    	return rankPopulation[k];
    }
    
//-Boltzman
    static double[] BoltzmanlSelection(double[][] population) {
    	double temperature = temperature();
        double totalFitness = 0;
        
        // Calculate total Boltzmann-adjusted fitness
        for (double[] individual : population) {
            double boltzmannFitness = Math.exp(fitness(individual) / temperature);
            totalFitness += boltzmannFitness + MAX_DOUBLE;
        }

        // Roulette Wheel Selection
        int i = 0;
        double sum = (Math.exp(fitness(population[i]) / temperature) + MAX_DOUBLE) / totalFitness;
        double r = rand.nextDouble();

        while (sum < r && i < population.length - 1) {
            i = i + 1;
            double boltzmannFitness = Math.exp(fitness(population[i]) / temperature);
            double probSelection = (boltzmannFitness + MAX_DOUBLE) / totalFitness;
            sum += probSelection;
        }

        return population[i];
    }
//-Elitism
    static double[] elitismSelection(double[][] population) {
    	
    	double[] best = findBestSolution(population);
    	return best;
    }

   
// ================================================================================================== Reproduction ==================================================================================================
// Write’s linear operator approach.
    static double[] writeLinearOperator(double[][] population) {
    	int choice = rand.nextInt(0, 4);
    	double[] first = randomSelection(population);
    	double[] second = randomSelection(population);
    	double[] child1 = new double[first.length];
    	double[] child2 = new double[first.length];
    	double[] child3 = new double[first.length];
    	double[] child4 = new double[first.length];
    	double[] third = elitismSelection(population);
    	for(int i = 0; i < first.length; i++)
    	{
    		child1[i] = first[i] + second[i];
    		child2[i] = 1.5*first[i] -0.5* second[i];
    		child3[i] = -1.5*first[i] + 0.5*second[i];
    		// New innovation
    		child4[i] = third[i] + 0.5*first[i] + 0.5*second[i];
    	}
    	
    	if(choice == 0)
    		return child1;
    	if(choice == 1)
    		return child2;
    	if(choice == 2)
        		return child3;
    	return child4;
    }
// Write’s directional heuristic crossover:
    static double[] writeHeuristicOperator(double[][] population) {
    	Random rand = new Random();
    	double u = rand.nextDouble();
    	double[] first = randomSelection(population);
    	double[] second = randomSelection(population);
    	double[] child = new double[DIMENSIONS];
    	if(fitness(first) > fitness(second))
    	{
    		for(int i = 0; i < 1; i++)
        	{
    			//System.out.print(i);
        		child[i] = u * (first[i] - second[i]) + first[i];
        	}
    	}
    	else
    	{
    		for(int i = 0; i < first.length; i++)
        	{
        		child[i] = u * (second[i] - first[i]) + second[i];
        	}
    	}
    	return child;
    }
// Michalewicz’s Arithmetic Crossover
    static double[] MichalewiczArithmetic(double[][] population) {
    	// We crate a small tournament, but do a wighted average, insted of chosing the best
    	int tournamentSize = 5;
    	double[][] tournamentPopulation = new double[POPULATION_SIZE][DIMENSIONS];
    	Random rand = new Random();
    	double[] child = new double[1];
    	for(int i = 0; i<tournamentSize; i++)
    	{
    		int j = rand.nextInt(POPULATION_SIZE); 
            double[] randIndividual = population[j];
            tournamentPopulation[i] = randIndividual;
    	}
    	double[] array = new double[tournamentSize];
    	double sum = 0;
    	
    	for(int i =0; i < tournamentSize; i ++)
    	{
    		array[i] = rand.nextDouble();
    		sum += array[i];
    	}
    	
    	for(int i = 0; i < tournamentSize; i++)
    	{
    		array[i] = array[i]/sum;
    	}
    	
    	double j = 0;
    	
    	for(int i = 0; i < tournamentSize; i++)
    	{
    		j += array[i]*tournamentPopulation[i][0];
    	}
    	
    	child[0] =  j;
    	return child;
    }
// Michalewicz’s Two-Parent Geometrical Crossover
    static double[] MichalewiczGeometric(double[][] population) {
    	// We crate a small tournament, but do a wighted average, insted of chosing the best
    	int tournamentSize = 5;
    	double[][] tournamentPopulation = new double[POPULATION_SIZE][DIMENSIONS];
    	double[] child = new double[1];
    	for(int i = 0; i<tournamentSize; i++)
    	{
    		int j = rand.nextInt(POPULATION_SIZE); 
            double[] randIndividual = population[j];
            tournamentPopulation[i] = randIndividual;
    	}
    	double[] array = new double[tournamentSize];
    	double sum = 0;
    	
    	for(int i =0; i < tournamentSize; i ++)
    	{
    		array[i] = rand.nextDouble();
    		sum += array[i];
    	}
    	
    	for(int i = 0; i < tournamentSize; i++)
    	{
    		array[i] = array[i]/sum;
    	}
    	
    	double j = 0;
    	for(int i = 0; i < tournamentSize; i++)
    	{
    		j = j * Math.pow(tournamentPopulation[i][0], array[i]);
    	}
    	
    	child[0] =  j;
    	//System.out.println(j);
    	return child;
    }
// Gaussian mutation / Headless chicken
    static double[] headlessChicken(double[][] population, double[] parent) {
    	Random rand = new Random();
    
    	double[] child = new double[DIMENSIONS];
    	double u = rand.nextDouble(-10000, 10000);
    	
    	child[0] = u;
    	
    	double[] ret = new double[DIMENSIONS];
    	
    	for(int i = 0; i < child.length; i++)
    	{
    		ret[i] = (child[i] + parent[i]);
    	}
    	return ret;
    }

// ================================================================================================== MUTATION ==================================================================================================
// Self-Adaption (Lognormal)
    static double selfAdaptionLogNormal(int generation) {
    	/*;
    	double t = 1/(Math.sqrt(2));
    	
    	double a = rand.nextDouble();
    	double b = rand.nextDouble();
    	
    	ss0 = (ss0*Math.pow(Math.E, (t*a + t*b)))/MAX_GENERATIONS;
    	double ret = Math.pow(a, b);*/
    	//System.out.println("Gen"+generation);
    	Random rand = new Random();
    	
    	double j = 100.0/(generation+1);
    	double x = rand.nextDouble(-1.0,1.0);
    	return x*j;
    }
    

// ================================================================================================== REPLACEMENT STRATEGY ==================================================================================================
//Replace worst
    static double[][] replaceWorst(double[][] population, double[] solution, int i) {
    	
    	double[] worst = findWorstSolution(population);
    	i = findWorstSolutionIndex(population);
    	population[i] = solution;
        return population;
    }
    
  //Replace worst
    static double[][] replaceBest(double[][] population, double[] solution) {
    	
    	double[] best = findBestSolution(population);
    	int bestIndex = findBestSolutionIndex(population);
    	population[bestIndex] = solution;
        return population;
    }
//Replace random
    static double[][] replaceRandom(double[][] population, double[] solution, int k) {
    	
    	//Random rand = new Random();
    	k = rand.nextInt(population.length);
    	double[] ret = population[k];
    	population[k] = solution;
        return population;
    }
    
//Parent child replacement
    static double[][] parentChild(double[][] population, double[] rep, int k)
    {
    	if (fitness(rep) > fitness(population[k])) { // Compare fitness for maximizing
            population[k] = rep;
        }
    	return population;
    }
    
//Elitist
    static double[][] elitistReplace(double[][] population, double[] solution) {
    	
    	int bestIndex = findBestSolutionIndex(population);
    	int randInt = rand.nextInt(0,bestIndex);
    	int randInt1 = rand.nextInt(bestIndex+1,population.length);
    	
    	int thirdint = rand.nextInt(2);
    	
    	if(thirdint  == 1)
    	{
    		population[randInt] = solution;
    	}
    	else
    	{
    		population[randInt1] = solution;
    	}
    	//best solution does not get replaced
        return population;
    }

  //Kill tournament
    static double[][] killTournament(double[][] population, double[] solution) {
    	
    	double[][] tournPopulation = new double[5][1];
    	
    	for(int r =0; r < tournPopulation.length; r++)
    	{
    		tournPopulation[r] = randomSelection(population);
    	}

    	int worstIndex = findWorstSolutionIndex(tournPopulation);
    	population[worstIndex] = solution;
    	
        return population;
    }
    
  //Kill tournament2
    static double[][] killTournament2(double[][] population, double[] solution) {
    	
    	double[][] tournPopulation = new double[5][1];
    	
    	for(int r =0; r < tournPopulation.length; r++)
    	{
    		tournPopulation[r] = randomSelection(population);
    	}
    	
    	//kill tornament where we replace the best
    	int bestIndex = findBestSolutionIndex(tournPopulation);
    	population[bestIndex] = solution;
    	
        return population;
    }


// ================================================================================================== CROSSOVER ==================================================================================================

 // Method for binary crossover
    static double[] binaryCrossover(double[] current, double[] mutant, double Cr) {
        int D = current.length; // Dimension of the vectors which is one 
        double[] trial = new double[D]; // Trial vector
        int jrand = rand.nextInt(D); // Ensuring at least one component is taken from mutant vector
        
        for (int i = 0; i < D; i++) {
            if (rand.nextDouble() < Cr || i == jrand) {
                trial[i] = mutant[i]; // Take component from mutant
            } else {
                trial[i] = current[i]; // Take component from target
            }
        }
        return trial;
    }

    // Method for exponential crossover
    static double[] exponentialCrossover(double[] current, double[] mutant, double Cr) {
        int D = current.length; // Dimension of the vectors
        double[] trial = new double[D]; // Trial vector
        int i = rand.nextInt(D); // Random starting point for the crossover
        int L = 0; // Counter for the number of genes copied from the mutant
        
        // Copy at least one gene from mutant
        do {
            trial[i] = mutant[i]; // Copy gene from mutant to trial
            i = (i + 1) % 1; // Move to the next gene (circular)
            L++;
        } while (rand.nextDouble() < Cr && L < D);

        // Fill the rest with the target vector values
        for (int j = 0; j < D; j++) {
            if (trial[j] == 0.0) {
                trial[j] = current[j]; // Copy gene from target if it wasn't taken from mutant
            }
        }
        return trial;
    }
    
// ================================================================================================== STOPPING CONDITIONS ==================================================================================================
    // When acceptable solution is found: (-1.088598, 3.96196)
    static boolean stopCondition_AS(double[][] population) {
    	double[] x = randomSelection(population);
    	if(Math.abs((Math.abs(x[0])- Math.abs(-1.088598)))<0.1)
   		{
   			blnContinue = false;
   		}
 
        return blnContinue;
    }
    
 // When objective function deviation reaches zero
    static boolean stopCondition_OFD(double[][] population) {
    	double sum_diff = 0;
    	double sum = 0;
    	double sigma = 0;
    	for (int i = 0; i < population.length; i++) {
    		sum += population[i][0];
    	}
    	
    	double mean = sum/population.length;
    	
    	for (int i = 0; i < population.length; i++) {
    		sum_diff += Math.pow(population[i][0]-mean, 2);
    	}
    	
    	sigma = Math.sqrt(sum_diff/population.length);
    	
    	if(sigma < 0.1)
    	{
    		blnContinue = false;
    	}
 
        return blnContinue;
    }
    
 
	
    
}
