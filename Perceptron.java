import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.stream.Collectors;

import java.util.Collections;
import java.io.FileWriter;
import java.io.PrintWriter;

public class Perceptron {
    double[] weights;

    public static void main(String[] args) throws FileNotFoundException {
        Scanner input = new Scanner(new File("C:\\Users\\Felix\\Documents\\Uni\\Artificial Intelligence\\Single Perceptron\\test_data.txt"));
//        Scanner input = new Scanner(System.in);
        Perceptron p = new Perceptron(3);
        p.initialize();
        List<String> trainingDataset = new ArrayList<>();
        List<String> dataset = new ArrayList<>();
        
        String nextLine = input.nextLine();
        while(input.hasNextLine() && !nextLine.equals("0,0,0")) {
                trainingDataset.add(nextLine);
                nextLine = input.nextLine();  
        }
        while(input.hasNextLine()) {
            nextLine = input.nextLine();
            dataset.add(nextLine);
        };

        //p.printStringList(dataset);
        // normalize the datset:
        List<List<Double>> normTestData = newNormalize(dataset);
        List<List<Double>> normTrainData = newNormalize(trainingDataset);
        
        
        List<Double> errorProp = p.learn(normTrainData);
        System.out.print("Weight 1: " + p.weights[0] + "\n");
        System.out.print("Weight 2: "  +p.weights[1] + "\n");
        System.out.print("Bias:     " + p.weights[2] + "\n");
        System.out.print("\n");
        // do enhamcement:
        // calculate first error:
        double firstError = 0;
        for (int i = 0; i < errorProp.size(); i++){
            firstError += errorProp.get(i);
        }
        
        //List<Double> errorSums = p.enhancement(10, normTrainData, firstError);
        /*
        int iterationStep;
        for (int i = 0; i < errorSums.size(); i++){
            iterationStep = i + 1;
            System.out.print("Error of Iteration " + iterationStep + ": " + errorSums.get(i) + "\n");
        }*/

        List<Double> results = p.output(normTestData);
        
        //p.printList(results);

        // now compare and see if your results are good:
        Scanner resultFile = new Scanner(new File("C:\\Users\\Felix\\Documents\\Uni\\Artificial Intelligence\\Single Perceptron\\result_data.txt"));
        List<String> resultData = new ArrayList<>();
        
        nextLine = resultFile.nextLine();
        while(resultFile.hasNextLine()) {
            resultData.add(nextLine);
            nextLine = resultFile.nextLine();
        };
        resultData.add(nextLine); // add the last line

        List<Double> realResults = parseResults(resultData);
        
        //p.printList(realResults);
        // now compare:
        double errorCounter = results.size();
        for (int i = 0; i < results.size(); i++){
            if (!realResults.get(i).equals(results.get(i))) errorCounter--;
        }

        double percentage = errorCounter / results.size() * 100;
        System.out.print("Your perceptron recognized " + errorCounter +
        " out of " + results.size() + " tests! \nSuccess rate: " + percentage  + " %");
    }
    
	// passt!
    public Perceptron(int n) {
        weights = new double[n];		// wieso machst du denn hier einmal als uebergabe 2 und dann erstellst du aber 3? hab ich mal geaendert ;)
    }

    public static List<Double> parseResults (List<String> test_results) {
        List<Double> input = new ArrayList<>();
        for(int i = 0; i < test_results.size(); i++) {
			input.add(Double.parseDouble(test_results.get(i)));
        }
        return input;
    }

    public static List<List<Double>> noNormalize (List<String> dataset) {
        List<List<Double>> input = new ArrayList<List<Double>>();
        for(int i = 0; i < dataset.size(); i++) {
            String[] line = dataset.get(i).split(",");
			input.add(Arrays.stream(line).map(Double::valueOf).collect(Collectors.toCollection(ArrayList::new)));
        }
        return input;
    }

    
    public static List<List<Double>> newNormalize (List<String> dataset) {
        List<List<Double>> input = new ArrayList<List<Double>>();
        for(int i = 0; i < dataset.size(); i++) {
            String[] line = dataset.get(i).split(",");
			input.add(Arrays.stream(line).map(Double::valueOf).collect(Collectors.toCollection(ArrayList::new)));
        }
        double maximum = 0;
        double minimum = 0;
        double firstEntry = 0;
        double secondEntry = 0;

        for(int i = 0; i < input.size(); i++){
            firstEntry = input.get(i).get(0);
            if (firstEntry > maximum) maximum = firstEntry;
            if (firstEntry < minimum) minimum = firstEntry;
            secondEntry = input.get(i).get(1);
            if (secondEntry > maximum) maximum = secondEntry;
            if (secondEntry < minimum) minimum = secondEntry;   
        }
        for (int i = 0; i < input.size(); i++) {
            input.get(i).set(0, 2*(input.get(i).get(0) - minimum) / (maximum - minimum) - 1);
            input.get(i).set(1, 2*(input.get(i).get(1) - minimum) / (maximum - minimum) - 1);
        }
        return input;
    }

    // print array list to console:
    public void printList(List<Double> dataset) {
        System.out.print("-------------------------------------------\n");
        for (int i = 0; i < dataset.size(); i++) {
            System.out.print("List entry " + i + " value: [" + dataset.get(i) + "]\n");
        }
        System.out.print("-------------------------------------------");
    }
    public void printStringList(List<String> listToPrint){
        System.out.print("-------------------------------------------\n");
        for (int i = 0; i < listToPrint.size(); i++){
            System.out.print("List entry " + i + " value: [" + listToPrint.get(i) + "]\n");
        }
        System.out.print("-------------------------------------------");
    }

    public static void printListList(List<List<Double>> listToPrint){
        List<Double> currentList;
        System.out.print("-------------------------------------------\n");
        for (int i = 0; i < listToPrint.size(); i++){
            currentList = listToPrint.get(i);
            System.out.print("List entry " + i + " value: [");
            for (int j = 0; j < 2; j++){
                System.out.print(currentList.get(j) + "   ");
            }
            System.out.print("]\n"); 
        }
        System.out.print("-------------------------------------------\n");
    }

    // Enhancement function
    public List<Double> enhancement(int n, List<List<Double>> trainingData, double firstError){
        List<Double> currentErrorList = new ArrayList<Double>();
        List<Double> errorSum = new ArrayList<Double>();
        double currentErrorSum = 0;
        errorSum.add(firstError);
        for (int i = 0; i < n; i++){
            Collections.shuffle(trainingData);
            currentErrorList = learn(trainingData);
            for (int j = 0; j < trainingData.size(); j++){
                currentErrorSum += currentErrorList.get(j);
            }
            errorSum.add(currentErrorSum);
            currentErrorSum = 0;
        }
        return errorSum;
    }

	// passt!
    public List<List<Double>> normalize (List<String> dataset) {
        List<List<Double>> input = new ArrayList<List<Double>>();
        for(int i = 0; i < dataset.size(); i++) {
            String[] line = dataset.get(i).split(",");
			input.add(Arrays.stream(line).map(Double::valueOf).collect(Collectors.toCollection(ArrayList::new)));
        }

        // calculate vector sum
        double maxVal = 0;
        double maxValSquared = 0;

        for (int i = 0; i < input.size(); i++){
            maxVal += input.get(i).get(0) + input.get(i).get(1);
            maxValSquared += Math.pow(input.get(i).get(0), 2) + Math.pow(input.get(i).get(1), 2);
        }
        // calculate mean
        double mean = maxVal / (input.size()*2);
		double var = maxValSquared / (input.size()*2);

        // calculate standard deviation
        double stddev = Math.sqrt(var - Math.pow(mean, 2));

        // normalize values
        for (int i = 0; i < input.size(); i++) {
            input.get(i).set(0, (input.get(i).get(0) - mean) / stddev);
            input.get(i).set(1, (input.get(i).get(1) - mean) / stddev);
        }
        return input;
    } 
    
    public List<Double> learn(List<List<Double>> trainingDataset) {
        double dW_i = 0;
        double out = 0;
        double net = 0;
        List<Double> errorList = new ArrayList<>();
        
        for(int i = 0; i < trainingDataset.size(); i++){
            List<Double> pattern = trainingDataset.get(i);	// get the current line
            double target = pattern.get(2);		// get the current target
            pattern.set(2,1.0);		// set third entry to 1 in order to use a third weight to simulate bias
            out = feedForward(pattern);  // calculate the output of the current line tanh(current_line)
            //System.out.print(out + "\n");
           
			// calculate deltaW: also ich hab jetzt hier mal noch den error mit 1/2 multipliziert, weil man das so macht, dann faellt das 2 weg
			// und dann hab ich das net mal rausgeloescht, weil da machst du irgendwie was doppelt, das brauchst du eigentlich gar nicht.
			// error_i(t) = 0.5 * (target_i(t) - out_i(t))^2
			// delta_error_i(t)/delta_out_i(t) = -(target - out)
			// delta_out_i(t)/delta_w_i(t) = P_i(t) * 1 / (cosh^2(out_i(t))
             
            // calculate argument:
            
            for (int j = 0; j < 3; j++){
                net += pattern.get(j) * weights[j];
            }
            //System.out.print(net + "\n");
            //System.out.print(pattern + "\n");
            for (int j = 0; j < 3; j++){ 
                dW_i =  - 0.0002  * (target - out) * derivativeTanh(net) * pattern.get(j);
                //System.out.print(dW_i + "\n");
                weights[j] -= dW_i;
            }     
            // print error after each iteration
            double newOut = feedForward(pattern);
            //System.out.print(out + " || " + newOut + " || " + target + "\n");
            double error = .5 * Math.pow((target - out),2);
            errorList.add(error);
            //System.out.print(error + "\n");
            net = 0;
        }
        
        return errorList;   
    }
    
	// passt!
    public double derivativeTanh (double x) {
        return 1/(Math.pow(Math.cosh(x), 2));
    }
    
	// passt!
    public double feedForward(List<Double> pattern) {
        double out = 0;
        for (int i = 0; i < 2; i++) {						// hab hier mal auf 2 geaendert, macht mehr Sinn in dem Fall, da es ja immer 2 sind in deiner Aufgabe
            out += pattern.get(i) * weights[i];
        }
        out += weights[2];
        out = Math.tanh(out);
        if(out >= 0){
            return 1;
        }
        else{
            return -1;
        }
        //return Math.tanh(out);
    }
    
	// passt!
    public void initialize() {
        
        Random r = new Random();
        double rangeMin = -1000;
        double rangeMax = +1000;
        for(int i = 0; i < 3; i++) {			// wie oben....
            weights[i] = (rangeMin + (rangeMax - rangeMin) * r.nextDouble()) / 100000;
        }
        //weights[0] = -0.003;
        //weights[1] = +0.002;
        //weights[2] = -0.001;       
    }
    
    public List<Double> output(List<List<Double>> dataset) {
        List<Double> resultList = new ArrayList<>();
       
        for (int i = 0; i < dataset.size(); i++) {
            List<Double> pattern = dataset.get(i);
            //printList(pattern);
            pattern.add(1.0);
            //printList(pattern);
            double out = feedForward(pattern); 
            /*if (out > 0.0) {
                System.out.println("+" + Math.signum(out));
            }
            else {
                System.out.println(Math.signum(out));
            } */
            //System.out.print(out + "\n");
            resultList.add(out);
        }
        return resultList;
        
    }
}
