import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.Formatter;
import java.util.HashMap;
import java.util.InputMismatchException;
import java.util.Map;
import java.util.Scanner;

public class Parser {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String folder = "blarge3bsv3dyn50";
		Scanner scanner;
		Map<Double, Double[]> fields = new HashMap<Double, Double[]>();
		Double[] baseStations = new Double[6];
		try {
			
			scanner = new Scanner(new FileInputStream("batch9/" + folder + "/" + folder + "MmWaveSinrTime.txt"));
			while(scanner.hasNext()) {
				double time = scanner.nextDouble();
				scanner.next();
				int baseStation = scanner.nextInt();
				try {
					double sinr = scanner.nextDouble();
					if(!fields.containsKey(time)) {
						Double[] toPut = new Double[6];
						toPut[baseStation-2] = sinr;
						fields.put(time, toPut);
					} else {
						fields.get(time)[baseStation-2] = sinr;
					}
				} catch (InputMismatchException e) {
					scanner.next();
				}
			
			}
			Formatter fileMaker = new Formatter(new FileOutputStream("batch9/" + folder + "/" + folder + "excelreadycsv.txt"));
			for(Map.Entry<Double, Double[]> entry: fields.entrySet()) {
				double[] lastBase = new double[6];
				fileMaker.format("%f", entry.getKey());
				Double[] array = entry.getValue();
				for(int i = 0; i < 6; i++) {
					double curSinr;
					if(!(array[i] == null)) {
						curSinr = array[i].doubleValue();
					} else {
						curSinr = 0.0;
					}
					
					if(curSinr == 0.0) {
						curSinr = lastBase[i];
					}
					lastBase[i] = curSinr;
					fileMaker.format(", %f", curSinr);
		
				}
				fileMaker.format("%n");
				fileMaker.flush();
			}
			fileMaker.close();
			System.out.print("Done");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}

}
