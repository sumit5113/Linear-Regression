import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * 
 */

/**
 * This method is used to calculate the linear regression using the formulae
 * given in the data sets.
 * 
 * @author sumit
 * 
 */
public class LinearRegression {

	/**
	 * This is a starter method of the linear regression process. We are
	 * implementing brute force approach for finding the different parameters
	 * measures for the obtained best fit curve data sets.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		// get the data points from the file
		PointSet pointSet = getCoordinateSystem();

		// perform the computation of best fit line
		if (pointSet != null) {
			performBruteForceComputation(pointSet);
			perfomComputation(pointSet);
		}

	}

	/**
	 * This method will perform the brute force approach for finding the best
	 * fit curves.
	 * 
	 * @param pointSet
	 */
	private static void performBruteForceComputation(PointSet pointSet) {
		// 1.find the minimum and maximum y value
		// 2.find the median/mean of x and y value
		// 3.create a y = m*b_1+b_0 line ->L1(b_1,b_0,x)
		// 4.loop for b_0 -> [minY,maxX]
		// 4.1. for each value of b_0, we can check whether sum_of_residual
		// squares for different value of sum(|L(b_1,b_0,x[i])-y[i]|^2) for i 1
		// to
		// n, we have minimum value print the final minimum Value.

		// data section initialize

		double XMedian = findMean(pointSet.getXCoordinate());
		double YMedian = findMean(pointSet.getYCoordinate());

		double minY = Collections.min(pointSet.getYCoordinate());
		double maxY = Collections.max(pointSet.getYCoordinate());

		double minX = Collections.min(pointSet.getXCoordinate());
		double maxX = Collections.max(pointSet.getXCoordinate());

		double beta_0 = 0.0;// intercept
		double beta_1 = 0.0;// slope

		double slope = 0.0;

		// getting extreme ends of beta_0 i.e [a,b]
		double beta_0_high = (maxY - getSlope(XMedian, YMedian, minX, maxY)
				* minX);
		double beta_0_low = (minY - getSlope(XMedian, YMedian, minX, minY)
				* minX);

		double beta_1_low = getSlope(0, beta_0_low, maxX, minY);
		double beta_1_high = getSlope(XMedian, YMedian, 0,beta_0_high);

		double yDiff = 0.0;
		int length = pointSet.getXCoordinate().size();
		List<Double> XCoord = pointSet.getXCoordinate();
		List<Double> YCoord = pointSet.getYCoordinate();
		double sum = Double.MAX_VALUE;

		// this is Line 1 -> which varies in y-intercept
		for (double d = beta_0_low; d < beta_0_high; d = d + .0001) {
			// equation of new line
			// slope between (0,d) and (XMedain,Ymedian)
			slope = getSlope(0, d, XMedian, YMedian);

			double tempSum = 0.0;

			if (Double.isInfinite(slope)) {
				// vertical line, compute only for (Xmedian-x[i])^2
				System.out.println("slope turns out to be infinite::");
				continue;
			}
			// computing the vertical sum
			for (int index = 0; index < length; index++) {
				yDiff = (YCoord.get(index) - (d + slope * XCoord.get(index)));
				tempSum = tempSum + yDiff * yDiff;
			}
			if (tempSum < sum) {
				sum = tempSum;
				beta_0 = d;
				beta_1 = slope;
			}
		}

		// this is line 2 which varies in slope
		slope = getSlope(XMedian, YMedian, minX, minY);
		double xIntercept = -(YMedian - slope * XMedian) / slope;
		double xMaxIntercept = -(YMedian - getSlope(XMedian, YMedian, maxX,
				minY) * XMedian)
				/ getSlope(XMedian, YMedian, maxX, minY);
		
		beta_1_high=Math.max(Math.atan(getSlope(XMedian,YMedian, maxX,minY)),Math.atan(beta_1_high));
		
		double c = 0.0d;

		while (xIntercept <= xMaxIntercept) {

			slope = getSlope(XMedian, YMedian, xIntercept, 0);

			
			double tempSum = 0.0;

			if (Double.isInfinite(slope)) {
				// vertical line, compute only for (Xmedian-x[i])^2
				System.out.println("slope turns out to be infinite::");
				continue;
			}

			c = YMedian - slope * XMedian;

			for (int index = 0; index < length; index++) {
				yDiff = (YCoord.get(index) - (c + slope * XCoord.get(index)));
				// yDiff = (XCoord.get(index) - ( (YCoord.get(index)-c)/slope));
				tempSum = tempSum + yDiff * yDiff;
			}
			if (tempSum < sum) {
				sum = tempSum;
				beta_0 = YMedian - slope * XMedian;// mx+C
				beta_1 = slope;
			}

			xIntercept += 0.00001;
		}

		System.out.println("=============Brute Force Analysis================");
		System.out.println("Range Of Beta_0=[" + beta_0_low + "," + beta_0_high
				+ "]");
		System.out.println("Range Of Beta_1=[" + beta_1_low + "," + Math.tan(beta_1_high)
				+ "]");
		System.out.println("");
		System.out.println("brute force calculation ::intercept :: " + beta_0
				+ ", slope:: " + beta_1);
		System.out.println("Best fitting line:: y = " + (beta_0) + " +("
				+ beta_1 + ") * x ");
		System.out.println("Residual sum of squares : " + Math.sqrt(sum));
	}

	/**
	 * Method to obtain the slope of two pint.
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return
	 */
	private static double getSlope(double x1, double y1, double x2, double y2) {
		return (y1 - y2) / (x1 - x2);

	}

	/**
	 * Method to compute the mean of the two line 
	 * @param numericData
	 * @return
	 */
	private static double findMean(List<Double> numericData) {
		double XMean = 0.0;
		int length = numericData.size();
		for (int index = 0; index < length; index++) {
			XMean = XMean + numericData.get(index);

		}
		return XMean / length;
	}

	/**
	 * 
	 * @param pointSet
	 */
	private static void perfomComputation(PointSet pointSet) {
		// I.Pass for mean computation , with n points O(n)
		// 1.calculate the mean of X and Y coordinate
		double XMean = 0.0;
		double YMean = 0.0;
		double sumOfX2 = 0.0;
		int length = pointSet.getXCoordinate().size();
		List<Double> XCoord = pointSet.getXCoordinate();
		List<Double> YCoord = pointSet.getYCoordinate();
		// System.out.println(YCoord);
		for (int index = 0; index < length; index++) {
			XMean = XMean + XCoord.get(index);
			YMean = YMean + YCoord.get(index);
			sumOfX2 = sumOfX2 + (XCoord.get(index) * XCoord.get(index));
		}

		XMean = XMean / length;
		YMean = YMean / length;
		// System.out.println(YMean);

		// II.pass calculation beta_0 and beta_1
		// 2.calculate xxProduct,YYProduct,XYProduct
		double XXProduct = 0.0, YYProduct = 0.0, XYProduct = 0.0;

		for (int index = 0; index < length; index++) {
			XXProduct += (XCoord.get(index) - XMean)
					* (XCoord.get(index) - XMean);
			YYProduct += (YCoord.get(index) - YMean)
					* (YCoord.get(index) - YMean);
			XYProduct += (XCoord.get(index) - XMean)
					* (YCoord.get(index) - YMean);
		}

		// 3.calculating beta_0 and beta_1, coefficients
		double beta_1 = XYProduct / XXProduct; // coefficient of x , i.e
												// regression coefficient
		double beta_0 = YMean - XMean * beta_1; //

		// III. measuring the goodness or accuracy of the best fitting line
		// 4.calculating residual sum of squares, regression sum of squares,
		// coefficient of determination,stand error in beta_0,beta_1,
		double fit = 0.0;
		// double R_square = 0.0;
		double residual_sum_squares = 0.0;
		double regression_sum_suquare = 0.0;
		// double stnd_err_beta_0 = 0.0;
		// double stnd_err_beta_1 = 0.0;
		double diff = length - 2;

		double residual;
		double min_residual = Double.MAX_VALUE;
		double max_residual = Double.MIN_VALUE;
		double residual_median = 0.0;
		double residual_1_quantile = 0.0;
		double residual_3_quantile = 0.0;
		List<Double> residuals = new ArrayList<>();

		for (int index = 0; index < length; index++) {
			fit = (XCoord.get(index) * beta_1 + beta_0);
			residual = -(fit - YCoord.get(index));
			residual_sum_squares += ((fit - YCoord.get(index)) * (fit - YCoord
					.get(index)));
			regression_sum_suquare += ((fit - YMean) * (fit - YMean));

			min_residual = Math.min(min_residual, residual);
			max_residual = Math.max(max_residual, residual);
			residuals.add(residual);
		}

		Collections.sort(residuals);
		// 1st quartile computation,median,3rd quartile
		residual_median = findMedianInSortedList(residuals);
		residual_1_quantile = findMedianInSortedList(residuals.subList(0,
				residuals.size() / 2));
		residual_3_quantile = findMedianInSortedList(residuals.subList(
				residuals.size() / 2, residuals.size()));

		// analysis result print
		System.out.println("=======Statistical Analysis===============");
		System.out.println("Best fitting line:: y = " + (beta_0) + " +("
				+ beta_1 + ") * x ");
		System.out.println("Residuals parameter :: ");
		System.out.println("min ::" + min_residual + ", median::"
				+ residual_median + ", max::" + max_residual + ", 1QT  :"
				+ residual_1_quantile + ", 3QT :" + residual_3_quantile);
		System.out.println("Coefficient of determination R^2 :: "
				+ (regression_sum_suquare / YYProduct));
		System.out.println("Residual sum of squares:: "
				+ Math.sqrt(residual_sum_squares));
		System.out
				.println("Residual Standard error in computation of best fit line:: "
						+ Math.sqrt(residual_sum_squares / diff)
						+ " on "
						+ diff + " of freedom");
		System.out.println("Standard error in beta_1 computation ::"
				+ Math.sqrt((residual_sum_squares / diff) / XXProduct));
		System.out
				.println("Standard error in beta_0 computation ::"
						+ Math.sqrt((sumOfX2 * (residual_sum_squares / diff) / (XXProduct * length))));

	}

	private static double findMedianInSortedList(List<Double> residuals) {
		int size = residuals.size();
		double median = 0.0;
		if (size % 2 == 0) {
			median = (residuals.get(size / 2 - 1) + residuals.get(size / 2)) / 2;

		} else {
			median = residuals.get(size / 2);
		}

		return median;
	}

	/**
	 * Read the coordinate system from the input file and prepare the coordinate
	 * system.
	 * 
	 * @return
	 */
	private static PointSet getCoordinateSystem() {
		try {
			String[] inputs = usageCommandArguments();
			List<String> allPairPoints = Files.readAllLines(
					Paths.get(inputs[0]), Charset.defaultCharset());
			// split all strings
			if (inputs[1].length() == 0
					|| inputs[1].replaceAll("[\\s]", "").length() == 0) {
				inputs[1] = "[\\s]+";
			}

			List<Double> X = new ArrayList<>();
			List<Double> Y = new ArrayList<>();

			// get the x and y coordinate values from the passed file name
			String[] text;
			for (String temp : allPairPoints) {
				text = temp.split(inputs[1]);
				if (text.length != 2) {
					throw new IllegalStateException(
							"Input file conatins invalid data: [" + temp + "]");
				}
				X.add(Double.valueOf(text[0]));
				Y.add(Math.log10(Double.valueOf(text[1])));
			}

			return new PointSet(X, Y);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("Error has occured:: [details ="
					+ e.getMessage() + "]");
			e.printStackTrace();
		}

		return null;
	}

	private static String[] usageCommandArguments() throws IOException {
		String[] input = new String[2];
		BufferedReader bfrdReader = new BufferedReader(new InputStreamReader(
				System.in));
		while (true) {
			System.out.println("**********Menu*********");
			System.out
					.println("Following are the input options: \n[1. File Name [F]]\n[2.Exit [E]]");
			System.out.println("Please enter the choice");
			System.out.print("Input Choice:");

			input[0] = bfrdReader.readLine();

			switch (input[0]) {
			case "f":
			case "F":
				System.out.print("Enter File Name (file path):");
				input[0] = bfrdReader.readLine();
				System.out
						.print("Enter data separator [if white space hit enter only]:");
				input[1] = bfrdReader.readLine();
				break;
			case "e":
			case "E":
				bfrdReader.close();
				System.exit(0);
				break;
			default:
				System.out.println("Invalid option");
				continue;
			}
			bfrdReader.close();
			return input;
		}
	}

	/**
	 * This class will hold the point sets, XY coordinate system
	 * 
	 * @author sumit
	 * 
	 */
	static class PointSet {
		List<Double> XCoordinate;

		List<Double> getXCoordinate() {
			return XCoordinate;
		}

		void setXCoordinate(List<Double> xCoordinate) {
			XCoordinate = xCoordinate;
		}

		List<Double> getYCoordinate() {
			return YCoordinate;
		}

		void setYCoordinate(List<Double> yCoordinate) {
			YCoordinate = yCoordinate;
		}

		List<Double> YCoordinate;

		PointSet(List<Double> pXCordinate, List<Double> pYCordinate) {
			if (pXCordinate == null || pYCordinate == null
					|| pXCordinate.size() != pYCordinate.size())
				throw new IllegalArgumentException(
						"Invalid Point Set Parameters");

			this.XCoordinate = pXCordinate;
			this.YCoordinate = pYCordinate;
		}

	}

}
