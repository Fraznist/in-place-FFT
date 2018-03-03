package main;

import java.util.ArrayList;

import org.apache.commons.math3.complex.Complex;

public class Polynomial {
	private static final double err = 0.000000000001;

	ArrayList<Complex> array;
	polyStrat representation;
	int degree = 0;
	
	// Constructs the Polynomial from a list of coefficients
	// index i represents the coeff for x^i.
	public Polynomial(double[] coeffs) {
		array = new ArrayList<Complex>();
		for (int i = 0; i < coeffs.length; i++) {	// copy the array into the list
			array.add(new Complex(coeffs[i]));
			if (coeffs[i] != 0) degree = i;	// track the degree of the polynomial
		}
		representation = new Coefficient();
	}
	public Polynomial(ArrayList<Complex> clist, int deg) {
		array = clist;
		degree = deg;
		representation = new Sample();
	}
	
	public Complex get(int i) {
		return array.get(i);
	}
	
	public void set(int i, Complex c) {
		array.set(i, c);
	}
	
	public Complex eval(Complex c) {
		return representation.eval(this, c);
	}
	
	public Polynomial add(Polynomial p) {
		return representation.add(this, p);
	}
	
	public void addInPlace(Polynomial p) {
		representation.addInPlace(this, p);
	}
	
	public Polynomial mult(Polynomial p) {
		return representation.mult(this, p);
	}
	
	public void multInPlace(Polynomial p) {
		representation.multInPlace(this, p);
	}
	
	private Polynomial changeRepresentation(int targetSamples) {
		return representation.changeRepresentation(this, targetSamples);
	}
	
	public String toString() {
		return representation.string(this);
	}
	
	public abstract interface polyStrat {
		public Complex eval(Polynomial p, Complex c);
		public Polynomial add(Polynomial p1, Polynomial p2);
		public void addInPlace(Polynomial p1, Polynomial p2);
		public Polynomial mult(Polynomial p1, Polynomial p2);
		public void multInPlace(Polynomial p1, Polynomial p2);
		public Polynomial changeRepresentation(Polynomial p, int targetSamples);
		public String string(Polynomial p);
		
		// For a list of size 2^n, represent each index of array with a n bit number
		// swap contents of index i and reverse(i)
		// Turns out by swapping elements this way we get the same order with
		// what we would get by dividing the polynomial into even and odd parts recursively
		public default void orderArrayForFFT(Polynomial p) {
			int size = p.array.size();
			int subset = 0;
			while (size > 1) {
				size *= 0.5;
				subset++;
			}
			for (int i = 0; i < p.array.size(); i++) {
				int x = reversePartially(i, subset);
				if (x > i) {
					Complex temp = p.get(i);
					p.set(i, p.get(x));
					p.set(x, temp);
				}
			}
		}
		
		// reverse the bit values of x, treating it as if it consists of subset bits
		public default int reversePartially(int x, int subset) {
			int y = 0;
			for (int i = 0; i < subset; i++) {
				y += ((x&1)<<subset-i-1);
				x >>= 1;
			}
			return y;
		}
		
		public default Complex roundComplex(Complex c) {
			if (Math.abs(Math.round(c.getReal()) - c.getReal()) < err)
				c = new Complex(Math.round(c.getReal()), c.getImaginary());
			if (Math.abs(Math.round(c.getImaginary()) - c.getImaginary()) < err)
				c = new Complex(c.getReal(), Math.round(c.getImaginary()));
			return c;
		}	
	}
	
	public class Coefficient implements polyStrat {

		@Override
		public Complex eval(Polynomial p, Complex c) {
			Complex xpow = Complex.ONE;
			Complex sum = Complex.ZERO;
			for (int i = 0; i <= degree; i++) {
				sum = sum.add(xpow.multiply(p.get(i)));
				xpow = xpow.multiply(c);
			}
			return sum;
		}

		@Override
		public Polynomial add(Polynomial p1, Polynomial p2) {
			if (p1.degree > p2.degree) 
				return addHelper(p2, p1);
			else 
				return addHelper(p1, p2);
		}
		
		@Override
		public void addInPlace(Polynomial p1, Polynomial p2) {
			if (p1.degree > p2.degree) 
				addHelperInPlace(p2, p1, p1);
			else 
				addHelperInPlace(p1, p2, p1);	
		}

		@Override
		public Polynomial mult(Polynomial p1, Polynomial p2) {
			// TODO Auto-generated method stub
			int n = 1;
			int multDeg = p1.degree + p2.degree;
			while (n < multDeg + 1) // Determine size of arraylist = 2^n
				n *= 2;
			p1.changeRepresentation(n);
			p2.changeRepresentation(n);
			
			Polynomial mult = p1.mult(p2);
			
			p1.changeRepresentation(0);
			p2.changeRepresentation(0);
			
			return mult.changeRepresentation(0);
		}
		
		@Override
		public void multInPlace(Polynomial p1, Polynomial p2) {
			int n = 1;
			int multDeg = p1.degree + p2.degree;
			while (n < multDeg + 1) // Determine size of arraylist = 2^n
				n *= 2;
			p1.changeRepresentation(n);
			p2.changeRepresentation(n);
			
			p1.multInPlace(p2);
			
			p1.changeRepresentation(0);
			p2.changeRepresentation(0);
		}

		@Override
		public String string(Polynomial p) {
			String sum = "";
			for (int i = 0; i < p.degree; i++)
				sum += printComplex(p.get(i)) + "x^" + i + " + ";
			sum += printComplex(p.get(p.degree)) + "x^" + p.degree;
			return sum;
		}
		
		// Convert contents of coefficient array to Sample representation,
		// aim to sample at "targetSamples" different points. 
		public Polynomial changeRepresentation(Polynomial p, int targetSamples) {
			// Divide phase of the divide-conquer-merge trifecta
			// dummy 0.0 coeffs are added to the end as necessary to have size = 2^n
			// This helps for FFT.
			for (int i = p.array.size(); i < targetSamples; i++) 
				p.array.add(Complex.ZERO);
			
			// Reorder coefficients
			orderArrayForFFT(p);
			
			// Conquer phase, take samples of every degree 0 polynomial
			// at x=1, which is the 1st Roots of unity of x=1
			// Since all polynomials are degree 0, any coefficient is already
			// equal to its respective samples at any point, no need to do anything.
			
			// Merge phase, thanks to the reordering at the divide phase,
			// Samples of the same x value of the two merging polynomials
			// will always be neighbors, and overwrite their own values. Magiq.
			
			// The resulting list of sample values will also be sorted ascending
			// with respect to the angle of the polar representation of the value
			// For instance {0, 45, 90, 135, 180, 225, 270, 315} for input size 8
			// Where index 2 represents the sample value for x = S(90).
			
			// Thinking about the recursion tree for the FFT algorithm,
			// Each iteration of the outer loop (i) represents a level of the recursion tree
			for (int i = 1; i < p.array.size(); i *= 2) {	
				Complex delta = new Complex(Math.cos(Math.toRadians(180 / i)), Math.sin(Math.toRadians(180 / i)));
				
				// Each iteration of loop j represents 2 merging nodes of the recursion tree
				// i is the size (degree + 1) of a single polynomial at this level.
				for (int j = 0; j < p.array.size(); j += 2*i) {	
					Complex samplePointX = Complex.ONE;
					Complex samplePointRootX = new Complex(-1);
					
					// Finally, loop k actually merges j'th and 1+j'th polynomials,
					// Both polynomials of size i
					for (int k = j; k < j + i; k++) {
						Complex even = p.get(k);
						Complex odd = p.get(k + i);
						
						p.set(k, even.add(odd.multiply(samplePointX)));
						p.set(k + i, even.add(odd.multiply(samplePointRootX)));
						
						samplePointX = samplePointX.multiply(delta);
						samplePointRootX = samplePointRootX.multiply(delta);
					}
				}
			}
			
			// Finally, for my leisure, round values that are too close to an integer,
			// which must be caused due to floating point arithmetics
			for (int i = 0; i < p.array.size(); i++) {
				p.set(i, roundComplex(p.get(i)));
			}
			
			p.representation = new Sample();
			
			return p;
		}
		
		private Polynomial addHelper(Polynomial small, Polynomial large) {
			double[] d = new double[large.degree + 1];
			for (int i = 0; i <= small.degree; i++) 
				d[i] = small.get(i).getReal() + large.get(i).getReal();
			for (int i = small.degree + 1; i <= large.degree; i++) 
				d[i] = large.get(i).getReal();
			return new Polynomial(d);
		}
		
		private void addHelperInPlace(Polynomial small, Polynomial large, Polynomial target) {
			for (int i = 0; i <= small.degree; i++) 
				target.set(i, small.get(i).add(large.get(i)));
			for (int i = small.degree + 1; i <= large.degree; i++) 
				target.array.add(large.get(i));
			target.degree = large.degree;
		}
		
		private String printComplex(Complex c) {
			if (c.getImaginary() == 0.0)
				return c.getReal() + "";
			else
				return c.toString();
		}
		
	}
	
	public class Sample implements polyStrat {

		@Override
		public Complex eval(Polynomial p, Complex c) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public Polynomial add(Polynomial p1, Polynomial p2) {			
			if (p1.array.size() == p2.array.size()) {
				ArrayList<Complex> list = new ArrayList<Complex>();
				for (int i = 0; i < p1.array.size(); i++) 
					list.add(p1.get(i).add(p2.get(i)));
				if (p1.degree > p2.degree)
					return new Polynomial(list, p1.degree);
				else 
					return new Polynomial(list, p2.degree);
			}
			else return null;
		}
		
		@Override
		public void addInPlace(Polynomial p1, Polynomial p2) {
			for (int i = 0; i < p1.array.size(); i++) 
				p1.set(i, p1.get(i).add(p2.get(i)));
			p1.degree = Math.max(p1.degree, p2.degree);
		}

		@Override
		public Polynomial mult(Polynomial p1, Polynomial p2) {
			// multiply each sample, the arguments should have the same 
			// arraylist size, which should be sufficient to hold at least
			// p1.degree + p2.degree + 1 samples. This will be handled
			// by the FFT operation.
			ArrayList<Complex> list = new ArrayList<Complex>();
			for (int i = 0; i < p1.array.size(); i++) 
				list.add(p1.get(i).multiply(p2.get(i)));
			return new Polynomial(list, p1.degree + p2.degree);
		}
		
		public void multInPlace(Polynomial p1, Polynomial p2) {
			for (int i = 0; i < p1.array.size(); i++) 
				p1.array.set(i, p1.get(i).multiply(p2.get(i)));
			
			p1.degree += p2.degree;
		}
		
		// arg targetSamples isn't needed, why am I even using strategy pattern
		// It's too late to change.
		public Polynomial changeRepresentation(Polynomial p, int targetSampes) {
			// Divide phase
			// Analogous to the complex merge phase of Coefficient.toSample(),
			// the inverse of toSample() needs to have a complex divide phase
			
			// i is half the size of a polynomial at this level,
			// i.e. i is the size of the odd or the even portion
			for (int i = (int) (p.array.size() * 0.5); i >= 1; i *= 0.5) {	
				Complex delta = new Complex(Math.cos(Math.toRadians(180 / i)), Math.sin(Math.toRadians(180 / i)));
				
				// Each iteration of loop j represents a dividing node of the recursion tree
				// i is half the size (degree + 1) of a single polynomial at this level.
				for (int j = 0; j < p.array.size(); j += 2*i) {	
					Complex x1 = Complex.ONE;
					Complex x2 = new Complex(-1);
					
					// Finally, loop k actually divides j'th polynomial,
					// Both polynomials of size i
					for (int k = j; k < j + i; k++) {
						Complex s1 = p.get(k);
						Complex s2 = p.get(k + i);
						
						Complex odd = s1.subtract(s2).divide(x1.subtract(x2));
						Complex even = s1.subtract(x1.multiply(odd));
						
						p.set(k, even);
						p.set(k + i, odd);
						
						x1 = x1.multiply(delta);
						x2 = x2.multiply(delta);
					}
				}
			}
			
			// Conquer phase
			// At this point, we have n polynomials of degree 0 in sample representation
			// Convert them to coefficient representation
			// i.e. do nothing, same reason with toSample() conquer phase.
			
			// Merge phase
			// Invert the special ordering done in toSample().
			orderArrayForFFT(p);
			
			// Due to using 2^n samples, the array probably has more elements
			// than degree + 1. Indices bigger than degree hold values that are equal to 0
			// They are obviously useless as coefficients. We remove them.
			for (int i = p.array.size() -1; i > p.degree ; i--) 
				// removing last element each time is more efficient.
				p.array.remove(i);
			
			// Finally, for my leisure, round values that are too close to an integer,
			// which must be caused due to floating point arithmetics
			for (int i = 0; i < p.array.size(); i++) {
				p.set(i, roundComplex(p.get(i)));
			}
			
			p.representation = new Coefficient();
			
			return p;
		}

		@Override
		public String string(Polynomial p) {
			String sum = "Degree: " + p.degree + "\n{";
			for (int i = 0; i < p.array.size() -1; i++)
				sum += p.get(i) + ", ";
			sum += p.get(p.array.size() -1) + "}";
			return sum;
		}
		
	}
}
