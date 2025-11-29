import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}
    /**
     * Test the PolynomFromPoints function - one test for each possible outcome (either ^2 or ^3)
     */
    @Test
    public void testPolynomFromTwoPoints() {
        // Line passing through points (1,3) and (3,7)
        double[] xx = {1, 3};
        double[] yy = {3, 7};

        double[] p = Ex1.PolynomFromPoints(xx, yy);

        // Expected line: f(x) = 2x + 1
        assertEquals(1, p[0], 1.0/1000);  // b
        assertEquals(2, p[1], 1.0/1000);  // a

        // Check values at sample points
        assertEquals(3, Ex1.f(p, 1), 1.0/1000);
        assertEquals(7, Ex1.f(p, 3), 1.0/1000);
    }
    @Test
    void testPolynomFromThreePoints() {
        // Points from the quadratic f(x) = 2x^2 - 3x + 5
        double[] xx = {0, 1, 2};
        double[] yy = {5, 2*1*1 - 3*1 + 5, 2*4 - 6 + 5}; // corresponding f(x) coordinates for each x in xx

        double[] p = Ex1.PolynomFromPoints(xx, yy);

        // Expected coefficients: c=5, b=-3, a=2
        assertEquals(5,  p[0], 1.0/1000);
        assertEquals(-3, p[1], 1.0/1000);
        assertEquals(2,  p[2], 1.0/1000);

        // Check that f(p, x) matches original function
        for (int x = 0; x <= 2; x++) {
            double expected = 2*x*x - 3*x + 5;
            assertEquals(expected, Ex1.f(p, x), 1.0/1000);
        }
    }
    //============================================================
    /* testing functions for each "is equal" scenario */
    @Test
    void testEqualSimpleLinear() {
        double[] p1 = {2, 3};      // 3x + 2
        double[] p2 = {2, 3};      // 3x + 2
        assertTrue(Ex1.equals(p1, p2));
    }

    @Test
    void testEqualDifferentLengthsButSamePolynomial() {
        double[] p1 = {2, 3};        // 3x + 2
        double[] p2 = {2, 3, 0, 0};  // 3x + 2 (same function)
        assertTrue(Ex1.equals(p1, p2));
    }

    @Test
    void testNotEqualDifferentCoefficients() {
        double[] p1 = {2, 3};      // 3x + 2
        double[] p2 = {2, 4};      // 4x + 2
        assertFalse(Ex1.equals(p1, p2));
    }

    @Test
    void testEqualQuadratic() {
        double[] p1 = {5, -3, 2};       // 2x^2 - 3x + 5
        double[] p2 = {5, -3, 2};       // same
        assertTrue(Ex1.equals(p1, p2));
    }

    @Test
    void testNotEqualQuadratic() {
        double[] p1 = {5, -3, 2};       // 2x^2 - 3x + 5
        double[] p2 = {5, -3, 2.001};   // slightly different a
        assertFalse(Ex1.equals(p1, p2));
    }

    @Test
    void testEqualWithinEps() {
        double[] p1 = {2, 3};         // 3x + 2
        double[] p2 = {2.0005, 3.0004}; // differences < EPS = 0.001
        assertTrue(Ex1.equals(p1, p2));
    }

    @Test
    void testNotEqualBeyondEps() {
        double[] p1 = {2, 3};
        double[] p2 = {2.002, 3}; // difference in c0 = 0.002 > EPS
        assertFalse(Ex1.equals(p1, p2));
    }

    @Test
    void testDifferentDegreesButNotSameFunction() {
        double[] p1 = {1};         // 1
        double[] p2 = {1, 0, 0, 1}; // x^3 + 1
        assertFalse(Ex1.equals(p1, p2));
    }

    @Test
    void testNullCases() {
        assertFalse(Ex1.equals(null, new double[]{1}));
        assertFalse(Ex1.equals(new double[]{1}, null));
        assertFalse(Ex1.equals(null, null));
    }
    // Testers for function sameValue
    @Test
    void testSameValueSimpleLinear() {
        // p1(x) = x
        double[] p1 = {0, 1};

        // p2(x) = 3 - x
        double[] p2 = {3, -1};

        // They cross at x = 1.5
        double x = Ex1.sameValue(p1, p2, 0, 3, 1.0/1000);

        assertEquals(1.5, x, 1.0/1000);
    }

    @Test
    void testSameValueQuadraticAndLine() {
        // p1(x) = x^2
        double[] p1 = {0, 0, 1};

        // p2(x) = x
        double[] p2 = {0, 1};

        // Solve x^2 = 2x + 1 => (x-1)^2 = 0 => intersection at x = 1
        double x = Ex1.sameValue(p1, p2, 0, 3, 1.0/1000);

        assertEquals(0.0, x, 1.0/1000);
    }

    @Test
    void testSameValueIntersectionAtLeftSide() {
        // p1(x) = 2x
        double[] p1 = {0, 2};

        // p2(x) = x
        double[] p2 = {0, 1};

        // They meet at x = 0
        double x = Ex1.sameValue(p1, p2, 0, 2, 1.0/1000);

        assertEquals(0.0, x, 1.0/1000);
    }

    @Test
    void testSameValueAlternateFunctions() {
        double[] p1 = {-3, 1}; // x-3
        double[] p2 = {0};     // 0

        double x = Ex1.sameValue(p1, p2, 1, 5, 1.0/1000);

        assertEquals(3.0, x, 1.0/1000);
    }
    // ==================================
    // test area()
    @Test
    void testAreaSimpleTriangle() {
        double[] p1 = {0, 1};  // x
        double[] p2 = {0};     // 0

        double area = Ex1.area(p1, p2, 0, 1, 1000);

        assertEquals(0.5, area, 0.01);
    }
    @Test
    void testAreaConstantDifference() {
        double[] p1 = {5};
        double[] p2 = {2};

        double area = Ex1.area(p1, p2, 0, 4, 100);

        assertEquals(12.0, area, 0.01);
    }
    @Test
    void testAreaQuadraticVsLine() {
        double[] p1 = {0,0,1}; // x^2
        double[] p2 = {0,1};   // x

        double area = Ex1.area(p1, p2, 0, 1, 2000);

        assertEquals(1.0/6.0, area, 0.01);
    }

    // ===============================
    // test add()
    @Test
    void testAddSimple() {
        // p1(x) = 2x + 1
        double[] p1 = {1, 2};
        // p2(x) = 4x + 3
        double[] p2 = {3, 4};

        double[] expected = {4, 6};
        assertArrayEquals(expected, Ex1.add(p1, p2));
    }

    @Test
    void testAddDifferentLengths() {
        // p1(x) = 3x^2 + 2x + 1
        double[] p1 = {1, 2, 3};
        // p2(x) = 6x + 5
        double[] p2 = {5, 6};

        double[] expected = {6, 8, 3};
        assertArrayEquals(expected, Ex1.add(p1, p2));
    }

    @Test
    void testAddZeroPolynomial() {
        double[] p1 = {2, -3, 5};
        double[] p2 = {0}; // ZERO polynomial

        double[] expected = {2, -3, 5};
        assertArrayEquals(expected, Ex1.add(p1, p2));
    }

    @Test
    void testAddOpposites() {
        // p1(x) = x^2 - 2x + 5
        double[] p1 = {5, -2, 1};
        // p2(x) = -x^2 + 2x - 5
        double[] p2 = {-5, 2, -1};

        // Should return ZERO polynomial {0}
        double[] expected = {0};
        assertArrayEquals(expected, Ex1.add(p1, p2));
    }

    @Test
    void testAddCreatesTrailingZerosButTrimsThem() {
        double[] p1 = {1, 0, 0};  // actual value = 1
        double[] p2 = {2};        // actual value = 2

        double[] expected = {3};
        assertArrayEquals(expected, Ex1.add(p1, p2));
    }

    @Test
    void testAddWithNegativeCoefficients() {
        double[] p1 = {3, -4, 2};    // 2x^2 - 4x + 3
        double[] p2 = {-1, 5, -2};   // -2x^2 + 5x - 1

        double[] expected = {2, 1, 0}; // 0 is highest-degree term
        double[] trimmedExpected = {2, 1}; // Must trim zero

        assertArrayEquals(trimmedExpected, Ex1.add(p1, p2));
    }
    //===================================================
    // test mul()
    @Test
    void testMulSimple() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        double[] p1 = {1, 2};
        double[] p2 = {3, 4};

        double[] expected = {3, 10, 8};
        assertArrayEquals(expected, Ex1.mul(p1, p2));
    }

    @Test
    void testMulDifferentLengths() {
        // (1 + 2x + 3x^2) * (5 + 6x)
        // = 5 + 16x + 27x^2 + 18x^3
        double[] p1 = {1, 2, 3};
        double[] p2 = {5, 6};

        double[] expected = {5, 16, 27, 18};
        assertArrayEquals(expected, Ex1.mul(p1, p2));
    }
    @Test
    void testMulByOne() {
        // Multiply by 1 (identity polynomial)
        double[] p1 = {5, -2, 7};
        double[] p2 = {1};

        double[] expected = {5, -2, 7};
        assertArrayEquals(expected, Ex1.mul(p1, p2));
    }

    @Test
    void testMulTrailingZeros() {
        double[] p1 = {1, 0, 0};
        double[] p2 = {2};

        double[] expected = {2};
        assertArrayEquals(expected, Ex1.mul(p1, p2));
    }

    //==========================================
    //derivative
    @Test
    void testDerivativeSimple() {
        // f(x) = 3x^2 + 2x + 1
        double[] p = {1, 2, 3};

        // f'(x) = 6x + 2  → {2, 6}
        double[] expected = {2, 6};

        assertArrayEquals(expected, Ex1.derivative(p));
    }
    @Test
    void testDerivativeConstant() {
        double[] p = {5}; // f(x) = 5
        assertArrayEquals(Ex1.ZERO, Ex1.derivative(p));
    }
    @Test
    void testDerivativeZero() {
        assertArrayEquals(Ex1.ZERO, Ex1.derivative(Ex1.ZERO));
    }
    @Test
    void testDerivativeHighDegree() {
        // f(x) = 4x^4 - 3x^3 + 2x + 7
        double[] p = {7, 2, 0, -3, 4};

        // f'(x) = 16x^3 - 9x^2 + 2  → {2, 0, -9, 16}
        double[] expected = {2, 0, -9, 16};

        assertArrayEquals(expected, Ex1.derivative(p));
    }
    @Test
    void testDerivativeTrailingZeros() {
        // f(x) = 0x^3 + 0x^2 + 5x + 0
        double[] p = {0, 5, 0, 0};

        // f'(x) = 5 → {5}
        assertArrayEquals(new double[]{5}, Ex1.derivative(p));
    }
}



