package com.synstorm.Space;

import com.google.common.util.concurrent.AtomicDoubleArray;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class SignalsTest {

    final double delta = 1e-10;


/**
 * Excluding cubes with edge < 3 because of 3d->1d conversion algorithm
 */
//    /**
//     * Cube 1x1x1
//     * Source in the only one cell
//     * All diffusion is leaked
//     */
//    @Test
//    public void centerPositionOneSignalSmallCube() {
//        double initialSignal0 = 100; double signalDiffRate0 = 0.1; int signalSource0 = coord(0,0,0,1);
//        int dim = 1;
//        Space s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});
//
//        for (int i = 0; i < 100; i++) {
//            double expectedLeak = chainsum(initialSignal0, signalDiffRate0, i);
//            assertSignal(s, 0,
//                    initialSignal0,
//                    expectedLeak,
//                    Collections.singletonMap(signalSource0, initialSignal0 - expectedLeak));
//            s.calcDiffusion(0);
//            s.updateDiffusion(0);
//        }
//    }
//
//    /**
//     * Cube 2x2x2
//     * Source at {0, 0, 0}
//     * All diffusion is leaked after 1 step
//     */
//    @Test
//    public void cube2() {
//        int dim=2;
//        double initialSignal0 = 100; double signalDiffRate0 = 0.1; int signalSource0 = coord(0,0,0,dim);
//        Space s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});
//
//        assertSignal(s, 0, initialSignal0, 0,
//                new double[]{initialSignal0, 0, 0, 0, 0, 0, 0, 0});
//
//        s.calcDiffusion(0);
//        s.updateDiffusion(0);
//
//        double ax = initialSignal0*signalDiffRate0/26;
//        assertSignal(s, 0, initialSignal0, 19*ax,
//                new double[]{initialSignal0*(1-signalDiffRate0), ax, ax, ax, ax, ax, ax, ax});
//    }


    /**
     * SECTION CUBE SIDE 3
     * SINGLE STEP
     */

    /**
     * Cube 3x3x3
     * Source in the center {1, 1, 1}
     */
    @Test
    public void cube3center() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(1, 1, 1, 3);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, x, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = initialSignal0 * signalDiffRate0 / 26;
        x *= 1 - signalDiffRate0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                d, d, d,  d, d, d,  d, d, d,
                d, d, d,  d, x, d,  d, d, d,
                d, d, d,  d, d, d,  d, d, d});
    }

    /**
     * Cube 3x3x3
     * Source in the corner {0, 0, 0}
     */
    @Test
    public void cube3corner() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(0, 0, 0, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                x, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = initialSignal0 * signalDiffRate0 / 26;
        x *= 1-signalDiffRate0;
        assertSignal(s, 0, initialSignal0, 19*d, new double[]{
                x, d, 0,  d, d, 0,  0, 0, 0,
                d, d, 0,  d, d, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});
    }

    /**
     * Cube 3x3x3
     * Source in the corner {2, 2, 2}
     */
    @Test
    public void cube3corner2() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(2, 2, 2, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, x});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = initialSignal0 * signalDiffRate0 / 26;
        x *= 1-signalDiffRate0;
        assertSignal(s, 0, initialSignal0, 19*d, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, d, d,  0, d, d,
                0, 0, 0,  0, d, d,  0, d, x});
    }

    /**
     * Cube 3x3x3
     * Source in the side center {1, 1, 0}
     */
    @Test
    public void cube3side() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(1, 1, 0, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0,  0, x, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = initialSignal0 * signalDiffRate0 / 26;
        x *= 1-signalDiffRate0;
        assertSignal(s, 0, initialSignal0, 9*d, new double[]{
                d, d, d,  d, x, d,  d, d, d,
                d, d, d,  d, d, d,  d, d, d,
                0, 0, 0,  0, 0, 0,  0, 0, 0});
    }

    /**
     * Cube 3x3x3
     * Source at the edge {1, 0, 0}
     */
    @Test
    public void cube3edge() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(1, 0, 0, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, x, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = initialSignal0 * signalDiffRate0 / 26;
        x *= 1-signalDiffRate0;
        assertSignal(s, 0, initialSignal0, 15*d, new double[]{
                d, x, d,  d, d, d,  0, 0, 0,
                d, d, d,  d, d, d,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});
    }

    /**
     * SECTION CUBE SIDE 3
     * Multiple steps
     */

    @Test
    public void cube3centerLong() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(1, 1, 1, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, x, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = x * signalDiffRate0 / 26;
        x *= (1-signalDiffRate0);
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                d, d, d,  d, d, d,  d, d, d,
                d, d, d,  d, x, d,  d, d, d,
                d, d, d,  d, d, d,  d, d, d});

        diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double l = s.getLeaked(0);
        double dd = (d*signalDiffRate0)/26;
        assertEquals(l, dd * (8 * 19 + 12 * 15 + 6 * 9), delta); // here we only count number of cells outside the cube
    }

    /**
     * SECTION CUBE SIDE 3
     * Multiple steps
     */

    @Test
    public void cube3centerLongDiffRate03() {
        int dim = 3;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.3;
        int signalSource0 = coord(1, 1, 1, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, x, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = x * signalDiffRate0 / 26;
        x *= (1-signalDiffRate0);
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                d, d, d,  d, d, d,  d, d, d,
                d, d, d,  d, x, d,  d, d, d,
                d, d, d,  d, d, d,  d, d, d});

        diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double l = s.getLeaked(0);
        double dd = (d*signalDiffRate0)/26;
        assertEquals(l, dd * (8 * 19 + 12 * 15 + 6 * 9), delta); // here we only count number of cells outside the cube
    }

    /**
     * SECTION CUBE SIDE 4
     * Multiple steps
     */

    @Test
    public void cube4centerLong() {
        int dim = 4;
        double initialSignal0 = 100;
        double signalDiffRate0 = 0.1;
        int signalSource0 = coord(1, 1, 1, dim);
        SpaceGrid s = init(dim, new double[][] {{initialSignal0}}, new int[][] {{signalSource0}}, new double[] {signalDiffRate0});

        double x = initialSignal0;
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                0, 0, 0, 0,  0, x, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
        });

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double d = x * signalDiffRate0 / 26;
        x *= (1-signalDiffRate0);
        assertSignal(s, 0, initialSignal0, 0, new double[]{
                d, d, d, 0,  d, d, d, 0,  d, d, d, 0,  0, 0, 0, 0,
                d, d, d, 0,  d, x, d, 0,  d, d, d, 0,  0, 0, 0, 0,
                d, d, d, 0,  d, d, d, 0,  d, d, d, 0,  0, 0, 0, 0,
                0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
        });

        diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double l = s.getLeaked(0);
        double dd = (d*signalDiffRate0)/26;
        assertEquals(l, dd * (1 * 19 + 6 * 15 + 12 * 9), delta); // here we only count number of cells outside the cube
    }

    /**
     * Cube 3x3x3
     * Sources are in the corners {0, 0, 0}, {2, 2, 0} {2, 2, 2}
     */
    @Test
    public void sources2() {
        int dim = 3;
        double xi = 101;
        double xd = 0.1;
        int signalSource0 = coord(0, 0, 0, dim);
        double yi = 33;
        double yd = 0.8;
        int signalSource1 = coord(1, 0, 0, dim);

        SpaceGrid s = init(dim, new double[][] {{xi}, {yi}}, new int[][] {{signalSource0}, {signalSource1}}, new double[] {xd, yd});

        assertSignal(s, 0, xi, 0, new double[]{
                xi, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        assertSignal(s, 1, yi, 0, new double[]{
                0, yi, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});
    }

    /**
     * Cube 3x3x3
     * three signal sources
     */
    @Test
    public void sources3() {
        int dim = 3;
        double xi = 100; double xd = 0.1; int signalSource0 = coord(0, 0, 0, dim);
        double yi = 100; double yd = 0.1; int signalSource1 = coord(1, 0, 0, dim);
        double zi = 100; double zd = 0.1; int signalSource2 = coord(1, 1, 0, dim);

        SpaceGrid s = init(dim, new double[][] {{xi}, {yi}, {zi}},
                new int[][] {{signalSource0}, {signalSource1}, {signalSource2}}, new double[] {xd, yd, zd});

        assertSignal(s, 0, xi, 0, new double[]{
                xi, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        assertSignal(s, 1, yi, 0, new double[]{
                0, yi, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        assertSignal(s, 2, zi, 0, new double[]{
                0, 0, 0,  0, zi, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});
    }


/**
     * Cube 3x3x3
     * Three different signal sources with leak
     */

    @Test
    public void cube3corner3Signals3Sources() {
        int dim = 3;
        double xi = 101;
        double xd = 0.1;
        int signalSource0 = coord(0, 0, 0, dim);
        double yi = 33;
        double yd = 0.8;
        int signalSource1 = coord(2, 2, 0, dim);
        double zi = 666;
        double zd = 0.5;
        int signalSource2 = coord(2, 2, 2, dim);

        SpaceGrid s = init(dim, new double[][] {{xi}, {yi}, {zi}},
                new int[][] {{signalSource0}, {signalSource1}, {signalSource2}},
                new double[] {xd, yd, zd});

        assertSignal(s, 0, xi, 0, new double[]{
                xi, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        assertSignal(s, 1, yi, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, yi,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        assertSignal(s, 2, zi, 0, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, zi});

        AtomicDoubleArray diff0 = s.calcDiffusion(0);
        s.updateDiffusion(0, diff0);
        AtomicDoubleArray diff1 = s.calcDiffusion(1);
        s.updateDiffusion(1, diff1);
        AtomicDoubleArray diff2 = s.calcDiffusion(2);
        s.updateDiffusion(2, diff2);

        double xj = xi * (1-xd); double dx = xi * xd/26;
        double yj = yi * (1-yd); double dy = yi * yd/26;
        double zj = zi * (1-zd); double dz = zi * zd/26;

        assertSignal(s, 0, xi, 19*dx, new double[]{
                xj, dx, 0,  dx, dx, 0,  0, 0, 0,
                dx, dx, 0,  dx, dx, 0,  0, 0, 0,
                0, 0, 0,    0, 0, 0,    0, 0, 0});

        assertSignal(s, 1, yi, 19*dy, new double[]{
                0, 0, 0,  0, dy, dy,  0, dy, yj,
                0, 0, 0,  0, dy, dy,  0, dy, dy,
                0, 0, 0,  0, 0, 0,    0, 0, 0});

        assertSignal(s, 2, zi, 19*dz, new double[]{
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, dz, dz,  0, dz, dz,
                0, 0, 0,  0, dz, dz,  0, dz, zj});
    }

    /**
     * Cube 3x3x3
     * Two sources of same signal
     */
    @Test
    public void sourcesSameSignal() {
        int dim = 3;
        double xi = 101;
        double dr = 0.1;
        int signalSource0 = coord(0, 0, 0, dim);
        double yi = 33;
        int signalSource1 = coord(1, 0, 0, dim);

        SpaceGrid s = init(dim, new double[][] {{xi, yi}}, new int[][] {{signalSource0, signalSource1}}, new double[] {dr});

        assertSignal(s, 0, xi+yi, 0, new double[]{
                xi, yi, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double xj = xi * (1-dr); double dx = xi * dr/26;
        double yj = yi * (1-dr); double dy = yi * dr/26;

        assertSignal(s, 0, xi+yi, 19*dx+15*dy, new double[]{
                xj+dy, yj+dx, dy,  dy+dx, dy+dx, dy,  0, 0, 0,
                dy+dx, dy+dx, dy,  dy+dx, dy+dx, dy,  0, 0, 0,
                0,     0,     0,   0,     0,     0,   0, 0, 0});

    }

    /**
     * Cube 3x3x3
     * Two sources of same signal in the same coord
     */
    @Test
    public void sourcesSameSignalSameCoord() {
        int dim = 3;
        double xi = 101;
        double dr = 0.1;
        int signalSource0 = coord(0, 0, 0, dim);
        double yi = 33;
        int signalSource1 = coord(0, 0, 0, dim);

        SpaceGrid s = init(dim, new double[][] {{xi, yi}}, new int[][] {{signalSource0, signalSource1}}, new double[] {dr});

        assertSignal(s, 0, xi+yi, 0, new double[]{
                xi+yi, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

        AtomicDoubleArray diff = s.calcDiffusion(0);
        s.updateDiffusion(0, diff);

        double z = (xi+yi) * (1-dr); double d = (xi+yi) * dr/26;

        assertSignal(s, 0, xi+yi, 19*d, new double[]{
                z, d, 0,  d, d, 0,  0, 0, 0,
                d, d, 0,  d, d, 0,  0, 0, 0,
                0, 0, 0,  0, 0, 0,  0, 0, 0});

    }


    private double totalConcentration(SpaceGrid s, final int signal) {
        double result = 0;
        int sizeCube = s.getSize() * s.getSize() * s.getSize();
        for (int i = 0; i < sizeCube; i++) {
            result += s.getConcentration(signal, i);
        }
        return result;
    }

    private void assertSignal(SpaceGrid s,
                              int signal,
                              double initAmount,
                              double expectedLeak,
                              double[] expectedPoints)
    {
        Map<Integer,Double> expectedCoords = new HashMap<>();
        for (int i = 0; i < expectedPoints.length; i++) {
            expectedCoords.put(i, expectedPoints[i]);
        }
        assertSignal(s, signal, initAmount, expectedLeak, expectedCoords);
    }

    private void assertSignal(SpaceGrid s,
                              int signal,
                              double initAmount,
                              double expectedLeak,
                              Map<Integer, Double> expectedPoints)
    {
        assertEquals(initAmount, totalConcentration(s, signal) + s.getLeaked(signal), delta);
        for (Integer coord: expectedPoints.keySet()) {
            assertEquals(expectedPoints.get(coord), s.getConcentration(signal, coord), delta);
        }
        assertEquals(expectedLeak, s.getLeaked(signal), delta);
    }

    private double chainsum(double x, double a, int n) {
        double res = 0;
        double rem = x;
        for (int i = 1; i <= n; i++) {
            double leaked = a * rem;
            res += leaked;
            rem -= leaked;
        }
        return res;
    }

    private int coord(int x, int y, int z, int d) {
        return x + y*d + z*d*d;
    }

    private SpaceGrid init(int dim, double[][] sourceConc, int[][] sourceCoord, double[] diffRate) {
        SpaceGrid s = new SpaceGrid(dim, diffRate);
        int objid=0;
        for (int i = 0; i < sourceConc.length; i++) {
            for (int j = 0; j < sourceConc[i].length; j++) {
                s.addObject(objid, sourceCoord[i][j], 1.);
                s.updateConcentration(objid, i, sourceConc[i][j]);
                objid++;
            }
//            s.updateDiffusion(i);
        }
        return s;
    }

    private double[] concArray(SpaceGrid s, int signal) {
        int coords = s.getSize()*s.getSize()*s.getSize();
        double[] res = new double[coords];
        for (int i = 0; i < coords; i++) {
            res[i] = s.getConcentration(0, i);
        }
        return res;
    }
}