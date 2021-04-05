package com.synstorm.Space;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Longs;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SpaceGridTest {
    private SpaceGrid sp = new SpaceGrid(100, new double[]{0.1,0.1});

    @Test
    public void updateConcentrations() {
    }

    private byte[][] createObjectsForPaths (int[][] paths, SpaceGrid sp) {
        final byte[][] tracks = new byte[paths.length][];
        int objCounter = 0;
        for (int i = 0; i < paths.length; i++) {
            byte[] tracksOut = new byte[1];
            byte[] objectsOut = new byte[1];
            for (int j = 0; j < paths[i].length; j++) {
                ++tracksOut[0];
                ++objectsOut[0];
                final int coordinate = paths[i][j];
                final long id = sp.addObject(objCounter++, coordinate, 1.);

                tracksOut = Bytes.concat(tracksOut, Ints.toByteArray(coordinate));
                objectsOut = Bytes.concat(objectsOut, new byte[]{1});
                objectsOut = Bytes.concat(objectsOut, Longs.toByteArray(id));
            }

            tracks[i] = Bytes.concat(tracksOut, objectsOut);
        }

        return tracks;
    }

    private void doPathTest(int[][] paths, SpaceGrid sp, int expectedAllowed) {
        byte[][] tracks = createObjectsForPaths(paths, sp);
        int[] allowedTracks = sp.resolveCollisions(tracks);
        assertEquals(expectedAllowed, allowedTracks.length);
    }

    // flat coordinate is calculated as 'size * (x * size + y) + z'
    // so for space with edge of 100 coordinate (1,0,0) will be 10000 and
    // coordinate (35, 35, 35) will be 353535 and so on

    // tracks should be straight without any changes of vector
    // they defines by start coordinate, vector and length
    // so if we have start in coordinate (1,0,0), vector (0,0,1) and length 4
    // we will get track {10000, 10001, 10002, 10003}
    int A=10000, B=10001, C=10002, D=10003, X=0, Y=1, Z=2;
    int[] trackABCD = new int[]{ A, B, C, D };
    int[] trackA = new int[]{ A };
    int[] trackD = new int[]{ D };
    int[] trackBC = new int[]{ 1,1,0, 1,2,0, 1,1,2, 1,0,3 };

    @Test
    public void oneTrackNoCollisions () {
        doPathTest(new int[][] { {A} }, sp, 1);
        doPathTest(new int[][] { {A, B} }, sp, 1);
    }

    @Test
    public void emptyTracks () {
        doPathTest(new int[][] { {} }, sp, 1);
        doPathTest(new int[][] { {A}, {} }, sp, 1);
        doPathTest(new int[][] { {A}, {}, {} }, sp, 1);
        doPathTest(new int[][] { {}, {A}, {} }, sp, 1);
    }

    @Test
    public void oneTrackOneCollision () {
        doPathTest(new int[][] { {A, B, A} }, sp, 1);
    }

    @Test
    public void twoTracksNoCollisions () {
        doPathTest(new int[][] { {A, B, C}, {D} }, sp, 2);
        doPathTest(new int[][] { {A, B}, {C, D} }, sp, 2);
    }

    @Test
    public void twoTracksOneCollision () {
        // Single collision
        doPathTest(new int[][] { {A, B, C, D}, {A} }, sp, 1);
        doPathTest(new int[][] { {A, B, C, D}, {D} }, sp, 1);
        doPathTest(new int[][] { {A, B, C}, {C, D} }, sp, 1);

        // Multiple collisions
        doPathTest(new int[][] { {A, B, C, D}, {A, D} }, sp, 1);

        // Subset collision
        doPathTest(new int[][] { {A, B, C, D}, {B, C} }, sp, 1);
        doPathTest(new int[][] { {B, C}, {A, B, C, D} }, sp, 1);

        // Identical (reverse)
        doPathTest(new int[][] { {A, B, C}, {A, B, C} }, sp, 1);
        doPathTest(new int[][] { {A, B, C}, {C, B, A} }, sp, 1);
    }

    @Test
    public void threeTracksOneCollision () {
        // Single collision
        doPathTest(new int[][] { {A, B, C, D}, {A}, {X,Y,Z} }, sp, 2);
        doPathTest(new int[][] { {A, B, C, D}, {D}, {X,Y,Z} }, sp, 2);
        doPathTest(new int[][] { {A, B, C}, {C, D}, {X,Y,Z} }, sp, 2);

        // Multiple collisions
        doPathTest(new int[][] { {A, B, C, D}, {A, D}, {X,Y,Z} }, sp, 2);

        // Subset collision
        doPathTest(new int[][] { {A, B, C, D}, {B, C}, {X,Y,Z} }, sp, 2);
        doPathTest(new int[][] { {B, C}, {A, B, C, D}, {X,Y,Z} }, sp, 2);

        // Identical (reverse)
        doPathTest(new int[][] { {A, B, C}, {A, B, C}, {X,Y,Z} }, sp, 2);
        doPathTest(new int[][] { {A, B, C}, {C, B, A}, {X,Y,Z} }, sp, 2);
    }

    @Test
    public void threeTracksTwoCollision () {
        // Single collision
        doPathTest(new int[][] { {A, B, C, D}, {A}, {A,X,Y,Z} }, sp, 1);
        doPathTest(new int[][] { {A, B, C, D}, {D}, {A,X,Y,Z} }, sp, 2);
        doPathTest(new int[][] { {A, B, C}, {C, D}, {B,X,Y,Z} }, sp, 2);

        // Multiple collisions
        doPathTest(new int[][] { {A, B, C, D}, {A, D}, {B,C} }, sp, 2);

        // Subset collision
        doPathTest(new int[][] { {A, B, C, D}, {B, C}, {A,D} }, sp, 2);
        doPathTest(new int[][] { {B, C}, {A, B, C, D}, {X,A,B,C,D,Z} }, sp, 1);

        // Identical (reverse)
        doPathTest(new int[][] { {A, B, C}, {A, B, C}, {A,B,C} }, sp, 1);
        doPathTest(new int[][] { {A, B, C}, {C, B, A}, {C,B,A} }, sp, 1);
    }

    @Test
    public void closed () {
        doPathTest(new int[][] { {A, B}, {B, C}, {C,A} }, sp, 1);
        doPathTest(new int[][] { {A, B}, {B, C}, {C,D}, {D,A} }, sp, 2);
        doPathTest(new int[][] { {A, B}, {B, C}, {C,D}, {D,X}, {X,A} }, sp, 2);
    }
}