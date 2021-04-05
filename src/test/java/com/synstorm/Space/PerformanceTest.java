package com.synstorm.Space;

import com.google.common.util.concurrent.AtomicDoubleArray;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class PerformanceTest {

    @Test
    public void measureTime() {

        // test params for grid search
        int[] signalSourceOpts = { 1, 10, 100, 1000, 10000, 100000 };
        int[] signalsOpts = { 1, 10/*, 100*/ };
        int[] sizeOpts = { 10, 25, 50/*, 100, 200*/ };
        int steps = 100;
        int trialsPerConfig = 10;

        System.out.println(
                "\"Signal sources\", "+
                "\"Signals\", "+
                "\"Space size\", "+
                "\"Steps\", "+
                "\"Time\""
        );
        for (int signalSourceN: signalSourceOpts) {
            for (int signalsN: signalsOpts) {
                for (int size: sizeOpts) {
                    for (int trial=0; trial<trialsPerConfig; trial++) {
                        SpaceGrid s = getSpace(signalSourceN, signalsN, size);
                        long start = System.nanoTime();
                        runSimulation(steps, signalsN, s);
                        long stop = System.nanoTime();

                        System.out.println(
                                Integer.toString(signalSourceN) + ", " +
                                        Integer.toString(signalsN) + ", " +
                                        Integer.toString(size) + ", " +
                                        Integer.toString(steps) + ", " +
                                        Long.toString(stop - start)
                        );
                    }
                }
            }
        }
    }

    @Test
    public void measureMemory() {

        // test params for grid search
        int[] signalSourceOpts = { /*1, 10, 100, 1000, 10000, */100000 };
        int[] signalsOpts = { /*1, */10/*, 100*/ };
        int[] sizeOpts = { /*10, 25, */50/*, 100, 200*/ };
        int steps = 1000;

        System.out.println(
                "\"Signal sources\", "+
                "\"Signals\", "+
                "\"Space size\", "+
                "\"Steps\", "+
                "\"Memory\""
        );
        Runtime runtime = Runtime.getRuntime();
        for (int signalSourceN: signalSourceOpts) {
            for (int signalsN: signalsOpts) {
                for (int size: sizeOpts) {
                    SpaceGrid s = getSpace(signalSourceN, signalsN, size);

                    for (int i=0; i<steps; i++) {
                        for (int sig = 0; sig < signalsN; sig++) {
                            AtomicDoubleArray diff = s.calcDiffusion(sig);
                            s.updateDiffusion(sig, diff);
                        }
                        System.out.println(
                                Integer.toString(signalSourceN) + ", " +
                                Integer.toString(signalsN) + ", " +
                                Integer.toString(size) + ", " +
                                Integer.toString(steps) + ", " +
                                Long.toString(runtime.totalMemory() - runtime.freeMemory())
                        );

                    }


                }
            }
        }
    }

    private void runSimulation(int steps, int signalsN, SpaceGrid s) {
        for (int i=0; i<steps; i++) {
            for (int sig=0; sig<signalsN; sig++) {
                AtomicDoubleArray diff = s.calcDiffusion(sig);
                s.updateDiffusion(sig, diff);
            }
        }
    }

    @NotNull
    private SpaceGrid getSpace(int signalSourceN, int signalsN, int size) {
        double[] diffRates = new double[signalsN];
        Arrays.fill(diffRates, 0.1);
        SpaceGrid s = new SpaceGrid(size, diffRates);
        int signalId = 0;
        int signalStep = Math.floorDiv(signalSourceN, signalsN);

        for (int signalSource=0; signalSource<Math.min(signalSourceN, size*size*size); signalSource++) {
            if (signalSource>0 && Math.floorMod(signalSource, signalStep) == 0) {
//                s.updateDiffusion(signalId);
                signalId++;
            }
            s.addObject(signalSource, signalSource, 1.);
            s.updateConcentration(signalSource, signalId, 100);
        }
        return s;
    }

}