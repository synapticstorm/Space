package com.synstorm.Space;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Longs;
import com.google.common.util.concurrent.AtomicDouble;
import com.google.common.util.concurrent.AtomicDoubleArray;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

import java.util.Comparator;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

public class SpaceGrid {
    //region Fields
    private final int size;
    private final int sizeSquare;
    private final int sizeCube;
    private final CellState[] cube;
    private final double[] diffRates;
    private final double[] leaked;
    private final AtomicDoubleArray[] concentrations;

    // Map for border cells restrictions
    private final TIntObjectHashMap<TIntHashSet> borders;

    // Flat vectors indices of neighborhood per each cell
    private final int[] nbh = new int[26];

    // Flat vectors indices of full neighborhood including coordinate of interest per each cell
    private final int[] fnbh = new int[27];

    // Map for <object, coordinate>
    private final ConcurrentHashMap<Long, Integer> objIdCoordinate;

    // Map for <object, random>
    private final ConcurrentHashMap<Long, Random> objRandoms;

    private final int intSize = Integer.BYTES;
    private final int longSize = Long.BYTES;
    private final int doubleSize = Double.BYTES;
    //endregion

    //region Inner classes
    private static class CellState {
        private final ConcurrentHashMap<Long, AtomicDouble> objects;
        private final AtomicDouble occupiedVolume;

        public CellState() {
            objects = new ConcurrentHashMap<>();
            occupiedVolume = new AtomicDouble(0.);
        }

        public void putObject(final long obj, final double volume) {
            objects.put(obj, new AtomicDouble(volume));
            occupiedVolume.addAndGet(volume);
        }

        public double removeObject(final long obj) {
            final double volume = objects.remove(obj).get();
            occupiedVolume.addAndGet(-volume);
            return volume;
        }

        public long[] getMovingCandidates(final double volume) {
            final double[] cumulativeVolume = new double[]{0.};
            return objects.entrySet().stream()
                    .sorted(Comparator.comparingDouble(e -> e.getValue().get()))
                    .filter(e -> (cumulativeVolume[0] += e.getValue().get()) <= volume)
                    .mapToLong(Map.Entry::getKey)
                    .toArray();
        }

        public double getSumVolume(final long[] objIds) {
            return LongStream.of(objIds)
                    .mapToDouble(i -> objects.get(i).get())
                    .sum();
        }
    }
    //endregion

    //region Constructors
    public SpaceGrid(int spaceSize, @NotNull double[] signalsDiffRate) {
        // Vectors of neighborhood per each cell
        final int[][] nbhD = {
                {-1,-1,-1}, {-1,-1,0}, {-1,-1,1},   {-1,0,-1}, {-1,0,0}, {-1,0,1},   {-1,1,-1}, {-1,1,0}, {-1,1,1},
                { 0,-1,-1}, { 0,-1,0}, { 0,-1,1},   { 0,0,-1},           { 0,0,1},   { 0,1,-1}, { 0,1,0}, { 0,1,1},
                { 1,-1,-1}, { 1,-1,0}, { 1,-1,1},   { 1,0,-1}, { 1,0,0}, { 1,0,1},   { 1,1,-1}, { 1,1,0}, { 1,1,1} };

        // Vectors of full neighborhood including coordinate of interest per each cell
        final int[][] fnbhD = {
                {-1,-1,-1}, {-1,-1,0}, {-1,-1,1},   {-1,0,-1}, {-1,0,0}, {-1,0,1},   {-1,1,-1}, {-1,1,0}, {-1,1,1},
                { 0,-1,-1}, { 0,-1,0}, { 0,-1,1},   { 0,0,-1}, { 0,0,0}, { 0,0,1},   { 0,1,-1}, { 0,1,0}, { 0,1,1},
                { 1,-1,-1}, { 1,-1,0}, { 1,-1,1},   { 1,0,-1}, { 1,0,0}, { 1,0,1},   { 1,1,-1}, { 1,1,0}, { 1,1,1} };

        final int signalCount = signalsDiffRate.length;
        size = spaceSize;
        sizeSquare = size * size;
        sizeCube = size * size * size;
        cube = new CellState[sizeCube];
        diffRates = new double[signalCount];
        leaked = new double[signalCount];
        concentrations = new AtomicDoubleArray[signalCount];
        System.arraycopy(signalsDiffRate, 0, diffRates, 0, signalsDiffRate.length);

        for (int i = 0; i < nbhD.length; i++)
            nbh[i] = getFlatCoordinate(nbhD[i][0], nbhD[i][1], nbhD[i][2]);

        for (int i = 0; i < fnbhD.length; i++)
            fnbh[i] = getFlatCoordinate(fnbhD[i][0], fnbhD[i][1], fnbhD[i][2]);

        // border for lower x
        final int[] x0 = new int[]{nbh[0], nbh[1], nbh[2], nbh[3], nbh[4], nbh[5], nbh[6], nbh[7], nbh[8]};
        // border for higher x
        final int[] xS = new int[]{nbh[17], nbh[18], nbh[19], nbh[20], nbh[21], nbh[22], nbh[23], nbh[24], nbh[25]};

        // border for lower y
        final int[] y0 = new int[]{nbh[0], nbh[1], nbh[2], nbh[9], nbh[10], nbh[11], nbh[17], nbh[18], nbh[19]};
        // border for higher y
        final int[] yS = new int[]{nbh[6], nbh[7], nbh[8], nbh[14], nbh[15], nbh[16], nbh[23], nbh[24], nbh[25]};

        // border for lower z
        final int[] z0 = new int[]{nbh[0], nbh[3], nbh[6], nbh[9], nbh[12], nbh[14], nbh[17], nbh[20], nbh[23]};
        // border for higher z
        final int[] zS = new int[]{nbh[2], nbh[5], nbh[8], nbh[11], nbh[13], nbh[16], nbh[19], nbh[22], nbh[25]};

        borders = new TIntObjectHashMap<>();
        objIdCoordinate = new ConcurrentHashMap<>();
        objRandoms = new ConcurrentHashMap<>();

        for (int i = 0; i < signalCount; ++i)
            concentrations[i] = new AtomicDoubleArray(sizeCube);

        for (int i = 0; i < sizeCube; ++i)
            cube[i] = new CellState();

        for (int x = 0; x < size; ++x)
            for (int y = 0; y < size; ++y)
                for (int z = 0; z < size; ++z) {
                    final int coordinate = getFlatCoordinate(x, y, z);
                    cube[coordinate] = new CellState();
                    if (x == 0 || x == size - 1 || y == 0 || y == size - 1 || z == 0 || z == size - 1) {
                        final TIntHashSet restricted = new TIntHashSet();
                        borders.putIfAbsent(coordinate, restricted);
                        if (x == 0) restricted.addAll(x0);
                        if (x == size - 1) restricted.addAll(xS);
                        if (y == 0) restricted.addAll(y0);
                        if (y == size - 1) restricted.addAll(yS);
                        if (z == 0) restricted.addAll(z0);
                        if (z == size - 1) restricted.addAll(zS);
                    }
                }
    }
    //endregion

    //region Getters and Setters
    public int getSize() {
        return size;
    }

    public Map<Long, int[]> getObjects() {
        return objIdCoordinate.entrySet()
                .parallelStream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> restoreCoordinateFromFlat(e.getValue())));
    }
    //endregion

    //region Public Methods
    /**
     * @param objId for object in space
     * @param coordinate for object
     * @param volume for volume: 0 - voluminous, 1 - volumetric
     * @return added status
     */
    public long addObject(final long objId, final int coordinate, final double volume) {
        final CellState cState = cube[coordinate];
        cState.putObject(objId, volume);

        objIdCoordinate.put(objId, coordinate);
        objRandoms.put(objId, new Random(objId));
        return objId;
    }

    /**
     * @param objId for object in space
     */
    public long removeObject(final long objId) {
        final int idx = objIdCoordinate.remove(objId);
        final CellState cState = cube[idx];
        cState.removeObject(objId);

        objRandoms.remove(objId);
        objIdCoordinate.remove(objId);
        return objId;
    }

    /**
     * @param objId id of moving object
     * @param dstCoordinate coordinate to move
     * @return id of moved object
     */
    public long moveObject(final long objId, final int dstCoordinate) {
        final int srcCoordinate = Objects.requireNonNull(objIdCoordinate.put(objId, dstCoordinate));
        final CellState srcState = cube[srcCoordinate];
        final CellState dstState = cube[dstCoordinate];
        final double volume = srcState.removeObject(objId);
        dstState.putObject(objId, volume);

        return objId;
    }

    /**
     * @param objId of interest
     * @return coordinate of given object id
     */
    public int getCoordinate(final long objId) {
        return objIdCoordinate.get(objId);
    }

    /**
     * @param signal of interest
     * @return total amount of leaked signal
     */
    public double getLeaked(final int signal) {
        return leaked[signal];
    }

    /**
     * @param signal of interest
     * @param coordinate of interest
     * @return concentration of signal in coordinate
     */
    public double getConcentration(final int signal, final int coordinate) {
        return concentrations[signal].get(coordinate);
    }

    /**
     * @param objId for object in space
     * @param signal index
     * @param val for spread/gather value
     */
    public void updateConcentration(final long objId, final int signal, final double val) {
        // here we update concentrations after spread/gather of signal
        final int coordinate = objIdCoordinate.get(objId);
        concentrations[signal].addAndGet(coordinate, val);
    }

    /**
     * @param signal index
     */
    public long[] updateDiffusion(final int signal, final @NotNull AtomicDoubleArray diff) {
        leaked[signal] += diff.get(sizeCube);

        return IntStream.range(0, sizeCube)
                .parallel()
                .filter(i -> diff.get(i) != 0.)
                .peek(i -> concentrations[signal].addAndGet(i, diff.get(i)))
                .mapToObj(i -> cube[i].objects.keySet().stream())
                .flatMapToLong(i -> i.mapToLong(Long::longValue))
                .toArray();
    }

    /**
     * @param signal of interest
     */
    public AtomicDoubleArray calcDiffusion(final int signal) {
        final AtomicDoubleArray diffResult = new AtomicDoubleArray(sizeCube + 1);
        final double diffusionRate = diffRates[signal];

        IntStream.range(0, sizeCube)
                .parallel()
                .forEach(i -> {
                    final double diffusedAmount = concentrations[signal].get(i) * diffusionRate;    // total diffused amount
                    if (diffusedAmount == 0) return;                                                // if nothing to diffuse - skip

                    final double dr = diffusedAmount / 26;                                          // diffusion per cell pair
                    diffResult.addAndGet(i, -diffusedAmount);                                       // update current cell
                    final TIntHashSet restricted = borders.get(i);
                    for (int dims : nbh) {
                        if (restricted != null && restricted.contains(dims)) {
                            diffResult.addAndGet(sizeCube, dr);                                     // tracking diffusion outside borders
                            continue;
                        }

                        final int di = i + dims;                                                    // flat cells
                        diffResult.addAndGet(di, dr);                                               // update neighbors cells
                    }
                });

        return diffResult;
    }

    public int getFirstCoordinate(@NotNull final byte[] bTrack) {
        return extractTrack(bTrack)[0];
    }

    /**
     * @param bTrack to shift
     * @return array of moved object id's
     */
    @NotNull
    public long[] shiftTrack(@NotNull final byte[] bTrack) {
        final int[] track = extractTrack(bTrack);
        final int objectsCount = getObjectCountInTrack(bTrack);
        final long[] movedIds = new long[objectsCount];
        byte[] buff = new byte[longSize];
        int position = track.length * intSize + 2;
        int objectsCounter = 0;
        for (int i = 0; i < track.length - 1; i++) {
            int objCountPos = bTrack[position++];
            for (int j = 0; j < objCountPos; j++) {
                System.arraycopy(bTrack, position, buff, 0, longSize);
                final long objIdToMove = Longs.fromByteArray(buff);
                movedIds[objectsCounter++] = moveObject(objIdToMove, track[i + 1]);
                position += longSize;
            }
        }

        return movedIds;
    }

    /**
     * @param objId for object to start shift from
     * @return array of shifting coordinates or empty for impossible shift
     */
    @NotNull
    public byte[] getShiftTrackCandidate(final long objId, final double volume) {
        final int entryCoordinate = objIdCoordinate.get(objId);
        final int vector = getRandomVPVector(objId, entryCoordinate);
        final int stepShiftRestriction = 5;

        final double rnd = objRandoms.get(objId).nextDouble();
        byte[] objRandom = Longs.toByteArray(Double.doubleToRawLongBits(rnd));

        byte[] tracksOut = new byte[1];
        byte[] objectsOut = new byte[1];
        double currentVolume = volume;
        int currentTrackPosition = entryCoordinate;

        for (int i = 0; i <= stepShiftRestriction; i++) {
            currentTrackPosition += vector;
            final TIntHashSet restricted = borders.get(currentTrackPosition);
            if (restricted != null && restricted.contains(vector)
                    || i == stepShiftRestriction)
                return new byte[0];

            final CellState cellState = cube[currentTrackPosition];
            final long[] objToMove = cellState.getMovingCandidates(currentVolume);
            ++tracksOut[0];
            tracksOut = Bytes.concat(tracksOut, Ints.toByteArray(currentTrackPosition));

            final byte objectsCount = (byte) objToMove.length;
            if (objectsCount == 0)
                break;

            objectsOut[0] += objectsCount;
            currentVolume = cellState.getSumVolume(objToMove);
            objectsOut = Bytes.concat(objectsOut, new byte[]{objectsCount});
            for (long id : objToMove)
                objectsOut = Bytes.concat(objectsOut, Longs.toByteArray(id));
        }

        return Bytes.concat(tracksOut, objectsOut, objRandom);
    }

    /**
     * @param bTracks are 2 dimensional array with flat coordinates and objIds each
     * @return indexes of allowed tracks
     */
    @NotNull
    public int[] resolveCollisions(@NotNull final byte[][] bTracks) {
        final int tracksCount = bTracks.length;
        final int[][] collisionMatrix = new int[tracksCount][tracksCount];
        final int[] collisionCounter = new int[tracksCount];
        final int[] wVector = new int[tracksCount];                             // weighted vector of collision count
                                                                                // for each [track * length] of track
        final int[][] tracks = new int[tracksCount][];
        final double[] rnd = new double[tracksCount];
        IntStream.range(0, bTracks.length)
                .parallel()
                .forEach(i -> {
                    tracks[i] = extractTrack(bTracks[i]);
                    wVector[i] = getObjectCountInTrack(bTracks[i]);
                    rnd[i] = getGeneratedRandom(bTracks[i]);
                });

        // counting collisions here in parallel by track
        IntStream.range(0, tracks.length)
                .parallel()
                .forEach(i -> {
                    for (int j = 0; j < tracks.length; j++) {
                        if (i == j)
                            continue;
                        final int collision = checkTracksCollision(tracks[i], tracks[j]);
                        collisionMatrix[i][j] = collision;
                        collisionCounter[i] += collision;
                    }
                });

        // calculating weights for each tracks
        for (int i = 0; i < wVector.length; i++)
            wVector[i] *= collisionCounter[i];

        // weight comparison
        IntStream.range(0, collisionMatrix.length)
                .parallel()
                .forEach(i -> {
                    for (int j = 0; j < collisionMatrix.length; j++) {
                        if (i == j)
                            continue;

                        if (wVector[i] < wVector[j] || (wVector[i] == wVector[j] && rnd[i] < rnd[j]))
                            collisionMatrix[i][j] = 0;
                    }
                });

        // filtering tracks with 0 collisions and sorting them by coordinate value
        // for avoiding further concurrent objects creation
        return IntStream.range(0, collisionMatrix.length)
                .parallel()
                .filter(i -> IntStream.of(collisionMatrix[i]).allMatch(e -> e == 0))
                .boxed()
                .sorted(Comparator.comparingInt(i -> tracks[i][0]))
                .mapToInt(i -> i)
                .toArray();
    }
    //endregion

    //region Private Methods

    @Contract(pure = true)
    private int getFlatCoordinate(int x, int y, int z) {
        return size * (x * size + y) + z;
    }

    private int[] restoreCoordinateFromFlat(int flatCoordinate) {
        final int[] restored = new int[3];
        restored[2] = flatCoordinate % size;
        restored[1] = flatCoordinate / size % size;
        restored[0] = flatCoordinate / sizeSquare;
        return restored;
    }

    /**
     * @param trackA 1 dimensional array with flat coordinates
     * @param trackB 1 dimensional array with flat coordinates
     * @return 1 if there is a collision, 0 if not
     */
    @Contract(pure = true)
    private int checkTracksCollision(@NotNull final int[] trackA, @NotNull final int[] trackB) {
        for (int item : trackA)
            for (int value : trackB)
                if (item == value)
                    return 1;

        return 0;
    }

    // TODO: select coordinate with respect to neighbors free volumes
    private int getRandomVPVector(final long objId, final int coordinate) {
        final int[][] vectors = new int[2][26];
        int freeIdx = 0;
        int busyIdx = 0;

        final TIntHashSet restricted = borders.get(coordinate);

        for (int i = 0; i < nbh.length; i++) {
            if (restricted != null && restricted.contains(nbh[i]))
                continue;

            if (cube[coordinate + nbh[i]].occupiedVolume.get() == 0)
                vectors[0][freeIdx++] = i;
            else
                vectors[1][busyIdx++] = i;
        }

        return (freeIdx > 0) ? nbh[vectors[0][objRandoms.get(objId).nextInt(freeIdx)]]
                : nbh[vectors[1][objRandoms.get(objId).nextInt(busyIdx)]];
    }

    @NotNull
    private int[] extractTrack(@NotNull final byte[] track) {
        final int trackLength = track[0];
        final int[] result = new int[trackLength];
        int position = 1;
        byte[] buff = new byte[intSize];

        for (int i = 0; i < trackLength; i++) {
            System.arraycopy(track, position, buff, 0, intSize);
            result[i] = Ints.fromByteArray(buff);
            position += intSize;
        }

        return result;
    }

    private int getObjectCountInTrack(@NotNull final byte[] track) {
        return track[1 + track[0] * intSize];
    }

    private double getGeneratedRandom(@NotNull final byte[] track) {
        byte[] buff = new byte[doubleSize];
        System.arraycopy(track, track.length - doubleSize, buff, 0, doubleSize);
        return Double.longBitsToDouble(Longs.fromByteArray(buff));
    }
    //endregion
}