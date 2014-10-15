package org.esa.s3tbx.c2rl8;

import org.esa.beam.util.io.CsvReader;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.List;

import static org.junit.Assert.*;


/**
 * Tests are performed against the breadboard landsat8_trans_nn_20141008_1.sce
 * To get the correct results from the breadboard the value 'nAlpha' in the file nnhs.sce m,ust be set to 100000.
 */
public class C2RL8AlgoTest {

    private double[] radiance_add = {-60.81587, -62.27619, -57.38698, -48.39193, -29.61345};
    private double[] radiance_mult = {0.012163, 0.012455, 0.011477, 0.0096784, 0.0059227};
    private double[] reflectance_add = {-0.1, -0.1, -0.1, -0.1, -0.1};
    private double[] reflectance_mult = {2.05e-5, 2.05e-5, 2.05e-5, 2.05e-5, 2.05e-5};
    private double[] lam_L8_5 = {440.0, 480.0, 560.0, 655.0, 865.0};

    private C2RL8Algo algo;

    @Before
    public void setUp() throws Exception {
        NeuralNet atmoNet = C2RL8Operator.loadNNFromResource(C2RL8Operator.DEFAULT_ATMO_NN_RESOURCE);
        NeuralNet iopNet = C2RL8Operator.loadNNFromResource(C2RL8Operator.DEFAULT_IOP_NN_RESOURCE);
        algo = new C2RL8Algo(radiance_add, radiance_mult, reflectance_add, reflectance_mult, lam_L8_5,
                             1, 3751, 16101.0 / 2, 153.34929033, 90 - 53.84615741, 1013, 350, atmoNet, iopNet);
    }

    @Test
    public void testComputeReflectances() throws Exception {
        double[] radiances = new double[]{53.94367218017578, 42.87092590332031, 25.57164192199707, 12.07862377166748, 2.6770551204681396};

        double[] reflectances = C2RL8Algo.computeReflectances(radiances, radiance_add, radiance_mult, reflectance_add, reflectance_mult);

        assertEquals(0.093420259, reflectances[0], 1.0e-8);
        assertEquals(0.073064301, reflectances[1], 1.0e-8);
        assertEquals(0.048179119, reflectances[2], 1.0e-8);
        assertEquals(0.028083810, reflectances[3], 1.0e-8);
        assertEquals(0.011765808, reflectances[4], 1.0e-8);
    }

    @Test
    public void testComputeViewGeometry() throws Exception {
        C2RL8Algo.GeometryAngels geometryAngels = algo.computeViewGeometry(6138.5, 2176.5, 54.39365005493164);
        assertEquals(1839.0, geometryAngels.ab, 1.0e-6);
        assertEquals(1.5990311, geometryAngels.view_zenith, 1.0e-6);
        assertEquals(0.0279047, geometryAngels.sin_view, 1.0e-6);
        assertEquals(0.9764083, geometryAngels.sin_azi_diff, 1.0e-6);
        assertEquals(-0.0060255, geometryAngels.x, 1.0e-6);
        assertEquals(0.0272464, geometryAngels.y, 1.0e-6);
        assertEquals(0.9996105, geometryAngels.z, 1.0e-6);
    }


    @Test
    public void testConvertToa_Tosa() throws Exception {
        double[] radiances = new double[]{53.94367218017578, 42.87092590332031, 25.57164192199707, 12.07862377166748, 2.6770551204681396};

        double[] reflectances = C2RL8Algo.computeReflectances(radiances, radiance_add, radiance_mult, reflectance_add, reflectance_mult);
        C2RL8Algo.GeometryAngels geometryAngels = algo.computeViewGeometry(6138.5, 2176.5, 54.39365005493164);

        double[] r_tosa = algo.convertToa_Tosa(reflectances, geometryAngels);

        double[] expectedTosa = {0.093628872, 0.074265467, 0.052198667, 0.029279485, 0.011766503};
        assertArrayEquals(expectedTosa, r_tosa, 1.0e-8);

        }

    @Test
    public void testTransect() throws IOException {

        List<double[]> inputRecords = getInputRecords();
        List<double[]> expectedRecords = getExpectedResultRecords();

        for (int i = 0; i < inputRecords.size(); i++) {
            validateRecord(inputRecords, expectedRecords, i);
        }

    }

    private void validateRecord(List<double[]> inputRecords, List<double[]> expectedRecords, int i) {
        double[] record = inputRecords.get(i);
        double[] expected = expectedRecords.get(i);
        int x = (int) Math.floor(record[0]);
        int y = (int) Math.floor(record[1]);
        double lat = record[3];
        double lon = record[2];
        double[] radiance_1_5 = {record[4], record[5], record[6], record[7], record[8]};
        C2RL8Result result = algo.runAlgo(x, y, lat, lon, radiance_1_5);
        assertEquals(expected[0], Math.exp(result.log_iops[0]), 1.0e-2);
        assertEquals(expected[1], Math.exp(result.log_iops[1]), 1.0e-2);
        assertEquals(expected[2], Math.exp(result.log_iops[2]), 1.0e-2);
        assertEquals(expected[3], Math.exp(result.log_iops[3]), 1.0e-2);
        assertEquals(expected[4], Math.exp(result.log_iops[4]), 1.0e-2);
        assertEquals(expected[5], Math.exp(result.log_rw[0]), 1.0e-5);
        assertEquals(expected[6], Math.exp(result.log_rw[1]), 1.0e-5);
        assertEquals(expected[7], Math.exp(result.log_rw[2]), 1.0e-5);
        assertEquals(expected[8], Math.exp(result.log_rw[3]), 1.0e-5);
        assertEquals(expected[9], Math.exp(result.log_rw[4]), 1.0e-5);
    }

    private List<double[]> getInputRecords() throws IOException {
        return getRecords("LC81970222013202LGN00_subset_helgoland_red_geometry_Mask2.txt");
    }

    private List<double[]> getExpectedResultRecords() throws IOException {
        return getRecords("c2rl8_breadboard_20141008_1_results.txt");
    }

    private List<double[]> getRecords(String resourceName) throws IOException {
        Reader streamReader = new InputStreamReader(getClass().getResourceAsStream(resourceName));
        List<double[]> records;
        try (CsvReader reader = new CsvReader(streamReader, new char[]{'\t'}, true, "#")) {
            String[] header = reader.readRecord();
            records = reader.readDoubleRecords();
        }
        return records;
    }
}