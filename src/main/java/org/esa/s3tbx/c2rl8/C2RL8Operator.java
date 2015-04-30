package org.esa.s3tbx.c2rl8;

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.pointop.PixelOperator;
import org.esa.beam.framework.gpf.pointop.ProductConfigurer;
import org.esa.beam.framework.gpf.pointop.Sample;
import org.esa.beam.framework.gpf.pointop.SampleConfigurer;
import org.esa.beam.framework.gpf.pointop.WritableSample;
import org.esa.beam.framework.ui.BooleanExpressionConverter;
import org.esa.beam.jai.ResolutionLevel;
import org.esa.beam.jai.VirtualBandOpImage;
import org.esa.beam.util.ProductUtils;
import org.esa.beam.util.StringUtils;

import java.awt.Rectangle;
import java.awt.image.Raster;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;


// todo (mp) - Add flags band and check for OOR of inputs and outputs of the NNs.
// todo (mp) - add min/max values of NN inputs and outputs to metadata

@OperatorMetadata(alias = "C2RL8", version = "0.1-SNAPSHOT",
                  authors = "Roland Doerffer, Marco Peters (Brockmann Consult)",
                  category = "Optical Processing/Thematic Water Processing",
                  copyright = "Copyright (C) 2014 by Brockmann Consult",
                  description = "Performs atmospheric correction and IOP retrieval on Landsat 8 data products.")
public class C2RL8Operator extends PixelOperator {

    static final String DEFAULT_ATMO_NN_RESOURCE = "/auxdata/nets/landsat8_20140914/atmo_invers31_08p5_1_65_75_fl_20140905/log_rw_l85/51x97x27_88584.8.net";
    static final String DEFAULT_IOP_NN_RESOURCE = "/auxdata/nets/landsat8_20140914/s2l5l8_iop_20140708/inv_l85_logrw_logiop_20140708_p5_fl/47x97x17_23343.3.net";
    private static final String[] EXPECTED_BANDNAMES = new String[]{"coastal_aerosol", "blue", "green", "red", "near_infrared"};

    private final BandConfig[] WATER_OUT_BANDS = new BandConfig[]{
            new BandConfig("conc_apig", "m^-1", "Pigment absorption coefficient"),
            new BandConfig("conc_adet", "m^-1", "Pigment absorption"),
            new BandConfig("conc_agelb", "m^-1", "Yellow substance absorption coefficient"),
            new BandConfig("conc_bpart", "m^-1", ""),
            new BandConfig("conc_bwit", "m^-1", "Backscattering of suspended particulate matter"),
            new BandConfig("tsm", "g m^-3", "Total suspended matter dry weight concentration",
                           "(conc_bpart + conc_bwit) * 1.7"),
            new BandConfig("atot", "m^-1", "Total absorption coefficient of all water constituents",
                           "conc_apig + conc_adet + conc_agelb"),
            new BandConfig("chl", "mg m^-3", "Chlorophyll concentration",
                           "pow(conc_apig, 1.04) * 20.0")
    };

    @SourceProduct(label = "Landsat 8 product", description = "The landsat 8 source product.")
    private Product source;

    @Parameter(label = "Valid pixel expression", defaultValue = "water_confidence_mid OR water_confidence_high",
               converter = BooleanExpressionConverter.class)
    private String validPixelExpression;

    @Parameter(defaultValue = "1013.25", unit = "hPa", interval = "(500, 1500)")
    private double pressure;

    @Parameter(defaultValue = "350.0", unit = "DU", interval = "(100, 1000)")
    private double ozone;

    @Parameter(label = "Use alternative neural net for AC", defaultValue = "false")
    private boolean useAlternativeAtmoNet;

    // todo (mp) - the initial enable state of the components not correct
    // todo (mp) - textboxes should not expand till infinity if file path is selected
    @Parameter(label = "Alternative neural net file for AC")
    private File alternativeAtmoNetFile;

    @Parameter(label = "Use alternative neural net for IOP retrieval", defaultValue = "false")
    private boolean useAlternativeIOPNet;

    // todo (mp) - the initial enable state of the components not correct
    // todo (mp) - textboxes should not expand till infinity if file path is selected
    @Parameter(label = "Alternative neural net file for IOP retrieval")
    private File alternativeIOPNetFile;

    @Parameter(label = "Include atmospherically corrected reflectances", defaultValue = "true")
    private boolean outputReflectances;

    private C2RL8Algo algo;
    private VirtualBandOpImage validMask;

    @Override
    protected void computePixel(int x, int y, Sample[] sourceSamples, WritableSample[] targetSamples) {

        if (isPixelValid(x, y)) {
            double[] radiance_1_5 = new double[sourceSamples.length];
            for (int i = 0; i < sourceSamples.length; i++) {
                radiance_1_5[i] = sourceSamples[i].getDouble();
            }

            GeoPos geoPos = source.getGeoCoding().getGeoPos(new PixelPos(x + 0.5f, y + 0.5f), null);
            C2RL8Result result = algo.runAlgo(x, y, geoPos.getLat(), geoPos.getLon(), radiance_1_5);

            int targetSamplesIndex = 0;
            for (int i = 0; i < result.log_iops.length; i++, targetSamplesIndex++) {
                targetSamples[i].set(Math.exp(result.log_iops[i]));
            }

            targetSamplesIndex += 3; // TODO - 3 virtual bands should not be in the targetSamples array

            if (outputReflectances) {
                for (int i = 0; i < result.log_rw.length; i++) {
                    targetSamples[targetSamplesIndex + i].set(Math.exp(result.log_rw[i]));
                }
            }
        } else {
            for (WritableSample targetSample : targetSamples) {
                if (targetSample.getNode() != null) {
                    try {
                        targetSample.set(Float.NaN);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }
        }

    }

    private boolean isPixelValid(int x, int y) {
        if (validMask != null) {
            Raster validMaskData = validMask.getData(new Rectangle(x, y, 1, 1));
            return validMaskData.getSample(x, y, 0) != 0;
        }else {
            return true;
        }
    }


    @Override
    protected void prepareInputs() throws OperatorException {
        super.prepareInputs();
        validateSource();
        if (!source.isCompatibleBandArithmeticExpression(validPixelExpression)) {
            throw new OperatorException("The given valid pixel expression can not be used with the given source product.");
        }


        MetadataElement metadataRoot = source.getMetadataRoot();
        MetadataElement l1MetadataFile = metadataRoot.getElement("L1_METADATA_FILE");
        MetadataElement imageAttributes = l1MetadataFile.getElement("IMAGE_ATTRIBUTES");
        double sun_azimuth = imageAttributes.getAttribute("SUN_AZIMUTH").getData().getElemDouble();
        double sun_elevation = imageAttributes.getAttribute("SUN_ELEVATION").getData().getElemDouble();
        double sun_zenith = 90 - sun_elevation;


        MetadataElement radiometricRescaling = l1MetadataFile.getElement("RADIOMETRIC_RESCALING");
        double[] radiance_add = new double[EXPECTED_BANDNAMES.length];
        double[] radiance_mult = new double[EXPECTED_BANDNAMES.length];
        double[] reflectance_add = new double[EXPECTED_BANDNAMES.length];
        double[] reflectance_mult = new double[EXPECTED_BANDNAMES.length];
        double[] lam_L8_5 = new double[EXPECTED_BANDNAMES.length];
        for (int i = 0; i < EXPECTED_BANDNAMES.length; i++) {
            Band sourceBand = source.getBand(EXPECTED_BANDNAMES[i]);
            radiance_add[i] = sourceBand.getScalingOffset();
            radiance_mult[i] = sourceBand.getScalingFactor();
            reflectance_add[i] = radiometricRescaling.getAttributeDouble(String.format("REFLECTANCE_ADD_BAND_%d", i + 1));
            reflectance_mult[i] = radiometricRescaling.getAttributeDouble(String.format("REFLECTANCE_MULT_BAND_%d", i + 1));
            lam_L8_5[i] = sourceBand.getSpectralWavelength();
        }

        double subsampling_x = 1;
//        double subsampling_y = 1;
        double x_off = 0;
//        double y_off = 0;
//        double subregion_width = source.getSceneRasterWidth();
//        double subregion_height = source.getSceneRasterHeight();
        MetadataElement history = metadataRoot.getElement("history");
        if (history != null) {
            MetadataElement subsetInfo = history.getElement("SubsetInfo");
            if (subsetInfo != null) {
                subsampling_x = subsetInfo.getAttributeInt("SubSampling.x", 1);
//                subsampling_y = subsetInfo.getAttributeInt("SubSampling.y", 1);
                x_off = subsetInfo.getAttributeInt("SubRegion.x", 0);
//                y_off = subsetInfo.getAttributeInt("SubRegion.y", 0);
//                subregion_width = subsetInfo.getAttributeInt("SubRegion.width", source.getSceneRasterWidth());
//                subregion_height = subsetInfo.getAttributeInt("SubRegion.height", source.getSceneRasterHeight());
            }
        }

        double x_c = metadataRoot.getProduct().getSceneRasterWidth() / 2;
//        double y_c = metadataRoot.getProduct().getSceneRasterHeight() / 2;

        NeuralNet atmoNeuralNet = getNeuralNet(useAlternativeAtmoNet, alternativeAtmoNetFile, DEFAULT_ATMO_NN_RESOURCE);
        NeuralNet iopNeuralNet = getNeuralNet(useAlternativeIOPNet, alternativeIOPNetFile, DEFAULT_IOP_NN_RESOURCE);

        algo = new C2RL8Algo(radiance_add, radiance_mult, reflectance_add, reflectance_mult, lam_L8_5, subsampling_x, x_off, x_c,
                             sun_azimuth, sun_zenith, pressure, ozone, atmoNeuralNet, iopNeuralNet);

        if (StringUtils.isNotNullAndNotEmpty(validPixelExpression)) {
            validMask = VirtualBandOpImage.createMask(validPixelExpression, source, ResolutionLevel.MAXRES);
        }
    }

    @Override
    protected void configureSourceSamples(SampleConfigurer sampleConfigurer) throws OperatorException {
        for (int i = 0; i < EXPECTED_BANDNAMES.length; i++) {
            sampleConfigurer.defineSample(i, EXPECTED_BANDNAMES[i]);
        }
    }


    @Override
    protected void configureTargetProduct(ProductConfigurer productConfigurer) {
        super.configureTargetProduct(productConfigurer);
        productConfigurer.copyMetadata();

        Product targetProduct = productConfigurer.getTargetProduct();
        int rasterWidth = targetProduct.getSceneRasterWidth();
        int rasterHeight = targetProduct.getSceneRasterHeight();
        for (BandConfig bandConfig : WATER_OUT_BANDS) {
            Band targetBand;
            if (bandConfig.asVirtual) {
                targetBand = new VirtualBand(bandConfig.name, ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight, bandConfig.expression);
                targetProduct.addBand(targetBand);
            } else {
                targetBand = targetProduct.addBand(bandConfig.name, ProductData.TYPE_FLOAT32);
            }
            targetBand.setUnit(bandConfig.unit);
            targetBand.setDescription(bandConfig.description);
            targetBand.setNoDataValue(bandConfig.noData);
            targetBand.setNoDataValueUsed(true);
        }

        if (outputReflectances) {
            for (int i = 0; i < EXPECTED_BANDNAMES.length; i++) {
                Band sourceBand = source.getBand(EXPECTED_BANDNAMES[i]);
                String description = String.format("Water leaving radiance reflectance at %f.2 nm", sourceBand.getSpectralWavelength());
                Band targetBand = productConfigurer.addBand(String.format("reflec_%d", i + 1), ProductData.TYPE_FLOAT32);
                ProductUtils.copySpectralBandProperties(sourceBand, targetBand);
                targetBand.setUnit("sr^-1");
                targetBand.setDescription(description);
                targetBand.setNoDataValue(Float.NaN);
                targetBand.setNoDataValueUsed(true);
            }
            targetProduct.setAutoGrouping("reflec");
        }

        ProductUtils.copyFlagBands(source, targetProduct, true);

    }

    @Override
    protected void configureTargetSamples(SampleConfigurer sampleConfigurer) throws OperatorException {
        for (int i = 0; i < WATER_OUT_BANDS.length; i++) {
            BandConfig bandConfig = WATER_OUT_BANDS[i];
            if (!bandConfig.asVirtual) {
                sampleConfigurer.defineSample(i, bandConfig.name);
            }
        }
        if (outputReflectances) {
            for (int i = 0; i < EXPECTED_BANDNAMES.length; i++) {
                sampleConfigurer.defineSample(WATER_OUT_BANDS.length + i, String.format("reflec_%d", i + 1));
            }
        }
    }

    private void validateSource() {
        if (source.getGeoCoding() == null) {
            throw new OperatorException("The source product must be geo-referenced.");
        }
        for (String expectedBandname : EXPECTED_BANDNAMES) {
            if (!source.containsBand(expectedBandname)) {
                throw new OperatorException(String.format("The band '%s' is missing. " +
                                                          "The source product is expected to be a Landsat 8 data product.", expectedBandname));
            }
        }
    }

    private NeuralNet getNeuralNet(boolean useAlternativeNet, File alternativeFile, String defaultNnResource) {
        if (useAlternativeNet) {
            try (InputStream nnStream = new FileInputStream(alternativeFile)) {
                return NeuralNet.read(nnStream);
            } catch (IOException e) {
                throw new OperatorException(String.format("Could not load neural net file.%n'%s'", alternativeAtmoNetFile), e);
            }
        } else {
            return loadNNFromResource(defaultNnResource);

        }
    }

    static NeuralNet loadNNFromResource(String defaultNnResource) {
        try (InputStream nnStream = C2RL8Operator.class.getResourceAsStream(defaultNnResource)) {
            if (nnStream != null) {
                return NeuralNet.read(nnStream);
            } else {
                throw new OperatorException(String.format("Could not load default neural net.%n'%s'", defaultNnResource));
            }
        } catch (IOException e) {
            throw new OperatorException(String.format("Could not load neural net file.%n'%s'", defaultNnResource), e);
        }
    }

    public static class Spi extends OperatorSpi {

        public Spi() {
            super(C2RL8Operator.class);
        }
    }

    private class BandConfig {

        String name;
        String description;
        String unit;
        String expression;
        boolean asVirtual;
        float noData;

        private BandConfig(String name, String unit, String description) {
            this(name, unit, description, null);
        }

        private BandConfig(String name, String unit, String description, String expression) {
            this.name = name;
            this.unit = unit;
            this.description = description;
            this.expression = expression;
            this.asVirtual = expression != null;
            this.noData = Float.NaN;
        }
    }
}
