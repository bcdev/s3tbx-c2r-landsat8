package org.esa.s3tbx.c2rl8;

/**
 * @author Marco Peters
 */
public class C2RL8Algo {

    private static final double[] ABSORB_OZONE = {2.82e-003, 2.076e-002, 1.022e-001, 5.313e-002, 0.0}; // has to be adjusted
    private final double[] radiance_add;
    private final double[] radiance_mult;
    private final double[] reflectance_add;
    private final double[] reflectance_mult;
    private final double[] lam_L8_5;
    private final double pressure;
    private final double ozone;
    private final NeuralNet atmoNeuralNet;
    private final NeuralNet iopNeuralNet;

    private final double sun_azimuth;
    private final double sun_zenith;

    private final double subsampling_x;
    private final double x_off;
    private final double x_c;

    private final double zenith_wink_fak;
    private final double temperature;
    private final double salinity;


    public C2RL8Algo(double[] radiance_add, double[] radiance_mult, double[] reflectance_add, double[] reflectance_mult, double[] lam_L8_5,
                     double subsampling_x, double x_off, double x_c, double sun_azimuth, double sun_zenith, double pressure, double ozone,
                     NeuralNet atmoNeuralNet, NeuralNet iopNeuralNet) {
        this.radiance_add = radiance_add;
        this.radiance_mult = radiance_mult;
        this.reflectance_add = reflectance_add;
        this.reflectance_mult = reflectance_mult;
        this.lam_L8_5 = lam_L8_5;
        this.subsampling_x = subsampling_x;
        this.x_off = x_off;
        this.x_c = x_c;
        this.pressure = pressure;
        this.ozone = ozone;
        this.sun_azimuth = sun_azimuth;
        this.sun_zenith = sun_zenith;
        this.atmoNeuralNet = atmoNeuralNet;
        this.iopNeuralNet = iopNeuralNet;
        this.zenith_wink_fak = 7.0 / this.x_c;

        // todo (mp)- load from waterradiance-auxdata per pixel
        this.temperature = 15.0;
        this.salinity = 30.0;

    }

    public C2RL8Result runAlgo(int x, int y, double lat, double lon, double[] radiance_1_5) {
        double[] reflectance_1_5 = computeReflectances(radiance_1_5, radiance_add, radiance_mult, reflectance_add, reflectance_mult);
        GeometryAngels geomAngels = computeViewGeometry(x + 0.5, y + 0.5, lat);
        double[] r_tosa = convertToa_Tosa(reflectance_1_5, geomAngels);
        double[] nn_in = new double[11];
        double[] atmoAuxdata = {sun_zenith, geomAngels.x, geomAngels.y, geomAngels.z, temperature, salinity};
        System.arraycopy(atmoAuxdata, 0, nn_in, 0, atmoAuxdata.length);
        for (int i = 0; i < r_tosa.length; i++) {
            nn_in[atmoAuxdata.length + i] = Math.log(r_tosa[i]);
        }

        C2RL8Result result = new C2RL8Result();
        // NNcompute rw from rtosa
        result.log_rw = atmoNeuralNet.get().calc(nn_in);

        // NN compute IOPs from rw
        double[] iopAuxdata = {sun_zenith, geomAngels.view_zenith, geomAngels.azi_diff_deg, temperature, salinity};
        double[] nn_in_inv = new double[10];
        System.arraycopy(iopAuxdata, 0, nn_in_inv, 0, iopAuxdata.length);
        System.arraycopy(result.log_rw, 0, nn_in_inv, iopAuxdata.length, result.log_rw.length);
        result.log_iops = iopNeuralNet.get().calc(nn_in_inv);
        return result;
    }

    public double[] convertToa_Tosa(double[] reflectance_1_5, GeometryAngels geomAngels) {
        double[] r_tosa = correctForOzone(reflectance_1_5, geomAngels);

        // pressure correction
        double rayl_mass_toa_tosa = pressure - 1013.2;
//        double rayl_rel_mass_toa_tosa =(pressure - 1013.2) / 1013.2;

        // calculate phase function for rayleigh path radiance
        double cos_scat_ang = -geomAngels.cos_view * geomAngels.cos_sun - geomAngels.sin_view * geomAngels.sin_sun * geomAngels.cos_azi_diff; // scattering angle without fresnel reflection
//        double cos_gamma_plus = cos_view * cos_sun - sin_view * sin_sun * cos_azi_diff; // for fresnel reflection
//        double phase_rayl_plus = 0.75 * (1.0 + cos_gamma_plus * cos_gamma_plus);
        double phase_rayl_min = 0.75 * (1.0 + cos_scat_ang * cos_scat_ang);

        for (int i = 0; i < lam_L8_5.length; i++) {
            double tau_rayl_standard = 0.008735 * Math.pow((lam_L8_5[i] / 1000.0), -4.08);  // lam in µm
            double tau_rayl_toa_tosa = tau_rayl_standard * rayl_mass_toa_tosa / 1013.2;
            double r_rayl_toa_tosa = geomAngels.cos_sun * tau_rayl_toa_tosa * phase_rayl_min / (4 * Math.PI) * (1.0 / geomAngels.cos_view) * Math.PI;
            double trans_rayl_press = Math.exp(-tau_rayl_toa_tosa * (1.0 / geomAngels.cos_view + 1.0 / geomAngels.cos_sun));
//            double trans_rayl_pressd = Math.exp(-tau_rayl_toa_tosa * (1.0 / cos_sun));
//            double trans_rayl_pressu = Math.exp(-tau_rayl_toa_tosa * (1.0 / cos_view));
            r_tosa[i] = r_tosa[i] / trans_rayl_press;
            r_tosa[i] = r_tosa[i] - r_rayl_toa_tosa;
        }

        return r_tosa;
    }

    private double[] correctForOzone(double[] reflectance_1_5, GeometryAngels geomAngels) {
        double[] r_tosa = new double[reflectance_1_5.length];
        double model_ozone = 0;
        for (int i = 0; i < r_tosa.length; i++) {
            double temp = -(ABSORB_OZONE[i] * ozone / 1000 - model_ozone);
            double trans_ozoned5 = Math.exp(temp / geomAngels.cos_sun);
            double trans_ozoneu5 = Math.exp(temp / geomAngels.cos_view);
            double trans_ozone5 = trans_ozoned5 * trans_ozoneu5;
            r_tosa[i] = reflectance_1_5[i] / trans_ozone5;
        }
        return r_tosa;
    }

    public static double[] computeReflectances(double[] radiance_1_5, double[] radiance_add_1_5, double[] radiance_mult_1_5,
                                               double[] reflectance_add, double[] reflectance_mult) {
        double[] refls = new double[radiance_1_5.length];
        for (int i = 0; i < refls.length; i++) {
            double count = (radiance_1_5[i] - radiance_add_1_5[i]) / radiance_mult_1_5[i];
            refls[i] = count * reflectance_mult[i] + reflectance_add[i];
        }
        return refls;
    }

    public GeometryAngels computeViewGeometry(double xb, double yb, double latitude) {
        GeometryAngels geomAngels = new GeometryAngels();
        double a = Math.toRadians(latitude);
        double alpha = 98.2; // Landsat 8 inclination
        double alpha_rad = Math.toRadians(alpha);
        double cos_alpha = Math.cos(alpha_rad);
        double sin_beta = cos_alpha / (Math.sin(Math.PI / 2 - a));
        double beta_rad = Math.asin(sin_beta);
        double beta_deg = Math.toDegrees(beta_rad);

//        double bb = y_c - (yb * subsampling_y + y_off);
        geomAngels.ab = xb * subsampling_x + x_off - x_c;
//        double cb = Math.sqrt(Math.pow(ab, 2) + Math.pow(bb, 2));
//        alpha_rad = Math.acos(bb / cb);
//        double db = Math.asin(beta_rad + alpha_rad) * cb;
        geomAngels.view_zenith = geomAngels.ab * zenith_wink_fak;
        double view_azimuth;
        // viewing left of ascending path, definition when standing on pixel looking to sensor
        if (geomAngels.ab < 0.0) {
            view_azimuth = beta_deg + 90.0;
        } else {
            // right of image
            view_azimuth = beta_deg + 270;
        }

        geomAngels.cos_sun = Math.cos(Math.toRadians(sun_zenith));
        geomAngels.cos_view = Math.cos(Math.toRadians(geomAngels.view_zenith));
        geomAngels.sin_sun = Math.sin(Math.toRadians(sun_zenith));
        geomAngels.sin_view = Math.sin(Math.toRadians(geomAngels.view_zenith));

        geomAngels.cos_azi_diff = Math.cos(Math.toRadians(view_azimuth - sun_azimuth));
        double azi_diff_rad = Math.acos(Math.cos(Math.toRadians(view_azimuth - sun_azimuth)));
        geomAngels.sin_azi_diff = Math.sin(azi_diff_rad);
        geomAngels.azi_diff_deg = Math.toDegrees(azi_diff_rad);

        geomAngels.x = geomAngels.sin_view * geomAngels.cos_azi_diff;
        geomAngels.y = geomAngels.sin_view * geomAngels.sin_azi_diff;
        geomAngels.z = geomAngels.cos_view;

        return geomAngels;
    }

    static class GeometryAngels {

        double ab;
        double view_zenith;
        double cos_sun;
        double sin_sun;
        double cos_view;
        double sin_view;
        double cos_azi_diff;
        double sin_azi_diff;
        double azi_diff_deg;
        double x;
        double y;
        double z;
    }
}
