// meris L1 transekte
// RD 20140219

cdir='C:\Users\Marco\Desktop\landsat8_nn_breadboard_20141002\';

execnn=0;
if execnn ==0 then
    cd 'C:\Users\Marco\Desktop\landsat8_nn_breadboard_20141002\';
    exec('./nnhs.sce');
    exec('./nnhs_ff.sce');
end
execnn=1;

lam_L8_5 = [440.0	480.0	560.0	655.0	865.0];
radiance_mult =	[0.012163,0.012455,0.011477, 0.0096784, 0.0059227];
radiance_add  = [-60.81587, -62.27619, -57.38698, -48.39193, -29.61345];

reflectance_mult = 2.05e-5;
reflectance_add  = -0.1;	

panchromatic_lines =	16301;		
panchromatic_samples =	16101;		
reflective_lines =	8151;		
reflective_samples =	8051;	
thermal_lines =	8151;		
thermal_samples =	8051;

subsampling_x =	1;	
subsampling_y =	1;		
subregion_x =	3751;
subregion_y =	6812;		
x_off = subregion_x;
y_off = subregion_y;		

subregion_width =	6696;		
subregion_height =	5273;

center_image_x = (panchromatic_samples./2);			
center_image_y = (panchromatic_lines./2);
x_c=center_image_x;
y_c=center_image_y;	

zenith_wink_fak = 7/x_c;

absorb_ozon=[2.82e-003, 2.076e-002, 1.022e-001, 5.313e-002, 0.0]; // has to be adjusted


//netpath='G:\backup_notebook_d_20130823\projekte\nn_hs';
inv_ac_nn1 =nnhs('.\landsat8_20140914\atmo_invers31_08p5_1_65_75_fl_20140905\log_rw_l85\51x97x27_88584.8.net');
inv_iop_nn1=nnhs('.\landsat8_20140914\s2l5l8_iop_20140708\inv_l85_logrw_logiop_20140708_p5_fl\47x97x17_23343.3.net');


deg2rad=%pi/180.0;
rad2deg=180.0/%pi;

temperature=15.0;
salinity=30.0;

sun_azi=153.34929033;	
sun_elevation=53.84615741;
sun_zeni = 90-sun_elevation;

// read transect data

fpin=mopen('LC81970222013202LGN00_subset_helgoland_red_geometry_Mask2.txt');

// jump over 6 header lines
head=mgetl(fpin,7);

// read transect line by line
longi_a=[];
iop_nn_a=[];
rw_a_nn=[];
r_toa_a=[];

npix=2483-7;
//npix=20;
for ipix=1:npix,
    data=mfscanf(16,fpin,'%f');
    xb=data(1);
    yb=data(2);
    longi=data(3);
    lati =data(4);
    radiance_1_5=data(5:9);
    radiance_6_11=data(10:15);
    flags=data(16);
    
    longi_a(ipix)=longi;

    // rescale radiance to reflectance
    count_1_5 = ((radiance_1_5' - radiance_add)./radiance_mult)
    reflectance_1_5 = ((radiance_1_5' - radiance_add)./radiance_mult).*reflectance_mult+reflectance_add;
    r_toa_a(ipix,:)=reflectance_1_5;
    
//    log_rtosa=log(reflectance_1_5);
    
    a=lati*deg2rad
    alpha= 98.2 // Landsat 8 inclination
    alpha_rad=alpha*deg2rad
    cos_alpha=cos(alpha_rad)
    sin_beta=cos_alpha./(sin(%pi./2-a))
    beta_rad=asin(sin_beta)
    beta_deg=beta_rad*rad2deg

    bb = y_c-(yb*subsampling_y + y_off);
    ab= xb*subsampling_x + x_off - x_c;
    cb =sqrt(ab.^2 +bb.^2);
    alpha_rad = acos(bb./cb);
    db=asin(beta_rad+alpha_rad)*cb;
    view_zeni= ab*zenith_wink_fak;

    // viewin left of ascendng path, definition when standing on pixe looking to sensor
    if ab < 0.0 then
        view_azi=beta_deg + 90.0;
    else
        // right of image
        view_azi=beta_deg + 270;
    end
    cos_sun=cos(sun_zeni*deg2rad);
    cos_view=cos(view_zeni*deg2rad);
    sin_sun=sin(sun_zeni*deg2rad);
    sin_view=sin(view_zeni*deg2rad);

    cos_azi_diff=cos((view_azi-sun_azi)*deg2rad);
    azi_diff_rad=acos(cos((view_azi-sun_azi)*deg2rad)); 
    sin_azi_diff=sin(azi_diff_rad);
    azi_diff_deg=azi_diff_rad*rad2deg;

    x=sin_view*cos_azi_diff;
    y=sin_view*sin_azi_diff;
    z=cos_view;

// toa_tosa corrections

        r_tosa=reflectance_1_5;

        //*** ozone correction ***/
        ozone = 350;
        atm_press=1013;
        model_ozone=0;

        trans_ozoned5=exp(-(absorb_ozon.*ozone./1000-model_ozone)./cos_sun);
        trans_ozoneu5=exp(-(absorb_ozon.*ozone./1000-model_ozone)./cos_view);
        trans_ozone5 =trans_ozoned5 .*trans_ozoneu5;

        r_tosa=r_tosa./trans_ozone5;


        // pressure correction

        rayl_mass_toa_tosa = atm_press - 1013.2;
        rayl_rel_mass_toa_tosa =(atm_press - 1013.2) / 1013.2;

        //* calculate phase function for rayleigh path radiance*/
        cos_scat_ang = -cos_view * cos_sun - sin_view * sin_sun * cos_azi_diff; // scattering angle without fresnel reflection
        cos_gamma_plus = cos_view * cos_sun - sin_view * sin_sun * cos_azi_diff; // for fresnel reflection
        phase_rayl_plus = 0.75 * (1.0 + cos_gamma_plus * cos_gamma_plus);
        phase_rayl_min = 0.75 * (1.0 + cos_scat_ang * cos_scat_ang);

        tau_rayl_standard = 0.008735 * (lam_L8_5 ./ 1000.0).^(-4.08);//* lam in µm */
        tau_rayl_toa_tosa = tau_rayl_standard * rayl_mass_toa_tosa/1013.2;

        r_rayl_toa_tosa   = cos_sun*tau_rayl_toa_tosa .* phase_rayl_min ./ (4 .* %pi) .* (1.0/cos_view ).*%pi;

        trans_rayl_press=exp(-tau_rayl_toa_tosa.*(1.0/cos_view + 1.0/ cos_sun));
        trans_rayl_pressd=exp(-tau_rayl_toa_tosa.*(1.0/ cos_sun));
        trans_rayl_pressu=exp(-tau_rayl_toa_tosa.*(1.0/cos_view));


        r_tosa = r_tosa ./ trans_rayl_press;

        r_tosa = r_tosa - r_rayl_toa_tosa;
        
        r_tosa_a(ipix,:)=r_tosa;
        
        log_r_tosa=log(r_tosa);


// end toa_tosa corrections

    nn_in=[sun_zeni,x,y,z,temperature,salinity,log_r_tosa];
    // NNcompute rw from rtosa
    log_rw=nnhs_ff(inv_ac_nn1,nn_in);
    rw=exp(log_rw);
    rw_a_nn(ipix,:)=rw;
    // NN compute IOPs from rw
    nn_in_inv=[sun_zeni view_zeni azi_diff_deg temperature salinity log_rw];
    log_result_nn1=nnhs_ff(inv_iop_nn1,nn_in_inv);
    iop_nn_a(ipix,:)=exp(log_result_nn1);
    
end
mclose(fpin);

hed='Test of Landsat 8 NN';
scf();
plot(longi_a, rw_a_nn(:,3),'+')
xtitle(hed,'longitude [deg]','rw_nn 560 nm [dl]');

scf();
plot(longi_a,iop_nn_a(:,4:5),'+');
xtitle(hed,'longitude [deg]','atot [m-1]');
legend(['b_part' 'b_wit']);

tsm=(iop_nn_a(:,4)+iop_nn_a(:,5)).*1.7;
scf();
plot(longi_a,tsm,'+');
xtitle(hed,'longitude [deg]','tsm [g m-3]');

// total a
a_tot=iop_nn_a(:,1)+iop_nn_a(:,2)+iop_nn_a(:,3);
scf();
plot(longi_a,a_tot,'+');
xtitle(hed,'longitude [deg]','atot [m-1]');

// chlorophyll
chl=iop_nn_a(:,1).^1.04 .* 20.0;
scf();
plot(longi_a,chl,'+');
xtitle(hed,'longitude [deg]','chl [mg m-3]');

scf();
plot(longi_a, r_toa_a(:,1), '+',longi_a, rw_a_nn(:,1),'+');
xtitle(hed,'longitude [deg]','r [dl]');
legend(['r_toa 440 nm' 'r_w 440 nm']);

