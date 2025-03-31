%Based on http://pulseq.github.io/writeEpiDiffusionRS.html by Prof. Maxim Zaitsev
%and the Pulseq team
% this is an experimental high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping
% it further features diffusion weighting using the standard
% Stejskal-Tanner scheme
%
% IMPORTANT NOTICE: be aware, that this sequence potentially uses very
% strong gradient that may overload your scanner!

%%%%SETTINGS%%%%
close all
seq=mr.Sequence();         % Create a new sequence object
fov=300e-3; Nx=130; Ny=130;  % Define FOV and resolution
thickness=(8e-3);            % slice thickness
distfac = (3e-3);
Nslices=3;                 % number of slices
bFactor=350; % s/mm^2   % b value
bFactorRef = 50; % s/mm^2 b reference b value
Repdelay = 0; % set at 0 if trigger signal is used, else desired repetition delay
TE = 50e-3; % echo time
small_delta_calc = 0.0043; %desired diffusion lobe duration (siemens definition 1 ramp + flat time), this value should be tuned
                           %based on desired TE until PNS limit is reached for a given sequence. b value calculations are
                           %automatically adjusted.
Ndirections = 3; % number of diffusion directions
pe_enable=1;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os=1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
readoutTime = 6.6e-4; % this controls the readout bandwidth
partFourierFactor = 0.5;    % partial Fourier factor: 1: full sampling 0: start with ky=0
filename = 'C:\dummydir\dummy.seq'; % desired output filepath
tRFex=3e-3;
tRFref=3e-3;
num_avgs = 1; %number of averages
longseq = "no";  %suppresses graphical output to avoid crashes, useful for long sequences
zerodir = "no";  %enable if a sequence without diff gradients is desired (untested feature)
% Set system limits
lims = mr.opts('MaxGrad',200,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89); % hardware limits

gradrate = 120; % defines separate slew rate for diffusion gradients
gradrate = gradrate*42576000; %proper units by gamma multiplication

STIR = "yes"; %set to "yes" if STIR fatsat should be used instead of saturation
STIRdelay = 220e-3; %in ms 205-225 ms range for 3T (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4359893/) Del Grande F et al. Fat-suppression techniques for 3-T MR imaging of the musculoskeletal system. Radiographics. 2014 Jan-Feb;34(1):217-33. doi: 10.1148/rg.341135130
ReducedFOV = "yes"; % shifts slice selection gradient of the 90 degree pulse to freq encoding direction to achieve inner volume excitation
partoforiginalFOV = 150e-3; % size of the desired inner FOV

%Inner Volume Implementation based on:
%https://pubs.rsna.org/doi/epdf/10.1148/radiology.156.3.4023236 Inner volume MR imaging: technical concepts and their application. D A Feinberg et al. Radiology 1985 156:3, 743-747 
%https://www.sciencedirect.com/science/article/pii/S0730725X13001884 Christopher J. Wargo et al. A comparison and evaluation of reduced-FOV methods for multi-slice 7T human imaging, Magnetic Resonance Imaging, Volume 31, Issue 8, 2013, Pages 1349-1359, ISSN 0730-725X, https://doi.org/10.1016/j.mri.2013.05.003.

NoFatSat = "no"; %disables fat saturation if set to "yes"

GRAPPA = "yes"; % enables GRAPPA
AccelFac = 3; % GRAPPA acceleration factor
ACSamount = 0; %number of central kspace lines used for autocalibration, kept at 0 since ACS lines are taken in a separate scan

NeoACStemp = ceil(Ny-(Ny/2)*partFourierFactor); %calculates number of ACS lines used in separate scan
NeoACS = ceil(NeoACStemp/AccelFac)+1;

multibreath = "yes";  %set to "yes" if multiple sequence files for different breathholds are used, recommended setting
dirsperbreath = 3;  %must be same number as Ndirections, gives the amount of diffusion directions in a given breathhold
part = 1; %index of the current breathhold number, i.e. 9 diffusion directions split across 3 breathholds will have disperbreath 3 and part indexes 1,2,3

directions = load("directions.mat"); %loads diffusion directions based on golden means
                                     %(Feng, L. Golden-Angle Radial MRI: Basics, Advances, and Applications. J Magn Reson Imaging, 56: 45-62.DOI: 10.1002/jmri.28187)  
directions = directions.coords;      %any direction schemes may be inserted here instead


%%%%%END OF SETTINGS%%%%%



% Create fat-sat pulse
B0= 2.89; % 1.5 2.89 3.0
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*lims.gamma;
sat_freq = -407.0; %siemens value from idea
rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gz_fs = mr.makeTrapezoid('z',lims,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',lims,'Duration',tRFex,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);
if isequal(ReducedFOV,"yes")
    gzReph.area = 0;
end
% Create 180 degree slice refocusing pulse and gradients
[rf180, gz180] = mr.makeSincPulse(pi,'system',lims,'Duration',tRFref,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');
[~, gzr_t, gzr_a]=mr.makeExtendedTrapezoidArea('z',gz180.amplitude,0,-gzReph.area+0.5*gz180.amplitude*gz180.fallTime,lims);
gz180n=mr.makeExtendedTrapezoid('z','system',lims,'times',[0 gz180.riseTime gz180.riseTime+gz180.flatTime+gzr_t]+gz180.delay, 'amplitudes', [0 gz180.amplitude gzr_a]);



% Create 180 degree STIR pulse and gradient
[rfSTIR, gzSTIR, gzRephSTIR] = mr.makeSincPulse(pi,'system',lims,'Duration',tRFex,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% Create ReducedFOV gradient



[rfdummy, gzy, gzyReph] = mr.makeSincPulse(pi/2,'system',lims,'Duration',tRFex,...
    'SliceThickness',partoforiginalFOV,'apodization',0.5,'timeBwProduct',4);
gzy.channel = 'y';   %was formerly y!
gzyReph.channel = 'y';   %was formerly y!

% define the output trigger to play out with every slice excitatuion

trig = mr.makeTrigger('physio1','duration', 1000e-6);  % possible channels: 'osc0','osc1','ext1'


if isequal(zerodir,"yes")
    small_delta_calc = 0;
end




% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
blip_durdouble = ceil(2*sqrt((deltak*AccelFac)/lims.maxSlew)/10e-6/2)*10e-6*2;
% the split code below fails if this really makes a trapezoid instead of a triangle...
gy = mr.makeTrapezoid('y',lims,'Area',-deltak,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
gydouble = mr.makeTrapezoid('y',lims,'Area',-deltak*AccelFac,'Duration',blip_durdouble);

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;

kWidthneo = NeoACS*deltak;
gxneo = mr.makeTrapezoid('x',lims,'Area',kWidthneo+extra_area,'duration',readoutTime+blip_dur);
actual_area=gxneo.area-gxneo.amplitude/gxneo.riseTime*blip_dur/2*blip_dur/2/2-gxneo.amplitude/gxneo.fallTime*blip_dur/2*blip_dur/2/2;
gxneo.amplitude=gxneo.amplitude/actual_area*kWidth;
gxneo.area = gxneo.amplitude*(gxneo.flatTime + gxneo.riseTime/2 + gxneo.fallTime/2);
gxneo.flatArea = gxneo.amplitude*gxneo.flatTime;




% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO

% split the blip into two halves and produnce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
gy_partsdouble = mr.splitGradientAt(gydouble, blip_durdouble/2, lims);
[gy_blipup, gy_blipdown]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
[gy_blipupdouble, gy_blipdowndouble]=mr.align('right',gy_partsdouble(1),'left',gy_partsdouble(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);
gy_blipdownupdouble = mr.addGradients({gy_blipdowndouble, gy_blipupdouble}, lims);
gy_longshort = mr.addGradients({gy_blipdowndouble, gy_blipup}, lims);
gy_shortlong = mr.addGradients({gy_blipdown, gy_blipupdouble}, lims);
% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

% phase encoding and partial Fourier

Ny_pre=round(partFourierFactor*Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
Ny_meas=Ny_pre+Ny_post;

% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2);
gxPreminus = mr.makeTrapezoid('x',lims,'Area',gx.area/2);
gyPredummy = mr.makeTrapezoid('y',lims,'Area',Ny_pre*deltak);
gyPre = mr.makeTrapezoid('y',lims,'Area',Ny_pre*deltak);
[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
[gxPreminus,gyPredummy]=mr.align('right',gxPreminus,'left',gyPredummy);
% relax the PE prepahser to reduce stimulation
gyPre = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
gyPre.amplitude=gyPre.amplitude*pe_enable;
gyPredummy = mr.makeTrapezoid('y',lims,'Area',gyPredummy.area,'Duration',mr.calcDuration(gxPreminus,gyPredummy));
gyPredummy.amplitude=gyPredummy.amplitude*pe_enable;





% Calculate delay times
durationToCenter = ((Ny_pre)+0.5)*mr.calcDuration(gx);
if isequal(GRAPPA,"yes")
    durationToCenter = ((((Ny_pre-ACSamount/2)/AccelFac)+((ACSamount/2)))+0.5)*mr.calcDuration(gx);
end
Ny_pre_orig = Ny_pre;
gxPreneo = mr.makeTrapezoid('x',lims,'Area',-gxneo.area/2);
linesuntilzeroline = ceil((durationToCenter/mr.calcDuration(gx))-0.5);
gyPreneo = mr.makeTrapezoid('y',lims,'Area',(linesuntilzeroline)*deltak);
gyPreneo.amplitude=gyPreneo.amplitude*pe_enable;


if ACSamount == 0
    Ny_pre = 0;
end

rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
delayTE1=ceil((TE/2 - mr.calcDuration(rf,gz) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
delayTE2tmp=ceil((TE/2 - mr.calcDuration(rf180,gz180n) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
assert(delayTE1>=0);
gxPre.delay=0;
gyPre.delay=0;
delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
[gxPreminus,gyPredummy]=mr.align('right',gxPreminus,'left',gyPredummy);
assert(delayTE2>=0);

% diffusion weithting calculation
% delayTE2 is our window for small_delta
% delayTE1+delayTE2-delayTE2 is our big delta
% we anticipate that we will use the maximum gradient amplitude, so we need
% to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
small_delta=delayTE2-ceil(lims.maxGrad/gradrate/lims.gradRasterTime)*lims.gradRasterTime;


big_delta=delayTE1+mr.calcDuration(rf180,gz180n);
g=sqrt(bFactor*1e6/bFactCalc(1,(1.4)*small_delta,(1.4)*big_delta)); % for now it looks too large!
gr=ceil(g/gradrate/lims.gradRasterTime)*lims.gradRasterTime;
gro = ceil(g/lims.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
grdiff = gr-gro;

%diffusion gradient poles

 gDiffzero=mr.makeTrapezoid('z','amplitude',0,'riseTime',gr,'flatTime',small_delta-gr,'system',lims);
 mingDiffzero = mr.makeTrapezoid('z','amplitude',0,'riseTime',2*gr,'flatTime',2*(small_delta-gr),'system',lims);

% rt = gDiff.riseTime;
rt = ceil(lims.maxGrad/gradrate/lims.gradRasterTime)*lims.gradRasterTime;
calc_rt = rt; %manual value

if isequal(zerodir,"yes")
    rt = 0;
end

assert(small_delta_calc>=rt);

mocompscaling = (17*small_delta_calc)/6 + (17*rt)/6 + ((small_delta_calc + rt)*(49*small_delta_calc + 97*rt))^(1/2)/6;
mocompscaling = mocompscaling/(big_delta);
calc_big_delta = mocompscaling*big_delta;
calc_big_delta = round(calc_big_delta/lims.gradRasterTime)*lims.gradRasterTime;
calc_small_delta = small_delta_calc;
gamma = 42.577478518;

calc_b = bFactor*(1e6);
calc_bzero = bFactorRef*(1e6);



%Stejskal-Tanner equation derived in full, Kuchel et al., 2012, Concepts Magn. Reson., 40A: 205-214.  DOI: 10.1002/cmr.a.21241
calc_g = sqrt(calc_b/(((2*pi)^2)*(-calc_small_delta^3+2*calc_small_delta^2*calc_rt+calc_big_delta*calc_small_delta^2-1.5*calc_small_delta*calc_rt^2+(3/10)*calc_rt^3)));
calc_gzero = sqrt(calc_bzero/(((2*pi)^2)*(-calc_small_delta^3+2*calc_small_delta^2*calc_rt+calc_big_delta*calc_small_delta^2-1.5*calc_small_delta*calc_rt^2+(3/10)*calc_rt^3)));


assert(calc_g<=lims.maxGrad)
assert(calc_gzero<=lims.maxGrad)


g = calc_g; %original g not used beyond this point
gr = calc_rt;

%Golden Ratio Diffusion Directions
%based on Feng, L. Golden-Angle Radial MRI: Basics, Advances, and Applications. J Magn Reson Imaging, 56: 45-62.DOI: 10.1002/jmri.28187

g1 = 0.4656;
g2 = 0.6823;
r = g; %maintain direction gradient magnitude
coords = zeros(3,Ndirections);
for m = 1:Ndirections

    phi = 2*pi*mod(m*g2,1);
    theta = (pi/2) - acos(mod(m*g1,1));
    [x,y,z] = sph2cart(phi,theta,r);
    coords(:,m) = [x,y,z];
end
coords = coords./r;

if isequal(multibreath,"yes")
    coords = directions(:,dirsperbreath*(part-1)+1:dirsperbreath*part);
end


gDiffzero=mr.makeTrapezoid('z','amplitude',calc_gzero,'riseTime',gr,'flatTime',round(calc_small_delta-gr,4),'system',lims);
mingDiffzero = mr.makeTrapezoid('z','amplitude',-calc_gzero,'riseTime',2*gr,'flatTime',2*round(calc_small_delta-gr,4),'system',lims);
gDiffzerox=mr.makeTrapezoid('x','amplitude',calc_gzero,'riseTime',gr,'flatTime',round(calc_small_delta-gr,4),'system',lims);
mingDiffzerox = mr.makeTrapezoid('x','amplitude',-calc_gzero,'riseTime',2*gr,'flatTime',2*round(calc_small_delta-gr,4),'system',lims);
gDiffzeroy=mr.makeTrapezoid('y','amplitude',calc_gzero,'riseTime',gr,'flatTime',round(calc_small_delta-gr,4),'system',lims);
mingDiffzeroy = mr.makeTrapezoid('y','amplitude',-calc_gzero,'riseTime',2*gr,'flatTime',2*round(calc_small_delta-gr,4),'system',lims);
truegDiffzero=mr.makeTrapezoid('z','amplitude',0,'riseTime',gr,'flatTime',round(calc_small_delta-gr,4),'system',lims);
truemingDiffzero = mr.makeTrapezoid('z','amplitude',0,'riseTime',2*gr,'flatTime',2*round(calc_small_delta-gr,4),'system',lims);

delayTE1 = delayTE1 - mr.calcDuration(gDiffzero);
delayTE2 = delayTE2 - mr.calcDuration(mingDiffzero);
actualdelta = mr.calcDuration(gDiffzerox)+mr.calcDuration(rf180,gz180n)+delayTE1;

deltadjust = actualdelta-calc_big_delta;
assert(deltadjust>=0);
delayTE1 = delayTE1 - deltadjust;
assert(delayTE1>=0);

delayTE1 = (delayTE1/lims.gradRasterTime)*lims.gradRasterTime;


gxamplitude_orig = gx.amplitude;

% Define sequence blocks
for s=1:Nslices
        if ~isequal(STIR,"yes") && isequal(NoFatSat,"no")
            seq.addBlock(trig);
            seq.addBlock(rf_fs,gz_fs);
        end
        rfSTIR.freqOffset=gzSTIR.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        rf.freqOffset = 0;
        rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        if isequal(STIR,"yes") && isequal(NoFatSat,"no")
            seq.addBlock(trig);
            seq.addBlock(mr.makeDelay(STIRdelay),rfSTIR,gzSTIR);
        end
        if ~isequal(ReducedFOV,"yes")
            if isequal(NoFatSat,"no")
                seq.addBlock(rf,gz);
            else
                seq.addBlock(trig);
                seq.addBlock(rf,gz);
            end
        else
            if isequal(NoFatSat,"no")
                seq.addBlock(rf,gzy);
                seq.addBlock(gzyReph);
            else
                seq.addBlock(trig);
                seq.addBlock(rf,gzy);
                seq.addBlock(gzyReph);
            end
        end
        seq.addBlock(mr.makeDelay(deltadjust));
        seq.addBlock(truegDiffzero);
        seq.addBlock(mr.makeDelay(delayTE1),truemingDiffzero);
        rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        seq.addBlock(rf180,gz180n);
        seq.addBlock(truemingDiffzero);
        seq.addBlock(mr.makeDelay(delayTE2),truegDiffzero);
        if isequal(GRAPPA,"yes")
            gx.amplitude = gxamplitude_orig;
        end
        seq.addBlock(gxPre);
        %Full Reference Scan
        for i=1:Ny_meas

        seq.addBlock(gx,adc); % Read the first line of k-space with a single half-blip at the end
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end
        if Repdelay > 0
            seq.addBlock(mr.makeDelay(Repdelay))
        end
end
for s=1:Nslices
        if ~isequal(STIR,"yes") && isequal(NoFatSat,"no")
            seq.addBlock(trig);
            seq.addBlock(rf_fs,gz_fs);
        end
        rfSTIR.freqOffset=gzSTIR.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        rf.freqOffset = 0;
        rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        if isequal(STIR,"yes") && isequal(NoFatSat,"no")
            seq.addBlock(trig);
            seq.addBlock(mr.makeDelay(STIRdelay),rfSTIR,gzSTIR);
        end
        if ~isequal(ReducedFOV,"yes")
            if isequal(NoFatSat,"no")
                seq.addBlock(rf,gz);
            else
                seq.addBlock(trig);
                seq.addBlock(rf,gz);
            end
        else
            if isequal(NoFatSat,"no")
                seq.addBlock(rf,gzy);
                seq.addBlock(gzyReph);
            else
                seq.addBlock(trig);
                seq.addBlock(rf,gzy);
                seq.addBlock(gzyReph);
            end
        end

        seq.addBlock(mr.makeDelay(deltadjust));
        seq.addBlock(gDiffzero);
        seq.addBlock(mr.makeDelay(delayTE1),mingDiffzero);
        rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
        seq.addBlock(rf180,gz180n);
        seq.addBlock(mingDiffzero);
        seq.addBlock(mr.makeDelay(delayTE2),gDiffzero);
        seq.addBlock(gxPre,gyPreneo);
        testinc = 0;
        for i=1:NeoACS
            if i==1
                seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
            elseif i==NeoACS
                seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
            elseif isequal(GRAPPA,"yes")
                if mod(i,1)~=0 && ~ismember(i,[(Ny_pre_orig-(NeoACS/2)):(Ny_pre_orig+(NeoACS/2))])
                    continue
                else
                    if ismember(i,[(Ny_pre_orig-(NeoACS/2)):(Ny_pre_orig+(NeoACS/2))])
                        if i == Ny_pre_orig-(NeoACS/2)
                            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        elseif i == Ny_pre_orig+(NeoACS/2)
                            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        else
                            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    else
                        seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    end
                end
            else
                seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
            end
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
            testinc = testinc +1;
        end
        if Repdelay > 0
            seq.addBlock(mr.makeDelay(Repdelay))
        end
end

for u = 1:num_avgs
    for s=1:Nslices
        for k = 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isequal(GRAPPA,"yes")
                gx.amplitude = gxamplitude_orig;
            end
            if ~isequal(STIR,"yes") && isequal(NoFatSat,"no")
                seq.addBlock(trig);
                seq.addBlock(rf_fs,gz_fs);
            end
            rfSTIR.freqOffset=gzSTIR.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            rf.freqOffset = 0;
            rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            if isequal(STIR,"yes") && isequal(NoFatSat,"no")
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(STIRdelay),rfSTIR,gzSTIR);
            end
            if ~isequal(ReducedFOV,"yes")
                if isequal(NoFatSat,"no")
                    seq.addBlock(rf,gz);
                else
                    seq.addBlock(trig);
                    seq.addBlock(rf,gz);
                end
            else
                if isequal(NoFatSat,"no")
                    seq.addBlock(rf,gzy);
                    seq.addBlock(gzyReph);
                else
                    seq.addBlock(trig);
                    seq.addBlock(rf,gzy);
                    seq.addBlock(gzyReph);
                end
            end

            seq.addBlock(mr.makeDelay(deltadjust));
            seq.addBlock(gDiffzero);
            seq.addBlock(mr.makeDelay(delayTE1),mingDiffzero);
            rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            seq.addBlock(rf180,gz180n);
            seq.addBlock(mingDiffzero);
            seq.addBlock(mr.makeDelay(delayTE2),gDiffzero);
            seq.addBlock(gxPre,gyPre);
            testinc = 0;
            for i=1:Ny_meas
                if i==1
                    seq.addBlock(gx,gy_blipupdouble,adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    seq.addBlock(gx,gy_blipdowndouble,adc); % Read the last line of k-space with a single half-blip at the beginning
                elseif isequal(GRAPPA,"yes")
                    if mod(i,AccelFac)~=0 && ~ismember(i,[(Ny_pre-(ACSamount/2)):(Ny_pre+(ACSamount/2))])
                        continue
                    else
                        if ismember(i,[(Ny_pre-(ACSamount/2)):(Ny_pre+(ACSamount/2))])
                            if i == Ny_pre-(ACSamount/2)
                                seq.addBlock(gx,gy_longshort,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            elseif i == Ny_pre+(ACSamount/2)
                                seq.addBlock(gx,gy_shortlong,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            else
                                seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            end
                        else
                            seq.addBlock(gx,gy_blipdownupdouble,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                else
                    seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                testinc = testinc +1;
            end
            if Repdelay > 0
                seq.addBlock(mr.makeDelay(Repdelay))
            end
        end
    end
    for k = 1:Ndirections
        for s=1:Nslices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %define custom diffusion gradients
            tempcoords = coords(:,k);
            gDiffz=mr.makeTrapezoid('z','amplitude',g*tempcoords(3),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            mingDiffz = mr.makeTrapezoid('z','amplitude',-g*tempcoords(3),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);
            gDiffx=mr.makeTrapezoid('x','amplitude',g*tempcoords(1),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            mingDiffx = mr.makeTrapezoid('x','amplitude',-g*tempcoords(1),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);
            gDiffy=mr.makeTrapezoid('y','amplitude',g*tempcoords(2),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            mingDiffy = mr.makeTrapezoid('y','amplitude',-g*tempcoords(2),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);

            neggDiffz=mr.makeTrapezoid('z','amplitude',-g*tempcoords(3),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            negmingDiffz = mr.makeTrapezoid('z','amplitude',g*tempcoords(3),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);
            neggDiffx=mr.makeTrapezoid('x','amplitude',-g*tempcoords(1),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            negmingDiffx = mr.makeTrapezoid('x','amplitude',g*tempcoords(1),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);
            neggDiffy=mr.makeTrapezoid('y','amplitude',-g*tempcoords(2),'riseTime',rt,'flatTime',round(calc_small_delta-gr,4),'system',lims);
            negmingDiffy = mr.makeTrapezoid('y','amplitude',g*tempcoords(2),'riseTime',rt*2,'flatTime',round(calc_small_delta-gr,4)*2,'system',lims);

            if isequal(GRAPPA,"yes")
                gx.amplitude = gxamplitude_orig;
            end
            if ~isequal(STIR,"yes") && isequal(NoFatSat,"no")
                seq.addBlock(trig);
                seq.addBlock(rf_fs,gz_fs);
            end
            rfSTIR.freqOffset=gzSTIR.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            %rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
            rf.freqOffset = 0;
            rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            if isequal(STIR,"yes") && isequal(NoFatSat,"no")
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(STIRdelay),rfSTIR,gzSTIR);
            end
            if ~isequal(ReducedFOV,"yes")
                if isequal(NoFatSat,"no")
                    seq.addBlock(rf,gz);
                else
                    seq.addBlock(trig);
                    seq.addBlock(rf,gz);
                end
            else
                if isequal(NoFatSat,"no")
                    seq.addBlock(rf,gzy);
                    seq.addBlock(gzyReph);
                else
                    seq.addBlock(trig);
                    seq.addBlock(rf,gzy);
                    seq.addBlock(gzyReph);
                end
            end
            seq.addBlock(mr.makeDelay(deltadjust));
            seq.addBlock(gDiffz,gDiffx,gDiffy);
            seq.addBlock(mr.makeDelay(delayTE1),mingDiffz,mingDiffx,mingDiffy);
            rf180.freqOffset=gz180.amplitude*(thickness+distfac)*(s-1-(Nslices-1)/2);
            seq.addBlock(rf180,gz180n);%%%%CHANGE 1
            seq.addBlock(mingDiffz,mingDiffx,mingDiffy);
            seq.addBlock(mr.makeDelay(delayTE2),gDiffz,gDiffx,gDiffy);
            if isequal(GRAPPA,"yes")
                gx.amplitude = gxamplitude_orig;
            end
            seq.addBlock(gxPre,gyPre);
            testinc2 = 0;
             for i=1:Ny_meas
                if i==1
                    seq.addBlock(gx,gy_blipupdouble,adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    seq.addBlock(gx,gy_blipdowndouble,adc); % Read the last line of k-space with a single half-blip at the beginning
                elseif isequal(GRAPPA,"yes")
                    if mod(i,AccelFac)~=0 && ~ismember(i,[(Ny_pre-(ACSamount/2)):(Ny_pre+(ACSamount/2))])
                        continue
                    else
                        if ismember(i,[(Ny_pre-(ACSamount/2)):(Ny_pre+(ACSamount/2))])
                            if i == Ny_pre-(ACSamount/2)
                                seq.addBlock(gx,gy_longshort,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            elseif i == Ny_pre+(ACSamount/2)
                                seq.addBlock(gx,gy_shortlong,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            else
                                seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                            end
                        else
                            seq.addBlock(gx,gy_blipdownupdouble,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                else
                    seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                testinc2 = testinc2 + 1;
            end
            if Repdelay > 0
                seq.addBlock(mr.makeDelay(Repdelay))
            end
        end
    end
end
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
if ~isequal(longseq,"yes")
    seq.plot();             % Plot sequence waveforms
end


seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi');

%seq.write('first_order_mc_diff_SE.seq');
seq.write(filename);
if ~isequal(longseq,"yes")
    figure(2)
    plot3(coords(1,:),coords(2,:),coords(3,:),'o');
end
sz = size(coords);
for i = 1:sz(2)
    tempcoords = coords(:,i);
        line([0,tempcoords(1)],[0,tempcoords(2)],[0,tempcoords(3)],'Color','red');
end



%% plot k-space diagrams


% k-space trajectory calculation

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-space
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

rep = seq.testReport;
fprintf([rep{:}]); % as for January 2019 TR calculation fails for fat-sat

function b=bFactCalc(g, delta, DELTA)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39–65 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end