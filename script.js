// script.js

// --- DOM Elements --- (No changes)
const sliders = {
    omega0: document.getElementById('omega0'),
    delta_omega: document.getElementById('delta_omega'),
    phi0: document.getElementById('phi0'),
    phi1: document.getElementById('phi1'),
    phi2: document.getElementById('phi2'),
    phi3: document.getElementById('phi3'),
};
const sliderVals = {
    omega0_val: document.getElementById('omega0_val'),
    delta_omega_val: document.getElementById('delta_omega_val'),
    phi0_val: document.getElementById('phi0_val'),
    phi1_val: document.getElementById('phi1_val'),
    phi2_val: document.getElementById('phi2_val'),
    phi3_val: document.getElementById('phi3_val'),
};
const canvasElements = {
    spectrumAndPhase: document.getElementById('spectrumChart'),
    individualCosines: document.getElementById('individualCosinesChart'),
    sumCosines: document.getElementById('sumCosinesChart'),
    sumIntensity: document.getElementById('sumIntensityChart'),
};
const showPeakConnectorCheckbox = document.getElementById('showPeakConnector');
const showSpectralPhaseCheckbox = document.getElementById('showSpectralPhase');
const showEnvelopeECheckbox = document.getElementById('showEnvelopeE');
const showEnvelopeICheckbox = document.getElementById('showEnvelopeI');
const gaussianModeRadio = document.getElementById('gaussianModeRadio');
const customModeRadio = document.getElementById('customModeRadio');
const gaussianControlsDiv = document.getElementById('gaussianControls');
const resetCustomSpectrumButton = document.getElementById('resetCustomSpectrumButton');
const autoOmega0CustomCheckbox = document.getElementById('autoOmega0Custom');
const omega0SliderGroup = document.getElementById('omega0SliderGroup');

// --- Chart Objects --- (No changes)
let charts = {};

// --- Global State --- (No changes)
let currentSpectrumMode = 'gaussian';
let customSpectrumDataPoints = [];
let isInitialSwitchToCustom = false;

// --- Global Constants --- (No changes)
const MIN_OMEGA_PHYSICAL = 0.01;
const NUM_T_POINTS = 401;
const T_MIN = -20;
const T_MAX = 20;
const BASE_NUM_COS_WAVES = 3;
const WAVES_TO_ADD_PER_INCREMENT = 2;
const NUM_DELTA_OMEGA_INCREMENTS = 10;
const INDIVIDUAL_COS_Y_OFFSET = 2.5;
const CUSTOM_SPECTRUM_POINT_DENSITY = 1.0;
const INDIVIDUAL_COSINE_MIN_AMPLITUDE = 0.05;

// --- Helper function: linspace --- (No changes)
function linspace(start, end, num) {
    const arr = [];
    if (num <= 0) return arr;
    if (num === 1) return [start];
    const step = (end - start) / (num - 1);
    for (let i = 0; i < num; i++) {
        arr.push(start + step * i);
    }
    return arr;
}

const timeArray = linspace(T_MIN, T_MAX, NUM_T_POINTS);
const NUM_OMEGA_PLOT_POINTS = 201;
const OMEGA_AXIS_MIN_DEFAULT = 0;
let OMEGA_AXIS_MAX_DEFAULT;
let omegaPlotArray;


// --- Catmull-Rom Spline Interpolation ---
function catmullRomInterpolate(p0, p1, p2, p3, t) {
    const t2 = t * t;
    const t3 = t2 * t;
    return 0.5 * (
        (2 * p1) +
        (-p0 + p2) * t +
        (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
        (-p0 + 3 * p1 - 3 * p2 + p3) * t3
    );
}

// --- Calculation Functions ---
function calculateSpectrumGaussian(omega, omega0, delta_omega_fwhm) { // No changes
    if (omega < MIN_OMEGA_PHYSICAL && omega0 < MIN_OMEGA_PHYSICAL && delta_omega_fwhm < 1e-5) return 0;
    if (omega < 0) return 0;
    if (delta_omega_fwhm <= 1e-6) return (Math.abs(omega - omega0) < 1e-6 && omega >= 0) ? 1 : 0;
    const factor = 2 * Math.sqrt(2 * Math.log(2));
    const sigma_effective = delta_omega_fwhm / factor;
    if (sigma_effective < 1e-9) return (Math.abs(omega - omega0) < 1e-6 && omega >= 0) ? 1 : 0;
    const exponent = -0.5 * Math.pow((omega - omega0) / sigma_effective, 2);
    return Math.exp(exponent);
}

// MODIFIED getCustomSpectrumValue to use Catmull-Rom
function getCustomSpectrumValue(omega) {
    if (!customSpectrumDataPoints || customSpectrumDataPoints.length === 0) return 0;

    const n = customSpectrumDataPoints.length;

    // Handle edge cases: omega outside the range of draggable points
    if (omega <= customSpectrumDataPoints[0].x) return Math.max(0, customSpectrumDataPoints[0].y);
    if (omega >= customSpectrumDataPoints[n - 1].x) return Math.max(0, customSpectrumDataPoints[n - 1].y);

    // Find the segment where omega lies for interpolation
    let segmentIndex = -1;
    for (let i = 0; i < n - 1; i++) {
        if (omega >= customSpectrumDataPoints[i].x && omega <= customSpectrumDataPoints[i + 1].x) {
            segmentIndex = i;
            break;
        }
    }

    if (segmentIndex === -1) { // Should be caught by edge cases, but as a fallback
        return Math.max(0, customSpectrumDataPoints[n - 1].y);
    }

    // For Catmull-Rom, we need 4 points: P0, P1, P2, P3
    // P1 and P2 are the points defining the current segment.
    // P0 is the point before P1, P3 is the point after P2.

    const p1 = customSpectrumDataPoints[segmentIndex];
    const p2 = customSpectrumDataPoints[segmentIndex + 1];

    // Handle endpoints for P0 and P3 by extrapolating or duplicating
    const p0 = (segmentIndex > 0) ? customSpectrumDataPoints[segmentIndex - 1] : {x: p1.x - (p2.x - p1.x), y: p1.y}; // Reflect or duplicate
    const p3 = (segmentIndex < n - 2) ? customSpectrumDataPoints[segmentIndex + 2] : {x: p2.x + (p2.x - p1.x), y: p2.y}; // Reflect or duplicate
    
    // Normalize omega (t) within the segment [P1.x, P2.x] to be [0, 1]
    // Important: Catmull-Rom t is usually for parameteric curves. For f(x), t is (x - x1) / (x2 - x1)
    let t;
    if (Math.abs(p2.x - p1.x) < 1e-9) { // Points are identical in x
        t = 0; // or 0.5 or 1, result should be p1.y or p2.y
    } else {
        t = (omega - p1.x) / (p2.x - p1.x);
    }
    t = Math.max(0, Math.min(1, t)); // Clamp t to [0, 1]

    const interpolatedY = catmullRomInterpolate(p0.y, p1.y, p2.y, p3.y, t);
    return Math.max(0, Math.min(1.1, interpolatedY)); // Clamp to valid amplitude range (0 to chart max)
}


function calculateSpectrum(omega, omega0_gaussian_center, delta_omega_fwhm_gaussian) {
    return currentSpectrumMode === 'custom' ? getCustomSpectrumValue(omega) : calculateSpectrumGaussian(omega, omega0_gaussian_center, delta_omega_fwhm_gaussian);
}

function calculateSpectralPhase(omega, omega0_ref, c0, c1, c2, c3) { // No changes
    const dw = omega - omega0_ref;
    return c0 + c1 * dw + c2 * Math.pow(dw, 2) + c3 * Math.pow(dw, 3);
}

function calculateEffectiveCustomFWHMAndCenter() { // No changes to its internal logic
    let sum_s = 0;
    let sum_os = 0;
    let sum_o2s = 0;

    // Use a denser sampling of the *smooth* custom spectrum for centroid calculation
    const numSamplesForCentroid = 100; // More samples for better accuracy
    const step = (OMEGA_AXIS_MAX_DEFAULT - OMEGA_AXIS_MIN_DEFAULT) / (numSamplesForCentroid -1);
    let prevOmega = OMEGA_AXIS_MIN_DEFAULT;
    let prevY = getCustomSpectrumValue(prevOmega); // Use smooth value

    for (let i = 1; i < numSamplesForCentroid; i++) {
        const currentOmega = OMEGA_AXIS_MIN_DEFAULT + i * step;
        const currentY = getCustomSpectrumValue(currentOmega); // Use smooth value
        
        const dx = currentOmega - prevOmega;
        const avg_y = (prevY + currentY) / 2;
        const avg_x = (prevOmega + currentOmega) / 2;

        sum_s += avg_y * dx;
        sum_os += avg_x * avg_y * dx;
        sum_o2s += avg_x * avg_x * avg_y * dx;

        prevOmega = currentOmega;
        prevY = currentY;
    }
    
    let centroid = parseFloat(sliders.omega0.value);
    if (sum_s > 1e-6) {
        centroid = sum_os / sum_s;
    } else if (customSpectrumDataPoints.length > 0) { // Fallback if spectrum is all zero
        let minX = Infinity, maxX = -Infinity, count = 0;
        customSpectrumDataPoints.forEach(p => {
            if (p.y > 1e-6) { // Consider points that were originally non-zero
                minX = Math.min(minX, p.x);
                maxX = Math.max(maxX, p.x);
                count++;
            }
        });
        if (count > 0) centroid = (minX + maxX) / 2;
        else if (customSpectrumDataPoints.length > 0) centroid = customSpectrumDataPoints[Math.floor(customSpectrumDataPoints.length / 2)].x; // Geometric center of all points
    }


    let s_max = 0;
    // Find s_max from the smooth curve
    for (let i=0; i < omegaPlotArray.length; i++) {
        const yVal = getCustomSpectrumValue(omegaPlotArray[i]);
        if (yVal > s_max) s_max = yVal;
    }


    if (s_max < 1e-3) return { center: centroid, fwhm: parseFloat(sliders.delta_omega.min) };

    const halfMax = s_max / 2;
    let omega1 = -1, omega2 = -1;

    // Find FWHM points on the smooth curve
    // Scan from left for omega1
    for (let i = 0; i < omegaPlotArray.length - 1; i++) {
        const x1 = omegaPlotArray[i];
        const y1 = getCustomSpectrumValue(x1);
        const x2 = omegaPlotArray[i+1];
        const y2 = getCustomSpectrumValue(x2);

        if ((y1 >= halfMax && y2 < halfMax) || (y1 < halfMax && y2 >= halfMax)) {
            if (Math.abs(y2 - y1) > 1e-6) {
                omega1 = x1 + (halfMax - y1) * (x2 - x1) / (y2 - y1);
            } else {
                omega1 = (x1 + x2) / 2;
            }
            break;
        }
    }
     if (omega1 === -1 && getCustomSpectrumValue(omegaPlotArray[0]) >= halfMax) omega1 = omegaPlotArray[0];


    // Scan from right for omega2
    for (let i = omegaPlotArray.length - 1; i > 0; i--) {
        const x1 = omegaPlotArray[i-1];
        const y1 = getCustomSpectrumValue(x1);
        const x2 = omegaPlotArray[i];
        const y2 = getCustomSpectrumValue(x2);
        if ((y1 < halfMax && y2 >= halfMax) || (y1 >= halfMax && y2 < halfMax)) {
            if (Math.abs(y2 - y1) > 1e-6) {
                omega2 = x1 + (halfMax - y1) * (x2 - x1) / (y2 - y1);
            } else {
                omega2 = (x1+x2)/2;
            }
            break;
        }
    }
    if (omega2 === -1 && getCustomSpectrumValue(omegaPlotArray[omegaPlotArray.length-1]) >= halfMax) omega2 = omegaPlotArray[omegaPlotArray.length-1];


    let fwhm;
    if (omega1 !== -1 && omega2 !== -1 && omega2 > omega1) {
        fwhm = Math.max(0.01, omega2 - omega1);
    } else if (sum_s > 1e-6) {
        const variance = (sum_o2s / sum_s) - Math.pow(centroid, 2);
        fwhm = (variance > 0) ? Math.max(0.01, 2.355 * Math.sqrt(variance)) : parseFloat(sliders.delta_omega.min);
    } else {
        fwhm = parseFloat(sliders.delta_omega.min);
    }
    return { center: centroid, fwhm: fwhm };
}

// --- Update and Plotting Logic ---
function updatePlots() {
    let omega0_slider_val_manual = parseFloat(sliders.omega0.value);
    let delta_omega_fwhm_gaussian = parseFloat(sliders.delta_omega.value);

    if (omega0_slider_val_manual < MIN_OMEGA_PHYSICAL) { 
        omega0_slider_val_manual = MIN_OMEGA_PHYSICAL; 
    } else if (omega0_slider_val_manual > parseFloat(sliders.omega0.max)) {
        omega0_slider_val_manual = parseFloat(sliders.omega0.max);
    }

    let max_delta_omega_fwhm_allowed = 2 * (omega0_slider_val_manual - MIN_OMEGA_PHYSICAL);
    if (max_delta_omega_fwhm_allowed < 0) max_delta_omega_fwhm_allowed = 0;
    if (delta_omega_fwhm_gaussian > max_delta_omega_fwhm_allowed) { delta_omega_fwhm_gaussian = max_delta_omega_fwhm_allowed; sliders.delta_omega.value = delta_omega_fwhm_gaussian; }
    if (delta_omega_fwhm_gaussian < parseFloat(sliders.delta_omega.min)) { delta_omega_fwhm_gaussian = parseFloat(sliders.delta_omega.min); sliders.delta_omega.value = delta_omega_fwhm_gaussian; }
    if (delta_omega_fwhm_gaussian < 0) { delta_omega_fwhm_gaussian = 0; sliders.delta_omega.value = delta_omega_fwhm_gaussian; }

    let omega0_phase_reference;
    let omega0_effective_center;
    let delta_omega_effective_fwhm;

    const customParams = calculateEffectiveCustomFWHMAndCenter();

    if (currentSpectrumMode === 'custom') {
        omega0_effective_center = customParams.center;
        delta_omega_effective_fwhm = customParams.fwhm;
        if (autoOmega0CustomCheckbox.checked) {
            omega0_phase_reference = customParams.center;
            sliders.omega0.value = omega0_phase_reference.toFixed(1); 
            sliderVals.omega0_val.value = omega0_phase_reference.toFixed(1);
        } else {
            omega0_phase_reference = omega0_slider_val_manual;
        }
    } else {
        omega0_phase_reference = omega0_slider_val_manual;
        omega0_effective_center = omega0_slider_val_manual;
        delta_omega_effective_fwhm = delta_omega_fwhm_gaussian;
    }
     if (omega0_phase_reference < MIN_OMEGA_PHYSICAL) omega0_phase_reference = MIN_OMEGA_PHYSICAL;


    const phi0_coeff = parseFloat(sliders.phi0.value), phi1_coeff = parseFloat(sliders.phi1.value), phi2_coeff = parseFloat(sliders.phi2.value), phi3_coeff = parseFloat(sliders.phi3.value);
    const showPeakConnector = showPeakConnectorCheckbox.checked, showSpectralPhase = showSpectralPhaseCheckbox.checked, showEnvE = showEnvelopeECheckbox.checked, showEnvI = showEnvelopeICheckbox.checked;

    sliderVals.omega0_val.value = omega0_phase_reference.toFixed(1);
    if (currentSpectrumMode === 'custom' && autoOmega0CustomCheckbox.checked) {
        sliders.omega0.value = omega0_phase_reference;
    }

    sliderVals.delta_omega_val.value = delta_omega_fwhm_gaussian.toFixed(1);
    sliderVals.phi0_val.value = phi0_coeff.toFixed(2); sliderVals.phi1_val.value = phi1_coeff.toFixed(2); sliderVals.phi2_val.value = phi2_coeff.toFixed(2); sliderVals.phi3_val.value = phi3_coeff.toFixed(2);

    const spectrumData = omegaPlotArray.map(om => ({ x: om, y: calculateSpectrum(om, omega0_phase_reference, delta_omega_fwhm_gaussian) }));
    const phaseData = omegaPlotArray.map(om => ({ x: om, y: calculateSpectralPhase(om, omega0_phase_reference, phi0_coeff, phi1_coeff, phi2_coeff, phi3_coeff) }));

    const deltaOmegaSliderMin = parseFloat(sliders.delta_omega.min);
    const deltaOmegaSliderMax = parseFloat(sliders.delta_omega.max);
    const deltaOmegaSliderRange = deltaOmegaSliderMax - deltaOmegaSliderMin;
    let numCosWavesBaseSampling;
    if (deltaOmegaSliderRange <= 1e-6) {
        numCosWavesBaseSampling = BASE_NUM_COS_WAVES;
    } else {
        const fwhmForProportion = delta_omega_effective_fwhm;
        const fractionOfRange = (fwhmForProportion - deltaOmegaSliderMin) / deltaOmegaSliderRange;
        const normalizedFraction = Math.max(0, Math.min(1, fractionOfRange));
        const currentIncrement = Math.floor(normalizedFraction * NUM_DELTA_OMEGA_INCREMENTS);
        const effectiveIncrement = Math.min(currentIncrement, NUM_DELTA_OMEGA_INCREMENTS - (NUM_DELTA_OMEGA_INCREMENTS > 0 ? 1 : 0));
        numCosWavesBaseSampling = BASE_NUM_COS_WAVES + (effectiveIncrement * WAVES_TO_ADD_PER_INCREMENT);
    }
    if (numCosWavesBaseSampling < 0) numCosWavesBaseSampling = 0;

    const individualCosinesDatasets = [];
    const sumCosinesDataY = new Array(NUM_T_POINTS).fill(0);
    const sumSinesDataY = new Array(NUM_T_POINTS).fill(0);
    let cosineOmegaValuesToPlot = [];

    if (currentSpectrumMode === 'custom') {
        for (let om_k = OMEGA_AXIS_MIN_DEFAULT; om_k <= OMEGA_AXIS_MAX_DEFAULT + 1e-9; om_k += 0.5) {
            if (Math.abs(om_k) < 1e-9) {
                const A_k_zero = getCustomSpectrumValue(0); // Uses smooth interpolation
                if (A_k_zero >= INDIVIDUAL_COSINE_MIN_AMPLITUDE) {
                    cosineOmegaValuesToPlot.push({ omega: 0, amplitude: A_k_zero });
                }
                continue;
            }
            if (om_k < MIN_OMEGA_PHYSICAL) continue;
            const A_k = getCustomSpectrumValue(om_k); // Uses smooth interpolation
            if (A_k >= INDIVIDUAL_COSINE_MIN_AMPLITUDE) {
                cosineOmegaValuesToPlot.push({ omega: om_k, amplitude: A_k });
            }
        }
    } else {
        let omega_sampling_min_gauss = omega0_effective_center - delta_omega_effective_fwhm / 2;
        let omega_sampling_max_gauss = omega0_effective_center + delta_omega_effective_fwhm / 2;
        if (omega_sampling_min_gauss < MIN_OMEGA_PHYSICAL) { omega_sampling_min_gauss = MIN_OMEGA_PHYSICAL; if (omega_sampling_max_gauss < omega_sampling_min_gauss) omega_sampling_max_gauss = omega_sampling_min_gauss; }
        if (omega_sampling_max_gauss <= omega_sampling_min_gauss) {
            if (delta_omega_effective_fwhm < 1e-5) {
                omega_sampling_min_gauss = Math.max(MIN_OMEGA_PHYSICAL, omega0_effective_center - 0.1);
                omega_sampling_max_gauss = omega0_effective_center + 0.1;
            } else {
                omega_sampling_min_gauss = Math.max(MIN_OMEGA_PHYSICAL, omega0_effective_center - 0.05);
                omega_sampling_max_gauss = omega0_effective_center + 0.05;
            }
            if (omega_sampling_min_gauss < MIN_OMEGA_PHYSICAL) omega_sampling_min_gauss = MIN_OMEGA_PHYSICAL;
            if (omega_sampling_max_gauss <= omega_sampling_min_gauss) omega_sampling_max_gauss = omega_sampling_min_gauss + 0.01;
        }

        let tempGaussianOmegas;
        if (numCosWavesBaseSampling === 0 || omega_sampling_max_gauss <= omega_sampling_min_gauss || delta_omega_effective_fwhm < 1e-6) {
            const centerFreqForSampling = Math.max(omega0_effective_center, omega_sampling_min_gauss);
            tempGaussianOmegas = (numCosWavesBaseSampling > 0 && centerFreqForSampling >=MIN_OMEGA_PHYSICAL) ? new Array(numCosWavesBaseSampling).fill(centerFreqForSampling) : [];
        } else if (numCosWavesBaseSampling === 1) {
             const om_cand = Math.max(omega0_effective_center, omega_sampling_min_gauss);
             tempGaussianOmegas = (om_cand >= MIN_OMEGA_PHYSICAL || Math.abs(om_cand) < 1e-9) ? [om_cand] : [];
        } else {
            tempGaussianOmegas = linspace(omega_sampling_min_gauss, omega_sampling_max_gauss, numCosWavesBaseSampling);
        }
        tempGaussianOmegas.forEach(om_k => {
            if (om_k >= MIN_OMEGA_PHYSICAL || Math.abs(om_k) < 1e-9) {
                const A_k = calculateSpectrumGaussian(om_k, omega0_phase_reference, delta_omega_fwhm_gaussian);
                if (A_k > 1e-5) {
                    cosineOmegaValuesToPlot.push({ omega: om_k, amplitude: A_k });
                }
            }
        });
    }

    const peakConnectorData = [];
    let plottedCosCount = 0;

    for (let i = 0; i < cosineOmegaValuesToPlot.length; i++) {
        const om_k = cosineOmegaValuesToPlot[i].omega;
        const A_k = cosineOmegaValuesToPlot[i].amplitude;
        const phi_k = calculateSpectralPhase(om_k, omega0_phase_reference, phi0_coeff, phi1_coeff, phi2_coeff, phi3_coeff);

        timeArray.forEach((t, t_idx) => {
            // MODIFICATION: Changed argument for cosine and sine to match common optics convention
            const arg = phi_k - om_k * t; 
            sumCosinesDataY[t_idx] += A_k * Math.cos(arg);
            sumSinesDataY[t_idx] += A_k * Math.sin(arg); // Envelope calculation should use the same argument for consistency
        });
        
        if (Math.abs(om_k) < 1e-9) continue; // Skip DC component for individual plotting if desired, or handle separately

        if (showPeakConnector) {
            // MODIFICATION: Adjusted t_peak_k for the new convention
            let t_peak_k = (Math.abs(om_k) > 1e-6) ? phi_k / om_k : 0; 
            t_peak_k = Math.max(T_MIN, Math.min(T_MAX, t_peak_k));
            peakConnectorData.push({ x: t_peak_k, y: A_k - plottedCosCount * INDIVIDUAL_COS_Y_OFFSET });
        }
        // MODIFICATION: Changed argument for individual cosine waves
        const waveYValues = timeArray.map(t => A_k * Math.cos(phi_k - om_k * t)); 
        individualCosinesDatasets.push({
            label: `ω=${om_k.toFixed(2)}, A=${A_k.toFixed(2)}`,
            data: timeArray.map((t, t_idx) => ({ x: t, y: waveYValues[t_idx] - plottedCosCount * INDIVIDUAL_COS_Y_OFFSET })),
            borderColor: `hsl(${(plottedCosCount * (360 / (Math.max(1, cosineOmegaValuesToPlot.filter(p => Math.abs(p.omega)>1e-9).length) ))) % 360}, 70%, 50%)`,
            borderWidth: 1, fill: false,
        });
        plottedCosCount++;
    }

    if (showPeakConnector && peakConnectorData.length > 0) {
        individualCosinesDatasets.push({ label: 'Peak Connector', data: peakConnectorData, borderColor: '#000000', borderWidth: 2, fill: false, showLine: true, pointRadius: 3, pointBackgroundColor: '#000000', order: 99 });
    }

    const sumCosinesData = timeArray.map((t, idx) => ({ x: t, y: sumCosinesDataY[idx] }));
    const sumIntensityData = timeArray.map((t, idx) => ({ x: t, y: Math.pow(sumCosinesDataY[idx], 2) }));
    // Envelope calculation uses sumCosinesDataY and sumSinesDataY which are now based on the new 'arg'
    const envelopeEData = timeArray.map((t, idx) => ({ x: t, y: Math.sqrt(Math.pow(sumCosinesDataY[idx], 2) + Math.pow(sumSinesDataY[idx], 2)) }));
    const envelopeIData = envelopeEData.map(point => ({ x: point.x, y: Math.pow(point.y, 2) }));

    charts.spectrumAndPhase.data.datasets[0].data = spectrumData;
    charts.spectrumAndPhase.data.datasets[1].data = phaseData;
    charts.spectrumAndPhase.data.datasets[1].hidden = !showSpectralPhase;
    charts.spectrumAndPhase.data.datasets[2].data = customSpectrumDataPoints;
    charts.spectrumAndPhase.data.datasets[2].hidden = (currentSpectrumMode !== 'custom');
    charts.spectrumAndPhase.options.scales.x.min = OMEGA_AXIS_MIN_DEFAULT;
    charts.spectrumAndPhase.options.scales.x.max = OMEGA_AXIS_MAX_DEFAULT;
    charts.spectrumAndPhase.update('none');
    charts.individualCosines.data.datasets = individualCosinesDatasets;
    charts.individualCosines.update('none');
    charts.sumCosines.data.datasets[0].data = sumCosinesData;
    charts.sumCosines.data.datasets[1].data = envelopeEData;
    charts.sumCosines.data.datasets[1].hidden = !showEnvE;
    charts.sumCosines.update('none');
    charts.sumIntensity.data.datasets[0].data = sumIntensityData;
    charts.sumIntensity.data.datasets[1].data = envelopeIData;
    charts.sumIntensity.data.datasets[1].hidden = !showEnvI;
    charts.sumIntensity.update('none');
}

// --- Custom Spectrum Point Management --- (No changes)
function initializeOrUpdateCustomSpectrumPoints(isInitialSwitch = false) {
    const previousCustomPoints = currentSpectrumMode === 'custom' ? [...customSpectrumDataPoints] : [];
    customSpectrumDataPoints = [];
    const density = CUSTOM_SPECTRUM_POINT_DENSITY;
    if (density <= 0) return;
    const omegaStep = 1.0 / density;
    let currentOmega = OMEGA_AXIS_MIN_DEFAULT;
    let safetyCounter = 0;
    while (currentOmega <= OMEGA_AXIS_MAX_DEFAULT + 1e-9 && safetyCounter < 1000) {
        let yVal = 0.5;
        if (isInitialSwitch) {
            const omega0_gauss = parseFloat(sliders.omega0.value), delta_omega_gauss = parseFloat(sliders.delta_omega.value);
            yVal = calculateSpectrumGaussian(currentOmega, omega0_gauss, delta_omega_gauss);
        } else if (previousCustomPoints.length > 0) yVal = getCustomSpectrumValueInterpolate(currentOmega, previousCustomPoints); // This is linear for point re-init
        customSpectrumDataPoints.push({ x: currentOmega, y: Math.max(0, Math.min(1.1, yVal)) });
        currentOmega += omegaStep; safetyCounter++;
    }
    if (customSpectrumDataPoints.length > 0 && Math.abs(customSpectrumDataPoints[customSpectrumDataPoints.length - 1].x - OMEGA_AXIS_MAX_DEFAULT) > 1e-3 && customSpectrumDataPoints[customSpectrumDataPoints.length - 1].x < OMEGA_AXIS_MAX_DEFAULT) {
        let yVal = 0.5;
        if (isInitialSwitch) {
            const omega0_gauss = parseFloat(sliders.omega0.value), delta_omega_gauss = parseFloat(sliders.delta_omega.value);
            yVal = calculateSpectrumGaussian(OMEGA_AXIS_MAX_DEFAULT, omega0_gauss, delta_omega_gauss);
        } else if (previousCustomPoints.length > 0) yVal = getCustomSpectrumValueInterpolate(OMEGA_AXIS_MAX_DEFAULT, previousCustomPoints); // Linear for point re-init
        customSpectrumDataPoints.push({ x: OMEGA_AXIS_MAX_DEFAULT, y: Math.max(0, Math.min(1.1, yVal)) });
    }
    customSpectrumDataPoints = customSpectrumDataPoints.filter(p => p.x <= OMEGA_AXIS_MAX_DEFAULT + 1e-9);
    customSpectrumDataPoints = customSpectrumDataPoints.filter((point, index, self) => index === self.findIndex((p) => Math.abs(p.x - point.x) < 1e-9));
    customSpectrumDataPoints.sort((a, b) => a.x - b.x);
    if (charts.spectrumAndPhase) charts.spectrumAndPhase.data.datasets[2].data = customSpectrumDataPoints;
}

// This function is used when re-initializing points after density change (which is now fixed, so less critical)
// It remains linear for simplicity in that specific context. The main getCustomSpectrumValue is now Catmull-Rom.
function getCustomSpectrumValueInterpolate(omega, pointsArray) {
    if (!pointsArray || pointsArray.length === 0) return 0.5;
    const firstPoint = pointsArray[0], lastPoint = pointsArray[pointsArray.length - 1];
    if (omega <= firstPoint.x) return Math.max(0, firstPoint.y);
    if (omega >= lastPoint.x) return Math.max(0, lastPoint.y);
    for (let i = 0; i < pointsArray.length - 1; i++) {
        const p1 = pointsArray[i], p2 = pointsArray[i + 1];
        if (omega >= p1.x && omega <= p2.x) {
            if (Math.abs(p1.x - p2.x) < 1e-9) return Math.max(0, p1.y);
            const t = (omega - p1.x) / (p2.x - p1.x);
            return Math.max(0, p1.y + t * (p2.y - p1.y));
        }
    }
    return Math.max(0, lastPoint.y);
}


// --- Chart Initialization ---
function initCharts() { // MODIFIED: Removed tension from S(omega) dataset
    const AXIS_TITLE_FONT_SIZE = 13, TICK_LABEL_FONT_SIZE = 11, LEGEND_LABEL_FONT_SIZE = 12;
    const baseChartOptions = {
        animation: { duration: 0 }, responsive: true, maintainAspectRatio: false,
        scales: {
            x: { type: 'linear', position: 'bottom', title: { display: true, font: { size: AXIS_TITLE_FONT_SIZE } }, ticks: { font: { size: TICK_LABEL_FONT_SIZE } } },
            y: { title: { display: true, font: { size: AXIS_TITLE_FONT_SIZE } }, ticks: { font: { size: TICK_LABEL_FONT_SIZE } } }
        },
        elements: { line: { tension: 0.0 }, point: { radius: 0 } }, // Default to no tension globally
        plugins: { legend: { labels: { font: { size: LEGEND_LABEL_FONT_SIZE } } }, tooltip: { enabled: false } }
    };

    charts.spectrumAndPhase = new Chart(canvasElements.spectrumAndPhase.getContext('2d'), {
        type: 'line', data: {
            datasets: [
                { label: 'S(ω)', data: [], borderColor: '#3498db', borderWidth: 2, fill: false, yAxisID: 'yS' /* tension: 0.4 REMOVED */ },
                { label: 'Φ(ω)', data: [], borderColor: '#2ecc71', borderWidth: 2, fill: false, yAxisID: 'yPhi', tension: 0.1 }, // Phase can still have slight visual tension
                { data: [], borderColor: '#e74c3c', backgroundColor: '#e74c3c', borderWidth: 1, fill: false, yAxisID: 'yS', showLine: false, pointRadius: 5, pointHoverRadius: 7, hidden: true, order: -1 }
            ]
        }, options: {
            ...baseChartOptions,
            scales: {
                x: { ...baseChartOptions.scales.x, min: OMEGA_AXIS_MIN_DEFAULT, max: OMEGA_AXIS_MAX_DEFAULT, title: { ...baseChartOptions.scales.x.title, text: 'ω (Frequency)' } },
                yS: { type: 'linear', position: 'left', title: { display: true, text: 'Amplitude S(ω)', font: { size: AXIS_TITLE_FONT_SIZE } }, min: 0, max: 1.1, grid: { drawOnChartArea: true }, ticks: { font: { size: TICK_LABEL_FONT_SIZE } } },
                yPhi: { type: 'linear', position: 'right', title: { display: true, text: 'Phase Φ(ω) (rad)', font: { size: AXIS_TITLE_FONT_SIZE } }, grid: { drawOnChartArea: false }, ticks: { font: { size: TICK_LABEL_FONT_SIZE } } }
            },
            plugins: { ...baseChartOptions.plugins, legend: { ...baseChartOptions.plugins.legend, onClick: null, labels: { ...baseChartOptions.plugins.legend.labels, filter: function(legendItem) { return legendItem.datasetIndex !== 2; } } } }
        }
    });
    charts.individualCosines = new Chart(canvasElements.individualCosines.getContext('2d'), { type: 'line', data: { datasets: [] }, options: { ...baseChartOptions, elements: {...baseChartOptions.elements, line: {tension: 0.1}}, plugins: { ...baseChartOptions.plugins, legend: { ...baseChartOptions.plugins.legend, display: true, position: 'top', onClick: Chart.defaults.plugins.legend.onClick } }, scales: { x: { ...baseChartOptions.scales.x, title: { ...baseChartOptions.scales.x.title, text: 't (Time)' } }, y: { ...baseChartOptions.scales.y, title: { ...baseChartOptions.scales.y.title, text: 'Amplitude (Offset)' } } } } });
    charts.sumCosines = new Chart(canvasElements.sumCosines.getContext('2d'), { type: 'line', data: { datasets: [{ label: 'Summed Pulse E(t)', data: [], borderColor: '#e74c3c', borderWidth: 2, fill: false, tension: 0.1 }, { label: 'E(t) Envelope', data: [], borderColor: '#f1a7a0', borderWidth: 1.5, borderDash: [5, 5], fill: false, tension: 0.1 }] }, options: { ...baseChartOptions, plugins: { ...baseChartOptions.plugins, legend: { ...baseChartOptions.plugins.legend, onClick: Chart.defaults.plugins.legend.onClick } }, scales: { x: { ...baseChartOptions.scales.x, title: { ...baseChartOptions.scales.x.title, text: 't (Time)' } }, y: { ...baseChartOptions.scales.y, title: { ...baseChartOptions.scales.y.title, text: 'Amplitude E(t)' } } } } });
    charts.sumIntensity = new Chart(canvasElements.sumIntensity.getContext('2d'), { type: 'line', data: { datasets: [{ label: 'Intensity I(t)', data: [], borderColor: '#9b59b6', borderWidth: 2, fill: false, tension: 0.1 }, { label: 'I(t) Envelope', data: [], borderColor: '#cda5d8', borderWidth: 1.5, borderDash: [5, 5], fill: false, tension: 0.1 }] }, options: { ...baseChartOptions, plugins: { ...baseChartOptions.plugins, legend: { ...baseChartOptions.plugins.legend, onClick: Chart.defaults.plugins.legend.onClick } }, scales: { x: { ...baseChartOptions.scales.x, title: { ...baseChartOptions.scales.x.title, text: 't (Time)' } }, y: { ...baseChartOptions.scales.y, title: { ...baseChartOptions.scales.y.title, text: 'Intensity I(t) (a.u.)' }, min: 0 } } } } )};


// --- UI Control Logic --- (No changes)
function updateOmega0ControlsCustomMode() {
    const isCustom = currentSpectrumMode === 'custom';
    const isAuto = autoOmega0CustomCheckbox.checked;
    autoOmega0CustomCheckbox.style.display = isCustom ? 'inline-block' : 'none';
    const autoOmega0Label = document.querySelector('label[for="autoOmega0Custom"]');
    if (autoOmega0Label) autoOmega0Label.style.display = isCustom ? 'inline-block' : 'none';
    if (isCustom && isAuto) {
        sliders.omega0.disabled = true;
        sliderVals.omega0_val.disabled = true;
        omega0SliderGroup.classList.add('disabled-look');
    } else {
        sliders.omega0.disabled = false;
        sliderVals.omega0_val.disabled = false;
        omega0SliderGroup.classList.remove('disabled-look');
    }
    document.getElementById('omega0_custom_note').style.display = isCustom ? 'block' : 'none';
}

// --- Event Listeners --- (No changes)
for (const key in sliders) {
    if (sliders.hasOwnProperty(key) && sliders[key]) {
        sliders[key].addEventListener('input', () => {
            if (key === 'omega0' && currentSpectrumMode === 'custom' && !autoOmega0CustomCheckbox.checked) {
                 sliderVals.omega0_val.value = sliders.omega0.value;
            }
            updatePlots();
        });
    }
}
function setupNumberInputListener(sliderKey, precision) {
    const numberInput = sliderVals[sliderKey + "_val"], rangeSlider = sliders[sliderKey];
    if (!numberInput || !rangeSlider) return;
    numberInput.addEventListener('change', () => {
        let value = parseFloat(numberInput.value);
        const min = parseFloat(rangeSlider.min), max = parseFloat(rangeSlider.max);
        if (isNaN(value)) value = parseFloat(rangeSlider.value); else value = Math.max(min, Math.min(max, value));
        rangeSlider.value = value;
        updatePlots();
    });
}
showPeakConnectorCheckbox.addEventListener('change', updatePlots);
showSpectralPhaseCheckbox.addEventListener('change', updatePlots);
showEnvelopeECheckbox.addEventListener('change', updatePlots);
showEnvelopeICheckbox.addEventListener('change', updatePlots);
gaussianModeRadio.addEventListener('change', () => {
    if (gaussianModeRadio.checked) {
        currentSpectrumMode = 'gaussian';
        gaussianControlsDiv.style.display = 'block';
        resetCustomSpectrumButton.style.display = 'none';
        customSpectrumHelpText.style.display = 'none';
        if (charts.spectrumAndPhase) charts.spectrumAndPhase.data.datasets[2].hidden = true;
        canvasElements.spectrumAndPhase.style.cursor = 'default';
        updateOmega0ControlsCustomMode();
        updatePlots();
    }
});
customModeRadio.addEventListener('change', () => {
    if (customModeRadio.checked) {
        isInitialSwitchToCustom = (currentSpectrumMode !== 'custom');
        currentSpectrumMode = 'custom';
        gaussianControlsDiv.style.display = 'none';
        resetCustomSpectrumButton.style.display = 'block';
        customSpectrumHelpText.style.display = 'block';
        initializeOrUpdateCustomSpectrumPoints(isInitialSwitchToCustom);
        if (charts.spectrumAndPhase) charts.spectrumAndPhase.data.datasets[2].hidden = false;
        isInitialSwitchToCustom = false;
        updateOmega0ControlsCustomMode();
        updatePlots();
    }
});
autoOmega0CustomCheckbox.addEventListener('change', () => {
    updateOmega0ControlsCustomMode();
    updatePlots();
});
resetCustomSpectrumButton.addEventListener('click', () => {
    if (currentSpectrumMode === 'custom') {
        customSpectrumDataPoints.forEach(point => { point.y = 0.0; });
        updatePlots();
    }
});
let draggedPointIndex = -1;
let isDragging = false;
let dragStartYPixelOffset = 0;
canvasElements.spectrumAndPhase.addEventListener('mousedown', (evt) => {
    if (currentSpectrumMode !== 'custom' || !charts.spectrumAndPhase) return;
    const chart = charts.spectrumAndPhase;
    const xScale = chart.scales.x, yScale = chart.scales.yS;
    const rect = canvasElements.spectrumAndPhase.getBoundingClientRect();
    const xPixel = evt.clientX - rect.left, yPixel = evt.clientY - rect.top;
    let newDraggedPointIndex = -1, minHorizontalPixelDiff = Infinity;
    const xClickPixelTolerance = 15;
    customSpectrumDataPoints.forEach((point, index) => {
        const pointXPixel = xScale.getPixelForValue(point.x);
        const horizontalPixelDiff = Math.abs(xPixel - pointXPixel);
        if (horizontalPixelDiff <= xClickPixelTolerance && horizontalPixelDiff < minHorizontalPixelDiff) {
            minHorizontalPixelDiff = horizontalPixelDiff; newDraggedPointIndex = index;
        }
    });
    if (newDraggedPointIndex !== -1) {
        if (Math.abs(customSpectrumDataPoints[newDraggedPointIndex].x) < 1e-9) return;
        draggedPointIndex = newDraggedPointIndex; isDragging = true;
        chart.options.animation = false; canvasElements.spectrumAndPhase.style.cursor = 'grabbing';
        const draggedPointYValue = customSpectrumDataPoints[draggedPointIndex].y;
        const draggedPointYPixel = yScale.getPixelForValue(draggedPointYValue);
        dragStartYPixelOffset = yPixel - draggedPointYPixel;
    }
});
canvasElements.spectrumAndPhase.addEventListener('mousemove', (evt) => {
    if (!charts.spectrumAndPhase) return;
    const chart = charts.spectrumAndPhase;
    if (currentSpectrumMode !== 'custom') { canvasElements.spectrumAndPhase.style.cursor = 'default'; return; }
    if (isDragging && draggedPointIndex !== -1) {
        const rect = canvasElements.spectrumAndPhase.getBoundingClientRect();
        const currentMouseYPixel = evt.clientY - rect.top;
        const yScale = chart.scales.yS;
        const targetPointYPixel = currentMouseYPixel - dragStartYPixelOffset;
        let newYValue = yScale.getValueForPixel(targetPointYPixel);
        newYValue = Math.max(yScale.min, Math.min(yScale.max, newYValue));
        customSpectrumDataPoints[draggedPointIndex].y = newYValue;
        const spectrumDataLine = omegaPlotArray.map(om => ({ x: om, y: getCustomSpectrumValue(om) })); // Recalculate with smooth getCustomSpectrumValue
        chart.data.datasets[0].data = spectrumDataLine; chart.data.datasets[2].data = customSpectrumDataPoints;
        chart.update('none');
    } else {
        const xScale = chart.scales.x, yScale = chart.scales.yS;
        const rect = canvasElements.spectrumAndPhase.getBoundingClientRect();
        const xPixel = evt.clientX - rect.left, yPixel = evt.clientY - rect.top;
        let onDraggablePoint = false; const hoverPixelTolerance = 10;
        customSpectrumDataPoints.forEach((point) => {
            if (Math.abs(point.x) < 1e-9) return;
            const pointXPixel = xScale.getPixelForValue(point.x), pointYPixel = yScale.getPixelForValue(point.y);
            if (Math.abs(xPixel - pointXPixel) <= hoverPixelTolerance && Math.abs(yPixel - pointYPixel) <= hoverPixelTolerance + 5) onDraggablePoint = true;
        });
        canvasElements.spectrumAndPhase.style.cursor = onDraggablePoint ? 'grab' : 'default';
    }
});
let lastMouseX = 0, lastMouseY = 0;
document.addEventListener('mousemove', (e) => { lastMouseX = e.clientX; lastMouseY = e.clientY; });
function handleDragEnd() {
    if (isDragging) {
        isDragging = false; draggedPointIndex = -1; dragStartYPixelOffset = 0;
        if (charts.spectrumAndPhase) charts.spectrumAndPhase.options.animation = { duration: 0 };
        const ev = new MouseEvent('mousemove', { clientX: lastMouseX, clientY: lastMouseY, bubbles: true, cancelable: true });
        canvasElements.spectrumAndPhase.dispatchEvent(ev);
        updatePlots();
    }
}
canvasElements.spectrumAndPhase.addEventListener('mouseup', handleDragEnd);
canvasElements.spectrumAndPhase.addEventListener('mouseout', (evt) => {
    if (isDragging && (evt.relatedTarget === null || evt.relatedTarget.nodeName === 'HTML')) handleDragEnd();
    else if (currentSpectrumMode === 'custom' && !isDragging) canvasElements.spectrumAndPhase.style.cursor = 'default';
});

// --- Initial Setup --- (No changes)
window.onload = () => {
    OMEGA_AXIS_MAX_DEFAULT = parseFloat(sliders.omega0.max) + parseFloat(sliders.delta_omega.max) + 2;
    omegaPlotArray = linspace(OMEGA_AXIS_MIN_DEFAULT, OMEGA_AXIS_MAX_DEFAULT, NUM_OMEGA_PLOT_POINTS);
    initCharts();
    Object.keys(sliders).forEach(key => {
        if (sliders[key]) setupNumberInputListener(key, key.startsWith('phi') ? 2 : 1);
    });
    updateOmega0ControlsCustomMode();
    initializeOrUpdateCustomSpectrumPoints(true);
    updatePlots();
};
