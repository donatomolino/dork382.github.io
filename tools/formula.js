function ANFC() // Aggregate Noise Figure calculator
{
    F1 = Number(document.calculator.number1.value);
    G1 = Number(document.calculator.number2.value);
    F2 = Number(document.calculator.number3.value);
    G2 = Number(document.calculator.number4.value);
    F3 = Number(document.calculator.number5.value);
    G3 = Number(document.calculator.number6.value);
    F4 = Number(document.calculator.number7.value);
    G4 = Number(document.calculator.number8.value);
    F1 = Math.pow(10, (F1 / 10));
    G1 = Math.pow(10, (G1 / 10));
    F2 = Math.pow(10, (F2 / 10));
    G2 = Math.pow(10, (G2 / 10));
    F3 = Math.pow(10, (F3 / 10));
    G3 = Math.pow(10, (G3 / 10));
    F4 = Math.pow(10, (F4 / 10));
    G4 = Math.pow(10, (G4 / 10));
    OA_NF = F1 + ((F2 - 1) / G1) + ((F3 - 1) / (G1 * G2)) + ((F4 - 1) / (G1 * G2 * G3));
    noise_figure_dB = 10 * Math.log10(OA_NF);
    document.calculator.NF.value = OA_NF;
    document.calculator.NF_dB.value = noise_figure_dB;
}

function CWCF() {
    a = Number(document.calculatorCWCF.number1.value);
    c1 = ((1.8412 * 3 * Math.pow(10, 8)) / (2 * Math.PI * a));
    c = (c1 / Math.pow(10, 9));
    document.calculatorCWCF.wavelength.value = c;
}

function DTD() {
    a = Number(document.calculatorDTD.number1.value);
    c = a - 30;
    document.calculatorDTD.dBW.value = c;
}

function DTW() {
    a = Number(document.calculatorDTW.number1.value);
    c = Math.log10(a * 1000) * 10;
    document.calculatorDTW.PdBm.value = c;
}

function WDT() {
    a1 = Number(document.calculatorWDT.number2.value);
    c1 = a1 / 10;
    c2 = Math.pow(10, c1) / 1000;
    document.calculatorWDT.PWatt.value = c2;
}


function MLIC() {
    Er = Number(document.calculatorMLIC.number1.value);
    H = Number(document.calculatorMLIC.number2.value);
    W = Number(document.calculatorMLIC.number3.value);
    if ((W / H) < 1) {
        m1 = W / H;
        P1 = Math.pow((1 + (12 / m1)), -0.5) + (0.04 * (Math.pow((1 - m1), 2)));
        c1 = ((Er + 1) / 2) + (((Er - 1) / 2) * P1);
        P2 = Math.log((8 / m1) + (0.25 * m1));
        c2 = ((60 / (Math.pow(c1, 0.5))) * P2);
        //c1=10;
        //c2=20;
    }
    else {
        m1 = W / H;
        P1 = Math.pow((1 + (12 / m1)), -0.5);
        c1 = ((Er + 1) / 2) + (((Er - 1) / 2) * P1);
        P2 = m1 + 1.393 + (0.666 * (Math.log(m1 + 1.444)));
        c2 = ((120 * Math.PI) / (Math.pow(c1, 0.5) * P2));
        //c1=10;
        //c2=20;
    }
    document.calculatorMLIC.effperm.value = c1;
    document.calculatorMLIC.impedance.value = c2;
}

function PAC() {
    z = Number(document.calculatorPAC.number1.value);
    att_dB = Number(document.calculatorPAC.number2.value);
    a = att_dB / 10;
    k = Math.pow(10, a);
    c1 = z * ((k - 1) / (k + 1));
    c2 = 2 * z * ((k) / (Math.pow(k, 2) - 1));
    document.calculatorPAC.R12_val.value = c1;
    document.calculatorPAC.R3_val.value = c2;
}

function PDC() {
    a = Number(document.calculatorPDC.number1.value);
    c = -10 * Math.log10(1 / a);
    document.calculatorPDC.PATHLOSS.value = c;
}

function RWC() {
    a = Number(document.calculatorRWC.number1.value);
    b = Number(document.calculatorRWC.number2.value);
    m = Number(document.calculatorRWC.number3.value);
    n = Number(document.calculatorRWC.number4.value);
    Er = Number(document.calculatorRWC.number5.value);
    Mu = Number(document.calculatorRWC.number6.value);
    pi_val = Math.PI;
    p1 = Math.pow(((m * pi_val) / a), 2);
    p2 = Math.pow(((n * pi_val) / b), 2);
    p = p1 + p2;
    c1 = Math.pow(p, 0.5);
    c2 = Math.pow((Mu * Er), 0.5);
    c = (1 / ((2 * Math.PI * c2))) * c1;
    c = c / Math.pow(10, 9);
    document.calculatorRWC.frequency.value = c;
}

function RLTVC() {
    a = Number(document.calculatorRLTVC.number1.value);
    x1 = Math.pow(10, (a / 20));
    c = ((1 + x1) / (x1 - 1));
    document.calculatorRLTVC.VSWR.value = c;
}

function SDC() {
    Freq = Number(document.calculatorSDC.number1.value);
    sigma = Number(document.calculatorSDC.number2.value);
    c1 = Math.PI * 4 * Math.PI * Math.pow(10, -7) * sigma * Freq * Math.pow(10, 9);
    c = (1 / Math.pow(c1, 0.5));
    document.calculatorSDC.skindepth.value = c;
}

function RFE() {
    E = Number(document.calculatorRFE.number1.value);
    sigma = Number(document.calculatorRFE.number2.value);
    md = Number(document.calculatorRFE.number3.value);
    c1 = ((sigma * (Math.pow(E, 2))) / md);
    c2 = ((Math.pow(E, 2)) / 377);
    document.calculatorRFE.SAR.value = c1;
    document.calculatorRFE.PD.value = c2;
}

function SIC() {
    W = Number(document.calculatorSIC.number1.value);
    h = Number(document.calculatorSIC.number2.value);
    Er = Number(document.calculatorSIC.number3.value);
    if ((W / h) > 0.35) {
        m1 = W / h;
        c1 = (((30 * Math.PI) / Math.sqrt(Er)) * (1 / (m1 + 0.441)));
        c1 = Math.round(c1);
    }
    else {
        m1 = ((W / h) - (0.35 - (Math.pow(m1, 2))));
        c1 = (((30 * Math.PI) / Math.sqrt(Er)) * (1 / (m1 + 0.441)));
        c1 = Math.round(c1);
    }
    document.calculatorSIC.striplineimpedance.value = c1;
}

function TAC() {
    z = Number(document.calculatorTAC.number1.value);
    att_dB = Number(document.calculatorTAC.number2.value);
    a = att_dB / 10;
    k = Math.pow(10, a);
    c1 = z * ((k + 1) / (k - 1));
    c2 = z * ((Math.pow(k, 2) - 1) / (2 * k));
    document.calculatorTAC.R13_val.value = c1;
    document.calculatorTAC.R2_val.value = c2;
}

function RFF() {
    ht = Number(document.calculatorRFF.number1.value);
    hr = Number(document.calculatorRFF.number2.value);
    current = Number(document.calculatorRFF.number3.value);
    d = Number(document.calculatorRFF.number4.value);
    f = Number(document.calculatorRFF.number5.value);
    lambda = ((3 * Math.pow(10, 8)) / (f * Math.pow(10, 9)));
    pi_val = Math.PI;
    c = ((120 * pi_val * ht * hr * current) / (lambda * d));
    document.calculatorRFF.fieldstrength.value = c;
}

function RFT() {
    zi = Number(document.calculatorRFT.number1.value);
    zl = Number(document.calculatorRFT.number2.value);
    c1 = zi / zl;
    c2 = Math.pow(c1, 0.5);
    c2 = 1 / c2;
    document.calculatorRFT.NsbyNp.value = c2;
}

function TEM() {
    Fr = Number(document.calculatorTEM.number1.value);
    Er = Number(document.calculatorTEM.number2.value);
    c = (300 / (Fr * (Math.pow(Er, 0.5))));
    document.calculatorTEM.wavelength.value = c;
}

function RWB() {
    a = Number(document.calculatorRWB.number1.value);
    b = Number(document.calculatorRWB.number2.value);
    f = Number(document.calculatorRWB.number3.value);
    lambda = (((3 * Math.pow(10, 8)) / (f * Math.pow(10, 9))) * 100);
    m2 = (lambda / (2 * a));
    m1 = (1 - (Math.pow(m2, 2)));
    c = (597 * a * b * (Math.pow(m1, 0.5)));
    document.calculatorRWB.power.value = c;
}

function PR() {
    len = Number(document.calculatorPR.number1.value);
    width = Number(document.calculatorPR.number2.value);
    thick = Number(document.calculatorPR.number3.value);
    res = Number(document.calculatorPR.number4.value);
    len = len * Math.pow(10, -3);
    width = width * Math.pow(10, -3);
    thick = thick * Math.pow(10, -6);
    c = ((len * res) / (width * thick))
    document.calculatorPR.resistance.value = c;
}

function CSI() {
    N = Number(document.calculatorCSI.number1.value);
    S = Number(document.calculatorCSI.number2.value);
    Di = Number(document.calculatorCSI.number3.value);
    Do = Number(document.calculatorCSI.number4.value);
    c = (0.03125 * (Math.pow(N, 2)) * (Math.pow(Do, 2)));
    c = (c / (Math.pow(10, 6)));
    document.calculatorCSI.inductance.value = c;
}

function IC() {
    Er = Number(document.calculatorIC.number1.value);
    n = Number(document.calculatorIC.number2.value);
    len = Number(document.calculatorIC.number3.value);
    W = Number(document.calculatorIC.number4.value);
    c1 = ((Er + 1) / W);
    c = (c1 * len * (((n - 3) * 0.089) + 0.10));
    document.calculatorIC.capacitance.value = c;
}

function RK() {
    F = Number(document.calculatorRK.number1.value);
    Vo = Number(document.calculatorRK.number2.value);
    Io = Number(document.calculatorRK.number3.value);
    N = Number(document.calculatorRK.number4.value);
    L = Number(document.calculatorRK.number5.value);
    c1 = ((0.398 * Vo * Io * (Math.pow(10, -3))) / N);
    c2 = (6.74 * (Math.pow(10, -6)) * (F * Math.pow(10, 9)) * (L) * ((Math.pow(Vo, 0.5)) / N)) - Vo;
    c3 = ((0.398 / N) * 100);
    document.calculatorRK.Power.value = c1;
    document.calculatorRK.repellervoltage.value = c2;
    document.calculatorRK.efficiency.value = c3;
}

function MC() {
    a = Number(document.calculatorMC.number1.value);
    b = Number(document.calculatorMC.number2.value);
    B = Number(document.calculatorMC.number3.value);
    Fr = Number(document.calculatorMC.number4.value);
    N = Number(document.calculatorMC.number5.value);
    m = (9.11 * (Math.pow(10, -31)));
    e = (1.6 * (Math.pow(10, -19)));
    m1 = (1 - ((Math.pow(a, 2)) / (Math.pow(b, 2))));
    c1 = (((e * (Math.pow(b, 2)) * (Math.pow(B, 2))) / (8 * m)) * (Math.pow(m1, 2)));
    c2 = ((2 * (Math.PI) * Fr * (Math.pow(10, 9)) * B) / N) * ((Math.pow(b, 2)) - (Math.pow(a, 2)));
    c1 = c1 / 1000;
    c2 = c2 / 1000;
    document.calculatorMC.Hull_potential.value = c1;
    document.calculatorMC.Hartree_potential.value = c2;
}

function RCR() {
    b = Number(document.calculatorRCR.number1.value);
    a = Number(document.calculatorRCR.number2.value);
    d = Number(document.calculatorRCR.number3.value);
    sigma = Number(document.calculatorRCR.number4.value);
    m = 1;
    n = 0;
    p = 1;
    m1 = Math.pow((m / a), 2) + Math.pow((n / b), 2) + Math.pow((p / d), 2);
    m2 = Math.pow(m1, 0.5);
    c1 = (((3 * (Math.pow(10, 10))) / 2) * m2);
    c1 = (c1 / (Math.pow(10, 9)));
    Rs1 = (((Math.PI) * c1 * (Math.pow(10, 9)) * 4 * (Math.PI) * (Math.pow(10, -7))) / sigma);
    Rs = Math.pow(Rs1, 0.5);
    m3 = 1 + (a / (2 * b));
    c2 = (((1.11 * 120 * (Math.PI)) / Rs) * (1 / m3));
    c3 = ((c1 * (Math.pow(10, 9))) / c2);
    document.calculatorRCR.frequency.value = c1;
    document.calculatorRCR.Q.value = c2;
    document.calculatorRCR.bandwidth.value = c3;
}

function CCR() {
    a = Number(document.calculatorCCR.number1.value);
    d = Number(document.calculatorCCR.number2.value);
    sigma = Number(document.calculatorCCR.number3.value);
    n = 0;
    p = 1;
    q = 0;
    m = 2.405;
    m1 = Math.pow((m / a), 2)
    m2 = Math.pow(m1, 0.5);
    c1 = (((3 * (Math.pow(10, 10))) / 2) * m2);
    c1 = (c1 / (Math.pow(10, 9)));
    Rs1 = (((Math.PI) * c1 * (Math.pow(10, 9)) * 4 * (Math.PI) * (Math.pow(10, -7))) / sigma);
    Rs = Math.pow(Rs1, 0.5);
    m3 = (1 + (a / d));
    c2 = (((1.202 * 377) / Rs) * (1 / m3));
    c3 = ((c1 * (Math.pow(10, 9))) / c2);
    document.calculatorCCR.frequency.value = c1;
    document.calculatorCCR.resistance.value = Rs;
    document.calculatorCCR.Q.value = c2;
    document.calculatorCCR.bandwidth.value = c3;
}

function RCC() {
    ZL = Number(document.calculatorRCC.number1.value);
    ZS = Number(document.calculatorRCC.number2.value);
    VI = Number(document.calculatorRCC.number3.value);
    c1 = ((ZL - ZS) / (ZL + ZS));
    c2 = (c1 * VI);
    document.calculatorRCC.reflection_coefficient.value = c1;
    document.calculatorRCC.reflected_voltage.value = c2;
}

function PWC() {
    Er = Number(document.calculatorPWC.number1.value);
    Mr = Number(document.calculatorPWC.number2.value);
    F = Number(document.calculatorPWC.number3.value);
    m1 = Er * Mr;
    c1 = ((3 * (Math.pow(10, 8))) / (Math.pow(m1, 0.5)));
    c2 = (c1 / (F * (Math.pow(10, 9))));
    m2 = Mr / Er;
    c3 = (377 * (Math.pow(m2, 0.5)));
    document.calculatorPWC.phasevelocity.value = c1;
    document.calculatorPWC.wavelength.value = c2;
    document.calculatorPWC.waveimpedance.value = c3;
}

function SC() {
    Er = Number(document.calculatorSC.number1.value);
    Fr = Number(document.calculatorSC.number2.value);
    c1 = ((Er + 1) / 2);
    lambda = ((3 * (Math.pow(10, 8))) / (Fr * (Math.pow(10, 9))));
    c2 = (lambda / (Math.pow(c1, 0.5)));
    document.calculatorSC.effective_dielectric_constant.value = c1;
    document.calculatorSC.guide_wavelength.value = c2;
}

function TDR() {
    Ve = Number(document.calculatorTDR.number1.value);
    Td = Number(document.calculatorTDR.number2.value);
    c1 = (((3 * Math.pow(10, 8)) * Ve * Td));
    document.calculatorTDR.TDR_length.value = c1;
}

function TWT() {
    F = Number(document.calculatorTWT.number1.value);
    K = Number(document.calculatorTWT.number2.value);
    V = Number(document.calculatorTWT.number3.value);
    I = Number(document.calculatorTWT.number4.value);
    Vo = Number(document.calculatorTWT.number5.value);
    L = Number(document.calculatorTWT.number6.value);
    Pi = 3.14159;
    F = (F * (Math.pow(10, 9)));
    I = I * Math.pow(10, -3);
    V = V * Math.pow(10, 3);
    m1 = ((I * K) / (4 * V));
    m2 = Math.pow(m1, 0.3333);
    m3 = (L * 2 * Pi * F) / (2 * Pi * Vo);
    c1 = 47.3 * m3 * m2;
    c2 = -9.45 + c1;
    document.calculatorTWT.TWT_Gain.value = c2;
}

function VDC() {
    C = Number(document.calculatorVDC.number1.value);
    Vb = Number(document.calculatorVDC.number2.value);
    V = Number(document.calculatorVDC.number3.value);
    Rs = Number(document.calculatorVDC.number4.value);
    F = Number(document.calculatorVDC.number5.value);
    m = Number(document.calculatorVDC.number6.value);
    K = 1;
    m1 = (Vb - V);
    c1 = ((C * K) / (Math.pow(m1, m)));
    //lambda=( (3* (Math.pow(10,8)))/(Fr*(Math.pow(10,9))) );
    c2 = (1 / (2 * (Math.PI) * Rs * c1 * (Math.pow(10, -12))));
    c2 = (c2 / (Math.pow(10, 6)));
    c3 = (c2 / F);
    document.calculatorVDC.diode_capacitance.value = c1;
    document.calculatorVDC.diode_cutoff_frequency.value = c2;
    document.calculatorVDC.quality_factor.value = c3;
}

function TDC() {
    R = Number(document.calculatorTDC.number1.value);
    Cj = Number(document.calculatorTDC.number2.value);
    Rs = Number(document.calculatorTDC.number3.value);
    Ls = Number(document.calculatorTDC.number4.value);
    Cj = (Cj * (Math.pow(10, -12)));
    Ls = (Ls * (Math.pow(10, -9)));
    m1 = ((R / Rs) - 1);
    c1 = ((1 / (2 * (Math.PI) * R * Cj)) * (Math.pow(m1, 0.5)));
    m2 = R * Cj;
    m3 = ((1 / (Ls * Cj)) - (1 / (Math.pow(m2, 2))));
    c2 = ((1 / (2 * 3.14)) * (Math.pow(m3, 0.5)));
    c1 = (c1 / (Math.pow(10, 6)));
    document.calculatorTDC.resistive_frequency.value = c1;
    c2 = (c2 / (Math.pow(10, 6)));
    document.calculatorTDC.resonant_frequency.value = c2;
}

function DCC() {
    P1 = Number(document.calculatorDCC.number1.value);
    P2 = Number(document.calculatorDCC.number2.value);
    P3 = Number(document.calculatorDCC.number3.value);
    P4 = Number(document.calculatorDCC.number4.value);
    c1 = (-10 * (Math.log10(P4 / P1)));
    c2 = (-10 * (Math.log10(1 - (P4 / P1))));
    c3 = (-10 * (Math.log10((P2 / P1) + (P3 / P1) + (P4 / P1))));
    c4 = (-10 * (Math.log10((P3 / P4))));
    document.calculatorDCC.coupling.value = c1;
    document.calculatorDCC.coupling_loss.value = c2;
    document.calculatorDCC.Insertion_loss.value = c3;
    document.calculatorDCC.directivity.value = c4;
}

function PPM() {
    Fr = Number(document.calculatorPPM.number1.value);
    PPM = Number(document.calculatorPPM.number2.value);
    c1 = (((Fr * (Math.pow(10, 9))) * PPM) / (Math.pow(10, 6)));
    document.calculatorPPM.freq.value = c1;
}

function MDSS() {
    B = Number(document.calculatorMDS.number1.value);
    NF = Number(document.calculatorMDS.number2.value);
    // -174 is the noise floor in dBm/Hz at 290 kelvin and with 1Hz bandwidth
    c1 = (-174 + (10 * Math.log10(B)) + NF);
    document.calculatorMDS.MDS.value = c1;
}

function RCS() {
    P = Number(document.calculatorRCS.number1.value);
    c = ((1 + P) / (1 - P));
    document.calculatorRCS.SWR.value = c;
}

function SFDRR() {
    IIP3 = Number(document.calculatorSFDR.number1.value);
    MDS1 = Number(document.calculatorSFDR.number2.value);
    c1 = ((2 / 3) * (IIP3 - MDS1));
    document.calculatorSFDR.SFDR.value = c1;
}

function RAT() {
    alpha = Number(document.calculatorRAT.number1.value);
    a1 = - (alpha / 4.34);
    c1 = (290 * (1 - Math.exp(a1)));
    document.calculatorRAT.rain_absorption_temperature.value = c1;
}

function MM() {
    N = Number(document.calculatorMM.number1.value);
    T = Number(document.calculatorMM.number2.value);
    F = Number(document.calculatorMM.number3.value);
    c1 = T / F;
    c2 = T / N;
    document.calculatorMM.MTBF.value = c1;
    document.calculatorMM.MTTF.value = c2;
}

function AEZ() {
    LAE = Number(document.calculatorAEZ.number1.value);
    LOE = Number(document.calculatorAEZ.number2.value);
    LOS = Number(document.calculatorAEZ.number3.value);
    L = LOE - LOS;
    L = -(L);
    m1 = (Math.cos(LAE * Math.PI / 180) * Math.cos(L * Math.PI / 180)) - 0.151;
    m2 = 1 - (Math.cos(LAE * Math.PI / 180) * Math.cos(LAE * Math.PI / 180) * Math.cos(L * Math.PI / 180) * Math.cos(L * Math.PI / 180));
    c1 = Math.atan(m1 / Math.pow(m2, 0.5));
    c1 = c1 * (180 / Math.PI);
    m3 = Math.tan(L * Math.PI / 180) / Math.sin(LAE * Math.PI / 180);
    c2 = Math.atan(m3);
    c2 = c2 * 180 / Math.PI;
    if ((LAE < 0) && (LOE < 0) && (L < 0)) {
        c3 = c2;
    }
    else if ((LAE < 0) && (L > 0)) {
        c3 = 360 - c2;
    }
    else if ((LAE > 0) && (L < 0)) {
        c3 = 180 + c2;
    }
    else if ((LAE > 0) && (L > 0)) {
        c3 = 180 - c2;
    }
    else {
        c3 = 0;
    }
    document.calculatorAEZ.elevation.value = c1;
    document.calculatorAEZ.azimuth.value = c3;
}

function ATP() {
    TOI = Number(document.calculatorATP.number1.value);
    PIN = Number(document.calculatorATP.number2.value);
    G = Number(document.calculatorATP.number3.value);
    m1 = (-((TOI - PIN - G) / 10));
    c1 = (13.2 * (Math.pow(10, m1)));
    document.calculatorATP.AM_to_PM_conversion.value = c1;
}

function PUC() {
    a = Number(document.calculatorPUC.number1.value);
    c = 10 * Math.log10(a * 1e6);
    document.calculatorPUC.PdBmicroWatt.value = c;
}

function RFP() {
    f = Number(document.calculatorRFP.number1.value);
    PN = Number(document.calculatorRFP.number2.value);
    p1 = PN / 10;
    p2 = 2 * Math.pow(10, p1);
    c1 = Math.pow(p2, 0.5);
    c2 = c1 / (2 * Math.PI * f * 1e6);
    document.calculatorRFP.phasejitter.value = c1;
    document.calculatorRFP.jitter.value = c2;
}

function ACPRR() {
    Pmax = Number(document.calculatorACPR.number1.value);
    Pavg = Number(document.calculatorACPR.number2.value);
    Pin = Number(document.calculatorACPR.number3.value);
    IP3 = Number(document.calculatorACPR.number4.value);
    Crestfactor = Pmax / Pavg;
    c1 = -20.75 + (1.6 * Crestfactor) + (2 * (Pin - IP3));
    document.calculatorACPR.ACPR.value = c1;
}

function CIP() {
    G1 = Number(document.calculatorCIP.number1.value);
    IP3_1 = Number(document.calculatorCIP.number2.value);
    G2 = Number(document.calculatorCIP.number3.value);
    IP3_2 = Number(document.calculatorCIP.number4.value);
    G3 = Number(document.calculatorCIP.number5.value);
    IP3_3 = Number(document.calculatorCIP.number6.value);
    G4 = Number(document.calculatorCIP.number7.value);
    IP3_4 = Number(document.calculatorCIP.number8.value);
    G5 = Number(document.calculatorCIP.number9.value);
    IP3_5 = Number(document.calculatorCIP.number10.value);
    G1 = Math.pow(10, (G1 / 10));
    IP3_1 = Math.pow(10, (IP3_1 / 10));
    G2 = Math.pow(10, (G2 / 10));
    IP3_2 = Math.pow(10, (IP3_2 / 10));
    G3 = Math.pow(10, (G3 / 10));
    IP3_3 = Math.pow(10, (IP3_3 / 10));
    G4 = Math.pow(10, (G4 / 10));
    IP3_4 = Math.pow(10, (IP3_4 / 10));
    G5 = Math.pow(10, (G5 / 10));
    IP3_5 = Math.pow(10, (IP3_5 / 10));
    AIP3 = (1 / IP3_1) + (G1 / IP3_2) + ((G1 * G2) / IP3_3) + ((G1 * G2 * G3) / IP3_4) + ((G1 * G2 * G3 * G4) / IP3_5);
    AIP3 = 1 / AIP3;
    AIP3_dBm = 10 * Math.log10(AIP3);
    document.calculatorCIP.cascaded_IP3_mw.value = AIP3;
    document.calculatorCIP.cascaded_IP3_dBm.value = AIP3_dBm;
}

function RFA() {
    Pi = Number(document.calculatorRFA.number1.value);
    Po = Number(document.calculatorRFA.number2.value);
    Pdc = Number(document.calculatorRFA.number3.value);
    Amplifier_PAE = ((Po - Pi) / (Pdc)) * 100;
    document.calculatorRFA.PAE.value = Amplifier_PAE;
}

function TWL() {
    Er = Number(document.calculatorTWL.number1.value);
    Mr = Number(document.calculatorTWL.number2.value);
    s = Number(document.calculatorTWL.number3.value);
    d = Number(document.calculatorTWL.number4.value);
    E0 = 8.854 * Math.pow(10, -12);
    M0 = 4 * Math.PI * Math.pow(10, -7);
    E = E0 * Er;
    Mu = M0 * Mr;
    L = (Mu / Math.PI) * Math.acosh((s / d));
    C = (Math.PI * E) / (Math.acosh(s / d))
    z = Math.pow((L / C), 0.5);
    document.calculatorTWL.twin_wire_line_inductance.value = L;
    document.calculatorTWL.twin_wire_line_capacitance.value = C;
    document.calculatorTWL.twin_wire_line_impedance.value = z;
}

function SINAD() {
    SINAD_in = Number(document.calculatorSINAD.number1.value);
    ENOB_out = (SINAD_in - 1.76) / 6.02;
    document.calculatorSINAD.ENOB.value = ENOB_out;
}

function SNR() {
    SNR_in = Number(document.calculatorSNR.number1.value);
    //SNR_in= Math.pow(10,(SNR_in/10));
    resolution_out = (SNR_in - 1.76) / 6.02;
    document.calculatorSNR.RESOLUTION.value = resolution_out;
}

function QC() {
    Rs = Number(document.calculatorQC.number1.value);
    Cs = Number(document.calculatorQC.number2.value);
    Cp = Number(document.calculatorQC.number3.value);
    Ls = Number(document.calculatorQC.number4.value);
    Cs = Cs * Math.pow(10, -15);
    Cp = Cp * Math.pow(10, -12);
    Ls = Ls * Math.pow(10, -3);
    var1 = Math.pow((Ls * Cs), 0.5);
    var2 = Ls * ((Cs * Cp) / (Cs + Cp))
    c1 = 1 / (2 * Math.PI * var1);
    c2 = 1 / (2 * Math.PI * Math.pow(var2, 0.5));
    c3 = 1 / (2 * Math.PI * c1 * Rs * Cs);
    c4 = Cp / Cs;
    document.calculatorQC.Fs.value = c1;
    document.calculatorQC.Fp.value = c2;
    document.calculatorQC.Q.value = c3;
    document.calculatorQC.pulling_factor.value = c4;
}

function FF() {
    Vrms = Number(document.calculatorFF.number1.value);
    Vdc = Number(document.calculatorFF.number2.value);
    Vpeak = Number(document.calculatorFF.number3.value);
    c1 = Vrms / Vdc;
    c2 = Vpeak / Vrms;
    document.calculatorFF.form_factor.value = c1;
    document.calculatorFF.crest_factor.value = c2;
}

function EL() {
    Freq = Number(document.calculatorEL.number1.value);
    Electrial_Degree = Number(document.calculatorEL.number2.value);
    lambda = 3 * Math.pow(10, 8) / Freq;
    c = (Electrial_Degree / 360) * lambda;
    document.calculatorEL.electrical_length_meters.value = c;
}


function DBU() {
    a = Number(document.calculatorDBU.number1.value);
    fac = a / 20;
    c = Math.pow(10, fac);
    c = c / 1000;
    document.calculatorDBU.mV_per_meter.value = c;
}

function DMS() {
    a1 = Number(document.calculatorDMS.number1.value);
    a2 = Number(document.calculatorDMS.number2.value);
    a3 = Number(document.calculatorDMS.number3.value);
    c = (a1 + (a2 / 60) + (a3 / 3600));
    document.calculatorDMS.decimal_degree.value = c;
}

function ADA() {
    Hb = Number(document.calculatorADA.number1.value);
    Hr = Number(document.calculatorADA.number2.value);
    D = Number(document.calculatorADA.number3.value);
    val = ((Hb - Hr)) / (D * 5280);
    c = Math.atan(val) * (180 / Math.PI);
    document.calculatorADA.downltilt_angle.value = c;
}

function KVA() {
    V1 = Number(document.calculatorKVA.number1.value);
    I1 = Number(document.calculatorKVA.number2.value);
    c1 = V1 * I1 / 1000;
    document.calculatorKVA.single_phase_transformer_KVA.value = c1;
}

function OAS() {
    Fr = Number(document.calculatorOAS.number1.value);
    V = Number(document.calculatorOAS.number2.value);
    c1 = 2 * Math.PI * Fr * V / Math.pow(10, 6);
    document.calculatorOAS.op_amp_slew_rate.value = c1;
}

function RTG() {
    RCS = Number(document.calculatorRAD.number1.value);
    Fr = Number(document.calculatorRAD.number2.value);
    //Lambda= ( (3*Math.pow(10,8))/Fr );
    //Lambda1 = 3*Math.pow(10,8)/Lambda;
    P1 = 10 * Math.log10(RCS);
    P2 = 20 * Math.log10(Fr);
    c1 = P1 + P2 - 38.54;
    document.calculatorRAD.Radar_target_gain_factor.value = c1;
}

function RAD() {
    Radius = Number(document.calculatorRADD.number1.value);
    c1 = Math.PI * Math.pow(Radius, 2);
    document.calculatorRADD.Radar_RCS.value = c1;
}

function WA() {
    Fr = Number(document.calculatorWA.number1.value);
    Len = Number(document.calculatorWA.number2.value);
    Dia = Number(document.calculatorWA.number3.value);
    c1 = (17 * (Len / 12) / ((Math.log10((24 * (Len / 12)) / Dia) - 1) * (1 - (Fr * (Len / 12) / 234) * (Fr * (Len / 12) / 234))));
    re1 = c1;
    c2 = ((1 / ((2 * 3.14159265359 * Fr * 1000000) * (2 * 3.14159265359 * Fr * 1000000) * re1 * 0.000000000001)) / 0.000000001);
    c3 = (300 / Fr);
    c3 = c3 / 4;
    re2 = c3;
    c4 = ((Len * 0.0254) / re2) * 100;
    c5 = ((Math.pow(((Len / 12) / 984 * Fr * 360), 2)) / 312);
    document.calculatorWA.whip_antenna_capacitance.value = c1;
    document.calculatorWA.whip_antenna_inductance.value = c2;
    document.calculatorWA.whip_antenna_quarter_wavelength.value = c3;
    document.calculatorWA.whip_antenna_length.value = c4;
    document.calculatorWA.whip_antenna_radiation_resistance.value = c5;
}


function MW() {
    Er = Number(document.calculatorMW.number1.value);
    H = Number(document.calculatorMW.number2.value);
    t = Number(document.calculatorMW.number3.value);
    z = Number(document.calculatorMW.number4.value);
    c5 = z * (Math.pow((Er + 1.41), 0.5) / 87);
    c5 = Math.exp(c5);
    W = ((5.98 * H) / c5) - t;
    W = W / 0.8;
    c1 = W;
    document.calculatorMW.microstrip_width.value = c1;
}

function MMB() {
    W = Number(document.calculatorMMB.number1.value);
    H = Number(document.calculatorMMB.number2.value);
    c1 = W * Math.pow(2, 0.5);
    D = c1;
    X = D * (0.52 + 0.65 * Math.exp(-1.35 * (W / H)));
    c2 = X;
    A = (X - D / 2) * Math.pow(2, 0.5);
    c3 = A;
    document.calculatorMMB.output_D.value = c1;
    document.calculatorMMB.output_X.value = c2;
    document.calculatorMMB.output_A.value = c3;
}

function EMI() {
    Er = Number(document.calculatorEMI.number1.value);
    w = Number(document.calculatorEMI.number2.value);
    t = Number(document.calculatorEMI.number3.value);
    h = Number(document.calculatorEMI.number4.value);
    h1 = Number(document.calculatorEMI.number5.value);
    Erp = Er * (1 - Math.exp(-1.55 * h1 / h));
    Zo = 60 / Math.sqrt(Erp) * Math.log(5.98 * h / (0.8 * w + t));
    c1 = Zo;
    document.calculatorEMI.characteristic_impedance.value = c1;
}

function DMI() {
    Er = Number(document.calculatorDMI.number1.value);
    w = Number(document.calculatorDMI.number2.value);
    d = Number(document.calculatorDMI.number3.value);
    t = Number(document.calculatorDMI.number4.value);
    h = Number(document.calculatorDMI.number5.value);
    Zo = 87 / Math.sqrt(Er + 1.41) * Math.log(5.98 * h / (0.8 * w + t));
    Zd = 2 * Zo * (1 - 0.48 / Math.exp(0.96 * (d / h)));
    c1 = Zo;
    c2 = Zd;
    document.calculatorDMI.characteristic_impedance.value = c1;
    document.calculatorDMI.differential_impedance.value = c2;
}


function AGC() {
    antenna_efficiency = Number(document.calculatorAGC.number1.value);
    antenna_diameter = Number(document.calculatorAGC.number2.value);
    antenna_frequency = Number(document.calculatorAGC.number3.value);
    antenna_aperture = Math.PI * Math.pow(antenna_diameter, 2) / 4;
    antenna_wavelength = ((3 * Math.pow(10, 8)) / (antenna_frequency * Math.pow(10, 9)));
    var c = (antenna_efficiency * 4 * Math.PI * antenna_aperture) / (Math.pow(antenna_wavelength, 2));
    c1 = 10 * Math.log10(c);
    c1 = Math.round(c1);
    document.calculatorAGC.Gain.value = c1;
}


function ARC() {
    antenna_gain = Number(document.calculator.number1.value);
    system_temp = Number(document.calculator.number2.value);
    var c = antenna_gain - 10 * Math.log10(system_temp);
    document.calculator.GT.value = c;
}


function CCI() {
    D = Number(document.calculatorCCI.number1.value);
    d1 = Number(document.calculatorCCI.number2.value);
    Eps = Number(document.calculatorCCI.number3.value);
    a2 = Math.pow(Eps, 0.5);
    c1 = Math.log10(D / d1);
    c = ((138 * c1) / a2);
    mu = 4 * Math.PI * Math.pow(10, -7);
    E0 = 8.854 * Math.pow(10, -12);
    c2 = ((mu / (2 * Math.PI)) * (Math.log(D / d1)));
    c3 = ((2 * Math.PI * E0 * Eps) / (Math.log(D / d1)));
    document.calculatorCCI.coaxialimpedance.value = c;
    document.calculatorCCI.inductance.value = c2;
    document.calculatorCCI.capacitance.value = c3;
}


function NTNS() {
    a = Number(document.calculatorNTNS.number1.value);
    c = 10 * Math.log10(1 + (a / 290));
    document.calculatorNTNS.noisefigure.value = c;
}

function RRC() {

    Pt = Number(document.calculatorRRC.number1.value);
    G = Number(document.calculatorRRC.number2.value);
    Rc = Number(document.calculatorRRC.number3.value);
    Ae = Number(document.calculatorRRC.number4.value);
    Pmin = Number(document.calculatorRRC.number5.value);
    //p1 = (Pt/10);
    //Pt = ( Math.pow(10, p1 ) /1000 ) ;
    //p2 = (Pmin/10);
    //Pmin= ( Math.pow(10, p2 ) /1000 ) ;
    G = Math.pow(10, (G / 10));
    a1 = 4 * Math.PI;
    a2 = Math.pow(a1, 2);
    c1 = ((Pt * G * Rc * Ae) / (a2 * Pmin));
    c1 = Math.pow(c1, 0.25);
    document.calculatorRRC.ANF.value = c1;
}


function RFC() {
    inductance = Number(document.calculatorRFC.number1.value);
    capacitance = Number(document.calculatorRFC.number2.value);
    c = inductance * capacitance * Math.pow(10, -6);
    var c = 1 / (2 * Math.PI * Math.pow(c, 0.5));
    document.calculatorRFC.Fr.value = c;
}


function WTFC() {
    a = Number(document.calculatorWTFC.number1.value);
    // b=Number(document.calculator.number2.value);
    c = (3 * Math.pow(10, 8)) / (a * Math.pow(10, 9));
    document.calculatorWTFC.wavelength.value = c;
}

function RCNR() {
    GbyT = Number(document.calculatorRCNR.number1.value);
    EIRP = Number(document.calculatorRCNR.number2.value);
    PLoss = Number(document.calculatorRCNR.number3.value);
    margins = Number(document.calculatorRCNR.number4.value);
    //k=1.38*Math.pow(10,-23); log10(k)=228.6;
    c = EIRP - margins - PLoss + GbyT + 228.6;
    document.calculatorRCNR.CNratio.value = c;
}


function RFLB() {
    f = Number(document.calculatorRFLB.number1.value);
    d = Number(document.calculatorRFLB.number2.value);
    Gt = Number(document.calculatorRFLB.number3.value);
    Gr = Number(document.calculatorRFLB.number4.value);
    Pt = Number(document.calculatorRFLB.number5.value);
    lambda = ((3 * Math.pow(10, 8)) / (f * Math.pow(10, 9)));
    c1 = -20 * Math.log10(lambda) + 20 * Math.log10(d) + 21.98;
    c2 = Pt - c1 + Gt + Gr;
    document.calculatorRFLB.pathloss.value = c1;
    document.calculatorRFLB.receiverpower.value = c2;
}


function PDA() {
    Fr = Number(document.calculatorPDA.number1.value);
    D = Number(document.calculatorPDA.number2.value);
    lambda = ((3 * (Math.pow(10, 8))) / (Fr * (Math.pow(10, 9))));
    m1 = D / lambda;
    c1 = 6 * Math.pow(m1, 2);
    c1 = (10 * (Math.log10(c1)));
    c2 = (60 * (lambda / D));
    c3 = (0.6 * ((Math.PI * (Math.pow(D, 2))) / 4));
    document.calculatorPDA.Gain.value = c1;
    document.calculatorPDA.half_power_beamwidth.value = c2;
    document.calculatorPDA.effective_aperture.value = c3;
}


function PRF() {
    PRF = Number(document.calculatorPRF.number1.value);
    c1 = ((3 * (Math.pow(10, 8))) / (2 * PRF));
    document.calculatorPRF.range.value = c1;
}


function RRR() {
    T = Number(document.calculatorRRR.number1.value);
    c1 = (((3 * (Math.pow(10, 8))) * (T * (Math.pow(10, -6)))) / 2);
    document.calculatorRRR.resolution.value = c1;
}

function ANF() {
    d = Number(document.calculatorANF.number1.value);
    Fr = Number(document.calculatorANF.number2.value);
    lambda = ((3 * (Math.pow(10, 8))) / (Fr * (Math.pow(10, 9))));
    c1 = (2 * (Math.pow(d, 2))) / lambda;
    fac = Math.pow(d, 3) / lambda;
    c2 = 0.62 * Math.pow(fac, 0.5);
    document.calculatorANF.reactive_nearfield_distance.value = c2;
    document.calculatorANF.antenna_nearfield_distance.value = c1;
}

function HA() {


    Fr = Number(document.calculatorHA.number1.value);
    b = Number(document.calculatorHA.number2.value);
    a = Number(document.calculatorHA.number3.value);
    Area = (b * a);
    lambda = ((3 * (Math.pow(10, 8))) / (Fr * (Math.pow(10, 9))));
    c1 = ((10 * Area) / (Math.pow(lambda, 2)));
    c1 = (10 * Math.log10(c1));
    c2 = ((51 * lambda) / b);
    c3 = ((70 * lambda) / a);
    document.calculatorHA.Gain.value = c1;
    document.calculatorHA.vertical_beamwidth.value = c2;
    document.calculatorHA.horizontal_beamwidth.value = c3;

}


function BPF() {
    Fc = Number(document.calculatorBPF.number1.value);
    BW3 = Number(document.calculatorBPF.number2.value);
    BW60 = Number(document.calculatorBPF.number3.value);
    c1 = (BW60 / BW3);
    c2 = ((Fc * (Math.pow(10, 3))) / BW3);
    document.calculatorBPF.BPF_shape_factor.value = c1;
    document.calculatorBPF.BPF_quality_factor.value = c2;
}


function RHC() {
    h = Number(document.calculatorRHC.number1.value);
    c1 = 1.42 * Math.pow(h, 0.5);
    document.calculatorRHC.radiohorizon.value = c1;
}

function EAG() {


    c = Number(document.calculatorEAG.number1.value);
    calls_lost = Number(document.calculatorEAG.number2.value);
    h = Number(document.calculatorEAG.number3.value);
    //Busy hour gives value of T equal to 60
    T = 60;
    c1 = ((c * h) / T);
    c2 = (calls_lost / c);
    document.calculatorEAG.Erlang.value = c1;
    document.calculatorEAG.GoS.value = c2;

}

function EVM() {
    EVMdB = Number(document.calculatorEVM.number1.value);
    m1 = (EVMdB / 20);
    c1 = (100 * (Math.pow(10, m1)));
    document.calculatorEVM.EVMrms.value = c1;
}

function GSM() {
    N = Number(document.calculatorGSM.number1.value);
    if (((N >= 259) && (N <= 293))) {
        c1 = 450.6 + 0.2 * (N - 259);
        c2 = c1 + 10;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 306) && (N <= 340))) {
        c1 = 479 + 0.2 * (N - 306);
        c2 = c1 + 10;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 438) && (N <= 511))) {
        c1 = 747.2 + 0.2 * (N - 438);
        c2 = c1 + 30;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 128) && (N <= 251))) {
        c1 = 824.2 + 0.2 * (N - 128);
        c2 = c1 + 45;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 1) && (N <= 124))) {
        c1 = 890 + 0.2 * N;
        c2 = c1 + 45;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 975) && (N <= 1023))) {
        c1 = 890 + 0.2 * (N - 1024);
        c2 = c1 + 45;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 940) && (N <= 974))) {
        c1 = 890 + 0.2 * (N - 1024);
        c2 = c1 + 45;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else if (((N >= 512) && (N <= 810))) {
        //i=1;
        c1 = 1710.2 + 0.2 * (N - 512);
        c2 = c1 + 95;
        str2 = " for DCS1800;";
        str3 = " for PCS1900";
        c3 = 1850.2 + 0.2 * (N - 512);
        c4 = c1 + 80;
        str1 = c1 + str2 + c3 + str3;
        str4 = c2 + str2 + c4 + str3
        document.calculatorGSM.Ful.value = str1;
        document.calculatorGSM.Fdl.value = str4;
    }
    else if (((N >= 811) && (N <= 885))) {
        c1 = 1710.2 + 0.2 * (N - 512);
        c2 = c1 + 95;
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
    else {
        c1 = "INVALID INPUT";
        c2 = "INVALID INPUT";
        document.calculatorGSM.Ful.value = c1;
        document.calculatorGSM.Fdl.value = c2;
    }
}

function UMTS() {
    Ndl = Number(document.calculatorUMTS.number1.value);
    if (((Ndl >= 10562) && (Ndl <= 10838))) {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = Ndl - 10562;
        Nul = 9612 + diff;
    }
    else if (((Ndl >= 9662) && (Ndl <= 9938))) {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = Ndl - 9662;
        Nul = 9262 + diff;
    }
    else if (((Ndl >= 1162) && (Ndl <= 1513))) {
        FDL_offset = 1575;
        FUL_offset = 1525;
        diff = Ndl - 1162;
        Nul = 937 + diff;
    }
    else if (((Ndl >= 1537) && (Ndl <= 1738))) {
        FDL_offset = 1805;
        FUL_offset = 1450;
        diff = Ndl - 1537;
        Nul = 1312 + diff;
    }
    else if (((Ndl >= 4357) && (Ndl <= 4458))) {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = Ndl - 4357;
        Nul = 4132 + diff;
    }
    else if (((Ndl >= 4387) && (Ndl <= 4413))) {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = Ndl - 4387;
        Nul = 4162 + diff;
    }
    else if (((Ndl >= 2237) && (Ndl <= 2563))) {
        FDL_offset = 2175;
        FUL_offset = 2100;
        diff = Ndl - 2237;
        Nul = 2012 + diff;
    }
    else if (((Ndl >= 2937) && (Ndl <= 3088))) {
        FDL_offset = 340;
        FUL_offset = 340;
        diff = Ndl - 2937;
        Nul = 2712 + diff;
    }
    else if (((Ndl >= 9237) && (Ndl <= 9387))) {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = Ndl - 9237;
        Nul = 8762 + diff;
    }
    else if (((Ndl >= 3112) && (Ndl <= 3388))) {
        FDL_offset = 1490;
        FUL_offset = 1135;
        diff = Ndl - 3112;
        Nul = 2887 + diff;
    }
    else if (((Ndl >= 3712) && (Ndl <= 3787))) {
        FDL_offset = 736;
        FUL_offset = 733;
        diff = Ndl - 3712;
        Nul = 3487 + diff;
    }
    else if (((Ndl >= 3842) && (Ndl <= 3903))) {
        FDL_offset = -37;
        FUL_offset = -22;
        diff = Ndl - 3842;
        Nul = 3617 + diff;
    }
    else if (((Ndl >= 4017) && (Ndl <= 4043))) {
        FDL_offset = -55;
        FUL_offset = 21;
        diff = Ndl - 4017;
        Nul = 3792 + diff;
    }
    else if (((Ndl >= 4117) && (Ndl <= 4143))) {
        FDL_offset = -63;
        FUL_offset = 12;
        diff = Ndl - 4117;
        Nul = 3892 + diff;
    }
    else if (((Ndl >= 712) && (Ndl <= 763))) {
        FDL_offset = 735;
        FUL_offset = 770;
        diff = Ndl - 712;
        Nul = 312 + diff;
    }
    else if (((Ndl >= 4512) && (Ndl <= 4638))) {
        FDL_offset = -109;
        FUL_offset = -23;
        diff = Ndl - 4512;
        Nul = 4287 + diff;
    }
    else if (((Ndl >= 862) && (Ndl <= 912))) {
        FDL_offset = 1326;
        FUL_offset = 1358;
        diff = Ndl - 862;
        Nul = 462 + diff;
    }
    else if (((Ndl >= 4662) && (Ndl <= 5038))) {
        FDL_offset = 2580;
        FUL_offset = 2525;
        diff = Ndl - 4662;
        Nul = 4437 + diff;
    }
    else if (((Ndl >= 5112) && (Ndl <= 5413))) {
        FDL_offset = 910;
        FUL_offset = 875;
        diff = Ndl - 5112;
        Nul = 4887 + diff;
    }
    else if (((Ndl >= 5762) && (Ndl <= 5913))) {
        FDL_offset = -291;
        FUL_offset = -291;
        diff = Ndl - 5762;
        Nul = 5537 + diff;
    }
    else {
        FDL_offset = 0;
        FUL_offset = 0;
        diff = 0;
        Nul = 0;
    }
    c1 = (FDL_offset + (0.2 * Ndl));
    c2 = (FUL_offset + (0.2 * Nul));
    c3 = Nul;
    document.calculatorUMTS.Fdl.value = c1;
    document.calculatorUMTS.Ful.value = c2;
    document.calculatorUMTS.UL_UARFCN.value = c3;
}


function SSRC() {
    h = Number(document.calculatorSSRC.number1.value);
    B = Number(document.calculatorSSRC.number2.value);
    Theta = Number(document.calculatorSSRC.number3.value);
    r = h + B;
    r = r * 1000;
    B = B * 1000;
    Epsilon = 90 - Theta;
    E = Math.cos(Epsilon);
    p1 = Math.pow((B * E), 2);
    p2 = Math.pow(r, 2);
    p3 = Math.pow(B, 2);
    p4 = B * Math.cos(Epsilon);
    c1 = p1 + p2 - p3;
    c1 = Math.pow(c1, 0.5);
    c1 = c1 - p4;
    c1 = c1 / 1000;
    document.calculatorSSRC.slant_range.value = c1;
}


function GS() {
    d = Number(document.calculatorGS.number1.value);
    go = 398600.5;
    re = 6378;
    m1 = (go / d);
    c1 = Math.pow(m1, 0.5);
    c2 = ((2 * (Math.PI) * (Math.pow(d, 1.5))) / (Math.pow(go, 0.5)));
    c3 = ((Math.pow(go, 0.5)) / (Math.pow(d, 1.5)));
    c4 = (go / (Math.pow(go, 2)));
    document.calculatorGS.velocity.value = c1;
    document.calculatorGS.orbit_period.value = c2;
    document.calculatorGS.angular_velocity.value = c3;
    document.calculatorGS.acceleration.value = c4;
}


function ESSLB() {
    F = Number(document.calculatorESSLB.number1.value);
    Dia = Number(document.calculatorESSLB.number2.value);
    Pt = Number(document.calculatorESSLB.number3.value);
    d = Number(document.calculatorESSLB.number4.value);
    eff = Number(document.calculatorESSLB.number5.value);
    //F=F*Math.pow(10,9);
    c1 = 92.4 + 20 * Math.log10(F) + 20 * Math.log10(d);
    c2 = 20.4 + 20 * Math.log10(F) + 20 * Math.log10(Dia) + 10 * Math.log10(eff);
    Pt = 10 * Math.log10(Pt);
    c3 = Pt + c2;
    c4 = c3 - c1;
    document.calculatorESSLB.path_loss.value = c1;
    document.calculatorESSLB.antenna_gain.value = c2;
    document.calculatorESSLB.EIRP.value = c3;
    document.calculatorESSLB.satellite_receive_power.value = c4;
}


function TNPC() {
    T = Number(document.calculatorTNPC.number1.value);
    R = Number(document.calculatorTNPC.number2.value);
    B = Number(document.calculatorTNPC.number3.value);
    K = 1.3807e-23;
    //K=1.3806503e-23
    c1 = 10 * Math.log10(1000 * K * T * B);
    c2 = Math.pow((4 * K * T * B * R), 0.5) * 1e6;
    document.calculatorTNPC.noisepower.value = c1;
    document.calculatorTNPC.noisevoltage.value = c2;
}

function AEC() {
    antenna_diameter = Number(document.calculatorAEC.number1.value);
    antenna_gain = Number(document.calculatorAEC.number2.value);
    antenna_frequency = Number(document.calculatorAEC.number3.value);
    //antenna_wavelength=(Math.pow(3,10)/(antenna_frequency*Math.pow(10,6)));
    c1 = 100 * Math.pow(10, ((-20 * Math.log10(antenna_frequency / 1000) - 20 * Math.log10(antenna_diameter) - 20.4 + antenna_gain) / 10));
    document.calculatorAEC.antenna_efficiency.value = c1;
}


function EBNO() {
    EbbyN0 = Number(document.calculatorEBNO.number1.value);
    bit_rate = Number(document.calculatorEBNO.number2.value);
    BW = Number(document.calculatorEBNO.number3.value);
    //k=1.38*Math.pow(10,-23); log10(k)=228.6;
    CbyN = EbbyN0 + 10 * Math.log10(bit_rate / BW);
    document.calculatorEBNO.CNR.value = CbyN;
}


function FVFPC() {
    floating_number = Number(document.calculatorFVFPC.number3.value);
    Q_format = Number(document.calculatorFVFPC.number4.value);
    //k=1.38*Math.pow(10,-23); log10(k)=228.6;
    fixed_number = floating_number * Math.pow(2, Q_format);
    fixed_number = Math.round(fixed_number);
    document.calculatorFVFPC.fixed.value = fixed_number;
}


function TTC() {
    a = Number(document.calculatorTTC.number1.value);
    c = a * (Math.pow(10, 12));
    document.calculatorTTC.freq_Hz.value = c;
}

function DAC() {
    Lphy = Number(document.calculatorDAC.number1.value);
    Fr = Number(document.calculatorDAC.number2.value);
    antenna_wavelength = (3 * Math.pow(10, 8) / (Fr * Math.pow(10, 6)));
    Leff1 = (2 * Lphy) / Math.PI;
    Leff2 = Lphy / 2;
    Leff3 = Lphy;
    val1 = Math.PI;
    c1 = 80 * Math.pow(val1, 2) * Math.pow((Leff1 / antenna_wavelength), 2);
    c2 = 80 * Math.pow(val1, 2) * Math.pow((Leff2 / antenna_wavelength), 2);
    c3 = 80 * Math.pow(val1, 2) * Math.pow((Leff3 / antenna_wavelength), 2);
    document.calculatorDAC.Radiation_resistance1.value = c1;
    document.calculatorDAC.Radiation_resistance2.value = c2;
    document.calculatorDAC.Radiation_resistance3.value = c3;
}

function CCCIC() {
    b = Number(document.calculatorCCCIC.number1.value);
    a = Number(document.calculatorCCCIC.number2.value);
    epsilon = Number(document.calculator.number3.value);
    Mu = Number(document.calculator.number4.value);
    epsilon1 = 8.85419 * Math.pow(10, -12) * epsilon;
    c1 = (2 * Math.PI * epsilon1) / (Math.log(b / a));
    Mu1 = 4 * Math.PI * Math.pow(10, -7) * Mu;
    c2 = (Mu1 / (2 * Math.PI)) * (Math.log(b / a));
    document.calculatorCCCIC.capacitance.value = c1;
    document.calculatorCCCIC.inductance.value = c2;
}

function CCC() {

    BW = Number(document.calculatorCCC.number1.value);
    SNR = Number(document.calculatorCCC.number2.value);
    c = 1 + SNR;
    c1 = BW * Math.log2(c);
    document.calculatorCCC.CC.value = c1;

}

function DFC() {
    v = Number(document.calculatorDFC.number1.value);
    Fr = Number(document.calculatorDFC.number2.value);
    lambda = (3 * Math.pow(10, 8) / (Fr * Math.pow(10, 6)));
    c1 = (2 * v) / lambda;
    document.calculatorDFC.Doppler_Frequency.value = c1;
}


function RBSC() {
    Fr = Number(document.calculatorRBSC.number1.value);
    Ts = Number(document.calculatorRBSC.number2.value);
    lambda = ((3 * Math.pow(10, 8)) / (Fr * Math.pow(10, 9)));
    c1 = lambda / (2 * Ts * Math.pow(10, -6));
    document.calculatorRBSC.Blind_Speed.value = c1;
}

function TVCN() {

    N = Number(document.calculatorTVCN.number1.value);
    if (((N >= 2) && (N <= 4))) {
        var1 = 42;
        c = (N * 6. + var1) + ' - ' + (N * 6. + var1 + 6);
        document.calculatorTVCN.TV_channel_Band.value = c;
    }
    else if (((N >= 5) && (N <= 6))) {
        var1 = 46;
        c = (N * 6. + var1) + ' - ' + (N * 6. + var1 + 6);
        document.calculatorTVCN.TV_channel_Band.value = c;
    }
    else if (((N >= 7) && (N <= 13))) {
        var1 = 132;
        c = (N * 6. + var1) + ' - ' + (N * 6. + var1 + 6);
        document.calculatorTVCN.TV_channel_Band.value = c;
    }
    else if (((N >= 14) && (N <= 69))) {
        var1 = 386;
        c = (N * 6. + var1) + ' - ' + (N * 6. + var1 + 6);
        document.calculatorTVCN.TV_channel_Band.value = c;
    }
    else if (((N >= 70) && (N <= 83))) {
        var1 = 386;
        c = (N * 6. + var1) + ' - ' + (N * 6. + var1 + 6);
        document.calculatorTVCN.TV_channel_Band.value = c;
    }
    else {
        c = "INVALID INPUT";
        document.calculatorTVCN.TV_channel_Band.value = c;
    }


}

function FMCN() {


    a = Number(document.calculatorFMCN.number1.value);
    c = ((((a - 200) * 0.2) + 87.9) * 10) / 10;
    document.calculatorFMCN.FM_freq.value = c;

}