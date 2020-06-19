from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np


# Size of variable arrays:
size_algebraic = 72 
size_states = 41
size_constants = 73


def init_consts():
    constants = [0.0] * size_constants 
    states = [0.0] * size_states
    
    # Table 8. Initial conditions
    states[0] = -82.4202  # V; Membrane potential:mV
    states[1] = 0.115001  # Cai; Myoplasmic Ca+ concentration:uM
    states[2] = 0.115001  # Cass; Subspace SR Ca+ concentration:uM
    states[3] = 1299.5  # CaJSR; JSR Ca+ concentration:uM
    states[4] = 1299.5  # CaNSR; NSR Ca+ concentration:uM
    states[5] = 0  # PRyR; RyR modulation factor
    states[6] = 11.2684  # LTRPNCa; Concentration Ca+ bound low-affinity troponin-binding sites:uM
    states[7] = 125.29  # HTEPNCa; Concentration Ca+ bound high-affinity troponin-binding sites:uM
    states[8] = 0.149102e-4  # PO1; %Fraction of RyR channels in sate POl
    states[9] = 0.951726e-10  # PO2; %Fraction of RyR channels in sate POl
    states[10] = 0.16774e-3  # PC2; Fraction of RyR channels in sate PC2
    states[11] = 0.930308e-18  # O; L-type Ca+ channel conducting state
    states[12] = 0.124216e-3  # C2; L-type Ca+ channel closed state
    states[13] = 0.578679e-8  # C3; L-type Ca+ channel closed state
    states[14] = 0.119816e-12  # C4; L-type Ca+ channel closed state
    states[15] = 0.497923e-18  # I1; L-type Ca+ channel inacticated state
    states[16] = 0.345847e-13  # I2; L-type Ca+ channel inacticated state
    states[17] = 0.185106e-13  # I3; L-type Ca+ channel inacticated state
    states[18] = 14237.1  # Nai; Myoplasmic Na+ concentration:uM
    states[19] = 0.713483e-6  # ONa; Open state of fast Na+ channel
    states[20] = 0.279132e-3  # CNa1; Closed state of fast Na+ channel
    states[21] = 0.020752  # CNa2; Closed state of fast Na+ channel
    states[22] = 0.673345e-6  # I1Na; Slow inactivated state 1 of fast Na+ channel
    states[23] = 0.155787e-8  # I2Na; Slow inactivated state 2 of fast Na+ channel
    states[24] = 0.153176e-3  # IFNa; Fast inactivated state of fast Na+ channel
    states[25] = 0.0113879  # ICNa2; Cloesd-inactivated state of fast Na+ channel
    states[26] = 0.34278  # ICNa3; Cloesd-inactivated state of fast Na+ channel
    states[27] = 143720  # Ki; Myoplasmic K+ concentration:uM
    states[28] = 0.265563e-2  # ato_f; Gating variable for transient outward K+ current
    states[29] = 0.999977  # ito_f; Gating variable for transient outward K+ current
    states[30] = 0.417069e-3  # ato_s; Gating variable for transient outward K+ current
    states[31] = 0.998543  # ito_s; Gating variable for transient outward K+ current
    states[32] = 0.262753e-3  # nKs; Gating variable for slow delayed-rectifier K+ current
    states[33] = 0.417069e-3  # aur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    states[34] = 0.998543  # iur; Gating variable for ultrarapidly activating delayed-rectifier K+ current
    states[35] = 0.417069e-3  # aKss; Gating variable for noninactivating steady-state K+ current
    states[36] = 1  # iKss; Gating variable for noninactivating steady-state K+ current
    states[37] = 0.175298e-3  # OK; mERG channel open state
    states[38] = 0.992513e-3  # CK1; mERG channel closed state
    states[39] = 0.641229e-3  # CK2; mERG channel closed state
    states[40] = 0.319129e-4  # IK; mERG channel inactivated state

    # Table 2. Cell geometry parameters
    constants[1] = 25.84e-6  # Vmyo; Myoplasmic volume:ul
    constants[2] = 0.12e-6  # VJS; Juncional SR volume:ul
    constants[3] = 2.098e-6  # VNSR; Network SR volume:ul
    constants[4] = 1.485e-9  # Vss; Subspace volume:ul
    constants[5] = 1.534e-4  # Acap; Capacitive membrane area:cm^2

    # Table 3. Extracellular ion concentrations
    constants[6] = 5400  # % Ko; Exracellular K+ concentration:uM
    constants[7] = 140000  # % Nao; Exracellular Na+ concentration:uM
    constants[8] = 1800  # % Cao; Exracellular Ca+ concentration:uM

    # Table 4. SR paramters
    constants[25] = 4.5  # nu1; Maximum RyR channel Ca+ permeability:ms^-1
    constants[27] = 1.74e-5  # nu2; Ca+ leak rate constant from the NSR:ms^-1
    constants[26] = 20  # tautr; Time constant for transfer from NSR to JSR:ms
    constants[28] = 8  # tauxfer; Time constant for transfer from subspace to myoplasm:ms
    constants[29] = 0.45  # nu3; SR Ca+-ATPase maximum pump rate:uM/ms
    constants[30] = 0.5  # Km_up; Half-saturation constant for SR Ca+-ATPase pump:uM
    constants[34] = 0.006075  # ka1; RyR Pc1-Po1 rate constant:uM^-4/ms
    constants[35] = 0.07125  # ka2; RyR Po1-Pc1 rate constant:ms^-1
    constants[36] = 0.00405  # kb1; RyR Po1-Po2 rate constant:uM^-3/ms
    constants[37] = 0.965  # kb2; RyR Po2-Po1 rate constant:ms^-1
    constants[38] = 0.009  # kc1; RyR Po1-Pc2 rate constant:ms^-1
    constants[39] = 0.0008  # kc2; RyR Pc2-Po1 rate constant:ms^-1
    constants[40] = 3  # m; RyR Ca+ cooperativity parameter Po1-Po2
    constants[41] = 4  # n; RyR Ca+ cooperativity parameter Pc1-Po1

    # Table 5. L-type Ca+ channel parameters
    constants[43] = 0.1729  # GCaL; Specific maximum conductivity for L-type Ca+ channel:mS/uF
    constants[42] = 63  # ECa; Reversal potential for L-type Ca+ channel:mV
    constants[45] = 0.23324  # Kpc_max; Maximum time constant for Ca+-induced inactivation:ms^-1
    constants[46] = 20  # Kpc_half; Half-saturation constant for Ca+-induced inactivation:uM
    constants[44] = 0.0005  # Kpcb; Voltage-insensitive rate constant for inactivation:ms^-1
    constants[33] = 7  # ICaL_max; Normalization constant for L-type Ca+ current:pA/pF 

    # Table 6. Buffering parameters
    constants[31] = 70  # LTRPN_tot; Total myoplasmic troponin low-affinity site concentration:uM
    constants[32] = 140  # HTRPN_tot; Total myoplasmic troponin high-affinity site concentration:uM
    constants[21] = 0.00237  # khtrpn1; Ca+ on rate constant for troponin high-affinity sites:uM^-1/ms
    constants[22] = 3.2e-5  # khtrpn2; Ca+ off rate constant for troponin high-affinity sites:ms^-1
    constants[23] = 0.0327  # kltrpn1; Ca+ on rate constant for troponin low-affinity sites:uM^-1/ms
    constants[24] = 0.0196  # kltrpn2; Ca+ off rate constant for troponin low-affinity sites:ms^-1
    constants[17] = 50  # CMDN_tot; Total myoplasmic calmodulin concentration:uM
    constants[18] = 15000  # CSQN_tot; Total junctional SR calsequestrin concentration:uM
    constants[19] = 0.238  # Km_CMDN; Ca+ half-saturaion constant for calmodulin:uM
    constants[20] = 800  # Km_CSQN; Ca+ half-saturaion constant for calsequestrin:uM

    # Table 7. Membrane current parameters
    constants[11] = 96.5  # F; Faraday constant:C/mmol
    constants[10] = 298  # T; Absolute temperature:K
    constants[9] = 8.314  # R; Ideal gas constant:J*mol^-1*K^-1
    constants[49] = 292.8  # kNaCa; Scaling factor of Na+/Ca+ exchange:pA/pF
    constants[50] = 87500  # Km_Na; Na+ half-saturation constant for Na+/Ca+ exchange:uM
    constants[51] = 1380  #  Km_Ca; Ca+ half-saturation constant for Na+/Ca+ exchange:uM
    constants[52] = 0.1  # ksat; Na+/Ca+ exchange saaturation factor at very negative potentials
    constants[53] = 0.35  # yita; Contals voltage dependence of Na+/Ca+ exchange
    constants[65] = 0.88  # INaK_max; Maximum Na+/K+ exchange current:pA/pF 
    constants[66] = 21000  # Km_Nai; Na+ half-saturation constant for Na+/K+ exchange currant:uM
    constants[67] = 1500  # Km_Ko; K+ half-saturation constant for Na+/K+ exchange currant:uM
    constants[47] = 1  # Ip_Ca_max; Maximun Ca+ pump current:pA/pF
    constants[48] = 0.5  # Km_p_Ca; Ca+ half-saturation constant for Ca+pump current:uM
    constants[54] = 0.000367  # GCab; Maximun background Ca+ current conductance:mS/uF
    constants[55] = 13  # GNa; Maximun fast Na+ current conductance:mS/uF
    constants[56] = 0.0026  # GNab; Maximun background Na+ current conductance:mS/uF
    constants[59] = 0.00575  # GKs; Maximum slow delayed-rectifier K+ current conductance:mS/uF

    # Table 7. Parameters for apex
    constants[57] = 0.4067  # GKtof; Maximum transient outward K+ current conductance(apex):mS/uF
    constants[58] = 0  # GKtos; Maximum transient outward K+ current conductance(apex):mS/uF
    constants[60] = 0.16  # GKur; Maximum ultrarapidly delayed-rectifier K+ current conductance(apex):mS/uF
    constants[61] = 0.05  # GKss; Maximum noninactivating steady-state K+ current conductance(apex):mS/uF

    # Table 7. Membrane current parameters
    constants[62] = 0.078  # GKr; Maximum rapid delayed-rectifier K+ current conductance:mS/uF
    constants[64] = 0.023761  # kf; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    constants[63] = 0.036778  # kb; Rate constant for rapid delayed-rectifier K+ current:ms^-1
    constants[68] = 10  # GCl_Ca; Maximum Ca+-activated Cl- current conductance:mS/uF
    constants[70] = 10  # Km_Cl; Half-saturaon constant for Ca+-activated Cl- current:uM
    constants[69] = -40  # ECl; Reversal potential for Ca+-activated Cl- current:mV
    
    # Parameters not in the SimuLink model
    constants[0] = 1
    constants[12] = 20
    constants[13] = 100000
    constants[14] = 71.43
    constants[15] = 0.5
    constants[16] = -80
    constants[71] = (1.00000/7.00000)*(np.exp(constants[7]/67300.0)-1.00000)
    constants[72] = 0.00000

    return (states, constants)

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return np.select(cases[0::2],cases[1::2])

def voltage_clamp(t, holding_p, holding_t, P1, P1_t, P2):
    v = np.piecewise(t, 
                     [t < holding_p, holding_t <= t <= P1_t, t > P1_t],
                     [holding_p, P1, P2])

    return v

def compute_rates(voi, states, constants, holding_p, holding_t, P1, P1_t, P2):
    rates = [0.0] * size_states 
    algebraic = [0.0] * size_algebraic

    # externally applied voltage (voltage clamp)
    algebraic[71] = voltage_clamp(voi, holding_p, holding_t, P1, P1_t, P2)

    # A94; iKss if CONSTANTS(:,72)=0  or A110; sigma if CONSTANTS(:,72)=A110
    rates[36] = constants[72]
    # A16; LTRPNCa
    rates[6] = constants[23]*states[1]*(constants[31]-states[6])-constants[24]*states[6]
    # A17; HTEPNCa
    rates[7] = constants[21]*states[1]*(constants[32]-states[7])-constants[22]*states[7]
    # A20; PO2
    rates[9] = constants[36]*(np.power(states[2], constants[40]))*states[8]-constants[37]*states[9]
    # A21; PC2
    rates[10] = constants[38]*states[8]-constants[39]*states[10]
    # A19; PC1
    algebraic[1] = 1.00000-(states[10]+states[8]+states[9])
    # A18; PO1
    rates[8] = (constants[34]*(np.power(states[2], constants[41]))*algebraic[1]+constants[37]*states[9]+constants[39]*states[10])-(constants[35]*states[8]+constants[36]*(np.power(states[2], constants[40]))*states[8]+constants[38]*states[8])
    # A71; alpha_a
    algebraic[4] = 0.180640*np.exp(0.0357700*(states[0]+30.0000))
    # A72; beta_a
    algebraic[14] = 0.395600*np.exp(-0.0623700*(states[0]+30.0000))
    # A69; ato_f
    rates[28] = algebraic[4]*(1.00000-states[28])-algebraic[14]*states[28]
    # A73; alpha_i
    algebraic[5] = (0.000152000*np.exp(-(states[0]+13.5000)/7.00000))/(0.00670830*np.exp(-(states[0]+33.5000)/7.00000)+1.00000)
    # A74; beta_i
    algebraic[15] = (0.000950000*np.exp((states[0]+33.5000)/7.00000))/(0.0513350*np.exp((states[0]+33.5000)/7.00000)+1.00000)
    # A70; ito_f
    rates[29] = algebraic[5]*(1.00000-states[29])-algebraic[15]*states[29]
    # A78; a_ss
    algebraic[6] = 1.00000/(1.00000+np.exp(-(states[0]+22.5000)/7.70000))
    # A80; tau_tas
    algebraic[16] = 0.493000*np.exp(-0.0629000*states[0])+2.05800
    # A76; ato_s
    rates[30] = (algebraic[6]-states[30])/algebraic[16]
    # A79; i_ss
    algebraic[7] = 1.00000/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    # A81; tau_tis
    algebraic[17] = 270.000+1050.00/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    # A77; ito_s
    rates[31] = (algebraic[7]-states[31])/algebraic[17]
    # A85; alpha_n
    algebraic[8] = (4.81333e-06*(states[0]+26.5000))/(1.00000-np.exp(-0.128000*(states[0]+26.5000)))
    # A86; beta_n
    algebraic[18] = 9.53333e-05*np.exp(-0.0380000*(states[0]+26.5000))
    # A84; nKs
    rates[32] = algebraic[8]*(1.00000-states[32])-algebraic[18]*states[32]
    # A90; tau_aur
    algebraic[19] = 0.493000*np.exp(-0.0629000*states[0])+2.05800
    # A88; aur
    rates[33] = (algebraic[6]-states[33])/algebraic[19]
    # A91; tau_iur
    algebraic[20] = 1200.00-170.000/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    # A89; iur
    rates[34] = (algebraic[7]-states[34])/algebraic[20]
    # A95; tau_Kss
    algebraic[21] = 39.3000*np.exp(-0.0862000*states[0])+13.1700
    # A93; aKss
    rates[35] = (algebraic[6]-states[35])/algebraic[21]
    # A104; alpha_a1
    algebraic[10] = 0.0137330*np.exp(0.0381980*states[0])
    # A105; beta_a1
    algebraic[23] = 6.89000e-05*np.exp(-0.0417800*states[0])
    # A99; CK2
    rates[39] = (constants[64]*states[38]+algebraic[23]*states[37])-(constants[63]*states[39]+algebraic[10]*states[39])
    # A24; C1
    algebraic[2] = 1.00000-(states[11]+states[12]+states[13]+states[14]+states[15]+states[16]+states[17])
    # A31; alpha
    algebraic[12] = (0.400000*np.exp((states[0]+12.0000)/10.0000)*((1.00000+0.700000*np.exp(-(np.power(states[0]+40.0000, 2.00000))/10.0000))-0.750000*np.exp(-(np.power(states[0]+20.0000, 2.00000))/400.000)))/(1.00000+0.120000*np.exp((states[0]+12.0000)/10.0000))
    # A32; beta
    algebraic[25] = 0.0500000*np.exp(-(states[0]+12.0000)/13.0000)
    # A25; C2
    rates[12] = (4.00000*algebraic[12]*algebraic[2]+2.00000*algebraic[25]*states[13])-(algebraic[25]*states[12]+3.00000*algebraic[12]*states[12])
    # A26; C3
    rates[13] = (3.00000*algebraic[12]*states[12]+3.00000*algebraic[25]*states[14])-(2.00000*algebraic[25]*states[13]+2.00000*algebraic[12]*states[13])
    # A97; CK0
    algebraic[9] = 1.00000-(states[38]+states[39]+states[37]+states[40])
    # A102; alpha_a0
    algebraic[22] = 0.0223480*np.exp(0.0117600*states[0])
    # A103; beta_a0
    algebraic[27] = 0.0470020*np.exp(-0.0631000*states[0])
    # A98; CK1
    rates[38] = (algebraic[22]*algebraic[9]+constants[63]*states[39])-(algebraic[27]*states[38]+constants[64]*states[38])
    # A106; alpha_i
    algebraic[28] = 0.0908210*np.exp(0.0233910*(states[0]+5.00000))
    # A107; beta_i
    algebraic[32] = 0.00649700*np.exp(-0.0326800*(states[0]+5.00000))
    # A100; OK
    rates[37] = (algebraic[10]*states[39]+algebraic[32]*states[40])-(algebraic[23]*states[37]+algebraic[28]*states[37])
    # A101; IK
    rates[40] = algebraic[28]*states[37]-algebraic[32]*states[40]
    # A33; gamma
    algebraic[30] = (constants[45]*states[2])/(constants[46]+states[2])
    # A34; K_pcf
    algebraic[34] = 13.0000*(1.00000-np.exp(-(np.power(states[0]+14.5000, 2.00000))/100.000))
    # A23; O
    rates[11] = (algebraic[12]*states[14]+constants[44]*states[15]+0.00100000*(algebraic[12]*states[16]-algebraic[34]*states[11]))-(4.00000*algebraic[25]*states[11]+algebraic[30]*states[11])
    # A27; C4
    rates[14] = (2.00000*algebraic[12]*states[13]+4.00000*algebraic[25]*states[11]+0.0100000*(4.00000*constants[44]*algebraic[25]*states[15]-algebraic[12]*algebraic[30]*states[14])+0.00200000*(4.00000*algebraic[25]*states[16]-algebraic[34]*states[14])+4.00000*algebraic[25]*constants[44]*states[17])-(3.00000*algebraic[25]*states[14]+algebraic[12]*states[14]+1.00000*algebraic[30]*algebraic[34]*states[14])
    # A28; I1
    rates[15] = (algebraic[30]*states[11]+0.00100000*(algebraic[12]*states[17]-algebraic[34]*states[15])+0.0100000*(algebraic[12]*algebraic[30]*states[14]-4.00000*algebraic[25]*algebraic[34]*states[15]))-constants[44]*states[15]
    # A29; I2
    rates[16] = (0.00100000*(algebraic[34]*states[11]-algebraic[12]*states[16])+constants[44]*states[17]+0.00200000*(algebraic[34]*states[14]-4.00000*algebraic[25]*states[16]))-algebraic[30]*states[16]
    # A30; I3
    rates[17] = (0.00100000*(algebraic[34]*states[15]-algebraic[12]*states[17])+algebraic[30]*states[16]+1.00000*algebraic[30]*algebraic[34]*states[14])-(4.00000*algebraic[25]*constants[44]*states[17]+constants[44]*states[17])
    # A8; B_JSR
    algebraic[29] = np.power(1.00000+(constants[18]*constants[20])/(np.power(constants[20]+states[3], 2.00000)), -1.00000)
    # A9; J_rel
    algebraic[33] = constants[25]*(states[8]+states[9])*(states[3]-states[2])*states[5]
    # A10; J_tr
    algebraic[36] = (states[4]-states[3])/constants[26]
    # A4; CaJSR
    rates[3] = algebraic[29]*(algebraic[36]-algebraic[33])
    # A12; J_leak
    algebraic[40] = constants[27]*(states[4]-states[1])
    # A13; J_up
    algebraic[42] = (constants[29]*(np.power(states[1], 2.00000)))/(np.power(constants[30], 2.00000)+np.power(states[1], 2.00000))
    # A5; CaNSR
    rates[4] = ((algebraic[42]-algebraic[40])*constants[1])/constants[3]-(algebraic[36]*constants[2])/constants[3]
    # A42; CNa3
    algebraic[3] = 1.00000-(states[19]+states[20]+states[21]+states[24]+states[22]+states[23]+states[25]+states[26])
    # A51; alpha_Na11
    algebraic[13] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/17.0000)+0.200000*np.exp(-(states[0]+2.50000)/150.000))
    # A54; beta_Na11
    algebraic[35] = 0.191700*np.exp(-(states[0]+2.50000)/20.3000)
    # A52; alpha_Na12
    algebraic[26] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/15.0000)+0.230000*np.exp(-(states[0]+2.50000)/150.000))
    # A55; beta_Na12
    algebraic[37] = 0.200000*np.exp(-(states[0]-2.50000)/20.3000)
    # A57; alpha_Na3
    algebraic[41] = 7.00000e-07*np.exp(-(states[0]+7.00000)/7.70000)
    # A58; beta_Na3
    algebraic[43] = 0.00840000+2.00000e-05*(states[0]+7.00000)
    # A43; CNa2
    rates[21] = (algebraic[13]*algebraic[3]+algebraic[37]*states[20]+algebraic[41]*states[25])-(algebraic[35]*states[21]+algebraic[26]*states[21]+algebraic[43]*states[21])
    # A53; alpha_Na13
    algebraic[31] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/12.0000)+0.250000*np.exp(-(states[0]+2.50000)/150.000))
    # A56; beta_Na13
    algebraic[39] = 0.220000*np.exp(-(states[0]-7.50000)/20.3000)
    # A44; CNa1
    rates[20] = (algebraic[26]*states[21]+algebraic[39]*states[19]+algebraic[41]*states[24])-(algebraic[37]*states[20]+algebraic[31]*states[20]+algebraic[43]*states[20])
    # A49; ICNa2
    rates[25] = (algebraic[13]*states[26]+algebraic[37]*states[24]+algebraic[43]*states[21])-(algebraic[35]*states[25]+algebraic[26]*states[25]+algebraic[41]*states[25])
    # A50; ICNa3
    rates[26] = (algebraic[35]*states[25]+algebraic[43]*algebraic[3])-(algebraic[13]*states[26]+algebraic[41]*states[26])
    # A22; I_CaL
    algebraic[46] = constants[43]*states[11]*(states[0]-constants[42])
    # A7; B_ss
    algebraic[24] = np.power(1.00000+(constants[17]*constants[19])/(np.power(constants[19]+states[2], 2.00000)), -1.00000)
    # A11; J_xfer
    algebraic[38] = (states[2]-states[1])/constants[28]
    # A3; Cass
    rates[2] = algebraic[24]*((algebraic[33]*constants[2])/constants[4]-((algebraic[38]*constants[1])/constants[4]+(algebraic[46]*constants[5]*constants[0])/(2.00000*constants[4]*constants[11])))
    # A15; PRyR
    rates[5] = -0.0400000*states[5]-((0.100000*algebraic[46])/constants[33])*np.exp(-(np.power(states[0]-5.00000, 2.00000))/648.000)
    # A59; alpha_Na2
    algebraic[45] = 1.00000/(0.188495*np.exp(-(states[0]+7.00000)/16.6000)+0.393956)
    # A60; beta_Na2
    algebraic[47] = (algebraic[31]*algebraic[45]*algebraic[41])/(algebraic[39]*algebraic[43])
    # A45; ONa
    rates[19] = (algebraic[31]*states[20]+algebraic[47]*states[24])-(algebraic[39]*states[19]+algebraic[45]*states[19])
    # A61; alpha_Na4
    algebraic[49] = algebraic[45]/1000.00
    # A62; beta_Na4
    algebraic[51] = algebraic[41]
    # A46; IFNa
    rates[24] = (algebraic[45]*states[19]+algebraic[43]*states[20]+algebraic[51]*states[22]+algebraic[26]*states[25])-(algebraic[47]*states[24]+algebraic[41]*states[24]+algebraic[49]*states[24]+algebraic[37]*states[24])
    # A35; I_p(Ca)
    algebraic[48] = (constants[47]*(np.power(states[1], 2.00000)))/(np.power(constants[48], 2.00000)+np.power(states[1], 2.00000))
    # A36; I_NaCa
    algebraic[50] = ((((((constants[49]*1.00000)/(np.power(constants[50], 3.00000)+np.power(constants[7], 3.00000)))*1.00000)/(constants[51]+constants[8]))*1.00000)/(1.00000+constants[52]*np.exp(((constants[53]-1.00000)*states[0]*constants[11])/(constants[9]*constants[10]))))*(np.exp((constants[53]*states[0]*constants[11])/(constants[9]*constants[10]))*(np.power(states[18], 3.00000))*constants[8]-np.exp(((constants[53]-1.00000)*states[0]*constants[11])/(constants[9]*constants[10]))*(np.power(constants[7], 3.00000))*states[1])
    # A38; E_CaN
    algebraic[52] = ((constants[9]*constants[10])/(2.00000*constants[11]))*np.log(constants[8]/states[1])
    # A37; I_Cab
    algebraic[54] = constants[54]*(states[0]-algebraic[52])
    # A6; Bi
    algebraic[11] = np.power(1.00000+(constants[17]*constants[19])/(np.power(constants[19]+states[1], 2.00000)), -1.00000)
    # A14; J_trpn
    algebraic[44] = (constants[21]*states[1]*(constants[32]-states[7])+constants[23]*states[1]*(constants[31]-states[6]))-(constants[22]*states[7]+constants[24]*states[6])
    # A2; Cai
    rates[1] = algebraic[11]*((algebraic[40]+algebraic[38])-(algebraic[42]+algebraic[44]+(((algebraic[54]+algebraic[48])-2.00000*algebraic[50])*constants[5]*constants[0])/(2.00000*constants[1]*constants[11])))
    # A63; alpha_Na5
    algebraic[53] = algebraic[45]/95000.0
    # A64; beta_Na5
    algebraic[55] = algebraic[41]/50.0000
    # A47; I1Na
    rates[22] = (algebraic[49]*states[24]+algebraic[55]*states[23])-(algebraic[51]*states[22]+algebraic[53]*states[22])
    # A48; I2Na
    rates[23] = algebraic[53]*states[22]-algebraic[55]*states[23]
    # A41; E_Na
    algebraic[56] = ((constants[9]*constants[10])/constants[11])*np.log((0.900000*constants[7]+0.100000*constants[6])/(0.900000*states[18]+0.100000*states[27]))
    # A40; I_Na
    algebraic[57] = constants[55]*states[19]*(states[0]-algebraic[56])
    # A65; I_Nab
    algebraic[58] = constants[56]*(states[0]-algebraic[56])
    # A109; f_NaK
    algebraic[67] = 1.00000/(1.00000+0.124500*np.exp((-0.100000*states[0]*constants[11])/(constants[9]*constants[10]))+0.0365000*constants[71]*np.exp((-states[0]*constants[11])/(constants[9]*constants[10])))
    # A108; I_NaK
    algebraic[68] = (((constants[65]*algebraic[67]*1.00000)/(1.00000+np.power(constants[66]/states[18], 1.50000)))*constants[6])/(constants[6]+constants[67])
    # A39; Nai
    rates[18] = (-(algebraic[57]+algebraic[58]+3.00000*algebraic[68]+3.00000*algebraic[50])*constants[5]*constants[0])/(constants[1]*constants[11])
    # A68; E_K
    algebraic[59] = ((constants[9]*constants[10])/constants[11])*np.log(constants[6]/states[27])
    # A67; I_Kto,f  
    algebraic[60] = constants[57]*(np.power(states[28], 3.00000))*states[29]*(states[0]-algebraic[59])
    # A75; I_Kto,s
    algebraic[61] = constants[58]*states[30]*states[31]*(states[0]-algebraic[59])
    # A82; I_K1
    algebraic[62] = (((0.293800*constants[6])/(constants[6]+210.000))*(states[0]-algebraic[59]))/(1.00000+np.exp(0.0896000*(states[0]-algebraic[59])))
    # A83; I_Ks
    algebraic[63] = constants[59]*(np.power(states[32], 2.00000))*(states[0]-algebraic[59])
    # A87; I_kUR
    algebraic[64] = constants[60]*states[33]*states[34]*(states[0]-algebraic[59])
    # A92; I_Kss
    algebraic[65] = constants[61]*states[35]*states[36]*(states[0]-algebraic[59])
    # A96; I_Kr
    algebraic[66] = constants[62]*states[37]*(states[0]-((constants[9]*constants[10])/constants[11])*np.log((0.980000*constants[6]+0.0200000*constants[7])/(0.980000*states[27]+0.0200000*states[18])))
    # A66; Ki
    rates[27] = (-((algebraic[60]+algebraic[61]+algebraic[62]+algebraic[63]+algebraic[65]+algebraic[64]+algebraic[66])-2.00000*algebraic[68])*constants[5]*constants[0])/(constants[1]*constants[11])
    # I_stim
    algebraic[0] = custom_piecewise([np.greater_equal(voi , constants[12]) & np.less_equal(voi , constants[13]) & np.less_equal((voi-constants[12])-np.floor((voi-constants[12])/constants[14])*constants[14] , constants[15]), constants[16] , True, 0.00000])
    # A112; O_ClCa
    algebraic[69] = 0.200000/(1.00000+np.exp(-(states[0]-46.7000)/7.80000))
    # A111; i_ClCa
    algebraic[70] = ((constants[68]*algebraic[69]*states[1])/(states[1]+constants[70]))*(states[0]-constants[69])
    # A1; V
    rates[0] = 0
    
    return rates

def compute_algebraic(constants, states, voi, holding_p, holding_t, P1, P1_t, P2):
    algebraic = np.array([[0.0] * len(voi)] * size_algebraic)
    states = np.array(states)
    voi = np.array(voi)
    algebraic[1] = 1.00000-(states[10]+states[8]+states[9])
    algebraic[4] = 0.180640*np.exp(0.0357700*(states[0]+30.0000))
    algebraic[14] = 0.395600*np.exp(-0.0623700*(states[0]+30.0000))
    algebraic[5] = (0.000152000*np.exp(-(states[0]+13.5000)/7.00000))/(0.00670830*np.exp(-(states[0]+33.5000)/7.00000)+1.00000)
    algebraic[15] = (0.000950000*np.exp((states[0]+33.5000)/7.00000))/(0.0513350*np.exp((states[0]+33.5000)/7.00000)+1.00000)
    algebraic[6] = 1.00000/(1.00000+np.exp(-(states[0]+22.5000)/7.70000))
    algebraic[16] = 0.493000*np.exp(-0.0629000*states[0])+2.05800
    algebraic[7] = 1.00000/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    algebraic[17] = 270.000+1050.00/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    algebraic[8] = (4.81333e-06*(states[0]+26.5000))/(1.00000-np.exp(-0.128000*(states[0]+26.5000)))
    algebraic[18] = 9.53333e-05*np.exp(-0.0380000*(states[0]+26.5000))
    algebraic[19] = 0.493000*np.exp(-0.0629000*states[0])+2.05800
    algebraic[20] = 1200.00-170.000/(1.00000+np.exp((states[0]+45.2000)/5.70000))
    algebraic[21] = 39.3000*np.exp(-0.0862000*states[0])+13.1700
    algebraic[10] = 0.0137330*np.exp(0.0381980*states[0])
    algebraic[23] = 6.89000e-05*np.exp(-0.0417800*states[0])
    algebraic[2] = 1.00000-(states[11]+states[12]+states[13]+states[14]+states[15]+states[16]+states[17])
    algebraic[12] = (0.400000*np.exp((states[0]+12.0000)/10.0000)*((1.00000+0.700000*np.exp(-(np.power(states[0]+40.0000, 2.00000))/10.0000))-0.750000*np.exp(-(np.power(states[0]+20.0000, 2.00000))/400.000)))/(1.00000+0.120000*np.exp((states[0]+12.0000)/10.0000))
    algebraic[25] = 0.0500000*np.exp(-(states[0]+12.0000)/13.0000)
    algebraic[9] = 1.00000-(states[38]+states[39]+states[37]+states[40])
    algebraic[22] = 0.0223480*np.exp(0.0117600*states[0])
    algebraic[27] = 0.0470020*np.exp(-0.0631000*states[0])
    algebraic[28] = 0.0908210*np.exp(0.0233910*(states[0]+5.00000))
    algebraic[32] = 0.00649700*np.exp(-0.0326800*(states[0]+5.00000))
    algebraic[30] = (constants[45]*states[2])/(constants[46]+states[2])
    algebraic[34] = 13.0000*(1.00000-np.exp(-(np.power(states[0]+14.5000, 2.00000))/100.000))
    algebraic[29] = np.power(1.00000+(constants[18]*constants[20])/(np.power(constants[20]+states[3], 2.00000)), -1.00000)
    algebraic[33] = constants[25]*(states[8]+states[9])*(states[3]-states[2])*states[5]
    algebraic[36] = (states[4]-states[3])/constants[26]
    algebraic[40] = constants[27]*(states[4]-states[1])
    algebraic[42] = (constants[29]*(np.power(states[1], 2.00000)))/(np.power(constants[30], 2.00000)+np.power(states[1], 2.00000))
    algebraic[3] = 1.00000-(states[19]+states[20]+states[21]+states[24]+states[22]+states[23]+states[25]+states[26])
    algebraic[13] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/17.0000)+0.200000*np.exp(-(states[0]+2.50000)/150.000))
    algebraic[35] = 0.191700*np.exp(-(states[0]+2.50000)/20.3000)
    algebraic[26] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/15.0000)+0.230000*np.exp(-(states[0]+2.50000)/150.000))
    algebraic[37] = 0.200000*np.exp(-(states[0]-2.50000)/20.3000)
    algebraic[41] = 7.00000e-07*np.exp(-(states[0]+7.00000)/7.70000)
    algebraic[43] = 0.00840000+2.00000e-05*(states[0]+7.00000)
    algebraic[31] = 3.80200/(0.102700*np.exp(-(states[0]+2.50000)/12.0000)+0.250000*np.exp(-(states[0]+2.50000)/150.000))
    algebraic[39] = 0.220000*np.exp(-(states[0]-7.50000)/20.3000)
    algebraic[46] = constants[43]*states[11]*(states[0]-constants[42])
    algebraic[24] = np.power(1.00000+(constants[17]*constants[19])/(np.power(constants[19]+states[2], 2.00000)), -1.00000)
    algebraic[38] = (states[2]-states[1])/constants[28]
    algebraic[45] = 1.00000/(0.188495*np.exp(-(states[0]+7.00000)/16.6000)+0.393956)
    algebraic[47] = (algebraic[31]*algebraic[45]*algebraic[41])/(algebraic[39]*algebraic[43])
    algebraic[49] = algebraic[45]/1000.00
    algebraic[51] = algebraic[41]
    algebraic[48] = (constants[47]*(np.power(states[1], 2.00000)))/(np.power(constants[48], 2.00000)+np.power(states[1], 2.00000))
    algebraic[50] = ((((((constants[49]*1.00000)/(np.power(constants[50], 3.00000)+np.power(constants[7], 3.00000)))*1.00000)/(constants[51]+constants[8]))*1.00000)/(1.00000+constants[52]*np.exp(((constants[53]-1.00000)*states[0]*constants[11])/(constants[9]*constants[10]))))*(np.exp((constants[53]*states[0]*constants[11])/(constants[9]*constants[10]))*(np.power(states[18], 3.00000))*constants[8]-np.exp(((constants[53]-1.00000)*states[0]*constants[11])/(constants[9]*constants[10]))*(np.power(constants[7], 3.00000))*states[1])
    algebraic[52] = ((constants[9]*constants[10])/(2.00000*constants[11]))*np.log(constants[8]/states[1])
    algebraic[54] = constants[54]*(states[0]-algebraic[52])
    algebraic[11] = np.power(1.00000+(constants[17]*constants[19])/(np.power(constants[19]+states[1], 2.00000)), -1.00000)
    algebraic[44] = (constants[21]*states[1]*(constants[32]-states[7])+constants[23]*states[1]*(constants[31]-states[6]))-(constants[22]*states[7]+constants[24]*states[6])
    algebraic[53] = algebraic[45]/95000.0
    algebraic[55] = algebraic[41]/50.0000
    algebraic[56] = ((constants[9]*constants[10])/constants[11])*np.log((0.900000*constants[7]+0.100000*constants[6])/(0.900000*states[18]+0.100000*states[27]))
    algebraic[57] = constants[55]*states[19]*(states[0]-algebraic[56])
    algebraic[58] = constants[56]*(states[0]-algebraic[56])
    algebraic[67] = 1.00000/(1.00000+0.124500*np.exp((-0.100000*states[0]*constants[11])/(constants[9]*constants[10]))+0.0365000*constants[71]*np.exp((-states[0]*constants[11])/(constants[9]*constants[10])))
    algebraic[68] = (((constants[65]*algebraic[67]*1.00000)/(1.00000+np.power(constants[66]/states[18], 1.50000)))*constants[6])/(constants[6]+constants[67])
    algebraic[59] = ((constants[9]*constants[10])/constants[11])*np.log(constants[6]/states[27])
    algebraic[60] = constants[57]*(np.power(states[28], 3.00000))*states[29]*(states[0]-algebraic[59])
    algebraic[61] = constants[58]*states[30]*states[31]*(states[0]-algebraic[59])
    algebraic[62] = (((0.293800*constants[6])/(constants[6]+210.000))*(states[0]-algebraic[59]))/(1.00000+np.exp(0.0896000*(states[0]-algebraic[59])))
    algebraic[63] = constants[59]*(np.power(states[32], 2.00000))*(states[0]-algebraic[59])
    algebraic[64] = constants[60]*states[33]*states[34]*(states[0]-algebraic[59])
    algebraic[65] = constants[61]*states[35]*states[36]*(states[0]-algebraic[59])
    algebraic[66] = constants[62]*states[37]*(states[0]-((constants[9]*constants[10])/constants[11])*np.log((0.980000*constants[6]+0.0200000*constants[7])/(0.980000*states[27]+0.0200000*states[18])))
    algebraic[0] = custom_piecewise([np.greater_equal(voi , constants[12]) & np.less_equal(voi , constants[13]) & np.less_equal((voi-constants[12])-np.floor((voi-constants[12])/constants[14])*constants[14] , constants[15]), constants[16] , True, 0.00000])
    algebraic[69] = 0.200000/(1.00000+np.exp(-(states[0]-46.7000)/7.80000))
    algebraic[70] = ((constants[68]*algebraic[69]*states[1])/(states[1]+constants[70]))*(states[0]-constants[69])
    algebraic[71] = algebraic[71] = voltage_clamp(voi, holding_p, holding_t, P1, P1_t, P2)

    return algebraic

def solve_model(holding_p, holding_t, P1, P1_t, P2):
    """Solve model with ODE solver"""
    
    # Initialise constants and state variables
    (init_states, constants) = init_consts()
    
    # Set timespan to solve over
    voi = np.linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(compute_rates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants, holding_p, holding_t, P1, P1_t, P2)
    
    # Solve model
    states = np.array([[0.0] * len(voi)] * size_states)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = compute_algebraic(constants, states, voi,  holding_p, holding_t, P1, P1_t, P2)
    
    return (voi, states, algebraic)


if __name__ == "__main__":
    holding_p = -140  # mV
    holding_t = 10  # msec
    P1s = np.arange(-130, 50, 10)  # mV
    P1_t = 80  # msec
    P2 = -20  # mV
    P2_t = 100  # msec

    voi_container = list()
    S_container = list()
    A_container = list()

    for P1 in P1s:
        (voi, S, A) = solve_model(holding_p, holding_t, P1, P1_t, P2)
        voi_container.append(voi)
        S_container.append(S)
        A_container.append(A)
