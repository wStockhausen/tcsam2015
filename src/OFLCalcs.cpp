#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "OFLCalcs.hpp"

using namespace tcsam;

////////////////////////////////////////////////////////////////////////////////
//Catch_Calculator
////////////////////////////////////////////////////////////////////////////////
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 * @param npFsh - number of fisheries
 */
Catch_Calculator::Catch_Calculator(int npZBs, int npFsh){
    nZBs = npZBs;
    nFsh = npFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    maxF = 1.0;//default scale
    ct_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    rm_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dm_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
}

/**
 * Find max fishing mortality rate for target fishery (f=1)
 * 
 * Modifies: maxF
 * 
 * @return max fishing mortality rate
 */
double Catch_Calculator::findMaxTargetCaptureRate(void){
    maxF = wts::max(capF_fmsz(1));
    return maxF;
}

/**
 * Calculate catch abundance (ct_msz, rm_fmsz, dm_fmsz) and post-fisheries
 * abundance based on target fishing mortality rate 'dirF' and initial population 
 * abundance n_msz.
 * 
 * Modifies: ct_msz, rm_fmsz, dm_fmsz
 * 
 * @param dirF - directed fishery fishing mortality rate
 * @param n_msz - pre-fisheries population size
 * 
 * @return np_msz - post-fisheries population size
 * 
 */
d3_array Catch_Calculator::calcCatch(double dirF, d3_array& n_msz){
    if (debug) cout<<"starting Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    double ratF = dirF/maxF;//target fishery (f=1) scaling ratio
    dvector totFM(1,nZBs);
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            totFM += elem_prod(ratF*capF_fmsz(1,m,s),
                               retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(capF_fmsz(f,m,s),
                                   retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            np_msz(m,s) = elem_prod(exp(-totFM),n_msz(m,s));//survival after all fisheries
            ct_msz(m,s) = n_msz(m,s)-np_msz(m,s);           //total catch, all fisheries
            for (int f=1;f<=nFsh;f++){
                rm_fmsz(f,m,s) = elem_prod(elem_div(elem_prod(capF_fmsz(f,m,s),retF_fmsz(f,m,s)),totFM),
                                              ct_msz(m,s));//retained catch mortality
                dm_fmsz(f,m,s) = elem_prod(elem_div(elem_prod(capF_fmsz(f,m,s),hm_f(f)*(1.0-retF_fmsz(f,m,s))),totFM),
                                              ct_msz(m,s));//discard catch mortality
            }//f
        }//m
    }//s
    if (debug) cout<<"finished Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    return np_msz;
}

/**
 * Calculates probabilities of surviving fisheries, given directed 
 * fishing capture rate 'dirF'.
 * 
 * @param dirF - capture rate for target fishery (f=1)
 * 
 * @return d3_array of survival probabilities S_msz.
 */
d3_array Catch_Calculator::calcSurvival(double dirF){
    double ratF = dirF/maxF;//target fishery (f=1) scaling ratio
    dvector totFM(1,nZBs);
    d3_array S_msz(1,nMSs,1,nSCs,1,nZBs);
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            totFM += elem_prod(ratF*capF_fmsz(1,m,s),
                              retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(dirF*capF_fmsz(f,m,s),
                                  retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            S_msz(m,s) = exp(-totFM);                //survival of fisheries
        }//m
    }//s
    return S_msz;
}
/**
 * Set sex-specific fishery capture rates.
 * 
 * @param x - sex to select
 * @param capF_fxmsz - input capture rates
 */
void Catch_Calculator::setCaptureRates(int x, d5_array& capF_fxmsz){
    capF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    capF_fmsz.initialize();
    for (int f=1;f<=nFsh;f++){
        for (int m=1;m<=nMSs;m++){
            capF_fmsz(f,m) = capF_fxmsz(f,x,m);
        }//m
    }//f
}

/**
 * Set sex-specific retention functions.
 * 
 * @param x - sex to select
 * @param retF_fxmsz - input retention functions
 */
void Catch_Calculator::setRetentionFcns(int x, d5_array& retF_fxmsz){
    retF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    retF_fmsz.initialize();
    for (int f=1;f<=nFsh;f++){
        for (int m=1;m<=nMSs;m++){
            retF_fmsz(f,m) = retF_fxmsz(f,x,m);
        }//m
    }//f
}

/**
 * Set handling mortality rates, by fishery.
 * 
 * @param pHM_f - dvector of handling mortality rates
 */
void Catch_Calculator::setHandlingMortality(dvector& pHM_f){
    hm_f.allocate(1,nFsh);
    hm_f.initialize();
    hm_f = pHM_f;
}

////////////////////////////////////////////////////////////////////////////////
//OFL_Calculator
////////////////////////////////////////////////////////////////////////////////
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 */
OFL_Calculator::OFL_Calculator(int Tier, int npZBs){
    nSXs = tcsam::nSXs;
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    nZBs = npZBs;
    
    tier = Tier;
    
    //other constants
    alpha = 0.1; 
    beta  = 0.25;
}

/**
 * Calculate Bmsy, given longterm (male) recruitment.
 * 
 * @param R - longterm (male) recruitment (units??)
 * 
 * @return value for Bmsy (units??)
 */
double OFL_Calculator::calcBmsy(double R){
    if (debug) cout<<"starting double OFL_Calculator::calcBmsy(double R)"<<endl;
    double Bmsy = 0.0;
    switch (tier){
        case 3: {
            Bmsy = pT3C->calcBmsy(R);
            break;
        }
        default:
            cout<<"Tier "<<tier<<" OFL calculations not yet implemented!"<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
            break;
    }
    if (debug) cout<<"finished double OFL_Calculator::calcBmsy(double R)"<<endl;
    return Bmsy;
}

/**
 * Calculate Fmsy, given longterm (male) recruitment.
 * 
 * @param R - longterm (male) recruitment (units??)
 * 
 * @return value for Fmsy
 */
double OFL_Calculator::calcFmsy(double R){
    if (debug) cout<<"starting double OFL_Calculator::calcFmsy(double R)"<<endl;
    double Fmsy = 0.0;
    switch (tier){
        case 3: {
            Fmsy = pT3C->calcFmsy(R);
            break;
        }
        default:
            cout<<"Tier "<<tier<<" OFL calculations not yet implemented!"<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
            break;
    }
    if (debug) cout<<"finished double OFL_Calculator::calcFmsy(double R)"<<endl;
    return Fmsy;
}

/**
 * Use the harvest control rule to calculate an Fofl based on "current" MMB and Bmsy.
 * 
 * @param currMMB - "current" (projected to Feb. 15) MMB
 * @param Bmsy - Bmsy (based on Tier)
 * @param Fmsy - Fmsy (based on Tier)
 * 
 * @return Fofl rate using harvest control rule 
 */
double OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy){
    if (debug) cout<<"starting OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy)"<<endl;
    double Fofl = 0.0;
    double ratio  = currMMB/Bmsy;
    if (ratio < beta){Fofl = 0.0;} else
    if (ratio < 1.0){
        Fofl = Fmsy*(ratio-alpha)/(1-alpha);
    } else Fofl = Fmsy;
    
    if (debug) cout<<"finished OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy)"<<endl;
    return Fofl;    
}

/**
 * Calculate Fofl, given longterm (male) recruitment.
 * 
 * @param R - longterm (male) recruitment (units??)
 * 
 * @return value for Fofl
 */
double OFL_Calculator::calcFofl(double Bmsy, double Fmsy, d3_array& n_msz){
    if (debug) cout<<"starting double OFL_Calculator::calcFofl(double R)"<<endl;
    //start with guess for Fofl based on currMMB
    //TODO: need to set directed fishery F to Fmsy in pPrjM
    double currMMB = pPrjM->calcMatureBiomass(Fmsy,n_msz);
    double Fofl = calcHCR(currMMB,Bmsy,Fmsy);
    //now iterate until Fofl yields currMMB
    double Foflp = 0.0;
    double criF = 0.001;
    double delF = std::numeric_limits<double>::infinity();
    while (sfabs(delF)>criF){
        //calculate currMMB based on Fofl
        //TODO: need to set directed fishery F to Fofl in pPrjM
        currMMB = pPrjM->calcMatureBiomass(Fofl,n_msz);
        //update Fofl based on currMMB
        Foflp = Fofl;
        Fofl  = calcHCR(currMMB,Bmsy,Fmsy);
        //check convergence
        delF = Fofl - Foflp;
    }
    if (debug) cout<<"finished double OFL_Calculator::calcFofl(double R)"<<endl;
    return Fofl;
}

double OFL_Calculator::calcOFL(double Fofl, d4_array& n_xmsz){
    pPrjM->calcCatch(Fofl,n_xmsz(  MALE));
    
    pPrjF->calcCatch(Fofl,n_xmsz(FEMALE));
}

////////////////////////////////////////////////////////////////////////////////
//Tier3_Calculator
////////////////////////////////////////////////////////////////////////////////
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 * @param npFsh - number of fisheries
 */
Tier3_Calculator::Tier3_Calculator(int npZBs, int npFsh){
    nZBs = npZBs;
    nFsh = npFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    //set other constants
    XX = 0.35;//SPR rate for MSY calculations
}

//-------------------------------------------------------------------------------------
double Tier3_Calculator::calcMMB(d3_array& n_msz){
    if (debug) cout<<"starting double Tier3_Calculator::calcMMB()"<<endl;
    double mmb = 0.0; 
    for (int s=1;s<=nSCs;s++) mmb += n_msz(MATURE,s)*w_mz(MATURE);//dot product here
    if (debug) cout<<"finished double Tier3_Calculator::calcMMB()"<<endl;
    return mmb;
}

//-------------------------------------------------------------------------------------
//calculate sex-specific equilibrium size distribution
//non-differentiable version
d3_array Tier3_Calculator::calcEqNatZ(dvector& R_z,d3_array& S1_msz, 
                                      dmatrix& Th_sz, d3_array& T_szz, d3_array& S2_msz){
    if (debug) cout<<"starting d3_array calcEqNatZ()"<<endl;

    //the equilibrium solution
    d3_array n_msz(1,nMSs,1,nSCs,1,nZBs);
    
    //create an identity matrix
    dmatrix I = identity_matrix(1,nZBs);
    
    //--calc the state transition matrices
    int i = tcsam::IMMATURE; 
    int m = tcsam::MATURE;
    int n = tcsam::NEW_SHELL;
    int o = tcsam::OLD_SHELL;
    //immature new shell crab
    dmatrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dmatrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dmatrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|pre-molt size,new shell, molting)
    dmatrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dmatrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dmatrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dmatrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dmatrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|pre-molt size,old shell, molting)
    dmatrix Ph_io = identity_matrix(1,nZBs);//pr(molt|pre-molt size, old shell) [assumed all immature crab molt]
    dmatrix S1_io = wts::diag(S1_msz(i,o)); //pr(survival|size) for immature old shell crab before molting occurs
    //mature new shell crab
    dmatrix Tr_mn = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab undergoing terminal molt (same as non-terminal molt)
    dmatrix Tr_mo = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab undergoing terminal molt (same as non-terminal molt)
    dmatrix S2_mn = wts::diag(S2_msz(m,n)); //pr(survival|size) for mature new shell crab after molting occurs
    dmatrix S1_mn = wts::diag(S1_msz(m,n)); //pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    //mature old shell crab
    dmatrix S2_mo = wts::diag(S2_msz(m,o)); //pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    dmatrix S1_mo = wts::diag(S1_msz(m,o)); //pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    
    //full state transition matrices
    dmatrix lA = S2_in * Tr_in * (I-Th_in) * Ph_in * S1_in;//imm, new -> imm, new
    dmatrix lB = S2_in * Tr_io * (I-Th_io) * Ph_io * S1_io;//imm, old -> imm, new
    dmatrix lC = S2_io * (I-Ph_in) * S1_in;                //imm, new -> imm, old
    dmatrix lD = S2_io * (I-Ph_io) * S1_io;                //imm, old -> imm, old
    dmatrix lE = S2_mn * Tr_mn * Th_in * Ph_in * S1_in;    //imm, new -> mat, new (terminal molt)
    dmatrix lF = S2_mn * Tr_mo * Th_io * Ph_io * S1_io;    //imm, old -> mat, new (terminal molt)
    dmatrix lG = S2_mo * S1_mn;                            //mat, new -> mat, old
    dmatrix lH = S2_mo * S1_mo;                            //mat, old -> mat, old
    //--done calculating transition matrices
    
    //calculate inverses of matrix quantities
    dmatrix iM1 = inv(I - lD);
    dmatrix iM2 = inv(I - lA - lB * iM1 * lC);
    dmatrix iM3 = inv(I - lH);
    
    //the equilibrium solution is
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
        
    if (debug) cout<<"finished d3_array calcEqNatZ()"<<endl;
    return n_msz;
}

d3_array Tier3_Calculator::calcEqNatZF0(double R){
    if (debug) cout<<"starting d3_array calcEqNatZF0(double R)"<<endl;

    d3_array S1_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array S2_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n_msz(1,nMSs,1,nSCs,1,nZBs);  //equilibrium size distribution
    
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            S1_msz(m,s) = exp(-M_msz(m,s)*dtM);      //survival until molting/growth/mating
            S2_msz(m,s) = exp(-M_msz(m,s)*(1.0-dtM));//survival after molting/growth/mating
        }//m
    }//s
    dvector R_zp = R*R_z;
    n_msz = calcEqNatZ(R_zp, S1_msz, Th_sz, T_szz, S2_msz);
    
    if (debug) cout<<"finished d3_array Tier3_Calculator::calcEqNatZF0(double R)"<<endl;
    return n_msz;
}

d3_array Tier3_Calculator::calcEqNatZFM(double R, double dirF){
    if (debug) cout<<"starting d3_array Tier3_Calculator::calcEqNatZFM(double R)"<<endl;
    dvector totF(1,nZBs); //total F-at-size
    d3_array S1_msz(1,nMSs,1,nSCs,1,nZBs);//survival until molting/mating
    d3_array S2_msz(1,nMSs,1,nSCs,1,nZBs);//survival after molting/mating
    d3_array n_msz(1,nMSs,1,nSCs,1,nZBs); //equilibrium size distribution
    
    if (dtF<=dtM){
        //fisheries occur BEFORE molting/growth/maturity 
        for (int s=1;s<=nSCs;s++){
            for (int m=1;m<=nMSs;m++){ 
                S1_msz(m,s)  = exp(-M_msz(m,s)*dtF);          //survival until fisheries
                totF.initialize();
                totF += elem_prod(dirF*capF_fmsz(1,m,s),
                                  retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
                for (int f=2;f<=nFsh;f++)
                    totF += elem_prod(dirF*capF_fmsz(f,m,s),
                                      retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
                S1_msz(m,s) = elem_prod(S1_msz(m,s),exp(-totF));                //survival of fisheries
                S1_msz(m,s) = elem_prod(S1_msz(m,s),exp(-M_msz(m,s)*(dtM-dtF)));//survival from fisheries to molting/growth/mating
                S2_msz(m,s) = exp(-M_msz(m,s)*(1.0-dtM));                       //survival after molting/growth/mating
            }//m
        }//s
    } else {
        //fisheries occur AFTER molting/growth/maturity 
        for (int s=1;s<=nSCs;s++){
            for (int m=1;m<=nMSs;m++){ 
                S1_msz(m,s)  = exp(-M_msz(m,s)*dtM);          //survival until molting/growth/mating
                totF.initialize();
                totF += elem_prod(dirF*capF_fmsz(1,m,s),
                                  retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
                for (int f=2;f<=nFsh;f++)
                    totF += elem_prod(dirF*capF_fmsz(f,m,s),
                                      retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
                S2_msz(m,s) = exp(-M_msz(m,s)*(dtF-dtM));                       //survival from molting/growth/mating until fisheries
                S2_msz(m,s) = elem_prod(S2_msz(m,s),exp(-totF));                //survival of fisheries
                S2_msz(m,s) = elem_prod(S2_msz(m,s),exp(-M_msz(m,s)*(1.0-dtF)));//survival after fisheries
            }//m
        }//s
    }
    dvector R_zp = R*R_z;
    n_msz = calcEqNatZ(R_zp, S1_msz, Th_sz, T_szz, S2_msz);
    
    if (debug) cout<<"finished Ter3_Calculator::calcEqNatZFM(double R)"<<endl;
    return n_msz; 
}

/**
 * Calculate Fmsy for a Tier 3 stock.
 * For Tier 3 stocks, Fmsy = F35% is the directed fishing capture/mortality rate
 * that results in an equilibrium population size that yields Bmsy. For Tier 3 stocks,
 * Bmsy = 0.35*B100, where B100 is unfished MMB (mature male biomass at time of mating).
 * NOTE: By setting XX to some number other than
 * 
 * @param R - recruitment level for calculations
 * 
 * @return value for Fmsy = F35%.
 */
double Tier3_Calculator::calcFmsy(double R){
    if (debug) cout<<"starting Tier3_Calculator::calcFmsy(double R)"<<endl;
    //calculate unfished mmb
    double mmb100 = calcB100(R);
    
    //From initial guess for FXX, iterate to improve FXX    
    double dF   = 0.0001;//"delta" F for derivative calculation
    double FXX  = 0.5;   //initial guess for FXX
    double dFXX = 1.0;   //"delta" FXX to update (initial value is a dummy)
    double mmbp, mmb, mmbm, dMMBdF, XXp;
    int i=0;
    while ((i<20)&&(sfabs(dFXX)>0.00001)){
        mmbp = calcEqMMBatF(R,FXX+dF);
        mmb  = calcEqMMBatF(R,FXX);
        mmbm = calcEqMMBatF(R,FXX-dF);
        dMMBdF = 0.5*(mmbp-mmbm)/dF;//derivative of mmb wrto F
        XXp   = mmb/mmb100;           //ratio of mmb for current FXX relative to unfished 
        dFXX  = (XX*mmb100 - mmb)/dMMBdF;
        FXX  += dFXX;
        if (debug) {
            cout<<"--iteration = "<<++i<<endl;
            cout<<"----mmbXX = "<<mmb<<". mmb100 = "<<mmb100<<endl;
            cout<<"----XXp   = "<<XXp<<". dFXX   = "<<dFXX<<endl;
            cout<<"----FXX   = "<<FXX<<endl;
        }
    }//i loop
    
    if (debug) cout<<"finished Tier3_Calculator::calcFmsy(double R)"<<endl;
    return FXX;//Fmsy for Tier 3 stock
}

/**
 * Calculate B100, which is MMB for an unfished stock.
 * 
 * @param R - longterm (male) recruitment (units??)
 * 
 * @return value for B100 (units??)
 */
double Tier3_Calculator::calcB100(double R){
    if (debug) cout<<"finished d3_array calcB100(double R)"<<endl;
    d3_array n_msz = calcEqNatZF0(R);
    double B100 = calcMMB(n_msz);
    if (debug) cout<<"finished double Tier3_Calculator::calcB100(double R)"<<endl;
    return B100;
}

/**
 * Calculate Bmsy, given longterm (male) recruitment.
 * For tier 3 stocks, Bmsy = 0.35*B100 where B100 is MMB for an unfished stock.
 * 
 * @param R - longterm (male) recruitment (units??)
 * 
 * @return value for Bmsy (units??)
 */
double Tier3_Calculator::calcBmsy(double R){
    if (debug) cout<<"finished d3_array calcBmsy(double R)"<<endl;
    double Bmsy = XX*calcB100(R);
    if (debug) cout<<"finished double Tier3_Calculator::calcBmsy(double R)"<<endl;
    return Bmsy;
}

/**
 * Calculate equilibrium MMB when the directed fishing capture/mortality rate 
 * is dirF and longterm male recruitment is R.
 * 
 * @param R - longterm male recruitment
 * @param dirF - directed fishing capture/mortality rate
 * 
 * @return equilibrium MMB 
 */
double Tier3_Calculator::calcEqMMBatF(double R, double dirF){
    if (debug) cout<<"starting Tier3_Calculator::calcEqMMBatF(double R, double dirF)"<<endl;
    //calculate equilibrium size distribution
    d3_array n_msz = calcEqNatZFM(R,dirF);
    
    //calculate mmb
    double mmb = calcMMB(n_msz);
    if (debug) cout<<"finished Tier3_Calculator::calcEqMMBatF(double R, double dirF)"<<endl;
    return mmb;
}
////////////////////////////////////////////////////////////////////////////////
//PopProjector
////////////////////////////////////////////////////////////////////////////////
/**
 * Constructor
 * 
 * @param npZBs - number of size bins
 */
PopProjector::PopProjector(int npZBs){
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    nZBs = npZBs;
}

/**
 * Project sex-specific component of population ahead one year, WITHOUT recruitment.
 * Also calculates:
 *      spB - spawning biomass at mating time
 *      ct_msz - total fishing mortality (abundance)
 *      rm_fmsz - retained mortality, by fishery
 *      dm_fmsz - discard mortality, by fishery
 * 
 * @param n_msz - initial numbers-at-maturity state/shell condition/size
 * @param dirF - multiplier on fishing mortality rate in directed fishery
 * 
 * @return final numbers-at-maturity state/shell condition/size, WITHOUT recruitment
 */
d3_array PopProjector::project(double dirF, d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::project()"<<endl;
    d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n2_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n3_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n4_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n5_msz(1,nMSs,1,nSCs,1,nZBs);
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        //apply natural mortality BEFORE fisheries
        n1_msz = applyNM(dtF, n_msz);
        //apply fisheries
        n2_msz = applyFM(dirF, n1_msz);
        //apply natural mortality after fisheries but before molting/growth
        if (dtF==dtM){
            n3_msz = n2_msz;
        } else {
            n3_msz = applyNM(dtM-dtF, n2_msz);
        }        
        //apply molting/growth
        n4_msz = applyMG(n3_msz);
        //apply natural mortality AFTER molting/growth
        n5_msz = applyNM(1.0-dtM, n4_msz);
        
        //calculate spawning biomass from pre-molting/growth abundance
        n3_msz(tcsam::IMMATURE) = 0.0;//need to set immature abundance to 0
        spB = calcBiomass(n3_msz);
    } else { //fisheries occur AFTER molting/growth/maturity 
        //apply natural mortality BEFORE molting/growth
        n1_msz = applyNM(dtM, n_msz);        
        //apply molting/growth
        n2_msz = applyMG(n1_msz);
        //apply natural mortality after molting/growth but before fisheries
        if (dtF==dtM){
            n3_msz = n2_msz;
        } else {
            n3_msz = applyNM(dtF-dtM, n2_msz);
        }
        //apply fisheries
        n4_msz = applyFM(dirF, n3_msz);
        //apply natural mortality AFTER fisheries
        n5_msz = applyNM(1.0-dtF, n4_msz);
        
        //calculate spawning biomass from pre-molting/growth abundance
        n1_msz(tcsam::IMMATURE) = 0.0;//need to set immature abundance to 0
        spB = calcBiomass(n1_msz);
    }
    
    if (debug) cout<<"finished PopProjector::project()"<<endl;   
    return n5_msz;
}

/**
 * Calculate mature biomass-at-mating by projecting population forward to mating time.
 * 
 * @param n_msz - initial (July 1) numbers-at-maturity state/shell condition/size
 * 
 * @return MMB-at-mating
 */
double PopProjector::calcMatureBiomass(double dirF, d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::project()"<<endl;
    d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        //apply natural mortality BEFORE fisheries
        n1_msz = applyNM(dtF,n_msz);
        //apply fisheries
        d3_array n2_msz = applyFM(dirF, n1_msz);
        //apply natural mortality after fisheries but before molting/growth
        d3_array n3_msz;
        if (dtF==dtM){
            n3_msz = n2_msz;
        } else {
            n3_msz = applyNM(dtM-dtF,n2_msz);
        }
        
        //calculate MMB
        n3_msz(tcsam::IMMATURE) = 0.0;//need to set immature abundance to 0
        spB = calcBiomass(n3_msz);
    } else { //fisheries occur AFTER molting/growth/maturity 
        //apply natural mortality BEFORE molting/growth
        n1_msz = applyNM(dtM,n_msz);
        
        //calculate MMB
        n1_msz(tcsam::IMMATURE) = 0.0;//need to set immature abundance to 0
        spB = calcBiomass(n1_msz);
    }
    
    if (debug) cout<<"finished PopProjector::calcMMB()"<<endl;   
    return spB;
}

/** Apply natural mortality to population */
d3_array PopProjector::applyNM(double dt, d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::applyNM()"<<endl;
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            np_msz(m,s)  = elem_prod(exp(-M_msz(m,s)*dt),n_msz(m,s));
        }//m
    }//s
    if (debug) cout<<"finished PopProjector::applyNM()"<<endl;
    return np_msz;
}

/** Apply fishing mortality to population */
/**
 * 
 * Also calculates:
 *      spB - spawning biomass at mating time
 *      ct_msz - total fishing mortality (abundance)
 *      rm_fmsz - retained mortality, by fishery
 *      dm_fmsz - discard mortality, by fishery
 * @param dirF
 * @param n_msz
 * @return 
 */
d3_array PopProjector::applyFM(double dirF, d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::applyFM()"<<endl;
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    dvector totFM(1,nZBs);
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            totFM += elem_prod(dirF*capF_fmsz(1,m,s),
                               retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(capF_fmsz(f,m,s),
                                   retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            np_msz(m,s) = elem_prod(exp(-totFM),n_msz(m,s));//survival after all fisheries
            ct_msz(m,s) = n_msz(m,s)-np_msz(m,s);           //total catch, all fisheries
            for (int f=1;f<=nFsh;f++){
                rm_fmsz(f,m,s) = elem_prod(elem_div(elem_prod(capF_fmsz(f,m,s),retF_fmsz(f,m,s)),totFM),
                                              ct_msz(m,s));//retained catch mortality
                dm_fmsz(f,m,s) = elem_prod(elem_div(elem_prod(capF_fmsz(f,m,s),hm_f(f)*(1.0-retF_fmsz(f,m,s))),totFM),
                                              ct_msz(m,s));//discard catch mortality
            }//f
        }//m
    }//s
    if (debug) cout<<"finished PopProjector::applyFM()"<<endl;
    return np_msz;
}

/** Apply molting/growth to population */
d3_array PopProjector::applyMG(d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::applyMG()"<<endl;
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    np_msz(IMMATURE,NEW_SHELL) = T_szz(NEW_SHELL)*elem_prod(1.0-Th_sz(NEW_SHELL),n_msz(IMMATURE,NEW_SHELL));
    np_msz(IMMATURE,OLD_SHELL) = 0.0;
    np_msz(MATURE,NEW_SHELL)   = T_szz(NEW_SHELL)*elem_prod(    Th_sz(NEW_SHELL),n_msz(IMMATURE,NEW_SHELL));
    np_msz(MATURE,OLD_SHELL)   = n_msz(MATURE,NEW_SHELL)+n_msz(MATURE,OLD_SHELL);
    if (debug) cout<<"finished PopProjector::applyMG()"<<endl;
    return np_msz;
}

/** Calculate biomass */
double PopProjector::calcBiomass(d3_array& n_msz){
    if (debug) cout<<"starting PopProjector::calcBiomass()"<<endl;
    double bio = 0.0;
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            bio += w_mz(m)*n_msz(m,s);
        }//m
    }//s
    if (debug) cout<<"finished PopProjector::calcBiomass()"<<endl;
    return bio;
}
