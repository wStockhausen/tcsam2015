    #include <math.h>
    #include <time.h>
    #include <admodel.h>
    #include "TCSAM.hpp"
    
    adstring model  = "TCSAM2015";
    adstring modVer = "01.01"; 
    
    time_t start,finish;
    
    //model objects
    ModelConfiguration*  ptrMC; //ptr to model configuration object
    ModelParametersInfo* ptrMPI;//ptr to model parameters info object
    ModelOptions*        ptrMOs;//ptr to model options object
    ModelDatasets*       ptrMDS;//ptr to model datasets object
    ModelDatasets*       ptrSimMDS;//ptr to simulated model datasets object
    
    //file streams and filenames
    std::ofstream mcmc;        //stream for mcmc output
    
    //filenames
    adstring fnMCMC = "TCSAM2015.MCMC.R";
    adstring fnConfigFile;//configuration file
    adstring fnResultFile;//results file name 
    adstring fnPin;       //pin file
    
    //runtime flags
    int resample  = 0;//use resampling for initial parameter values (default is 0=false)
    int opModMode = 0;//run as operating model, no fitting
    int usePin    = 0;//flag to initialize parameter values using a pin file
    
    int iSeed = 999;//default random number generator seed
    random_number_generator rng(iSeed);//random number generator
    //debug flags
    int debugModelConfig     = 0;
    int debugModelDatasets   = 0;
    int debugModelParamsInfo = 0;
    int debugModelParams     = 0;
    
    int debugDATA_SECTION    = 0;
    int debugPARAMS_SECTION  = 0;
    int debugPRELIM_CALCS    = 0;
    int debugPROC_SECTION    = 0;
    int debugREPORT_SECTION  = 0;
    
    int showActiveParams = 0;    
    int debugRunModel    = 0;    
    int debugObjFun      = 0;
    
    int debugMCMC = 0;
    
    int dbgCalcProcs = 10;
    int dbgObjFun = 20;
    int dbgPriors = tcsam::dbgPriors;
    int dbgPopDy  = 70;
    int dbgApply  = 80;
    int dbgDevs   = 90;
    int dbgAll    = tcsam::dbgAll;
    
    int nSXs   = tcsam::nSXs;
    int MALE   = tcsam::MALE;
    int FEMALE = tcsam::FEMALE;
    int ANY_SX = tcsam::ANY_SEX;
    
    int nMSs     = tcsam::nMSs;
    int IMMATURE = tcsam::IMMATURE;
    int MATURE   = tcsam::MATURE;
    int ANY_MS   = tcsam::ANY_MATURITY;
    
    int nSCs      = tcsam::nSCs;
    int NEW_SHELL = tcsam::NEW_SHELL;
    int OLD_SHELL = tcsam::OLD_SHELL;
    int ANY_SC    = tcsam::ANY_SHELL;
    
    double smlVal = 0.00001;//small value to keep things > 0
    
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <TCSAM2015.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
    rpt::echo<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    rpt::echo<<"#Starting DATA_SECTION"<<endl;
    cout<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    cout<<"#Starting DATA_SECTION"<<endl;
    int on = 0;
    int flg = 0;
    //configFile
    fnConfigFile = "TCSAM2015_ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        flg=1;
    }
    //resultsFile
    fnResultFile = "TCSAM2015_ResultFile";//default results file name
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resultsFile"))>-1) {
        fnResultFile = ad_comm::argv[on+1];
        rpt::echo<<"#results file changed to '"<<fnResultFile<<"'"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
    }
    //opModMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-opModMode"))>-1) {
        opModMode=1;
        rpt::echo<<"#operating model mode turned ON"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        rpt::echo<<"#debugModelDatasets turned ON"<<endl;
        flg = 1;
    }
    //debugModelParamsInfo
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParamsInfo"))>-1) {
        debugModelParamsInfo=1;
        rpt::echo<<"#debugModelParamsInfo turned ON"<<endl;
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        if (on+1<argc) {
            iSeed=atoi(ad_comm::argv[on+1]);
            rng.reinitialize(iSeed);
            rpt::echo<<"#Resampling for initial parameter values using "<<iSeed<<endl;
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter iSeed for random number generator: ";
            cin>>iSeed;
        }
        rng.reinitialize(iSeed);
        rpt::echo<<"#Resampling for initial parameter values using "<<iSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        flg = 1;
    }
    //debugModelParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParams"))>-1) {
        debugModelParams=1;
        cout<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#debugModelParams turned ON"<<endl;
        flg = 1;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        rpt::echo<<"#debugDATA_SECTION turned ON"<<endl;
        flg = 1;
    }
    //debugPARAMS_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPARAMS_SECTION"))>-1) {
        debugPARAMS_SECTION=1;
        rpt::echo<<"#debugPARAMS_SECTION turned ON"<<endl;
        flg = 1;
    }
    //debugPRELIM_CALCS
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPRELIM_CALCS"))>-1) {
        debugPRELIM_CALCS=1;
        rpt::echo<<"debugPRELIM_CALCS turned ON"<<endl;
        flg = 1;
    }
    //debugPROC_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPROC_SECTION"))>-1) {
        debugPROC_SECTION=1;
        rpt::echo<<"#debugPROC_SECTION turned ON"<<endl;
        flg = 1;
    }
    //debugREPORT_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugREPORT_SECTION"))>-1) {
        debugREPORT_SECTION=1;
        rpt::echo<<"#debugREPORT_SECTION turned ON"<<endl;
        flg = 1;
    }
    //debugRunModel
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugRunModel")>-1) {
        debugRunModel=1;
        rpt::echo<<"#debugRunModel turned ON"<<endl;
        flg = 1;
    }
    //debugObjFun
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugObjFun")>-1) {
        debugObjFun=1;
        rpt::echo<<"#debugObjFun turned ON"<<endl;
        flg = 1;
    }
    //showActiveParams
    if (option_match(ad_comm::argc,ad_comm::argv,"-showActiveParams")>-1) {
        showActiveParams=1;
        rpt::echo<<"#showActiveParams turned ON"<<endl;
        flg = 1;
    }
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
        flg = 1;
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading configuration file '"<<fnConfigFile<<"'"<<endl;
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    rpt::echo<<"#----finished model configuration---"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelConfiguration-----------------"<<endl;
        cout<<(*ptrMC);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model configuration---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    mxYrp1 = mxYr+1;
    nFsh   = ptrMC->nFsh;
    nSrv   = ptrMC->nSrv;
    nZBs   = ptrMC->nZBs;
  zBs.allocate(1,nZBs);
zBs  = ptrMC->zMidPts;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading parameters info file '"<<ptrMC->fnMPI<<"'"<<endl;
    if (debugModelParamsInfo) ModelParametersInfo::debug=1;
    ptrMPI = new ModelParametersInfo(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMPI);
    ptrMPI->read(*(ad_comm::global_datafile));
    if (debugModelParamsInfo) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelParamsInfo;
        if (debugModelParamsInfo<0) exit(1);
        ModelParametersInfo::debug=debugModelParamsInfo;
    }
    rpt::echo<<"#----finished model parameters info---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelParametersInfo-----------------"<<endl;
        cout<<(*ptrMPI);
        cout<<"#----finished model parameters info---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file '"<<ptrMC->fnMDS<<"'"<<endl;
    if (debugModelDatasets) {
        BioData::debug=1;
        FisheryData::debug=1;
        SurveyData::debug=1;
    }
    ptrMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrMDS->read(*(ad_comm::global_datafile));
    if (debugModelDatasets) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelDatasets;
        if (debugModelDatasets<0) exit(1);
        ModelDatasets::debug=debugModelDatasets;
        BioData::debug=debugModelDatasets;
        SurveyData::debug=debugModelDatasets;
    }
    rpt::echo<<"#----finished model datasets---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelDatasets-----------------"<<endl;
        cout<<(*ptrMDS);
        cout<<"#----finished model datasets---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file again to create SimMDS object '"<<ptrMC->fnMDS<<"'"<<endl;
    ptrSimMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrSimMDS->read(*(ad_comm::global_datafile));
    rpt::echo<<"#finished SimMDS object"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading model options file '"<<ptrMC->fnMOs<<"'"<<endl;
    if (debugModelParamsInfo) ModelOptions::debug=1;
    ptrMOs = new ModelOptions(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMOs);
    ptrMOs->read(*(ad_comm::global_datafile));
    if (debugModelParamsInfo) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelParamsInfo;
        if (debugModelParamsInfo<0) exit(1);
        ModelOptions::debug=debugModelParamsInfo;
    }
    rpt::echo<<"#----finished model options---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelOptions-----------------"<<endl;
        cout<<(*ptrMOs);
        cout<<"#----finished model options---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
  mapD2MFsh.allocate(1,nFsh);
  mapM2DFsh.allocate(1,nFsh);
    {
     int idx;
     for (int f=1;f<=nFsh;f++){
         idx = wts::which(ptrMDS->ppFsh[f-1]->name,ptrMC->lblsFsh);
         mapD2MFsh(f)   = idx;//map from fishery data object f to model fishery idx
         mapM2DFsh(idx) = f;  //map from model fishery idx to fishery data object f
     }
     rpt::echo<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
     cout<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
    }
  mapD2MSrv.allocate(1,nSrv);
  mapM2DSrv.allocate(1,nSrv);
    {
     int idx;
     for (int v=1;v<=nSrv;v++){
         idx = wts::which(ptrMDS->ppSrv[v-1]->name,ptrMC->lblsSrv);
         mapD2MSrv(v)   = idx;//map from survey data object v to model survey idx
         mapM2DSrv(idx) = v;  //map from model survey idx to survey data object v
     }
     rpt::echo<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
     cout<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
    }
  phsLnR.allocate();
  lbLnR.allocate();
  ubLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnR,npLnR,lbLnR,ubLnR,phsLnR,rpt::echo);
  phsLnRCV.allocate();
  lbLnRCV.allocate();
  ubLnRCV.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRCV,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV,rpt::echo);
  phsLgtRX.allocate();
  lbLgtRX.allocate();
  ubLgtRX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLgtRX,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX,rpt::echo);
  phsLnRa.allocate();
  lbLnRa.allocate();
  ubLnRa.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRa,npLnRa,lbLnRa,ubLnRa,phsLnRa,rpt::echo);
  phsLnRb.allocate();
  lbLnRb.allocate();
  ubLnRb.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRb,npLnRb,lbLnRb,ubLnRb,phsLnRb,rpt::echo);
  mniDevsLnR.allocate();
  mxiDevsLnR.allocate();
  idxsDevsLnR.allocate();
  lbDevsLnR.allocate();
  ubDevsLnR.allocate();
  phsDevsLnR.allocate();
tcsam::setParameterInfo(ptrMPI->ptrRec->pDevsLnR,npDevsLnR,mniDevsLnR,mxiDevsLnR,idxsDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,rpt::echo);
  phsLnM.allocate();
  lbLnM.allocate();
  ubLnM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnM,npLnM,lbLnM,ubLnM,phsLnM,rpt::echo);
  phsLnDMT.allocate();
  lbLnDMT.allocate();
  ubLnDMT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMT,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT,rpt::echo);
  phsLnDMX.allocate();
  lbLnDMX.allocate();
  ubLnDMX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMX,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX,rpt::echo);
  phsLnDMM.allocate();
  lbLnDMM.allocate();
  ubLnDMM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMM,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM,rpt::echo);
  phsLnDMXM.allocate();
  lbLnDMXM.allocate();
  ubLnDMXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMXM,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM,rpt::echo);
zMref = ptrMPI->ptrNM->zRef;
  mniLgtPrMat.allocate();
  mxiLgtPrMat.allocate();
  idxsLgtPrMat.allocate();
  lbLgtPrMat.allocate();
  ubLgtPrMat.allocate();
  phsLgtPrMat.allocate();
tcsam::setParameterInfo(ptrMPI->ptrMat->pLgtPrMat,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,idxsLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,rpt::echo);
  phsLnGrA.allocate();
  lbLnGrA.allocate();
  ubLnGrA.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrA,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA,rpt::echo);
  phsLnGrB.allocate();
  lbLnGrB.allocate();
  ubLnGrB.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrB,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB,rpt::echo);
  phsLnGrBeta.allocate();
  lbLnGrBeta.allocate();
  ubLnGrBeta.allocate();
tcsam::setParameterInfo(ptrMPI->ptrGr->pLnGrBeta,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta,rpt::echo);
nSel = ptrMPI->ptrSel->nPCs;//number of selectivity functions defined
  phsS1.allocate();
  lbS1.allocate();
  ubS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS1,npS1,lbS1,ubS1,phsS1,rpt::echo);
  phsS2.allocate();
  lbS2.allocate();
  ubS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS2,npS2,lbS2,ubS2,phsS2,rpt::echo);
  phsS3.allocate();
  lbS3.allocate();
  ubS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS3,npS3,lbS3,ubS3,phsS3,rpt::echo);
  phsS4.allocate();
  lbS4.allocate();
  ubS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS4,npS4,lbS4,ubS4,phsS4,rpt::echo);
  phsS5.allocate();
  lbS5.allocate();
  ubS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS5,npS5,lbS5,ubS5,phsS5,rpt::echo);
  phsS6.allocate();
  lbS6.allocate();
  ubS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pS6,npS6,lbS6,ubS6,phsS6,rpt::echo);
  mniDevsS1.allocate();
  mxiDevsS1.allocate();
  idxsDevsS1.allocate();
  lbDevsS1.allocate();
  ubDevsS1.allocate();
  phsDevsS1.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS1,npDevsS1,mniDevsS1,mxiDevsS1,idxsDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,rpt::echo);
  mniDevsS2.allocate();
  mxiDevsS2.allocate();
  idxsDevsS2.allocate();
  lbDevsS2.allocate();
  ubDevsS2.allocate();
  phsDevsS2.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS2,npDevsS2,mniDevsS2,mxiDevsS2,idxsDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,rpt::echo);
  mniDevsS3.allocate();
  mxiDevsS3.allocate();
  idxsDevsS3.allocate();
  lbDevsS3.allocate();
  ubDevsS3.allocate();
  phsDevsS3.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS3,npDevsS3,mniDevsS3,mxiDevsS3,idxsDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,rpt::echo);
  mniDevsS4.allocate();
  mxiDevsS4.allocate();
  idxsDevsS4.allocate();
  lbDevsS4.allocate();
  ubDevsS4.allocate();
  phsDevsS4.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS4,npDevsS4,mniDevsS4,mxiDevsS4,idxsDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,rpt::echo);
  mniDevsS5.allocate();
  mxiDevsS5.allocate();
  idxsDevsS5.allocate();
  lbDevsS5.allocate();
  ubDevsS5.allocate();
  phsDevsS5.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS5,npDevsS5,mniDevsS5,mxiDevsS5,idxsDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,rpt::echo);
  mniDevsS6.allocate();
  mxiDevsS6.allocate();
  idxsDevsS6.allocate();
  lbDevsS6.allocate();
  ubDevsS6.allocate();
  phsDevsS6.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS6,npDevsS6,mniDevsS6,mxiDevsS6,idxsDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,rpt::echo);
  phsHM.allocate();
  lbHM.allocate();
  ubHM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pHM,npHM,lbHM,ubHM,phsHM,rpt::echo);
  phsLnC.allocate();
  lbLnC.allocate();
  ubLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnC,npLnC,lbLnC,ubLnC,phsLnC,rpt::echo);
  phsLnDCT.allocate();
  lbLnDCT.allocate();
  ubLnDCT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCT,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT,rpt::echo);
  phsLnDCX.allocate();
  lbLnDCX.allocate();
  ubLnDCX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCX,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX,rpt::echo);
  phsLnDCM.allocate();
  lbLnDCM.allocate();
  ubLnDCM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCM,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM,rpt::echo);
  phsLnDCXM.allocate();
  lbLnDCXM.allocate();
  ubLnDCXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCXM,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM,rpt::echo);
  mniDevsLnC.allocate();
  mxiDevsLnC.allocate();
  idxsDevsLnC.allocate();
  lbDevsLnC.allocate();
  ubDevsLnC.allocate();
  phsDevsLnC.allocate();
tcsam::setParameterInfo(ptrMPI->ptrFsh->pDevsLnC,npDevsLnC,mniDevsLnC,mxiDevsLnC,idxsDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,rpt::echo);
  phsLnQ.allocate();
  lbLnQ.allocate();
  ubLnQ.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnQ,npLnQ,lbLnQ,ubLnQ,phsLnQ,rpt::echo);
  phsLnDQT.allocate();
  lbLnDQT.allocate();
  ubLnDQT.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQT,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT,rpt::echo);
  phsLnDQX.allocate();
  lbLnDQX.allocate();
  ubLnDQX.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQX,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX,rpt::echo);
  phsLnDQM.allocate();
  lbLnDQM.allocate();
  ubLnDQM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQM,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM,rpt::echo);
  phsLnDQXM.allocate();
  lbLnDQXM.allocate();
  ubLnDQXM.allocate();
tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQXM,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM,rpt::echo);
  dtF.allocate(mnYr,mxYr);
dtF = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
  dtM.allocate(mnYr,mxYr);
dtM = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
  optsFcAvg.allocate(1,nFsh);
optsFcAvg = ptrMOs->optsFcAvg;
  avgEff.allocate(1,nFsh);
npcRec = ptrMPI->ptrRec->nPCs;
npcNM = ptrMPI->ptrNM->nPCs;
npcMat = ptrMPI->ptrMat->nPCs;
npcGr = ptrMPI->ptrGr->nPCs;
npcSel = ptrMPI->ptrSel->nPCs;
npcFsh = ptrMPI->ptrFsh->nPCs;
npcSrv = ptrMPI->ptrSrv->nPCs;
    rpt::echo<<"#finished DATA_SECTION"<<endl;
    cout<<"#finished DATA_SECTION"<<endl;
}

void model_parameters::initializationfunction(void)
{
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
rpt::echo<<"#Starting PARAMETER_SECTION"<<endl;
cout<<"#Starting PARAMETER_SECTION"<<endl;
  pLnR.allocate(1,npLnR,lbLnR,ubLnR,phsLnR,"pLnR");
  pLnRCV.allocate(1,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV,"pLnRCV");
  pLgtRX.allocate(1,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX,"pLgtRX");
  pLnRa.allocate(1,npLnRa,lbLnRa,ubLnRa,phsLnRa,"pLnRa");
  pLnRb.allocate(1,npLnRb,lbLnRb,ubLnRb,phsLnRb,"pLnRb");
  pDevsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,"pDevsLnR");
  devsLnR.allocate(1,npDevsLnR,mniDevsLnR,mxiDevsLnR+1,"devsLnR");
  #ifndef NO_AD_INITIALIZE
    devsLnR.initialize();
  #endif
  pLnM.allocate(1,npLnM,lbLnM,ubLnM,phsLnM,"pLnM");
  pLnDMT.allocate(1,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT,"pLnDMT");
  pLnDMX.allocate(1,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX,"pLnDMX");
  pLnDMM.allocate(1,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM,"pLnDMM");
  pLnDMXM.allocate(1,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM,"pLnDMXM");
  pLnGrA.allocate(1,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA,"pLnGrA");
  pLnGrB.allocate(1,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB,"pLnGrB");
  pLnGrBeta.allocate(1,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta,"pLnGrBeta");
  pLgtPrMat.allocate(1,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,"pLgtPrMat");
  pS1.allocate(1,npS1,lbS1,ubS1,phsS1,"pS1");
  pS2.allocate(1,npS2,lbS2,ubS2,phsS2,"pS2");
  pS3.allocate(1,npS3,lbS3,ubS3,phsS3,"pS3");
  pS4.allocate(1,npS4,lbS4,ubS4,phsS4,"pS4");
  pS5.allocate(1,npS5,lbS5,ubS5,phsS5,"pS5");
  pS6.allocate(1,npS6,lbS6,ubS6,phsS6,"pS6");
  pDevsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,"pDevsS1");
  pDevsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,"pDevsS2");
  pDevsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,"pDevsS3");
  pDevsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,"pDevsS4");
  pDevsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,"pDevsS5");
  pDevsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,"pDevsS6");
  devsS1.allocate(1,npDevsS1,mniDevsS1,mxiDevsS1+1,"devsS1");
  #ifndef NO_AD_INITIALIZE
    devsS1.initialize();
  #endif
  devsS2.allocate(1,npDevsS2,mniDevsS2,mxiDevsS2+1,"devsS2");
  #ifndef NO_AD_INITIALIZE
    devsS2.initialize();
  #endif
  devsS3.allocate(1,npDevsS3,mniDevsS3,mxiDevsS3+1,"devsS3");
  #ifndef NO_AD_INITIALIZE
    devsS3.initialize();
  #endif
  devsS4.allocate(1,npDevsS4,mniDevsS4,mxiDevsS4+1,"devsS4");
  #ifndef NO_AD_INITIALIZE
    devsS4.initialize();
  #endif
  devsS5.allocate(1,npDevsS5,mniDevsS5,mxiDevsS5+1,"devsS5");
  #ifndef NO_AD_INITIALIZE
    devsS5.initialize();
  #endif
  devsS6.allocate(1,npDevsS6,mniDevsS6,mxiDevsS6+1,"devsS6");
  #ifndef NO_AD_INITIALIZE
    devsS6.initialize();
  #endif
  pHM.allocate(1,npHM,lbHM,ubHM,phsHM,"pHM");
  pLnC.allocate(1,npLnC,lbLnC,ubLnC,phsLnC,"pLnC");
  pLnDCT.allocate(1,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT,"pLnDCT");
  pLnDCX.allocate(1,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX,"pLnDCX");
  pLnDCM.allocate(1,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM,"pLnDCM");
  pLnDCXM.allocate(1,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM,"pLnDCXM");
  pDevsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,"pDevsLnC");
  devsLnC.allocate(1,npDevsLnC,mniDevsLnC,mxiDevsLnC+1,"devsLnC");
  #ifndef NO_AD_INITIALIZE
    devsLnC.initialize();
  #endif
  pLnQ.allocate(1,npLnQ,lbLnQ,ubLnQ,phsLnQ,"pLnQ");
  pLnDQT.allocate(1,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT,"pLnDQT");
  pLnDQX.allocate(1,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX,"pLnDQX");
  pLnDQM.allocate(1,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM,"pLnDQM");
  pLnDQXM.allocate(1,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM,"pLnDQXM");
  objFun.allocate("objFun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  spb_yx.allocate(mnYr,mxYr,1,nSXs,"spb_yx");
  #ifndef NO_AD_INITIALIZE
    spb_yx.initialize();
  #endif
  n_yxmsz.allocate(mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_yxmsz");
  #ifndef NO_AD_INITIALIZE
    n_yxmsz.initialize();
  #endif
  R_y.allocate(mnYr,mxYr,"R_y");
  #ifndef NO_AD_INITIALIZE
    R_y.initialize();
  #endif
  Rx_c.allocate(1,npcRec,"Rx_c");
  #ifndef NO_AD_INITIALIZE
    Rx_c.initialize();
  #endif
  R_yx.allocate(mnYr,mxYr,1,nSXs,"R_yx");
  #ifndef NO_AD_INITIALIZE
    R_yx.initialize();
  #endif
  R_cz.allocate(1,npcRec,1,nZBs,"R_cz");
  #ifndef NO_AD_INITIALIZE
    R_cz.initialize();
  #endif
  R_yz.allocate(mnYr,mxYr,1,nZBs,"R_yz");
  #ifndef NO_AD_INITIALIZE
    R_yz.initialize();
  #endif
  zscrDevsLnR.allocate(mnYr,mxYr,"zscrDevsLnR");
  #ifndef NO_AD_INITIALIZE
    zscrDevsLnR.initialize();
  #endif
  M_cxm.allocate(1,npcNM,1,nSXs,1,nMSs,"M_cxm");
  #ifndef NO_AD_INITIALIZE
    M_cxm.initialize();
  #endif
  M_yxmz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nZBs,"M_yxmz");
  #ifndef NO_AD_INITIALIZE
    M_yxmz.initialize();
  #endif
  prMat_cz.allocate(1,npcMat,1,nZBs,"prMat_cz");
  #ifndef NO_AD_INITIALIZE
    prMat_cz.initialize();
  #endif
  prMat_yxz.allocate(mnYr,mxYr,1,nSXs,1,nZBs,"prMat_yxz");
  #ifndef NO_AD_INITIALIZE
    prMat_yxz.initialize();
  #endif
  prGr_czz.allocate(1,npcGr,1,nZBs,1,nZBs,"prGr_czz");
  #ifndef NO_AD_INITIALIZE
    prGr_czz.initialize();
  #endif
  prGr_yxmzz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nZBs,1,nZBs,"prGr_yxmzz");
  #ifndef NO_AD_INITIALIZE
    prGr_yxmzz.initialize();
  #endif
  sel_cz.allocate(1,npcSel,1,nZBs,"sel_cz");
  #ifndef NO_AD_INITIALIZE
    sel_cz.initialize();
  #endif
  sel_iyz.allocate(1,nSel,mnYr,mxYr+1,1,nZBs,"sel_iyz");
  #ifndef NO_AD_INITIALIZE
    sel_iyz.initialize();
  #endif
  avgFc.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,"avgFc");
  #ifndef NO_AD_INITIALIZE
    avgFc.initialize();
  #endif
  avgRatioFc2Eff.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,"avgRatioFc2Eff");
  #ifndef NO_AD_INITIALIZE
    avgRatioFc2Eff.initialize();
  #endif
  Fhm_fy.allocate(1,nFsh,mnYr,mxYr,"Fhm_fy");
  #ifndef NO_AD_INITIALIZE
    Fhm_fy.initialize();
  #endif
  Fc_fyxms.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,"Fc_fyxms");
  #ifndef NO_AD_INITIALIZE
    Fc_fyxms.initialize();
  #endif
  Fc_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"Fc_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    Fc_fyxmsz.initialize();
  #endif
  Fr_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"Fr_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    Fr_fyxmsz.initialize();
  #endif
  Fd_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"Fd_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    Fd_fyxmsz.initialize();
  #endif
  Ft_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"Ft_yxmsz");
  #ifndef NO_AD_INITIALIZE
    Ft_yxmsz.initialize();
  #endif
  c_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"c_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    c_fyxmsz.initialize();
  #endif
  r_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"r_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    r_fyxmsz.initialize();
  #endif
  d_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"d_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    d_fyxmsz.initialize();
  #endif
  spb_vyx.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,"spb_vyx");
  #ifndef NO_AD_INITIALIZE
    spb_vyx.initialize();
  #endif
  q_vyxms.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,"q_vyxms");
  #ifndef NO_AD_INITIALIZE
    q_vyxms.initialize();
  #endif
  q_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"q_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    q_vyxmsz.initialize();
  #endif
  n_vyxmsz.allocate(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"n_vyxmsz");
  #ifndef NO_AD_INITIALIZE
    n_vyxmsz.initialize();
  #endif
  fPenRecDevs.allocate(1,npDevsLnR,"fPenRecDevs");
  #ifndef NO_AD_INITIALIZE
    fPenRecDevs.initialize();
  #endif
  fPenSmoothLgtPrMat.allocate(1,npLgtPrMat,"fPenSmoothLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenSmoothLgtPrMat.initialize();
  #endif
  fPenNonDecLgtPrMat.allocate(1,npLgtPrMat,"fPenNonDecLgtPrMat");
  #ifndef NO_AD_INITIALIZE
    fPenNonDecLgtPrMat.initialize();
  #endif
  fPenDevsS1.allocate(1,npDevsS1,"fPenDevsS1");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS1.initialize();
  #endif
  fPenDevsS2.allocate(1,npDevsS2,"fPenDevsS2");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS2.initialize();
  #endif
  fPenDevsS3.allocate(1,npDevsS3,"fPenDevsS3");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS3.initialize();
  #endif
  fPenDevsS4.allocate(1,npDevsS4,"fPenDevsS4");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS4.initialize();
  #endif
  fPenDevsS5.allocate(1,npDevsS5,"fPenDevsS5");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS5.initialize();
  #endif
  fPenDevsS6.allocate(1,npDevsS6,"fPenDevsS6");
  #ifndef NO_AD_INITIALIZE
    fPenDevsS6.initialize();
  #endif
  fPenDevsLnC.allocate(1,npDevsLnC,"fPenDevsLnC");
  #ifndef NO_AD_INITIALIZE
    fPenDevsLnC.initialize();
  #endif
  nllRecDevs.allocate("nllRecDevs");
  #ifndef NO_AD_INITIALIZE
  nllRecDevs.initialize();
  #endif
  sdrLnR_y.allocate(mnYr,mxYr,"sdrLnR_y");
  sdrSpB_xy.allocate(1,nSXs,mnYr+5,mxYr,"sdrSpB_xy");
cout<<"#finished PARAMETER_SECTION"<<endl;
rpt::echo<<"#finished PARAMETER_SECTION"<<endl;
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
    rpt::echo<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    cout<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    int debug=1;
        
    {cout<<"writing data to R"<<endl;
     ofstream echo1; echo1.open("ModelData.R", ios::trunc);
     ReportToR_Data(echo1,0,cout);
    }
    
    {cout<<"writing parameters info to R"<<endl;
     ofstream echo1; echo1.open("ModelParametersInfo.R", ios::trunc);
     ptrMPI->writeToR(echo1);
    }
    
    //calculate average effort for model fisheries over specified time periods
    avgEff = 0.0;
    for (int f=1;f<=nFsh;f++){//fishery data object
        if (ptrMDS->ppFsh[f-1]->ptrEff){
            IndexRange* pir = ptrMDS->ppFsh[f-1]->ptrEff->ptrAvgIR;
            int fm = mapD2MFsh(f);//index of corresponding model fishery
            int mny = max(mnYr,pir->getMin());
            int mxy = min(mxYr,pir->getMax());
            avgEff(fm) = mean(ptrMDS->ppFsh[f-1]->ptrEff->eff_y(mny,mxy));
        }
    }
    if (debug) rpt::echo<<"avgEff = "<<avgEff<<endl;
    
    //set initial values for all parameters
    //recruitment parameters
    setInitVals(ptrMPI->ptrRec->pLnR,    pLnR,    0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRCV,  pLnRCV,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLgtRX,  pLgtRX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRa,   pLnRa,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRb,   pLnRb,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,0,rpt::echo);
    //natural mortality parameters
    setInitVals(ptrMPI->ptrNM->pLnM,   pLnM,   0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMT, pLnDMT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMX, pLnDMX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMM, pLnDMM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMXM,pLnDMXM,0,rpt::echo);
    //growth parameters
    setInitVals(ptrMPI->ptrGr->pLnGrA,   pLnGrA,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGr->pLnGrB,   pLnGrB,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGr->pLnGrBeta,pLnGrBeta,0,rpt::echo);
    //maturity parameters
    setInitVals(ptrMPI->ptrMat->pLgtPrMat,pLgtPrMat,0,rpt::echo);
    //selectivity parameters
    setInitVals(ptrMPI->ptrSel->pS1, pS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS2, pS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS3, pS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS4, pS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS5, pS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS6, pS6,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS1, pDevsS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS2, pDevsS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS3, pDevsS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS4, pDevsS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS5, pDevsS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS6, pDevsS6,0,rpt::echo);
    //fully-selected fishing capture rate parameters
    setInitVals(ptrMPI->ptrFsh->pHM,     pHM,     0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnC,    pLnC,    0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCT,  pLnDCT,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCX,  pLnDCX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCM,  pLnDCM,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCXM, pLnDCXM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,0,rpt::echo);
    //survey catchability parameters
    setInitVals(ptrMPI->ptrSrv->pLnQ,   pLnQ,   0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQT, pLnDQT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQX, pLnDQX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQM, pLnDQM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQXM,pLnDQXM,0,rpt::echo);
    cout<<"testing setAllDevs()"<<endl;
    setAllDevs(0,rpt::echo);
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")<0) {
        cout<<"testing calcRecruitment():"<<endl;
        calcRecruitment(0,rpt::echo);
        cout<<"testing calcNatMort():"<<endl;
        calcNatMort(0,rpt::echo);
        cout<<"testing calcGrowth():"<<endl;
        calcGrowth(0,rpt::echo);
        cout<<"testing calcMaturity():"<<endl;
        calcMaturity(0,rpt::echo);
        cout<<"testing calcSelectivities():"<<endl;
        calcSelectivities(dbgCalcProcs,rpt::echo);
        cout<<"testing calcFisheryFs():"<<endl;
        calcFisheryFs(dbgCalcProcs,rpt::echo);
        cout<<"testing calcSurveyQs():"<<endl;
        calcSurveyQs(dbgCalcProcs,cout);
        cout<<"testing runPopDyMod():"<<endl;
        runPopDyMod(0,cout);
        rpt::echo<<"n_yxm:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                   rpt::echo<<y<<cc;
                   rpt::echo<<tcsam::getSexType(x)<<cc;
                   rpt::echo<<tcsam::getMaturityType(m)<<cc;
                   rpt::echo<<sum(n_yxmsz(y,x,m))<<endl;
                }
            }
        }
        rpt::echo<<"n_yxmsz:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                       rpt::echo<<y<<cc;
                       rpt::echo<<tcsam::getSexType(x)<<cc;
                       rpt::echo<<tcsam::getMaturityType(m)<<cc;
                       rpt::echo<<tcsam::getShellType(s)<<cc;
                       rpt::echo<<n_yxmsz(y,x,m,s)<<endl;
                    }
                }
            }
        }
        cout<<"Testing calcObjFun()"<<endl;
        calcObjFun(1,rpt::echo);
        {cout<<"writing model results to R"<<endl;
            ofstream echo1; echo1.open("ModelRes0.R", ios::trunc);
            ReportToR(echo1,1,cout);
        }
        {cout<<"writing model sim data to file"<<endl;
            createSimData(0,cout);
            ofstream echo1; echo1.open("ModelSimData0.dat", ios::trunc);
            writeSimData(echo1,0,cout);
        }
        cout<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
        rpt::echo<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
    } else {
        writeMCMCHeader();
        cout<<"MCEVAL is on"<<endl;
    }
    
    
}

void model_parameters::userfunction(void)
{
  objFun =0.0;
    objFun.initialize();
    runPopDyMod(0,rpt::echo);
    calcObjFun(0,rpt::echo);
    
    if (sd_phase()){
        sdrLnR_y = log(R_y);
        for (int x=1;x<=nSXs;x++){
            for (int y=mnYr+ptrMDS->ptrBio->recLag; y<=mxYr; y++){
                sdrSpB_xy(x,y) = spb_yx(y,x);
            }
        }
    }
    
    if (mceval_phase()){
        updateMPI(0, cout);
        writeMCMCtoR(mcmc);
    }
}

void model_parameters::writeMCMCHeader(void)
{
    mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::trunc);
    mcmc<<"mcmc=list("<<endl;
    mcmc.close();
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,NumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedNumberVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,BoundedVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ostream& mcmc,DevsVectorVectorInfo* ptr)
{
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
}

void model_parameters::writeMCMCtoR(ofstream& mcmc)
{
    mcmc.open((char *) fnMCMC, ofstream::out|ofstream::app);
    mcmc<<"list(objFun="<<objFun<<cc<<endl;
    //write parameter values
        //recruitment values
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnR);   mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRCV); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLgtRX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRa);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRb);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pDevsLnR);  mcmc<<cc<<endl;
        //natural mortality parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMXM); mcmc<<cc<<endl;
        //growth parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrA); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrB); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGr->pLnGrBeta); mcmc<<cc<<endl;
        //maturity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrMat->pLgtPrMat); mcmc<<cc<<endl;
        //selectivity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS6); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS6); mcmc<<cc<<endl;
        //fully-selected fishing capture rate parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pHM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCXM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDevsLnC); mcmc<<cc<<endl;
        //survey catchability parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnQ);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQT);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQX);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQM);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQXM); mcmc<<cc<<endl;
    
        //write other quantities
        mcmc<<"R_y="; wts::writeToR(mcmc,value(R_y)); mcmc<<cc<<endl;
        mcmc<<"spb_xy="; wts::writeToR(mcmc,trans(value(spb_yx)),ptrMC->csvSXs,ptrMC->csvYrs); //mcmc<<cc<<endl;
        
    mcmc<<")"<<cc<<endl;
    mcmc.close();
    
}

void model_parameters::createSimData(int debug, ostream& cout)
{
    if (debug)cout<<"simulating model results as data"<<endl;
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vc_fyxmsz = wts::value(c_fyxmsz);
    d6_array vr_fyxmsz = wts::value(r_fyxmsz);
    for (int v=1;v<=nSrv;v++) {
        if (debug) cout<<"survey "<<v<<endl;
        (ptrSimMDS->ppSrv[v-1])->replaceCatchData(vn_vyxmsz(v),ptrSimMDS->ptrBio->wAtZ_xmz);
    }
    for (int f=1;f<=nFsh;f++) {
        if (debug) cout<<"fishery f: "<<f<<endl;
        (ptrSimMDS->ppFsh[f-1])->replaceCatchData(vc_fyxmsz(f),vr_fyxmsz(f),ptrSimMDS->ptrBio->wAtZ_xmz);
    }
    if (debug) cout<<"finished simulating model results as data"<<endl;
     
}

void model_parameters::writeSimData(ostream& os, int debug, ostream& cout)
{
    if (debug)cout<<"writing model results as data"<<endl;
    for (int v=1;v<=nSrv;v++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSimMDS->ppSrv[v-1]))<<endl;
    }
    //     cout<<4<<endl;
    for (int f=1;f<=nFsh;f++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSimMDS->ppFsh[f-1]))<<endl;
    }
    if (debug) cout<<"finished writing model results as data"<<endl;
     
}

void model_parameters::setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();
        for (int i=1;i<=np;i++) p(i)=vls(i);
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"p   = "<<p<<std::endl;
            std::cout<<"vls = "<<vls<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=vls.indexmin();j<=vls.indexmax();j++) p(i,j)=vls(j);
            if (debug>=dbgAll){
                std::cout<<"p(i)   = "<<p(i)<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
            }
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals "<<p(i).label()<<":"<<tb<<p(i)<<endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=vls.indexmin();j<=(vls.indexmax()-1);j++) p(i,j)=vls(j);
            if (debug>=dbgAll){
                std::cout<<"p(i)   = "<<p(i)<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
            }
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals "<<p(i).label()<<":"<<tb<<p(i)<<endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
}

void model_parameters::setAllDevs(int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"starting setAllDevs()"<<endl;
    tcsam::setDevs(devsLnR, pDevsLnR,debug,cout);
    tcsam::setDevs(devsS1, pDevsS1,debug,cout);
    tcsam::setDevs(devsS2, pDevsS2,debug,cout);
    tcsam::setDevs(devsS3, pDevsS3,debug,cout);
    tcsam::setDevs(devsS4, pDevsS4,debug,cout);
    tcsam::setDevs(devsS5, pDevsS5,debug,cout);
    tcsam::setDevs(devsS6, pDevsS6,debug,cout);
    
    tcsam::setDevs(devsLnC, pDevsLnC,debug,cout);
    if (debug>=dbgAll) cout<<"finished setAllDevs()"<<endl;
    
}

void model_parameters::runPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting runPopDyMod()"<<endl;
    //initialize population model
    initPopDyMod(debug, cout);
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,debug,cout);
        runPopDyModOneYear(y,debug,cout);        
    }
    doSurveys(mxYr+1,debug,cout);//do final surveys
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyMod()"<<endl;
    
}

void model_parameters::initPopDyMod(int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting initPopDyMod()"<<endl;
    
    spb_yx.initialize();
    n_yxmsz.initialize();
       
    setAllDevs(debug,cout);//set devs vectors
    
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcMaturity(debug,cout);   //calculate maturity ogives
    
    calcSelectivities(debug,cout);
    calcFisheryFs(debug,cout);
    calcSurveyQs(debug,cout);
    
    if (debug>=dbgPopDy) cout<<"finished initPopDyMod()"<<endl;
}

void model_parameters::doSurveys(int y,int debug,ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"starting doSurveys("<<y<<")"<<endl;
    for (int v=1;v<=nSrv;v++){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    n_vyxmsz(v,y,x,m,s) = elem_prod(q_vyxmsz(v,y,x,m,s),n_yxmsz(y,x,m,s));
                }
            }
        }
    }
    for (int v=1;v<=nSrv;v++){
        spb_vyx(v,y) = calcSpB(n_vyxmsz(v,y),y,debug,cout);
    }
    if (debug>=dbgPopDy) cout<<"finished doSurveys("<<y<<")"<<endl;
}

void model_parameters::runPopDyModOneYear(int yr, int debug, ostream& cout)
{
    if (debug>=dbgPopDy) cout<<"Starting runPopDyModOneYear("<<yr<<")"<<endl;
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n2_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n3_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n4_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n5_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    
    if (dtF(yr)<=dtM(yr)){//fishery occurs BEFORE molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs BEFORE molting/growth/maturity"<<endl;
        //apply natural mortality before fisheries
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtF(yr),debug,cout);
        //conduct fisheries
        n2_xmsz = applyFshMort(n1_xmsz,yr,debug,cout);
        //apply natural mortality after fisheries to molting/growth/maturity
        if (dtF(yr)==dtM(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtM(yr)-dtF(yr),debug,cout);
        }
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spb_yx(yr) = calcSpB(n3_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n4_xmsz = applyMGM(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtM(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtM(yr),debug,cout);
        }
    } else {              //fishery occurs AFTER molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs AFTER molting/growth/maturity"<<endl;
        //apply natural mortality before molting/growth/maturity
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtM(yr),debug,cout);
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spb_yx(yr) = calcSpB(n1_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n2_xmsz = applyMGM(n1_xmsz,yr,debug,cout);
        //apply natural mortality after molting/growth/maturity to fisheries
        if (dtM(yr)==dtF(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtF(yr)-dtM(yr),debug,cout);
        }
        //conduct fisheries
        n4_xmsz = applyFshMort(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtF(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtF(yr),debug,cout);
        }
    }
    
    //advance surviving individuals to next year
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n_yxmsz(yr+1,x,m,s) = n5_xmsz(x,m,s);
            }
        }
    }
    //add in recruits
    for (int x=1;x<=nSXs;x++) n_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) += R_y(yr)*R_yx(yr,x)*R_yz(yr);
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyModOneYear("<<yr<<")"<<endl;
    
}

dvar_vector model_parameters::calcSpB(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector spb(1,nSXs); spb.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int s=1;s<=nSCs;s++) spb(x) += n0_xmsz(x,MATURE,s)*ptrMDS->ptrBio->wAtZ_xmz(x,MATURE);//dot product here
    }
    if (debug>dbgApply) cout<<"finished calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return spb;
    
}

dvar4_array model_parameters::applyNatMort(dvar4_array& n0_xmsz, int y, double dt, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n1_xmsz(x,m,s) = elem_prod(exp(-M_yxmz(y,x,m)*dt),n0_xmsz(x,m,s));
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

dvar4_array model_parameters::applyFshMort(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector tm_z(1,nZBs);
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                Ft_yxmsz(y,x,m,s) = 0.0;//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) Ft_yxmsz(y,x,m,s) += Fr_fyxmsz(f,y,x,m,s)+Fd_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(mfexp(-Ft_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheriess
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);                                //numbers killed by all fisheries
                for (int f=1;f<=nFsh;f++){                   
                    c_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(Fc_fyxmsz(f,y,x,m,s),Ft_yxmsz(y,x,m,s)),tm_z);//numbers captured in fishery f
                    r_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(Fr_fyxmsz(f,y,x,m,s),Ft_yxmsz(y,x,m,s)),tm_z);//retained mortality in fishery f (numbers)
                    d_fyxmsz(f,y,x,m,s) = c_fyxmsz(f,y,x,m,s)-r_fyxmsz(f,y,x,m,s);//discarded catch (NOT mortality) in fishery f (numbers)                    
                }
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

dvar4_array model_parameters::applyMGM(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
{
    if (debug>dbgApply) cout<<"starting applyMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        n1_xmsz(x,IMMATURE,NEW_SHELL) = prGr_yxmzz(y,x,IMMATURE)*elem_prod(1.0-prMat_yxz(y,x),n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,IMMATURE,OLD_SHELL) = 0.0;
        n1_xmsz(x,MATURE,NEW_SHELL)   = prGr_yxmzz(y,x,MATURE)  *elem_prod(    prMat_yxz(y,x),n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,MATURE,OLD_SHELL)   = n0_xmsz(x,MATURE,NEW_SHELL)+n0_xmsz(x,MATURE,OLD_SHELL);
    }
    if (debug>dbgApply) cout<<"finished applyNatMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
}

void model_parameters::calcRecruitment(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcRecruitment()"<<endl;
    RecruitmentInfo* ptrRI = ptrMPI->ptrRec;
    
    R_y.initialize();
    Rx_c.initialize();
    R_yx.initialize();
    R_cz.initialize();
    R_yz.initialize();
    zscrDevsLnR.initialize();
    
    int k; int y;
    dvector dzs = zBs+(zBs[2]-zBs[1])/2.0-zBs[1];
    for (int pc=1;pc<=ptrRI->nPCs;pc++){
        ivector pids = ptrRI->getPCIDs(pc);
        k=ptrRI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvariable mnLnR    = pLnR(pids[k++]);
        dvariable lnRCV    = pLnRCV(pids[k++]);
        dvariable lgtRX    = pLgtRX(pids[k++]);
        dvariable lnRa     = pLnRa(pids[k++]);
        dvariable lnRb     = pLnRb(pids[k++]);
        if (debug>dbgCalcProcs){
            cout<<"pids  = "<<pids<<endl;
            cout<<"mnLnR = "<<mnLnR<<endl;
            cout<<"lnRCV = "<<lnRCV<<endl;
            cout<<"lgtRX = "<<lgtRX<<endl;
            cout<<"lnRa  = "<<lnRa<<endl;
            cout<<"lnRb  = "<<lnRb<<endl;
        }
        int useDevs = pids[k]; k++;
        dvariable mdR;
        dvar_vector dvsLnR;
        ivector idxDevsLnR;
        if (useDevs) {
            dvsLnR     = devsLnR(useDevs);
            idxDevsLnR = idxsDevsLnR(useDevs);
            if (debug>dbgCalcProcs) {
                cout<<"lims(dvsLnR) = "<<dvsLnR.indexmin()<<cc<<dvsLnR.indexmax()<<endl;
                cout<<"idx(dvsLnR) = "<<idxDevsLnR<<endl;
                cout<<"dvsLnR = "<<dvsLnR<<endl;
            }
        } else {
            mdR = mfexp(mnLnR);
        }
        
        Rx_c(pc) = 1.0/(1.0+mfexp(-lgtRX));
        R_cz(pc) = elem_prod(pow(dzs,mfexp(lnRa-lnRb)-1.0),mfexp(-dzs/mfexp(lnRb)));
        R_cz(pc) /= sum(R_cz(pc));//normalize to sum to 1
        imatrix idxs = ptrRI->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if (debug>dbgCalcProcs) cout<<"y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            if (useDevs){
                R_y(y) = mfexp(mnLnR+dvsLnR[idxDevsLnR[y]]);
            } else {
                R_y(y) = mdR;
            }
            if (debug>dbgCalcProcs) cout<<R_y(y)<<tb;
            R_yx(y,MALE)   = Rx_c(pc);
            if (debug>dbgCalcProcs) cout<<R_yx(y,MALE)<<endl;
            R_yx(y,FEMALE) = 1.0-R_yx(y,MALE);
            
            R_yz(y) = R_cz(pc);
            if (debug>dbgCalcProcs) cout<<R_yz(y)<<endl;
            
            zscrDevsLnR(y) = dvsLnR(idxDevsLnR(y))/sqrt(log(1.0+mfexp(2.0*lnRCV)));//standardized ln-scale rec devs
        }
    }
    
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz = "<<endl<<R_yz<<endl;
        cout<<"zscr = "<<zscrDevsLnR<<endl;
        cout<<"finished calcRecruitment()"<<endl;
    }
}

void model_parameters::calcNatMort(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcNatMort()"<<endl;
    
    NaturalMortalityInfo* ptrNM = ptrMPI->ptrNM;
    
    dvar_matrix lnM(1,nSXs,1,nMSs);
    dvar3_array M_xmz(1,nSXs,1,nMSs,1,nZBs);
    
    M_cxm.initialize();
    M_yxmz.initialize();
    int y; 
    for (int pc=1;pc<=ptrNM->nPCs;pc++){
        lnM.initialize();
        ivector pids = ptrNM->getPCIDs(pc);
        int k=ptrNM->nIVs+1;//1st parameter variable column
        //add in base (ln-scale) natural mortality (mature males)
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x) += pLnM(pids[k]);}   k++;
        //add in main temporal offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x) += pLnDMT(pids[k]);} k++;
        //add in female offset
        if (pids[k]) {lnM(FEMALE) += pLnDMX(pids[k]);}                      k++;
        //add in immature offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x,IMMATURE) += pLnDMM(pids[k]);} k++;
        //add in offset immature females for stanza
        if (pids[k]) {lnM(FEMALE,IMMATURE) += pLnDMXM(pids[k]);}            k++; //advance k to zScaling in pids
        
        //convert from ln-scale to arithmetic scale
        M_cxm(pc) = mfexp(lnM);
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"lnM:"<<endl<<lnM<<endl;
            cout<<"pc: "<<pc<<tb<<"M_xm:"<<endl<<M_cxm(pc)<<endl;
        }
        
        //add in size-scaling, if requested
        M_xmz.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                if (pids[k]&&(current_phase()>=pids[k])) {
                    M_xmz(x,m) = M_cxm(pc,x,m)*(zMref/zBs);//factor in size dependence
                } else {
                    M_xmz(x,m) = M_cxm(pc,x,m);//no size dependence
                }
            }
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrNM->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//only model index for natural mortality is year
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    M_yxmz(y,x,m) = M_xmz(x,m);
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"Finished calcNatMort()"<<endl;
    
}

void model_parameters::calcMaturity(int debug, ostream& cout)
{
    if (debug>dbgCalcProcs) cout<<"starting calcMaturity()"<<endl;
    MaturityInfo* ptrMI = ptrMPI->ptrMat;
    
    prMat_cz.initialize();
    prMat_yxz.initialize();
    
    int k; int y; int x;
    for (int pc=1;pc<=ptrMI->nPCs;pc++){
        ivector pids = ptrMI->getPCIDs(pc);
        k=ptrMI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvar_vector lgtPrMat = pLgtPrMat(pids[k++]);
        int vmn = lgtPrMat.indexmin();
        int vmx = lgtPrMat.indexmax();
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<vmn<<", mx = "<<vmx<<endl;
            cout<<"lgtPrMat = "<<lgtPrMat<<endl;
        }
        prMat_cz(pc) = 1.0;//default is 1
        prMat_cz(pc)(vmn,vmx) = 1.0/(1.0+mfexp(-lgtPrMat));
            
        imatrix idxs = ptrMI->getModelIndices(pc);
        if (debug) cout<<"maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            x = idxs(idx,2);
            if (debug) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(x)<<endl;
            prMat_yxz(y,x) = prMat_cz(pc);//note: this change made a difference, but not sure why!
        }
    }
    
    if (debug>dbgCalcProcs) cout<<"finished calcMaturity()"<<endl;
}

void model_parameters::calcGrowth(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcGrowth()"<<endl;
    
    GrowthInfo* ptrGrI = ptrMPI->ptrGr;
    
    dvariable grA;
    dvariable grB;
    dvariable grBeta;
    
    prGr_czz.initialize();
    prGr_yxmzz.initialize();
    
    dvar_matrix prGr_zz(1,nZBs,1,nZBs);
    int y; int x;
    for (int pc=1;pc<=ptrGrI->nPCs;pc++){
        ivector pids = ptrGrI->getPCIDs(pc);
        int k=ptrGrI->nIVs+1;//1st parameter column
        grA = mfexp(pLnGrA(pids[k])); k++; //"a" coefficient for mean growth
        grB = mfexp(pLnGrB(pids[k])); k++; //"b" coefficient for mean growth
        grBeta = mfexp(pLnGrBeta(pids[k])); k++; //shape factor for gamma function growth transition
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<". grB:"<<tb<<grB<<". grBeta:"<<grBeta<<endl;
        }
        
        //compute growth transition matrix for this pc
        prGr_zz.initialize();
        dvar_vector mnZ = mfexp(grA)*pow(zBs,grB);//mean size after growth from zBs
        dvar_vector alZ = (mnZ-zBs)/grBeta;//scaled mean growth increment from zBs
        for (int z=1;z<nZBs;z++){//pre-molt growth bin
            dvar_vector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
            if (debug) cout<<"dZs: "<<dZs.indexmin()<<":"<<dZs.indexmax()<<endl;
            dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),mfexp(-dZs/grBeta)); //pr(dZ|z)
            if (debug) cout<<"prs: "<<prs.indexmin()<<":"<<prs.indexmax()<<endl;
            if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
            if (debug) cout<<prs<<endl;
            prs = prs/sum(prs);//normalize to sum to 1
            if (debug) cout<<prs<<endl;
            prGr_zz(z)(z,nZBs) = prs;
        }
        prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., "to") z's so n+ = prGr_zz*n
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrGrI->getModelIndices(pc);
        if (debug) cout<<"growth indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1); //year index
            x = idxs(idx,2); //sex index
            for (int m=1;m<=nMSs;m++){
                for (int z=1;z<=nZBs;z++){
                    prGr_yxmzz(y,x,m,z) = prGr_czz(pc,z);
                }
            }
        }
    }
    
    if (debug>dbgCalcProcs) cout<<"finished calcGrowth()"<<endl;
}

void model_parameters::calcSelectivities(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSelectivities()"<<endl;
    
    SelectivityInfo* ptrSel = ptrMPI->ptrSel;
    double fsZ;             //fully selected size
    int idSel;              //selectivity function id
    int idxFSZ = ptrSel->nIVs+ptrSel->nPVs+1;//index for fsZ in pids vector below
    
    ivector mniSelDevs(1,6);//min indices of devs vectors
    ivector mxiSelDevs(1,6);//max indices of devs vectors
    dvar_vector params(1,6);//vector for number_vector params
        
    sel_cz.initialize();//selectivities w/out deviations
    sel_iyz.initialize();//selectivity array
    int y;
    for (int pc=1;pc<=ptrSel->nPCs;pc++){
        params.initialize();
        ivector pids = ptrSel->getPCIDs(pc);
        //extract the number parameters
        int k=ptrSel->nIVs+1;//1st parameter variable column
        if (pids[k]) {params[1] = pS1(pids[k]);}   k++;
        if (pids[k]) {params[2] = pS2(pids[k]);}   k++;
        if (pids[k]) {params[3] = pS3(pids[k]);}   k++;
        if (pids[k]) {params[4] = pS4(pids[k]);}   k++;
        if (pids[k]) {params[5] = pS5(pids[k]);}   k++;
        if (pids[k]) {params[6] = pS6(pids[k]);}   k++;
        if (debug>dbgCalcProcs) {
            cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
            cout<<tb<<"params:"<<tb<<params<<endl;
        }
        
        int useDevsS1=pids[k++];
        dvar_vector dvsS1; ivector idxDevsS1;
        if (useDevsS1){
            dvsS1 = devsS1(useDevsS1);
            idxDevsS1 = idxsDevsS1(useDevsS1);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS1) = "<<idxDevsS1<<endl;
                cout<<"dvsS1      = "<<dvsS1<<endl;
            }
        }
        int useDevsS2=pids[k++];
        dvar_vector dvsS2; ivector idxDevsS2;
        if (useDevsS2){
            dvsS2 = devsS2(useDevsS2);
            idxDevsS2 = idxsDevsS2(useDevsS2);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS2) = "<<idxDevsS2<<endl;
                cout<<"dvsS2      = "<<dvsS2<<endl;
            }
        }
        int useDevsS3=pids[k++];
        dvar_vector dvsS3; ivector idxDevsS3;
        if (useDevsS3){
            dvsS3 = devsS3(useDevsS3);
            idxDevsS3 = idxsDevsS3(useDevsS3);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS3) = "<<idxDevsS3<<endl;
                cout<<"dvsS3      = "<<dvsS3<<endl;
            }
        }
        int useDevsS4=pids[k++];
        dvar_vector dvsS4; ivector idxDevsS4;
        if (useDevsS4){
            dvsS4 = devsS4(useDevsS4);
            idxDevsS4 = idxsDevsS4(useDevsS4);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS4) = "<<idxDevsS4<<endl;
                cout<<"dvsS4      = "<<dvsS4<<endl;
            }
        }
        int useDevsS5=pids[k++];
        dvar_vector dvsS5; ivector idxDevsS5;
        if (useDevsS5){
            dvsS5 = devsS5(useDevsS5);
            idxDevsS5 = idxsDevsS5(useDevsS5);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS5) = "<<idxDevsS5<<endl;
                cout<<"dvsS5      = "<<dvsS5<<endl;
            }
        }
        int useDevsS6=pids[k++];
        dvar_vector dvsS6; ivector idxDevsS6;
        if (useDevsS6){
            dvsS6 = devsS6(useDevsS6);
            idxDevsS6 = idxsDevsS6(useDevsS6);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS6) = "<<idxDevsS6<<endl;
                cout<<"dvsS6      = "<<dvsS6<<endl;
            }
        }
        fsZ   = pids[idxFSZ];
        idSel = pids[idxFSZ+1];
        if (debug>dbgCalcProcs) cout<<tb<<"fsZ: "<<fsZ<<tb<<"idSel"<<tb<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<endl;;
        sel_cz(pc) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
            
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSel->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//year
            k=ptrSel->nIVs+1+6;//1st devs vector variable column
            if (useDevsS1) {params[1] += devsS1(useDevsS1,idxDevsS1[y]);}
            if (useDevsS2) {params[2] += devsS2(useDevsS2,idxDevsS2[y]);}
            if (useDevsS3) {params[3] += devsS3(useDevsS3,idxDevsS3[y]);}
            if (useDevsS4) {params[4] += devsS4(useDevsS4,idxDevsS4[y]);}
            if (useDevsS5) {params[5] += devsS5(useDevsS5,idxDevsS5[y]);}
            if (useDevsS6) {params[6] += devsS6(useDevsS6,idxDevsS6[y]);}
            sel_iyz(pc,y) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
            if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"sel: "<<sel_iyz(pc,y)<<endl;
        }
    }
    if (debug>dbgCalcProcs) cout<<"finished calcSelectivities()"<<endl;
}

void model_parameters::calcFisheryFs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcFisheryFs()"<<endl;
    
    FisheriesInfo* ptrFsh = ptrMPI->ptrFsh;
    
    dvariable hm;                   //handling mortality
    dvar_matrix lnC(1,nSXs,1,nMSs); //ln-scale capture rate
    dvar_matrix C_xm(1,nSXs,1,nMSs);//arithmetic-scale capture rate
    
    Fhm_fy.initialize();   //handling mortality
    Fc_fyxms.initialize(); //fully-selected capture rate
    Fc_fyxmsz.initialize();//size-specific capture rate
    Fr_fyxmsz.initialize();//retention rate
    Fd_fyxmsz.initialize();//discard mortality rate
    Ft_yxmsz.initialize(); //total mortality rate
    
    /******************************************************\n
     * Fully-selected annual capture rates are calculated  \n
     * using 2 approaches:                                 \n
     *  1. from parameter values (if useER=0 below)        \n
     *  2. based on effort and the average ratio between   \n
     *     effort and capture rates over some time period  \n
     *     in the fishery        (if useER=1 below)        \n
     * Consequently, calculating all the above quantities  \n
     * requires 2 passes through parameter combinations.   \n
    ******************************************************/
    int idxER = ptrFsh->idxUseER;//index into pids below for flag to use effort ratio
    int y; int f; int x; int idSel; int idRet; int useER; int useDevs;
    //Pass 1: calculations based on parameter values
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids: "<<pids<<endl;
        useER = pids[idxER];//flag to use effort ratio
        if (!useER){//calculate capture rates from parameters
            lnC.initialize();
            C_xm.initialize();
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);}                                   k++;
            //set base (ln-scale) capture rate (mature males)
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x) += pLnC(pids[k]);}   k++;
            //add in main temporal offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x) += pLnDCT(pids[k]);} k++;
            //add in female offset
            if (pids[k]) {lnC(FEMALE) += pLnDCX(pids[k]);}                      k++;
            //add in immature offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x,IMMATURE) += pLnDCM(pids[k]);} k++;
            //add in offset immature females for stanza
            if (pids[k]) {lnC(FEMALE,IMMATURE) += pLnDCXM(pids[k]);}            k++; 
            //extract devs vector
            useDevs = pids[k]; k++;
            dvar_vector dvsLnC;             
            ivector idxDevsLnC;
            if (useDevs) {
                dvsLnC     = devsLnC(useDevs);
                idxDevsLnC = idxsDevsLnC(useDevs);
                if (debug>dbgCalcProcs){
                    cout<<"y   idx    devsLnC"<<endl;
                    for (int i=idxDevsLnC.indexmin();i<=idxDevsLnC.indexmax();i++) {
                        cout<<i<<tb<<idxDevsLnC(i)<<tb;
                        if (idxDevsLnC(i)) cout<<dvsLnC[idxDevsLnC(i)];
                        cout<<endl;
                    }
                }
            } else {
                C_xm = mfexp(lnC);
            }
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];//selectivity function id
            idRet = pids[k++];//retention function id
        
            //convert from ln-scale to arithmetic scale
            if (debug>dbgCalcProcs){
                cout<<"pc: "<<pc<<". idSel = "<<idSel<<". idRet = "<<idRet<<"."<<endl;
                cout<<tb<<tb<<"lnC:"<<endl<<lnC<<endl;
                if (useDevs) {
                    cout<<tb<<tb<<"dvsLnC["<<dvsLnC.indexmin()<<cc<<dvsLnC.indexmax()<<"] = "<<dvsLnC<<endl;
                } else {
                    cout<<tb<<tb<<"C_xm:"<<endl<<C_xm<<endl;
                }
            }
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                x = idxs(idx,3);//sex
                if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                if (useDevs) C_xm = mfexp(lnC+dvsLnC[idxDevsLnC[y]]);//recalculate C_xm w/ devs
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        Fc_fyxms(f,y,x,m,s)  = C_xm(x,m);                 //fully-selected capture rate
                        Fc_fyxmsz(f,y,x,m,s) = C_xm(x,m)*sel_iyz(idSel,y);//size-specific capture rate
                        if (idRet){//fishery has retention
                            Fr_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         Fc_fyxmsz(f,y,x,m,s));//retention mortality
                            Fd_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),Fc_fyxmsz(f,y,x,m,s));//discard mortality
                        } else {//discard only
                            Fd_fyxmsz(f,y,x,m,s) = hm*Fc_fyxmsz(f,y,x,m,s);//discard mortality
                        }
                    }
                }
            }
        }//useER=FALSE
    }
    if (debug) cout<<"finished pass 1"<<endl;
    
    //calculate ratio of average capture rate to effort
    if (debug>dbgCalcProcs) cout<<"calculating avgRatioFc2Eff"<<endl;
    dvariable tot;
    for (int f=1;f<=nFsh;f++){//model fishery objects
        int fd = mapM2DFsh(f);//index of corresponding fishery data object
        if (ptrMDS->ppFsh[fd-1]->ptrEff){
            IndexRange* pir = ptrMDS->ppFsh[fd-1]->ptrEff->ptrAvgIR;
            int mny = max(mnYr,pir->getMin());//adjust for model min
            int mxy = min(mxYr,pir->getMax());//adjust for model max
            if (debug>dbgCalcProcs) cout<<"f,mny,mxy = "<<f<<tb<<mny<<tb<<mxy<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        tot.initialize();
                        switch (optsFcAvg(f)){
                            case 1:
                                for (int y=mny;y<=mxy;y++) tot += Fc_fyxms(f,y,x,m,s);
                                avgFc(f,x,m,s) = tot/(mxy-mny+1); break;
                            case 2:
                                for (int y=mny;y<=mxy;y++) tot += 1.0-mfexp(-Fc_fyxms(f,y,x,m,s));
                                avgFc(f,x,m,s) = tot/(mxy-mny+1); break;
                            case 3:
                                for (int y=mny;y<=mxy;y++) tot += mean(Fc_fyxmsz(f,y,x,m,s));
                                avgFc(f,x,m,s) = tot/(mxy-mny+1); break;
                            default:
                                cout<<"optsFcAvg("<<f<<") = "<<optsFcAvg(f)<<" to calculate average Fc is invalid."<<endl;
                                cout<<"Aborting..."<<endl;
                                exit(-1);
                        }
                        avgRatioFc2Eff(f,x,m,s) = avgFc(f,x,m,s)/avgEff(f);
                    }
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"calculated avgRatioFc2Eff"<<endl;
    
    //Pass 2: calculations based on effort and average effort:capture rate ratios
    int fd; double eff;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        useER = pids[idxER];//flag to use effort ratio
        if (useER){//calculate capture rates from parameters
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);} k++;
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];   //selectivity function id
            idRet = pids[k++];   //retention function id
            
            if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<". hm = "<<hm<<". idSel = "<<idSel<<". idRet = "<<idRet<<". Using ER"<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                x = idxs(idx,3);//sex
                fd = mapM2DFsh(f);//index of corresponding fishery data object
                eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                if (debug>dbgCalcProcs) cout<<"f,y,x,eff = "<<f<<tb<<y<<tb<<x<<tb<<eff<<endl;
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        //fully-selected capture rate
                        switch(optsFcAvg(f)) {
                            case 1:
                                Fc_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                            case 2:
                                Fc_fyxms(f,y,x,m,s) = -log(1.0-avgRatioFc2Eff(f,x,m,s)*eff); break;
                            case 3:
                                Fc_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                        }
                        Fc_fyxmsz(f,y,x,m,s) = Fc_fyxms(f,y,x,m,s)*sel_iyz(idSel,y);//size-specific capture rate
                        if (idRet){//fishery has retention
                            Fr_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         Fc_fyxmsz(f,y,x,m,s));//retention mortality
                            Fd_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),Fc_fyxmsz(f,y,x,m,s));//discard mortality
                        } else {//discard only
                            Fd_fyxmsz(f,y,x,m,s) = hm*Fc_fyxmsz(f,y,x,m,s);//discard mortality
                        }
                    }
                }
            }
        }//useER=TRUE
    }
    if (debug>dbgCalcProcs) cout<<"finished pass 2."<<endl;
    if (debug>dbgCalcProcs) cout<<"finished calcFisheryFs()"<<endl;
}

void model_parameters::calcSurveyQs(int debug, ostream& cout)
{
    if(debug>dbgCalcProcs) cout<<"Starting calcSurveyQs()"<<endl;
    
    SurveysInfo* ptrSrv = ptrMPI->ptrSrv;
    
    dvar_matrix lnQ(1,nSXs,1,nMSs);
    dvar_matrix Q_xm(1,nSXs,1,nMSs);
    
    q_vyxms.initialize();
    q_vyxmsz.initialize();
    int y; int v; int x; int idSel;
    for (int pc=1;pc<=ptrSrv->nPCs;pc++){
        lnQ.initialize();
        Q_xm.initialize();
        ivector pids = ptrSrv->getPCIDs(pc);
        int k=ptrSrv->nIVs+1;//1st parameter variable column
        //add in base (ln-scale) catchability (mature males)
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x) += pLnQ(pids[k]);}   k++;
        //add in main temporal offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x) += pLnDQT(pids[k]);} k++;
        //add in female offset
        if (pids[k]) {lnQ(FEMALE) += pLnDQX(pids[k]);}                      k++;
        //add in immature offsets
        if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x,IMMATURE) += pLnDQM(pids[k]);} k++;
        //add in offset immature females for stanza
        if (pids[k]) {lnQ(FEMALE,IMMATURE) += pLnDQXM(pids[k]);}            k++; 
        
        idSel = pids[k];//selectivity function id
        
        //convert from ln-scale to arithmetic scale
        Q_xm = mfexp(lnQ);
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"lnQ:"<<endl<<lnQ<<endl;
            cout<<"pc: "<<pc<<tb<<"Q_xm:"<<endl<<Q_xm<<endl;
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSrv->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            v = idxs(idx,1);//survey
            y = idxs(idx,2);//year
            x = idxs(idx,3);//sex
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    q_vyxms(v,y,x,m,s)  = Q_xm(x,m);
                    q_vyxmsz(v,y,x,m,s) = Q_xm(x,m)*sel_iyz(idSel,y);
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
    
}

void model_parameters::calcPenalties(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    //smoothness penalties on maturity parameters (NOT maturity ogives)
    double penWgtLgtPrMat = 1.0;//TODO: read in value from input file
    fPenSmoothLgtPrMat.initialize();
    cout<<"fPenSmooth = ";
    for (int i=1;i<=npLgtPrMat;i++){
        dvar_vector v; v = 1.0*pLgtPrMat(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        cout<<fPenSmoothLgtPrMat(i)<<tb;
        objFun += penWgtLgtPrMat*fPenSmoothLgtPrMat(i);
    }
    cout<<endl;
    //non-decreasing penalties on maturity parameters (NOT maturity ogives)
    double penWgtNonDecLgtPrMat = 1.0;//TODO: read in value from input file
    fPenNonDecLgtPrMat.initialize();
    cout<<"fPenNonDec = ";
    for (int i=1;i<=npLgtPrMat;i++){
        dvar_vector v; v = calc1stDiffs(pLgtPrMat(i));
        for (int iv=v.indexmin();iv<=v.indexmax();iv++){
            posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
        }
        cout<<fPenNonDecLgtPrMat(i)<<tb;
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
    }
    cout<<endl;
    
    if (debug>=dbgObjFun) cout<<"Finished calcPenalties()"<<endl;
}

dvar_vector model_parameters::calc1stDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector cp; cp = d;
    dvar_vector r(mn,mx-1);
    r = cp(mn+1,mx).shift(mn)-cp(mn,mx-1);
    RETURN_ARRAYS_DECREMENT();
    return r;
}

dvar_vector model_parameters::calc2ndDiffs(const dvar_vector& d)
{
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector r = calc1stDiffs(calc1stDiffs(d));
    RETURN_ARRAYS_DECREMENT();
    return r;
}

void model_parameters::calcNLLCompRec(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Starting calcNLLCompRec"<<endl;
    //recruitment devs
    double nllWgtRecDevs = 1.0;//TODO: read in from input file (as vector?))
    nllRecDevs.initialize();
    nllRecDevs = 0.5*norm2(zscrDevsLnR);//
    objFun += nllWgtRecDevs*nllRecDevs;
    if (debug<0){
        cout<<"list(type='normal',nll="<<nllRecDevs<<cc; wts::writeToR(cout,value(zscrDevsLnR)); cout<<")";
    }
    if (debug>=dbgObjFun) cout<<"Finished calcNLLCompRec"<<endl;
}

void model_parameters::calcObjFun(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Starting calcObjFun"<<endl;
    //objective function penalties
    calcPenalties(debug,cout);
    //prior likelihoods
    calcAllPriors(debug,cout);
    //recruitment component
    calcNLLCompRec(debug,cout);
    
    //data components
    calcNLLs_Fisheries(debug,cout);
    calcNLLs_Surveys(debug,cout);
    
    if (debug>=dbgObjFun) cout<<"Finished calcObjFun"<<endl;
    
}

dvariable model_parameters::calcNorm2NLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcNorm2NLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y]);
        }
    }
    nll += 0.5*norm2(zscr);
    if (debug<0){
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(type='norm2',nll="<<nll<<cc; 
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),yrs); cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNorm2NLL()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return nll;
    
}

dvariable model_parameters::calcNormalNLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
   if (debug>=dbgAll) cout<<"Starting calcNormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y])/stdv[i];
            nll += log(stdv[i]);
        }
    }
    nll += 0.5*norm2(zscr);
    if (debug<0){
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(type='normal',nll="<<nll<<cc;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),yrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcNormalNLL()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return nll;
    
}

dvariable model_parameters::calcLognormalNLL(dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcLognormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (log(obs[i]+smlVal)-log(mod[y]+smlVal))/stdv[i];
            nll += log(stdv[i]);
        }
    }
    nll += 0.5*norm2(zscr);
    if (debug<0){
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(type='lognormal',nll="<<nll<<cc; 
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),yrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcLognormalNLL()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return nll;
    
}

dvar_vector model_parameters::calcNLLs_CatchAbundance(AggregateCatchData* ptrA, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcNLLs_CatchAbundance()"<<endl;
    dvar_vector nlls;
    dvar_matrix tA_xy;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    if (debug<0) cout<<"list(type='"<<tcsam::getFitType(ptrA->optFit)<<"',fits=list("<<endl;
    if (ptrA->optFit==tcsam::FIT_BY_SEX){
        nlls.allocate(1,nSXs); nlls.initialize();
        tA_xy.allocate(1,nSXs,mny,mxy); tA_xy.initialize();
        ivector yrs = ptrA->yrs;
        for (int x=1;x<=nSXs;x++){
            for (int y=mny;y<=mxy;y++) tA_xy(x,y) = sum(mA_yxmsz(y,x));//sum by sex
            if (debug<0) cout<<tcsam::getSexType(x)<<"=";
            nlls(x) = calcNLL(ptrA->llType, tA_xy(x), ptrA->C_xy(x), ptrA->stdv_xy(x), ptrA->yrs, debug, cout); 
            if (debug<0) cout<<","<<endl;
        }
        if (debug<0) cout<<"NULL)";
    } else 
    if (ptrA->optFit==tcsam::FIT_BY_TOTAL){
        nlls.allocate(ANY_SX,ANY_SX); nlls.initialize();
        tA_xy.allocate(ANY_SX,ANY_SX,mny,mxy);
        tA_xy.initialize();
        for (int y=mny;y<=mxy;y++) tA_xy(ANY_SX,y) = sum(mA_yxmsz(y));//sum over sexes
        if (debug<0) cout<<tcsam::getSexType(ANY_SX)<<"=";
        nlls(ANY_SX) = calcNLL(ptrA->llType, tA_xy(ANY_SX), ptrA->C_xy(ANY_SX), ptrA->stdv_xy(ANY_SX), ptrA->yrs, debug, cout);                
        if (debug<0) cout<<")";
    } else {
        std::cout<<"Calling calcNLLs_CatchAbundance with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrA->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<")";
    if (debug>0) cout<<"nlls = "<<nlls<<endl;
    ptrA->saveNLLs(nlls);
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_CatchAbundance()"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return nlls;
}

dvar_vector model_parameters::calcNLLs_CatchBiomass(AggregateCatchData* ptrB, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcNLLs_CatchBiomass()"<<endl;
    dvar_vector nlls;
    dvar_matrix tB_xy;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    if (debug<0) cout<<"list(type='"<<tcsam::getFitType(ptrB->optFit)<<"',fits=list("<<endl;
    if (ptrB->optFit==tcsam::FIT_BY_SEX){
        nlls.allocate(1,nSXs); nlls.initialize();
        tB_xy.allocate(1,nSXs,mny,mxy); tB_xy.initialize();
        ivector yrs = ptrB->yrs;
        //calc catch biomass by sex
        for (int x=1;x<=nSXs;x++){
            for (int y=mny;y<=mxy;y++) {
                for (int m=1;m<=nMSs;m++) {
                    //sum by sex over maturity states, shell conditions, sizes
                    for (int s=1;s<=nSCs;s++) tB_xy(x,y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                }
            }
            if (debug<0) cout<<tcsam::getSexType(x)<<"=";
            nlls(x) = calcNLL(ptrB->llType, tB_xy(x), ptrB->C_xy(x), ptrB->stdv_xy(x), ptrB->yrs, debug, cout); 
            if (debug<0) cout<<","<<endl;
        }
        if (debug<0) cout<<"NULL)";
    } else 
    if (ptrB->optFit==tcsam::FIT_BY_TOTAL){
        nlls.allocate(ANY_SX,ANY_SX); nlls.initialize();
        tB_xy.allocate(ANY_SX,ANY_SX,mny,mxy);
        tB_xy.initialize();
        //calc catch biomass over sexes
        for (int x=1;x<=nSXs;x++){
            for (int y=mny;y<=mxy;y++) {
                for (int m=1;m<=nMSs;m++) {
                    for (int s=1;s<=nSCs;s++) tB_xy(ANY_SX,y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                }
            }
        }
        if (debug<0) cout<<tcsam::getSexType(ANY_SX)<<"=";
        nlls(ANY_SX) = calcNLL(ptrB->llType, tB_xy(ANY_SX), ptrB->C_xy(ANY_SX), ptrB->stdv_xy(ANY_SX), ptrB->yrs, debug, cout);                
        if (debug<0) cout<<")";
    } else {
        std::cout<<"Calling calcNLLs_CatchBiomass with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrB->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<")";
    if (debug>0) cout<<"nlls = "<<nlls<<endl;
    ptrB->saveNLLs(nlls);
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_CatchBiomass()"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return nlls;
}

dvariable model_parameters::calcMultinomialNLL(dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcMultinomialNLL()"<<endl;
    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
    if (debug<0){
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(type='multinomial',nll="<<nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        cout<<"nlls=";  wts::writeToR(cout,nlls, ptrMC->csvZBs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,ptrMC->csvZBs); cout<<endl;
        cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcMultinomialNLL()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return nll;
    
}

dvariable model_parameters::calcNLL(int llType, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    dvariable nll; nll.initialize();
    switch (llType){
        case tcsam::LL_NONE:
            break;
        case tcsam::LL_LOGNORMAL:
            nll = calcLognormalNLL(mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORMAL:
            nll = calcNormalNLL(mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORM2:
            nll = calcNorm2NLL(mod,obs,stdv,yrs,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(1)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }    
    RETURN_ARRAYS_DECREMENT();
    return nll;
}

dvariable model_parameters::calcNLL(int llType, dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    dvariable nll; nll.initialize();
    switch (llType){
        case tcsam::LL_NONE:
            break;
        case tcsam::LL_MULTINOMIAL:
            nll = calcMultinomialNLL(mod,obs,ss,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(2)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }
    
    RETURN_ARRAYS_DECREMENT();
    return nll;
}

dvar_matrix model_parameters::calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgAll) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    ivector yrs = ptrZFD->yrs;
    dvar_matrix nlls;
    int y;
    double ss;
    dvariable nT;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    if (ptrZFD->optFit==tcsam::FIT_BY_TOTAL){
        nlls.allocate(ANY_SX,ANY_SX,1,yrs.size()); nlls.initialize();
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_SEX){
        nlls.allocate(1,nSXs,1,yrs.size()); nlls.initialize();
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_SEX_EXTENDED){
        nlls.allocate(1,nSXs,1,yrs.size()); nlls.initialize();
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_SEX_MAT_EXTENDED){
        nlls.allocate(1,nSXs*nMSs,1,yrs.size()); nlls.initialize();
    } else
    {
        std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<"list("<<endl;
    for (int iy=1;iy<=yrs.size();iy++) {
        y = yrs[iy];
        if (debug>0) cout<<"y = "<<y<<endl;
        if ((mny<=y)&&(y<=mxy)) {
            if (ptrZFD->optFit==tcsam::FIT_BY_TOTAL){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    dvector     oP_z(1,nZBs); oP_z.initialize();//observed size comp.
                    dvar_vector mP_z(1,nZBs); mP_z.initialize();//model size comp.
                    for (int x=1;x<=(nSXs+1);x++){
                        for (int m=1;m<=(nMSs+1);m++) {
                            for (int s=1;s<=(nSCs+1);s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                    }
                    oP_z /= sum(oP_z);
                    if (debug>0){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++)mP_z += mA_yxmsz(y,x,m,s);
                        }
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                    if (debug<0) {
                        cout<<"'"<<y<<"'=list(";
                        cout<<"fit.type="<<ptrZFD->optFit<<cc;
                        cout<<"sx='ANY_SEX'"<<cc;
                        cout<<"ms='ANY_MATURITY'"<<cc;
                        cout<<"sc='ANY_SHELL_CONDITION'"<<cc;
                        cout<<"fit=";
                    }
                    nlls(ANY_SX,iy) = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                    if (debug<0) cout<<")"<<cc<<endl;
                }
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_SEX){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        dvector     oP_z(1,nZBs); oP_z.initialize();//observed size comp.
                        dvar_vector mP_z(1,nZBs); mP_z.initialize();//model size comp.
                        for (int m=1;m<=(nMSs+1);m++) {
                            for (int s=1;s<=(nSCs+1);s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        oP_z /= sum(oP_z);
                        if (debug>0){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type="<<ptrZFD->optFit<<cc;
                            cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"ms='ANY_MATURITY'"<<cc;
                            cout<<"sc='ANY_SHELL_CONDITION'"<<cc;
                            cout<<"fit=";
                        }
                        nlls(x,iy) = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_SEX_EXTENDED){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    dvector     oP_z(1,nSXs*nZBs); oP_z.initialize();//observed size comp.
                    dvar_vector mP_z(1,nSXs*nZBs); mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        for (int m=1;m<=(nMSs+1);m++) {
                            for (int s=1;s<=(nSCs+1);s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                        }
                    }//x
                    oP_z /= sum(oP_z);
                    if (debug>0){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = mnz+nZBs-1;
                        dvar_vector mPt = mP_z(mnz,mxz);
                        dvector oPt = oP_z(mnz,mxz);
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type="<<ptrZFD->optFit<<cc;
                            cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"ms='ANY_MATURITY'"<<cc;
                            cout<<"sc='ANY_SHELL_CONDITION'"<<cc;
                            cout<<"fit=";
                        }
                        nlls(x,iy) = calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                }//nT>0
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_SEX_MAT_EXTENDED){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    dvector     oP_z(1,nSXs*nMSs*nZBs); oP_z.initialize();//observed size comp.
                    dvar_vector mP_z(1,nSXs*nMSs*nZBs); mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=(nMSs+1);m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            for (int s=1;s<=(nSCs+1);s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (debug>0){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            if (m<=nMSs) {for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);}
                        }//m
                    }//x
                    oP_z /= sum(oP_z);
                    if (debug>0){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type="<<ptrZFD->optFit<<cc;
                                cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"ms='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"sc='ANY_SHELL_CONDITION'"<<cc;
                                cout<<"fit=";
                            }
                            nlls(m+(x-1)*nMSs,iy) = calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                }//nT>0
            }//FIT_BY_SEX_MAT_EXTENDED
        } //if ((mny<=y)&&(y<=mxy))
    } //loop over iy
    if (debug<0) cout<<"NULL)";
    if (debug>0) cout<<"nlls:"<<endl<<nlls<<endl;
    ptrZFD->saveNLLs(nlls);
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_CatchNatZ()"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return nlls;
}

void model_parameters::calcNLLs_Fisheries(int debug, ostream& cout)
{
    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Fisheries()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug>0) cout<<"calculating NLLs for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<ptrMC->lblsFsh[f]<<"=list("<<endl;
        FisheryData* ptrObs = ptrMDS->ppFsh[f-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasN && ptrObs->ptrRCD->ptrN->optFit){
                if (debug>0) cout<<"---retained catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrRCD->ptrN,r_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>0) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrRCD->ptrB,r_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>0) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                dvar_matrix nlls = calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,r_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasN && ptrObs->ptrTCD->ptrN->optFit){
                if (debug>0) cout<<"---total catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrTCD->ptrN,c_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>0) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrTCD->ptrB,c_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>0) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                dvar_matrix nlls = calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,c_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasN && ptrObs->ptrDCD->ptrN->optFit){
                if (debug>0) cout<<"---discard catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrDCD->ptrN,d_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>0) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrDCD->ptrB,d_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>0) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                dvar_matrix nlls = calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,d_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Fisheries()"<<endl;
}

void model_parameters::calcNLLs_Surveys(int debug, ostream& cout)
{
    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Surveys()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug>0) cout<<"calculating NLLs for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<ptrMC->lblsSrv[v]<<"=list("<<endl;
        SurveyData* ptrObs = ptrMDS->ppSrv[v-1];
        if (ptrObs->hasN && ptrObs->ptrN->optFit){
            if (debug>0) cout<<"---survey abundance"<<endl;
            if (debug<0) cout<<"abundance="<<endl;
            dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrN,n_vyxmsz(v),debug,cout);
            objFun += sum(nlls);//TODO: incorporate likelihood weights
            if (debug<0) cout<<","<<endl;
        }
        if (ptrObs->hasB && ptrObs->ptrB->optFit){
            if (debug>0) cout<<"---survey biomass"<<endl;
            if (debug<0) cout<<"biomass="<<endl;
            dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrB,n_vyxmsz(v),debug,cout);
            objFun += sum(nlls);//TODO: incorporate likelihood weights
            if (debug<0) cout<<","<<endl;
        }
        if (ptrObs->hasZFD && ptrObs->ptrZFD->optFit){
            if (debug>0) cout<<"---survey size frequencies"<<endl;
            if (debug<0) cout<<"n.at.z="<<endl;
            dvar_matrix nlls = calcNLLs_CatchNatZ(ptrObs->ptrZFD,n_vyxmsz(v),debug,cout);
            objFun += sum(nlls);//TODO: incorporate likelihood weights
            if (debug<0) cout<<","<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Surveys()"<<endl;
}

void model_parameters::calcAllPriors(int debug, ostream& cout)
{
    if (debug>=dbgPriors) cout<<"Starting calcAllPriors()"<<endl;
    //recruitment parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,  pLnR,  debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRCV,pLnRCV,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLgtRX,pLgtRX,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRa, pLnRa, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRb, pLnRb, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout);
   
    //natural mortality parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnM,   pLnM,   debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMT, pLnDMT, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMX, pLnDMX, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMM, pLnDMM, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMXM,pLnDMXM,debug,cout);
    
    //growth parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrA,   pLnGrA,   debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrB,   pLnGrB,   debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrBeta,pLnGrBeta,debug,cout);
    
    //maturity parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrMat->pLgtPrMat,pLgtPrMat,debug,cout);
    
    //selectivity parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS1,pS1,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS2,pS2,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS3,pS3,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS4,pS4,debug,cout);     
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS5,pS5,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS6,pS6,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS1,devsS1,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS2,devsS2,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS3,devsS3,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS4,devsS4,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS5,devsS5,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS6,devsS6,debug,cout);
    
    //fishing mortality parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC,    pLnC,   debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCT,  pLnDCT, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCX,  pLnDCX, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCM,  pLnDCM, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCXM, pLnDCXM,debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout);
   
    //survey catchability parameters
    tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnQ,    pLnQ,   debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQT,  pLnDQT, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQX,  pLnDQX, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQM,  pLnDQM, debug,cout);
    tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQXM, pLnDQXM,debug,cout);
    
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;
}

void model_parameters::ReportToR_Data(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Data(...)"<<endl;
    ptrMDS->writeToR(os,"data",0);
    if (debug) cout<<"Finished ReportToR_Data(...)"<<endl;
}

void model_parameters::ReportToR_NLLs(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_NLLs(...)"<<endl;
    os<<"ofcs=list("<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_NLLs(...)"<<endl;
}

void model_parameters::ReportToR_Params(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Params(...)"<<endl;
    ptrMPI->writeToR(os);
    if (debug) cout<<"Finished ReportToR_Params(...)"<<endl;
}

void model_parameters::ReportToR_PopQuants(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_PopQuants(...)"<<endl;
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    d5_array n_xmsyz = tcsam::rearrangeYXMSZtoXMSYZ(vn_yxmsz);
    
    os<<"pop.quants=list("<<endl;
    os<<"R.y=";  wts::writeToR(os,value(R_y), ptrMC->csvYrs);               os<<cc<<endl;
    os<<"Rx.y=";  wts::writeToR(os,trans(value(R_yx))(MALE), ptrMC->csvYrs);os<<cc<<endl;
    os<<"Rx.c="; wts::writeToR(os,value(Rx_c),adstring("1:"+str(npcRec)));  os<<cc<<endl;
    os<<"R.cz="; wts::writeToR(os,value(R_cz),adstring("1:"+str(npcRec)),ptrMC->csvZBs); os<<cc<<endl;
    os<<"M.cxm=";   wts::writeToR(os,value(M_cxm),   adstring("1:"+str(npcNM)),ptrMC->csvSXs,ptrMC->csvMSs); os<<cc<<endl;
    os<<"prMat.cz=";wts::writeToR(os,value(prMat_cz),adstring("1:"+str(npcMat)),ptrMC->csvZBs); os<<cc<<endl;
    os<<"prGr.czz=";wts::writeToR(os,value(prGr_czz),adstring("1:"+str(npcGr)),ptrMC->csvZBs,ptrMC->csvZBs); os<<cc<<endl;
    os<<"spb.yx=";  wts::writeToR(os,value(spb_yx),ptrMC->csvYrs,ptrMC->csvSXs); os<<cc<<endl;
    os<<"n.xmsyz="; wts::writeToR(os,n_xmsyz,ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrsP1,ptrMC->csvZBs); os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_PopQuants(...)"<<endl;
}

void model_parameters::ReportToR_SelFuncs(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_SelFuncs(...)"<<endl;
    os<<"sel.funcs=list("<<endl;
    os<<"sel_cz=";  wts::writeToR(os,value(sel_cz),  adstring("1:"+str(npcSel)),ptrMC->csvZBs);os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_SelFuncs(...)"<<endl;
}

void model_parameters::ReportToR_FshQuants(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_FshQuants(...)"<<endl;
    d5_array vFc_fyxms = wts::value(Fc_fyxms);
    d5_array Fc_fxmsy = tcsam::rearrangeIYXMStoIXMSY(vFc_fyxms);
    d6_array vFc_fyxmsz = wts::value(Fc_fyxmsz);
    d6_array Fc_fxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vFc_fyxmsz);
    d6_array vFr_fyxmsz = wts::value(Fr_fyxmsz);
    d6_array Fr_fxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vFr_fyxmsz);
    
    d6_array vc_fyxmsz = wts::value(c_fyxmsz);
    d3_array c_fxy = tcsam::calcIXYfromIYXMSZ(vc_fyxmsz);
    d3_array cb_fxy = tcsam::calcIXYfromIYXMSZ(vc_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array c_fxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vc_fyxmsz);
    
    d6_array vr_fyxmsz = wts::value(r_fyxmsz);
    d3_array r_fxy = tcsam::calcIXYfromIYXMSZ(vr_fyxmsz);
    d3_array rb_fxy = tcsam::calcIXYfromIYXMSZ(vr_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array r_fxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vr_fyxmsz);
   
    
    os<<"fisheries=list("<<endl;
    for (int f=1;f<=nFsh;f++){
        os<<ptrMC->lblsFsh[f]<<"=list("<<endl;
        os<<"tot=list("<<endl;
            os<<"n.xy="; wts::writeToR(os,c_fxy(f), ptrMC->csvSXs,ptrMC->csvYrs); os<<cc<<endl;
            os<<"b.xy="; wts::writeToR(os,cb_fxy(f),ptrMC->csvSXs,ptrMC->csvYrs); os<<cc<<endl;
            os<<"n.xmsyz="; wts::writeToR(os,c_fxmsyz(f),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrs,ptrMC->csvZBs); os<<cc<<endl;
            os<<"Fc.xmsy="; wts::writeToR(os,Fc_fxmsy(f),  ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrs); os<<cc<<endl;
            os<<"F.xmsyz="; wts::writeToR(os,Fc_fxmsyz(f),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrs,ptrMC->csvZBs); 
        os<<")"<<cc<<endl;
        if (sum(r_fxy(f))>0){
            os<<"ret=list("<<endl;
                os<<"n.xy="; wts::writeToR(os,r_fxy(f), ptrMC->csvSXs,ptrMC->csvYrs); os<<cc<<endl;
                os<<"b.xy="; wts::writeToR(os,rb_fxy(f),ptrMC->csvSXs,ptrMC->csvYrs); os<<cc<<endl;
                os<<"n.xmsyz="; wts::writeToR(os,r_fxmsyz(f),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrs,ptrMC->csvZBs); os<<cc<<endl;
                os<<"F.xmsyz="; wts::writeToR(os,Fr_fxmsyz(f),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrs,ptrMC->csvZBs); 
            os<<")"<<cc<<endl;
        }
        os<<"NULL),"<<endl;    
    }
    os<<"NULL)";
    if (debug) cout<<"Finished ReportToR_FshQuants(...)"<<endl;
}

void model_parameters::ReportToR_SrvQuants(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_SrvQuants(...)"<<endl;
    d5_array vq_vyxms = wts::value(q_vyxms);
    d5_array q_vxmsy = tcsam::rearrangeIYXMStoIXMSY(vq_vyxms);
    d6_array vq_vyxmsz = wts::value(q_vyxmsz);
    d6_array q_vxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vq_vyxmsz);
    
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d3_array n_vxy = tcsam::calcIXYfromIYXMSZ(vn_vyxmsz);
    d3_array b_vxy = tcsam::calcIXYfromIYXMSZ(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array n_vxmsyz = tcsam::rearrangeIYXMSZtoIXMSYZ(vn_vyxmsz);
    
    os<<"surveys=list("<<endl;
    for (int v=1;v<=nSrv;v++){
        os<<ptrMC->lblsSrv[v]<<"=list("<<endl;
        os<<"q.xmsy=";  wts::writeToR(os,q_vxmsy(v), ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrsP1); os<<cc<<endl;
        os<<"q.xmsyz="; wts::writeToR(os,q_vxmsyz(v),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrsP1,ptrMC->csvZBs); os<<cc<<endl;
        os<<"n.xy=";    wts::writeToR(os,n_vxy(v),  ptrMC->csvSXs,ptrMC->csvYrsP1); os<<cc<<endl;
        os<<"b.xy=";    wts::writeToR(os,b_vxy(v),  ptrMC->csvSXs,ptrMC->csvYrsP1); os<<cc<<endl;
        os<<"spb.yx=";  wts::writeToR(os,value(spb_vyx(v)),ptrMC->csvYrsP1,ptrMC->csvSXs); os<<cc<<endl;
        os<<"n.xmsyz="; wts::writeToR(os,n_vxmsyz(v),ptrMC->csvSXs,ptrMC->csvMSs,ptrMC->csvSCs,ptrMC->csvYrsP1,ptrMC->csvZBs); os<<endl;
        os<<")"<<cc<<endl;
    }
    os<<"NULL)";
    if (debug) cout<<"Finished ReportToR_SrvQuants(...)"<<endl;
}

void model_parameters::ReportToR_ModelFits(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    os<<"model.fits=list("<<endl;
    os<<"fisheries="; calcNLLs_Fisheries(-1,os); os<<","<<endl; //recalc and write results to os
    os<<"surveys="; calcNLLs_Surveys(-1,os); os<<endl; //recalc and write results to os
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelFits(...)"<<endl;
}

void model_parameters::updateMPI(int debug, ostream& cout)
{
    if (debug) cout<<"Starting updateMPI(...)"<<endl;
    //recruitment parameters
    ptrMPI->ptrRec->pLnR->setFinalVals(pLnR);
    ptrMPI->ptrRec->pLnRCV->setFinalVals(pLnRCV);
    ptrMPI->ptrRec->pLgtRX->setFinalVals(pLgtRX);
    ptrMPI->ptrRec->pLnRa->setFinalVals(pLnRa);
    ptrMPI->ptrRec->pLnRb->setFinalVals(pLnRb);
    //cout<<"setting final vals for pDevsLnR"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) (*ptrMPI->ptrRec->pDevsLnR)[p]->setFinalVals(pDevsLnR(p));
     
    //natural mortality parameters
    ptrMPI->ptrNM->pLnM->setFinalVals(pLnM);
    ptrMPI->ptrNM->pLnDMT->setFinalVals(pLnDMT);
    ptrMPI->ptrNM->pLnDMX->setFinalVals(pLnDMX);
    ptrMPI->ptrNM->pLnDMM->setFinalVals(pLnDMM);
    ptrMPI->ptrNM->pLnDMXM->setFinalVals(pLnDMXM);
    
    //growth parameters
    ptrMPI->ptrGr->pLnGrA->setFinalVals(pLnGrA);
    ptrMPI->ptrGr->pLnGrB->setFinalVals(pLnGrB);
    ptrMPI->ptrGr->pLnGrBeta->setFinalVals(pLnGrBeta);
    
    //maturity parameters
    //cout<<"setting final vals for pLgtPrMat"<<endl;
    for (int p=1;p<=npLgtPrMat;p++) (*ptrMPI->ptrMat->pLgtPrMat)[p]->setFinalVals(pLgtPrMat(p));
    
    //selectivity parameters
    ptrMPI->ptrSel->pS1->setFinalVals(pS1);
    ptrMPI->ptrSel->pS2->setFinalVals(pS2);
    ptrMPI->ptrSel->pS3->setFinalVals(pS3);
    ptrMPI->ptrSel->pS4->setFinalVals(pS4);
    ptrMPI->ptrSel->pS5->setFinalVals(pS5);
    ptrMPI->ptrSel->pS6->setFinalVals(pS6);
    //cout<<"setting final vals for pDevsS1"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS1->getSize();p++) (*ptrMPI->ptrSel->pDevsS1)[p]->setFinalVals(pDevsS1(p));
    //cout<<"setting final vals for pDevsS2"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS2->getSize();p++) (*ptrMPI->ptrSel->pDevsS2)[p]->setFinalVals(pDevsS2(p));
    //cout<<"setting final vals for pDevsS3"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS3->getSize();p++) (*ptrMPI->ptrSel->pDevsS3)[p]->setFinalVals(pDevsS3(p));
    //cout<<"setting final vals for pDevsS4"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS4->getSize();p++) (*ptrMPI->ptrSel->pDevsS4)[p]->setFinalVals(pDevsS4(p));
    //cout<<"setting final vals for pDevsS5"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS5->getSize();p++) (*ptrMPI->ptrSel->pDevsS5)[p]->setFinalVals(pDevsS5(p));
    //cout<<"setting final vals for pDevsS6"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS6->getSize();p++) (*ptrMPI->ptrSel->pDevsS6)[p]->setFinalVals(pDevsS6(p));
     
    //fully-selected fishing capture rate parameters
    ptrMPI->ptrFsh->pHM->setFinalVals(pHM);
    ptrMPI->ptrFsh->pLnC->setFinalVals(pLnC);
    ptrMPI->ptrFsh->pLnDCT->setFinalVals(pLnDCT);
    ptrMPI->ptrFsh->pLnDCX->setFinalVals(pLnDCX);
    ptrMPI->ptrFsh->pLnDCM->setFinalVals(pLnDCM);
    ptrMPI->ptrFsh->pLnDCXM->setFinalVals(pLnDCXM);
    //cout<<"setting final vals for pDevsLnC"<<endl;
    for (int p=1;p<=ptrMPI->ptrFsh->pDevsLnC->getSize();p++) (*ptrMPI->ptrFsh->pDevsLnC)[p]->setFinalVals(pDevsLnC(p));
    
    //survey catchability parameters
    ptrMPI->ptrSrv->pLnQ->setFinalVals(pLnQ);
    ptrMPI->ptrSrv->pLnDQT->setFinalVals(pLnDQT);
    ptrMPI->ptrSrv->pLnDQX->setFinalVals(pLnDQX);
    ptrMPI->ptrSrv->pLnDQM->setFinalVals(pLnDQM);
    ptrMPI->ptrSrv->pLnDQXM->setFinalVals(pLnDQXM);
    
    if (debug) cout<<"Finished updateMPI(...)"<<endl;
}

void model_parameters::ReportToR(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR(...)"<<endl;
    updateMPI(debug,cout);
        
    os<<"res=list("<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        
        //model data
        ptrMDS->writeToR(os,"data",0); os<<","<<endl;
        
        //objective function values
        ReportToR_NLLs(os,debug,cout); os<<","<<endl;
        
        //parameter values
        ReportToR_Params(os,debug,cout); os<<","<<endl;
        //selectivity functions
        ReportToR_SelFuncs(os,debug,cout); os<<","<<endl;
        //population quantities
        ReportToR_PopQuants(os,debug,cout); os<<","<<endl;
        //fishery quantities
        ReportToR_FshQuants(os,debug,cout); os<<","<<endl;
        //survey quantities 
        ReportToR_SrvQuants(os,debug,cout); os<<","<<endl;
        //model fit quantities
        ReportToR_ModelFits(os,debug,cout); os<<","<<endl;
        
        //simulated model data
        createSimData(debug, cout);
        ptrSimMDS->writeToR(os,"sim.data",0); 
        os<<endl;
    os<<")"<<endl;
    if (debug) cout<<"Finished ReportToR(...)"<<endl;
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
        
    //write active parameters to rpt::echo
    rpt::echo<<"Finished phase "<<current_phase()<<endl;
    
    //write report as R file
    ReportToR(report,0,rpt::echo);
}

void model_parameters::between_phases_calculations(void)
{
}

void model_parameters::final_calcs()
{
    {cout<<"writing model sim data to file"<<endl;
        ofstream echo1; echo1.open("ModelSimData.dat", ios::trunc);
        writeSimData(echo1,0,cout);
    }
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::app);
        mcmc<<"NULL)"<<endl;
        mcmc.close();
    }
    
    long hour,minute,second;
    double elapsed_time;
    
    time(&finish); 
    elapsed_time = difftime(finish,start);
    
    hour = long(elapsed_time)/3600;
    minute = long(elapsed_time)%3600/60;
    second = (long(elapsed_time)%3600)%60;
    cout << endl << endl << "Starting time: " << ctime(&start);
    cout << "Finishing time: " << ctime(&finish);
    cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
    
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{1000,5000,5000,5000,5000,5000,10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1,1,.01,.001,.001,.001,1e-3,1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 1000000000; //must be smaller than 2,147,483,647
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(40000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(4000);
  time(&start);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
