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
    
    //runtime flags (0=false)
    int jitter     = 0;//use jittering for initial parameter values
    int resample   = 0;//use resampling for initial parameter values
    int opModMode  = 0;//run as operating model, no fitting
    int usePin     = 0;//flag to initialize parameter values using a pin file
    int doRetro    = 0;//flag to facilitate a retrospective model run
    int fitSimData = 0;//flag to fit model to simulated data calculated in the PRELIMINARY_CALCs section
    
    int yRetro = 0; //number of years to decrement for retrospective model run
    int iSeed =  0;//default random number generator seed
    random_number_generator rng(-1);//random number generator
    int iSimDataSeed = 0;
    random_number_generator rngSimData(-1);//random number generator for data simulation
    
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
    
    int nSXs    = tcsam::nSXs;
    int MALE    = tcsam::MALE;
    int FEMALE  = tcsam::FEMALE;
    int ALL_SXs = tcsam::ALL_SXs;
    
    int nMSs     = tcsam::nMSs;
    int IMMATURE = tcsam::IMMATURE;
    int MATURE   = tcsam::MATURE;
    int ALL_MSs  = tcsam::ALL_MSs;
    
    int nSCs      = tcsam::nSCs;
    int NEW_SHELL = tcsam::NEW_SHELL;
    int OLD_SHELL = tcsam::OLD_SHELL;
    int ALL_SCs   = tcsam::ALL_SCs;
    
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
    rpt::echo<<"#------Reading command line options---------"<<endl;
    //configFile
    fnConfigFile = "TCSAM2015_ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //resultsFile
    fnResultFile = "TCSAM2015_ResultFile";//default results file name
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resultsFile"))>-1) {
        fnResultFile = ad_comm::argv[on+1];
        rpt::echo<<"#results file changed to '"<<fnResultFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //opModMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-opModMode"))>-1) {
        opModMode=1;
        rpt::echo<<"#operating model mode turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        rpt::echo<<"#debugModelDatasets turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParamsInfo
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParamsInfo"))>-1) {
        debugModelParamsInfo=1;
        rpt::echo<<"#debugModelParamsInfo turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //doRetro
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-doRetro"))>-1) {
        doRetro=1;
        cout<<"#doRetro turned ON"<<endl;
        rpt::echo<<"#doRetro turned ON"<<endl;
        if (on+1<argc) {
            yRetro=atoi(ad_comm::argv[on+1]);
            cout<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
            rpt::echo<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //fitSimData
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-fitSimData"))>-1) {
        fitSimData=1;
        if (on+1<argc) {
            iSimDataSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed (0 -> deterministic) for data simulation: ";
            cin>>iSimDataSeed;
        }
        if (iSimDataSeed) rng.reinitialize(iSimDataSeed);
        rpt::echo<<"#Simulating data to fit using "<<iSimDataSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //seed
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-seed"))>-1) {
        if (on+1<argc) {
            iSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed for jittering/resampling: ";
            cin>>iSeed;
        }
        rng.reinitialize(iSeed);
        rpt::echo<<"#Random number seed set to "<<iSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //jitter
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-jitter"))>-1) {
        jitter=1;
        rpt::echo<<"#Jittering for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        rpt::echo<<"#Resampling for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParams"))>-1) {
        debugModelParams=1;
        cout<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        rpt::echo<<"#debugDATA_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPARAMS_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPARAMS_SECTION"))>-1) {
        debugPARAMS_SECTION=1;
        rpt::echo<<"#debugPARAMS_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPRELIM_CALCS
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPRELIM_CALCS"))>-1) {
        debugPRELIM_CALCS=1;
        rpt::echo<<"debugPRELIM_CALCS turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPROC_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPROC_SECTION"))>-1) {
        debugPROC_SECTION=1;
        rpt::echo<<"#debugPROC_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugREPORT_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugREPORT_SECTION"))>-1) {
        debugREPORT_SECTION=1;
        rpt::echo<<"#debugREPORT_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugRunModel
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugRunModel")>-1) {
        debugRunModel=1;
        rpt::echo<<"#debugRunModel turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugObjFun
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugObjFun")>-1) {
        debugObjFun=1;
        rpt::echo<<"#debugObjFun turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //showActiveParams
    if (option_match(ad_comm::argc,ad_comm::argv,"-showActiveParams")>-1) {
        showActiveParams=1;
        rpt::echo<<"#showActiveParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading configuration file '"<<fnConfigFile<<"'"<<endl;
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    if (doRetro){mxYr = mxYr-yRetro; ptrMC->setMaxModelYear(mxYr);}
    if (jitter)   {ptrMC->jitter=1;}
    if (resample) {ptrMC->resample = 1;}
    
    rpt::echo<<"#------------------ModelConfiguration-----------------"<<endl;
    rpt::echo<<(*ptrMC);
    rpt::echo<<"#-----------------------------------"<<endl;
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
    rpt::echo<<"---SimMDS object after reading datasets---"<<endl;
    rpt::echo<<(*ptrSimMDS);
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
  nmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"nmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    nmN_yxmsz.initialize();
  #endif
  tmN_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmN_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmN_yxmsz.initialize();
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
  stdvDevsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"stdvDevsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    stdvDevsLnR_cy.initialize();
  #endif
  zscrDevsLnR_cy.allocate(1,npcRec,mnYr,mxYr,"zscrDevsLnR_cy");
  #ifndef NO_AD_INITIALIZE
    zscrDevsLnR_cy.initialize();
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
  cF_fyxms.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,"cF_fyxms");
  #ifndef NO_AD_INITIALIZE
    cF_fyxms.initialize();
  #endif
  cF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cF_fyxmsz.initialize();
  #endif
  rmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmF_fyxmsz.initialize();
  #endif
  dmF_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmF_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmF_fyxmsz.initialize();
  #endif
  tmF_yxmsz.allocate(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"tmF_yxmsz");
  #ifndef NO_AD_INITIALIZE
    tmF_yxmsz.initialize();
  #endif
  cN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"cN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    cN_fyxmsz.initialize();
  #endif
  dN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dN_fyxmsz.initialize();
  #endif
  rmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"rmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    rmN_fyxmsz.initialize();
  #endif
  dmN_fyxmsz.allocate(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs,"dmN_fyxmsz");
  #ifndef NO_AD_INITIALIZE
    dmN_fyxmsz.initialize();
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
  nllRecDevs.allocate(1,npcRec,"nllRecDevs");
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
    
    //set initial values for all parameters
    if (usePin) {
        rpt::echo<<"NOTE: setting initial values for parameters using pin file"<<endl;
    } else {
        rpt::echo<<"NOTE: setting initial values for parameters using setInitVals(...)"<<endl;
        setInitVals();
    }
    cout<<"testing setAllDevs()"<<endl;
    setAllDevs(0,rpt::echo);
        
    {cout<<"writing data to R"<<endl;
     ofstream echo1; echo1.open("ModelData.R", ios::trunc);
     ReportToR_Data(echo1,0,cout);
    }
    
    {cout<<"writing parameters info to R"<<endl;
     ofstream echo1; echo1.open("ModelParametersInfo.R", ios::trunc);
     ptrMPI->writeToR(echo1);
    }
    
    //calculate average effort for fisheries over specified time periods
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
    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")<0) {
        cout<<"testing calcRecruitment():"<<endl;
        calcRecruitment(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcNatMort():"<<endl;
        calcNatMort(0,rpt::echo);
        rpt::echo<<"testing calcGrowth():"<<endl;
        calcGrowth(0,rpt::echo);
        rpt::echo<<"testing calcMaturity():"<<endl;
        calcMaturity(0,rpt::echo);
        rpt::echo<<"testing calcSelectivities():"<<endl;
        calcSelectivities(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcFisheryFs():"<<endl;
        calcFisheryFs(dbgCalcProcs+1,rpt::echo);
        rpt::echo<<"testing calcSurveyQs():"<<endl;
        calcSurveyQs(dbgCalcProcs+1,cout);
        rpt::echo<<"testing runPopDyMod():"<<endl;
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
        
        if (fitSimData){
            cout<<"creating sim data to fit in model"<<endl;
            createSimData(1,rpt::echo,iSimDataSeed,ptrMDS);//stochasstic if iSimDataSeed<>0
            {cout<<"re-writing data to R"<<endl;
             ofstream echo1; echo1.open("ModelData.R", ios::trunc);
             ReportToR_Data(echo1,0,cout);
            }
        }
        
        cout<<"Testing calcObjFun()"<<endl;
        calcObjFun(1,rpt::echo);
        {cout<<"writing model results to R"<<endl;
            ofstream echo1; echo1.open("ModelRes0.R", ios::trunc);
            ReportToR(echo1,1,cout);
        }
        
        {cout<<"writing model sim data to file"<<endl;
            createSimData(1,rpt::echo,0,ptrSimMDS);//deterministic
            ofstream echo1; echo1.open("ModelSimData0.dat", ios::trunc);
            writeSimData(echo1,0,rpt::echo,ptrSimMDS);
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

void model_parameters::setInitVals(void)
{
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
        ivector bnds = wts::getBounds(spb_yx);
        adstring dsxs = tcsamDims::getSXsForR(bnds[3],bnds[4]);
        mcmc<<"spb_xy="; wts::writeToR(mcmc,trans(value(spb_yx)),dsxs,ptrMC->dimYrsToR); //mcmc<<cc<<endl;
        
    mcmc<<")"<<cc<<endl;
    mcmc.close();
    
}

void model_parameters::createSimData(int debug, ostream& cout, int iSimDataSeed, ModelDatasets* ptrSim)
{
    if (debug)cout<<"simulating model results as data"<<endl;
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vcN_fyxmsz = wts::value(cN_fyxmsz);
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    for (int v=1;v<=nSrv;v++) {
        if (debug) cout<<"survey "<<v<<endl;
        (ptrSim->ppSrv[v-1])->replaceCatchData(iSimDataSeed,rngSimData,vn_vyxmsz(v),ptrSim->ptrBio->wAtZ_xmz);
    }
    for (int f=1;f<=nFsh;f++) {
        if (debug) cout<<"fishery f: "<<f<<endl;
        (ptrSim->ppFsh[f-1])->replaceCatchData(iSimDataSeed,rngSimData,vcN_fyxmsz(f),vrmN_fyxmsz(f),ptrSim->ptrBio->wAtZ_xmz);
    }
    if (debug) cout<<"finished simulating model results as data"<<endl;
     
}

void model_parameters::writeSimData(ostream& os, int debug, ostream& cout, ModelDatasets* ptrSim)
{
    if (debug)cout<<"writing model results as data"<<endl;
    for (int v=1;v<=nSrv;v++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppSrv[v-1]))<<endl;
    }
    //     cout<<4<<endl;
    for (int f=1;f<=nFsh;f++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppFsh[f-1]))<<endl;
    }
    if (debug) cout<<"finished writing model results as data"<<endl;
     
}

void model_parameters::setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = pI->getInitVals();//defaults are initial values
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            NumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);//assign initial value based on resampling prior pdf
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
     
}

void model_parameters::setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
{
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = 0.5*(pI->getUpperBounds()+pI->getLowerBounds());//defaults are midpoints of ranges
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            BoundedNumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<"jittering "<<p(i).label()<<endl;
                p(i) = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//only done if parameter phase > 0
            } else 
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
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
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            BoundedVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
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
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            DevsVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
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
    nmN_yxmsz.initialize();
    tmN_yxmsz.initialize();
       
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
        //apply natural mortality from fisheries to molting/growth/maturity
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
        //apply natural mortality from molting/growth/maturity to fisheries
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
                n1_xmsz(x,m,s) = elem_prod(exp(-M_yxmz(y,x,m)*dt),n0_xmsz(x,m,s));//survivors
                nmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
                tmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
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
    dvar_vector tm_z(1,nZBs);//total mortality (numbers) by size
    dvar_vector tvF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvector     tdF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                tmF_yxmsz(y,x,m,s) = 0.0;//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) tmF_yxmsz(y,x,m,s) += rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(mfexp(-tmF_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheries
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);  //numbers killed by all fisheries
                tmN_yxmsz(y,x,m,s) += tm_z;            //add in numbers killed by all fisheries to total killed
                
                //calculate fishing rate components (need to ensure NOT dividing by 0)
                tdF_z = value(tmF_yxmsz(y,x,m,s));
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),tmF_yxmsz(y,x,m,s)) + 
                                  wts::isEQ(tdF_z,0.0);
                for (int f=1;f<=nFsh;f++){                   
                    cN_fyxmsz(f,y,x,m,s)  = elem_prod(elem_div( cF_fyxmsz(f,y,x,m,s),tvF_z),tm_z);//numbers captured in fishery f
                    rmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(rmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z);//retained mortality in fishery f (numbers)
                    dmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(dmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z);//discards mortality in fishery f (numbers)
                    dN_fyxmsz(f,y,x,m,s)  = cN_fyxmsz(f,y,x,m,s)-rmN_fyxmsz(f,y,x,m,s);//discarded catch (NOT mortality) in fishery f (numbers)                    
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
    stdvDevsLnR_cy.initialize();
    zscrDevsLnR_cy.initialize();
    
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
            if ((mnYr<=y)&&(y<=mxYr)){
                if (debug>dbgCalcProcs) cout<<"y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
                if (useDevs){
                    R_y(y) = mfexp(mnLnR+dvsLnR[idxDevsLnR[y]]);
                } else {
                    R_y(y) = mdR;
                }
                if (debug>dbgCalcProcs) cout<<R_y(y)<<tb;
                R_yx(y,MALE)   = Rx_c(pc);
                if (debug>dbgCalcProcs) cout<<R_yx(y,MALE)<<endl;
                if (FEMALE<=nSXs) R_yx(y,FEMALE) = 1.0-R_yx(y,MALE);
                R_yz(y) = R_cz(pc);
                if (debug>dbgCalcProcs) cout<<R_yz(y)<<endl;
                stdvDevsLnR_cy(pc,y) = sqrt(log(1.0+mfexp(2.0*lnRCV)));           //ln-scale std dev
                zscrDevsLnR_cy(pc,y) = dvsLnR(idxDevsLnR(y))/stdvDevsLnR_cy(pc,y);//standardized ln-scale rec devs
            } else {
                if (debug>dbgCalcProcs) cout<<"skipping y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            }
        }//idx
    }//pc
    
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz = "<<endl<<R_yz<<endl;
        cout<<"zscr = "<<zscrDevsLnR_cy<<endl;
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
        if (FEMALE<=nSXs){
            //add in female offset
            if (pids[k]) {lnM(FEMALE) += pLnDMX(pids[k]);}                      k++;
            //add in immature offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnM(x,IMMATURE) += pLnDMM(pids[k]);} k++;
            //add in offset immature females for stanza
            if (pids[k]) {lnM(FEMALE,IMMATURE) += pLnDMXM(pids[k]);}            k++; //advance k to zScaling in pids
        }
        
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
            if ((mnYr<=y)&&(y<=mxYr)){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        M_yxmz(y,x,m) = M_xmz(x,m);
                    }
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
        if (debug>dbgCalcProcs) cout<<"maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2);
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(x)<<endl;
                prMat_yxz(y,x) = prMat_cz(pc);//note: this change made a difference, but not sure why!
            }
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
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2); //sex index
                for (int m=1;m<=nMSs;m++){
                    for (int z=1;z<=nZBs;z++){
                        prGr_yxmzz(y,x,m,z) = prGr_czz(pc,z);
                    }
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
            if ((mnYr<=y)&&(y<=mxYr+1)){
                k=ptrSel->nIVs+1+6;//1st devs vector variable column
                if (useDevsS1) {params[1] += devsS1(useDevsS1,idxDevsS1[y]);}
                if (useDevsS2) {params[2] += devsS2(useDevsS2,idxDevsS2[y]);}
                if (useDevsS3) {params[3] += devsS3(useDevsS3,idxDevsS3[y]);}
                if (useDevsS4) {params[4] += devsS4(useDevsS4,idxDevsS4[y]);}
                if (useDevsS5) {params[5] += devsS5(useDevsS5,idxDevsS5[y]);}
                if (useDevsS6) {params[6] += devsS6(useDevsS6,idxDevsS6[y]);}
                sel_iyz(pc,y) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"sel: "<<sel_iyz(pc,y)<<endl;
            } else {
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
            }
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
    cF_fyxms.initialize(); //fully-selected capture rate
    cF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retention rate
    dmF_fyxmsz.initialize();//discard mortality rate
    tmF_yxmsz.initialize(); //total mortality rate
    
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
            if (FEMALE<=nSXs){
                //add in female offset
                if (pids[k]) {lnC(FEMALE) += pLnDCX(pids[k]);}                      k++;
                //add in immature offsets
                if (pids[k]) {for (int x=1;x<=nSXs;x++) lnC(x,IMMATURE) += pLnDCM(pids[k]);} k++;
                //add in offset immature females for stanza
                if (pids[k]) {lnC(FEMALE,IMMATURE) += pLnDCXM(pids[k]);}            k++; 
            }
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
                if ((mnYr<=y)&&(y<=mxYr)){
                    x = idxs(idx,3);//sex
                    if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                    if (useDevs) C_xm = mfexp(lnC+dvsLnC[idxDevsLnC[y]]);//recalculate C_xm w/ devs
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            cF_fyxms(f,y,x,m,s)  = C_xm(x,m);                 //fully-selected capture rate
                            cF_fyxmsz(f,y,x,m,s) = C_xm(x,m)*sel_iyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         cF_fyxmsz(f,y,x,m,s));//retention mortality
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),cF_fyxmsz(f,y,x,m,s));//discard mortality
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cF_fyxmsz(f,y,x,m,s);//discard mortality
                            }
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
                                for (int y=mny;y<=mxy;y++) tot += cF_fyxms(f,y,x,m,s);
                                avgFc(f,x,m,s) = tot/(mxy-mny+1); break;
                            case 2:
                                for (int y=mny;y<=mxy;y++) tot += 1.0-mfexp(-cF_fyxms(f,y,x,m,s));
                                avgFc(f,x,m,s) = tot/(mxy-mny+1); break;
                            case 3:
                                for (int y=mny;y<=mxy;y++) tot += mean(cF_fyxmsz(f,y,x,m,s));
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
                if ((mnYr<=y)&&(y<=mxYr)){
                    x = idxs(idx,3);//sex
                    fd = mapM2DFsh(f);//index of corresponding fishery data object
                    eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                    if (debug>dbgCalcProcs) cout<<"f,y,x,eff = "<<f<<tb<<y<<tb<<x<<tb<<eff<<endl;
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            //fully-selected capture rate
                            switch(optsFcAvg(f)) {
                                case 1:
                                    cF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                                case 2:
                                    cF_fyxms(f,y,x,m,s) = -log(1.0-avgRatioFc2Eff(f,x,m,s)*eff); break;
                                case 3:
                                    cF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                            }
                            cF_fyxmsz(f,y,x,m,s) = cF_fyxms(f,y,x,m,s)*sel_iyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_iyz(idRet,y),         cF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_iyz(idRet,y)),cF_fyxmsz(f,y,x,m,s));//discard mortality rate
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cF_fyxmsz(f,y,x,m,s);//discard mortality rate
                            }
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
        if (FEMALE<=nSXs){
            //add in female offset
            if (pids[k]) {lnQ(FEMALE) += pLnDQX(pids[k]);}                      k++;
            //add in immature offsets
            if (pids[k]) {for (int x=1;x<=nSXs;x++) lnQ(x,IMMATURE) += pLnDQM(pids[k]);} k++;
            //add in offset immature females for stanza
            if (pids[k]) {lnQ(FEMALE,IMMATURE) += pLnDQXM(pids[k]);}            k++; 
        }
        
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
            if ((mnYr<=y)&&(y<=mxYrp1)){
                x = idxs(idx,3);//sex
                if (y <= mxYr+1){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            q_vyxms(v,y,x,m,s)  = Q_xm(x,m);
                            q_vyxmsz(v,y,x,m,s) = Q_xm(x,m)*sel_iyz(idSel,y);
                        }
                    }
                }
            }
        }
    }
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
    
}

void model_parameters::calcPenalties(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    if (debug<0) cout<<"list("<<endl;//start list of penalties by category
    if (debug<0) cout<<tb<<"maturity=list("<<endl;//start of maturity penalties list
    //smoothness penalties on maturity parameters (NOT maturity ogives)
    double penWgtLgtPrMat = 1.0;//TODO: read in value from input file
    fPenSmoothLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = 1.0*pLgtPrMat(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtLgtPrMat*fPenSmoothLgtPrMat(i)<<"),"<<endl;
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = 1.0*pLgtPrMat(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtLgtPrMat*fPenSmoothLgtPrMat(i)<<")"<<endl;
    }
    if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list
    //non-decreasing penalties on maturity parameters (NOT maturity ogives)
    double penWgtNonDecLgtPrMat = 1.0;//TODO: read in value from input file
    fPenNonDecLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"nondecreasing=list(";//start of non-decreasing penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = calc1stDiffs(pLgtPrMat(i));
        for (int iv=v.indexmin();iv<=v.indexmax();iv++){
            posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<"),";
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = calc1stDiffs(pLgtPrMat(i));
        for (int iv=v.indexmin();iv<=v.indexmax();iv++){
            posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<")";
    }
    if (debug<0) cout<<tb<<tb<<")"<<endl;//end of non-decreasing penalties list    
    if (debug<0) cout<<tb<<")";//end of maturity penalties list
    
    if (debug<0) cout<<")";//end of penalties list
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

void model_parameters::calcNLLs_Recruitment(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Starting calcNLLs_Recruitment"<<endl;
    double nllWgtRecDevs = 1.0;//TODO: read in from input file (as vector?))
    nllRecDevs.initialize();
    if (debug<0) cout<<"list("<<endl;
    if (debug<0) cout<<tb<<"recDevs=list("<<endl;
    for (int pc=1;pc<npcRec;pc++){
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<cc<<endl;
        }
    }//pc
    {
        int pc = npcRec;
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<endl;
        }
    }//pc
    if (debug<0) cout<<tb<<")";//recDevs
    if (debug<0) cout<<")";
    if (debug>=dbgObjFun) cout<<"Finished calcNLLs_Recruitment"<<endl;
}

void model_parameters::calcObjFun(int debug, ostream& cout)
{
    if (debug>=dbgObjFun) cout<<"Starting calcObjFun"<<endl;
    //objective function penalties
    calcPenalties(debug,cout);
    //prior likelihoods
    calcAllPriors(debug,cout);
    //recruitment component
    calcNLLs_Recruitment(debug,cout);
    
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
        double wgt = 1.0;//TODO: implement likelihood weights
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='norm2',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc; 
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
        double wgt = 1.0;//TODO: implement likelihood weights
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc; 
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
        double wgt = 1.0;//TODO: implement likelihood weights
        adstring yrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='lognormal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc; 
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
    if (debug<0) cout<<"list(fit.type='"<<tcsam::getFitType(ptrA->optFit)<<"',fits=list("<<endl;
    if (ptrA->optFit==tcsam::FIT_BY_X){
        nlls.allocate(1,nSXs); nlls.initialize();
        tA_xy.allocate(1,nSXs,mny,mxy); tA_xy.initialize();
        ivector yrs = ptrA->yrs;
        for (int x=1;x<=nSXs;x++){
            for (int y=mny;y<=mxy;y++) tA_xy(x,y) = sum(mA_yxmsz(y,x));//sum by sex
            if (debug<0) cout<<tcsam::getSexType(x)<<"=";
            nlls(x) = calcNLL(ptrA->llType, tA_xy(x), ptrA->C_xy(x), ptrA->sd_xy(x), ptrA->yrs, debug, cout); 
            if (debug<0) cout<<","<<endl;
        }
        if (debug<0) cout<<"NULL)";
    } else 
    if (ptrA->optFit==tcsam::FIT_BY_TOT){
        nlls.allocate(ALL_SXs,ALL_SXs); nlls.initialize();
        tA_xy.allocate(ALL_SXs,ALL_SXs,mny,mxy);
        tA_xy.initialize();
        for (int y=mny;y<=mxy;y++) tA_xy(ALL_SXs,y) = sum(mA_yxmsz(y));//sum over sexes
        if (debug<0) cout<<tcsam::getSexType(ALL_SXs)<<"=";
        nlls(ALL_SXs) = calcNLL(ptrA->llType, tA_xy(ALL_SXs), ptrA->C_xy(ALL_SXs), ptrA->sd_xy(ALL_SXs), ptrA->yrs, debug, cout);                
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
    if (debug<0) cout<<"list(fit.type='"<<tcsam::getFitType(ptrB->optFit)<<"',fits=list("<<endl;
    if (ptrB->optFit==tcsam::FIT_BY_X){
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
            nlls(x) = calcNLL(ptrB->llType, tB_xy(x), ptrB->C_xy(x), ptrB->sd_xy(x), ptrB->yrs, debug, cout); 
            if (debug<0) cout<<","<<endl;
        }
        if (debug<0) cout<<"NULL)";
    } else 
    if (ptrB->optFit==tcsam::FIT_BY_TOT){
        nlls.allocate(ALL_SXs,ALL_SXs); nlls.initialize();
        tB_xy.allocate(ALL_SXs,ALL_SXs,mny,mxy);
        tB_xy.initialize();
        //calc catch biomass over sexes
        for (int x=1;x<=nSXs;x++){
            for (int y=mny;y<=mxy;y++) {
                for (int m=1;m<=nMSs;m++) {
                    for (int s=1;s<=nSCs;s++) tB_xy(ALL_SXs,y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                }
            }
        }
        if (debug<0) cout<<tcsam::getSexType(ALL_SXs)<<"=";
        nlls(ALL_SXs) = calcNLL(ptrB->llType, tB_xy(ALL_SXs), ptrB->C_xy(ALL_SXs), ptrB->sd_xy(ALL_SXs), ptrB->yrs, debug, cout);                
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
        double wgt = 1.0;//TODO: incorporate weights
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(nll.type='multinomial',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
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

void model_parameters::calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
{
    if (debug>=dbgAll) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    if (ptrZFD->optFit==tcsam::FIT_NONE) return;
    ivector yrs = ptrZFD->yrs;
    int y;
    double ss;
    dvariable nT;
    dvariable nlls;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvector     oP_z;//observed size comp.
    dvar_vector mP_z;//model size comp.
    if (ptrZFD->optFit==tcsam::FIT_BY_XE){
        oP_z.allocate(1,nSXs*nZBs);
        mP_z.allocate(1,nSXs*nZBs);
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XME){
        oP_z.allocate(1,nSXs*nMSs*nZBs);
        mP_z.allocate(1,nSXs*nMSs*nZBs);
    } else {
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
    }
    if (debug<0) cout<<"list("<<endl;
    for (int iy=1;iy<=yrs.size();iy++) {
        y = yrs[iy];
        if (debug>0) cout<<"y = "<<y<<endl;
        if ((mny<=y)&&(y<=mxy)) {
            if (ptrZFD->optFit==tcsam::FIT_BY_TOT){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=ALL_SXs;x++){
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
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
                        cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                        cout<<"yr="<<y<<cc;
                        cout<<"sx='ALL_SEX'"<<cc;
                        cout<<"ms='ALL_MATURITY'"<<cc;
                        cout<<"sc='ALL_SHELL_CONDITION'"<<cc;
                        cout<<"fit=";
                    }
                    nlls = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                    objFun += nlls;//TODO: add in likelihood weights
                    if (debug<0) cout<<")"<<cc<<endl;
                }
                //FIT_BY_TOT
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
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
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"yr="<<y<<cc;
                            cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"ms='ALL_MATURITY'"<<cc;
                            cout<<"sc='ALL_SHELL_CONDITION'"<<cc;
                            cout<<"fit=";
                        }
                        nlls = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                        objFun += nlls;//TODO: add in likelihood weights
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
                //FIT_BY_X
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XE){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
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
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"yr="<<y<<cc;
                            cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"ms='ALL_MATURITY'"<<cc;
                            cout<<"sc='ALL_SHELL_CONDITION'"<<cc;
                            cout<<"fit=";
                        }
                        nlls = calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                        objFun += nlls;//TODO: add in likelihood weights
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                }//nT>0
                //FIT_BY_XE
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_XM){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            oP_z /= sum(oP_z);
                            if (debug>0){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"yr="<<y<<cc;
                                cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"ms='"<<tcsam::getMaturityType(x)<<"'"<<cc;
                                cout<<"sc='ALL_SHELL_CONDITION'"<<cc;
                                cout<<"fit=";
                            }
                            nlls = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                            objFun += nlls;//TODO: add in likelihood weights
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XM
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XME){
                ss = 0;
                nT = sum(mA_yxmsz[y]);//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            for (int s=1;s<=ALL_SCs;s++) {
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
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"yr="<<y<<cc;
                                cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"ms='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"sc='ALL_SHELL_CONDITION'"<<cc;
                                cout<<"fit=";
                            }
                            nlls = calcNLL(ptrZFD->llType,mPt,oPt,ss,debug,cout);
                            objFun += nlls;//TODO: add in likelihood weights
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                    //FIT_BY_XME
                }//nT>0
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XS){
                for (int x=1;x<=nSXs;x++) {
                    for (int s=1;s<=nSCs;s++){
                        ss = 0;
                        nT.initialize();
                        for (int m=1;m<=ALL_MSs;m++) nT += sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            oP_z /= sum(oP_z);
                            if (debug>0){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int m=1;m<=nMSs;m++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"yr="<<y<<cc;
                                cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"ms='ALL_MATURITY'"<<cc;
                                cout<<"sc='"<<tcsam::getShellType(x)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            nlls = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                            objFun += nlls;//TODO: add in likelihood weights
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XS
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XMS){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) {
                            ss = 0;
                            nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                            if (value(nT)>0){
                                oP_z.initialize();//observed size comp.
                                mP_z.initialize();//model size comp.                            
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                oP_z /= sum(oP_z);
                                if (debug>0){
                                    cout<<"ss = "<<ss<<endl;
                                    cout<<"oP_Z = "<<oP_z<<endl;
                                }
                                mP_z += mA_yxmsz(y,x,m,s);
                                mP_z /= nT;//normalize model size comp
                                if (debug>0) cout<<"mP_z = "<<mP_z<<endl;
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"yr="<<y<<cc;
                                    cout<<"sx='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"ms='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"sc='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                nlls = calcNLL(ptrZFD->llType,mP_z,oP_z,ss,debug,cout);
                                objFun += nlls;//TODO: add in likelihood weights
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//nT>0
                        }//s
                    }//m
                }//x
                //FIT_BY_XMS
            } else 
            {
                std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
                std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
                std::cout<<"Aborting..."<<endl;
                exit(-1);
            }
        } //if ((mny<=y)&&(y<=mxy))
    } //loop over iy
    if (debug<0) cout<<"NULL)";
    if (debug>0) cout<<"nlls:"<<endl<<nlls<<endl;
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_CatchNatZ()"<<endl;
    }
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
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrRCD->ptrN,rmN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>0) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrRCD->ptrB,rmN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>0) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasN && ptrObs->ptrTCD->ptrN->optFit){
                if (debug>0) cout<<"---total catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrTCD->ptrN,cN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>0) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrTCD->ptrB,cN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>0) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasN && ptrObs->ptrDCD->ptrN->optFit){
                if (debug>0) cout<<"---discard catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                dvar_vector nlls = calcNLLs_CatchAbundance(ptrObs->ptrDCD->ptrN,dN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>0) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                dvar_vector nlls = calcNLLs_CatchBiomass(ptrObs->ptrDCD->ptrB,dN_fyxmsz(f),debug,cout);
                objFun += sum(nlls);//TODO: incorporate likelihood weights
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>0) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dN_fyxmsz(f),debug,cout);
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
            calcNLLs_CatchNatZ(ptrObs->ptrZFD,n_vyxmsz(v),debug,cout);
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
    if (debug<0) cout<<"list("<<endl;
    //recruitment parameters
    if (debug) cout<<tb<<"recruitment=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,  pLnR,  debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRCV,pLnRCV,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLgtRX,pLgtRX,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRa, pLnRa, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRb, pLnRb, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
   
    //natural mortality parameters
    if (debug) cout<<tb<<"'natural mortality'=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnM,   pLnM,   debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMT, pLnDMT, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMX, pLnDMX, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMM, pLnDMM, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMXM,pLnDMXM,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
    
    //growth parameters
    if (debug) cout<<tb<<"growth=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrA,   pLnGrA,   debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrB,   pLnGrB,   debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGr->pLnGrBeta,pLnGrBeta,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
    
    //maturity parameters
    if (debug) cout<<tb<<"maturity=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrMat->pLgtPrMat,pLgtPrMat,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
    
    //selectivity parameters
    if (debug) cout<<tb<<"'selectivity functions'=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS1,pS1,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS2,pS2,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS3,pS3,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS4,pS4,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS5,pS5,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS6,pS6,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS1,devsS1,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS2,devsS2,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS3,devsS3,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS4,devsS4,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS5,devsS5,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS6,devsS6,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
    
    //fishing mortality parameters
    if (debug) cout<<tb<<"fisheries=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC,    pLnC,   debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCT,  pLnDCT, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCX,  pLnDCX, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCM,  pLnDCM, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCXM, pLnDCXM,debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<cc<<endl;
   
    //survey catchability parameters
    if (debug) cout<<tb<<"surveys=list("<<endl;
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnQ,    pLnQ,   debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQT,  pLnDQT, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQX,  pLnDQX, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQM,  pLnDQM, debug,cout); if (debug){cout<<cc<<endl;}
    if (debug) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQXM, pLnDQXM,debug,cout); if (debug){cout<<endl;}
    if (debug) cout<<tb<<")"<<endl;
    
    if (debug<0) cout<<")"<<endl;
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;
}

void model_parameters::ReportToR_Data(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_Data(...)"<<endl;
    ptrMDS->writeToR(os,"data",0);
    if (debug) cout<<"Finished ReportToR_Data(...)"<<endl;
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
    //abundance
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    d5_array n_xmsyz = tcsam::rearrangeYXMSZtoXMSYZ(vn_yxmsz);
    //natural mortality (numbers)
    d5_array vnmN_yxmsz = wts::value(nmN_yxmsz);
    d5_array nmN_xmsyz = tcsam::rearrangeYXMSZtoXMSYZ(vnmN_yxmsz);
    //total mortality (numbers)
    d5_array vtmN_yxmsz = wts::value(tmN_yxmsz);
    d5_array tmN_xmsyz = tcsam::rearrangeYXMSZtoXMSYZ(vtmN_yxmsz);
        
    //total fishing mortality rates
    d5_array vtmF_yxmsz = wts::value(tmF_yxmsz);
    d5_array tmF_xmsyz  = tcsam::rearrangeYXMSZtoXMSYZ(vtmF_yxmsz);
    if (debug) cout<<"finished tmF_xmsyz"<<endl;
    
    os<<"pop.quants=list("<<endl;
    os<<"R.y=";       wts::writeToR(os,value(R_y),               ptrMC->dimYrsToR);os<<cc<<endl;
    os<<"Rx.y=";      wts::writeToR(os,trans(value(R_yx))(MALE), ptrMC->dimYrsToR);os<<cc<<endl;
    os<<"Rx.c=";      wts::writeToR(os,value(Rx_c),adstring("pc=1:"+str(npcRec))); os<<cc<<endl;
    os<<"R.cz=";      wts::writeToR(os,value(R_cz),adstring("pc=1:"+str(npcRec)),    ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"M.cxm=";     wts::writeToR(os,value(M_cxm),   adstring("pc=1:"+str(npcNM)), ptrMC->dimSXsToR,ptrMC->dimMSsToR); os<<cc<<endl;
    os<<"prMat.cz=";  wts::writeToR(os,value(prMat_cz),adstring("pc=1:"+str(npcMat)),ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"prGr.czz=";  wts::writeToR(os,value(prGr_czz),adstring("pc=1:"+str(npcGr)), ptrMC->dimZBsToR,ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"spb.yx=";    wts::writeToR(os,value(spb_yx),ptrMC->dimYrsToR,ptrMC->dimSXsToR); os<<cc<<endl;
    os<<"n.xmsyz=";   wts::writeToR(os,n_xmsyz,  ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsP1ToR,ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"nmN.xmsyz="; wts::writeToR(os,nmN_xmsyz,ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"tmN.xmsyz="; wts::writeToR(os,tmN_xmsyz,ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"tmF.xmsyz="; wts::writeToR(os,tmF_xmsyz,ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_PopQuants(...)"<<endl;
}

void model_parameters::ReportToR_SelFuncs(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_SelFuncs(...)"<<endl;
    os<<"sel.funcs=list("<<endl;
    os<<"sel_cz=";  wts::writeToR(os,value(sel_cz),  adstring("pc=1:"+str(npcSel)),ptrMC->dimZBsToR);os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_SelFuncs(...)"<<endl;
}

void model_parameters::ReportToR_FshQuants(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_FshQuants(...)"<<endl;
    //fishing capture rates (NOT mortality rates)
    d5_array vcF_fyxms = wts::value(cF_fyxms);//fully-selected
    d5_array cF_fxmsy  = tcsam::rearrangeIYXMStoIXMSY(vcF_fyxms);
    if (debug) cout<<"finished cF_fxmsy"<<endl;
    d6_array vFc_fyxmsz = wts::value(cF_fyxmsz);
    d6_array cF_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vFc_fyxmsz);
    if (debug) cout<<"finished cF_fxmsyz"<<endl;
    
    //retention mortality rates
    d6_array vrmF_fyxmsz = wts::value(rmF_fyxmsz);
    d6_array rmF_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vrmF_fyxmsz);
    if (debug) cout<<"finished rmF_fxmsyz"<<endl;
    
    //discard mortality rates
    d6_array vdmF_fyxmsz = wts::value(dmF_fyxmsz);
    d6_array dmF_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vdmF_fyxmsz);
    if (debug) cout<<"finished dmF_fxmsyz"<<endl;
    
    //numbers, biomass captured (NOT mortality)
    d6_array vcN_fyxmsz = wts::value(cN_fyxmsz);
    d3_array cN_fxy     = tcsam::calcIXYfromIYXMSZ(vcN_fyxmsz);
    d3_array cB_fxy     = tcsam::calcIXYfromIYXMSZ(vcN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array cN_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vcN_fyxmsz);
    if (debug) cout<<"finished cN_fxmsyz"<<endl;
    
    //numbers, biomass discards (NOT mortality)
    d6_array vdN_fyxmsz = wts::value(dN_fyxmsz);
    d3_array dN_fxy     = tcsam::calcIXYfromIYXMSZ(vdN_fyxmsz);
    d3_array dB_fxy     = tcsam::calcIXYfromIYXMSZ(vdN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array dN_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vdN_fyxmsz);
    if (debug) cout<<"finished dN_fxmsyz"<<endl;
    
    //numbers, biomass retained (mortality))
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    d3_array rmN_fxy     = tcsam::calcIXYfromIYXMSZ(vrmN_fyxmsz);
    d3_array rmB_fxy     = tcsam::calcIXYfromIYXMSZ(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array rmN_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vrmN_fyxmsz);
    if (debug) cout<<"finished rmN_fxmsyz"<<endl;
    
    //numbers, biomass discard mortality
    d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
    d3_array dmN_fxy     = tcsam::calcIXYfromIYXMSZ(vdmN_fyxmsz);
    d3_array dmB_fxy     = tcsam::calcIXYfromIYXMSZ(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d6_array dmN_fxmsyz  = tcsam::rearrangeIYXMSZtoIXMSYZ(vdmN_fyxmsz);
    if (debug) cout<<"finished dmN_fxmsyz"<<endl;   
    
    os<<"fisheries=list("<<endl;
    for (int f=1;f<=nFsh;f++){
        os<<ptrMC->lblsFsh[f]<<"=list("<<endl;
        os<<"cap=list("<<endl;
            os<<"n.xy=";   wts::writeToR(os,cN_fxy(f),   ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
            os<<"b.xy=";   wts::writeToR(os,cB_fxy(f),   ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
            os<<"n.xmsyz=";wts::writeToR(os,cN_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<cc<<endl;
            os<<"F.xmsy="; wts::writeToR(os,cF_fxmsy(f), ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR); os<<cc<<endl;
            os<<"F.xmsyz=";wts::writeToR(os,cF_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); 
        os<<")"<<cc<<endl;
        os<<"dm=list("<<endl;
            os<<"n.xy="; wts::writeToR(os,dmN_fxy(f),      ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
            os<<"b.xy="; wts::writeToR(os,dmB_fxy(f),      ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
            os<<"n.xmsyz="; wts::writeToR(os,dmN_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<cc<<endl;
            os<<"F.xmsyz="; wts::writeToR(os,dmF_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); 
        os<<")"<<cc<<endl;
        if (sum(rmN_fxy(f))>0){
            os<<"rm=list("<<endl;
                os<<"n.xy="; wts::writeToR(os,rmN_fxy(f),      ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
                os<<"b.xy="; wts::writeToR(os,rmB_fxy(f),      ptrMC->dimSXsToR,ptrMC->dimYrsToR); os<<cc<<endl;
                os<<"n.xmsyz="; wts::writeToR(os,rmN_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); os<<cc<<endl;
                os<<"F.xmsyz="; wts::writeToR(os,rmF_fxmsyz(f),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsToR,ptrMC->dimZBsToR); 
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
        os<<"q.xmsy=";  wts::writeToR(os,q_vxmsy(v), ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsP1ToR); os<<cc<<endl;
        os<<"q.xmsyz="; wts::writeToR(os,q_vxmsyz(v),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsP1ToR,ptrMC->dimZBsToR); os<<cc<<endl;
        os<<"n.xy=";    wts::writeToR(os,n_vxy(v),  ptrMC->dimSXsToR,ptrMC->dimYrsP1ToR); os<<cc<<endl;
        os<<"b.xy=";    wts::writeToR(os,b_vxy(v),  ptrMC->dimSXsToR,ptrMC->dimYrsP1ToR); os<<cc<<endl;
        os<<"spb.yx=";  wts::writeToR(os,value(spb_vyx(v)),ptrMC->dimYrsP1ToR,ptrMC->dimSXsToR); os<<cc<<endl;
        os<<"n.xmsyz="; wts::writeToR(os,n_vxmsyz(v),ptrMC->dimSXsToR,ptrMC->dimMSsToR,ptrMC->dimSCsToR,ptrMC->dimYrsP1ToR,ptrMC->dimZBsToR); os<<endl;
        os<<")"<<cc<<endl;
    }
    os<<"NULL)";
    if (debug) cout<<"Finished ReportToR_SrvQuants(...)"<<endl;
}

void model_parameters::ReportToR_ModelFits(ostream& os, int debug, ostream& cout)
{
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    //recalc objective function components and and write results to os
    os<<"model.fits=list("<<endl;
        os<<tb<<"penalties="; calcPenalties(-1,os);      os<<cc<<endl;
        os<<tb<<"priors=";    calcAllPriors(-1,os);      os<<cc<<endl;
        os<<tb<<"components=list("<<endl;
            os<<tb<<tb<<"recruitment="; calcNLLs_Recruitment(-1,os); os<<endl;
        os<<tb<<")"<<cc<<endl;
        os<<tb<<"fisheries="; calcNLLs_Fisheries(-1,os);  os<<cc<<endl; 
        os<<tb<<"surveys=";   calcNLLs_Surveys(-1,os);    os<<endl;  
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
        createSimData(debug, cout, 0, ptrSimMDS);//deterministic
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
        writeSimData(echo1,0,cout,ptrSimMDS);
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
