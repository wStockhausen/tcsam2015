#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "SummaryFunctions.hpp"

d3_array tcsam::calcIXYfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array n_ixy(mni,mxi,mnx,mxx,mny,mxy);
    n_ixy.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                n_ixy(i,x,y) = sum(n_iyxmsz(i,y,x));
            }
        }
    }
    return n_ixy;
}

d3_array tcsam::calcIXYfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array b_ixy(mni,mxi,mnx,mxx,mny,mxy);
    b_ixy.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) b_ixy(i,x,y) += n_iyxmsz(i,y,x,m,s)*w_xmz(x,m);//dot product here
                }
            }
        }
    }
    return b_ixy;
}

d4_array tcsam::calcIXYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d4_array n_ixyz(mni,mxi,mnx,mxx,mny,mxy,mnz,mxz);
    n_ixyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixyz(i,x,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixyz;
}

d5_array tcsam::calcIXMYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d5_array n_ixmyz(mni,mxi,mnx,mxx,mnm,mxm,mny,mxy,mnz,mxz);
    n_ixmyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixmyz(i,x,m,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixmyz;
}

d5_array tcsam::calcIXSYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d5_array n_ixsyz(mni,mxi,mnx,mxx,mns,mxs,mny,mxy,mnz,mxz);
    n_ixsyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixsyz(i,x,s,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixsyz;
}

d5_array tcsam::rearrangeYXMSZtoXMSYZ(d5_array& n_yxmsz){
    ivector bnds = wts::getBounds(n_yxmsz);
    int mny = bnds( 1); int mxy = bnds( 2);
    int mnx = bnds( 3); int mxx = bnds( 4);
    int mnm = bnds( 5); int mxm = bnds( 6);
    int mns = bnds( 7); int mxs = bnds( 8);
    int mnz = bnds( 9); int mxz = bnds(10);
    d5_array n_xmsyz(mnx,mxx,mnm,mxm,mns,mxs,mny,mxy,mnz,mxz);
    n_xmsyz.initialize();
    for (int x=mnx;x<=mxx;x++){
        for (int y=mny;y<=mxy;y++){
            for (int m=mnm;m<=mxm;m++) {
                for (int s=mns;s<=mxs;s++) n_xmsyz(x,m,s,y) = n_yxmsz(y,x,m,s);
            }
        }
    }
    return n_xmsyz;
}

d5_array tcsam::rearrangeIYXMStoIXMSY(d5_array& n_iyxms){
    ivector bnds = wts::getBounds(n_iyxms);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    d5_array n_ixmsy(mni,mxi,mnx,mxx,mnm,mxm,mns,mxs,mny,mxy);
    n_ixmsy.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++) {
                    for (int s=mns;s<=mxs;s++) n_ixmsy(i,x,m,s,y) = n_iyxms(i,y,x,m,s);
                }
            }
        }
    }
    return n_ixmsy;
}

d6_array tcsam::rearrangeIYXMSZtoIXMSYZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d6_array n_ixmsyz(mni,mxi,mnx,mxx,mnm,mxm,mns,mxs,mny,mxy,mnz,mxz);
    n_ixmsyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++) {
                    for (int s=mns;s<=mxs;s++) n_ixmsyz(i,x,m,s,y) = n_iyxmsz(i,y,x,m,s);
                }
            }
        }
    }
    return n_ixmsyz;
}

/**
 * Extract (possibly summary) vector from 5d array.
 * 
 * @param x - sex index             (tcsam::ALL_SXs yields sum over sex)
 * @param m - maturity index        (tcsam::ALL_MSs yields sum over maturity state)
 * @param s - shell condition index (tcsam::ALL_SCs yields sum over shell condition)
 * @param y - year index
 * @param n_xmsyz - d5_array from which to extract vector at z's
 * 
 * @return extracted vector (indices consistent with z's)
 */
dvector tcsam::extractFromXMSYZ(int x, int m, int s, int y, d5_array& n_xmsyz){
    ivector bnds = wts::getBounds(n_xmsyz);
    dvector n_z(bnds(9,10));//dimension for z index
    int xmn, xmx;
    xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::nSXs;}
    int mmn, mmx;
    mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::nMSs;}
    int smn, smx;
    smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::nSCs;}
    n_z.initialize();
    for (int xp=xmn;xp<=xmx;xp++){
        for (int mp=mmn;mp<=mmx;mp++){
            for (int sp=smn;sp<=smx;sp++){
                n_z += n_xmsyz(xp,mp,sp,y);
            }
        }
    }
    return(n_z);
}

/**
 * Extract (possibly summary) vector from 5d array w/ indices
 * y,x,m,s,z.
 * 
 * @param y - year index
 * @param x - sex index             (tcsam::ALL_SXs yields sum over sex)
 * @param m - maturity index        (tcsam::ALL_MSs yields sum over maturity state)
 * @param s - shell condition index (tcsam::ALL_SCs yields sum over shell condition)
 * @param n_yxmsz - d5_array from which to extract vector at z's
 * 
 * @return extracted vector (indices consistent with z's)
 */
dvector tcsam::extractFromYXMSZ(int y, int x, int m, int s, d5_array& n_yxmsz){
    ivector bnds = wts::getBounds(n_yxmsz);
//    cout<<"in extractFromYXMSZ"<<endl;
    dvector n_z(bnds(9),bnds(10));//dimension for z index
//    cout<<n_z.indexmin()<<" "<<n_z.indexmax()<<endl;
//    cout<<y<<" "<<x<<" "<<m<<" "<<" "<<s<<endl;
    int xmn, xmx;
    xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::nSXs;}
    int mmn, mmx;
    mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::nMSs;}
    int smn, smx;
    smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::nSCs;}
    n_z.initialize();
    for (int xp=xmn;xp<=xmx;xp++){
        for (int mp=mmn;mp<=mmx;mp++){
            for (int sp=smn;sp<=smx;sp++){
//                cout<<xp<<" "<<mp<<" "<<sp<<endl;
                n_z += n_yxmsz(y,xp,mp,sp);
            }
        }
    }
//    cout<<"finished extractFromYXMSZ"<<endl;
    return(n_z);
}
