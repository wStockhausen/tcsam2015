#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "SummaryFunctions.hpp"

/**
 * Calculate marginal sums over msz by ixy.
 * @param n_iyxmsz - array to calculate sums on
 * @return d3_array with dims ixy.
 */
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

/**
 * Calculate marginal sums over msz by ixy, weighted by w_xmz.
 * @param n_iyxmsz - array to calculate sums on
 * @param w_xmz - weight-at-xmz array
 * @return d3_array with dims ixy.
 */
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

/**
 * Calculate marginal sums over ms by ixyz.
 * @param n_iyxmsz - array to calculate sums on
 * @return d4_array with dims ixyz.
 */
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

/**
 * Calculate marginal sums over s by ixmyz.
 * @param n_iyxmsz - array to calculate sums on
 * @return d5_array with dims ixmyz.
 */
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

/**
 * Sum iyxmsz array over maturity states and rearrange indices to ixsyz.
 * 
 * @param n_iyxmsz - d6_array w/ indices iyxmsz to sum
 * 
 * @return d5_array with indices ixsyz
 */
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
//    cout<<"In rearrangeYXMSZtoXMSYZ(...)"<<endl;
    ivector perm(1,5);//{4,1,2,3,5};
    perm[1]=4;perm[2]=1;perm[3]=2;perm[4]=3;perm[5]=5;
    d5_array n_xmsyz = wts::permuteDims(perm,n_yxmsz);
//    cout<<"Finished rearrangeYXMSZtoXMSYZ(...)"<<endl;
    return n_xmsyz;
}

/**
 * Rearrange indices for a d5_array from IYXMS to IXMSY.
 * 
 * @param n_iyxms - d5_array with indices iyxms
 * 
 * @return d5_array with indices ixmsy
 */
d5_array tcsam::rearrangeIYXMStoIXMSY(d5_array& n_iyxms){
//    cout<<"starting rearrangeIYXMStoIXMSY(...)"<<endl;
    ivector perm(1,5);//{1,5,2,3,4};
    perm[1]=1;perm[2]=5;perm[3]=2;perm[4]=3;perm[5]=4;
    d5_array n_ixmsy = wts::permuteDims(perm,n_iyxms);
//    cout<<"finished rearrangeIYXMStoIXMSY(...)"<<endl;
    return n_ixmsy;
}

/**
 * Rearrange indices for a d6_array from IYXMSZ to IXMSYZ.
 * 
 * @param n_iyxmsz - d6_array with indices iyxmsz
 * 
 * @return d6_array with indices ixmsyz
 */
d6_array tcsam::rearrangeIYXMSZtoIXMSYZ(d6_array& n_iyxmsz){
//    cout<<"starting rearrangeIYXMSZtoIXMSYZ(...)"<<endl;
    ivector perm(1,6);//{1,5,2,3,4,6};
    perm[1]=1;perm[2]=5;perm[3]=2;perm[4]=3;perm[5]=4;perm[6]=6;
    d6_array n_ixmsyz = wts::permuteDims(perm,n_iyxmsz);
//    cout<<"finished rearrangeIYXMSZtoIXMSYZ(...)"<<endl;
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
