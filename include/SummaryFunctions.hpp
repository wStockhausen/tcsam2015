/* 
 * File:   SummaryFunctions.hpp
 * Author: WilliamStockhausen
 *
 * Created on April 4, 2014, 3:39 PM
 */

#ifndef SUMMARYFUNCTIONS_HPP
#define	SUMMARYFUNCTIONS_HPP

namespace tcsam {
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param  n_ij - dmatrix
     * @param  w    - weighting vector
     * @return n_i  - dvector
     */
    dvector sumOverLastDim(const dmatrix& n_ij, dvector& w);
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param  n_ijk - d3_array
     * @param  w    - weighting vector
     * @return n_ij  - dmatrix
     */
    dmatrix sumOverLastDim(const d3_array& n_ijk, dvector& w);
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param n_ijkl - d4_array
     * @param  w    - weighting vector
     * @return n_ijk - d3_array
     */
    d3_array sumOverLastDim(const d4_array& n_ijkl, dvector& w);
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param  n_ijklm - d5_array
     * @param  w    - weighting vector
     * @return n_ijkl  - d4_array
     */
    d4_array sumOverLastDim(const d5_array& n_ijklm, dvector& w);
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param  n_ijklmn - d6_array
     * @param  w    - weighting vector
     * @return n_ijklm  - d4_array
     */
    d5_array sumOverLastDim(const d6_array& n_ijklmn, dvector& w);
    
    /**
     * Compute marginal weighted sums over last dimension.
     * 
     * @param  n_ijklmno - d7_array
     * @param  w         - weighting vector
     * @return n_ijklmn  - d6_array
     */
    d6_array sumOverLastDim(const d7_array& n_ijklmno, dvector& w);
    
    /**
     * Calculate the sub-array n_iyx from full array n_iyxmsz by
     * summing over indices m, s and z.
     * @param n_iyxmsz
     * @return d3_array of resulting sums
     */
    d3_array calcIYXfromIYXMSZ(d6_array& n_iyxmsz);
    
    /**
     * Calculate the biomass-weighted sub-array b_ixy 
     * from full abundance array n_iyxmsz by summing over 
     * indices m, s and z weighted by weight-at-sex/maturity/size.
     * @param n_iyxmsz - full abundance array (1000's)
     * @param w_xmz    - weight (kg) at sex/maturity/size.
     * @return - biomass for index i, year, sex in mt.
     */
    d3_array calcIYXfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz);
    
    /**
     * Calculate the biomass-weighted sub-array b_yxms 
     * from full abundance array n_yxmsz by summing over 
     * indices m, s and z weighted by weight-at-sex/maturity/size.
     * @param n_yxmsz - full abundance array (1000's)
     * @param w_xmz   - weight (kg) at sex/maturity/size.
     * @return - biomass for index i, year, sex in mt.
     */
    d4_array calcYXMSfromYXMSZ(d5_array& n_yxmsz, d3_array w_xmz);
    
    /**
     * Calculate the biomass-weighted sub-array b_iyxms 
     * from full abundance array n_iyxmsz by summing over 
     * indices m, s and z weighted by weight-at-sex/maturity/size.
     * @param n_iyxmsz - full abundance array (1000's)
     * @param w_xmz    - weight (kg) at sex/maturity/size.
     * @return - biomass for index i, year, sex in mt.
     */
    d5_array calcIYXMSfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz);
    
    /**
     * Calculate the sub-array n_ixy from full array n_iyxmsz by
     * summing over indices m, s and z.
     * @param n_iyxmsz
     * @return d3_array of resulting sums
     */
    d3_array calcIXYfromIYXMSZ(d6_array& n_iyxmsz);
    
    /**
     * Calculate the biomass-weighted sub-array b_ixy 
     * from full abundance array n_iyxmsz by summing over 
     * indices m, s and z weighted by weight-at-sex/maturity/size.
     * @param n_iyxmsz - full abundance array (1000's)
     * @param w_xmz    - weight (kg) at sex/maturity/size.
     * @return - biomass for index i, sex, year in mt.
     */
    d3_array calcIXYfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz);
    
    /**
     * Calculate the sub-array n_ixyz from full array n_iyxmsz by
     * summing over indices m and s.
     * @param n_iyxmsz
     * @return 
     */
    d4_array calcIXYZfromIYXMSZ(d6_array& n_iyxmsz);
        
    /**
     * Calculate the sub-array n_ixmyz from full array n_iyxmsz by
     * summing over index s.
     * @param n_iyxmsz
     * @return 
     */
    d5_array calcIXMYZfromIYXMSZ(d6_array& n_iyxmsz);
    
    /**
     * Calculate the sub-array n_ixsyz from full array n_iyxmsz by
     * summing over index m.
     * @param n_iyxmsz
     * @return 
     */
    d5_array calcIXSYZfromIYXMSZ(d6_array& n_iyxmsz);
    
    /**
     * Rearrange the input array n_yxmsz to n_xmsyz.
     * 
     * @param d5_array n_yxmsz
     * @return d5_array n_xmsyz
     */
    d5_array rearrangeYXMSZtoXMSYZ(d5_array& n_yxmsz);
    
    /**
     * Rearrange the input array n_iyxms to n_ixmsy.
     * 
     * @param d5_array n_iyxms
     * @return d5_array n_ixmsy
     */
    d5_array rearrangeIYXMStoIXMSY(d5_array& n_iyxms);
    
    /**
     * Rearrange the input array n_iyxmsz to n_ixmsyz.
     * 
     * @param d6_array n_iyxmsz
     * @return d6_array n_ixmsyz
     */
    d6_array rearrangeIYXMSZtoIXMSYZ(d6_array& n_iyxmsz);
    
    /**
     * Extract (possibly summary) vector from 5d array w/ indices
     * x,m,s,y,z.
     * 
    * @param x
    * @param m
    * @param s
    * @param y
    * @param n_xmsyz
     * @return extracted vector (indices consistent with z)
     */
    dvector extractFromXMSYZ(int x, int m, int s, int y, d5_array& n_xmsyz);
    
    /**
     * Extract (possibly summary) value from 4d array w/ indices
     * y,x,m,s.
     * 
    * @param y
    * @param x
    * @param m
    * @param s
    * @param n_yxms
     * @return extracted value
     */
    double extractFromYXMS(int y, int x, int m, int s, d4_array& n_yxms);
    
    /**
     * Extract (possibly summary) vector from 5d array w/ indices
     * y,x,m,s,z.
     * 
    * @param y
    * @param x
    * @param m
    * @param s
    * @param n_yxmsz
     * @return extracted vector (indices consistent with z)
     */
    dvector extractFromYXMSZ(int y, int x, int m, int s, d5_array& n_yxmsz);
}
#endif	/* SUMMARYFUNCTIONS_HPP */

