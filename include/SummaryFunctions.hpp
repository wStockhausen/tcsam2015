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
}
#endif	/* SUMMARYFUNCTIONS_HPP */

