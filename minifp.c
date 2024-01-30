
#include <stdio.h>
#include <stdlib.h>
#include "common_structs.h"
#include "common_definitions.h"
#include "common_functions.h"
#include "minifp.h"

#define NAN 0x7ff
#define POS_INFINITY 0x3c0
#define NEG_INFINITY 0x7c0
#define POS_ZERO 0x000
#define BIAS 7
#define DENORM_E -6


//helper function - determines if minifp is nan
int is_nan(minifp_s val){
    //checks that exp = 15 and frac > 0
    return (((val & (0xf << 6)) >> 6 == 0xf) && ((val & 0x3f) > 0));
}

//helper function - determines if minifp is infinity
int is_infinity(minifp_s val){
    
    return (val == POS_INFINITY || val == NEG_INFINITY);
}

//helper function - determines if minifp is denormalized
int is_denorm(minifp_s val){
    //checks that exp = 0
    return ((val & (0xf << 6)) >> 6 == 0);
}

//helper function - determines if minifp is negative
int is_negative(minifp_s val){
    //checks sign bit
    return (val & (1 << 10)) >> 10;
}

//helper function - determines if minifp is zero
int is_zero(minifp_s val){
    //checks that exp and frac are 0
    return ((val & 0x3ff) == 0);
}

/* helper function for arithmetic functions
*  extracts mantissas from minifps, matches E values, then shifts mantissas until both fracs are zero
*/
void float_to_base10(minifp_s val1, minifp_s val2, unsigned int * M1, unsigned int * M2, int * E){

    //determine values - denormalized or normalized
    int E1 = (is_denorm(val1)) ? DENORM_E : ((val1 & (0xf << 6)) >> 6) - BIAS;
    int E2 = (is_denorm(val2)) ? DENORM_E : ((val2 & (0xf << 6)) >> 6) - BIAS;
    unsigned int whole1 = (is_denorm(val1)) ? 0 : 1;
    unsigned int whole2 = (is_denorm(val2)) ? 0 : 1;
    unsigned int frac1 = (val1 & 0x3f) << 26;
    unsigned int frac2 = (val2 & 0x3f) << 26;

    //shift val1 so E1 and E2 match
    while(E1 != E2){
        if(E1 < E2){
            int lsb = (1 & whole1);
            whole1 >>= 1;
            frac1 >>= 1;
            frac1 += (lsb << 31);
            E1++;
        }
        else{
            int msb = ((1 << 31) & frac1) >> 31;
            frac1 <<= 1;
            whole1 <<= 1;
            whole1 += msb;
            E1--;
        }
    }

    //shift until both fracs are 0
    int result_E = E1;
    while(frac1 != 0 || frac2 != 0){
        //shift first mantissa
        int msb1 = ((1 << 31) & frac1) >> 31;
        frac1 <<= 1;
        whole1 <<= 1;
        whole1 += msb1;
        //shift second mantissa
        int msb2 = ((1 << 31) & frac2) >> 31;
        frac2 <<= 1;
        whole2 <<= 1;
        whole2 += msb2;
        //decrement E
        result_E--;
    }

    //update final values
    *M1 = whole1;
    *M2 = whole2;
    *E = result_E;
}

/* helper function for arithmetic functions
*  converts whole number mantissa, E, and S to minifp
*/
minifp_s base10_to_float(unsigned int M, int E, int S){

    //if result was zero
    if(M == 0) return (POS_ZERO ^ S);

    unsigned int whole = M;
    unsigned int frac = 0;
    int temp_E = E;

    //shift until mantissa has leading one
    while(whole > 1){
        int lsb = (1 & whole);
        whole >>= 1;
        frac >>= 1;
        frac += (lsb << 31);
        E++;
    }

    //determine exp
    int exp = (BIAS + E);
    //if result overflowed
    if(exp >= 15) return (POS_INFINITY ^ S);

    //if exp is <= 0, use denormalized process
    if(exp <= 0){
        //reset whole, fraction, E, and exp
        whole = M;
        frac = 0;
        exp = 0;
        E = temp_E;
        //shift mantissa until correct
        while(E != DENORM_E){
            int lsb = (1 & whole);
            frac >>= 1;
            whole >>= 1;
            frac += (lsb << 31);
            E++;
        }
        
    }

    //encode sign, frac, and exp
    minifp_s result = S + ((exp & 0xF) << 6) + (frac >> 26);

    return result;
}

/* toMiniFP - Converts a Number Struct (whole and fraction parts) into a MiniFP Value
 *  - number is managed by MUAN, DO NOT FREE number.
 *    - You may change the contents of number, but do not free it.
 *  - Follow the project documentation for this function.
 * Return the MiniFP Value or any legal MiniFP NaN representation if number is NULL.
 */
minifp_s toMiniFP(Number_s *number) {

    //if number is null, return NaN
    if(!number) return NAN;

    //if NaN flag is set, return NaN
    if(number->is_nan) return NAN;

    //if infinity flag is set or whole number is too large, return infinity
    if((number->whole >= 256 || number->is_infinity) && number->is_negative) return NEG_INFINITY;
    else if((number->whole >= 256 || number->is_infinity) && !number->is_negative) return POS_INFINITY;

    //if whole and fraction are zero, return zero
    if(!number->whole && !number->fraction) return 0;

    //set sign bit
    minifp_s minifp = (number->is_negative) ? 0x400 : 0x000;

    //create copies of values in case denormalization is needed
    unsigned int whole = number->whole;
    unsigned int fraction = number->fraction;
    int E = 0;

    //shift until mantissa is in range while adjusting E
    while(whole != 1){
        //if whole num is greater than 1
        if(whole > 1){
            //track lsb of whole before each shift and replace msb of fraction after shift
            int lsb = (1 & whole);
            whole >>= 1;
            fraction >>= 1;
            fraction += (lsb << 31);
            E++;
        }
        //if whole num is less than 1
        else{
            //track msb of fraction before each shift and replace lsb of whole after shift
            int msb = ((1 << 31) & fraction) >> 31;
            fraction <<= 1;
            whole <<= 1;
            whole += msb;
            E--;
        }
    }

    //determine exp
    int exp = (BIAS + E);

    //if exp is <= 0, use denormalized process
    if(exp <= 0){
        //reset whole, fraction, E, and exp
        whole = number->whole;
        fraction = number->fraction;
        E = 0;
        exp = 0;
        //shift mantissa until correct
        while(E != DENORM_E){
            int msb = ((1 << 31) & fraction) >> 31;
            fraction <<= 1;
            whole <<= 1;
            whole += msb;
            E--;
        }
        
    }

    //set exp bits
    minifp += (exp & 0xF) << 6;

    //set frac bits
    minifp += fraction >> 26;

    return minifp;
}

/* toNumber - Converts a MiniFP Value into a Number Struct (whole and fraction parts)
 *  - number is managed by MUAN, DO NOT FREE or re-Allocate number.
 *    - It is already allocated.  Do not call malloc/calloc for number.
 *  - Follow the project documentation for this function.
 *  If the conversion is successful, return 0. 
 *  - If number is NULL or there are any inconsistencies, return -1.
 */
int toNumber(Number_s *number, minifp_s value) {

    //if number is null
    if(!number) return -1;

    //if value is negative
    if(is_negative(value)){
        number->is_negative = 1;
    }

    //if value is infinity
    if(is_infinity(value)){
        number->is_infinity = 1;
        return 0;
    }

    //if value is NaN
    if(is_nan(value)){
        number->is_nan = 1;
        return 0;
    }

    int exp = (value & (0xF << 6)) >> 6;
    int E = 0;
    unsigned int whole = 0;
    unsigned int fraction = (value & 0x3F) << 26;

    //if exp is zero begin denorm process
    if(exp == 0){
        E = 1 - BIAS;
        whole = 0;
    }
    //otherwise use norm
    else{
        E = exp - BIAS;
        whole = 1;
    }

    //shift whole and fraction until E is zero
    while(E != 0){
        //shifting right
        if(E < 0){
            int lsb = (1 & whole);
            whole >>= 1;
            fraction >>= 1;
            fraction += (lsb << 31);
            E++;
        }
        //shifting left
        else{
            int msb = ((1 << 31) & fraction) >> 31;
            fraction <<= 1;
            whole <<= 1;
            whole += msb;
            E--;
        }
    }

    number->whole = whole;
    number->fraction = fraction;

    return 0;
}

/* mulMiniFP - Performs an operation on two miniFP values
 *  - Follow the project documentation for this function.
 * Return the resulting minifp_s value
 */
minifp_s mulMiniFP(minifp_s val1, minifp_s val2) {

    //SPECIAL CASES
    //determine sign by sign rules
    minifp_s sign = (is_negative(val1) ^ is_negative(val2)) ? 0x400 : 0x000;
    //NaN * anything = NaN
    if(is_nan(val1) || is_nan(val2)) return NAN;
    //infinity * 0 = NaN
    if((is_infinity(val1) && is_zero(val2)) || (is_infinity(val2) && is_zero(val1))) return NAN;
    //infinity * infinity = infinity
    if(is_infinity(val1) && is_infinity(val2)) return (POS_INFINITY ^ sign);
    //infinity * x = infinity
    if(is_infinity(val1) || is_infinity(val2)) return (POS_INFINITY ^ sign);
    //0 * x = 0
    if(is_zero(val1) || is_zero(val2)) return (POS_ZERO ^ sign);

    //REGULAR OPERATION
    int E = 0;
    unsigned int M1, M2 = 0;

    //shift mantissas to whole numbers and track current E
    //mantissas will be shifted to same E so E of product after multiplication will be E * 2
    float_to_base10(val1, val2, &M1, &M2, &E);
    E *= 2;

    //multiply shifted mantissas
    unsigned int product = M1 * M2;

    //create final minifp value using product, E, and sign bit
    minifp_s result = base10_to_float(product, E, sign);

    return result;
}

/* addMiniFP - Performs an operation on two miniFP values
 *  - Follow the project documentation for this function.
 * Return the resulting minifp_s value
 */
minifp_s addMiniFP(minifp_s val1, minifp_s val2) {

    //SPECIAL CASES
    //determine sign by sign rules
    minifp_s sign = (is_negative(val1) && is_negative(val2)) ? 0x400 : 0x000;
    //infinity - infinity = NaN
    if((val1 == POS_INFINITY && val2 == NEG_INFINITY) || (val1 == NEG_INFINITY && val2 == POS_INFINITY)) return NAN;
    //NaN + x = NaN
    if(is_nan(val1) || is_nan(val2)) return NAN;
    //infinity + infinity = infinity
    if(is_infinity(val1) && is_infinity(val2)) return (POS_INFINITY ^ sign);
    //x - x = 0
    if(val1 == negateMiniFP(val1)) return POS_ZERO;
    //infinity + x = infinity
    if(val1 == POS_INFINITY || val2 == POS_INFINITY) return POS_INFINITY;
    //-infinity + x = -infinity
    if(val1 == NEG_INFINITY || val2 == NEG_INFINITY) return NEG_INFINITY;
    //0 + 0 = 0
    if(is_zero(val1) && is_zero(val2)) return (POS_ZERO ^ sign);
    //0 + x = x
    if(is_zero(val1)) return val2;
    //x + 0 = x
    if(is_zero(val2)) return val1;

    //REGULAR OPERATION
    int E = 0;
    unsigned int M1, M2 = 0;
    
    //shift both mantissas until Es match, then shift until mantissas are whole numbers
    float_to_base10(val1, val2, &M1, &M2, &E);

    //determine addends and their parity
    int addend_1 = ((val1 & (1 << 10)) >> 10) ? -((int) M1) : (int) M1;
    int addend_2 = ((val2 & (1 << 10)) >> 10) ? -((int) M2) : (int) M2;
    int sum = addend_1 + addend_2;
    sign = (sum < 0) ? 0x400 : 0x000;

    //create final minifp value using sum, E, and sign bit
    minifp_s result = base10_to_float((unsigned int) (abs(sum)), E, sign);

    return result;
}

/* subMiniFP - Performs an operation on two miniFP values
 *  - Follow the project documentation for this function.
 * Return the resulting minifp_s value
 */
minifp_s subMiniFP(minifp_s val1, minifp_s val2) {

    return addMiniFP(val1, negateMiniFP(val2));
}

/* negateMiniFP - Negates a MiniFP Value.
 *  - Follow the project documentation for this function.
 * Return the resulting MiniFP Value
 */
minifp_s negateMiniFP(minifp_s value) {
  
    return (value ^ (1 << 10));
}
