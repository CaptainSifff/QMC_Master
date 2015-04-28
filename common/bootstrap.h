/***************************************************************************
 *   Copyright (C) 2009 by Florian Goth   *
 *   fgoth@wthp095   *
 *                                                                         *
 *   All rights reserved.                                                  *
 *                                                                         *
 *   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: *
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. *
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. *
 *     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. *
 *                                                                         *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR *
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR *
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, *
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,   *
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF *
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          *
 ***************************************************************************/
#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H
#include <iostream>
#include <cmath>//for sqrt
#include "errordata.h"

using namespace std;
namespace mc_analysis
{
/**
This functions calculates the statistical measure given by Func via a bootstrap estimate.
the template parameters are as follows:
Func   the Function-Object that represents the statistical measure we want to calculate from the data
Cont   the Container that is used for the Storage of the data
@param func the function we want to evaluate
@param data an array of Containers with the stored samples
@param numcont the number of arrays pointed to by data
*/
template <class Func, class Cont>
mc_analysis::errordata<typename Func::RetType> bootstrap(Func func, const Cont *const cont, unsigned int numcont)
{
    class Ziff
    {
    public:
        Ziff(unsigned int seed) : Random_Max(2147483648U), Random_nd(0)//forgotten in Hinrichsens Code ?
        {
            // we initialize via a "good" LCG random number generator
            class LCG
            {
                //A Linear Congruential Generator
            public:
                inline unsigned int rnd()
                {
                    return (x =(62089911 * x + 4349) % 2147483647u);// the familiar scheme for generating random numbers with a linear congruential one
                }
                LCG(unsigned int seed) : x(seed) {}
            private:
                unsigned int x;///< the value that is stored in between the generation of a random number
            } lcg(seed);
            for (unsigned int k = 0; k <= Random_M; ++k)//some warmup
                Random_ra[k] = lcg.rnd();
            //some warmup
            for (unsigned int k = 0; k < 100; ++k) rnd();
            return;
        }
        inline unsigned int rnd()
        {
            ++Random_nd;
            return (Random_ra[Random_nd & Random_M] =
                        Random_ra[(Random_nd - Random_A) & Random_M] ^
                        Random_ra[(Random_nd - Random_B) & Random_M] ^
                        Random_ra[(Random_nd - Random_C) & Random_M] ^
                        Random_ra[(Random_nd - Random_D) & Random_M]);
        }
    private:
        enum  //various constants. we store them in an enum
        {
            Random_A = 471,
            Random_B = 1586,
            Random_C = 6988,
            Random_D = 9689,
            Random_M = 16383,
        };
        const unsigned int Random_Max;///< the Maximum number that can be generated
        int Random_nd;
        int Random_ra[Random_M+1];///<the array of numbers the RNG is working on
    } prng(141799041);//from random.org. guaranteed to be completely random...
    typedef typename Func::RetType RetType;
    typedef typename Cont::value_type ElementType;
    const unsigned int len = cont[0].size();
    const unsigned int samplecnt = 500;//the number of bootstrapsamples that we generate for estimating the error
    unsigned int mempos = 0;
    ElementType** bootstrapsample;
    try
    {
        bootstrapsample = new ElementType*[len];//this contains the memory that is used for one bootstrapsample
        for (mempos = 0; mempos < len; ++mempos)
            bootstrapsample[mempos] = new ElementType[numcont];
    }
    catch (bad_alloc& ba)
    {
        if (mempos != 0)
        {
            for (unsigned int i = 0; i < mempos; ++i)
                delete [] bootstrapsample[i];
            delete [] bootstrapsample;
        }
        cerr<<"[Bootstrap] FAILED TO ALLOCATE MEMORY FOR ONE BOOTSTRAP SAMPLE!!"<<endl;
        throw;
    }
    RetType* res;
    try
    {
        res = new RetType[samplecnt];
    }
    catch (bad_alloc& ba)
    {
        cerr<<"[Bootstrap] FAILED TO ALLOCATE MEMORY FOR BOOTSTRAP MEASUREMENTS!!"<<endl;
        throw;
    }
    RetType mean(0);
    for ( unsigned int k = 0; k < samplecnt; ++k)//do something with every sample
    {
        for (unsigned int i = 0; i < len; ++i)
        {
            unsigned int pos = prng.rnd() % len;
            for (unsigned int j = 0; j < numcont; ++j)
                bootstrapsample[i][j] = cont[j][pos];
        }
//bootstrapsample now contains one permutated bootstrapsample
        res[k] = func(bootstrapsample, len);
        mean += res[k];
    }
    for (unsigned int k = 0; k < len; ++k)
        delete [] bootstrapsample[k];
    delete [] bootstrapsample;
    mean /= samplecnt;
    //calculate the variance
    RetType var(0);
    for (unsigned int k = 0; k < samplecnt; ++k)
        var += (res[k] - mean) * (res[k] - mean);
    var /= samplecnt - 1.0;
    delete [] res;
    return mc_analysis::errordata<RetType>(mean, std::sqrt(var));
}

/**
A Mean Object for calculating means with the BootstrapTool
*/
template <typename T>
class Mean
{
public:
    T operator() (T** data, unsigned int len)
    {
        T RetVal = 0;
        for (unsigned int l = 0; l < len; ++l)
            RetVal += data[l][0];
        return RetVal/static_cast<T>(len);
    }
    typedef T RetType;
};
}

#endif
