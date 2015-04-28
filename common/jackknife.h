/***************************************************************************
 *   Copyright (C) 2009, 2010 by Florian Goth   *
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
#ifndef JACKKNIFE_H
#include <cmath>
#include <omp.h>
#include "errordata.h"
#include "analysis_functions.h"
namespace mc_analysis
{
/**
this function calculates the expectation value of the samples stored in the container cont
@param cont an array of length len that stores the samples
@param len the length of the array
@return the average value of the stored samples
*/
template<typename T, typename IntType>
T expectationvalue(const T *const cont, IntType len)
{
    T retval(cont[0]);
    for (IntType k = 1; k < len; ++k)
        retval += cont[k];
    return retval/static_cast<typename ToScalar<T>::RetType>(len);
}

/**
This is a trait that maps between pointer and object like behaviour
*/
template <class Cont>
struct Access_Trait
{
    static inline const Cont& deref(const Cont& a)
    {
        return a;
    }
    typedef typename Cont::value_type ContainerElementType;
};

template <class Cont>
struct Access_Trait<Cont*>
{
    static inline const Cont& deref(const Cont *const& a)
    {
        const Cont& retval(*a);
        return retval;
    }
    typedef typename Cont::value_type ContainerElementType;
};

/**
This function calculates the function func from the given data and performs a jackknife analysis afterwards
The template parameters are as follows:
Func   the Function-Object that represents the Function we want to calculate from the data
Cont   the Container that is used for the Storage of the data
@param func the function we want to evaluate
@param data an array of Containers with the stored samples
@param numcont the number of arrays pointed to by data
*/
template <typename Func, typename Cont>
errordata<typename Func::res_t > jackknife(Func func, const Cont *const data, unsigned int numcont)
{
    typedef typename Access_Trait<Cont>::ContainerElementType ElementType;
    typedef typename Func::res_t RetType;
    unsigned int len = Access_Trait<Cont>::deref(data[0]).size();//data[0].size();//we assume all datasets have the same length
    ElementType mean[numcont];
    for (unsigned int n = 0; n < numcont; ++n)
        mean[n] = mc_analysis::mean(Access_Trait<Cont>::deref(data[n]));
    //mean now contains the averages of each dataset
    ElementType x_J[numcont];
    RetType* jackbin = new RetType[len];
    //calculate jackknife samples of the function
    for (unsigned int k = 0; k < len; ++k)
    {
        for (unsigned int j = 0; j < numcont; ++j)
        {
            x_J[j] = ElementType(1.0/len) *(ElementType(len) * mean[j] - Access_Trait<Cont>::deref(data[j])[k]);
        }
        jackbin[k] = func(x_J);
    }

    //calculate the jackknife average of the jackknife samples
    RetType jackmean = expectationvalue(jackbin, len);
    RetType fun_mean = func(mean);
    RetType fun_error(0);

    //calculate the error
    for (unsigned int k = 0; k < len; ++k)
    {
	RetType temp = jackmean - jackbin[k];
        fun_error += temp * temp;
    }
    delete [] jackbin;
    fun_error *= static_cast<double>(len - 1)/static_cast<double>(len);
    fun_error = std::sqrt(fun_error);
    return errordata<RetType>(fun_mean, fun_error, (fun_mean - jackmean) * RetType(len - 1));
}

/**
This function calculates the function func from the given data and performs a jackknife analysis afterwards
The template parameters are as follows:
Func   the Function-Object that represents the Function we want to calculate from the data
Cont   the Container that is used for the Storage of the data
@param func the function we want to evaluate
@param data an array of Containers with the stored samples
@param numcont the number of arrays pointed to by data
*/

template <typename Func, typename Cont>
errordata<typename Func::res_t > vecjackknife(Func func, const Cont& data, unsigned int numcont, bool covariance)
{
    typedef typename Cont::value_type::value_type ElementType;
    typedef typename Func::res_t RetType;
    auto len = data[0].size();//we assume all datasets have the same length. this is the number of functions not the number of bins!!!
    auto mean = mc_analysis::mean(data);
    //mean now contains the averages of each dataset
    auto x_J = mean;
    RetType* jackbin = new RetType[data.size()];
    auto norm = ElementType(1.0/data.size());
    //calculate jackknife samples of the function
    for (unsigned int k = 0; k < data.size(); ++k)
    {
        x_J = norm *(ElementType(data.size()) * mean - data[k]);
        jackbin[k] = func(x_J);
    }
    //calculate the jackknife average of the jackknife samples
    RetType jackmean = expectationvalue(jackbin, data.size());
    RetType fun_mean = func(mean);
    RetType fun_error = (jackmean - jackbin[0]) * (jackmean - jackbin[0]);
    //calculate the error
    for (unsigned int k = 1; k < data.size(); ++k)
    {
	RetType temp = jackmean - jackbin[k];
        fun_error += temp * temp;
    }
    typedef typename ToScalar<RetType>::RetType ScalarType;
    auto fac = typename ToScalar<RetType>::RetType(data.size() - 1);
    fun_error *= fac/static_cast<double>(data.size());
    fun_error = std::sqrt(fun_error);
    errordata<RetType> retval(fun_mean, fun_error, (fun_mean - jackmean) * fac);
    //let's do the covariance matrix
    if(covariance)
    {
#ifdef _OPENMP
    double start = omp_get_wtime();
#endif
    std::valarray<ScalarType> cov(numcont*numcont);
#pragma omp parallel for schedule(dynamic)
    for(uint y = 0; y < numcont; ++y)
    {
      ScalarType fmy = fun_mean[y];
      for(uint x = 0; x < numcont; ++x)
      {
	ScalarType fmx = fun_mean[x];
	ScalarType coventry(0);
	for(uint k = 0; k < data.size(); ++k)
	{
	   coventry += (jackbin[k][x] - fmx) * (jackbin[k][y] - fmy);
	}
	coventry *= (data.size() - 1.0)/data.size();
	cov[y * numcont + x] = coventry;
      }
    }
#ifdef _OPENMP
    std::cout<<"calculation of covariance took "<<omp_get_wtime() - start<<" s."<<std::endl;
#endif
    retval.setCov(cov);
    }
    delete [] jackbin;
    return retval;
}
}
#endif
