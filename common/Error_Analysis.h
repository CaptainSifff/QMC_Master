/***************************************************************************
 *   Copyright (C) 2009-2012 by Florian Goth   *
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
#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H
#include <string>
#include <fstream>
#include <type_traits>
#include "AverageSign.h"
#include "errordata.h"
#include "jackknife.h"
#include "ClientState.h"
#include <complex>
#include <valarray>
#include <algorithm>

#include <unistd.h>
/*
//forward declare Observable
template <class ObsT, class SignType>
class Observable;*/

//let's carve in stone that the sign is always stored last...
template <typename T>
class ReFunction
{
public:
    inline T operator()(std::complex<T>* dat)
    {
        return (dat[0].real()*dat[len-1].real() + dat[0].imag()*dat[len-1].imag())/(dat[len-1].real()*dat[len-1].real() + dat[len-1].imag()*dat[len-1].imag());
    }
    typedef T res_t;
    ReFunction(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class ImFunction
{
public:
    inline T operator() (std::complex<T>* dat)
    {
        return (dat[0].imag() * dat[len-1].real() - dat[len-1].imag()*dat[0].real())/(dat[len-1].real()*dat[len-1].real() + dat[len-1].imag() * dat[len-1].imag());
    }
    typedef T res_t;
    ImFunction(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class ReVecFunction
{
public:
    inline std::valarray<T> operator()(std::valarray<std::complex<T> >& dat)
    {
      auto denom = dat[len-1].real()*dat[len-1].real() + dat[len-1].imag()*dat[len-1].imag();
      std::valarray<T> retval(len);
      for(uint k = 0; k < (len - 1); ++k)
      {
	retval[k] = (dat[k].real()*dat[len-1].real() + dat[k].imag()*dat[len-1].imag())/denom;
      }
      return retval;
    }
    typedef std::valarray<T> res_t;
    ReVecFunction(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class ImVecFunction
{
public:
    inline std::valarray<T> operator() (std::valarray<std::complex<T> >& dat)
    {
      auto denom = dat[len-1].real()*dat[len-1].real() + dat[len-1].imag()*dat[len-1].imag();
      std::valarray<T> retval(len);
      for(uint k = 0; k < (len - 1); ++k)
	retval[k] = (dat[k].imag() * dat[len-1].real() - dat[len-1].imag()*dat[k].real())/denom;
      return retval;
    }
    typedef std::valarray<T> res_t;
    ImVecFunction(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class PlainSign
{
public:
    inline T operator() (T* dat)
    {
        return dat[0]/dat[len-1];
    }
    typedef T res_t;
    PlainSign(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class PlainVecSign
{
public:
    inline T operator() (T* dat)
    {
      std::valarray<T> retval(len);
      for(uint k = 0; k < len-1; ++k)
	retval[k] = dat[k]/dat[len-1];
      return retval;
    }
    inline std::valarray<T> operator() (std::valarray<T>& dat)
    {
      std::valarray<T> retval(len);
      for(uint k = 0; k < (len-1); ++k)
	retval[k] = dat[k]/dat[len-1];
      return retval;
    }
    typedef std::valarray<T> res_t;
    PlainVecSign(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class PlainSign<std::complex<T> >
{//this is supposed to be used if we mix a complex observable and a real sign
public:
    inline std::complex<T> operator() (std::complex<T>* dat)
    {
        return dat[0]/real(dat[len-1]);
    }
    typedef std::complex<T> res_t;
    PlainSign< std::complex< T > >(uint l) : len(l) {}
private:
  uint len;
};

template <typename T>
class PlainVecSign<std::complex<T> >
{//this is supposed to be used if we mix a complex observable and a real sign
public:
    inline std::valarray<std::complex<T> > operator() (std::complex<T>* dat)
    {
        std::valarray<std::complex<T> > retval(len);
	for(uint k = 0; k < len-1; ++k)
	  retval[k] = dat[k]/ dat[len-1];
        return retval;
    }
    typedef std::valarray<std::complex<T> > res_t;
    PlainVecSign< std::complex< T > >(uint l) : len(l) {}
private:
  uint len;
};

template <typename FPType>
inline mc_analysis::errordata<std::complex<FPType> > handlesign(FileVector<std::complex<FPType> >& a, FileVector<std::complex<FPType> >& b)
{
    typedef std::complex<FPType> SignType;
    std::valarray<std::complex<FPType> > cont[2];
    cont[0].resize(a.size());
    cont[1].resize(b.size());
    for (std::size_t k = 0; k < a.size(); ++k)
    {
        cont[0][k] = a[k];
        cont[1][k] = b[k];
    }
    mc_analysis::errordata<SignType> ed = mc_analysis::jackknife(ReFunction<FPType>(2), cont, 2);
    mc_analysis::errordata<SignType> edim = mc_analysis::jackknife(ImFunction<FPType>(2), cont, 2);
    return mc_analysis::errordata<std::complex<FPType> >(std::complex<FPType>(ed.get_mean().real(), edim.get_mean().real()), std::complex<FPType>(ed.get_error().real(), edim.get_error().real()), std::complex<FPType>(ed.get_bias().real(), edim.get_bias().real()));
}

template <typename FPType>
inline mc_analysis::errordata<std::complex<FPType> > handlesign(FileVector<std::complex<FPType> >& a, std::valarray<std::complex<FPType> >& b)
{
    typedef std::complex<FPType> SignType;
    std::valarray<std::complex<FPType> > cont[2];
    cont[0].resize(a.size());
    cont[1].resize(b.size());
    for (typename FileVector<std::complex<FPType> >::size_type k = 0; k < a.size(); ++k)
    {
        cont[0][k] = a[k];
        cont[1][k] = b[k];
    }
    mc_analysis::errordata<typename ReFunction<FPType>::res_t> ed = mc_analysis::jackknife(ReFunction<FPType>(2), (&cont[0]), 2);
    mc_analysis::errordata<typename ImFunction<FPType>::res_t> edim = mc_analysis::jackknife(ImFunction<FPType>(2), (&cont[0]), 2);
    return mc_analysis::errordata<std::complex<FPType> >(std::complex<FPType>(ed.get_mean(), edim.get_mean()), std::complex<FPType>(ed.get_error(), edim.get_error()), std::complex<FPType>(ed.get_bias(), edim.get_bias()));
}

template <typename FPType>
inline mc_analysis::errordata<std::complex<FPType> > handlesign(FileVector<std::complex<FPType> >& a, std::valarray<FPType>& b)
{//the twist is here that the observable is a complex quantity whereas the sign is real
    typedef FPType SignType;
    std::valarray<std::complex<FPType> > cont[2];
    cont[0].resize(a.size());
    cont[1].resize(b.size());
    for (typename FileVector<std::complex<FPType> >::size_type k = 0; k < a.size(); ++k)
    {
        cont[0][k] = a[k];
        cont[1][k] = b[k];
    }
    return mc_analysis::jackknife(PlainSign<std::complex<FPType> >(2), cont, 2);
}

template <typename FPType>
mc_analysis::errordata<std::complex<FPType> > handlesign(const std::valarray<std::complex<FPType> >& a, const std::valarray<FPType>& b)
{
    typedef std::complex<FPType> ObsType;
    std::valarray<ObsType> cont[2];
    cont[0].resize(a.size());
    cont[1].resize(b.size());
    for (std::size_t k = 0; k < a.size(); ++k)
    {
        cont[0][k] = a[k];
        cont[1][k] = b[k];
    }
    return mc_analysis::jackknife (PlainSign<ObsType>(2), cont, 2);
}

template <typename FPType>
mc_analysis::errordata<FPType> handlesign(const std::valarray<FPType>& a, const std::valarray<FPType>& b)
{
    typedef FPType SignType;
    const std::valarray<SignType> *const cont[2] = {&a, &b};
    return mc_analysis::jackknife (PlainSign<FPType>(2), cont, 2);
}

template <typename FPType>
inline mc_analysis::errordata<FPType> handlesign(FileVector<FPType>& a, std::valarray<FPType>& b)
{
    typedef FPType SignType;
    std::valarray<FPType> cont[2];
    cont[0].resize(a.size());
    cont[1].resize(b.size());
    for (typename FileVector<FPType>::size_type k = 0; k < a.size(); ++k)
    {
        cont[0][k] = a[k];
        cont[1][k] = b[k];
    }
    return mc_analysis::jackknife(PlainSign<FPType>(2), cont, 2);
}

template <typename FPType>
mc_analysis::errordata<std::complex<FPType> > handlesign(const std::valarray<std::complex<FPType> >& a, const std::valarray<std::complex<FPType> >& b)
{
    const std::valarray<std::complex<FPType> > *const cont[2] = {&a, &b};
    mc_analysis::errordata<typename ReFunction<FPType>::res_t> ed = mc_analysis::jackknife (ReFunction<FPType>(2), cont, 2);
    mc_analysis::errordata<typename ImFunction<FPType>::res_t> edim = mc_analysis::jackknife(ImFunction<FPType>(2), cont, 2);
    return mc_analysis::errordata<std::complex<FPType> >(std::complex<FPType>(ed.get_mean(), edim.get_mean()), std::complex<FPType>(ed.get_error(), edim.get_error()), std::complex<FPType>(ed.get_bias(), edim.get_bias()));
}

template <typename FPType>
mc_analysis::errordata<std::valarray<std::complex<FPType> > > handlesign(const std::valarray<std::valarray<std::complex<FPType> > >& a, const std::valarray<FPType>& b)
{
    std::valarray<std::complex<FPType> > sign(b.size());
    for(uint k = 0; k < b.size(); ++k)
      sign[k] = b[k];
    const std::valarray<std::complex<FPType> > ** cont = new const std::valarray<std::complex<FPType> > *[a.size() + 1];
    for(uint k = 0; k < a.size(); ++k)
      cont[k] = &(a[k]);
    cont[a.size()] = &sign;
    mc_analysis::errordata<typename PlainVecSign<std::complex<FPType> >::res_t> ed = mc_analysis::jackknife(PlainVecSign<std::complex<FPType> >(a.size() + 1), cont, a.size() + 1);
    return ed;
}

template <typename FPType>
mc_analysis::errordata<std::valarray<std::complex<FPType> > > handlesign(const std::valarray<std::valarray<std::complex<FPType> > >& a, bool covariance)
{
  uint funnr = a[0].size();
    mc_analysis::errordata<typename ReVecFunction<FPType>::res_t> ed = mc_analysis::vecjackknife(ReVecFunction<FPType>(funnr), a, funnr, covariance);
    mc_analysis::errordata<typename ImVecFunction<FPType>::res_t> edim = mc_analysis::vecjackknife(ImVecFunction<FPType >(funnr), a, funnr, covariance);
    std::valarray<std::complex<FPType> > retval(funnr-1);
    std::valarray<std::complex<FPType> > retvalerr(funnr-1);
    std::valarray<std::complex<FPType> > retvalbias(funnr-1);
    std::valarray<std::complex<FPType> > retvalcov((funnr - 1) * (funnr - 1));
    for(uint y = 0; y < (funnr-1); ++y)
    {
      retval[y] = std::complex<FPType>(ed.get_mean()[y], edim.get_mean()[y]);
      retvalerr[y] = std::complex<FPType>(ed.get_error()[y], edim.get_error()[y]);
      retvalbias[y] = std::complex<FPType>(ed.get_bias()[y], edim.get_bias()[y]);
      if(covariance)
      {
	for(uint j = 0; j < (funnr - 1); ++j)
	{
	  retvalcov[y*(funnr - 1) + j] = std::complex<FPType>(ed.getCov()[y*funnr + j], edim.getCov()[y*funnr + j]);
	}
      }
    }
    mc_analysis::errordata<std::valarray<std::complex<FPType> > > finalretval(retval, retvalerr, retvalbias);
    if(covariance)
      finalretval.setCov(retvalcov);
    return finalretval;
}
/**
 * This handlesign deals with the preparation of the covariance in the case of a real function.
 * The array "a" contains in the last row the sign.
 * @param a an array containing the function and in the last index the average sign
 * @param covariance Shall we calculate a covariance matrix
 * @return errors and covariances
 * */
template <typename FPType>
mc_analysis::errordata<std::valarray<FPType> > handlesign(const std::valarray<std::valarray<FPType> >& a, bool covariance)
{
    uint funnr = a[0].size();
    mc_analysis::errordata<typename PlainVecSign<FPType>::res_t> ed = mc_analysis::vecjackknife(PlainVecSign<FPType>(funnr), a, funnr, covariance);
    std::valarray<FPType> retval(funnr-1);
    std::valarray<FPType> retvalerr(funnr-1);
    std::valarray<FPType> retvalbias(funnr-1);
    std::valarray<FPType> retvalcov((funnr - 1) * (funnr - 1));
    for(uint y = 0; y < (funnr-1); ++y)
    {
      retval[y] = ed.get_mean()[y];
      retvalerr[y] = ed.get_error()[y];
      retvalbias[y] = ed.get_bias()[y];
      if(covariance)
      {
	for(uint j = 0; j < (funnr - 1); ++j)
	{
	  retvalcov[y*(funnr - 1) + j] = ed.getCov()[y*funnr + j];
	}
      }
    }
    mc_analysis::errordata<std::valarray<FPType> > finalretval(retval, retvalerr, retvalbias);
    if(covariance)
      finalretval.setCov(retvalcov);
    return finalretval;
}

template <typename FPType>
mc_analysis::errordata<std::valarray<std::complex<FPType> > > handlesign(const FileVector<std::valarray<std::complex<FPType> > >& a, const std::valarray<std::complex<FPType> >& b)
{
    const std::valarray<std::complex<FPType> > *const* cont = new const std::valarray<std::complex<FPType> > *const[a.size() + 1];
    for(uint k = 0; k < a.size(); ++k)
      cont[k] = &(a[k]);
    cont[a.size()] = &b;
    mc_analysis::errordata<typename ReVecFunction<FPType>::res_t> ed = mc_analysis::jackknife (ReFunction<FPType>(a.size() + 1), cont, a.size() + 1);
    mc_analysis::errordata<typename ImVecFunction<FPType>::res_t> edim = mc_analysis::jackknife(ImFunction<FPType>(a.size() + 1), cont, a.size() + 1);
    std::valarray<std::complex<FPType> > retval(a.size());
    std::valarray<std::complex<FPType> > retvalerr(a.size());
    std::valarray<std::complex<FPType> > retvalbias(a.size());
    for(uint k = 0; k < a.size(); ++k)
    {
      retval[k] = std::complex<FPType>(ed.get_mean()[k], edim.get_mean()[k]);
      retvalerr[k] = std::complex<FPType>(ed.get_error()[k], edim.get_error()[k]);
      retvalbias[k] = std::complex<FPType>(ed.get_bias()[k], edim.get_bias()[k]);
    }
    return mc_analysis::errordata<std::complex<FPType> >(retval, retvalerr, retvalbias);
}

template<class ObsType, typename SignType>
class Error_Analysis
{
public:
    /**
    @param p the path to which we store
    @param obs the observable
    @param avsign the AverageSign used for determining the value of the observable
    */
    static void analyze(const std::string& p, Observable<ObsType, SignType>& obs, AverageSign<SignType>& avsign);
private:
};

template<class ObsType, typename SignType>
void Error_Analysis<ObsType, SignType>::analyze(const std::string& p, Observable<ObsType, SignType>& obs, AverageSign<SignType>& avsign)
{
    std::ofstream out((p + obs.myname).c_str());
    mc_analysis::errordata<ObsType> ed = handlesign<typename TypetoFPType<SignType>::Ret >(obs.fv, avsign.signcache);
    ObsType mean = ed.get_mean();
    ObsType err = ed.get_error();
    out<<mean<<" +- "<<err<<"      Bias: "<<ed.get_bias()<<std::endl;
}

template<class T, typename SignType>
class Error_Analysis<std::valarray<T>, SignType>
{
public:
    /**
    @param p the path to which we store
    @param obs the Observable
    @param avsign the AverageSign used for determining the value of the observable
    */
    static void analyze(const std::string& p, Observable<std::valarray<T>, SignType>& obs, AverageSign<SignType>& avsign);
private:
};

template <typename FPType>
struct FPIndex
{
  static inline FPType call(uint j, FPType delta_s ) {return j * delta_s;}
};

template <typename FPType>
struct IntIndex
{
  static inline FPType call(uint j, FPType) {return j;}
};

/**
@param fileRe the file to which we write the realpart
@param fileIm the file to which we write the imaginarypart
@param slices the number of points the function has
@param delta_s the resolution of the function
*/
template <class Cont, class SignType, class FPType, class IndexTrait>
inline void writeFunctiontoFile(Cont& container, const valarray<SignType>& avs, std::ofstream& fileRe, std::ofstream& fileIm, unsigned int functionpoints, FPType delta_s, IndexTrait indextrait)
{//Cont is usually a file vector
    typedef typename Cont::value_type::value_type ScalarType;
    unsigned int nroffunctions = container.size();//get nr. of functions
    std::valarray<ScalarType> temp(nroffunctions);//allocate a valarray for one function
    for (unsigned int j = 0; j < functionpoints; ++j)//for all functionpoints
    {
        for (unsigned int k = 0; k < nroffunctions; ++k)//copy function points
            temp[k] = container[k][j];
        mc_analysis::errordata<ScalarType> edp(handlesign(temp, avs));//do sign - analysis
        //output data
	FPType idx = IndexTrait::call(j, delta_s);
        fileRe<<idx<<" "<<real(edp.get_mean())<<" "<<real(edp.get_error())<<std::endl;
        fileIm<<idx<<" "<<imag(edp.get_mean())<<" "<<imag(edp.get_error())<<std::endl;
    }
    return;
}

/**
@param file the file to which we write the realpart
@param slices the number of points the function has
@param delta_s the resolution of the function
*/
template <class Cont, class FPType>
inline void writeFunctiontoFile(Cont& container, const valarray<FPType>& avs, std::ofstream& file, unsigned int functionpoints, FPType delta_s)
{
    unsigned int nroffunctions = container.size();//get nr. of functions
    std::valarray<FPType> temp(nroffunctions);//allocate a valarray for one function

    for (unsigned int j = 0; j < functionpoints; ++j)//for all functionpoints
    {
        for (unsigned int k = 0; k < nroffunctions; ++k)//copy function points
            temp[k] = container[k][j];
        mc_analysis::errordata<FPType> edp(handlesign(temp, avs));//do sign - analysis
        //output data
        file<<j * delta_s<<" "<<edp.get_mean()<<" "<<edp.get_error()<<std::endl;
    }
    return;
}

template <class FPType>
struct Point
{
    Point(FPType a, unsigned int b, FPType c, FPType d) : t(a), x(b), y(c), dy(d) {}
    FPType t;
    unsigned int x;
    FPType y,dy;
    inline bool operator<(const Point& rhs) const
    {
        if (t >= rhs.t)
        {
            if (t == rhs.t)
            {
                if (x >= rhs.x)
                {
                    return false;
                }
                return true;
            }
            return false;
        }
        return true;
    }
};

/**
 * This gets called for every complex Function, e.g. a complex Green's function
@param fileRe the file to which we write the realpart
@param fileIm the file to which we write the imaginarypart
@param len the length of the Vector
@param slices the number of points the function has
@param delta_s the resolution of the function
*/
template <class Cont, class SignType, typename FPType>
static inline void writeVectorFunctiontoFile(Cont& container, unsigned int len, const valarray<SignType>& avsign, std::string& nameRe, std::string& nameIm, unsigned int slices, FPType delta_s, bool covariance)
{
    std::ofstream fileRe(nameRe.c_str());
    std::ofstream fileIm(nameIm.c_str());
    std::ofstream covfileRe((nameRe+"Cov").c_str());
    std::ofstream covfileIm((nameIm+"Cov").c_str());
    //the template parameter Cont is usually a filevector
    typedef typename Cont::value_type::value_type::value_type ScalarType;
    std::ofstream xdep[5];
    for (unsigned int i = 0; i < 5; ++i)
    {
        unsigned int dt = 2 << i;
        xdep[i].open((nameRe + "_" + toString(dt)).c_str());
    }
    container.sync();
    std::vector<Point<FPType> > points;
    points.reserve(4*slices);
    unsigned int nroffunctions = container.size();//get nr. of bins
    std::valarray<std::valarray<ScalarType> > temp(nroffunctions + 1);//allocate a valarray for one function
    for(uint j = 0; j < (nroffunctions+1); ++j)
      temp[j].resize(slices+1);
    for(uint k = 0; k < len; ++k) // for every point
    {
      for (uint j = 0; j < nroffunctions; ++j) //get every bin
      {
         auto ret = container(j, k);
         for(uint i = 0; i < ret.size(); ++i)
            temp[j][i] = ret[i];
	 temp[j][ret.size()] = avsign[j];//we copy the sign bins into the last index
      }
      mc_analysis::errordata<std::valarray<ScalarType> > edp(handlesign(temp, covariance));//do sign - analysis
      for(uint l = 0; l < slices; ++l)
      {
            //output data
            points.push_back(Point<FPType>(l * delta_s, k, real(edp.get_mean()[l]), real(edp.get_error()[l])));
            fileRe<<l * delta_s<<" "<<real(edp.get_mean()[l])<<" "<<real(edp.get_error()[l])<<std::endl;
            fileIm<<l * delta_s<<" "<<imag(edp.get_mean()[l])<<" "<<imag(edp.get_error()[l])<<std::endl;
      }
        //output xmgrace separators
        fileRe<<"&"<<std::endl;
        fileIm<<"&"<<std::endl;
    if(covariance)
    {
      for(uint k = 0; k < slices; ++k)
      {
	for(uint j = 0; j < slices; ++j)
	{
	  covfileRe<<real(edp.getCov()[k*slices + j])<<" ";
	  covfileIm<<imag(edp.getCov()[k*slices + j])<<" ";
	}
	covfileRe<<std::endl;
	covfileIm<<std::endl;
      }
    }
    }
//     for (unsigned int k = 0; k < len; ++k)//for every function that we stored
//     {
//         for (unsigned int j = 0; j < slices; ++j)//for all functionpoints
//         {
//             for (unsigned int l = 0; l < nroffunctions; ++l)//copy function points
//                 temp[l] = container(l,k)[j];
//             mc_analysis::errordata<FPType> edp(handlesign(temp, avsign));//do sign - analysis
//             //output data
//             points.push_back(Point<FPType>(j * delta_s, k, edp.get_mean(), edp.get_error()));
//             file<<j * delta_s<<" "<<edp.get_mean()<<" "<<edp.get_error()<<std::endl;
//         }
//         //output xmgrace separators
//         file<<"&"<<std::endl;
//     }
    //now let's generate the transformed data
    sort(points.begin(), points.end());
#pragma omp parallel for
    for(size_t i = 0; i < 5; ++i)
    {
       unsigned int inc = 2 << i;
       for(unsigned int k = 0; k < slices; k += inc)
       {
          for(unsigned int j = 0; j < len; ++j)
          {
             xdep[i]<<j<<" "<<points[k*len+j].y<<" "<<points[k*len+j].dy<<std::endl;
          }
          xdep[i]<<"&"<<std::endl;
       }
    }
    return;
}

// /**
// @param fileRe the file to which we write the realpart
// @param fileIm the file to which we write the imaginarypart
// @param len the length of the Vector
// @param slices the number of points the function has
// @param delta_s the resolution of the function
// */
// template <class Cont, class SignType, typename FPType>
// static inline void writeVectorFunctiontoFile(Cont& container, unsigned int len, const valarray<SignType>& avsign, std::string& nameRe, std::string& nameIm, unsigned int slices, FPType delta_s)
// {
//     std::ofstream fileRe(nameRe.c_str());
//     std::ofstream fileIm(nameIm.c_str());
//     //the template parameter Cont is usually a filevector
//     typedef typename Cont::value_type::value_type::value_type ScalarType;
//     std::ofstream xdep[5];
//     for (unsigned int i = 0; i < 5; ++i)
//     {
//         unsigned int dt = 2 << i;
//         xdep[i].open((nameRe + "_" + toString(dt)).c_str());
//     }
//     container.sync();
//     std::vector<Point<FPType> > points;
//     points.reserve(4*slices);
//     unsigned int nrbins = container.size();//get nr. of bins
//     uint availmem = sysconf (_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE);//some approximation to the amount of available memory
//     uint memoryperpoint = sizeof(ScalarType)*nrbins;
//     int possibleparallelism = std::min(availmem / memoryperpoint, slices);
//     std::valarray<ScalarType>* temp = NULL;
//     bool memcritical = (possibleparallelism == 0);
//     std::cout<<"Possible parallelism: "<<possibleparallelism<<std::endl;
//     if(!memcritical)
//     {
//       try
//       {
//       temp = new std::valarray<ScalarType>[possibleparallelism];
//       for(int j = 0; j < possibleparallelism; ++j)
// 	temp[j].resize(nrbins);
//       }
//       catch(std::bad_alloc& e)
//       {
// 	delete [] temp;
// 	memcritical = true;
//       }
//     }
//     if(memcritical)
//     {//memory is critical. let's try it anyway and see if the kernel frees memory. 
//       possibleparallelism = 1;
//       temp = new std::valarray<ScalarType>(nrbins);
//     }
// //    std::valarray<ScalarType> temp(nrbins);//allocate a valarray for one function    
//     for (unsigned int k = 0; k < len; ++k)//for every index of the function that we stored
//     {
//         for (unsigned int j = 0; j < slices; j += possibleparallelism)//for all functionpoints
//         {
// 	  int upperlimit = possibleparallelism;
// 	  if(j + possibleparallelism >= slices)
// 	    upperlimit = slices - j;//determine the remaining indices
//           for (unsigned int l = 0; l < nrbins; ++l)//copy function points
//             {
// 	      auto tempdatafromfile = container(l, k);
// 	      for(int j1 = 0; j1 < upperlimit; ++j1)
// 		temp[j1][l] = tempdatafromfile[j+j1];
// //	        temp[l] = container(l, k)[j];
//             }
//             for(int j1 = 0; j1 < upperlimit; ++j1)
// 	    {
//             mc_analysis::errordata<ScalarType> edp(handlesign(temp[j1], avsign));//do sign - analysis
//             //output data
//             points.push_back(Point<FPType>((j+j1) * delta_s, k, real(edp.get_mean()), real(edp.get_error())));
//             fileRe<<(j+j1) * delta_s<<" "<<real(edp.get_mean())<<" "<<real(edp.get_error())<<std::endl;
//             fileIm<<(j+j1) * delta_s<<" "<<imag(edp.get_mean())<<" "<<imag(edp.get_error())<<std::endl;
// 	    }
//         }
//         //output xmgrace separators
//         fileRe<<"&"<<std::endl;
//         fileIm<<"&"<<std::endl;
//     }
//     delete [] temp;
//     //now let's generate the transformed data
//     sort(points.begin(), points.end());
// #pragma omp parallel for schedule(dynamic)
//     for(size_t i = 0; i < 5; ++i)
//     {
//        unsigned int inc = 2 << i;
//        for(unsigned int k = 0; k < slices; k += inc)
//        {
//           for(unsigned int j = 0; j < len; ++j)
//           {
//              xdep[i]<<j<<" "<<points[k*len+j].y<<" "<<points[k*len+j].dy<<std::endl;
//           }
//           xdep[i]<<"&"<<std::endl;
//        }
//     }
//     return;
// }

/**
 * This gets called for every real function, e.g. The Green's function if it is real.
@param name the file to which we write the data
@param len the length of the Vector
@param slices the number of points the function has
@param delta_s the resolution of the function
*/
template <class Cont, typename FPType>
static inline void writeVectorFunctiontoFile(Cont& container, unsigned int len, const valarray<FPType>& avsign, std::string& name, unsigned int slices, FPType delta_s, bool covariance)
{
  //the template parameter Cont is usually a filevector
    typedef typename Cont::value_type::value_type::value_type ScalarType;
    std::ofstream file(name.c_str());
    std::ofstream covfile((name + "_Cov").c_str());
    std::ofstream xdep[5];
    for (unsigned int i = 0; i < 5; ++i)
    {
        unsigned int dt = 2 << i;
        xdep[i].open((name + "_" +toString(dt)).c_str());
    }
    container.sync();
    std::vector<Point<FPType> > points;
    points.reserve(4*slices);
    unsigned int nroffunctions = container.size();//get nr. of bins
    std::valarray<std::valarray<ScalarType> > temp(nroffunctions + 1);//allocate a valarray for one function + the sign
    for(uint j = 0; j < (nroffunctions + 1); ++j)
      temp[j].resize(slices + 1);
    for (unsigned int k = 0; k < len; ++k)//for every function that we stored
    {
      for(uint j = 0; j < nroffunctions; ++j) //get every bin
      {
	auto ret = container(j, k);
	for(uint i = 0; i < ret.size(); ++i)
	  temp[j][i] = ret[i];
	temp[j][ret.size()] = avsign[j];
      }
      mc_analysis::errordata< std::valarray< ScalarType > > edp(handlesign(temp, covariance));//do sign - analysis
      for(uint l = 0; l < slices; ++l)
      {
	//output data
	points.push_back(Point<FPType>(l * delta_s, k, edp.get_mean()[l], edp.get_error()[l]));
	file<<l * delta_s<<" "<<edp.get_mean()[l]<<" "<<edp.get_error()[l]<<std::endl;
      }
        //output xmgrace separators
      file<<"&"<<std::endl;
      if(covariance)
      {
	for(uint k = 0; k < slices; ++k)
	{
	  for(uint j = 0; j < slices; ++j)
	    covfile<<edp.getCov()[k*slices + j]<<" ";
	  covfile<<std::endl;
	}
      }
    }
    //now let's generate the transformed data
    sort(points.begin(), points.end());
#pragma omp parallel for
    for(size_t i = 0; i < 5; ++i)
    {
       unsigned int inc = 2 << i;
       for(unsigned int k = 0; k < slices; k += inc)
       {
          for(unsigned int j = 0; j < len; ++j)
          {
             xdep[i]<<j<<" "<<points[k*len+j].y<<" "<<points[k*len+j].dy<<std::endl;
          }
          xdep[i]<<"&"<<std::endl;
       }
    }
    return;
}


template<class T, typename SignType>
void Error_Analysis<std::valarray<T>, SignType>::analyze(const std::string& p, Observable<std::valarray<T>, SignType>& obs, AverageSign<SignType>& avsign)
{
    std::ofstream outRe((p + obs.myname + string("_Re")).c_str());
    std::ofstream outIm((p + obs.myname + string("_Im")).c_str());
    typedef typename remove_const<typename remove_reference<decltype(obs.delta_s)>::type>::type mytype;
    if(obs.myname.compare(0, 10, "KondoCloud") == 0)
      writeFunctiontoFile(obs.fv, avsign.signcache, outRe, outIm, obs.functionpoints, obs.delta_s, IntIndex<mytype>());
    else
      writeFunctiontoFile(obs.fv, avsign.signcache, outRe, outIm, obs.functionpoints, obs.delta_s, FPIndex<mytype>());
    return;
}

template<>
void Error_Analysis<std::valarray<double>, double>::analyze(const std::string& p, Observable<std::valarray<double>, double>& obs, AverageSign<double>& avsign)
{
    std::ofstream outRe((p + obs.myname).c_str());
    writeFunctiontoFile(obs.fv, avsign.signcache, outRe, obs.functionpoints, obs.delta_s);
    return;
}

template<class T, typename SignType>
class Error_Analysis<std::valarray<std::valarray<T> >, SignType>
{
public:
    /**
    @param p the path to which we store
    @param obs the Observable
    @param avsign the AverageSign used for determining the value of the observable
    */
    static void analyze(const std::string& p, Observable<std::valarray<std::valarray<T> >, SignType>& obs, AverageSign<SignType>& avsign);
private:
};

template<class T, typename SignType>
void Error_Analysis<std::valarray<std::valarray<T> >, SignType>::analyze(const std::string& p, Observable<std::valarray<std::valarray<T> >, SignType>& obs, AverageSign<SignType>& avsign)
{
    std::string outRe(p + obs.myname + string("_Re"));//create file for outputting the real part
    std::string outIm(p + obs.myname + string("_Im"));//create file for outputting the imaginary part
    writeVectorFunctiontoFile(obs.fv, obs.tensorindices, avsign.signcache, outRe, outIm, obs.functionpoints, obs.delta_s, obs.covariance);//well, write to files...
}

template<>
void Error_Analysis<std::valarray<std::valarray<double> >, double>::analyze(const std::string& p, Observable<std::valarray<std::valarray<double> >, double>& obs, AverageSign<double>& avsign)
{
    std::string outRe(p + obs.myname);//create file for outputting the real part
    writeVectorFunctiontoFile(obs.fv, obs.tensorindices, avsign.signcache, outRe, obs.functionpoints, obs.delta_s, obs.covariance);//well, write to files...
}
#endif
