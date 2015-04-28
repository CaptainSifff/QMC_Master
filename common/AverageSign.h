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
#ifndef AVERAGESIGN_H
#define AVERAGESIGN_H
#include "ParallelServer.h"
#include "Master.h"
#include "filevector.h"
#include <fstream>
#include <string>
#include <algorithm>
#include "bootstrap.h"

template <class T>
struct TypetoFPType
{
    typedef T Ret;
};

template <class T>
struct TypetoFPType<std::complex<T> >
{
    typedef T Ret;
};

template <class SignType>
class MeanRe
    {
    public:
        typename TypetoFPType<SignType>::Ret operator() (SignType** data, unsigned int len)
        {
            typename TypetoFPType<SignType>::Ret RetVal = 0;
            for (unsigned int l = 0; l < len; ++l)
                RetVal += data[l][0].real();
            return RetVal/static_cast<typename TypetoFPType<SignType>::Ret>(len);
        }
        typedef typename TypetoFPType<SignType>::Ret RetType;
    };

template <class SignType>
    class MeanIm
    {
    public:
        typename TypetoFPType<SignType>::Ret operator() (SignType** data, unsigned int len)
        {
            typename TypetoFPType<SignType>::Ret RetVal = 0;
            for (unsigned int l = 0; l < len; ++l)
                RetVal += data[l][0].imag();
            return RetVal/static_cast<typename TypetoFPType<SignType>::Ret>(len);
        }
        typedef typename TypetoFPType<SignType>::Ret RetType;
    };

template <class FPType>
class Mean
{
public:
        FPType operator() (FPType** data, unsigned int len)
        {
            FPType RetVal = 0;
            for (unsigned int l = 0; l < len; ++l)
                RetVal += data[l][0];
            return RetVal/static_cast<FPType>(len);
        }
typedef FPType RetType;
};

template <class SignType>
class AverageSign
{
public:
    template <class T, class S> friend class Error_Analysis;
    void commit();
    void cacheData(Master&);
    void analyzeData(const std::string&);
    AverageSign(const std::string&);
    unsigned int size()
    {
        return fv.size();
    }
    inline void sync() {fv.sync();}
private:
    SignType cache;
    FileVector<SignType> fv;
    std::valarray<SignType> signcache; 
};

template <class SignType>
void AverageSign<SignType>::analyzeData(const std::string& p)
{
    std::ofstream file((p + std::string("AverageSign.dat")).c_str());
    //load bins into memory
    signcache.resize(fv.size());
    unsigned int k = 0;
    for(typename FileVector<SignType>::const_iterator it = fv.begin(); it != fv.end(); ++it, ++k)
      signcache[k] = *it;
    mc_analysis::errordata<typename TypetoFPType<SignType>::Ret> meanre = mc_analysis::bootstrap(MeanRe<SignType>(), &signcache, 1);
    mc_analysis::errordata<typename TypetoFPType<SignType>::Ret> meanim = mc_analysis::bootstrap(MeanIm<SignType>(), &signcache, 1);
    mc_analysis::errordata<SignType> ed(0.0, 0.0);
    ed.set_mean(SignType(meanre.get_mean(), meanim.get_mean()));
    ed.set_error(SignType(meanre.get_error(), meanim.get_error()));
    file<<"AverageSign: "<<ed.get_mean()<<" +- "<<ed.get_error()<<endl;
}

template <>
void AverageSign<double>::analyzeData(const std::string& p)
{
    std::ofstream file((p + std::string("AverageSign.dat")).c_str());
    //load bins into memory
    signcache.resize(fv.size());
    unsigned int k = 0;
    for(FileVector<double>::const_iterator it = fv.begin(); it != fv.end(); ++it, ++k)
      signcache[k] = *it;
    mc_analysis::errordata<double> mean = mc_analysis::bootstrap(Mean<double>(), &signcache, 1);
    file<<"AverageSign: "<<mean.get_mean()<<" +- "<<mean.get_error()<<endl;
}

template <class ObsT>
AverageSign<ObsT>::AverageSign(const std::string& bp)
{
    fv.load(bp + std::string("AverageSign"));//FIXME: we load from binpath + name of observable
}

template <class SignType>
void AverageSign<SignType>::commit()
{
    fv.push_back(cache);
}

template <class SignType>
void AverageSign<SignType>::cacheData(Master& m)
{
    cache = m.recvfromanyclient<SignType>().msg;
}
#endif
