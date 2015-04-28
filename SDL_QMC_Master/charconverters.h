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
#ifndef CHARCONVERTERS_H
#define CHARCONVERTERS_H
#include <complex>
#include <stdint.h>//fixme C++1x has support for this...
#include <string>
#include <cstring>

template<class T>
class CH
{
public:
    inline CH(const T& d) : val(d) {}
    inline CH( const char *const arg)
    {
        for (unsigned int k = 0; k < sizeof(T); ++k)
            bytes[k] = arg[k];
    }
    operator T ()
    {
        return val;
    }
    inline char* accessbytes() throw()
    {
        return bytes;
    }
    inline char accessbytes(unsigned int k) throw()
    {
        return bytes[k];
    }
    inline std::size_t size() const throw()
    {
        return sizeof(T);
    }
private:
    union
    {
        T val;
        char bytes[sizeof(T)];
    };
};

template <class FPType>
class CH<std::complex<FPType> >
{
public:
    inline CH(const std::complex<FPType>& d)
    {
        val[0] = d.real();
        val[1] = d.imag();
    }
    inline CH(const char *const arg)
    {
        for (unsigned int k = 0; k < 2 * sizeof(FPType); ++k) bytes[k] = arg[k];
    }
    operator std::complex<FPType>()
    {
        return std::complex<FPType>(val[0], val[1]);
    }
    inline char* accessbytes() throw()
    {
        return bytes;
    }
    inline char accessbytes(unsigned int k) throw()
    {
        return bytes[k];
    }
    inline std::size_t size() const throw()
    {
        return sizeof(std::complex<FPType>);
    }
private:
    union
    {
        FPType val[2];
        char bytes[2*sizeof(FPType)];
    };
};

template<>
class CH<std::string>
{
public:
    inline CH(const std::string& d) : len(d.size()+1), data(new char[len])
    {
        data[len-1] = 0;
        d.copy(data, len);
    }
    inline CH(const char *const arg) : len(strlen(arg) + 1/*the terminating zero*/), data(new char[len])//to prevent usage of the implicit conversion via string
    {
        strcpy(data, arg);
    }
    operator std::string ()
    {
        return std::string(data);
    }
    inline char* accessbytes() throw()
    {
        return data;
    }
    inline char accessbytes(unsigned int k ) throw()
    {
        return data[k];
    }
    inline std::size_t size() const throw()
    {
        return len;
    }
    ~CH()
    {
        delete [] data;
    }
private:
    std::size_t len;
    char* data;
};
#endif
