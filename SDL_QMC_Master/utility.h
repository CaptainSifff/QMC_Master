/***************************************************************************
 *   Copyright (C) 2008 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef UTILITY_H
#define UTILITY_H
#include <string>
#include <sstream>
#include <iostream>
#include "stripped_SDL_Net.h"

using namespace std;
/**
@param ip a 4-byte IP in the hosts byte-Order
@return a string represantation of the IP
*/

inline string convertIPtoString(uint32_t ip)
{
    uint32_t ipaddr = SDL_SwapBE32(ip);
    stringstream st[4];
    st[0] << (ipaddr >> 24);
    st[1] << ( (ipaddr >> 16) & 0xff);
    st[2] << ((ipaddr >> 8) & 0xff);
    st[3] << ( ipaddr & 0xff);
    string str[4];
    for (unsigned int k = 0; k < 4 ; ++k)
        st[k]>>str[k];
    string RetVal = str[0] + string(".") + str[1] + string(".") + str[2] + string(".") + str[3];
    return RetVal;
}

/**
Converts a given String to an IP-address in uint format
@param work the string to convert.
@return an ip-address. Returns 0xFFFFFFFF if an error occured
*/

inline uint32_t StringtoIP(const string& work)
{
    //Assumption: IP looks like 127.0.0.1
    uint32_t addr[4] = {0};
    stringstream ss[4];
    size_t dot = 0;
    size_t nr = 0;
    for (unsigned int k = 0; k < 3; ++k)
    {
        nr = work.find_first_of("0123456789", dot);
        dot = work.find_first_of('.', dot+1);
        if ((nr == string::npos ) || (dot == string::npos)) return 0xFFFFFFFF;
        if (nr > dot) return 0xFFFFFFFF;
        ss[k] << work.substr(nr, dot);
        ss[k] >> addr[k];
        if (addr[k] > 255) return 0xFFFFFFFF;
    }

    nr = work.find_first_of("0123456789", dot);
    dot = work.find_first_of('.', dot);
    ss[3] << work.substr(nr, dot);
    ss[3] >> addr[3];
    if (addr[3] > 255) return 0xFFFFFFFF;

    std::cout<<addr[0]<<"."<<addr[1]<<"."<<addr[2]<<"."<<addr[3]<<std::endl;
    return SDL_SwapBE32((addr[0]<<24) | (addr[1] << 16) | (addr[2] << 8) | (addr[3]));
}

/**
Tries to resolve a given string to sth. network usable.
@param arg The string we try to resolve to an IP.
@param testport The port where to test
@return an IP Address on success, else 0xFFFFFFFF
*/

inline uint32_t resolveArgtoIP(const std::string& arg, uint16_t testport)
{
    uint32_t retval = StringtoIP(arg);
    if ( retval == 0xFFFFFFFF)
    {
        IPaddress ipaddr;
        int ret = SDLNet_ResolveHost(&ipaddr, arg.c_str(), testport);
        if (ret == -1)
            return 0xFFFFFFFF;
        retval = ipaddr.host;
    }
    return retval;
}
#endif
