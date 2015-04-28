/***************************************************************************
 *   Copyright (C) 2008,2009 by Florian Goth   *
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
#ifndef NETCOMMONBASE_H
#define NETCOMMONBASE_H
#include "stripped_SDL_Net.h"
#include "Client.h"
#include "NetError.h"
#include <cstring>
#include <iostream>

#define MAX_MTU 40000

class TimeOutError : public NetError
{
public:
    TimeOutError(std::string a, TCPsocket s) : NetError(a), sock(s){}
    TCPsocket sock;
private:
};

class SDL_Error : public SubSystemError
{
public:
    SDL_Error(std::string a) : SubSystemError(a) {}
};

// send a string buffer over a TCP socket with error checking
// returns 0 on any errors, length sent on success
inline unsigned int csend(TCPsocket sock, uint32_t len, const char* buf)
{
    // change endianness to network order
    len = SDL_SwapBE32(len);

    // send the length of the string
    uint32_t result = SDLNet_TCP_Send(sock,&len,sizeof(len));
    if (result < sizeof(len))
    {
        if (SDLNet_GetError() && strlen(SDLNet_GetError())) // sometimes blank!
            printf("SDLNet_TCP_Send: %s\n", SDLNet_GetError());
        return 0;
    }

    // revert to our local byte order
    len = SDL_SwapBE32(len);

    // send the buffer, with the NULL as well
    result = SDLNet_TCP_Send(sock,buf,len);
    if (result < len)
    {
        if (SDLNet_GetError() && strlen(SDLNet_GetError())) // sometimes blank!
            printf("SDLNet_TCP_Send: %s\n", SDLNet_GetError());
        return 0;
    }

    // return the length sent
    return result;
}

inline void blockrecv(TCPsocket sock, char *const data, const unsigned int len)
{
    // get the string buffer over the socket
    unsigned int pos = 0;
    unsigned int remlen = len;
    unsigned int delay = 2;
    do
    {
        unsigned int result = SDLNet_TCP_Recv(sock, &(data[pos]), remlen);
        //we have read result bytes
        pos += result;
        remlen -= result;
        if (remlen != 0)
        {
            SDL_Delay(delay);//wait two ms
            std::cout<<"waiting for data..."<<std::endl;
            if (result == 0)
            {
                delay <<= 1;
            }
            else
            {
                delay = 2;
            }
            if (delay > 500)//we didn't receive any data in about half a second -> the connection timed out
            {
                throw(TimeOutError("Connection timed out!", sock));
            }
        }
    }
    while (remlen != 0);
    return;
}

/**
@param sock the socket to use for receiving
@param buf A pointer to the beginning of an uninitialized array. The array gets allocated by recv.
@return the new array
*/
inline char* crecv(TCPsocket sock, char* data, uint32_t& len)
{
    char** buf = &data;
    static char *_buf;

    // allow for a NULL buf, use a static internal one...
    if (!buf)
        buf=&_buf;

    // free the old buffer
    if (*buf)
        delete [] *buf;
    *buf = NULL;

    // receive the length of the message
    char* ptr = reinterpret_cast<char*>(&len);
    blockrecv(sock, ptr, sizeof(len));

    // swap byte order to our local order
    len = SDL_SwapBE32(len);
    // check if anything is strange, like a zero length buffer
    if (!len)
        return NULL;

    // allocate the buffer memory
    *buf = new char[len];
    if (!(*buf))
        return NULL;
    try
    {
        blockrecv(sock, *buf, len);
    }
    catch (TimeOutError& e)
    {
        if (data == NULL) delete [] *buf;//we assume that the whole buffer is then invalid.
        throw;
    }
    // return the new buffer
    return *buf;
}

class Message
{
public:
    char* msg;
    unsigned int len;
    inline Message(char* a, unsigned int b) : msg(a), len(b) {}
    inline Message() : msg(NULL), len(0) {}
    inline void free();
    inline unsigned int alloc(unsigned int);
};

void Message::free()
{
    delete [] msg;
}

unsigned int Message::alloc(unsigned int l)
{
    msg = new char[l];
    len = l;
    return l;
}

template<typename T>
class NetCommonBase
{
public:
protected:
    static inline Message recv(T, Message&);
    static inline unsigned int send(T, const Message&);
private:
};

// send a string buffer over a TCP socket with error checking
// returns length sent on success
template <>
inline unsigned int NetCommonBase<TCPsocket>::send(TCPsocket sock, const Message& m)
{
    unsigned int sentbytes = 0;
    unsigned int firstslicepayload = MAX_MTU - 6;
    unsigned int nextslicepayload = MAX_MTU - 2;
    if(!m.msg)
        throw(NetError("[SDL_Common] Message pointer not set!"));
    if(!m.len)
        throw(NetError("[SDL_Common] Message has length 0! "));
    if (m.len < firstslicepayload)
    {
        //one packet is sufficient
        char* packet = new char[6 + m.len];
        SDLNet_Write32(m.len, &(packet[0]));
        SDLNet_Write16(1, &(packet[4]));
        memcpy(&packet[6], m.msg, m.len);
        csend(sock, 6 + m.len, packet);
        sentbytes += m.len;
        delete [] packet;
    }
    else
    {
        //we need to send more packets...
        char* packet;
        packet = new char[6 + firstslicepayload];
        uint16_t nrofpackets = 1;
        unsigned int remlen = m.len - firstslicepayload;
        while (remlen > nextslicepayload)
        {
            ++nrofpackets;
            remlen -= nextslicepayload;
        }
        ++nrofpackets;//account for the last non-full packet
        SDLNet_Write32(m.len, &(packet[0]));
        SDLNet_Write16(nrofpackets, &(packet[4]));
        memcpy(&packet[6], m.msg, firstslicepayload);
        csend(sock, 6 + firstslicepayload, packet);
        sentbytes += firstslicepayload;
        delete [] packet;
        const char* msgptr = &m.msg[firstslicepayload];
        //we've sent the first packet.
        packet = new char[2 + nextslicepayload];//the following packets have all the same length
        for (uint16_t k = 1; k < (nrofpackets -  1); ++k)
        {
            SDLNet_Write16(k, packet);
            memcpy(&packet[2], msgptr, nextslicepayload);
            csend(sock, 2 + nextslicepayload, packet);
            sentbytes += nextslicepayload;
            msgptr = &msgptr[nextslicepayload];
        }
        delete [] packet;
        packet = new char[2 + remlen];
	uint16_t arg = static_cast<uint16_t>(nrofpackets -  static_cast<uint16_t>(1u));
        SDLNet_Write16(arg, packet);
        memcpy(&packet[2], msgptr, remlen);
        csend(sock, remlen + 2, packet);
        sentbytes += remlen;
        delete [] packet;
    }
    return sentbytes;
}

/**
@param sock the socket to use for receiving
@param buf A pointer to the beginning of an uninitialized array. The array gets allocated by recv.
@return the new array
*/
template <>
inline Message NetCommonBase<TCPsocket>::recv(TCPsocket sock, Message& m)
{
    uint32_t packetlen = 0;
    char* packet = crecv(sock, NULL, packetlen);//receive the first packet and let's see what's in it
    if (packetlen < 2) throw(NetError("invalid Packet received!"));
    m.alloc(SDLNet_Read32(packet));
    unsigned int nrofpackets = SDLNet_Read16(&packet[4]);
    memcpy(m.msg, &packet[6], packetlen - 6);
    unsigned int recvbytes = packetlen - 6;
    //we're done if the packet contains only one slice
    delete [] packet;
    if (nrofpackets > 1)
    {
        char* msgptr = &m.msg[packetlen - 6];
        for (unsigned int k = 1; k <  nrofpackets; ++k)
        {
            packet = crecv(sock, NULL, packetlen);
            memcpy(msgptr, &packet[2], packetlen - 2);
            recvbytes += packetlen - 2;
            msgptr = &msgptr[packetlen -2];
            delete [] packet;
        }
    }
    // return the new buffer
    if (recvbytes != m.len) throw(NetError("Not enough bytes received!"));
    return m;
}
#endif
