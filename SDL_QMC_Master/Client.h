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
#ifndef CLIENT_H
#define CLIENT_H
#include <string>

typedef struct _TCPsocket *TCPsocket;

class Client
{
public:
    inline const std::string& getName() const;
    inline const TCPsocket& getSocket() const throw();
    inline Client(TCPsocket&, std::string&);
    inline Client() {}
    inline Client(const Client&);
    inline bool operator==(const Client& a) const throw();
    inline bool operator!=(const Client& a) const throw();
private:
    TCPsocket sock;
    std::string name;
};

Client::Client(const Client& rhs) : sock(rhs.sock), name(rhs.name)
{
}

Client::Client(TCPsocket& t, std::string& s) : sock(t), name(s)
{
}

const std::string& Client::getName() const
{
    return name;
}

const TCPsocket& Client::getSocket() const throw()
{
    return sock;
}

/**
Compares two Clients.
@return true if it is the same socket and the same Client-Name
*/
bool Client::operator==(const Client& a) const throw()
{
    return (a.sock == sock) && (a.name == name);
}

bool Client::operator!=(const Client& a) const throw()
{
    return !(*this == a);
}
#endif
