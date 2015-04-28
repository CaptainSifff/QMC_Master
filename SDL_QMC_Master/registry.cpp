/***************************************************************************
 *   Copyright (C) 2005-2009 by Florian Goth   *
 *   Captainsifff@gmx.de   *
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
#include "registry.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <dirent.h>
#include <unistd.h>
#include <cstring>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

using namespace std;

// General tool to strip spaces from both ends:
static inline std::string trim(const std::string& s)
{//from "Thinking in C++" by Bruce Eckel
    if(s.length() == 0)
        return s;
    std::string::size_type b = s.find_first_not_of(" \t");
    std::string::size_type e = s.find_last_not_of(" \t");
    if(b == string::npos) // No non-spaces
        return "";
    return std::string(s, b, e - b + 1);
}

void Block::push_back_Key(std::string arg)
{
    string::size_type pos = arg.find_first_of("=");
    string key(arg , 0 , pos );
    string value(arg , pos + 1, arg.size() - pos - 1);
    Block_Data[trim(key)] = trim(value);
    return;
}

Block::Block(std::string blockname, const std::vector<std::string>& block) : BlockName(blockname)
{//Now has an array of all it's arrays
    for ( unsigned int k = 0 ; k < block.size() ; ++ k )
        push_back_Key(block[k]);
    NrOfKeys = Block_Data.size();
    Keys.reserve( NrOfKeys );
    std::map<std::string, std::string>::iterator it;
    for(it = Block_Data.begin(); it != Block_Data.end(); it++)
        Keys.push_back( it->first );
    return;
}

inline pair<string , unsigned int > FindNextBracketPair(vector<string>& file , unsigned int pos ) throw(runtime_error)
{
    for (unsigned int k = pos; k < file.size(); ++k )
    {
        string::size_type first = file[k].find_first_of("[");
        if(  first != string::npos )
        {
            //Block Begin found!
            string::size_type last = file[k].find_last_of("]");
            if (last == string::npos)
                throw(runtime_error("Error in Parsing ConfigFile!! Closing Bracket not found"));
            string Blockname(file[k] , first+1 , last -first - 1 );
            return make_pair( Blockname , k );
        }
    }
    throw(runtime_error("No new Block found"));
}

cfile::cfile(std::string& file ) : filename(file)
{
    ifstream Config;
    vector<string> Zeilen;
    Config.open(file.c_str(), ios_base::in);
    if (!Config)
    {//the destructor of Config gets called because of stack unwinding
        throw(runtime_error(string("Couldn't open Config File : ") + file));
    }
    else
    {
        string line;
        while(getline(Config, line))
            if ( ((line = trim(line)).size() != 0) && (( line = trim(line)).at(0) != '#'))
                Zeilen.push_back( trim(line) ); // Add the line to the end
        //Now all Comments ( those starting with a '#' are removed
        //vector Zeilen contains Config File
    }
    if (std::find(Zeilen.begin(), Zeilen.end(), string("[END]")) == Zeilen.end())
    {
       throw(runtime_error(string("File not recognized as Config File! ") + file));
    }
    //Now Split up in Blocks
    pair<string, unsigned int > Ret("" , 0 );
    while( (Ret.first != "END") && (Ret.second < Zeilen.size()) )
    {
        Ret = FindNextBracketPair(Zeilen , Ret.second);
        unsigned int BlockBegin = Ret.second + 1;
        string BlockName = trim(Ret.first);

        Ret = FindNextBracketPair(Zeilen , BlockBegin );
        cfgfile.insert( make_pair(BlockName , Block(BlockName , vector<string> (Zeilen.begin() + BlockBegin , Zeilen.begin() + Ret.second) )) );
    }
    return;
}

#ifndef HAVE_ALPHASORT
inline int alphasort(const void* d1, const void* d2)
{//MINGW doesn't define this
    return(strcmp((*static_cast<struct dirent **>(const_cast<void*>(d1) ) )->d_name, (*static_cast<struct dirent **>(const_cast<void*>(d2) ) )->d_name));
}
#endif

#ifndef HAVE_SCANDIR
#include <unistd.h>
#include <string.h>
inline int scandir(const char* dirname, struct dirent *(*namelist[]), int (*select)(const struct dirent *), int (*dcomp)(const void* , const void*))
{//from here: http://www.tina-vision.net/doxygen/tina-libs/html/fileUtil__name_8c-source.html
    DIR *dirp;
    if ((dirp = opendir(dirname)) == NULL)
        return -1;
    // 256 should be big enough for a spool directory of any size * big enough for 128 jobs anyway.
    if ((*namelist = static_cast<struct dirent **> ( calloc(256, sizeof(struct dirent *) ) ) ) == NULL)
    {
        closedir(dirp);
        return -1;
    }
    struct dirent *res;
    const int tdirsize = sizeof(struct dirent);
    if ((res = static_cast<struct dirent *> ( malloc(tdirsize + 256) ) )== NULL )
    {
        closedir(dirp);
        return -1;
    }
    int i = 0;
    while (
#ifdef HAVE_READDIR_R
        readdir_r(dirp, res)
#else
        ( res = readdir(dirp) )//This might be not thread-safe...
#endif
        != NULL)
    {
        if (select(res))
        {
            if (((*namelist)[i] = static_cast<struct dirent *> (malloc(tdirsize + 256)) ) == NULL)
            {
                closedir(dirp);
                return -1;
            }
            memcpy((*namelist)[i], res, sizeof(res)+256);
            i++;
        }
    }

    if (dcomp != NULL)
        qsort(reinterpret_cast <char *> (&( (*namelist)[0] ) ), i, sizeof(struct dirent *), dcomp);
    free(res);
    closedir(dirp);
    return i;
}
#endif

static int
IsFile (const struct dirent *used)
{
    string A(used->d_name);
    if ( ( A == "." ) || ( A == ".." ) ) //Don't use these
        return 0;
    return 1;
}

RegistryDB::RegistryDB( const string& arg)
{
    Init(arg);
    return;
}

int RegistryDB::Init( const string& cfgDir )
{
    //The following Directory - Lister is taken directly from the glibc - Documentation
    struct dirent **eps = NULL;
    string* configs = NULL;
    int n = scandir( cfgDir.c_str() , &eps , IsFile , alphasort );
    //cout<<"Content of Configuration Directory  " << cfgDir<<" : "<<endl;
    if ( n > 0)
    {
        configs = new string[n];
        for (unsigned int cnt = 0 ; cnt < static_cast<unsigned int>(n) ; ++cnt )
        {
            configs[cnt] = eps[cnt]->d_name;
        }
    }
    else
    {
        free (eps);//Never ever insert delete [] here, scandir allocates the array with malloc so use free
        throw(Registry_Exception("Configuration Directory not found or empty!!"));
    }
    int success = chdir(cfgDir.c_str());//change working Directory to the Config Directory
    if(success != 0)
    {
        delete [] configs;
        for(unsigned int k = 0; k < static_cast<unsigned int>(n); ++k)
            free(eps[k]);
        free (eps); //Never ever insert delete [] here, scandir allocates the array with malloc so use free
        cout<<"Error while changing to directory "<<cfgDir<<" with error code: "<<success<<endl;
        throw(Registry_Exception(std::string("Error while changing directory: ") + cfgDir));
    }
    for ( unsigned int k = 0; k < static_cast<unsigned int>(n); ++k )
        Reg.insert( make_pair(configs[k] , cfile(configs[k] ) ) );
    //FIXME
    int errcode = chdir("..");//Get Back to the previous Directory
    delete [] configs;
    for(unsigned int k = 0; k < static_cast<unsigned int>(n); ++k)
        free(eps[k]);
    free (eps); //Never ever insert delete [] here, scandir allocates the array with malloc so use free
    if(errcode != 0)
    {
      throw(Registry_Exception(strerror(errcode)));
    }
    return 0;
}
