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
#ifndef FILEVECTOR_H
#define FILEVECTOR_H
#include <cstring>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <valarray>
#include <inttypes.h>
#include <errno.h>
#include <iterator>
#include <stdexcept>

#if ((__GLIBC__ > 1) && __GLIBC_MINOR__ > 1)
//we define the USELFS macro to signify that we are on a gnu libc. LFS support is available since 2.2
#define USELFS
//request LFS support
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#endif
/**
@file filevector.h This File contains the file-vector. It's supposed to provide an interface similar to std::vector, but store its content in a file.
The file that the FileVector uses stays open permanently. So take care not to overflow your OSs limit on the maximum number of open files.
*/

template <class T>
class FileVector;

namespace FileVector_Private
{
/**
@class FileVectorException This is the base class for all exceptions the filevector can throw. It is derived from std::runtime_error.
*/
class FileVectorException : public std::runtime_error
{
public:
    FileVectorException(const std::string& a) : std::runtime_error(a){}
};

/**
@class FileVectorReadError This exception signifies errors while reading.
*/
class FileVectorReadError : public FileVectorException
{
public:
    FileVectorReadError(const std::string& a) : FileVectorException(a) {}
};

/**
@class FileVectorWriteError This exception signifies errors while writing data.
*/
class FileVectorWriteError : public FileVectorException
{
public:
    FileVectorWriteError(const std::string& a) : FileVectorException(a) {}
};

/**
@class FileVectorOpenError This exception signifies errors that happened during acquiring the file handle.
*/
class FileVectorOpenError : public FileVectorException
{
public:
    FileVectorOpenError(const std::string& a) : FileVectorException(a) {}
};

/**
@class Lock This class represents the Lock that a FileVector holds onto a certain file. The FileVector will not do changes to a file that is used by another filevector
*/
class Lock
{
public:
    inline Lock(const std::string& fi) : filename(fi+ std::string(".lock"))
    {
        FILE* lockptr = fopen(filename.c_str(), "wb+");
        fclose(lockptr);
    }
    inline ~Lock()
    {
        if (remove(filename.c_str()) != 0)
            std::cerr<<"[Lock] An error happened: "<<errno<<std::endl;
    }
private:
    const std::string filename;///< the filename that is used for the lockfile
};

static inline bool islocked(const std::string& fi)
{
    FILE* test = fopen((fi + std::string(".lock")).c_str(), "rb");
    if (test != NULL)
        fclose(test);
    return test != NULL;
}

template <class T>
class FileVector_Iterator : std::iterator<std::input_iterator_tag, T>
{
public:
    inline FileVector_Iterator(FileVector<T>& fv) : file(fv.dump), pos(ftello(file)), isbinary(fv.isbinary)
    {
        if(ferror(file) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file for iteration."));
        if (isbinary)
        {
            if (fread(&elem, sizeof(T), 1, file) != 1) throw FileVectorReadError("[FileVector] couldn't read element!");
        }
        else
        {
            std::string temp;
            unsigned char t;
            while ((t = static_cast<char>(getc(file))) != '\n')//getc could return EOF too...
            {
                temp += t;
            }
            std::istringstream i(temp);
            if (!(i >> elem) )
            {
                std::string errmsg("[FileVector] couldn't decipher text: ");
                errmsg += temp;
                throw FileVectorException(errmsg);
            }
        }
    }
    /**
    Having multiple iterators work on the same filevector is possible, but the behaviour is undefined.
    */
    inline FileVector_Iterator(const FileVector_Iterator& mit) : file(mit.file), pos(mit.pos), isbinary(mit.isbinary), elem(mit.elem) {}
    /**
    Advance one position forward, prefix version
    */
    inline FileVector_Iterator& operator++()// equals ++it
    {
        if (isbinary)
        {
            size_t readelems = fread(&elem, sizeof(T), 1, file);
            if((readelems != 1) && (feof(file) == 0)) throw(FileVectorReadError("[FileVector] Error while advancing iterator one element!"));
        }
        else
        {
            std::string temp;
            unsigned char t;
            while (((t = static_cast<unsigned char>(getc(file))) != '\n') && (feof(file) == 0))
            {
                temp += t;
            }
            if (!temp.empty())
            {
                std::istringstream i(temp);
                if (!(i >> elem) )
                {
                    std::string errmsg("[FileVector] couldn't decipher text: ");
                    errmsg += temp;
                    throw FileVectorException(errmsg);
                }
            }
        }
        ++pos;
        return *this;
    }
    /**
    Advance one position forward, postfix version
    */
    inline FileVector_Iterator operator++(int) //equals it++
    {
        FileVector_Iterator retval(*this);//save state for later
        if (isbinary)
        {
            fread(&elem, sizeof(T), 1, file);
        }
        else
        {
            std::string temp;
            unsigned char t;
            while (((t = getc(file)) != '\n') && (feof(file) == 0))
            {
                temp += t;
            }
            if (!temp.empty())
            {
                std::istringstream i(temp);
                if (!(i >> elem) )
                {
                    std::string errmsg("[FileVector] couldn't decipher text: ");
                    errmsg += temp;
                    throw FileVectorException(errmsg);
                }
            }
        }
        pos++;
        return retval;
    }
    /**
    Checks if two iterators are equal with respect to the file they're pointing to and if they point to the same offset
    @param rhs a file vector iterator
    @return true if the iterators point to the same file and to the same offset, else false
    */
    inline bool operator==(const FileVector_Iterator& rhs)
    {
        return (file == rhs.file) && (pos == rhs.pos);
    }
    /**
    Checks if two iterators are unequal with respect to the file they're pointing to and if they point to the same offset
    @param rhs a file vector iterator
    @return false if the iterators point to the same file or to the same offset, else false
    */
    inline bool operator!=(const FileVector_Iterator& rhs)
    {
        return !((file == rhs.file) && (pos == rhs.pos));
    }
    /**
    dereference the iterator
    @return the element the iterator is pointing to
    */
    inline T operator*()
    {
        return elem;
    }
private:
    friend class FileVector<T>;
    inline FileVector_Iterator(FILE*& fp, uint32_t p, bool isb) : file(fp), pos(p), isbinary(isb)
    {
    }
    FILE*& file;///< a reference to the file pointer of the filevector
    uint32_t pos;///< the offset in units of the stored items
    const bool isbinary;
    T elem;
};
};

template<class T>
class FileVector
{
public:
    typedef T value_type;///< this typedef provides the type of the contained data. And makes this class STL-compatible
    typedef off_t size_type;///< this typedef provides the type of the numbers. And makes this class STL-compatible
    typedef T& reference;///< this typedef provides the type of the reference. And makes this class STL-compatible
    typedef const T& const_reference;
    friend class FileVector_Private::FileVector_Iterator<T>;///< Iterators are our friends!
    typedef FileVector_Private::FileVector_Iterator<T> const_iterator;
    /**
    Returns a FileVector_Iterator to the beginning of the file. Note that for that the internal file pointer gets set to the beginning.
    @return an iterator that points to the beginning of the file
    */
    inline const_iterator begin();
    /**
    Returns a FileVector_Iterator past the end of the file
    @return an iterator past the end of the file
    */
    inline const_iterator end();
    /**
    Default Constructor.
    @param isb this signifies if the filevector uses text files or binary files
    */
    inline FileVector(bool isb);
    /**
    Constructor that directly opens a file with the specified mode. If the File exists the file will be overwritten
    @param filename the FileName
    @param isbinary signifies whether the file is binary or not
    */
    inline FileVector(const std::string& filename, bool isbinary);
    /**
    Append an element.
    This appends an Element to the end of the file.
    @param arg the element to store
    */
    inline void push_back(const T& arg);
    /**
    Size of the FileVector
    Returns the numer of values that the file stores
    @return the number of members
    */
    inline size_type size() const throw();
    /**
    Access operator to the Elements stored in a file.
    @param n which element to get
    @return the element at position n
    */
    inline value_type operator[] ( size_type n );
    /**
    Open a file
    This function opens a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    */
    inline void open ( const char *const filename);
    /**
    Open a file
    This functions open a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    */
    inline void open ( const std::string& filename);
    /**
    Read from a file.
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read(const std::string& filename);
    /**
    Read from a file
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read ( const char *const filename);
    /**
    Load existing file.
    This functions makes the content of the files accessible, but also allows appending new data to the file.
    @param filename the filename of the file
    */
    inline void load(const std::string& filename);
    /**
    Load existing file.
    This functions makes the content of the files accessible, but also allows appending new data to the file.
    @param filename the filename of the file
    */
    inline void load ( const char *const filename);
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create(const std::string& filename);
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create( const char *const filename);
    /**
    The destructor. It deletes the lock-file
    */
    inline virtual ~FileVector();
    inline void sync(){fflush(dump);}
private:
    inline std::string toString(const T& Value) const;
    FILE* dump;///< the pointer to the file
    size_type size_;///< the length of the filevector
    const bool isbinary;///< a bool to signify whether a file is binary or not
    std::string filename;///< the filename of the file
    FileVector_Private::Lock* lock;///< the lockfile
};

template<class T>
typename FileVector<T>::const_iterator FileVector<T>::begin()
{
    rewind(dump);
    return const_iterator(*this);
}

template<class T>
typename FileVector<T>::const_iterator FileVector<T>::end()
{
    return const_iterator(dump, size_, isbinary);
}

template<class T>
FileVector<T>::~FileVector()
{
    if (dump != NULL)
    {
        fclose(dump);
        delete lock;
    }
}

template<class T>
std::string FileVector<T>::toString(const T& Value) const
{
    std::stringstream ss;
    ss << Value;
    return ss.str();
}

template<class T>
void FileVector<T>::create( const std::string& filename)
{
    this->create(filename.c_str());
}

template<class T>
void FileVector<T>::create( const char *const fi)
{
    filename = fi;
    dump = fopen(fi, "r");
    if (dump == NULL)
        std::cerr<<"File doesn't exist! Creating"<<std::endl;
    else
        fclose(dump);
    if (isbinary)
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
    else
        dump = fopen(fi, "w+");
    lock = new FileVector_Private::Lock(fi);
}

template<class T>
void FileVector<T>::load( const std::string& filename)
{
    this->load(filename.c_str());
}

template<class T>
void FileVector<T>::load( const char *const fi)
{
    filename = fi;
    if (FileVector_Private::islocked(fi))
    {
        std::cerr<<"[FileVector] Can't open file which is locked by another filevector! Aborting"<<std::endl;
        exit(-1);
    }
    dump = fopen(fi, "r");
    if (dump == NULL)
        std::cerr<<"[FileVector] File doesn't exist! Creating"<<std::endl;
    else
        fclose(dump);
    if (isbinary)
    {
        dump = fopen(fi, "ab+");//Create an empty file for both reading and writing. This is a binary file
        if (dump == NULL)
        {
            throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
        }
        off_t pos = ftello(dump);
        if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
        fseeko(dump, 0, SEEK_END);
        size_ = ftello(dump)/sizeof(T);
        if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
        fseeko(dump, pos, SEEK_SET);
    }
    else
    {
        dump = fopen(fi, "a+");
        std::cerr<<"[FileVector] ERROR! Loading of text Files not even implemented!"<<std::endl;
        exit(-1);
    }
    lock = new FileVector_Private::Lock(fi);
    return;
}

template<class T>
void FileVector<T>::read( const std::string& filename)
{
    this->read(filename.c_str());
}

template<class T>
void FileVector<T>::read( const char *const fi)
{
    filename = fi;
    if (isbinary)
    {
        dump = fopen(fi, "rb");
        if (dump != NULL)
        {
            fseeko(dump, 0, SEEK_END);
            size_ = ftello(dump)/sizeof(T);
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error reading file: ") + fi);
            fseeko(dump, 0, SEEK_SET);
        }
        else
        {
            std::string errmsg("[FileVector] File not found: ");
            errmsg += filename;
            throw FileVector_Private::FileVectorOpenError(errmsg);
        }
    }
    else
    {
        dump = fopen(fi, "r");
        std::cerr<<"[FileVector] ERROR! Reading of text Files not even implemented!"<<std::endl;
    }
    if (dump == NULL)
    {
        std::string errmsg("[FileVector] File not found: ");
        errmsg += filename;
        throw FileVector_Private::FileVectorOpenError(errmsg);
    }
    return;
}

template<class T>
void FileVector<T>::open( const std::string& filename)
{
    this->open(filename.c_str());
}

template<class T>
void FileVector<T>::open( const char *const fi)
{
    filename = fi;
    if (FileVector_Private::islocked(fi))
    {
        std::cerr<<"[FileVector] Can't open file which is locked by another filevector! Aborting"<<std::endl;
        exit(-1);
    }
    dump = fopen(fi, "rb");
    if (dump != NULL)
    {
        std::cerr<<"[FileVector] WARNING A FILE WITH THAT NAME: "<<filename<<" ALREADY EXISTS! OVERWRITING!"<<std::endl;
        fclose(dump);
    }
    if (isbinary)
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
    else
        dump = fopen(fi, "w+");
    if (dump == NULL)
    {
        throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
    }
    lock = new FileVector_Private::Lock(filename);
}

template<class T>
typename FileVector<T>::value_type FileVector<T>::operator[](const size_type n)
{
    value_type retval;
    if (isbinary)
    {
        fseeko(dump, static_cast<off_t>(n) * sizeof(T), SEEK_SET);
        size_t readelems = fread(&retval, sizeof(value_type), 1, dump);
        if (readelems != 1) throw FileVector_Private::FileVectorReadError("Error reading from file!");
    }
    else
    {
        //untested code-path
        off_t pos = ftello(dump);
        fseeko(dump, 0, SEEK_END);
        off_t length = ftello(dump);
        fseeko(dump, pos, SEEK_SET);
        char* buf = new char[length + 1];
        if(fread(buf, length, 1, dump) != 1) throw FileVector_Private::FileVectorReadError("[FileVector] couldn't read element!");
        buf[length] = 0;
        char* pch = strtok(buf, "\n");
        size_type len = 0;
        while ((len < n) && pch != NULL)
        {
            pch = strtok(NULL,"\n");
            ++len;
        }
        std::string key(strtok(NULL, "\n"));
        std::stringstream ss;
        ss << key;
        ss >> retval;
    }
    return retval;
}

template<class T>
typename FileVector<T>::size_type FileVector<T>::size() const throw()
{
    return size_;
}

template<class T>
void FileVector<T>::push_back(const T& arg)
{
    value_type arg2(arg);
    if (isbinary)
    {
        fseeko(dump, 0, SEEK_END);
        size_t writtenelems = fwrite(reinterpret_cast<char*>(&arg2), sizeof(T), 1, dump);
        if (writtenelems != 1)
            throw FileVector_Private::FileVectorWriteError("[FileVector] Error while writing data to file. File is perhaps read-only?");
    }
    else
    {
        std::string temp(toString(arg));
        temp += '\n';
        fseeko(dump, 0, SEEK_END);
        if(fwrite(temp.c_str(), temp.size(), 1, dump) != 1) throw FileVector_Private::FileVectorWriteError("[FileVector] Error while writing data to file. File is perhaps read-only?");
    }
    ++size_;
    return;
}

template<class T>
FileVector<T>::FileVector(bool isb = true) : dump(NULL), size_(0), isbinary(isb), lock(NULL)
{
}

template<class T>
FileVector<T>::FileVector(const std::string& filename, bool isb = true) : size_(0), dump(NULL), isbinary(isb)
{
    this->open(filename);
    return;
}

template<class T>
class FileVector<std::valarray<T> >
{
public:
    typedef std::valarray<T> value_type;///< this typedef provides the type of the contained data. And makes this function STL-compatible
    typedef off_t size_type;///< this typedef provides the type of the numbers. And makes this function STL-compatible
    typedef std::valarray<T>& reference;///< this typedef provides the type of the reference. And makes this function STL-compatible
    typedef const std::valarray<T>& const_reference;
    /**
    Default Constructor. Does nothing.
    */
    inline FileVector();
    inline virtual ~FileVector();
    /**
    Constructor that directly opens a file with the specified mode. If the File exists the file will be overwritten
    @param filename the FileName
    @param l the length of a single entry
    @param isbinary signifies whether the file is binary or not
    */
    inline FileVector(const std::string& filename, uint32_t l, bool isbinary);
    /**
    Appends an Element to the end of the file.
    @param arg the element to store
    */
    inline void push_back(const std::valarray<T>& arg);
    /**
    Returns the number of values that the file stores
    @return the number of members
    */
    inline size_type size() const throw();
    /**
    Access operator to the Elements stored in a file.
    @param n which element to get
    @return the element at position n
    */
    inline value_type operator[] ( size_type n );
    /**
    Provide array-like access to the contents of a bin function
    @param i the function stored in which bin
    @param k which point of the function
    */
    inline T operator()(unsigned int i, unsigned int k);
    /**
    This opens a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    @param l the length of a stored entry
    */
    inline void open ( const char *const filename, uint32_t l);
    /**
    This opens a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    @param l the length of a stored entry
    */
    inline void open ( const std::string& filename, uint32_t l);
    /**
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read(const std::string& filename);
    /**
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read ( const char *const filename);
    /**
    This function makes the content of the files accessible, but also allows appending new data to the file. If the file is not found it works as open().
    @param filename the filename of the file.
    @param slices the number of slices of the function, in case the file is not found. The data in the file has precedence.
    */
    inline void load(const std::string& filename, unsigned int slices = 0);
    /**
    This function makes the content of the files accessible, but also allows appending new data to the file. If the file is not found it works as open().
    @param filename the filename of the file.
    @param slices the number of slices of the function, in case the file is not found. The data in the file has precedence.
    */
    inline void load ( const char *const filename, unsigned int slices = 0);
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create(const std::string& filename) throw();
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create( const char *const filename) throw();
    inline void sync(){fflush(dump);}
private:
    /**
    An internal helper function for conversion of elements to strings
    */
    inline std::string toString(const value_type& Value) const;
    FILE* dump;///< a pointer to the FILE structure
    size_type size_;///< how many valarrays did we store
    const bool isbinary;///< a bool to determine whether the file is a binary file or not
    std::string filename;///< the bins are stored under this filename
    uint32_t len;///< this is the length of a single valarray
    static const size_type offs = sizeof(uint32_t);///< where starts the beginning of the file
    FileVector_Private::Lock* lock;
};

template<class T>
std::string FileVector<std::valarray<T> >::toString(const value_type& value) const
{
    std::stringstream ss;
    ss<<"{";
    for (unsigned int k = 0; k < value.size(); ++k)
        ss << value[k] <<",";
    ss<<"}";
    return ss.str();
}

template<class T>
void FileVector<std::valarray<T> >::create( const std::string& filename) throw()
{
    this->create(filename.c_str());
}

template<class T>
void FileVector<std::valarray<T> >::create( const char *const fi) throw()
{
    filename = fi;
#ifdef USELFS
    dump = fopen64(fi, "r");
#else
    dump = fopen(fi, "r");
#endif
    if (dump == NULL)
        std::cerr<<"[FileVector] File doesn't exist! Creating"<<std::endl;
    else
        fclose(dump);
    if (isbinary)
    {
#ifdef USELFS
        dump = fopen64(fi, "wb+");
#else
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
    }
    else
    {
#ifdef USELFS
        dump = fopen64(fi, "w+");
#else
        dump = fopen(fi, "w+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
    }
}

template<class T>
void FileVector<std::valarray<T> >::load( const std::string& filename, unsigned int slices)
{
    this->load(filename.c_str(), slices);
}

template<class T>
void FileVector<std::valarray<T> >::load( const char *const fi, unsigned int slices)
{
    filename = fi;
    bool fileexists = true;
#ifdef USELFS
    dump = fopen64(fi, "r");
#else
    dump = fopen(fi, "r");
#endif
    if (dump == NULL)
    {
        std::cerr<<"[FileVector] File doesn't exist! Creating"<<std::endl;
        fileexists = false;
    }
    else
        fclose(dump);
    if (FileVector_Private::islocked(fi))
    {
        std::cerr<<"[FileVector] Can't open file which is locked by another filevector! Aborting"<<std::endl;
        exit(-1);
    }
    if (isbinary)
    {
#ifdef USELFS
        dump = fopen64(fi, "ab+");
#else
        dump = fopen(fi, "ab+");//Create an empty file for both reading and writing to the end of the file(appending). This is a binary file
#endif
        if (dump == NULL)
        {
            throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
        }
        if (fileexists)
        {
            if (fread(&len, sizeof(uint32_t), 1, dump) != 1)
                throw(FileVector_Private::FileVectorReadError("[FileVector] Error while reading slice number!"));
            if ((slices != 0) && (len != slices))
            {
                std::cerr<<"[FileVector] Warning: file and expected slice number don't match!"<<std::endl;
            }
            
            
            
            uint binsize = sizeof(T) * len;
            off64_t pos = ftello(dump);//pos of the beginning of the datastream
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error loading file: ") + fi);
            fseeko(dump, 0, SEEK_END);
            size_ = static_cast<size_type>((ftello64(dump) - offs)/static_cast<off64_t>(binsize));//retrieve the number of stored bins
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error loading file: ") + fi);
            fseeko(dump, static_cast<off_t>(pos), SEEK_SET);
        }
        else
        {
            len = slices;
            fwrite(&len, sizeof(uint32_t), 1, dump);
            size_ = 0;
        }
    }
    else
    {
        dump = fopen(fi, "a+");
        std::cerr<<"[FileVector] ERROR! Loading of text files not supported!"<<std::endl;
        exit(-1);
    }
    lock = new FileVector_Private::Lock(fi);
    return;
}

template<class T>
void FileVector<std::valarray<T> >::read( const std::string& filename)
{
    this->read(filename.c_str());
}

template<class T>
void FileVector<std::valarray<T> >::read( const char *const fi)
{
    filename = fi;
    if (isbinary)
    {
#ifdef USELFS
        dump = fopen64(fi, "rb");
#else
        dump = fopen(fi, "rb");
#endif
        if (dump != NULL)
        {
            fread(&len, sizeof(uint32_t), 1, dump);
            off_t pos = ftello(dump);//pos of the beginning of the datastream
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error reading file: ") + fi);
            fseeko(dump, 0, SEEK_END);
            size_ = (ftello(dump) - offs)/(sizeof(T)*len);
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error reading file: ") + fi);
            fseek(dump, pos, SEEK_SET);
        }
        else
        {
            std::string errmsg("[FileVector] File not found: ");
            errmsg += std::string(fi);
            throw FileVector_Private::FileVectorOpenError(errmsg);
        }
    }
    else
    {
#ifdef USELFS
        dump = fopen64(fi, "r");
#else
        dump = fopen(fi, "r");
#endif
    }
    if (dump == NULL)
    {
        std::string errmsg("[FileVector] File not found: ");
        errmsg += std::string(fi);
        throw FileVector_Private::FileVectorOpenError(errmsg);
    }
    return;
}

template<class T>
void FileVector<std::valarray<T> >::open( const std::string& filename, uint32_t l)
{
    this->open(filename.c_str(), l);
}

template<class T>
void FileVector<std::valarray<T> >::open( const char *const fi, uint32_t l)
{
    filename = fi;
    len = l;
    dump = fopen(fi, "rb");
    if (dump != NULL)
    {
        std::cerr<<"[FileVector] WARNING A FILE WITH THAT NAME: "<<filename<<" ALREADY EXISTS! OVERWRITING!"<<std::endl;
        fclose(dump);
    }
    if (isbinary)
    {
#ifdef USELFS
        dump = fopen64(fi, "wb+");
#else
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
    }
    else
    {
#ifdef USELFS
        dump = fopen64(fi, "w+");
#else
        dump = fopen(fi, "w+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
    }
    if (dump == NULL)
    {
        throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
    }
    else
    {
        if(fwrite(&len, sizeof(uint32_t), 1, dump) != 1) throw FileVector_Private::FileVectorWriteError("[FileVector] Error writing element to file!");
    }
}

template<class T>
T FileVector<std::valarray<T> >::operator()(unsigned int i, unsigned int k)
{
    T retval;
    if (isbinary)
    {
        off_t offset = static_cast<off_t>(offs) + (static_cast<off_t>(i) * static_cast<off_t>(len)  + static_cast<off_t>(k)) * static_cast<off_t>(sizeof(T));
        fseeko(dump, offset, SEEK_SET);
        size_t readelems = fread(&retval, sizeof(T), 1, dump);
        if (readelems != 1)
        {
            std::cerr<<"[FileVector] Error reading file! "<<readelems<< " "<<len<<std::endl;
        }
    }
    else
    {
        //untested code-path
        off_t pos = ftello(dump);
        fseeko(dump, 0, SEEK_END);
        fseeko(dump, pos, SEEK_SET);
        std::cerr<<"[FileVector] random access operators in text file based valarrays not supported!"<<std::endl;
        exit(-1);
        /*        char* buf = new char[length + 1];
                fread(buf, length, 1, dump);
                buf[length] = 0;
                char* pch = strtok(buf, "\n");
                unsigned int len = 0;
                while ((len < n) && pch != NULL)
                {
                    pch = strtok(NULL,"\n");
                    ++len;
                }
                std::string key(strtok(NULL, "\n"));
                std::stringstream ss;
                ss << key;
                ss >>retval;*/
    }
    return retval;
}

template<class T>
std::valarray<T> FileVector<std::valarray<T> >::operator[](const size_type n)
{
    value_type retval(len);
    if (isbinary)
    {
        off_t offset = static_cast<off_t>(offs) + static_cast<off_t>(n) * static_cast<off_t>(len) * static_cast<off_t>(sizeof(T));
        fseeko(dump, offset, SEEK_SET);
        T *const tempmem(new T[len]);
        size_t readelems = fread(tempmem, sizeof(T), len, dump);
        if (readelems != len)
        {
            std::stringstream s("[FileVector] Error reading file! ");
            s<< readelems << " " << len;
            throw FileVector_Private::FileVectorReadError(s.str());
        }
        for (unsigned int k = 0; k < len; ++k)
        {
            retval[k] = tempmem[k];
        }
        delete [] tempmem;
    }
    else
    {
        //untested code-path
        off_t pos = ftello(dump);
        fseeko(dump, 0, SEEK_END);
//        long int length = ftell(dump);
        fseeko(dump, pos, SEEK_SET);
        std::cerr<<"[FileVector] random access operators in text file based valarrays not supported!"<<std::endl;
        exit(-1);
        /*        char* buf = new char[length + 1];
                fread(buf, length, 1, dump);
                buf[length] = 0;
                char* pch = strtok(buf, "\n");
                unsigned int len = 0;
                while ((len < n) && pch != NULL)
                {
                    pch = strtok(NULL,"\n");
                    ++len;
                }
                std::string key(strtok(NULL, "\n"));
                std::stringstream ss;
                ss << key;
                ss >>retval;*/
    }
    return retval;
}

template<class T>
typename FileVector<std::valarray<T> >::size_type FileVector<std::valarray<T> >::size() const throw()
{
    return size_;
}

template<class T>
void FileVector<std::valarray<T> >::push_back(const value_type& arg)
{
    if (isbinary)
    {
        T *const tempmem(new T[arg.size()]);
        for (unsigned int k = 0; k < arg.size(); ++k)
            tempmem[k] = arg[k];
        fseeko(dump, 0, SEEK_END);
        size_t writtenelems = fwrite(reinterpret_cast<char *const>(tempmem), sizeof(T), arg.size(), dump);
        delete [] tempmem;
        if (writtenelems != arg.size())
            throw FileVector_Private::FileVectorWriteError("[FileVector] Error while writing data to file. File is perhaps read-only?");
    }
    else
    {
        std::string temp(toString(arg));
        temp += '\n';
        fseeko(dump, 0, SEEK_END);
        fwrite(temp.c_str(), temp.size(), 1, dump);
    }
    ++size_;
    return;
}

template<class T>
FileVector<std::valarray<T> >::FileVector() : dump(NULL), size_(0), isbinary(true), lock(NULL)
{
}

template<class T>
FileVector<std::valarray<T> >::~FileVector()
{
    if (dump != NULL)
    {
        fclose(dump);
        delete lock;
    }
}

template<class T>
FileVector<std::valarray<T> >::FileVector(const std::string& filename, uint32_t l, bool isb = true) : size_(0), dump(NULL), isbinary(isb), len(l), lock(NULL)
{
    this->open(filename, l);
    return;
}

/**
 * Another instance of the filevector. This time for Vectorfunction. Hence every bin contains several functions.
 * Note that we removed the text mode...
 * */
template<class T>
class FileVector<std::valarray<std::valarray<T> > >
{
public:
    typedef std::valarray<T> Function;
    typedef std::valarray<Function> value_type;///< this typedef provides the type of the contained data. And makes this function STL-compatible
    typedef off_t size_type;///< this typedef provides the type of the numbers. And makes this function STL-compatible
    typedef std::valarray<Function>& reference;///< this typedef provides the type of the reference. And makes this function STL-compatible
    typedef const std::valarray<Function>& const_reference;
    /**
    Default Constructor. Does nothing.
    */
    inline FileVector();
    inline virtual ~FileVector();
    /**
    Constructor that directly opens a file with the specified mode. If the File exists the file will be overwritten
    @param filename the FileName
    @param l the length of a single entry
    */
    inline FileVector(const std::string& filename, uint32_t l);
    /**
    Appends an Element to the end of the file.
    @param arg the element to store
    */
    inline void push_back(const_reference arg);
    /**
    Returns the number of values that the file stores
    @return the number of members
    */
    inline size_type size() const throw();
    /**
    Access operator to the Elements stored in a file.
    @param n which element to get
    @return the element at position n
    */
    inline value_type operator[] ( size_type n );
    /** 
     * an access operator to access the i'th bin and the k'th function index
     * @param i which bin you access
     * @param k the function index that you want to index
     * @return the k'th function(a full valarray)  of the i'th bin
    */
    inline Function operator()(unsigned int i, unsigned int k);
    /**
    This opens a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    @param l the length of a stored entry
    */
    inline void open ( const char *const filename, uint32_t l, uint32_t nrofFunctions);
    /**
    This opens a file for reading and writing. Note, that if a file with the same name exists the content of the file is destroyed and overwritten
    @param filename the filename
    @param l the length of a stored entry
    */
    inline void open ( const std::string& filename, uint32_t l, uint32_t nrofFunctions);
    /**
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read(const std::string& filename);
    /**
    This functions opens a file only for reading.
    @param filename the filename
    */
    inline void read ( const char *const filename);
    /**
    This function makes the content of the files accessible, but also allows appending new data to the file. If the file is not found it works as open().
    @param filename the filename of the file.
    @param slices the number of slices of the function, in case the file is not found. The data in the file has precedence.
    */
    inline void load(const std::string& filename, unsigned int slices = 0, unsigned int nrofFunctions = 0);
    /**
    This function makes the content of the files accessible, but also allows appending new data to the file. If the file is not found it works as open().
    @param filename the filename of the file.
    @param slices the number of slices of the function, in case the file is not found. The data in the file has precedence.
    */
    inline void load ( const char *const filename, unsigned int slices = 0, unsigned int nrofFunctions = 0);
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create(const std::string& filename) throw();
    /**
    Try to create a new file with the given filename
    @param filename the filename of the file
    */
    inline void create( const char *const filename) throw();
    inline void sync(){fflush(dump);}
private:
    FILE* dump;///< a pointer to the FILE structure
    size_type size_;///< how many valarrays did we store
    std::string filename;///< the bins are stored under this filename
    uint32_t len;///< this is the length of a single valarray
    uint32_t functions;///< this is the number of functions that a single function makes up
    uint32_t binsize;///< the length of a bin in bytes
    static const size_type offs = 2*sizeof(uint32_t);///< where starts the beginning of the file
    FileVector_Private::Lock* lock;
    off64_t curpos;///<< here we store the current position. That way we get sometmes rid of an expensive sync that is incurred by fseek
};

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::create( const std::string& filename) throw()
{
    this->create(filename.c_str());
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::create( const char *const fi) throw()
{
    filename = fi;
#ifdef USELFS
    dump = fopen64(fi, "r");
#else
    dump = fopen(fi, "r");
#endif
    if (dump == NULL)
        std::cerr<<"[FileVector] File doesn't exist! Creating"<<std::endl;
    else
        fclose(dump);
#ifdef USELFS
        dump = fopen64(fi, "wb+");
#else
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::load( const std::string& filename, unsigned int slices, unsigned int nrofFunctions)
{
    this->load(filename.c_str(), slices, nrofFunctions);
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::load( const char *const fi, unsigned int slices, unsigned int nrofFunctions)
{
    filename = fi;
    bool fileexists = true;
#ifdef USELFS
    dump = fopen64(fi, "r");
#else
    dump = fopen(fi, "r");
#endif    
    if (dump == NULL)
    {
        std::cerr<<"[FileVector] File doesn't exist! Creating"<<std::endl;
        fileexists = false;
    }
    else
        fclose(dump);
    if (FileVector_Private::islocked(fi))
    {
        std::cerr<<"[FileVector] Can't open file which is locked by another filevector! Aborting"<<std::endl;
        exit(-1);
    }
#ifdef USELFS
        dump = fopen64(fi, "ab+");
#else
        dump = fopen(fi, "ab+");//Create an empty file for both reading and writing to the end of the file(appending). This is a binary file
#endif
        if (dump == NULL)
        {
            throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
        }
        if (fileexists)
        {
            if (fread(&len, sizeof(uint32_t), 1, dump) != 1)
                throw(FileVector_Private::FileVectorReadError("[FileVector] Error while reading slice number!"));
	    if (fread(&functions, sizeof(uint32_t), 1, dump) != 1)
                throw(FileVector_Private::FileVectorReadError("[FileVector] Error while reading number of functions!"));
            if ((slices != 0) && (len != slices))
            {
                std::cerr<<"[FileVector] Warning: file and expected slice number don't match!"<<std::endl;
            }
            if ((nrofFunctions != 0) && (functions != nrofFunctions))
            {
                std::cerr<<"[FileVector] Warning: file and expected function number don't match!"<<std::endl;
            }
            binsize = sizeof(T) * len * functions;
            off64_t pos = ftello(dump);//pos of the beginning of the datastream
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error loading file: ") + fi);
            fseeko(dump, 0, SEEK_END);
            size_ = static_cast<size_type>((ftello64(dump) - offs)/static_cast<off64_t>(binsize));//retrieve the number of stored bins
            if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error loading file: ") + fi);
            fseeko(dump, static_cast<off_t>(pos), SEEK_SET);
	    curpos = pos;
        }
        else
        {
            len = slices;
	    functions = nrofFunctions;
            fwrite(&len, sizeof(uint32_t), 1, dump);
	    fwrite(&functions, sizeof(uint32_t), 1, dump);
            size_ = 0;
	    binsize = sizeof(T) * len * functions;
        }
    lock = new FileVector_Private::Lock(fi);
    return;
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::read( const std::string& filename)
{
    this->read(filename.c_str());
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::read( const char *const fi)
{
    filename = fi;
#ifdef USELFS
    dump = fopen64(fi, "rb");
#else
    dump = fopen(fi, "rb");
#endif
    if (dump != NULL)
    {
       fread(&len, sizeof(uint32_t), 1, dump);
       fread(&functions, sizeof(uint32_t), 1, dump);
       off_t pos = ftello(dump);//pos of the beginning of the datastream
       if (pos == (off_t) -1) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error positioning file pointer ") + fi);
       if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error reading file: ") + fi);
       int seekerr = fseeko(dump, 0, SEEK_END);
       if (seekerr != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error positioning file pointer ") + fi);
       binsize = sizeof(T)*len*functions;
       off_t end;
#ifdef USELFS
       end = ftello64(dump);
#else
       end = ftello(dump);
#endif
       if (end == (off_t) -1) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error positioning file pointer ") + fi);
       size_ = (end - offs)/binsize;
       if(ferror(dump) != 0) throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error reading file: ") + fi);
       fseek(dump, pos, SEEK_SET);
    }
    else
    {
       throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] File not found: ") + fi);
    }
    if (dump == NULL)
    {
        throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] File not found: ") + fi);
    }
    return;
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::open( const std::string& filename, uint32_t l, uint32_t nrofFunctions)
{
    this->open(filename.c_str(), l, nrofFunctions);
}

template<class T>
void FileVector<std::valarray< std::valarray<T> > >::open( const char *const fi, uint32_t l, uint32_t nrofFunctions)
{
    filename = fi;
    len = l;
    functions = nrofFunctions;
    dump = fopen(fi, "rb");
    if (dump != NULL)
    {
        std::cerr<<"[FileVector] WARNING A FILE WITH THAT NAME: "<<filename<<" ALREADY EXISTS! OVERWRITING!"<<std::endl;
        fclose(dump);
    }
#ifdef USELFS
        dump = fopen64(fi, "wb+");
#else
        dump = fopen(fi, "wb+");//Create an empty file for both reading and writing. If a file with the same name already exists its content is erased and the file is treated as a new empty file. This is a binary file
#endif
    if (dump == NULL)
    {
        throw FileVector_Private::FileVectorOpenError(std::string("[FileVector] Error opening file: ") + fi);
    }
    else
    {
        if(fwrite(&len, sizeof(uint32_t), 1, dump) != 1) throw FileVector_Private::FileVectorWriteError("[FileVector] Error writing element to file!");
	if(fwrite(&functions, sizeof(uint32_t), 1, dump) != 1) throw FileVector_Private::FileVectorWriteError("[FileVector] Error writing element to file!");
    }
}

#include <unistd.h>

template<class T>
std::valarray<T> FileVector<std::valarray< std::valarray<T> > >::operator()(unsigned int i, unsigned int k)
{
  //we don't need to update curpos here, since the FILE* stream is not touched
  off64_t offset = static_cast<off64_t>(offs) + static_cast<off64_t>(i) * static_cast<off64_t>(binsize) + static_cast<off64_t>(k) * static_cast<off64_t>(len * sizeof(T));
  std::valarray<T> retval(len);
/*  fseeko(dump, offset, SEEK_SET);
  size_t readelems = fread(&(retval[0]), sizeof(T), len, dump);*/
  ssize_t readbytes = pread64(fileno(dump), &(retval[0]), len*sizeof(T), offset);//this led to a speed up of the error analysis of 30%
  if (readbytes != len*sizeof(T))
  {
    std::cout<<"[FileVector] pointer argument: "<<offset<<std::endl;
    std::cerr<<"[FileVector] Error reading file! "<<readbytes<< " "<<len<<std::endl;
  }
  return retval;
}

template<class T>
std::valarray<std::valarray<T> > FileVector<std::valarray<std::valarray<T> > >::operator[](const size_type n)
{
    value_type retval(functions);
    off_t offset = static_cast<off_t>(offs) + static_cast<off_t>(n) * static_cast<off_t>(binsize);
    if(curpos != offset)
    {
      fseeko(dump, offset, SEEK_SET);
      curpos = ftello(dump);
    }
    T *const tempmem(new T[len]);
    for (uint i = 0; i < functions; ++i)
    {
        retval[i].resize(len);
        size_t readelems = fread(tempmem, sizeof(T), len, dump);
        if (readelems != len)
        {
            std::stringstream s("[FileVector] Error reading file! ");
            s<< readelems << " " << len;
            throw FileVector_Private::FileVectorReadError(s.str());
        }
        for (unsigned int k = 0; k < len; ++k)
        {
            retval[i][k] = tempmem[k];
        }
    }
    delete [] tempmem;
    return retval;
}

template<class T>
typename FileVector<std::valarray<std::valarray<T> > >::size_type FileVector<std::valarray<std::valarray<T> > >::size() const throw()
{
    return size_;
}

template<class T>
void FileVector<std::valarray<std::valarray<T> > >::push_back(const value_type& arg)
{
  T *const tempmem(new T[arg[0].size()]);
  for(unsigned int i = 0; i < arg.size(); ++i)
  {
    for (unsigned int k = 0; k < arg[0].size(); ++k)
        tempmem[k] = arg[i][k];
    if(ftello(dump) != curpos)//avoids syncing as long as possible
      fseeko(dump, 0, SEEK_END);
    size_t writtenelems = fwrite(reinterpret_cast<char *const>(tempmem), sizeof(T), arg[i].size(), dump);
    if (writtenelems != arg[i].size())
    {
      delete [] tempmem;
      throw FileVector_Private::FileVectorWriteError("[FileVector] Error while writing data to file. File is perhaps read-only?");
    }
  }
  delete [] tempmem;
  ++size_;
  curpos = ftello(dump);
  return;
}

template<class T>
FileVector<std::valarray<std::valarray<T> > >::FileVector() : dump(NULL), size_(0), lock(NULL)
{
}

template<class T>
FileVector<std::valarray<std::valarray<T> > >::~FileVector()
{
    if (dump != NULL)
    {
        fclose(dump);
        delete lock;
    }
}

template<class T>
FileVector<std::valarray<std::valarray<T> > >::FileVector(const std::string& filename, uint32_t l) : size_(0), dump(NULL), len(l), lock(NULL)
{
    this->open(filename, l);
    return;
}
#endif
