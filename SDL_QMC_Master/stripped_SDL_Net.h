
/*
A version of the SDL_Net.h stripped of UDP support
*/

/*
    SDL_net:  An example cross-platform network library for use with SDL
    Copyright (C) 1997-2004 Sam Lantinga

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Sam Lantinga
    slouken@libsdl.org
*/

/* $Id: SDL_net.h 3281 2007-07-15 05:58:56Z slouken $ */

#ifndef _SDL_NET_H
#define _SDL_NET_H

//from SDL_stdinc.h
#include <SDL/SDL_config.h>


#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#if defined(STDC_HEADERS)
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
#else
# if defined(HAVE_STDLIB_H)
#  include <stdlib.h>
# elif defined(HAVE_MALLOC_H)
#  include <malloc.h>
# endif
# if defined(HAVE_STDDEF_H)
#  include <stddef.h>
# endif
# if defined(HAVE_STDARG_H)
#  include <stdarg.h>
# endif
#endif
#ifdef HAVE_STRING_H
# if !defined(STDC_HEADERS) && defined(HAVE_MEMORY_H)
#  include <memory.h>
# endif
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif
#if defined(HAVE_INTTYPES_H)
# include <inttypes.h>
#elif defined(HAVE_STDINT_H)
# include <stdint.h>
#endif
#ifdef HAVE_CTYPE_H
# include <ctype.h>
#endif
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#if defined(HAVE_ICONV) && defined(HAVE_ICONV_H)
# include <iconv.h>
#endif

/* The number of elements in an array */
#define SDL_arraysize(array)	(sizeof(array)/sizeof(array[0]))
#define SDL_TABLESIZE(table)	SDL_arraysize(table)

/* Use proper C++ casts when compiled as C++ to be compatible with the option
 -Wold-style-cast of GCC (and -Werror=old-style-cast in GCC 4.2 and above. */
#ifdef __cplusplus
#define SDL_reinterpret_cast(type, expression) reinterpret_cast<type>(expression)
#define SDL_static_cast(type, expression) static_cast<type>(expression)
#else
#define SDL_reinterpret_cast(type, expression) ((type)(expression))
#define SDL_static_cast(type, expression) ((type)(expression))
#endif

/* Basic data types */
typedef enum SDL_bool
{
    SDL_FALSE = 0,
    SDL_TRUE = 1
} SDL_bool;

/**
 * \typedef Sint8
 * \brief A signed 8-bit integer type.
 */
typedef int8_t Sint8;
/**
 * \typedef Uint8
 * \brief An unsigned 8-bit integer type.
 */
typedef uint8_t Uint8;
/**
 * \typedef Sint16
 * \brief A signed 16-bit integer type.
 */
typedef int16_t Sint16;
/**
 * \typedef Uint16
 * \brief An unsigned 16-bit integer type.
 */
typedef uint16_t Uint16;
/**
 * \typedef Sint32
 * \brief A signed 32-bit integer type.
 */
typedef int32_t Sint32;
/**
 * \typedef Uint32
 * \brief An unsigned 32-bit integer type.
 */
typedef uint32_t Uint32;

#ifdef SDL_HAS_64BIT_TYPE
/**
 * \typedef Sint64
 * \brief A signed 64-bit integer type.
 * \warning On platforms without any sort of 64-bit datatype, this is equivalent to Sint32!
 */
typedef int64_t Sint64;
/**
 * \typedef Uint64
 * \brief An unsigned 64-bit integer type.
 * \warning On platforms without any sort of 64-bit datatype, this is equivalent to Uint32!
 */
typedef uint64_t Uint64;
#else
/* This is really just a hack to prevent the compiler from complaining */
typedef Sint32 Sint64;
typedef Uint32 Uint64;
#endif

/* Make sure the types really have the right sizes */
#define SDL_COMPILE_TIME_ASSERT(name, x)               \
       typedef int SDL_dummy_ ## name[(x) * 2 - 1]
#ifndef DOXYGEN_SHOULD_IGNORE_THIS
SDL_COMPILE_TIME_ASSERT(uint8, sizeof(Uint8) == 1);
SDL_COMPILE_TIME_ASSERT(sint8, sizeof(Sint8) == 1);
SDL_COMPILE_TIME_ASSERT(uint16, sizeof(Uint16) == 2);
SDL_COMPILE_TIME_ASSERT(sint16, sizeof(Sint16) == 2);
SDL_COMPILE_TIME_ASSERT(uint32, sizeof(Uint32) == 4);
SDL_COMPILE_TIME_ASSERT(sint32, sizeof(Sint32) == 4);
#ifndef __NINTENDODS__          /* TODO: figure out why the following happens:
                                   include/SDL_stdinc.h:150: error: size of array 'SDL_dummy_uint64' is negative
                                   include/SDL_stdinc.h:151: error: size of array 'SDL_dummy_sint64' is negative */
SDL_COMPILE_TIME_ASSERT(uint64, sizeof(Uint64) == 8);
SDL_COMPILE_TIME_ASSERT(sint64, sizeof(Sint64) == 8);
#endif
#endif /* DOXYGEN_SHOULD_IGNORE_THIS */

/* Check to make sure enums are the size of ints, for structure packing.
   For both Watcom C/C++ and Borland C/C++ the compiler option that makes
   enums having the size of an int must be enabled.
   This is "-b" for Borland C/C++ and "-ei" for Watcom C/C++ (v11).
*/
/* Enable enums always int in CodeWarrior (for MPW use "-enum int") */
#ifdef __MWERKS__
#pragma enumsalwaysint on
#endif

#ifndef DOXYGEN_SHOULD_IGNORE_THIS
#ifndef __NINTENDODS__          /* TODO: include/SDL_stdinc.h:174: error: size of array 'SDL_dummy_enum' is negative */
typedef enum
{
    DUMMY_ENUM_VALUE
} SDL_DUMMY_ENUM;

SDL_COMPILE_TIME_ASSERT(enum, sizeof(SDL_DUMMY_ENUM) == sizeof(int));
#endif
#endif /* DOXYGEN_SHOULD_IGNORE_THIS */

#define SDL_LIL_ENDIAN	1234
#define SDL_BIG_ENDIAN	4321

#ifndef SDL_BYTEORDER	/* Not defined in SDL_config.h? */
#if defined(__hppa__) || \
    defined(__m68k__) || defined(mc68000) || defined(_M_M68K) || \
    (defined(__MIPS__) && defined(__MISPEB__)) || \
    defined(__ppc__) || defined(__POWERPC__) || defined(_M_PPC) || \
    defined(__sparc__)
#define SDL_BYTEORDER	SDL_BIG_ENDIAN
#else
#define SDL_BYTEORDER	SDL_LIL_ENDIAN
#endif
#endif /* !SDL_BYTEORDER */

//from begin_code.h

/* Some compilers use a special export keyword */
#ifndef DECLSPEC
# if defined(__BEOS__) || defined(__HAIKU__)
#  if defined(__GNUC__)
#   define DECLSPEC     __declspec(dllexport)
#  else
#   define DECLSPEC     __declspec(export)
#  endif
# elif defined(__WIN32__)
#  ifdef __BORLANDC__
#   ifdef BUILD_SDL
#    define DECLSPEC
#   else
#    define DECLSPEC    __declspec(dllimport)
#   endif
#  else
#   define DECLSPEC     __declspec(dllexport)
#  endif
# else
#  if defined(__GNUC__) && __GNUC__ >= 4
#   define DECLSPEC     __attribute__ ((visibility("default")))
#  else
#   define DECLSPEC
#  endif
# endif
#endif

/* By default SDL uses the C calling convention */
#ifndef SDLCALL
#if defined(__WIN32__) && !defined(__GNUC__)
#define SDLCALL __cdecl
#else
#define SDLCALL
#endif
#endif /* SDLCALL */

/* Removed DECLSPEC on Symbian OS because SDL cannot be a DLL in EPOC */
#ifdef __SYMBIAN32__
#undef DECLSPEC
#define DECLSPEC
#endif /* __SYMBIAN32__ */

/* Force structure packing at 4 byte alignment.
   This is necessary if the header is included in code which has structure
   packing set to an alternate value, say for loading structures from disk.
   The packing is reset to the previous value in close_code.h
 */
#if defined(_MSC_VER) || defined(__MWERKS__) || defined(__BORLANDC__)
#ifdef _MSC_VER
#pragma warning(disable: 4103)
#endif
#ifdef __BORLANDC__
#pragma nopackwarning
#endif
#pragma pack(push,4)
#endif /* Compiler needs structure packing set */

/* Set up compiler-specific options for inlining functions */
#ifndef SDL_INLINE_OKAY
#ifdef __GNUC__
#define SDL_INLINE_OKAY
#else
/* Add any special compiler-specific cases here */
#if defined(_MSC_VER) || defined(__BORLANDC__) || \
    defined(__DMC__) || defined(__SC__) || \
    defined(__WATCOMC__) || defined(__LCC__) || \
    defined(__DECC)
#ifndef __inline__
#define __inline__      __inline
#endif
#define SDL_INLINE_OKAY
#else
#if !defined(__MRC__) && !defined(_SGI_SOURCE)
#ifndef __inline__
#define __inline__ inline
#endif
#define SDL_INLINE_OKAY
#endif /* Not a funky compiler */
#endif /* Visual C++ */
#endif /* GNU C */
#endif /* SDL_INLINE_OKAY */

/* If inlining isn't supported, remove "__inline__", turning static
   inlined functions into static functions (resulting in code bloat
   in all files which include the offending header files)
*/
#ifndef SDL_INLINE_OKAY
#define __inline__
#endif

/* Apparently this is needed by several Windows compilers */
#if !defined(__MACH__)
#ifndef NULL
#ifdef __cplusplus
#define NULL 0
#else
#define NULL ((void *)0)
#endif
#endif /* NULL */
#endif /* ! Mac OS X - breaks precompiled headers */


/* Set up for C function definitions, even when using C++ */
#ifdef __cplusplus
extern "C" {
#endif

//from SDL.h

/* This function loads the SDL dynamically linked library and initializes 
 * the subsystems specified by 'flags' (and those satisfying dependencies)
 * Unless the SDL_INIT_NOPARACHUTE flag is set, it will install cleanup
 * signal handlers for some commonly ignored fatal signals (like SIGSEGV)
 */
extern DECLSPEC int SDLCALL SDL_Init(Uint32 flags);

/* This function cleans up all initialized subsystems and unloads the
 * dynamically linked library.  You should call it upon all exit conditions.
 */
extern DECLSPEC void SDLCALL SDL_Quit(void);

//from SDL_error.h
extern DECLSPEC void SDLCALL SDL_SetError(const char *fmt, ...);
extern DECLSPEC char * SDLCALL SDL_GetError(void);

//from SDL_timer.h

/* This is the OS scheduler timeslice, in milliseconds */
#define SDL_TIMESLICE		10

/* This is the maximum resolution of the SDL timer on all platforms */
#define TIMER_RESOLUTION	10	/* Experimentally determined */

/* Get the number of milliseconds since the SDL library initialization.
 * Note that this value wraps if the program runs for more than ~49 days.
 */ 
extern DECLSPEC Uint32 SDLCALL SDL_GetTicks(void);

/* Wait a specified number of milliseconds before returning */
extern DECLSPEC void SDLCALL SDL_Delay(Uint32 ms);

//from SDL_endian.h

/* Use inline functions for compilers that support them, and static
   functions for those that do not.  Because these functions become
   static for compilers that do not support inline functions, this
   header should only be included in files that actually use them.
*/
#if defined(__GNUC__) && defined(__i386__) && \
   !(__GNUC__ == 2 && __GNUC_MINOR__ <= 95 /* broken gcc version */)
static __inline__ Uint16 SDL_Swap16(Uint16 x)
{
	__asm__("xchgb %b0,%h0" : "=q" (x) :  "0" (x));
	return x;
}
#elif defined(__GNUC__) && defined(__x86_64__)
static __inline__ Uint16 SDL_Swap16(Uint16 x)
{
	__asm__("xchgb %b0,%h0" : "=Q" (x) :  "0" (x));
	return x;
}
#elif defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc__))
static __inline__ Uint16 SDL_Swap16(Uint16 x)
{
	Uint16 result;

	__asm__("rlwimi %0,%2,8,16,23" : "=&r" (result) : "0" (x >> 8), "r" (x));
	return result;
}
#elif defined(__GNUC__) && (defined(__M68000__) || defined(__M68020__))
static __inline__ Uint16 SDL_Swap16(Uint16 x)
{
	__asm__("rorw #8,%0" : "=d" (x) :  "0" (x) : "cc");
	return x;
}
#else
static __inline__ Uint16 SDL_Swap16(Uint16 x) {
	return((x<<8)|(x>>8));
}
#endif

#if defined(__GNUC__) && defined(__i386__) && \
   !(__GNUC__ == 2 && __GNUC_MINOR__ <= 95 /* broken gcc version */)
static __inline__ Uint32 SDL_Swap32(Uint32 x)
{
	__asm__("bswap %0" : "=r" (x) : "0" (x));
	return x;
}
#elif defined(__GNUC__) && defined(__x86_64__)
static __inline__ Uint32 SDL_Swap32(Uint32 x)
{
	__asm__("bswapl %0" : "=r" (x) : "0" (x));
	return x;
}
#elif defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc__))
static __inline__ Uint32 SDL_Swap32(Uint32 x)
{
	Uint32 result;

	__asm__("rlwimi %0,%2,24,16,23" : "=&r" (result) : "0" (x>>24), "r" (x));
	__asm__("rlwimi %0,%2,8,8,15"   : "=&r" (result) : "0" (result),    "r" (x));
	__asm__("rlwimi %0,%2,24,0,7"   : "=&r" (result) : "0" (result),    "r" (x));
	return result;
}
#elif defined(__GNUC__) && (defined(__M68000__) || defined(__M68020__))
static __inline__ Uint32 SDL_Swap32(Uint32 x)
{
	__asm__("rorw #8,%0\n\tswap %0\n\trorw #8,%0" : "=d" (x) :  "0" (x) : "cc");
	return x;
}
#else
static __inline__ Uint32 SDL_Swap32(Uint32 x) {
	return((x<<24)|((x<<8)&0x00FF0000)|((x>>8)&0x0000FF00)|(x>>24));
}
#endif

#ifdef SDL_HAS_64BIT_TYPE
#if defined(__GNUC__) && defined(__i386__) && \
   !(__GNUC__ == 2 && __GNUC_MINOR__ <= 95 /* broken gcc version */)
static __inline__ Uint64 SDL_Swap64(Uint64 x)
{
	union { 
		struct { Uint32 a,b; } s;
		Uint64 u;
	} v;
	v.u = x;
	__asm__("bswapl %0 ; bswapl %1 ; xchgl %0,%1" 
	        : "=r" (v.s.a), "=r" (v.s.b) 
	        : "0" (v.s.a), "1" (v.s.b)); 
	return v.u;
}
#elif defined(__GNUC__) && defined(__x86_64__)
static __inline__ Uint64 SDL_Swap64(Uint64 x)
{
	__asm__("bswapq %0" : "=r" (x) : "0" (x));
	return x;
}
#else
static __inline__ Uint64 SDL_Swap64(Uint64 x)
{
	Uint32 hi, lo;

	/* Separate into high and low 32-bit values and swap them */
	lo = (Uint32)(x&0xFFFFFFFF);
	x >>= 32;
	hi = (Uint32)(x&0xFFFFFFFF);
	x = SDL_Swap32(lo);
	x <<= 32;
	x |= SDL_Swap32(hi);
	return(x);
}
#endif
#else
/* This is mainly to keep compilers from complaining in SDL code.
   If there is no real 64-bit datatype, then compilers will complain about
   the fake 64-bit datatype that SDL provides when it compiles user code.
*/
#define SDL_Swap64(X)	(X)
#endif /* SDL_HAS_64BIT_TYPE */


/* Byteswap item from the specified endianness to the native endianness */
#if SDL_BYTEORDER == SDL_LIL_ENDIAN
#define SDL_SwapLE16(X)	(X)
#define SDL_SwapLE32(X)	(X)
#define SDL_SwapLE64(X)	(X)
#define SDL_SwapBE16(X)	SDL_Swap16(X)
#define SDL_SwapBE32(X)	SDL_Swap32(X)
#define SDL_SwapBE64(X)	SDL_Swap64(X)
#else
#define SDL_SwapLE16(X)	SDL_Swap16(X)
#define SDL_SwapLE32(X)	SDL_Swap32(X)
#define SDL_SwapLE64(X)	SDL_Swap64(X)
#define SDL_SwapBE16(X)	(X)
#define SDL_SwapBE32(X)	(X)
#define SDL_SwapBE64(X)	(X)
#endif

//from SDL_version.h

typedef struct SDL_version {
	Uint8 major;
	Uint8 minor;
	Uint8 patch;
} SDL_version;


/* Printable format: "%d.%d.%d", MAJOR, MINOR, PATCHLEVEL
*/
#define SDL_NET_MAJOR_VERSION	1
#define SDL_NET_MINOR_VERSION	2
#define SDL_NET_PATCHLEVEL	7

/* This macro can be used to fill a version structure with the compile-time
 * version of the SDL_net library.
 */
#define SDL_NET_VERSION(X)						\
{									\
	(X)->major = SDL_NET_MAJOR_VERSION;				\
	(X)->minor = SDL_NET_MINOR_VERSION;				\
	(X)->patch = SDL_NET_PATCHLEVEL;				\
}

#include <errno.h>

#ifdef macintosh
#ifndef USE_GUSI_SOCKETS
#define MACOS_OPENTRANSPORT
//#error Open Transport driver is broken
#endif
#endif /* macintosh */

/* Include system network headers */
#ifdef MACOS_OPENTRANSPORT
#include <OpenTransport.h>
#include <OpenTptInternet.h>
#else
#if defined(__WIN32__) || defined(WIN32)
#define __USE_W32_SOCKETS
#include <winsock.h>
#else /* UNIX */
#ifdef __OS2__
#include <types.h>
#include <sys/ioctl.h>
#endif
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <netinet/in.h>
#ifndef __BEOS__
#include <arpa/inet.h>
#endif
#ifdef linux /* FIXME: what other platforms have this? */
#include <netinet/tcp.h>
#endif
#include <netdb.h>
#include <sys/socket.h>
#endif /* WIN32 */
#endif /* Open Transport */

/* System-dependent definitions */
#ifdef MACOS_OPENTRANSPORT
//#define closesocket	OTCloseProvider
#define closesocket OTSndOrderlyDisconnect
#define SOCKET		EndpointRef
#define INVALID_SOCKET	kOTInvalidEndpointRef
#else
#ifndef __USE_W32_SOCKETS
#ifdef __OS2__
#define closesocket     soclose
#else  /* !__OS2__ */
#define closesocket	close
#endif /* __OS2__ */
#define SOCKET	int
#define INVALID_SOCKET	-1
#define SOCKET_ERROR	-1
#endif /* __USE_W32_SOCKETS */
#endif /* Open Transport */


/***********************************************************************/
/* Error reporting functions                                           */
/***********************************************************************/

/* We'll use SDL's functions for error reporting */
#define SDLNet_SetError	SDL_SetError
#define SDLNet_GetError	SDL_GetError



/* Initialize/Cleanup the network API
   SDL must be initialized before calls to functions in this library,
   because this library uses utility functions from the SDL library.
*/
extern int SDLNet_started;

#ifdef MACOS_OPENTRANSPORT

#include <Events.h>

extern DNSStatus dnsStatus;
extern Uint32 OTlocalhost;

/* We need a notifier for opening DNS.*/
/* ( 010311 masahiro minami<elsur@aaa.letter.co.jp>) */
static inline pascal void OpenDNSNotifier(
	void* context, OTEventCode code, OTResult result, void* cookie )
{
	switch( code )
	{
		case T_OPENCOMPLETE:
			// DNS is ready now.
			if( result == kOTNoError )
			{
				dnsStatus.dns = (InetSvcRef)cookie;
				dnsStatus.stat = dnsReady;
			}
			else
			{
				SDLNet_SetError("T_DNRSTRINGTOADDRCOMPLETE event returned an error");
				dnsStatus.dns = NULL;
				dnsStatus.stat = dnsError;
			}
			break;
		case T_DNRSTRINGTOADDRCOMPLETE:
			// DNR resolved the name to address
			// WORK IN PROGRESS (TODO )
			dnsStatus.stat = dnsResolved;
			break;
		default:
			if( result != kOTNoError )
				dnsStatus.stat = dnsError;
	}
	// Is there anything else to be done here ???
	// ( 010311 masahiro minami<elsur@aaa.letter.co.jp> )
	// (TODO)
}

/* Local functions for initializing and cleaning up the DNS resolver */
static inline int OpenDNS(void)
{
	int retval;
	OSStatus status;

	retval = 0;
	status = OTAsyncOpenInternetServices(
		kDefaultInternetServicesPath, 0, OpenDNSNotifier, NULL);
	if ( status == noErr ) {
		InetInterfaceInfo	info;
		
		dnsStatus.stat = dnsNotReady;
		
		while( dnsStatus.stat != dnsError && dnsStatus.dns == NULL)
		{
			// what's to be done ? Yield ? WaitNextEvent ? or what ?
			// ( 010311 masahiro minami<elsur@aaa.letter.co.jp> )
			//YieldToAnyThread();
		}
		/* Get the address of the local system -
		   What should it be if ethernet is off?
		 */
		OTInetGetInterfaceInfo(&info, kDefaultInetInterface);
		OTlocalhost = info.fAddress;
	} else {
		SDLNet_SetError("Unable to open DNS handle");
		retval = status;
	}
	
	return(retval);
}

static inline void CloseDNS(void)
{
	if ( dnsStatus.dns ) {
		OTCloseProvider(dnsStatus.dns);
		dnsStatus.dns = 0;
		dnsStatus.stat = dnsNotReady;
	}
	
	OTlocalhost = 0;
}


/* Initialize/Cleanup the network API */
static inline int  SDLNet_Init(void)
{
	OSStatus status;
	int retval;

	dnsStatus.stat = dnsNotReady;
	dnsStatus.dns = 0;


	retval = 0;
	if ( ! SDLNet_started ) {
		status = InitOpenTransport();
		if ( status == noErr ) {
			retval = OpenDNS();
			if ( retval < 0 ) {
				SDLNet_Quit();
			}
		} else {
			SDLNet_SetError("Unable to initialize Open Transport");
			retval = status;
		}
	}
	if ( retval == 0 ) {
		++SDLNet_started;
	}
	return(retval);
}

static inline void SDLNet_Quit(void)
{
	if ( SDLNet_started == 0 ) {
		return;
	}
	if ( --SDLNet_started == 0 ) {
		CloseDNS();
		CloseOpenTransport();
	}
}

#else /* !MACOS_OPENTRANSPORT */

#ifndef __USE_W32_SOCKETS
#include <signal.h>
#endif

/* Initialize/Cleanup the network API */
static inline int SDLNet_Init(void)
{
	if ( !SDLNet_started ) {
#ifdef __USE_W32_SOCKETS
		/* Start up the windows networking */
		WORD version_wanted = MAKEWORD(1,1);
		WSADATA wsaData;

		if ( WSAStartup(version_wanted, &wsaData) != 0 ) {
			SDLNet_SetError("Couldn't initialize Winsock 1.1\n");
			return(-1);
		}
#else
		/* SIGPIPE is generated when a remote socket is closed */
		void (*handler)(int);
		handler = signal(SIGPIPE, SIG_IGN);
		if ( handler != SIG_DFL ) {
			signal(SIGPIPE, handler);
		}
#endif
	}
	++SDLNet_started;
	return(0);
}

static inline void SDLNet_Quit(void)
{
	if ( SDLNet_started == 0 ) {
		return;
	}
	if ( --SDLNet_started == 0 ) {
#ifdef __USE_W32_SOCKETS
		/* Clean up windows networking */
		if ( WSACleanup() == SOCKET_ERROR ) {
			if ( WSAGetLastError() == WSAEINPROGRESS ) {
				WSACancelBlockingCall();
				WSACleanup();
			}
		}
#else
		/* Restore the SIGPIPE handler */
		void (*handler)(int);
		handler = signal(SIGPIPE, SIG_DFL);
		if ( handler != SIG_IGN ) {
			signal(SIGPIPE, handler);
		}
#endif
	}
}

#endif

/***********************************************************************/
/* IPv4 hostname resolution API                                        */
/***********************************************************************/

typedef struct {
	Uint32 host;			/* 32-bit IPv4 host address */
	Uint16 port;			/* 16-bit protocol port */
} IPaddress;

/* Resolve a host name and port to an IP address in network form.
   If the function succeeds, it will return 0.
   If the host couldn't be resolved, the host portion of the returned
   address will be INADDR_NONE, and the function will return -1.
   If 'host' is NULL, the resolved host will be set to INADDR_ANY.
 */
#ifndef INADDR_ANY
#define INADDR_ANY		0x00000000
#endif
#ifndef INADDR_NONE
#define INADDR_NONE		0xFFFFFFFF
#endif
#ifndef INADDR_BROADCAST
#define INADDR_BROADCAST	0xFFFFFFFF
#endif

#ifdef MACOS_OPENTRANSPORT

#include <Events.h>
/* Resolve a host name and port to an IP address in network form */
static inline int SDLNet_ResolveHost(IPaddress *address, const char *host, Uint16 port)
{
	int retval = 0;

	/* Perform the actual host resolution */
	if ( host == NULL ) {
		address->host = INADDR_ANY;
	} else {
/*		int a[4];

		address->host = INADDR_NONE;
		
		if ( sscanf(host, "%d.%d.%d.%d", a, a+1, a+2, a+3) == 4 ) {
			if ( !(a[0] & 0xFFFFFF00) && !(a[1] & 0xFFFFFF00) &&
			     !(a[2] & 0xFFFFFF00) && !(a[3] & 0xFFFFFF00) ) {
				address->host = ((a[0] << 24) |
				                 (a[1] << 16) |
				                 (a[2] <<  8) | a[3]);
				if ( address->host == 0x7F000001 ) {
					address->host = OTlocalhost;
				}
			}
		}
		
		if ( address->host == INADDR_NONE ) {*/
			InetHostInfo hinfo;
			
			/* Check for special case - localhost */
			if ( strcmp(host, "localhost") == 0 )
				return(SDLNet_ResolveHost(address, "127.0.0.1", port));

			/* Have OpenTransport resolve the hostname for us */
			retval = OTInetStringToAddress(dnsStatus.dns, (char *)host, &hinfo);
			if (retval == noErr) {
				while( dnsStatus.stat != dnsResolved )
					{WaitNextEvent(everyEvent, 0, 1, NULL );}
				address->host = hinfo.addrs[0];
			}
		//}
	}
	
	address->port = SDL_SwapBE16(port);

	/* Return the status */
	return(retval);
}

/* Resolve an ip address to a host name in canonical form.
   If the ip couldn't be resolved, this function returns NULL,
   otherwise a pointer to a static buffer containing the hostname
   is returned.  Note that this function is not thread-safe.
*/
/* MacOS implementation by Roy Wood
 */
static inline const char *SDLNet_ResolveIP(IPaddress *ip)
{
	if (ip != nil)
	{
	InetHost				theIP;
	static InetDomainName	theInetDomainName;
	OSStatus				theOSStatus;
	
		
		/*	Default result will be null string */
		
		theInetDomainName[0] = '\0';	
		
		
		/*	Do a reverse DNS lookup */
		
		theIP = ip->host;
		
		theOSStatus = OTInetAddressToName(dnsStatus.dns,theIP,theInetDomainName);
		
		/*	If successful, return the result */
			
		if (theOSStatus == kOTNoError)
		{
			while( dnsStatus.stat != dnsResolved )
				{ /*should we yield or what ? */ }
			return(theInetDomainName);
		}
	}
	
	SDLNet_SetError("Can't perform reverse DNS lookup");
	
	return(NULL);
}

#else /* !MACOS_OPENTRANSPORT */

#ifndef __USE_W32_SOCKETS
#include <signal.h>
#endif

/* Resolve a host name and port to an IP address in network form */
static inline int SDLNet_ResolveHost(IPaddress *address, const char *host, Uint16 port)
{
	int retval = 0;

	/* Perform the actual host resolution */
	if ( host == NULL ) {
		address->host = INADDR_ANY;
	} else {
		address->host = inet_addr(host);
		if ( address->host == INADDR_NONE ) {
			struct hostent *hp;

			hp = gethostbyname(host);
			if ( hp ) {
				memcpy(&address->host,hp->h_addr,hp->h_length);
			} else {
				retval = -1;
			}
		}
	}
	address->port = SDL_SwapBE16(port);

	/* Return the status */
	return(retval);
}

/* Resolve an ip address to a host name in canonical form.
   If the ip couldn't be resolved, this function returns NULL,
   otherwise a pointer to a static buffer containing the hostname
   is returned.  Note that this function is not thread-safe.
*/
/* Written by Miguel Angel Blanch.
 * Main Programmer of Arianne RPG.
 * http://come.to/arianne_rpg
 */
static inline const char *SDLNet_ResolveIP(IPaddress *ip)
{
	struct hostent *hp;

	hp = gethostbyaddr((char *)&ip->host, 4, AF_INET);
	if ( hp != NULL ) {
		return hp->h_name;
	}
  	return NULL;
}
#endif

/***********************************************************************/
/* TCP network API                                                     */
/***********************************************************************/

typedef struct _TCPsocket *TCPsocket;
#ifdef MACOS_OPENTRANSPORT

#include <Events.h>
#include <Threads.h>
#include <OpenTransport.h>
#include <OpenTptInternet.h>
#include <OTDebug.h>

struct _TCPsocket {
	int ready;
	SOCKET channel;
	
	// These are taken from GUSI interface.
	// I'm not sure if it's really necessary here yet
	// ( masahiro minami<elsur@aaa.letter.co.jp> )
	// ( 01/02/19 )
	OTEventCode curEvent;
	OTEventCode newEvent;
	OTEventCode event;
	OTEventCode curCompletion;
	OTEventCode newCompletion;
	OTEventCode completion;
	OSStatus	error;
	TEndpointInfo info;
	Boolean		readShutdown;
	Boolean		writeShutdown;
	Boolean		connected;
	OTConfigurationRef	config;		// Master configuration. you can clone this.
	TCPsocket	nextListener;
	// ( end of new members --- masahiro minami<elsur@aaa.letter.co.jp>
	
	IPaddress remoteAddress;
	IPaddress localAddress;
	int sflag;
	
	// Maybe we don't need this---- it's from original SDL_net
	// (masahiro minami<elsur@aaa.letter.co.jp>)
	// ( 01/02/20 )
	int rcvdPassConn;
};

// To be used in WaitNextEvent() here and there....
// (010311 masahiro minami<elsur@aaa.letter.co.jp>)
extern EventRecord macEvent;

#if TARGET_API_MAC_CARBON
/* for Carbon */
extern OTNotifyUPP notifier;
#endif

/* A helper function for Mac OpenTransport support*/
// This function is a complete copy from GUSI
// ( masahiro minami<elsur@aaa.letter.co.jp> )
// ( 01/02/19 )
static __inline__ Uint32 CompleteMask(OTEventCode code)	
{ 	
	return 1 << (code & 0x1F);
}

/* Create a new TCPsocket */
// Because TCPsocket structure gets bigger and bigger,
// I think we'd better have a constructor function and delete function.
// ( 01/02/25 masahiro minami<elsur@aaa.letter.co.jp> )
static inline TCPsocket AsyncTCPNewSocket()
{
	TCPsocket sock;
	
	sock = (TCPsocket)malloc(sizeof(*sock));
	if ( sock == NULL ) {
		SDLNet_SetError("Out of memory");
		return NULL;
	}
	
	sock->newEvent = 0;
	sock->event = 0;
	sock->curEvent = 0;
	sock->newCompletion = 0;
	sock->completion = 0;
	sock->curCompletion = 0;
	//sock->info = NULL;
	sock->readShutdown = sock->writeShutdown = sock->connected = false;
	sock->error = 0;
	sock->config = NULL;
	sock->nextListener = NULL;
	sock->sflag = 0;
	return sock;	
}

/* Retrieve OT event */
// This function is taken from GUSI interface.
// ( 01/02/19 masahiro minami<elsur@aaa.letter.co.jp> )
static inline void AsyncTCPPopEvent( TCPsocket sock )
{
	// Make sure OT calls are not interrupted
	// Not sure if we really need this.
	OTEnterNotifier( sock->channel );
	
	sock->event |= (sock->curEvent = sock->newEvent );
	sock->completion |= ( sock->curCompletion = sock->newCompletion );
	sock->newEvent = sock->newCompletion = 0;
	
	OTLeaveNotifier( sock->channel );
	
	if( sock->curEvent & T_UDERR)
	{
		// We just clear the error.
		// Should we feed this back to users ?
		// (TODO )
		OTRcvUDErr( sock->channel, NULL );
		
#ifdef DEBUG_NET
		printf("AsyncTCPPopEvent  T_UDERR recognized");
#endif
	}
	
	// Remote is disconnecting...
	if( sock->curEvent & ( T_DISCONNECT | T_ORDREL ))
	{
		sock->readShutdown = true;
	}
	
	if( sock->curEvent &T_CONNECT )
	{
		// Ignore the info of remote (second parameter).
		// Shoule we care ?
		// (TODO)
		OTRcvConnect( sock->channel, NULL );
		sock->connected = 1;
	}
	
	if( sock->curEvent & T_ORDREL )
	{
		OTRcvOrderlyDisconnect( sock->channel );
	}
	
	if( sock->curEvent & T_DISCONNECT )
	{
		OTRcvDisconnect( sock->channel, NULL );
	}
	
	// Do we need to ?
	// (masahiro minami<elsur@aaa.letter.co.jp>)
	//YieldToAnyThread();
}

/* Close a TCP network socket */
static inline void SDLNet_TCP_Close(TCPsocket sock)
{
	if ( sock != NULL ) {
		if ( sock->channel != INVALID_SOCKET ) {
			//closesocket(sock->channel);
			OTSndOrderlyDisconnect( sock->channel );
		}
		free(sock);
	}
}

/* Accept an incoming connection on the given server socket.
   The newly created socket is returned, or NULL if there was an error.
*/
static inline TCPsocket SDLNet_TCP_Accept(TCPsocket server)
{
	
	/* Only server sockets can accept */
	if ( ! server->sflag ) {
		SDLNet_SetError("Only server sockets can accept()");
		return(NULL);
	}
	server->ready = 0;

	/* Accept a new TCP connection on a server socket */
	{
		InetAddress peer;
		TCall peerinfo;
		TCPsocket sock = NULL;
		Boolean mustListen = false;
		OTResult err;
		
		memset(&peerinfo, 0, (sizeof peerinfo ));
		peerinfo.addr.buf = (Uint8 *) &peer;
		peerinfo.addr.maxlen = sizeof(peer);
		
		while( mustListen || !sock )
		{
			// OTListen
			// We do NOT block ---- right thing ? (TODO)
			err = OTListen( server->channel, &peerinfo );

			if( err )
				goto error_return;
			else
			{
				mustListen = false;
				sock = AsyncTCPNewSocket();
				if( ! sock )
					goto error_return;
			}
		}
		if( sock )
		{
			// OTAsyncOpenEndpoint
			server->error = OTAsyncOpenEndpoint( OTCloneConfiguration( server->config ),
				NULL, &(sock->info), (OTNotifyProcPtr)AsyncTCPNotifier, sock );
			AsyncTCPPopEvent( sock );
			while( !sock->error && !( sock->completion & CompleteMask( T_OPENCOMPLETE)))
			{
				AsyncTCPPopEvent( sock );
			}
			if( ! sock->channel )
			{
				mustListen = false;
				goto error_return;
			}
			
			// OTAccept
			server->completion &= ~(CompleteMask(T_ACCEPTCOMPLETE));
			server->error = OTAccept( server->channel, sock->channel, &peerinfo );
			AsyncTCPPopEvent( server );
			while( !(server->completion & CompleteMask(T_ACCEPTCOMPLETE)))
			{
				AsyncTCPPopEvent( server );
			}
			
			switch( server->error )
			{
				case kOTLookErr:
					switch( OTLook(server->channel ))
					{
						case T_LISTEN:
							mustListen = true;
							break;
						case T_DISCONNECT:
							goto error_return;
					}
					break;
				case 0:
					sock->nextListener = server->nextListener;
					server->nextListener = sock;
					sock->remoteAddress.host = peer.fHost;
					sock->remoteAddress.port = peer.fPort;
					return sock;
					// accept successful
					break;
				default:
					free( sock );
			}
		}
		sock->remoteAddress.host = peer.fHost;
		sock->remoteAddress.port = peer.fPort;
		sock->sflag = 0;
		sock->ready = 0;

		/* The socket is ready */
		return(sock);
	
	// Error; close the socket and return	
	error_return:
		SDLNet_TCP_Close(sock);
		return(NULL);
	}
}

/* Send 'len' bytes of 'data' over the non-server socket 'sock'
   This function returns the actual amount of data sent.  If the return value
   is less than the amount of data sent, then either the remote connection was
   closed, or an unknown socket error occurred.
*/
static inline int SDLNet_TCP_Send(TCPsocket sock, const void *datap, int len)
{
	const Uint8 *data = (const Uint8 *)datap;	/* For pointer arithmetic */
	int sent, left;

	/* Server sockets are for accepting connections only */
	if ( sock->sflag ) {
		SDLNet_SetError("Server sockets cannot send");
		return(-1);
	}

	/* Keep sending data until it's sent or an error occurs */
	left = len;
	sent = 0;
	errno = 0;
	do {
		len = OTSnd(sock->channel, (void *)data, left, 0);
		if (len == kOTFlowErr)
			len = 0;
		if ( len > 0 ) {
			sent += len;
			left -= len;
			data += len;
		}
		// Do we need to ?
		// ( masahiro minami<elsur@aaa.letter.co.jp> )
		// (TODO)
		//WaitNextEvent(everyEvent, &macEvent, 1, NULL);
		//AsyncTCPPopEvent(sock);
	} while ( (left > 0) && (len > 0) );

	return(sent);
}

/* Open a TCP network socket
   If 'remote' is NULL, this creates a local server socket on the given port,
   otherwise a TCP connection to the remote host and port is attempted.
   The newly created socket is returned, or NULL if there was an error.

   ( re-written by masahiro minami<elsur@aaa.letter.co.jp>
     Now endpoint is created in Async mode.
     01/02/20 )
*/
TCPsocket SDLNet_TCP_Open(IPaddress *ip)
{
	EndpointRef dummy = NULL;

	TCPsocket sock = AsyncTCPNewSocket();
	if( ! sock)
		return NULL;

	// Determin whether bind locally, or connect to remote
	if ( (ip->host != INADDR_NONE) && (ip->host != INADDR_ANY) )
	{
		// ######## Connect to remote
		OTResult stat;
		InetAddress inAddr;
		TBind bindReq;

		// Open endpoint
		sock->error = OTAsyncOpenEndpoint(
			OTCreateConfiguration(kTCPName), NULL, &(sock->info),
			(OTNotifyProcPtr)(AsyncTCPNotifier),
			sock );

		AsyncTCPPopEvent( sock );
		while( !sock->error && !( sock->completion & CompleteMask(T_OPENCOMPLETE)))
		{
			//SetThreadState( kCurrentThreadID, kReadyThreadState, kNoThreadID );
			//YieldToAnyThread();
			//WaitNextEvent(everyEvent, &macEvent, 1, NULL);
			AsyncTCPPopEvent( sock );
		}

		if( !sock->channel )
		{
			SDLNet_SetError("OTAsyncOpenEndpoint failed --- client socket could not be opened");
			goto error_return;
		}

		// Set blocking mode
		// I'm not sure if this is a good solution....
		// Check out Apple's sample code, OT Virtual Server
		// ( 010314 masahiro minami<elsur@aaa.letter.co.jp>)

		sock->error = OTSetBlocking( sock->channel );
		if( sock->error != kOTNoError )
		{
			SDLNet_SetError("OTSetBlocking() returned an error");
			goto error_return;
		}

		// Bind the socket
		OTInitInetAddress(&inAddr, 0, 0 );
		bindReq.addr.len = sizeof( InetAddress );
		bindReq.addr.buf = (unsigned char*)&inAddr;
		bindReq.qlen = 0;

		sock->error = OTBind( sock->channel, &bindReq, NULL );
		AsyncTCPPopEvent(sock);
		while( !sock->error && !( sock->completion & CompleteMask(T_BINDCOMPLETE)))
		{
			//YieldToAnyThread();
			//WaitNextEvent(everyEvent, &macEvent, 1, NULL);
			AsyncTCPPopEvent(sock);
		}


		switch( stat = OTGetEndpointState( sock->channel ))
		{
			InetAddress inAddr;
			TCall sndCall;
			OTResult res;

			case T_OUTCON:
				SDLNet_SetError("SDLNet_Open() failed -- T_OUTCON");
				goto error_return;
				break;
			case T_IDLE:
				sock->readShutdown = false;
				sock->writeShutdown = false;
				sock->event &=~T_CONNECT;

				OTMemzero(&sndCall, sizeof(TCall));
				OTInitInetAddress(&inAddr, ip->port, ip->host );
				sndCall.addr.len = sizeof(InetAddress);
				sndCall.addr.buf = (unsigned char*)&inAddr;
				sock->connected = 0;
				res = OTConnect( sock->channel, &sndCall, NULL );
				AsyncTCPPopEvent(sock);
				while( sock->error == kOTNoDataErr || !sock->connected )
					AsyncTCPPopEvent(sock);
				break;
			default:
				// What's to be done ? (TODO)
				SDLNet_SetError("SDLNet_TCP_Open() failed -- EndpointState not good");
				goto error_return;

		}
		if( !(sock->event & (T_CONNECT|T_DISCONNECT)))
			goto error_return;

		AsyncTCPPopEvent( sock );
		while( !(sock->event & (T_CONNECT|T_DISCONNECT)))
		{
			AsyncTCPPopEvent( sock );
		}
		// OTConnect successfull
		if( sock->event & T_CONNECT)
		{
			sock->remoteAddress.host = inAddr.fHost;
			sock->remoteAddress.port = inAddr.fPort;
			sock->sflag = false;
		}
		else
		{
			// OTConnect failed
			sock->event &= ~T_DISCONNECT;
			goto error_return;
		}
	}
	else
	{
		// ######## Bind locally
		TBind bindReq;
		InetAddress	inAddr;

	// First, get InetInterfaceInfo.
	// I don't search for all of them.
	// Does that matter ?

		sock->error = OTAsyncOpenEndpoint(
			OTCreateConfiguration("tilisten, tcp"), NULL, &(sock->info),
			(OTNotifyProcPtr)(AsyncTCPNotifier),
			sock);
		AsyncTCPPopEvent( sock );
		while( !sock->error && !( sock->completion & CompleteMask( T_OPENCOMPLETE)))
		{
			AsyncTCPPopEvent( sock );
		}

		if( ! sock->channel )
		{
			SDLNet_SetError("OTAsyncOpenEndpoint failed --- server socket could not be opened");
			goto error_return;
		}

		// Create a master OTConfiguration
		sock->config = OTCreateConfiguration(kTCPName);
		if( ! sock->config )
		{
			SDLNet_SetError("Could not create master OTConfiguration");
			goto error_return;
		}

		// Bind the socket
		OTInitInetAddress(&inAddr, ip->port, 0 );
		inAddr.fAddressType = AF_INET;
		bindReq.addr.len = sizeof( InetAddress );
		bindReq.addr.buf = (unsigned char*)&inAddr;
		bindReq.qlen = 35;	// This number is NOT well considered. (TODO)
		sock->localAddress.host = inAddr.fHost;
		sock->localAddress.port = inAddr.fPort;
		sock->sflag = true;
		
		sock->error = OTBind( sock->channel, &bindReq, NULL );
		AsyncTCPPopEvent(sock);
		while( !sock->error && !( sock->completion & CompleteMask(T_BINDCOMPLETE)))
		{
			AsyncTCPPopEvent(sock);
		}
		if( sock->error != kOTNoError )
		{
			SDLNet_SetError("Could not bind server socket");
			goto error_return;
		}
		
		if( dummy )
			OTCloseProvider( dummy );

	}
	
	sock->ready = 0;
	return sock;
	
	error_return:
	if( dummy )
		OTCloseProvider( dummy );
	SDLNet_TCP_Close( sock );
	return NULL;	
}

/* Receive up to 'maxlen' bytes of data over the non-server socket 'sock',
   and store them in the buffer pointed to by 'data'.
   This function returns the actual amount of data received.  If the return
   value is less than or equal to zero, then either the remote connection was
   closed, or an unknown socket error occurred.
*/
static inline int SDLNet_TCP_Recv(TCPsocket sock, void *data, int maxlen)
{
	int len = 0;
	OSStatus res;
	/* Server sockets are for accepting connections only */
	if ( sock->sflag ) {
		SDLNet_SetError("Server sockets cannot receive");
		return(-1);
	}

	do
	{
		res = OTRcv(sock->channel, data, maxlen-len, 0);
		if (res > 0) {
			len = res;
		}

#ifdef DEBUG_NET
		if ( res != kOTNoDataErr )
			printf("SDLNet_TCP_Recv received ; %d\n", res );
#endif
		
		AsyncTCPPopEvent(sock);
		if( res == kOTLookErr )
		{
			res = OTLook(sock->channel );
			continue;
		}
	} while ( (len == 0) && (res == kOTNoDataErr) );

	sock->ready = 0;
	if ( len == 0 ) { /* Open Transport error */
#ifdef DEBUG_NET
		printf("Open Transport error: %d\n", res);
#endif
		return(-1);
	}
	return(len);
}

/* Get the IP address of the remote system associated with the socket.
   If the socket is a server socket, this function returns NULL.
*/
IPaddress *SDLNet_TCP_GetPeerAddress(TCPsocket sock)
{
	if ( sock->sflag ) {
		return(NULL);
	}
	return(&sock->remoteAddress);
}
#else /* !MACOS_OPENTRANSPORT */

struct _TCPsocket {
	int ready;
	SOCKET channel;
	IPaddress remoteAddress;
	IPaddress localAddress;
	int sflag;
};

/* Close a TCP network socket  */
static inline void SDLNet_TCP_Close(TCPsocket sock)
{
	if ( sock != NULL ) {
		if ( sock->channel != INVALID_SOCKET ) {
			closesocket(sock->channel);
		}
		free(sock);
	}
}

/* Accept an incoming connection on the given server socket.
   The newly created socket is returned, or NULL if there was an error.
*/
static inline TCPsocket SDLNet_TCP_Accept(TCPsocket server)
{
	TCPsocket sock;
	struct sockaddr_in sock_addr;
#ifdef USE_GUSI_SOCKETS
        unsigned int sock_alen;
#elif defined(WIN32)
        int sock_alen;
#else
        socklen_t sock_alen; /*default to some POSIX behaviour*/
#endif

	/* Only server sockets can accept */
	if ( ! server->sflag ) {
		SDLNet_SetError("Only server sockets can accept()");
		return(NULL);
	}
	server->ready = 0;

	/* Allocate a TCP socket structure */
	sock = (TCPsocket)malloc(sizeof(*sock));
	if ( sock == NULL ) {
		SDLNet_SetError("Out of memory");
		goto error_return;
	}

	/* Accept a new TCP connection on a server socket */
	sock_alen = sizeof(sock_addr);
	sock->channel = accept(server->channel, (struct sockaddr *)&sock_addr, &sock_alen);
	if ( sock->channel == INVALID_SOCKET ) {
		SDLNet_SetError("accept() failed");
		goto error_return;
	}
#ifdef WIN32
	{
		/* passing a zero value, socket mode set to block on */
		unsigned long mode = 0;
		ioctlsocket (sock->channel, FIONBIO, &mode);
	}
#elif defined(O_NONBLOCK)
	{
		int flags = fcntl(sock->channel, F_GETFL, 0);
		fcntl(sock->channel, F_SETFL, flags & ~O_NONBLOCK);
	}
#endif /* WIN32 */
	sock->remoteAddress.host = sock_addr.sin_addr.s_addr;
	sock->remoteAddress.port = sock_addr.sin_port;

	sock->sflag = 0;
	sock->ready = 0;

	/* The socket is ready */
	return(sock);

error_return:
	SDLNet_TCP_Close(sock);
	return(NULL);
}

/* Open a TCP network socket
   If 'remote' is NULL, this creates a local server socket on the given port,
   otherwise a TCP connection to the remote host and port is attempted.
   The newly created socket is returned, or NULL if there was an error.
*/
static inline TCPsocket SDLNet_TCP_Open(IPaddress *ip)
{
	TCPsocket sock;
	struct sockaddr_in sock_addr;

	/* Allocate a TCP socket structure */
	sock = (TCPsocket)malloc(sizeof(*sock));
	if ( sock == NULL ) {
		SDLNet_SetError("Out of memory");
		goto error_return;
	}

	/* Open the socket */
	sock->channel = socket(AF_INET, SOCK_STREAM, 0);
	if ( sock->channel == INVALID_SOCKET ) {
		SDLNet_SetError("Couldn't create socket");
		goto error_return;
	}

	/* Connect to remote, or bind locally, as appropriate */
	if ( (ip->host != INADDR_NONE) && (ip->host != INADDR_ANY) ) {

	// #########  Connecting to remote
	
		memset(&sock_addr, 0, sizeof(sock_addr));
		sock_addr.sin_family = AF_INET;
		sock_addr.sin_addr.s_addr = ip->host;
		sock_addr.sin_port = ip->port;

		/* Connect to the remote host */
		if ( connect(sock->channel, (struct sockaddr *)&sock_addr,
				sizeof(sock_addr)) == SOCKET_ERROR ) {
			SDLNet_SetError("Couldn't connect to remote host");
			goto error_return;
		}
		sock->sflag = 0;
	} else {

	// ##########  Binding locally

		memset(&sock_addr, 0, sizeof(sock_addr));
		sock_addr.sin_family = AF_INET;
		sock_addr.sin_addr.s_addr = INADDR_ANY;
		sock_addr.sin_port = ip->port;

/*
 * Windows gets bad mojo with SO_REUSEADDR:
 * http://www.devolution.com/pipermail/sdl/2005-September/070491.html
 *   --ryan.
 */
#ifndef WIN32
		/* allow local address reuse */
		{ int yes = 1;
			setsockopt(sock->channel, SOL_SOCKET, SO_REUSEADDR, (char*)&yes, sizeof(yes));
		}
#endif

		/* Bind the socket for listening */
		if ( bind(sock->channel, (struct sockaddr *)&sock_addr,
				sizeof(sock_addr)) == SOCKET_ERROR ) {
			SDLNet_SetError("Couldn't bind to local port");
			goto error_return;
		}
		if ( listen(sock->channel, 5) == SOCKET_ERROR ) {
			SDLNet_SetError("Couldn't listen to local port");
			goto error_return;
		}

		/* Set the socket to non-blocking mode for accept() */
#if defined(__BEOS__) && defined(SO_NONBLOCK)
		/* On BeOS r5 there is O_NONBLOCK but it's for files only */
		{
			long b = 1;
			setsockopt(sock->channel, SOL_SOCKET, SO_NONBLOCK, &b, sizeof(b));
		}
#elif defined(O_NONBLOCK)
		{
			fcntl(sock->channel, F_SETFL, O_NONBLOCK);
		}
#elif defined(WIN32)
		{
			unsigned long mode = 1;
			ioctlsocket (sock->channel, FIONBIO, &mode);
		}
#elif defined(__OS2__)
		{
			int dontblock = 1;
			ioctl(sock->channel, FIONBIO, &dontblock);
		}
#else
#warning How do we set non-blocking mode on other operating systems?
#endif
		sock->sflag = 1;
	}
	sock->ready = 0;

#ifdef TCP_NODELAY
	/* Set the nodelay TCP option for real-time games */
	{ int yes = 1;
	setsockopt(sock->channel, IPPROTO_TCP, TCP_NODELAY, (char*)&yes, sizeof(yes));
	}
#endif /* TCP_NODELAY */

	/* Fill in the channel host address */
	sock->remoteAddress.host = sock_addr.sin_addr.s_addr;
	sock->remoteAddress.port = sock_addr.sin_port;

	/* The socket is ready */
	return(sock);

error_return:
	SDLNet_TCP_Close(sock);
	return(NULL);
}

/* Send 'len' bytes of 'data' over the non-server socket 'sock'
   This function returns the actual amount of data sent.  If the return value
   is less than the amount of data sent, then either the remote connection was
   closed, or an unknown socket error occurred.
*/
static inline int SDLNet_TCP_Send(TCPsocket sock, const void *datap, int len)
{
	const Uint8 *data = (const Uint8 *)datap;	/* For pointer arithmetic */
	int sent, left;

	/* Server sockets are for accepting connections only */
	if ( sock->sflag ) {
		SDLNet_SetError("Server sockets cannot send");
		return(-1);
	}

	/* Keep sending data until it's sent or an error occurs */
	left = len;
	sent = 0;
	errno = 0;
	do {
		len = send(sock->channel, (const char *) data, left, 0);//Flo: we use send instead of write because of portability
		if ( len > 0 ) {
			sent += len;
			left -= len;
			data += len;
		}
	} while ( (left > 0) && ((len > 0) || (errno == EINTR)) );
	return(sent);
}

/* Receive up to 'maxlen' bytes of data over the non-server socket 'sock',
   and store them in the buffer pointed to by 'data'.
   This function returns the actual amount of data received.  If the return
   value is less than or equal to zero, then either the remote connection was
   closed, or an unknown socket error occurred.
*/
static inline int SDLNet_TCP_Recv(TCPsocket sock, void *data, int maxlen)
{
	int len;

	/* Server sockets are for accepting connections only */
	if ( sock->sflag ) {
		SDLNet_SetError("Server sockets cannot receive");
		return(-1);
	}

	errno = 0;
	do {
		len = recv(sock->channel, (char *) data, maxlen, 0);//Flo: we use recv instead of read because of portability
	} while ( errno == EINTR );
	sock->ready = 0;
	return(len);
}

/* Get the IP address of the remote system associated with the socket.
   If the socket is a server socket, this function returns NULL.
*/
static inline IPaddress *SDLNet_TCP_GetPeerAddress(TCPsocket sock)
{
	if ( sock->sflag ) {
		return(NULL);
	}
	return(&sock->remoteAddress);
}
#endif /* MACOS_OPENTRANSPORT */

/***********************************************************************/
/* Hooks for checking sockets for available data                       */
/***********************************************************************/
struct SDLNet_Socket {
	int ready;
	SOCKET channel;
#ifdef MACOS_OPENTRANSPORT
	OTEventCode curEvent;
#endif
};

struct _SDLNet_SocketSet {
	int numsockets;
	int maxsockets;
	struct SDLNet_Socket **sockets;
};


typedef _SDLNet_SocketSet *SDLNet_SocketSet;

/* Any network socket can be safely cast to this socket type */
typedef struct {
	int ready;
} *SDLNet_GenericSocket;

/* Allocate a socket set for use with SDLNet_CheckSockets() 
   This returns a socket set for up to 'maxsockets' sockets, or NULL if 
   the function ran out of memory. 
 */
static inline SDLNet_SocketSet SDLNet_AllocSocketSet(int maxsockets)
{
	struct _SDLNet_SocketSet *set;
	int i;

	set = (struct _SDLNet_SocketSet *)malloc(sizeof(*set));
	if ( set != NULL ) {
		set->numsockets = 0;
		set->maxsockets = maxsockets;
		set->sockets = (struct SDLNet_Socket **)malloc
					(maxsockets*sizeof(*set->sockets));
		if ( set->sockets != NULL ) {
			for ( i=0; i<maxsockets; ++i ) {
				set->sockets[i] = NULL;
			}
		} else {
			free(set);
			set = NULL;
		}
	}
	return(set);
}

/* Add a socket to a set of sockets to be checked for available data */
#define SDLNet_TCP_AddSocket(set, sock) \
			SDLNet_AddSocket(set, (SDLNet_GenericSocket)sock)
#define SDLNet_UDP_AddSocket(set, sock) \
			SDLNet_AddSocket(set, (SDLNet_GenericSocket)sock)
/* Add a socket to a set of sockets to be checked for available data */
static inline int SDLNet_AddSocket(SDLNet_SocketSet set, SDLNet_GenericSocket sock)
{
	if ( sock != NULL ) {
		if ( set->numsockets == set->maxsockets ) {
			SDLNet_SetError("socketset is full");
			return(-1);
		}
		set->sockets[set->numsockets++] = (struct SDLNet_Socket *)sock;
	}
	return(set->numsockets);
}


/* Remove a socket from a set of sockets to be checked for available data */
#define SDLNet_TCP_DelSocket(set, sock) \
			SDLNet_DelSocket(set, (SDLNet_GenericSocket)sock)
#define SDLNet_UDP_DelSocket(set, sock) \
			SDLNet_DelSocket(set, (SDLNet_GenericSocket)sock)
static inline int SDLNet_DelSocket(SDLNet_SocketSet set, SDLNet_GenericSocket sock)
{
	int i;

	if ( sock != NULL ) {
		for ( i=0; i<set->numsockets; ++i ) {
			if ( set->sockets[i] == (struct SDLNet_Socket *)sock ) {
				break;
			}
		}
		if ( i == set->numsockets ) {
			SDLNet_SetError("socket not found in socketset");
			return(-1);
		}
		--set->numsockets;
		for ( ; i<set->numsockets; ++i ) {
			set->sockets[i] = set->sockets[i+1];
		}
	}
	return(set->numsockets);
}

/* This function checks to see if data is available for reading on the
   given set of sockets.  If 'timeout' is 0, it performs a quick poll,
   otherwise the function returns when either data is available for
   reading, or the timeout in milliseconds has elapsed, which ever occurs
   first.  This function returns the number of sockets ready for reading,
   or -1 if there was an error with the select() system call.
*/
#ifdef MACOS_OPENTRANSPORT
static inline int SDLNet_CheckSockets(SDLNet_SocketSet set, Uint32 timeout)
{
Uint32	stop;
int 	numReady;

	/* Loop, polling the network devices */
	
	stop = SDL_GetTicks() + timeout;
	
	do 
	{
	OTResult status;
	size_t	numBytes;
	int 	i;
		
		numReady = 0;
	
		for (i = set->numsockets-1;i >= 0;--i) 
		{
			status = OTLook( set->sockets[i]->channel );
			if( status > 0 )
			{
				switch( status )
				{
					case T_UDERR:
						OTRcvUDErr( set->sockets[i]->channel , nil);
						break;
					case T_DISCONNECT:
						OTRcvDisconnect( set->sockets[i]->channel, nil );
						break;
					case T_ORDREL:
						OTRcvOrderlyDisconnect(set->sockets[i]->channel );
						break;
					case T_CONNECT:
						OTRcvConnect( set->sockets[i]->channel, nil );
						break;
					
				
					default:
						set->sockets[i]->ready = 1;
						++numReady;
				}
			}
			else if( OTCountDataBytes(set->sockets[i]->channel, &numBytes ) != kOTNoDataErr )
			{
				set->sockets[i]->ready = 1;
				++numReady;
			}
			else
				set->sockets[i]->ready = 0;
		}
		
	} while (!numReady && (SDL_GetTicks() < stop));

	return(numReady);
}
#else
static inline int SDLNet_CheckSockets(SDLNet_SocketSet set, Uint32 timeout)
{
	int i;
	SOCKET maxfd;
	int retval;
	struct timeval tv;
	fd_set mask;

	/* Find the largest file descriptor */
	maxfd = 0;
	for ( i=set->numsockets-1; i>=0; --i ) {
		if ( set->sockets[i]->channel > maxfd ) {
			maxfd = set->sockets[i]->channel;
		}
	}
	/* Check the file descriptors for available data */
	do {
		errno = 0;

		/* Set up the mask of file descriptors */
		FD_ZERO(&mask);
		for ( i=set->numsockets-1; i>=0; --i ) {
			FD_SET(set->sockets[i]->channel, &mask);
		}

		/* Set up the timeout */
		tv.tv_sec = timeout/1000;
		tv.tv_usec = (timeout%1000)*1000;

		/* Look! */
		retval = select(maxfd+1, &mask, NULL, NULL, &tv);
	} while ( errno == EINTR );

	/* Mark all file descriptors ready that have data available */
	if ( retval > 0 ) {
		for ( i=set->numsockets-1; i>=0; --i ) {
			if ( FD_ISSET(set->sockets[i]->channel, &mask) ) {
				set->sockets[i]->ready = 1;
			}
		}
	}
	return(retval);
}
#endif /* MACOS_OPENTRANSPORT */

/* After calling SDLNet_CheckSockets(), you can use this function on a
   socket that was in the socket set, to find out if data is available
   for reading.
*/
#define SDLNet_SocketReady(sock) \
		((sock != NULL) && ((SDLNet_GenericSocket)sock)->ready)

/* Free a set of sockets allocated by SDL_NetAllocSocketSet() */
/* Free a set of sockets allocated by SDL_NetAllocSocketSet() */
static inline void SDLNet_FreeSocketSet(SDLNet_SocketSet set)
{
	if ( set ) {
		free(set->sockets);
		free(set);
	}
}

/***********************************************************************/
/* Platform-independent data conversion functions                      */
/***********************************************************************/

#if !SDL_DATA_ALIGNED /* function versions for binary compatibility */

/* Write a 16 bit value to network packet buffer */
#undef SDLNet_Write16
inline void   SDLNet_Write16(Uint16 value, void *areap)
{
	(*(Uint16 *)(areap) = SDL_SwapBE16(value));
}

/* Write a 32 bit value to network packet buffer */
#undef SDLNet_Write32
inline void   SDLNet_Write32(Uint32 value, void *areap)
{
	*(Uint32 *)(areap) = SDL_SwapBE32(value);
}

/* Read a 16 bit value from network packet buffer */
#undef SDLNet_Read16
inline Uint16 SDLNet_Read16(void *areap)
{
	return (SDL_SwapBE16(*(Uint16 *)(areap)));
}

/* Read a 32 bit value from network packet buffer */
#undef SDLNet_Read32
inline Uint32 SDLNet_Read32(void *areap)
{
	return (SDL_SwapBE32(*(Uint32 *)(areap)));
}

#endif /* !SDL_DATA_ALIGNED */

/* Inline macro functions to read/write network data */

/* Warning, some systems have data access alignment restrictions */
#if defined(sparc) || defined(mips)
#define SDL_DATA_ALIGNED	1
#endif
#ifndef SDL_DATA_ALIGNED
#define SDL_DATA_ALIGNED	0
#endif

/* Write a 16 bit value to network packet buffer */
#if !SDL_DATA_ALIGNED
#define SDLNet_Write16(value, areap)	\
	(*(Uint16 *)(areap) = SDL_SwapBE16(value))
#else
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
#define SDLNet_Write16(value, areap)	\
do 					\
{					\
	Uint8 *area = (Uint8 *)(areap);	\
	area[0] = (value >>  8) & 0xFF;	\
	area[1] =  value        & 0xFF;	\
} while ( 0 )
#else
#define SDLNet_Write16(value, areap)	\
do 					\
{					\
	Uint8 *area = (Uint8 *)(areap);	\
	area[1] = (value >>  8) & 0xFF;	\
	area[0] =  value        & 0xFF;	\
} while ( 0 )
#endif
#endif /* !SDL_DATA_ALIGNED */

/* Write a 32 bit value to network packet buffer */
#if !SDL_DATA_ALIGNED
#define SDLNet_Write32(value, areap) 	\
	*(Uint32 *)(areap) = SDL_SwapBE32(value);
#else
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
#define SDLNet_Write32(value, areap) 	\
do					\
{					\
	Uint8 *area = (Uint8 *)(areap);	\
	area[0] = (value >> 24) & 0xFF;	\
	area[1] = (value >> 16) & 0xFF;	\
	area[2] = (value >>  8) & 0xFF;	\
	area[3] =  value       & 0xFF;	\
} while ( 0 )
#else
#define SDLNet_Write32(value, areap) 	\
do					\
{					\
	Uint8 *area = (Uint8 *)(areap);	\
	area[3] = (value >> 24) & 0xFF;	\
	area[2] = (value >> 16) & 0xFF;	\
	area[1] = (value >>  8) & 0xFF;	\
	area[0] =  value       & 0xFF;	\
} while ( 0 )
#endif
#endif /* !SDL_DATA_ALIGNED */

/* Read a 16 bit value from network packet buffer */
#if !SDL_DATA_ALIGNED
#define SDLNet_Read16(areap) 		\
	(SDL_SwapBE16(*(Uint16 *)(areap)))
#else
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
#define SDLNet_Read16(areap) 		\
	((((Uint8 *)areap)[0] <<  8) | ((Uint8 *)areap)[1] <<  0)
#else
#define SDLNet_Read16(areap) 		\
	((((Uint8 *)areap)[1] <<  8) | ((Uint8 *)areap)[0] <<  0)
#endif
#endif /* !SDL_DATA_ALIGNED */

/* Read a 32 bit value from network packet buffer */
#if !SDL_DATA_ALIGNED
#define SDLNet_Read32(areap) 		\
	(SDL_SwapBE32(*(Uint32 *)(areap)))
#else
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
#define SDLNet_Read32(areap) 		\
	((((Uint8 *)areap)[0] << 24) | (((Uint8 *)areap)[1] << 16) | \
	 (((Uint8 *)areap)[2] <<  8) |  ((Uint8 *)areap)[3] <<  0)
#else
#define SDLNet_Read32(areap) 		\
	((((Uint8 *)areap)[3] << 24) | (((Uint8 *)areap)[2] << 16) | \
	 (((Uint8 *)areap)[1] <<  8) |  ((Uint8 *)areap)[0] <<  0)
#endif
#endif /* !SDL_DATA_ALIGNED */

#ifdef MACOS_OPENTRANSPORT
#endif
/* Ends C function definitions when using C++ */
#ifdef __cplusplus
}
#endif

//from close_code.h
#if defined(_MSC_VER) || defined(__MWERKS__) || defined(__WATCOMC__)  || defined(__BORLANDC__)
#ifdef __BORLANDC__
#pragma nopackwarning
#endif
#pragma pack(pop)
#endif /* Compiler needs structure packing set */

#endif /* _SDL_NET_H */

