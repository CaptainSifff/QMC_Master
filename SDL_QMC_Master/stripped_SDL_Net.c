#ifdef macintosh
#ifndef USE_GUSI_SOCKETS
#define MACOS_OPENTRANSPORT
//#error Open Transport driver is broken
#endif
#endif /* macintosh */

int SDLNet_started = 0;

#ifdef MACOS_OPENTRANSPORT

#include <Events.h>
#include <Threads.h>
#include <OpenTransport.h>
#include <OpenTptInternet.h>
//hopefully this header exists on MoacOS and provides the same functionality as the SDL...
#include <inttypes.h>

DNSStatus dnsStatus;
uint32_t OTlocalhost = 0;

// To be used in WaitNextEvent() here and there....
// (010311 masahiro minami<elsur@aaa.letter.co.jp>)
EventRecord macEvent;

#if TARGET_API_MAC_CARBON
/* for Carbon */
OTNotifyUPP notifier;
#endif

#endif
