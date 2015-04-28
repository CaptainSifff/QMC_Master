#ifndef MODELS_H
#define MODELS_H
#include <inttypes.h>
//typed enum. That's an C++0x feature, which is present since gcc-4.4. there's no need to give a special flag to gcc. it just emits a warning
enum Models :
uint32_t
{
    HUBBARD_CHAIN = 0,
    COLD_ATOMS = 1,
    SIAM = 2,
    IMAG_MODEL_FROM_FILE = 3,
    KONDO_IMP_TI = 4,
    RASHBA_CHAIN = 5,
    RASHBA_CHAIN_EXPONENTIAL = 6,
    RASHBA_CHAIN_POWER_LAW = 7
};
#endif
