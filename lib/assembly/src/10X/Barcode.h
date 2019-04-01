
#ifndef TENX_BARCODE_H
#define TENX_BARCODE_H

#include "CoreTools.h"

// actually, uint32_t should enough alerady
typedef int64_t barcode_t;

barcode_t createBarcodeIndex(String &barcode_str);

#endif
