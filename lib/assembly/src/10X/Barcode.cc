#include "10X/Barcode.h"

barcode_t createBarcodeIndex(String &barcode_str) {

    std::size_t pos;
    barcode_t barcode_int = stol(barcode_str, &pos) * 1537 * 1537;
    barcode_str = barcode_str.substr(pos + 1);
    barcode_int += stol(barcode_str, &pos) * 1537;
    barcode_str = barcode_str.substr(pos + 1);
    barcode_int += stol(barcode_str);
    
    return barcode_int;
}
