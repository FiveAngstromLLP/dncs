#include "libdncs.h"
#include <iostream>
#include <cstdlib> // For std::free

int main() {
    const char* seq = "YGGFM";
    const char* folder = "Result";
    const char* ff = "amber-fb15";
    size_t sample_count = 100;
    const char* method = "fold";

    int result = sample(seq, folder, ff, sample_count, method);
    if (result != 0) {
        std::cerr << "Sample function returned error code: " << result << std::endl;
        return 1; // Indicate error in main
    }

    return 0;
}
