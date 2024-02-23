# Rspt Library

Rspt is a C++ library designed to facilitate the compression and filtering of digital time-domain signal sequences. The signal processing algorithms used are targeting real-time processing during sampling and are optimized accordingly.

## Compression

The library's compression algorithms operate on fixed-size input data, typically those encountered during real-time signal sampling. Therefore, data sizes must be specified when creating the compression object, while subsequent calls to the compression function do not allow for data size specification.
During compression initialization, size specifications are not provided as a single number, as the efficiency of the compressors requires knowledge of the internal structures of the data to be compressed. Thus, under data sizes, we understand 3 pieces of information:
- *bytes_per_sample*: Indicates the resolution of the sampled data, in bytes. A typical ADC often provides data on 24 bits, in which case the value of *bytes_per_sample* should be 3.
- *nr_channels*: Number of data channels. This is the number of individual signals present in the input data, like different leads in an ECG signal.
- *nr_samples*: Number of samples to be compressed in each channel.

The total size in bytes is the product of these: *bytes_per_sample* x *nr_channels* x *nr_samples*. With knowledge of the internal structure of the data to be compressed, the compressor can optimize compression.

Currently, four types of compressors are implemented, which can be created with the appropriate factory functions (lib_rspt/signal_packer.h):
- [*hzr*](#lib_hzr): Lossless. [RLE](#abbreviations) + [Huffman](#abbreviations) coding.
- *xdelta hzr*: Lossless. Combination of delta encoding, offseting, xor encoding, and hzr compression.
- *dct*: Compression based on DCT transformation, with uniform quantization, combined with hzr compression.
- *hadamard*: Compression based on  [Hadamard-Walsh Transform (HWT)](#abbreviations), with uniform quantization, combined with hzr compression. See the [fwht library](#lib_fwht)

Interface and factory functions of the filters are provided in [lib_rspt/signal_packer.h](https://github.com/tamask1s/rspt/blob/main/lib_rspt/signal_packer.h) file.

## Filtering

Two types of digital filters are implemented: [IIR](#abbreviations) and [FIR](#abbreviations).
Designing the coefficients of the appropriate filters must be done with another tool, as it is not included in the Rspt library. Therefore, filter initialization is not based on filtering frequencies and sampling frequency, but directly on the provision of the filter's coefficients.

When using the filter, the code using the library does not need to preserve the previous input and output values, as this is done by the Rspt library. The filter() or filter_opt() functions will always return a filtered output, with the preservation of the history.
Interface and factory functions of the filters are provided in [lib_rspt/filter.h](https://github.com/tamask1s/rspt/blob/main/lib_rspt/filter.h) file.

### [IIR](#abbreviations):
- During initialization, two pointers to arrays of doubles - *n* and *d* must be provided. These are the numerator and the denominator of the digital [IIR](#abbreviations) filter coefficients.

### [FIR](#abbreviations):
- During initialization, a pointer to the array of doubles containing the kernel coefficients must be provided.

## Examples

### Compression

Simulating, compressing, and decompressing sinusoidal data with a lossless compression algorithm:

```cpp
#include <iostream>
#include <math.h>

#include "signal_packer.h"

int main()
{
    /** Simulate simple sinusoidal data. 1 Channels and 32 bit resolution, stored in an */
    /** array of int32_t. Total number of samples: 8192, total data size: 32768 Bytes. */
    const int bytes_per_sample = 4;
    const int nr_samples = 8192;
    const int nr_channels = 1;
    int32_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;

    /** Initialize packer. For the sake of efficiency, packers needs to know about the */
    /** internal structure of the native data. This is why not a simple [size] is */
    /**  given as an argument, but [bytes_per_sample, nr_channels, nr_samples] */
    i_signal_packer* c = i_signal_packer::new_xdelta_hzr(bytes_per_sample, nr_channels, nr_samples);

    /** allocate sufficient room for compressed data, then compress the data */
    size_t dst_max_len = nr_samples * nr_channels * bytes_per_sample * 2;
    unsigned char dst[dst_max_len];
    size_t compressed_size;
    c->compress((uint8_t*)data_stream, dst, dst_max_len, compressed_size);
    std::cout << "compressed_size: " << compressed_size;

    /** Allocate space for decompression, then decompress the compressed data. */
    size_t cmpr_size;
    unsigned char decdst[bytes_per_sample * nr_samples * nr_channels];
    c->decompress(dst, cmpr_size, decdst);
    std::cout << "  compressed len: " << cmpr_size << " compression CR = ";
    std::cout << (double)(nr_channels * bytes_per_sample * nr_samples) / cmpr_size << std::endl;
    return 0;
}
```

#### Result:

```cpp
compressed_size: 2022  compressed len: 2022 compression CR = 16.2057
```

#### Comparison of different packers:

xdelta_hzr: Compression ratio CR = 15.9 PRDN[%] = 0. PRDN is 0 as it is a lossless algorythm.

![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/compression_xdelta_hzr.png)

hadamard: Compression ratio CR = 52.7 PRDN[%] = 2.2. Better CR, but we have a non-0 PRDN

![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/compression_hadamard.png)

dct: Compression ratio CR = 142.4 PRDN[%] = 1.5. Even better CR, but note the artifacts on the beginning of the signal.

![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/compression_dct.png)

Legend:

![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/legend_.png)

For PRDN formula and different quality metrics, please check:

https://www.researchgate.net/figure/List-of-reconstructed-ECG-quality-assessment-tool_tbl2_269935665

DCT compression on real ECG data. CR ~= 24 PRDN ~= 3.5:

![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/compression_dct_ecg.png)

### Filtering

#### [IIR](#abbreviations)

Simulate a combination of 4Hz + 70Hz data

```cpp
/** Simulate sinusoidal data. 1 Channels and 32 bit resolution, stored in an */
/** array of int32_t. Signal with 2 sinusoids: 4Hz and 70Hz combined. */
#include <inttypes.h>
#include <math.h>
#include "iir_filter.h"

const int nr_samples = 8192;
const double sample_rate = 2000;
const double sr_pi = 3.14159265358979323846 * 2.0 / sample_rate;
int32_t data_stream[nr_samples];
for (int i = 0; i < nr_samples; ++i)
    data_stream[i] = sin(i * sr_pi* 4.0) * 1000.0 + sin(i * sr_pi * 70.0) * 1000.0;
```
Result:
![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/filtering_orig.png)

Filter it with a low-pass filter of 5Hz

```cpp
double n[] = {1.00000000000, -1.97778648378, 0.97803050849}; /// LP 5Hz @ 2kSps
double d[] = {0.00006100618, 0.00012201236, 0.00006100618};
i_filter* lp_filter = i_filter::new_iir(n, d, 3);
lp_filter->init_history_values(data_stream[0], sample_rate);
for (int i = 0; i < nr_samples; ++i)
    data_stream[i] = lp_filter->filter_opt(data_stream[i]);
```

Result: only the 4Hz component remains in the signal.
![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/filtering_lp.png)

Reset the data, then filter it with a high-pass filter of 50Hz.
```cpp
for (int i = 0; i < nr_samples; ++i)
    data_stream[i] = sin(i * sr_pi * 4.0) * 1000.0 + sin(i * sr_pi * 70.0) * 1000.0;

double n2[] = {1.00000000000, -1.77863177782, 0.80080264667}; /// HP 50Hz @ 2kSps
double d2[] = {0.89485860612, -1.78971721225, 0.89485860612};
i_filter* hp_filter = i_filter::new_iir(n2, d2, 3);
hp_filter->init_history_values(data_stream[0], sample_rate);
for (int i = 0; i < nr_samples; ++i)
    data_stream[i] = hp_filter->filter_opt(data_stream[i]);
```
Result: only the 70Hz component remains in the signal.
![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/filtering_hp.png)

Original + filtered signals:
![alt text](https://github.com/tamask1s/rspt/blob/main/lib_rspt_doc/filtering_hp_all.png)

## Abbreviations

- **Hadamard-Walsh Transform (HWT)**: [Hadamard-Walsh Transform](https://en.wikipedia.org/wiki/Hadamard_transform)
- **FIR**: [Finite Impulse Response (FIR) Filter](https://en.wikipedia.org/wiki/Finite_impulse_response)
- **IIR**: [Infinite Impulse Response (IIR) Filter](https://en.wikipedia.org/wiki/Infinite_impulse_response)
- **DCT**: [Discrete Cosine Transform (DCT)](https://en.wikipedia.org/wiki/Discrete_cosine_transform)
- **Delta**: [Delta encoding](https://en.wikipedia.org/wiki/Delta_encoding)
- **RLE**: [Run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding)
- **Huffmann**: [Huffmann coding](https://en.wikipedia.org/wiki/Huffman_coding)

## License

This library is licensed under the Apache 2 license, but it contains other libraries with different licenses.

### lib_fwht
 https://github.com/bvssvni/fwht
 
 FWHT - Fast Walsh-Hadamard Transform in C
 
 BSD license.
 
 by Sven Nilsen, 2012
 
 Please check it in the /lib_fwht directory

### lib_hzr 
 https://github.com/mbitsnbites/hzr
 
 https://gitlab.com/mbitsnbites/hzr
 
 hzr - A Huffman + RLE compression library.
 
 Copyright (C) 2016 Marcus Geelnard
 
 Permission is granted to anyone to use this software for any purpose, including commercial applications.
 
 Please check it in the /libhzr directory
