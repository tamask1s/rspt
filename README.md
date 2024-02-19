# Rspt library

The Rspt is a C++ library designed to facilitate the compression and filtering of digital time-domain signal sequences. The signal processing algorithms used target real-time processing during sampling and are optimized accordingly.

## Compression

The library's compression algorithms operate on fixed-size input data, typically those encountered during real-time signal sampling. Therefore, data sizes must be specified when creating the compression object, while subsequent calls to the compression function do not allow for data size specification.
During compression initialization, size specifications are not provided as a single number, as the efficiency of the compressors requires knowledge of the internal structures of the data to be compressed. Thus, under data sizes, we understand 3 pieces of information:
- BYTESPERSAMPLE: Indicates the resolution of the sampled data, in bytes. A typical ADC often provides data on 24 bits, in which case the value of BYTESPERSAMPLE should be 3.
- nr_channels: Number of data channels. This is the number of individual signals present in the input data, like different leads in an ECG signal.
- nr_samples: Number of samples to be compressed in each channel.
The total size in bytes is the product of these: BYTESPERSAMPLE x nr_channels x nr_samples. With knowledge of the internal structure of the data to be compressed, the compressor can optimize compression.

Currently, four types of compressors are implemented, which can be created with the appropriate factory functions (lib_rspt/signal_packer.h):
- hzr: Lossless. Simple RLE + Huffman coding.
- xdelta hzr: Lossless. Combination of delta encoding, offseting, xor encoding, and hzr compression.
- dct: Compression based on DCT transformation, with uniform quantization, combined with hzr compression.
- hadamard: Compression based on Hadamard-Welsh transformation, with uniform quantization, combined with hzr compression.

## Filtering

Two types of digital filters are implemented: IIR and FIR.
Designing the coefficients of the appropriate filters must be done with another tool, as it is not included in the Rspt library. Therefore, filter initialization is not based on filtering frequencies and sampling frequency, but directly on the provision of the filter's coefficients.

When using the filter, the code using the library does not need to preserve the previous input and output values, as this is done by the Rspt library. The filter() or filter_opt() function will always return a filtered output, with the preservation of the history.
Interface and factory functions of the filetrs are provided in lib_rspt/iir_filter.h file.

### IIR:
- During initialization, n and d denote the numerator and denominator of the digital IIR filter.

### FIR:
- During initialization, the content and size of the filter kernel must be provided.

## License

This library is licensed under the Apache 2 license, but it contains other libraries with different licenses.

### Directory: lib_fwht
 FWHT - Fast Walsh-Hadamard Transform in C
 BSD license.
 by Sven Nilsen, 2012
 http://www.cutoutpro.com
 Please check it in the /lib_fwht directory

### Directory: libhzr 
 hzr - A Huffman + RLE compression library.
 Copyright (C) 2016 Marcus Geelnard
 Permission is granted to anyone to use this software for any purpose, including commercial applications.
 Please check it in the /libhzr directory
