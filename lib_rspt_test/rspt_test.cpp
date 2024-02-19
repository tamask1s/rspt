#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <numeric>
#include <vector>
#include <functional>
#include <iostream>
#include <math.h>

/** lib_rspt interfaces */
#include "../lib_rspt/signal_packer.h"
#include "../lib_rspt/iir_filter.h"

/** Bellow includes are not necessary to use the rspt lib, but the functionalities can be used to easy data manipulation. */
#include "../lib_rspt/lib_signalpacker/utils.h"
#include "../lib_rspt/lib_zaxtensor/ZaxJsonParser.h"
#include "../lib_rspt/lib_zaxtensor/ZaxTensor.h"

using namespace std;

bool write_buffer_(const char* filename, const unsigned char* buffer, int len)
{
    std::ofstream file(filename, std::ios::binary | std::ios::trunc);
    if (!file.is_open())
        return false;

    if (!file.write((char*)buffer, len))
        return false;

    file.close();
    return true;
}

bool read_buffer_(const char* filename, std::vector<char>& buffer)
{
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file.is_open())
        return false;

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    buffer.resize(size);
    if (!file.read(buffer.data(), size))
        return false;

    file.close();
    return true;
}

void test_packer_(i_signal_packer* cmpr, int nr_samples_to_encode, int Channels, uint8_t* data_stream, int BYTESPERSAMPLE)
{
    write_buffer_("_original.bin", data_stream, nr_samples_to_encode * Channels * BYTESPERSAMPLE);
    /** Create a [Channels x nr_samples_to_encode] matrix to store the original data in an int32_t matrix */
    tensor_i32 orig;
    orig.resize(Channels, nr_samples_to_encode);
    convert_native_to_i32(orig.d2d, data_stream, nr_samples_to_encode, Channels, BYTESPERSAMPLE, false);

    /** allocate a space to compress the data in, then compress the input data */
    size_t dst_max_len = nr_samples_to_encode * Channels * BYTESPERSAMPLE * 2;
    size_t dst_len_;
    unsigned char dst[dst_max_len];
    cmpr->compress(data_stream, dst, dst_max_len, dst_len_);
    cout << "dst_len_: " << dst_len_ << endl;

    /** Allocat space for decompression, then decompress the compressed data. */
    size_t compressed_len;
    unsigned char decdst[dst_max_len];
    cmpr->decompress(dst, compressed_len, decdst);
    cout << "compressed len: " << compressed_len << "   COMPRESSION CR = " << (double)(Channels * BYTESPERSAMPLE * nr_samples_to_encode) / compressed_len << std::endl;

    /** Write the decoded data to the disc. In case of lossless compression this must match the input data */
    write_buffer_("_decoded.bin", decdst, nr_samples_to_encode * Channels * BYTESPERSAMPLE);

    /** Create a [Channels x nr_samples_to_encode] matrix to store the decoded data in an int32_t matrix */
    tensor_i32 enc;
    enc.resize(Channels, nr_samples_to_encode);
    convert_native_to_i32(enc.d2d, decdst, nr_samples_to_encode, Channels, BYTESPERSAMPLE, false);

    /** Calculate and print PRDN */
    double MSE = 0;
    double origg = 0;
    for (int j = 0; j < Channels; ++j)
    {
        int32_t mean = average_32(orig.d2d[j], nr_samples_to_encode);
        for (int i = 0; i < nr_samples_to_encode; ++i)
        {
            double tmp = orig.d2d[j][i] - enc.d2d[j][i];
            MSE += tmp * tmp;
            origg += (orig.d2d[j][i] - mean) * (orig.d2d[j][i] - mean);
        }
    }
    cout << "PRDN[%] " << sqrt(MSE / origg) * 100.0 << endl;
}

void test_data(uint8_t* data_stream, int BYTESPERSAMPLE, int nr_channels, int nr_samples, bool filter_data)
{
    if (filter_data)
    {
        /** We convert the native data to a [nr_channels x nr_samples] sized matrix. */
        tensor_i32 enc;
        enc.resize(nr_channels, nr_samples);
        convert_native_to_i32(enc.d2d, data_stream, nr_samples, nr_channels, BYTESPERSAMPLE, false);
        /** Bellow filter coefficients represent a bandpass butterworth filter of 0.4-200Hz @ 2000Sps */
        double n[] = {1.00000000000, -3.14332095199, 3.70064088865, -1.97083923944, 0.41351972908};
        double d[] = {0.06722876941, 0.00000000000, -0.13445753881, 0.00000000000, 0.06722876941};
        /** Create the filter object and use it on the data */
        i_filter* filter = i_filter::new_iir(n, d, 5);
        for (int j = 0; j < nr_channels; ++j)
        {
            filter->init_history_values(enc.d2d[j][0]);
            for (int i = 0; i < nr_samples; ++i)
                enc.d2d[j][i] = filter->filter_opt(enc.d2d[j][i]);
        }
        /** We convert the [nr_channels x nr_samples] sized matrix back to native data. */
        convert_i32_to_native(data_stream, enc.d2d, nr_samples, nr_channels, BYTESPERSAMPLE, false);
    }

    /** Initialize packer. For the sake of efficiency, packers needs to know about the structure of native data. */
    /** This is why not a simple [size] is given as an argument, but [BYTESPERSAMPLE, nr_channels, nr_samples] */
    cout << "\n******************************* Testing packers *******************************\nxdelta_hzr:" << endl;
    i_signal_packer* cmpr = i_signal_packer::new_xdelta_hzr(BYTESPERSAMPLE, nr_channels, nr_samples);
    test_packer_(cmpr, nr_samples, nr_channels, data_stream, BYTESPERSAMPLE);

    /** For transformation based compression methods we need to provide a number of samples of 2^N, therefore we truncate the data @ 16384 samples */
    cout << "hadamard:" << endl;
    cmpr = i_signal_packer::new_hadamard(BYTESPERSAMPLE, nr_channels, 16384);
    test_packer_(cmpr, 16384, nr_channels, data_stream, BYTESPERSAMPLE);

    /** For dct we only provide fewer samples to encode as it is slow, therefore we truncate the data @ 4096 samples */
    cout << "dct:" << endl;
    cmpr = i_signal_packer::new_dct(BYTESPERSAMPLE, nr_channels, 4096);
    test_packer_(cmpr, 4096, nr_channels, data_stream, BYTESPERSAMPLE);
}

void test_1()
{
    int BYTESPERSAMPLE = 3;
    std::vector<char> data_stream;
    int nr_channels = 3;
    /** "data_stream.bin": ECG data sampled @ 2000Sps, 3 Channels and 24 bit resolution, each sample stored in 3 bytes.
    * Total number of samples: 60.000, 20.000 Samples per channel.
    * Structure: B0: most significant byte, B2: less significant byte
    *    - [CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2][CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2] ...
    *    - see the function convert_native_to_i32() */
    if (read_buffer_("data_stream.bin", data_stream))
    {
        int nr_samples = data_stream.size() / (nr_channels * BYTESPERSAMPLE);
        test_data((uint8_t*)data_stream.data(), BYTESPERSAMPLE, nr_channels, nr_samples, true);
    }
}

void test_2()
{
    const int BYTESPERSAMPLE = 4;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 32 bit resolution, stored in an array of int32_t. */
    /** Total number of samples: 16384, total data size: 16384 * 4 Bytes. */
    int32_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;
    test_data((uint8_t*)data_stream, BYTESPERSAMPLE, nr_channels, nr_samples, true);
}

void test_3()
{
    const int BYTESPERSAMPLE = 2;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 16 bit resolution, stored in an array of int16_t. */
    /** Total number of samples: 16384, total data size: 32768 Bytes. */
    int16_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;
    test_data((uint8_t*)data_stream, BYTESPERSAMPLE, nr_channels, nr_samples, true);
}

void test_4()
{
    const int BYTESPERSAMPLE = 1;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 8 bit resolution, stored in an array of int8_t. */
    /** Total number of samples: 16384, total data size: 16384 Bytes. */
    int8_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 100.0;
    test_data((uint8_t*)data_stream, BYTESPERSAMPLE, nr_channels, nr_samples, true);
}

void test_5()
{
    /** Simulate simple sinusoidal data. 1 Channels and 32 bit resolution, stored in an array of int32_t. */
    /** Total number of samples: 8192, total data size: 32768 Bytes. */
    const int BYTESPERSAMPLE = 4;
    const int nr_samples = 8192;
    const int nr_channels = 1;
    int32_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;

    /** Initialize packer. For the sake of efficiency, packers needs to know about the structure of native data. */
    /** This is why not a simple [size] is given as an argument, but [BYTESPERSAMPLE, nr_channels, nr_samples] */
    i_signal_packer* cmpr = i_signal_packer::new_xdelta_hzr(BYTESPERSAMPLE, nr_channels, nr_samples);

    /** allocate sufficient room for compressed data, then compress the data */
    size_t dst_max_len = nr_samples * nr_channels * BYTESPERSAMPLE * 2;
    unsigned char dst[dst_max_len];
    size_t compressed_size;
    cmpr->compress((uint8_t*)data_stream, dst, dst_max_len, compressed_size);
    cout << "compressed_size: " << compressed_size;

    /** Allocate space for decompression, then decompress the compressed data. */
    size_t decompressed_size;
    unsigned char decdst[dst_max_len];
    cmpr->decompress(dst, decompressed_size, decdst);
    cout << "  compressed len: " << decompressed_size << " compression CR = " << (double)(nr_channels * BYTESPERSAMPLE * nr_samples) / decompressed_size << std::endl;
}

int main()
{
    test_1();
    test_2();
    test_3();
    test_4();
    test_5();
    return 0;
}
