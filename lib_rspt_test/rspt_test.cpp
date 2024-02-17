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

void test_packer_(i_signal_packer* cmpr, int nr_samples_to_encode, int Channels, std::vector<char>& data_stream, int BYTESPERSAMPLE)
{
    /** Create a [Channels x nr_samples_to_encode] matrix to store the original data in an int32_t matrix */
    tensor_i32 orig;
    orig.resize(Channels, nr_samples_to_encode);
    convert_i24native_to_i32(orig.d2d, (uint8_t*)data_stream.data(), nr_samples_to_encode, Channels, BYTESPERSAMPLE, false);

    /** allocate a space to compress the data in, then compress the input data */
    size_t dst_max_len = nr_samples_to_encode * Channels * BYTESPERSAMPLE * 2;
    size_t dst_len_;
    unsigned char dst[dst_max_len];
    cmpr->compress((uint8_t*)data_stream.data(), dst, dst_max_len, dst_len_);
    cout << "dst_len_: " << dst_len_ << endl;

    /** Allocat space for decompression, then decompress the compressed data. */
    size_t compressed_len;
    unsigned char decdst[dst_max_len];
    cmpr->decompress(dst, compressed_len, decdst);
    cout << "compressed len: " << compressed_len << "   COMPRESSION CR = " << (double)(Channels * BYTESPERSAMPLE * nr_samples_to_encode) / compressed_len << std::endl;

    /** Write the decoded data to the disc. In case of lossless compression this must match the input data */
    write_buffer_("decoded.bin", decdst, data_stream.size());

    /** Create a [Channels x nr_samples_to_encode] matrix to store the decoded data in an int32_t matrix */
    tensor_i32 enc;
    enc.resize(Channels, nr_samples_to_encode);
    convert_i24native_to_i32(enc.d2d, decdst, nr_samples_to_encode, Channels, BYTESPERSAMPLE, false);

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

void test_1()
{
    std::vector<char> data_stream;
    int BYTESPERSAMPLE = 3;
    int Channels = 3;
    bool filter_data = true;
    /** "data_stream.bin": ECG data sampled @ 2000Sps, 3 Channels and 24 bit resolution, each sample stored in 3 bytes.
    * Total number of samples: 60.000, 20.000 Samples per channel.
    * Structure: B0: most significant byte, B2: less significant byte
    *    - [CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2][CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2] ...
    *    - see the function convert_i24native_to_i32() */
    if (read_buffer_("data_stream.bin", data_stream))
    {
        int nrsamples = data_stream.size() / (Channels * BYTESPERSAMPLE);
        if (filter_data)
        {
            /** We convert the native data to a [Channels x nrsamples] sized matrix. */
            tensor_i32 enc;
            enc.resize(Channels, nrsamples);
            convert_i24native_to_i32(enc.d2d, (uint8_t*)data_stream.data(), nrsamples, Channels, BYTESPERSAMPLE, false);
            /** Bellow filter coefficients represent a bandpass butterworth filter of 0.4-200Hz @ 2000Sps */
            double n[] = {1.00000000000, -3.14332095199, 3.70064088865, -1.97083923944, 0.41351972908};
            double d[] = {0.06722876941, 0.00000000000, -0.13445753881, 0.00000000000, 0.06722876941};
            /** Create the filter object and use it on the data */
            i_filter* filter = i_filter::new_iir(n, d, 5);
            for (int j = 0; j < Channels; ++j)
            {
                filter->init_history_values(enc.d2d[j][0]);
                for (int i = 0; i < nrsamples; ++i)
                    enc.d2d[j][i] = filter->filter_opt(enc.d2d[j][i]);
            }
            /** We convert the [Channels x nrsamples] sized matrix back to native data. */
            convert_i32_to_i24native((uint8_t*)data_stream.data(), enc.d2d, nrsamples, Channels, BYTESPERSAMPLE, nrsamples, false);
        }

        /** For the sake of efficient compression, packers needs to know about the structure of native data.
        * this is why not a simple [size] is given as an argument for compression, but [BYTESPERSAMPLE, Channels, nrsamples] */
        cout << "Testing packers. xdelta_hzr:" << endl;
        i_signal_packer* cmpr = i_signal_packer::new_xdelta_hzr(BYTESPERSAMPLE, Channels, nrsamples);
        test_packer_(cmpr, nrsamples, Channels, data_stream, BYTESPERSAMPLE);

        /** For transformation based compression methods we need to provide a number of samples of 2^N */
        cout << "hadamard:" << endl;
        cmpr = i_signal_packer::new_hadamard(BYTESPERSAMPLE, Channels, 4096 * 4);
        test_packer_(cmpr, 4096 * 4, Channels, data_stream, BYTESPERSAMPLE);

        /** For dct we only provide fewer samples to encode as it is slow */
        cout << "dct:" << endl;
        cmpr = i_signal_packer::new_dct(BYTESPERSAMPLE, Channels, 4096 * 2);
        test_packer_(cmpr, 4096 * 2, Channels, data_stream, BYTESPERSAMPLE);
        cout << "EO Testing packers." << endl;
    }
}

int main()
{
    test_1();
    return 0;
}
