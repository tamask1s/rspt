#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <numeric>
#include <vector>
#include <functional>
#include <iostream>
#include <chrono>
#include <math.h>
#include <random>
#include <set>
#include <deque>

/** lib_rspt interfaces */
#include "../lib_rspt/signal_packer.h"
#include "../lib_rspt/filter.h"

/** Bellow includes are not necessary to use the rspt lib, but the functionalities can be used to easy data manipulation. */
#include "../lib_rspt/lib_signalpacker/utils.h"
#include "../lib_rspt/lib_zaxtensor/ZaxJsonParser.h"
#include "../lib_rspt/lib_zaxtensor/ZaxTensor.h"
#include "../lib_rspt/lib_ring_buffer/ring_buffers.h"
#include "../lib_rspt/lib_stat/rolling_window_median.h"

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

void test_packer_(i_signal_packer* cmpr, int nr_samples_to_encode, int Channels, uint8_t* data_stream, int bytes_per_sample)
{
    write_buffer_("_original.bin", data_stream, nr_samples_to_encode * Channels * bytes_per_sample);
    /** Create a [Channels x nr_samples_to_encode] matrix to store the original data in an int32_t matrix */
    tensor_i32 orig;
    orig.resize(Channels, nr_samples_to_encode);
    convert_native_to_i32(orig.d2d, data_stream, nr_samples_to_encode, Channels, bytes_per_sample, false);

    /** allocate a space to compress the data in, then compress the input data */
    size_t dst_max_len = nr_samples_to_encode * Channels * bytes_per_sample * 2;
    size_t dst_len_;
    unsigned char* dst = new unsigned char[dst_max_len];
    auto start = std::chrono::system_clock::now();
    cmpr->compress(data_stream, dst, dst_max_len, dst_len_);
    cout << "compression finished. compressed len: " << dst_len_ << " time elapsed: " << (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start)).count() << " ms." << endl;

    /** Allocat space for decompression, then decompress the compressed data. */
    size_t compressed_len;
    unsigned char* decdst = new unsigned char[dst_max_len];
    start = std::chrono::system_clock::now();
    if (cmpr->decompress(dst, compressed_len, decdst) != 0)
    {
        cout << "WARNING: decompression was not successful." << endl;
        compressed_len = dst_len_;
        memcpy(decdst, data_stream, nr_samples_to_encode * Channels * bytes_per_sample);
    }

    delete[] dst;
    cout << "decomp finished. " << " time:" << (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start)).count() << "ms. orig len: " << Channels * bytes_per_sample * nr_samples_to_encode << " compressed len: " << compressed_len << "   COMPRESSION CR = " << (double)(Channels * bytes_per_sample * nr_samples_to_encode) / compressed_len;

    /** Write the decoded data to the disc. In case of lossless compression this must match the input data */
    write_buffer_("_decoded.bin", decdst, nr_samples_to_encode * Channels * bytes_per_sample);

    /** Create a [Channels x nr_samples_to_encode] matrix to store the decoded data in an int32_t matrix */
    tensor_i32 enc;
    enc.resize(Channels, nr_samples_to_encode);
    convert_native_to_i32(enc.d2d, decdst, nr_samples_to_encode, Channels, bytes_per_sample, false);

    delete[] decdst;

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
    cout << " PRDN[%] = " << sqrt(MSE / origg) * 100.0 << endl;
}

void test_data(uint8_t* data_stream, int bytes_per_sample, int nr_channels, int nr_samples, bool filter_data, size_t nr_bytes_to_encode)
{
    if (filter_data)
    {
        /** We convert the native data to a [nr_channels x nr_samples] sized matrix. */
        tensor_i32 enc;
        enc.resize(nr_channels, nr_samples);
        convert_native_to_i32(enc.d2d, data_stream, nr_samples, nr_channels, bytes_per_sample, false);

        /** Bellow filter coefficients represent a bandpass butterworth filter of 0.4-200Hz @ 2000Sps */
        double n[] = {1.00000000000, -3.14332095199, 3.70064088865, -1.97083923944, 0.41351972908};
        double d[] = {0.06722876941, 0.00000000000, -0.13445753881, 0.00000000000, 0.06722876941};
        /** Create the filter object and use it on the data */
        i_filter* filter = i_filter::new_iir(n, d, 5);
        for (int j = 0; j < nr_channels; ++j)
        {
            filter->init_history_values(enc.d2d[j][0], 2000);
            for (int i = 0; i < nr_samples; ++i)
                enc.d2d[j][i] = filter->filter_opt(enc.d2d[j][i]);
        }
        /** We convert the [nr_channels x nr_samples] sized matrix back to native data. */
        convert_i32_to_native(data_stream, enc.d2d, nr_samples, nr_channels, bytes_per_sample, false);
    }

    /** Initialize packer. For the sake of efficiency, packers needs to know about the structure of native data. */
    /** This is why not a simple [size] is given as an argument, but [bytes_per_sample, nr_channels, nr_samples] */
    cout << "******************************* Testing packers *******************************\n----------------------xdelta_hzr:----------------------" << endl;
    i_signal_packer* cmpr = i_signal_packer::new_xdelta_hzr(bytes_per_sample, nr_channels, nr_samples, nr_bytes_to_encode);
    test_packer_(cmpr, nr_samples, nr_channels, data_stream, bytes_per_sample);

    /** For transformation based compression methods we need to provide a number of samples of 2^N, therefore we truncate the data @ 16384 samples */
    cout << "----------------------hadamard:----------------------" << endl;
    cmpr = i_signal_packer::new_hadamard(bytes_per_sample, nr_channels, 16384);
    test_packer_(cmpr, 16384, nr_channels, data_stream, bytes_per_sample);

    /** For dct we only provide fewer samples to encode as it is slow, therefore we truncate the data @ 4096 samples */
    cout << "----------------------dct:----------------------" << endl;
    cmpr = i_signal_packer::new_dct(bytes_per_sample, nr_channels, 4096);
    test_packer_(cmpr, 4096, nr_channels, data_stream, bytes_per_sample);

    // another compression method
//    cout << "----------------------Lala:----------------------" << endl;
//    cmpr = i_signal_packer::new_lala(bytes_per_sample, nr_channels, nr_samples);
//    test_packer_(cmpr, nr_samples, nr_channels, data_stream, bytes_per_sample);
}

void test_1()
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of ECG data sampled @ 2000Sps, 3 Channels and 24 bit res 20.000 Samples per channel." << endl;
    int bytes_per_sample = 3;
    std::vector<char> data_stream;
    int nr_channels = 3;
    /** "data_stream.bin": ECG data sampled @ 2000Sps, 3 Channels and 24 bit resolution, each sample stored in 3 bytes.
    * Total number of samples: 60.000, 20.000 Samples per channel.
    * Structure: B0: most significant byte, B2: less significant byte
    *    - [CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2][CH0 B0][CH0 B1][CH0 B2][CH1 B0][CH1 B1][CH1 B2][CH2 B0][CH2 B1][CH2 B2] ...
    *    - see the function convert_native_to_i32() */
    if (read_buffer_("data_stream.bin", data_stream))
    {
        int nr_samples = data_stream.size() / (nr_channels * bytes_per_sample);
        cout << "data_stream.bin" << " file loaded: " << nr_samples << endl;
        test_data((uint8_t*)data_stream.data(), bytes_per_sample, nr_channels, nr_samples, false, 3);
    }
}

void test_2()
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of a sine wave, 1 Channels and 32 bit res 16384 Samples." << endl;
    const int bytes_per_sample = 4;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 32 bit resolution, stored in an array of int32_t. */
    /** Total number of samples: 16384, total data size: 16384 * 4 Bytes. */
    int32_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;
    test_data((uint8_t*)data_stream, bytes_per_sample, nr_channels, nr_samples, true, 3);
}

void test_3()
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of a sine wave, 1 Channels and 16 bit res 16384 Samples." << endl;
    const int bytes_per_sample = 2;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 16 bit resolution, stored in an array of int16_t. */
    /** Total number of samples: 16384, total data size: 32768 Bytes. */
    int16_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 1000.0;
    test_data((uint8_t*)data_stream, bytes_per_sample, nr_channels, nr_samples, true, 3);
}

void test_4()
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of a sine wave, 1 Channels and 8 bit res 16384 Samples." << endl;
    const int bytes_per_sample = 1;
    const int nr_samples = 16384;
    const int nr_channels = 1;
    /** Simulate simple sinusoidal data. 1 Channels and 8 bit resolution, stored in an array of int8_t. */
    /** Total number of samples: 16384, total data size: 16384 Bytes. */
    int8_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i / 100.0) * 100.0;
    test_data((uint8_t*)data_stream, bytes_per_sample, nr_channels, nr_samples, true, 3);
}

void test_5()
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of a sine wave compressed with xdelta_hzr method only, 1 Channels and 32 bit res 8192 Samples." << endl;
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
    i_signal_packer* c = i_signal_packer::new_xdelta_hzr(bytes_per_sample, nr_channels, nr_samples, 3);

    /** allocate sufficient room for compressed data, then compress the data */
    size_t dst_max_len = nr_samples * nr_channels * bytes_per_sample * 2;
    unsigned char dst[dst_max_len];
    size_t compressed_size;
    c->compress((uint8_t*)data_stream, dst, dst_max_len, compressed_size);
    std::cout << "compressed_size: " << compressed_size;

    /** Allocate space for decompression, then decompress the compressed data. */
    size_t cmpr_size;
    unsigned char decdst[dst_max_len];
    c->decompress(dst, cmpr_size, decdst);
    std::cout << "\n\n Compressing simple sine wave\n  compressed len: " << cmpr_size << " compression CR = ";
    std::cout << (double)(nr_channels * bytes_per_sample * nr_samples) / cmpr_size << std::endl;
}

void test_6()
{
    /** cout << "\nFiltering of a simulated sine wave combined of a 4Hz and a 70Hz wave @2kSps, 1 Channels and 32 bit res 8192 Samples." << endl; */
    /** Simulate sinusoidal data. 1 Channels and 32 bit resolution, stored in an */
    /** array of int32_t. Signal with 2 sinusoids: 4Hz and 70Hz combined. */
    const int nr_samples = 8192;
    const double sample_rate = 2000;
    const double sr_pi = 3.14159265358979323846 * 2.0 / sample_rate;
    int32_t data_stream[nr_samples];
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i * sr_pi* 4.0) * 1000.0 + sin(i * sr_pi * 70.0) * 1000.0;

    double n[] = {1.00000000000, -1.97778648378, 0.97803050849}; /// LP 5Hz @ 2kSps
    double d[] = {0.00006100618, 0.00012201236, 0.00006100618};
    i_filter* lp_filter = i_filter::new_iir(n, d, 3);
    lp_filter->init_history_values(data_stream[0], sample_rate);
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = lp_filter->filter_opt(data_stream[i]);

    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = sin(i * sr_pi * 4.0) * 1000.0 + sin(i * sr_pi * 70.0) * 1000.0;

    double n2[] = {1.00000000000, -1.77863177782, 0.80080264667}; /// HP 50Hz @ 2kSps
    double d2[] = {0.89485860612, -1.78971721225, 0.89485860612};
    i_filter* hp_filter = i_filter::new_iir(n2, d2, 3);
    hp_filter->init_history_values(data_stream[0], sample_rate);
    for (int i = 0; i < nr_samples; ++i)
        data_stream[i] = hp_filter->filter_opt(data_stream[i]);
}

void test_7(const char* filename, size_t nr_bytes_to_encode)
{
    cout << "\n*******************************************************************************" << endl;
    cout << "Cpmression of ECG data sampled @ 1000Sps, 12 Channels and 32 bit res 1.801.625 Samples per channel." << endl;
    int bytes_per_sample = 4;
    std::vector<char> data_stream;
    int nr_channels = 12;
    /** "12_chan_32bit_1801625_samples.bin": ECG data sampled @ 1000Sps, 12 Channels and 32 bit resolution, each sample stored in 4 bytes.
    * Structure: B0: most significant byte, B2: less significant byte
    *    - [CH0 B0][CH0 B1][CH0 B2][CH0 B3][CH1 B0][CH1 B1][CH1 B2][CH1 B3] ... */
    if (read_buffer_(filename, data_stream))
    {
        int nr_samples = data_stream.size() / (nr_channels * bytes_per_sample);
        cout << filename << " file loaded: " << nr_samples << endl;
        test_data((uint8_t*)data_stream.data(), bytes_per_sample, nr_channels, nr_samples, false, nr_bytes_to_encode);
    }
}

//#include "../lib_rspt/lib_cl/reader.h"
//void convert_raw_to_bin(const char* filename)
//{
//    auto inputChannels = ReadFile2(filename);
//    size_t sampleCount = inputChannels[0].size();
//    size_t channelCount = inputChannels.size();
//    std::cout << filename << "file loaded: " << channelCount << " channels with " << sampleCount << " samples\n";
//    int32_t* flat = new int32_t[channelCount * sampleCount];
//    for (size_t j = 0; j < sampleCount; ++j)
//    {
//        for (size_t i = 0; i < channelCount; ++i)
//        {
//            flat[j * channelCount + i] = inputChannels[i][j];
//        }
//    }
//    string out_filename = std::to_string(channelCount) + "_chan_32bit_" + std::to_string(sampleCount) + "_samples_" + filename + ".bin";
//    cout << "writing: " << out_filename << endl;
//    write_buffer_(out_filename.c_str(), (uint8_t*)flat, sizeof(int32_t) * channelCount * sampleCount);
//    delete[] flat;
//}

template<const size_t k, const size_t n>
void rolling_window_median_tester(std::vector<double> expected_result)
{
    std::vector<double> data(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 100);

    for (size_t i = 0; i < n; ++i)
        data[i] = dis(gen);

    std::vector<double> nums = {1, 2, 3, 4, 5, 6, 7, 8, 4, 5, 6, 5, 4, 3, 2, 1, 1, 1, 1, 9};
    for (size_t i = 0; i < nums.size(); ++i)
        data[i] = nums[i];

    std::vector<double> result;
    result.reserve(n);

    rolling_window_median<double> rwm(k);

    auto start = std::chrono::high_resolution_clock::now();

    for (double num : data)
        result.push_back(rwm.insert(num));

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    cout << "result:   ";
    for (unsigned int i = 0; i < expected_result.size(); ++i)
        std::cout << result[i] << " ";

    cout << endl << "expected: ";
    for (unsigned int i = 0; i < expected_result.size(); ++i)
        std::cout << expected_result[i] << " ";

    bool equal = true;
    for (unsigned int i = 0; i < expected_result.size(); ++i)
        if (result[i] != expected_result[i])
            equal = false;

    cout << endl << "equal: " << (equal ? "true" : "false") << endl << endl;
}

void test_8_rolling_window_median()
{
    rolling_window_median_tester<5, 1000000>({1, 1.5, 2, 2.5, 3, 4, 5, 6, 6, 6, 6, 5, 5, 5, 4, 3, 2, 1, 1, 1});
    rolling_window_median_tester<6, 1000000>({1, 1.5, 2, 2.5, 3, 3.5, 4.5, 5.5, 5.5, 5.5, 6, 5.5, 5, 4.5, 4.5, 3.5, 2.5, 1.5, 1, 1});
    rolling_window_median_tester<7, 1000000>({1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 5, 5, 6, 6, 5, 5, 4, 4, 3, 2, 1, 1});
    rolling_window_median_tester<1500, 1000000>({1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 4, 4.5, 5, 5, 5, 4.5, 4, 4, 4, 4, 4, 4});

//    rolling_window_median_tester<3, 50000000>({});
//    rolling_window_median_tester<5, 50000000>({});
//    rolling_window_median_tester<10, 50000000>({});
//    rolling_window_median_tester<100, 50000000>({});
//    rolling_window_median_tester<1000, 50000000>({});
//    rolling_window_median_tester<10000, 50000000>({});
//    rolling_window_median_tester<100000, 50000000>({});
//    rolling_window_median_tester<177800, 50000000>({});
//    rolling_window_median_tester<316200, 50000000>({});
//    rolling_window_median_tester<562300, 50000000>({});
//    rolling_window_median_tester<1000000, 50000000>({});
//    rolling_window_median_tester<5000000, 50000000>({});
//    rolling_window_median_tester<10000000, 50000000>({});
//    rolling_window_median_tester<15000000, 50000000>({});
//    rolling_window_median_tester<20000000, 50000000>({});
//    rolling_window_median_tester<30000000, 50000000>({});

}

int main()
{
//    test_1();
//    test_2();
//    test_3();
//    test_4();
//    test_5();
//    test_6();
//    test_7("12_chan_32bit_34199_samples_r00000135fghd8.raw.bin", 1);
    test_8_rolling_window_median();
//    convert_raw_to_bin("r000000b520wf2.raw");
//    convert_raw_to_bin("r000000k54yy4m.raw");
//    convert_raw_to_bin("r000000v5qk36w.raw");
//    convert_raw_to_bin("r00000134vdjb6.raw");
//    convert_raw_to_bin("r00000135fghd8.raw");
//    convert_raw_to_bin("r000001b57e2n8.raw");

//    test_7("12_chan_32bit_115225_samples_r00000134vdjb6.raw.bin", 3);
//    test_7("12_chan_32bit_1801625_samples_r000001b57e2n8.raw.bin", 3);
//    test_7("12_chan_32bit_1801853_samples_r000000k54yy4m.raw.bin", 4);
//    test_7("12_chan_32bit_34199_samples_r00000135fghd8.raw.bin", 3);
//    test_7("12_chan_32bit_56120_samples_r000000b520wf2.raw.bin", 4);
//    test_7("12_chan_32bit_68300_samples_r000000v5qk36w.raw.bin", 4);

//    test_7("12_chan_32bit_115225_samples_r00000134vdjb6.raw.bin", 1);
//    test_7("12_chan_32bit_1801625_samples_r000001b57e2n8.raw.bin", 1);
//    test_7("12_chan_32bit_1801853_samples_r000000k54yy4m.raw.bin", 1);
//    test_7("12_chan_32bit_34199_samples_r00000135fghd8.raw.bin", 1);
//    test_7("12_chan_32bit_56120_samples_r000000b520wf2.raw.bin", 1);
//    test_7("12_chan_32bit_68300_samples_r000000v5qk36w.raw.bin", 1);

//    test_7("12_chan_32bit_115225_samples_r00000134vdjb6.raw.bin", 3);
//    test_7("12_chan_32bit_1801625_samples_r000001b57e2n8.raw.bin", 3);
//    test_7("12_chan_32bit_1801853_samples_r000000k54yy4m.raw.bin", 4);
//    test_7("12_chan_32bit_34199_samples_r00000135fghd8.raw.bin", 2);
//    test_7("12_chan_32bit_56120_samples_r000000b520wf2.raw.bin", 4);
//    test_7("12_chan_32bit_68300_samples_r000000v5qk36w.raw.bin", 4);
    return 0;
}
