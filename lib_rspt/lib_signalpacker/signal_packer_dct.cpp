/***************************************************************************
* Copyright 2024 Tamas Levente Kis - tamkis@gmail.com                      *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <numeric>
#include <vector>
#include <iostream>
#include <math.h>

#include "../signal_packer.h"
#include "../lib_zaxtensor/ZaxJsonParser.h"
#include "../lib_zaxtensor/ZaxTensor.h"
#include "signal_packer_base.h"
#include "utils.h"

using namespace std;

class signal_packer_dct: public signal_packer_int32_base
{
    const double PI = 3.14159265358979323846;
    const double quality = 128.0; /// 1 denotes the best possible quality
    size_t nr_of_channels_;
    size_t bytes_per_channel_;
    size_t nr_of_samples_in_each_channel_;
    tensor_f32 COSINES;
    tensor_f32 Cs;
    tensor_i32 dct;
    unsigned int nr_bytes_to_compress_ = 2;

public:
    signal_packer_dct(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel)
        : nr_of_channels_(nr_of_channels),
          bytes_per_channel_(bytes_per_channel),
          nr_of_samples_in_each_channel_(nr_of_samples_in_each_channel)
    {
        enc_.resize(nr_of_channels_, nr_of_samples_in_each_channel_);
        dct.resize(nr_of_channels_, nr_of_samples_in_each_channel_);
        serialized_.resize(nr_bytes_to_compress_, enc_.d2 * enc_.d1);
        init_cos_table(enc_.d2);
    }

    void init_cos_table(int n)
    {
        COSINES.resize(n, n);
        Cs.resize(n);
        double pi_n_2 = PI / (n * 2.0);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                COSINES.d2d[i][j] = cos(((i << 1) * j + j) * pi_n_2);
            if (i)
                Cs.d1d[i] = 1;
            else
                Cs.d1d[i] = 1 / sqrt(2);
        }
    }

    void transform_dct(int n, const int *src, int *dst)
    {
        double ratio1 = sqrt(2.0 / n);
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int x = 0; x < n; x++)
                sum += src[x] * COSINES.d2d[x][i];
            sum *= Cs.d1d[i] * ratio1 / quality;
            dst[i] = sum;
        }
    }

    void transform_idct(int n, const int *dct, int *idct)
    {
        double ratio1 = sqrt(2.0 / n);
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int x = 0; x < n; x++)
                sum += Cs.d1d[x] * dct[x] * COSINES.d2d[i][x];
            sum *= ratio1 * quality;
            idct[i] = sum;
        }
    }

    virtual void compress(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len)
    {
        convert_i24native_to_i32(enc_.d2d, src, nr_of_samples_in_each_channel_, nr_of_channels_, bytes_per_channel_, false);
        int32_t means[enc_.d1];
        for (int i = 0; i < enc_.d1; ++i)
        {
            means[i] = average_32(enc_.d2d[i], enc_.d2);
            offset_32(enc_.d2d[i], enc_.d2, -means[i]);
        }
        for (int j = 0; j < enc_.d1; ++j)
        {
            transform_dct(enc_.d2, enc_.d2d[j], dct.d2d[j]);
            for (int i = 0; i < enc_.d2; ++i)
                enc_.d2d[j][i] = dct.d2d[j][i];
        }
        delta_encode(enc_.data(), enc_.d1 * enc_.d2);
        offset_32(enc_.data(), enc_.m_2d->d1 * enc_.m_2d->d2, -128);
        xor_encode_32(enc_.data(), enc_.d1 * enc_.d2);
        unsigned char header[enc_.d1 * 3];
        for (int i = 0; i < enc_.d1; ++i)
        {
            header[i * 3 + 0] = means[i] & 0x000000FF;
            header[i * 3 + 1] = (means[i] & 0x0000FF00) >> 8;
            header[i * 3 + 2] = (means[i] & 0x00FF0000) >> 16;
        }
        compress_i32(src, dst, dst_max_len, dst_len, 1, nr_bytes_to_compress_, header, enc_.d1 * 3);
    }

    virtual void decompress(const unsigned char* src, size_t& src_len, unsigned char* dst)
    {
        uint8_t compression_method;
        unsigned char header[enc_.d1 * 3];
        decompress_i32(src, src_len, dst, compression_method, nr_bytes_to_compress_, header, enc_.d1 * 3);
        if (compression_method != 1)
            cout << "ERROR: compression method unsupported." << endl;
        xor_decode_32(enc_.data(), enc_.d1 * enc_.d2);
        offset_32(enc_.data(), enc_.m_2d->d1 * enc_.m_2d->d2, 128);
        delta_decode(enc_.data(), enc_.d1 * enc_.d2, 0);
        for (int j = 0; j < enc_.d1; ++j)
        {
            transform_idct(enc_.d2, enc_.d2d[j], dct.d2d[j]);
            for (int i = 0; i < enc_.d2; ++i)
                enc_.d2d[j][i] = dct.d2d[j][i];
        }
        int32_t means[enc_.d1];
        for (int i = 0; i < enc_.d1; ++i)
            means[i] = ((header[i * 3] | (header[i * 3 + 1] << 8) | (header[i * 3 + 2] << 16)) << 8) >> 8; /// << 8 >> 8 handles the case when bytes 1,2,3 are representing a signed integer.
        for (int i = 0; i < enc_.d1; ++i)
            offset_32(enc_.d2d[i], enc_.d2, means[i]);
        convert_i32_to_i24native(dst, enc_.d2d, nr_of_samples_in_each_channel_, nr_of_channels_, bytes_per_channel_, false);
    }

    virtual ~signal_packer_dct() = default;
};

i_signal_packer* i_signal_packer::new_dct(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel)
{
    return new signal_packer_dct(bytes_per_channel, nr_of_channels, nr_of_samples_in_each_channel);
}

void i_signal_packer::delete_dct(i_signal_packer* instance)
{
    delete((signal_packer_dct*)instance);
}
