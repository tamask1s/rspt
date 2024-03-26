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

#include "../signal_packer.h"
#include "../lib_zaxtensor/ZaxJsonParser.h"
#include "../lib_zaxtensor/ZaxTensor.h"
#include "signal_packer_base.h"
#include "utils.h"

using namespace std;

class signal_packer_xdelta_hzr: public signal_packer_int32_base
{
    size_t nr_of_channels_;
    size_t bytes_per_channel_;
    size_t nr_of_samples_in_each_channel_;
    unsigned int nr_bytes_to_compress_;

public:
    signal_packer_xdelta_hzr(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel, size_t nr_bytes_to_encode)
        : nr_of_channels_(nr_of_channels),
          bytes_per_channel_(bytes_per_channel),
          nr_of_samples_in_each_channel_(nr_of_samples_in_each_channel),
          nr_bytes_to_compress_(nr_bytes_to_encode)
    {
        enc_.resize(nr_of_channels_, nr_of_samples_in_each_channel_);
        serialized_.resize(nr_bytes_to_compress_, enc_.d2 * enc_.d1);
    }

    virtual void compress(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len)
    {
        convert_native_to_i32(enc_.d2d, src, nr_of_samples_in_each_channel_, nr_of_channels_, bytes_per_channel_, false);
        delta_encode(enc_.data(), enc_.d1 * enc_.d2);
        offset_32(enc_.data(), enc_.m_2d->d1 * enc_.m_2d->d2, -128);
        xor_encode_32(enc_.data(), enc_.d1 * enc_.d2);
        compress_i32(src, dst, dst_max_len, dst_len, 0, nr_bytes_to_compress_);
        unsigned char* decompressed = new unsigned char[nr_of_samples_in_each_channel_ * nr_of_channels_ * bytes_per_channel_];
        size_t src_len;
        decompress(dst, src_len, decompressed);
        if (memcmp(src, decompressed, nr_of_samples_in_each_channel_ * nr_of_channels_ * bytes_per_channel_))
        {
            delete[] decompressed;
            cout << "Compression needs one more byte to encode." << endl;
            nr_bytes_to_compress_++;
            serialized_.resize(nr_bytes_to_compress_, enc_.d2 * enc_.d1);
            compress(src, dst, dst_max_len, dst_len);
        }
        else
            delete[] decompressed;
    }

    virtual int decompress(const unsigned char* src, size_t& src_len, unsigned char* dst)
    {
        uint8_t compression_method;
        decompress_i32(src, src_len, dst, compression_method, nr_bytes_to_compress_);
        if (compression_method != 0)
            cout << "ERROR: compression method unsupported." << endl;
        xor_decode_32(enc_.data(), enc_.d1 * enc_.d2);
        offset_32(enc_.data(), enc_.m_2d->d1 * enc_.m_2d->d2, 128);
        delta_decode(enc_.data(), enc_.d1 * enc_.d2, 0);
        convert_i32_to_native(dst, enc_.d2d, nr_of_samples_in_each_channel_, nr_of_channels_, bytes_per_channel_, false);
        return 0;
    }

    virtual ~signal_packer_xdelta_hzr() = default;
};

i_signal_packer* i_signal_packer::new_xdelta_hzr(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel, size_t nr_bytes_to_encode)
{
    return new signal_packer_xdelta_hzr(bytes_per_channel, nr_of_channels, nr_of_samples_in_each_channel, nr_bytes_to_encode);
}

void i_signal_packer::delete_xdelta_hzr(i_signal_packer* instance)
{
    delete((signal_packer_xdelta_hzr*)instance);
}
