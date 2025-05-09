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
#include <climits>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <cstdio>
#include <map>
#include <vector>
#include <limits>
#include <numeric>

//#include "../lib_lz4/lz4.h"
//#include "../lib_lz4/lz4hc.h"
//#include "../lib_zzip/zzlib.h"
#include "../signal_packer.h"
#include "../lib_zaxtensor/ZaxJsonParser.h"
#include "../lib_zaxtensor/ZaxTensor.h"
#include "../lib_hzr/libhzr.h"
#include "signal_packer_base.h"

using namespace std;
#define CHUNK_MAX_SIZE_TYPE uint32_t

void signal_packer_int32_base::compress_i32(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len, uint8_t compression_method, unsigned int nr_bytes_to_compress, const unsigned char* header, size_t header_size)
{
    int bindx = -1;
    if (nr_bytes_to_compress == 4)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j)
            {
                serialized_.d2d[0][++bindx] = enc_.d2d[i][j] & 0x000000FF;
                serialized_.d2d[1][bindx] = (enc_.d2d[i][j] & 0x0000FF00) >> 8;
                serialized_.d2d[2][bindx] = (enc_.d2d[i][j] & 0x00FF0000) >> 16;
                serialized_.d2d[3][bindx] = (enc_.d2d[i][j] & 0xFF000000) >> 24;
            }
    else if (nr_bytes_to_compress == 3)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j)
            {
                serialized_.d2d[0][++bindx] = enc_.d2d[i][j] & 0x000000FF;
                serialized_.d2d[1][bindx] = (enc_.d2d[i][j] & 0x0000FF00) >> 8;
                serialized_.d2d[2][bindx] = (enc_.d2d[i][j] & 0x00FF0000) >> 16;
            }
    else if (nr_bytes_to_compress == 2)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j)
            {
                serialized_.d2d[0][++bindx] = enc_.d2d[i][j] & 0x000000FF;
                serialized_.d2d[1][bindx] = (enc_.d2d[i][j] & 0x0000FF00) >> 8;
            }
    else
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j)
                serialized_.d2d[0][++bindx] = enc_.d2d[i][j] & 0x000000FF;
    auto compress_bytes = [this](unsigned char*& dst, unsigned char* src, size_t& dst_max_len)
    {
        size_t out_len;
        hzr_encode(src, serialized_.d2, dst + sizeof(CHUNK_MAX_SIZE_TYPE), dst_max_len, &out_len);
        //out_len = LZ4_compress_default((const char*)src, (char*)dst + sizeof(CHUNK_MAX_SIZE_TYPE), serialized_.d2, dst_max_len);
        //out_len = LZ4_compress_HC((const char*)src, (char*)dst + sizeof(CHUNK_MAX_SIZE_TYPE), serialized_.d2, dst_max_len, 2);
        //memcpy(dst + sizeof(CHUNK_MAX_SIZE_TYPE), src, serialized_.d2);
        //out_len = ZzCompressBlock(dst + sizeof(CHUNK_MAX_SIZE_TYPE), serialized_.d2, 1, 1);
        //cout << "OUTLEN " << out_len << endl;
        *((CHUNK_MAX_SIZE_TYPE*)dst) = out_len;
        dst += sizeof(CHUNK_MAX_SIZE_TYPE) + out_len;
        dst_max_len -= sizeof(CHUNK_MAX_SIZE_TYPE) + out_len;
        return out_len;
    };
    *dst = compression_method;
    dst += 1;
    dst_len = 1;
    if (header && header_size)
    {
        memcpy(dst, header, header_size);
        dst += header_size;
        dst_len += header_size;
    }
    dst_max_len -= 1;
    dst_len += nr_bytes_to_compress * sizeof(CHUNK_MAX_SIZE_TYPE);
    for (unsigned int i = 0; i < nr_bytes_to_compress; ++i)
        dst_len += compress_bytes(dst, serialized_.d2d[i], dst_max_len);
}

void signal_packer_int32_base::decompress_i32(const unsigned char* src, size_t& src_len, unsigned char* dst, uint8_t& compression_method, unsigned int nr_bytes_to_compress, unsigned char* header, size_t header_size)
{
    const unsigned char* orig_src = src;
    compression_method = *src;
    auto decompress_bytes = [this](const unsigned char*& src, unsigned char* dst, size_t len)
    {
        size_t comp_len = *((CHUNK_MAX_SIZE_TYPE*)src);
        src += sizeof(CHUNK_MAX_SIZE_TYPE);
        hzr_decode(src, comp_len, dst, len);
        //int comp_len2 = LZ4_decompress_safe((char*)src, (char*)dst, comp_len, len);
        //memcpy(dst, src, comp_len);
        //int comp_len2 = ZzUncompressBlock(dst);
        src += comp_len;
    };
    src++;
    if (header && header_size)
    {
        memcpy(header, src, header_size);
        src += header_size;
    }
    serialized_.fill(0);
    for (unsigned int i = 0; i < nr_bytes_to_compress; ++i)
        decompress_bytes(src, serialized_.d2d[i], serialized_.d2);
    src_len = src - orig_src;
    int bindx = 0;
    if (nr_bytes_to_compress == 4)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j, ++bindx)
                enc_.d2d[i][j] = serialized_.d2d[0][bindx] | (serialized_.d2d[1][bindx] << 8) | (serialized_.d2d[2][bindx] << 16) | (serialized_.d2d[3][bindx] << 24);
    else if (nr_bytes_to_compress == 3)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j, ++bindx)
                enc_.d2d[i][j] = ((serialized_.d2d[0][bindx] | (serialized_.d2d[1][bindx] << 8) | (serialized_.d2d[2][bindx] << 16)) << 8) >> 8; /// << 8 >> 8 handles the case when bytes 1,2,3 are representing a signed integer.
    else if (nr_bytes_to_compress == 2)
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j, ++bindx)
                enc_.d2d[i][j] = ((serialized_.d2d[0][bindx] | (serialized_.d2d[1][bindx] << 8)) << 16) >> 16; /// << 16 >> 16 handles the case when bytes 1,2 are representing a signed integer.
    else
        for (int i = 0; i < enc_.d1; ++i)
            for (int j = 0; j < enc_.d2; ++j, ++bindx)
                enc_.d2d[i][j] = (serialized_.d2d[0][bindx] << 24) >> 24; /// << 24 >> 24 handles the case when byte 1 is representing a signed integer.
}

