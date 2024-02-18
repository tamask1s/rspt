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
#include <functional>
#include <math.h>

#include "utils.h"

using namespace std;

int32_t average_32(int32_t* arr, size_t len)
{
//    double average = arr[0];
//    for (size_t i = 1; i < len; ++i)
//        average = average * (i / (i + 1.0)) + (arr[i] / (i + 1.0));
    int64_t average = 0;
    for (size_t i = 0; i < len; ++i)
        average += arr[i];
    average /= len;
    return average;
}

double average_32(double * arr, size_t len)
{
    double average = 0;
    for (size_t i = 0; i < len; ++i)
        average += arr[i];
    average /= len;
    return average;
}

void convert_i32_to_native(uint8_t* native, int32_t** src, int nr_samples, int NRCHANNELS, int BYTESPERCHANNEL, bool reverse_byte_order)
{
    if (BYTESPERCHANNEL == 4)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0xFF000000) >> 24;
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x00FF0000) >> 16;
                    *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 3 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 3 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0xFF000000) >> 24;
                    *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x00FF0000) >> 16;
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
    }
    else if (BYTESPERCHANNEL == 3)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x00FF0000) >> 16;
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x00FF0000) >> 16;
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
    }
    else if (BYTESPERCHANNEL == 2)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x0000FF00) >> 8;
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
                }
    }
    else if (BYTESPERCHANNEL == 1)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL) = (src[i][s] & 0x000000FF);
    }
}

void convert_native_to_i32(int32_t** dst, const uint8_t* native, int nr_samples, int NRCHANNELS, int BYTESPERCHANNEL, bool reverse_byte_order)
{
    if (BYTESPERCHANNEL == 4)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    uint8_t tmp1 = *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp2 = *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp3 = *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp4 = *(native + 3 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint32_t tmp = ((tmp1 << 24) | (tmp2 << 16) | (tmp3 << 8) | tmp4);
                    dst[i][s] = *((int32_t*)&tmp);
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    dst[i][s] = *((int32_t*)(native + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL));
    }
    else if (BYTESPERCHANNEL == 3)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    uint8_t tmp1 = *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp2 = *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp3 = *(native + 2 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint32_t tmp = ((tmp1 << 16) | (tmp2 << 8) | tmp3) << 8;
                    dst[i][s] = *((int32_t*)&tmp) >> 8;
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    dst[i][s] = (*((int32_t*)(native + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL)) << 8) >> 8;
    }
    else if (BYTESPERCHANNEL == 2)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    uint8_t tmp1 = *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint8_t tmp2 = *(native + 1 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint32_t tmp = ((tmp1 << 8) | tmp2) << 16;
                    dst[i][s] = *((int32_t*)&tmp) >> 16;
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    dst[i][s] = (*((int32_t*)(native + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL)) << 16) >> 16;
    }
    else if (BYTESPERCHANNEL == 1)
    {
        if (reverse_byte_order)
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                {
                    uint8_t tmp1 = *(native + 0 + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL);
                    uint32_t tmp = tmp1 << 24;
                    dst[i][s] = *((int32_t*)&tmp) >> 24;
                }
        else
            for (int s = 0; s < nr_samples; ++s)
                for (uint16_t i = 0; i < NRCHANNELS; ++i)
                    dst[i][s] = (*((int32_t*)(native + NRCHANNELS * s * BYTESPERCHANNEL + i * BYTESPERCHANNEL)) << 24) >> 24;
    }
}

void delta_encode(int32_t* arr, size_t len)
{
    int32_t last = 0;
    for (unsigned int i = 0; i < len; ++i)
    {
        int32_t curr = arr[i];
        arr[i] = (curr - last);
        last = curr;
    }
}

void delta_decode(int32_t* arr, size_t len, int32_t min_val)
{
    int32_t last = 0;
    for (unsigned int i = 0; i < len; ++i)
    {
        int32_t delta = arr[i];
        arr[i] = (delta + last + min_val);
        last = arr[i];
    }
}

void offset_32(int32_t* arr, size_t len, int32_t val)
{
    for (unsigned int i = 0; i < len; ++i)
        arr[i] += val;
}

void xor_encode_32(int32_t* arr, size_t len)
{
    int32_t last = 0;
    for (unsigned int i = 0; i < len; ++i)
    {
        int32_t diff = last ^ arr[i];
        last = arr[i];
        arr[i] = diff;
    }
}

void xor_decode_32(int32_t* arr, size_t len)
{
    for (unsigned int i = 1; i < len; ++i)
        arr[i] = arr[i - 1] ^ arr[i];
}
