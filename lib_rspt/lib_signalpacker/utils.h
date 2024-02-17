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

void convert_i24native_to_i32(int32_t** dst, const uint8_t* native, int nr_samples, int NRCHANNELS = 3, int BYTESPERCHANNEL = 3, bool reverse_byte_order = true);
void convert_i32_to_i24native(uint8_t* native, int32_t** src, int nr_samples, int NRCHANNELS, int BYTESPERCHANNEL, int NR_SAMPLES_PER_PACKET_PER_CHANNEL, bool reverse_byte_order = true);
void delta_encode(int32_t* arr, size_t len);
void delta_decode(int32_t* arr, size_t len, int32_t min_val);
void offset_32(int32_t* arr, size_t len, int32_t val);
int32_t average_32(int32_t* arr, size_t len);
double average_32(double* arr, size_t len);
void xor_encode_32(int32_t* arr, size_t len);
void xor_decode_32(int32_t* arr, size_t len);
