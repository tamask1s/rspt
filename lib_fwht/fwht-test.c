//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <map>
//#include <numeric>
//#include <vector>
//#include <functional>
//#include <math.h>
//
//#include "FifoBuffer.h"
//#include "ZaxJsonParser.h"
//#include "ZaxTensor.h"
//
//#include "fwht.h"
//
////#define Err() do {printf("%i Error\n", __LINE__); exit(1);} while (0)
////
////void test(void)
////{
////    {
////        int a[] = {1,0,1,0,0,1,1,0};
////        int n = 8;
////        int b[n];
////        fwht_transform(n, a, b);
////        if (b[0] != 4) Err();
////        if (b[1] != 2) Err();
////        if (b[2] != 0) Err();
////        if (b[3] != -2) Err();
////        if (b[4] != 0) Err();
////        if (b[5] != 2) Err();
////        if (b[6] != 0) Err();
////        if (b[7] != 2) Err();
////        fwht_transform(n, b, a);
////        fwht_normalize(n, a);
////        if (a[0] != 1) Err();
////        if (a[1] != 0) Err();
////        if (a[2] != 1) Err();
////        if (a[3] != 0) Err();
////        if (a[4] != 0) Err();
////        if (a[5] != 1) Err();
////        if (a[6] != 1) Err();
////        if (a[7] != 0) Err();
////    }
////
////    {
////        int a[] = {1,0,1,0,0,1,1,0,1,1,1,0,1,0,0,0};
////        int n = sizeof(a)/sizeof(*a);
////        int b[n];
////        fwht_transform(n, a, b);
////        fwht_transform(n, b, a);
////        fwht_normalize(n, a);
////        int ans[] = {1,0,1,0,0,1,1,0,1,1,1,0,1,0,0,0};
////        int i;
////        for (i = 0; i < n; i++) {
////            if (a[i] != ans[i]) Err();
////        }
////    }
////    {
////        int a[] = {1,0,1,0,0,1,1,0};
////        int b[] = {1,0,0,0,0,1,1,0};
////        int n = 8;
////        double diff_a = fwht_sum_absolute_difference(n, a, b);
////        double diff_b = fwht_sum_absolute_difference(n, b, a);
////        if (diff_a != diff_b) Err();
////
////        int a2[n];
////        int b2[n];
////        fwht_transform(n, a, a2);
////        fwht_transform(n, b, b2);
////        double diff_a2 = fwht_sum_absolute_difference(n, a2, b2);
////        double diff_b2 = fwht_sum_absolute_difference(n, b2, a2);
////        if (diff_a2 != diff_b2) Err();
////
////        fwht_transform(n, a2, a);
////        fwht_transform(n, b2, b);
////        fwht_normalize(n, a);
////        fwht_normalize(n, b);
////
////        double diff_a3 = fwht_sum_absolute_difference(n, a, b);
////        double diff_b3 = fwht_sum_absolute_difference(n, b, a);
////
////        if (diff_a != diff_a3) Err();
////        if (diff_b != diff_b3) Err();
////    }
////}
//
//int main(int argc, char *argv[])
//{
//    tensor_ui8 serialized;
//    std::vector<char> saved_stream;
//    int BYTESPERSAMPLE = 3;
//    int Channels = 3;
//    tensor_i32 enc;
//    tensor_i32 orig;
//    tensor_i32 downsampled;
//    tensor_i32 resampled;
//    int downsample_ratio = 8;
//    int shift_bits = 3;
//
//    if (read_buffer("something.mdebin", saved_stream))
//    {
//        int nrsamples = saved_stream.size() / (Channels * BYTESPERSAMPLE);
//        enc.resize(Channels, nrsamples);
//        orig.resize(Channels, nrsamples);
//        resampled.resize(Channels, nrsamples);
//        downsampled.resize(Channels, nrsamples / downsample_ratio);
//
//        for (int j = 0; j < Channels; ++j)
//            for (int i = 0; i < nrsamples; ++i)
//            {
//                uint32_t tmp = *((uint32_t*)(saved_stream.data() + Channels * i * BYTESPERSAMPLE + j * BYTESPERSAMPLE));
//                tmp <<= 8;
//                tmp = (int32_t)tmp;
//                tmp >>= 8;
//                enc.d2d[j][i] = tmp;
//                orig.d2d[j][i] = enc.d2d[j][i];
//            }
//
//
//    }
//    return 0;
//}
