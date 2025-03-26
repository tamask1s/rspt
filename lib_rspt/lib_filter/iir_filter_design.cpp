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
#include <inttypes.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#include "../filter.h"

bool create_filter_iir_butterworth(vector<double>& n, vector<double>& d, filter_type type, int order, double sampling_rate, double cutoff_low)
{
    if(order != 2 || (type != low_pass && type != high_pass) || sampling_rate <= 0 || cutoff_low <= 0)
        return false;

    // Prewarpolt cutoff: K = tan(pi*fc/fs)
    double K = tan(M_PI * cutoff_low / sampling_rate);
    double K2 = K * K;
    double sqrt2 = sqrt(2.0);

    // Kozos nevezo egyutthatok (mindkét esetben)
    double a0 = 1.0 + sqrt2 * K + K2;
    double a1 = 2.0 * (K2 - 1.0);
    double a2 = 1.0 - sqrt2 * K + K2;

    if(type == low_pass)
    {
        // Low pass: a számlálo analog egyutthatok
        double b0 = K2;
        double b1 = 2.0 * K2;
        double b2 = K2;
        n = { b0 / a0, b1 / a0, b2 / a0 };
    }
    else // high_pass
    {
        // High pass: analog high pass prototípus: H(s)= s^2/(s^2+sqrt2*K*s+K2)
        // Bilineáris transzformácio után:
        double b0 = 1.0;
        double b1 = -2.0;
        double b2 = 1.0;
        n = { b0 / a0, b1 / a0, b2 / a0 };
    }
    d = { 1.0, a1 / a0, a2 / a0 };

    return true;
}

//Ok butterworth LP, HP
//{
//    // Csak a 2. rendu alul- és felulátereszto szuro támogatása
//    if(order != 2 || (type != low_pass && type != high_pass) || sampling_rate <= 0 || cutoff_low <= 0)
//        return false;
//
//    // Prewarpolt cutoff: K = tan(pi*fc/fs)
//    double K = tan(M_PI * cutoff_low / sampling_rate);
//    double K2 = K * K;
//    double sqrt2 = sqrt(2.0);
//
//    // Közös nevezo egyutthatok (mindkét esetben)
//    double a0 = 1.0 + sqrt2 * K + K2;
//    double a1 = 2.0 * (K2 - 1.0);
//    double a2 = 1.0 - sqrt2 * K + K2;
//
//    if(type == low_pass)
//    {
//        // Low pass: a számlálo analog egyutthatok
//        double b0 = K2;
//        double b1 = 2.0 * K2;
//        double b2 = K2;
//        n = { b0 / a0, b1 / a0, b2 / a0 };
//    }
//    else // high_pass
//    {
//        // High pass: analog high pass prototípus: H(s)= s^2/(s^2+sqrt2*K*s+K2)
//        // Bilineáris transzformácio után:
//        double b0 = 1.0;
//        double b1 = -2.0;
//        double b2 = 1.0;
//        n = { b0 / a0, b1 / a0, b2 / a0 };
//    }
//    d = { 1.0, a1 / a0, a2 / a0 };
//
//    return true;
//}

// Ok butterworth LP
//{
//    // Most csak a 2. rendu alulátereszto szurot támogatjuk
//    if(order != 2 || type != low_pass || sampling_rate <= 0 || cutoff_low <= 0)
//        return false;
//
//    // Prewarpált cutoff
//    double K = tan(M_PI * cutoff_low / sampling_rate);
//    double K2 = K * K;
//    double sqrt2 = sqrt(2.0);
//
//    // Analog prototípus egyutthatoi
//    double a0 = 1.0 + sqrt2 * K + K2;
//    double a1 = 2.0 * (K2 - 1.0);
//    double a2 = 1.0 - sqrt2 * K + K2;
//
//    double b0 = K2;
//    double b1 = 2.0 * K2;
//    double b2 = K2;
//
//    // Normalizálás: osztás az a0-val, hogy d[0] = 1
//    n = { b0 / a0, b1 / a0, b2 / a0 };
//    d = { 1.0, a1 / a0, a2 / a0 };
//
//    return true;
//}



// --- Polinom muveletek --- //

// Két polinom szorzása: mindkettot magasfoktól konstans tagig tárolva.
vector<double> poly_multiply(const vector<double>& p, const vector<double>& q)
{
    vector<double> result(p.size() + q.size() - 1, 0.0);
    for (size_t i = 0; i < p.size(); i++)
    {
        for (size_t j = 0; j < q.size(); j++)
        {
            result[i+j] += p[i] * q[j];
        }
    }
    return result;
}

// Két polinom összeadása. Feltételezzük, hogy a két vektor azonos hosszúságú.
vector<double> poly_add(const vector<double>& p, const vector<double>& q)
{
    size_t n = max(p.size(), q.size());
    vector<double> r(n, 0.0);
    // Az eltéréseket nullákkal kiegészítjük.
    size_t p_offset = n - p.size();
    size_t q_offset = n - q.size();
    for (size_t i = 0; i < n; i++)
    {
        double a = (i < p_offset) ? 0.0 : p[i - p_offset];
        double b = (i < q_offset) ? 0.0 : q[i - q_offset];
        r[i] = a + b;
    }
    return r;
}

vector<double> poly_scale(const vector<double>& p, double s)
{
    vector<double> r = p;
    for (auto &coef : r)
        coef *= s;
    return r;
}

// --- Készít egy polinomot (z - 1)^n vagy (z + 1)^n --- //
// A polinomokat magasfoktól (z^n) a konstans tagig tároljuk.
vector<double> poly_z_minus_1(int n)
{
    vector<double> poly(n+1, 0.0);
    // Együtthatók: binom(n,k)*(-1)^k, k=0..n, legmagasabb fok: k=0
    for (int k = 0; k <= n; k++)
    {
        double coeff = 1.0;
        for (int i = 1; i <= k; i++)
        {
            coeff *= double(n - i + 1) / i;
        }
        poly[k] = coeff * ((k % 2 == 0) ? 1.0 : -1.0);
    }
    return poly;
}

vector<double> poly_z_plus_1(int n)
{
    vector<double> poly(n+1, 0.0);
    for (int k = 0; k <= n; k++)
    {
        double coeff = 1.0;
        for (int i = 1; i <= k; i++)
        {
            coeff *= double(n - i + 1) / i;
        }
        poly[k] = coeff;
    }
    return poly;
}

// --- Fo függvény --- //
// A cél: 2. rendu lowpass prototípusból (Butterworth) sávpass szuro készítése
// Az analóg filter átviteli függvényét a frekvencia-transzformáció: s -> (s^2+W0^2)/(Bw*s)
// és a bilineáris transzformáció segítségével digitalizáljuk.
// A végso digitalizált polinomokat úgy kell eloállítani, hogy a kimenet egy 4. rendu filter (5 koefficiens),
// melynek értékei megfelelnek a scipy által adottaknak.
bool create_filter_iir_butterworth_bandpass(vector<double>& n, vector<double>& d, filter_type type, int order, double sampling_rate, double cutoff_low, double cutoff_high)
{
    // Csak a band_pass és a 2. rendu prototípust támogatjuk.
    if (order != 2 || type != band_pass || sampling_rate <= 0 || cutoff_low <= 0 || cutoff_high <= cutoff_low)
        return false;

    double T = 1.0 / sampling_rate;
    double k = 2.0 / T; // bilineáris transzformáció: s = k*(z-1)/(z+1)

    // Prewarping: analóg cutoffok (Omega = k*tan(pi*f/fs))
    double Omega1 = k * tan(M_PI * cutoff_low / sampling_rate);
    double Omega2 = k * tan(M_PI * cutoff_high / sampling_rate);
    double Bw = Omega2 - Omega1;
    double W0 = sqrt(Omega1 * Omega2);

    // Analóg Butterworth lowpass prototípus (2. rendu) átviteli függvény: 1/(s^2+sqrt2*s+1)
    // A frekvencia-transzformáció miatt az analóg bandpass filter átviteli függvénye:
    // H(s) = [Bw^2 * s^2] / [s^4 + sqrt2*Bw*s^3 + (2W0^2+Bw^2)*s^2 + sqrt2*Bw*W0^2*s + W0^4]
    // Az analóg nevezo együtthatói:
    double a4 = 1.0;
    double a3 = sqrt(2.0) * Bw;
    double a2 = 2.0 * W0 * W0 + Bw * Bw;
    double a1 = sqrt(2.0) * Bw * W0 * W0;
    double a0 = W0 * W0 * W0 * W0;
    // Az analóg számláló: N(s) = Bw^2 * s^2.
    double b2 = Bw * Bw;

    // Most végezzük el a bilineáris transzformációt az analóg polinomokra.
    // Az analóg polinomokat s helyére k*(z-1)/(z+1) behelyettesítjük, majd megszorozzuk (z+1)^4-el,
    // hogy egy polinomot kapjunk z-ben.
    // Nevezo: D(z) = a4*k^4*(z-1)^4 + a3*k^3*(z-1)^3*(z+1) + a2*k^2*(z-1)^2*(z+1)^2 + a1*k*(z-1)*(z+1)^3 + a0*(z+1)^4.
    // Számláló: N(z) = b2*k^2*(z-1)^2*(z+1)^2.

    // Számoljuk ki az egyes tagokat polinomként.
    vector<double> zm1_4 = poly_z_minus_1(4); // (z-1)^4, 5 elem: z^4,...,z^0
    vector<double> zp1_4 = poly_z_plus_1(4);  // (z+1)^4, 5 elem
    vector<double> zm1_3 = poly_z_minus_1(3); // 4 elem
    vector<double> zp1_1 = poly_z_plus_1(1);  // 2 elem
    vector<double> zm1_2 = poly_z_minus_1(2); // 3 elem
    vector<double> zp1_2 = poly_z_plus_1(2);  // 3 elem
    vector<double> zm1_1 = poly_z_minus_1(1); // 2 elem
    vector<double> zp1_3 = poly_z_plus_1(3);  // 4 elem

    vector<double> Term1 = poly_scale(zm1_4, a4 * pow(k,4)); // a4*k^4*(z-1)^4
    vector<double> Term2 = poly_scale(poly_multiply(zm1_3, zp1_1), a3 * pow(k,3)); // a3*k^3*(z-1)^3*(z+1)
    vector<double> Term3 = poly_scale(poly_multiply(zm1_2, zp1_2), a2 * pow(k,2)); // a2*k^2*(z-1)^2*(z+1)^2
    vector<double> Term4 = poly_scale(poly_multiply(zm1_1, zp1_3), a1 * k);         // a1*k*(z-1)*(z+1)^3
    vector<double> Term5 = poly_scale(zp1_4, a0);                                   // a0*(z+1)^4

    // Nevezo polinom: összeadjuk a Term1...Term5
    d = Term1;
    d = poly_add(d, Term2);
    d = poly_add(d, Term3);
    d = poly_add(d, Term4);
    d = poly_add(d, Term5);
    // d tartalmazza a D(z) együtthatóit magasfoktól (z^4) a konstans tagig.

    // Számláló polinom: N(z) = b2*k^2*(z-1)^2*(z+1)^2, de (z-1)^2*(z+1)^2 = (z^2-1)^2 = z^4 - 2z^2 + 1.
    n = poly_scale({1.0, 0.0, -2.0, 0.0, 1.0}, b2 * pow(k,2));

    // Normalizáljuk úgy, hogy a d legmagasabb fokú tagja 1 legyen.
    double norm = d[0];
    for (auto &coef : d) coef /= norm;
    for (auto &coef : n) coef /= norm;
    return true;
}

bool create_filter_iir(vector<double>& n, vector<double>& d, filter_kind kind, filter_type type, int order, double sampling_rate, double cutoff_low, double cutoff_high)
{
    if (type == low_pass || type == high_pass)
        return create_filter_iir_butterworth(n, d, type, order, sampling_rate, cutoff_low);
    else
        return create_filter_iir_butterworth_bandpass(n, d, type, order, sampling_rate, cutoff_low, cutoff_high);
}
